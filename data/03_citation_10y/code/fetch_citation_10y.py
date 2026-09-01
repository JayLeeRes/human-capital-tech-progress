"""Stage 03 — Compute 10-year citation count per paper via Scopus CitationOverview.

For each paper in publication_history.csv (T6 + non-T6, all venues):
  - Window = [publication_date, publication_date + 120 months)
  - If window end > today, truncate; the 'complete' flag captures this
  - Sum citations in window using Scopus yearly citation counts, pro-rated by
    months for the partial first/last calendar years

Inputs : ../input/publication_history.csv  (← stage 02 output, symlinked)
Outputs: ../output/citation_10y_cache.csv   (eid × year × count; resumable raw cache)
         ../output/citation_10y.csv          (eid → cited_by_10y, cited_by_10y_complete)

Env vars: SCOPUS_API_KEY (required), SCOPUS_INST_TOKEN (off-campus)

Caveats:
  - Scopus has expanded cited-reference coverage back to 1970, but pre-1996
    coverage depends on whether publishers supplied historical references. The
    `complete` flag below is time-based, not a guarantee of index completeness.
  - CitationOverview supplies calendar-year totals, not citation event dates.
    Partial first/last years are therefore estimated by day-weighting yearly totals.
"""
import csv
import os
import time
from pathlib import Path

import pandas as pd
from pybliometrics.exception import (
    Scopus401Error,
    Scopus403Error,
    Scopus404Error,
    Scopus429Error,
)
from pybliometrics.scopus import CitationOverview, init

ROOT = Path(__file__).resolve().parent.parent
INP = ROOT / "input" / "publication_history.csv"
CACHE = ROOT / "output" / "citation_10y_cache.csv"
OUT = ROOT / "output" / "citation_10y.csv"

WINDOW_MONTHS = 120
BATCH = 25
MAX_ATTEMPTS = 4
RETRY_BASE_SECONDS = 5
CACHE_CHUNK_SIZE = 500_000
TODAY = pd.Timestamp.today().normalize()
TOMORROW = TODAY + pd.Timedelta(days=1)
FATAL_API_ERRORS = (Scopus401Error, Scopus403Error, Scopus429Error)


def parse_pub_date(s):
    if pd.isna(s) or not str(s).strip():
        return None
    try:
        return pd.to_datetime(s)
    except Exception:
        return None


def compute_cited_10y(pub_date_str, yearly):
    """yearly: list of (year, count). Returns (sum, complete) or (None, None)."""
    pub_dt = parse_pub_date(pub_date_str)
    if pub_dt is None:
        return None, None
    target_end = pub_dt + pd.DateOffset(months=WINDOW_MONTHS)
    complete = 1 if target_end <= TOMORROW else 0
    end_dt = min(target_end, TOMORROW)

    total = 0.0
    for y, c in yearly:
        if y < pub_dt.year:
            continue
        elig_start = pub_dt if y == pub_dt.year else pd.Timestamp(y, 1, 1)
        elig_end = TOMORROW if y == TODAY.year else pd.Timestamp(y + 1, 1, 1)
        if elig_end <= elig_start:
            continue
        win_end = min(end_dt, elig_end)
        if win_end <= elig_start:
            continue
        win_days = (win_end - elig_start).days
        elig_days = (elig_end - elig_start).days
        if elig_days <= 0:
            continue
        total += c * (win_days / elig_days)
    return total, complete


def safe_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def _fetch_batch_once(eids, start_year, end_year):
    """Returns dict: eid -> list of (year, count).

    pybliometrics ≥4.x: identifier_type='eid' returns 404. Use 'scopus_id'
    after stripping the '2-s2.0-' prefix.
    """
    out = {eid: [] for eid in eids}
    scopus_ids = [eid.replace("2-s2.0-", "") for eid in eids]
    co = CitationOverview(
        identifier=scopus_ids,
        identifier_type="scopus_id",
        date=f"{start_year}-{end_year}",
        refresh=30,
    )
    cc = co.cc or []
    for idx, eid in enumerate(eids):
        if idx >= len(cc) or cc[idx] is None:
            continue
        for entry in cc[idx]:
            try:
                year, count = entry
                y = safe_int(year)
                c = safe_int(count)
                if y is not None and c is not None:
                    out[eid].append((y, c))
            except (TypeError, ValueError):
                continue
    return out


def fetch_batch(eids, start_year, end_year):
    """Fetch a batch with bounded retries for transient transport failures."""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            return _fetch_batch_once(eids, start_year, end_year)
        except FATAL_API_ERRORS:
            raise
        except Scopus404Error:
            return {eid: [] for eid in eids}
        except Exception:
            if attempt == MAX_ATTEMPTS:
                raise
            delay = RETRY_BASE_SECONDS * (2 ** (attempt - 1))
            print(
                f"  .. transient batch failure; retry {attempt}/{MAX_ATTEMPTS} "
                f"in {delay}s",
                flush=True,
            )
            time.sleep(delay)
    return {eid: [] for eid in eids}


def read_cached_eids(path):
    if not path.exists() or path.stat().st_size == 0:
        return set()
    done = set()
    for chunk in pd.read_csv(
        path,
        usecols=["eid"],
        dtype={"eid": "string"},
        chunksize=CACHE_CHUNK_SIZE,
    ):
        done.update(chunk["eid"].dropna().astype(str))
    return done


def compose_output(unique):
    """Compose citation totals in chunks to keep the 200MB cache memory-bounded."""
    meta = unique[["eid", "pub_date"]].copy()
    meta["pub_dt"] = pd.to_datetime(meta["pub_date"], errors="coerce")
    meta["target_end"] = meta["pub_dt"] + pd.DateOffset(months=WINDOW_MONTHS)
    meta["end_dt"] = meta["target_end"].clip(upper=TOMORROW)
    meta["complete"] = (meta["target_end"] <= TOMORROW).astype("int8")

    cached_eids = set()
    totals = pd.Series(dtype="float64")
    meta_join = meta[["eid", "pub_dt", "end_dt"]]

    if CACHE.exists() and CACHE.stat().st_size > 0:
        for cache in pd.read_csv(
            CACHE,
            dtype={"eid": "string"},
            chunksize=CACHE_CHUNK_SIZE,
        ):
            cache = cache.dropna(subset=["eid"])
            cache["eid"] = cache["eid"].astype(str)
            cached_eids.update(cache["eid"].unique())
            cache = cache[cache["year"] >= 0]
            if cache.empty:
                continue

            x = cache.merge(meta_join, on="eid", how="inner", validate="many_to_one")
            x = x[x["year"] >= x["pub_dt"].dt.year]
            if x.empty:
                continue
            year_start = pd.to_datetime(
                x["year"].astype("int64").astype(str) + "-01-01", errors="coerce"
            )
            year_end = year_start + pd.DateOffset(years=1)
            current_year = x["year"].eq(TODAY.year)
            year_end.loc[current_year] = TOMORROW

            elig_start = year_start.where(x["year"].ne(x["pub_dt"].dt.year), x["pub_dt"])
            win_end = year_end.where(year_end.le(x["end_dt"]), x["end_dt"])
            numer = (win_end - elig_start).dt.days.clip(lower=0)
            denom = (year_end - elig_start).dt.days
            weight = numer.div(denom.where(denom.gt(0))).fillna(0.0)
            x["weighted_count"] = pd.to_numeric(x["count"], errors="coerce").fillna(0) * weight
            part = x.groupby("eid", sort=False)["weighted_count"].sum()
            totals = totals.add(part, fill_value=0.0)

    out_df = meta[["eid", "complete"]].copy()
    out_df["cited_by_10y"] = out_df["eid"].map(totals)
    fetched = out_df["eid"].isin(cached_eids)
    # A successful fetch with no returned year rows is stored as a -1 sentinel.
    out_df.loc[fetched & out_df["cited_by_10y"].isna(), "cited_by_10y"] = 0.0
    out_df.loc[~fetched, "complete"] = pd.NA
    out_df = out_df[["eid", "cited_by_10y", "complete"]].rename(
        columns={"complete": "cited_by_10y_complete"}
    )
    out_df["cited_by_10y_complete"] = out_df[
        "cited_by_10y_complete"
    ].astype("Int8")

    tmp = OUT.with_suffix(OUT.suffix + ".tmp")
    out_df.to_csv(tmp, index=False)
    tmp.replace(OUT)
    return out_df


def main():
    api_key = os.environ.get("SCOPUS_API_KEY")
    if not api_key:
        raise RuntimeError("SCOPUS_API_KEY env var not set")
    inst_token = os.environ.get("SCOPUS_INST_TOKEN")
    init(keys=[api_key], inst_tokens=[inst_token] if inst_token else None)

    required = ["eid", "publication_year", "cover_date"]
    df = pd.read_csv(INP, usecols=required, low_memory=False)
    print(f"[load] publication_history rows: {len(df):,}", flush=True)

    # cover_date preferred (YYYY-MM-DD), fallback to publication_year + Jan 1
    unique = df.drop_duplicates("eid")[required].copy()
    unique = unique.dropna(subset=["eid"])
    unique["eid"] = unique["eid"].astype(str).str.strip()
    unique = unique[unique["eid"].ne("")]
    cover_date = pd.to_datetime(unique["cover_date"], errors="coerce")
    pub_year = pd.to_numeric(unique["publication_year"], errors="coerce").astype("Int64")
    year_date = pd.to_datetime(pub_year.astype("string") + "-01-01", errors="coerce")
    unique["pub_date"] = cover_date.fillna(year_date)
    unique = unique.dropna(subset=["pub_date"])
    unique["publication_year"] = unique["pub_date"].dt.year.astype("int64")
    print(f"[load] unique papers with pub_date: {len(unique):,}", flush=True)

    # Resume from cache
    done = read_cached_eids(CACHE)
    todo = unique[~unique["eid"].isin(done)].copy()
    todo = todo.sort_values(["publication_year", "eid"], kind="stable")
    print(f"[scope] cached={len(done):,}  todo={len(todo):,}", flush=True)

    if len(todo) > 0:
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        write_header = not CACHE.exists() or CACHE.stat().st_size == 0
        with CACHE.open("a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if write_header:
                w.writerow(["eid", "year", "count"])
            start_t = time.time()
            rows = todo.to_dict("records")
            for i in range(0, len(rows), BATCH):
                batch = rows[i:i + BATCH]
                eids = [r["eid"] for r in batch]
                pub_years = [int(r["publication_year"]) for r in batch]
                start_year = max(min(pub_years), 1900)
                end_year = min(max(pub_years) + 10, TODAY.year)
                if end_year < start_year:
                    end_year = start_year

                try:
                    result = fetch_batch(eids, start_year, end_year)
                    for eid in eids:
                        rows_for_eid = result.get(eid, [])
                        if not rows_for_eid:
                            w.writerow([eid, -1, 0])  # sentinel: fetched, no data
                        else:
                            for y, c in rows_for_eid:
                                w.writerow([eid, y, c])
                except FATAL_API_ERRORS as e:
                    f.flush()
                    print(f"  !! batch@{i}: {type(e).__name__}: {e}; stopping", flush=True)
                    raise
                except Exception as e:
                    # Failed batches are intentionally not cached. They remain
                    # in `todo` and will be retried on the next run.
                    print(f"  !! batch@{i}: {type(e).__name__}: {e}", flush=True)

                if ((i // BATCH) + 1) % 20 == 0 or i + BATCH >= len(rows):
                    f.flush()
                    done_n = i + len(batch)
                    elapsed = time.time() - start_t
                    rate = done_n / elapsed if elapsed > 0 else 0
                    eta = (len(rows) - done_n) / rate / 60 if rate > 0 else 0
                    print(f"  {done_n:,}/{len(rows):,}  rate={rate:.1f}/s  eta={eta:.1f}min",
                          flush=True)

    # Compose final output
    print("[compose] computing 10y window sums", flush=True)
    out_df = compose_output(unique)
    n_filled = out_df["cited_by_10y"].notna().sum()
    n_complete = (out_df["cited_by_10y_complete"] == 1).sum()
    print(f"[done] {OUT}  rows={len(out_df):,}  filled={n_filled:,}  complete={n_complete:,}",
          flush=True)


if __name__ == "__main__":
    main()
