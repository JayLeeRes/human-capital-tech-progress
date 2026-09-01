"""Stage 02 — Fetch full Scopus publication history for every T6 author.

Inputs : ../input/t6_authors.csv  (← stage 01 output)
Outputs: ../output/publication_history.csv  (author × work, resumable)

Scopus query: AU-ID(<author_id>)
Env vars: SCOPUS_API_KEY (required), SCOPUS_INST_TOKEN (off-campus only)
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
from pybliometrics.scopus import ScopusSearch, init

ROOT = Path(__file__).resolve().parent.parent
INP = ROOT / "input" / "t6_authors.csv"
OUT = ROOT / "output" / "publication_history.csv"
STATE = ROOT / "output" / "publication_history_completed_authors.csv"

FIELDS = [
    "author_id", "eid", "doi", "title",
    "publication_year", "cover_date", "subtype",
    "publication_name", "issn", "source_id",
    "cited_by_count", "author_count", "creator",
    "author_ids", "author_afids",
    "afid", "affilname", "affiliation_city", "affiliation_country",
]
MAX_ATTEMPTS = 4
RETRY_BASE_SECONDS = 5
FATAL_API_ERRORS = (Scopus401Error, Scopus403Error, Scopus429Error)


def row_from_result(author_id, r):
    cover = r.coverDate or ""
    return {
        "author_id":           author_id,
        "eid":                 r.eid,
        "doi":                 r.doi or "",
        "title":               (r.title or "").replace("\n", " ").strip(),
        "publication_year":    cover[:4],
        "cover_date":          cover,
        "subtype":             r.subtypeDescription or r.subtype or "",
        "publication_name":    r.publicationName or "",
        "issn":                r.issn or "",
        "source_id":           r.source_id or "",
        "cited_by_count":      r.citedby_count or 0,
        "author_count":        r.author_count or 0,
        "creator":             r.creator or "",
        # author_ids and author_afids are ';'-aligned -> per-author affiliation;
        # map each author's afid via the paper's afid/affilname/affiliation_city.
        "author_ids":          r.author_ids or "",
        "author_afids":        r.author_afids or "",
        "afid":                r.afid or "",
        "affilname":           r.affilname or "",
        "affiliation_city":    r.affiliation_city or "",
        "affiliation_country": r.affiliation_country or "",
    }


def search_author(author_id):
    """Query one author, retrying transient failures but not quota/auth errors."""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            search = ScopusSearch(
                f"AU-ID({author_id})", view="COMPLETE", refresh=30, verbose=False
            )
            return search.results or []
        except FATAL_API_ERRORS:
            raise
        except Scopus404Error:
            return []
        except Exception:
            if attempt == MAX_ATTEMPTS:
                raise
            delay = RETRY_BASE_SECONDS * (2 ** (attempt - 1))
            print(
                f"  .. {author_id}: transient failure; retry {attempt}/{MAX_ATTEMPTS} "
                f"in {delay}s",
                flush=True,
            )
            time.sleep(delay)
    return []


def read_completed_authors(path):
    """Read the existing author column in bounded memory for resume support."""
    if not path.exists() or path.stat().st_size == 0:
        return set()
    done = set()
    for chunk in pd.read_csv(
        path,
        usecols=["author_id"],
        dtype={"author_id": "string"},
        on_bad_lines="error",
        chunksize=250_000,
    ):
        done.update(chunk["author_id"].dropna().astype(str))
    return done


def initialize_state(existing_done, reset=False):
    """Create/migrate the explicit completion ledger used by future resumes."""
    if not reset and STATE.exists() and STATE.stat().st_size > 0:
        state = pd.read_csv(STATE, dtype={"author_id": "string"})
        if "author_id" not in state.columns:
            raise RuntimeError(f"Missing author_id column in {STATE}")
        return set(state["author_id"].dropna().astype(str))

    tmp = STATE.with_suffix(STATE.suffix + ".tmp")
    pd.DataFrame({"author_id": sorted(existing_done)}).to_csv(tmp, index=False)
    tmp.replace(STATE)
    return existing_done


def main():
    api_key = os.environ.get("SCOPUS_API_KEY")
    if not api_key:
        raise RuntimeError("SCOPUS_API_KEY env var not set")
    inst_token = os.environ.get("SCOPUS_INST_TOKEN")
    init(keys=[api_key], inst_tokens=[inst_token] if inst_token else None)

    authors_df = pd.read_csv(INP, usecols=["author_id"], dtype={"author_id": "string"})
    authors = authors_df["author_id"].dropna().str.strip()
    authors = authors[authors.ne("")].drop_duplicates().tolist()
    if not authors:
        raise RuntimeError(f"No valid author_id values found in {INP}")
    print(f"[load] {len(authors):,} unique T6 authors", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    inferred_done = read_completed_authors(OUT)
    done = initialize_state(inferred_done, reset=not OUT.exists())
    todo = [a for a in authors if a not in done]
    print(f"[scope] done={len(done):,}  todo={len(todo):,}", flush=True)

    write_header = not OUT.exists() or OUT.stat().st_size == 0
    state_header = not STATE.exists() or STATE.stat().st_size == 0
    with OUT.open("a", newline="", encoding="utf-8") as f, STATE.open(
        "a", newline="", encoding="utf-8"
    ) as state_f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        state_w = csv.writer(state_f)
        if write_header:
            w.writeheader()
        if state_header:
            state_w.writerow(["author_id"])
        start = time.time()
        for i, aid in enumerate(todo):
            try:
                # Materialize one author's rows before writing. This sharply
                # reduces the chance that an interrupted run leaves a partial
                # author that is then mistaken for complete on resume.
                author_rows = [row_from_result(aid, r) for r in search_author(aid)]
                w.writerows(author_rows)
                f.flush()
                state_w.writerow([aid])
                state_f.flush()
            except FATAL_API_ERRORS as e:
                print(f"  !! {aid}: {type(e).__name__}: {e}; stopping", flush=True)
                raise
            except Exception as e:
                print(f"  !! {aid}: {type(e).__name__}: {e}", flush=True)
            if (i + 1) % 50 == 0:
                rate = (i + 1) / (time.time() - start)
                eta = (len(todo) - i - 1) / rate / 60 if rate else 0
                print(f"  {i+1:,}/{len(todo):,}  rate={rate:.2f}/s  eta={eta:.1f}min", flush=True)

    print(f"[done] {OUT}", flush=True)


if __name__ == "__main__":
    main()
