"""Stage 02b -- Scopus AbstractRetrieval for each unique paper (resumable).

Adds, per paper, fields the Search API cannot give cleanly:
  abstract (full), refcount, subject areas (codes+names), author keywords,
  language, funding count, and EXACT per-author affiliation (authorgroup:
  author Scopus id -> affiliation id).

Input : unique eids from ../../03_citation_10y/output/citation_10y.csv
Output: ../output/abstracts.csv  (resumable)

Env: SCOPUS_API_KEY (required), SCOPUS_INST_TOKEN (for FULL view)

NOTE: one API call per paper. With ~216k papers this is a large, quota-heavy
run -- launch deliberately.
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
from pybliometrics.scopus import AbstractRetrieval, init

ROOT = Path(__file__).resolve().parent.parent
EID_SRC = ROOT.parent / "03_citation_10y" / "output" / "citation_10y.csv"
OUT = ROOT / "output" / "abstracts.csv"

FIELDS = [
    "eid", "abstract", "refcount", "subject_area_codes", "subject_areas",
    "authkeywords", "language", "n_funding", "author_affils", "fetch_status",
]
MAX_ATTEMPTS = 4
RETRY_BASE_SECONDS = 5
FATAL_API_ERRORS = (Scopus401Error, Scopus403Error, Scopus429Error)


def safe(obj, attr, default=None):
    try:
        return getattr(obj, attr)
    except Exception:
        return default


def row_from_abstract(eid, ab):
    sas = ab.subject_areas or []
    ag = ab.authorgroup or []
    # author Scopus id : affiliation id  (exact per-author affiliation)
    pairs = []
    for a in ag:
        auid = safe(a, "auid")
        afid = safe(a, "affiliation_id")
        if auid:
            pairs.append(f"{auid}:{afid or ''}")
    kw = ab.authkeywords
    return {
        "eid":                eid,
        "abstract":           (ab.abstract or "").replace("\n", " ").strip(),
        "refcount":           ab.refcount or "",
        "subject_area_codes": ";".join(str(s.code) for s in sas),
        "subject_areas":      ";".join(s.area for s in sas),
        "authkeywords":       ";".join(kw) if kw else "",
        "language":           ab.language or "",
        "n_funding":          len(ab.funding or []),
        "author_affils":      ";".join(pairs),
        "fetch_status":       "ok",
    }


def fetch_abstract(eid):
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            return AbstractRetrieval(eid, view="FULL", refresh=False)
        except FATAL_API_ERRORS:
            raise
        except Scopus404Error:
            return None
        except Exception:
            if attempt == MAX_ATTEMPTS:
                raise
            delay = RETRY_BASE_SECONDS * (2 ** (attempt - 1))
            print(
                f"  .. {eid}: transient failure; retry {attempt}/{MAX_ATTEMPTS} "
                f"in {delay}s",
                flush=True,
            )
            time.sleep(delay)
    return None


def read_done_eids(path):
    if not path.exists() or path.stat().st_size == 0:
        return set()
    done = set()
    for chunk in pd.read_csv(
        path,
        usecols=["eid"],
        dtype={"eid": "string"},
        chunksize=250_000,
    ):
        done.update(chunk["eid"].dropna().astype(str))
    return done


def main():
    api_key = os.environ.get("SCOPUS_API_KEY")
    if not api_key:
        raise RuntimeError("SCOPUS_API_KEY env var not set")
    tok = os.environ.get("SCOPUS_INST_TOKEN")
    init(keys=[api_key], inst_tokens=[tok] if tok else None)

    eid_df = pd.read_csv(EID_SRC, usecols=["eid"], dtype={"eid": "string"})
    eids = eid_df["eid"].dropna().str.strip()
    eids = eids[eids.ne("")].drop_duplicates().tolist()
    if not eids:
        raise RuntimeError(f"No valid EIDs found in {EID_SRC}")
    print(f"[load] {len(eids):,} unique eids", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists() and OUT.stat().st_size > 0:
        existing_fields = pd.read_csv(OUT, nrows=0).columns.tolist()
        if existing_fields != FIELDS:
            raise RuntimeError(
                f"Existing {OUT} has an incompatible schema; expected {FIELDS}"
            )
    done = read_done_eids(OUT)
    todo = [e for e in eids if e not in done]
    print(f"[scope] done={len(done):,}  todo={len(todo):,}", flush=True)

    write_header = not OUT.exists() or OUT.stat().st_size == 0
    with OUT.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            w.writeheader()
        start = time.time()
        for i, eid in enumerate(todo):
            try:
                ab = fetch_abstract(eid)
                if ab is None:
                    w.writerow({"eid": eid, "fetch_status": "not_found"})
                else:
                    w.writerow(row_from_abstract(eid, ab))
            except FATAL_API_ERRORS as e:
                print(f"  !! {eid}: {type(e).__name__}: {e}; stopping", flush=True)
                raise
            except Exception as e:
                print(f"  !! {eid}: {type(e).__name__}: {e}", flush=True)
            f.flush()
            if (i + 1) % 100 == 0:
                rate = (i + 1) / (time.time() - start)
                eta = (len(todo) - i - 1) / rate / 3600 if rate else 0
                print(f"  {i+1:,}/{len(todo):,}  rate={rate:.2f}/s  eta={eta:.1f}h",
                      flush=True)

    print(f"[done] {OUT}", flush=True)


if __name__ == "__main__":
    main()
