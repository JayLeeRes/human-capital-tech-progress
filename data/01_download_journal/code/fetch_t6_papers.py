"""Stage 01 — Fetch T6 journal papers (1970+) via Scopus and extract unique authors.

Inputs : none (T6 ISSNs hard-coded)
Outputs: ../output/t6_papers.csv   — paper × T6 (one row per paper)
         ../output/t6_authors.csv  — unique authors who appeared on ≥1 T6 paper

Scopus query: ISSN(X) AND PUBYEAR > 1969 AND (DOCTYPE(ar) OR DOCTYPE(re))
Env vars: SCOPUS_API_KEY (required), SCOPUS_INST_TOKEN (off-campus only)
"""
import os
from pathlib import Path

import pandas as pd
from pybliometrics.scopus import ScopusSearch, init

ROOT = Path(__file__).resolve().parent.parent
OUT_PAPERS = ROOT / "output" / "t6_papers.csv"
OUT_AUTHORS = ROOT / "output" / "t6_authors.csv"

T6 = {
    "AER":          "0002-8282",
    "JPE":          "0022-3808",
    "Econometrica": "0012-9682",
    "QJE":          "0033-5533",
    "RES":          "0034-6527",
    "REStat":       "0034-6535",
}
YEAR_START = 1970

PAPER_FIELDS = [
    "short_name", "issn", "eid", "doi", "title", "publication_year",
    "cover_date", "subtype", "publication_name", "source_id",
    "cited_by_count", "author_count", "author_ids", "author_names",
    "affilname", "affiliation_country",
]


def write_csv_atomic(df, path):
    """Write a CSV without exposing a partially-written final file."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(path)


def main():
    api_key = os.environ.get("SCOPUS_API_KEY")
    if not api_key:
        raise RuntimeError("SCOPUS_API_KEY env var not set")
    inst_token = os.environ.get("SCOPUS_INST_TOKEN")
    init(keys=[api_key], inst_tokens=[inst_token] if inst_token else None)

    OUT_PAPERS.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for short, issn in T6.items():
        query = f"ISSN({issn}) AND PUBYEAR > {YEAR_START - 1} AND (DOCTYPE(ar) OR DOCTYPE(re))"
        print(f"[{short}] {query}", flush=True)
        s = ScopusSearch(query, refresh=30, verbose=True)
        results = s.results or []
        for r in results:
            cover = r.coverDate or ""
            rows.append({
                "short_name":         short,
                "issn":               issn,
                "eid":                r.eid,
                "doi":                r.doi or "",
                "title":              (r.title or "").replace("\n", " ").strip(),
                "publication_year":   cover[:4],
                "cover_date":         cover,
                "subtype":            r.subtypeDescription or r.subtype or "",
                "publication_name":   r.publicationName or "",
                "source_id":          r.source_id or "",
                "cited_by_count":     r.citedby_count or 0,
                "author_count":       r.author_count or 0,
                "author_ids":         r.author_ids or "",       # ';'-separated
                "author_names":       r.author_names or "",     # ';'-separated
                "affilname":          r.affilname or "",
                "affiliation_country": r.affiliation_country or "",
            })
        print(f"[{short}] {len(results):,} papers", flush=True)

    df = pd.DataFrame(rows, columns=PAPER_FIELDS)
    if df.empty:
        raise RuntimeError("Scopus returned no T6 papers; refusing to overwrite outputs")
    if df["eid"].isna().any() or df["eid"].duplicated().any():
        raise RuntimeError("T6 paper EIDs are missing or duplicated")
    df = df.sort_values(["publication_year", "eid"], kind="stable")
    write_csv_atomic(df, OUT_PAPERS)
    print(f"[done] {OUT_PAPERS}  rows={len(df):,}", flush=True)

    # Unique authors who appeared on ≥1 T6 paper (with first T6 venue + year)
    auth_rows = []
    for r in df.itertuples(index=False):
        ids = str(r.author_ids or "").split(";")
        names = str(r.author_names or "").split(";")
        for i, aid in enumerate(ids):
            aid = aid.strip()
            if not aid:
                continue
            name = names[i].strip() if i < len(names) else ""
            try:
                year = int(r.publication_year) if r.publication_year else 9999
            except (TypeError, ValueError):
                year = 9999
            auth_rows.append({
                "author_id": aid,
                "author_name": name,
                "first_t6": r.short_name,
                "first_year": year,
            })
    aa = pd.DataFrame(auth_rows)
    if aa.empty:
        raise RuntimeError("No valid author IDs were returned")
    aa = (
        aa.sort_values(["author_id", "first_year"], kind="stable")
        .drop_duplicates(subset="author_id", keep="first")
    )
    write_csv_atomic(aa, OUT_AUTHORS)
    print(f"[done] {OUT_AUTHORS}  unique authors={len(aa):,}", flush=True)


if __name__ == "__main__":
    main()
