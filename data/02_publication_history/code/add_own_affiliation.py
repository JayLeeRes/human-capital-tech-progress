"""Stage 02 post-process -- add each paper's OWN-author affiliation.

Using the Scopus COMPLETE-view per-author affiliation ids, extract, for every
author x paper row, THAT author's own institution (not coauthors'). This is the
exact per-author affiliation that replaces the single-affiliation heuristic.

Mapping per row:
  author_ids and author_afids are ';'-aligned (one entry per author; an author's
  multiple affiliation ids are '-' joined). afid / affilname / affiliation_city
  are ';'-aligned (paper's distinct affiliations). For the target author_id,
  locate its position in author_ids, take its afid(s), map to name/city.

Input/Output: ../output/publication_history.csv  (adds own_affiliation, own_city)

Run (AFTER the COMPLETE fetch finishes): python add_own_affiliation.py
"""
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "output" / "publication_history.csv"
CHUNK_SIZE = 100_000


def own_affil(target, author_ids, author_afids, afid, affilname, city):
    target = target.strip()
    aids = [x.strip() for x in author_ids.split(";")]
    aafs = [x.strip() for x in author_afids.split(";")]
    af = [x.strip() for x in afid.split(";")]
    nm = [x.strip() for x in affilname.split(";")]
    ci = [x.strip() for x in city.split(";")]
    a2n = {a: n for a, n in zip(af, nm) if a}
    a2c = {a: c for a, c in zip(af, ci) if a}
    if target not in aids:
        return "", ""
    i = aids.index(target)
    if i >= len(aafs):
        return "", ""
    my = [x.strip() for x in aafs[i].split("-") if x.strip()]
    names = "; ".join(dict.fromkeys(a2n[x] for x in my if a2n.get(x)))
    cities = "; ".join(dict.fromkeys(a2c[x] for x in my if a2c.get(x)))
    return names, cities


def main():
    need = ["author_id", "author_ids", "author_afids", "afid",
            "affilname", "affiliation_city"]
    columns = pd.read_csv(OUT, nrows=0).columns.tolist()
    miss = [c for c in need if c not in columns]
    if miss:
        raise SystemExit(f"missing columns {miss} -- run the COMPLETE fetch first")

    tmp = OUT.with_suffix(OUT.suffix + ".tmp")
    total = 0
    covered = 0
    try:
        for chunk_no, df in enumerate(
            pd.read_csv(
                OUT,
                dtype=str,
                keep_default_na=False,
                chunksize=CHUNK_SIZE,
            )
        ):
            values = [
                own_affil(*row)
                for row in df[need].itertuples(index=False, name=None)
            ]
            df["own_affiliation"] = [x[0] for x in values]
            df["own_city"] = [x[1] for x in values]
            total += len(df)
            covered += df["own_affiliation"].ne("").sum()
            df.to_csv(
                tmp,
                mode="w" if chunk_no == 0 else "a",
                header=chunk_no == 0,
                index=False,
            )
    except Exception:
        tmp.unlink(missing_ok=True)
        raise

    if total == 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{OUT} contains no rows")
    tmp.replace(OUT)
    cov = covered / total * 100
    print(f"[own_affiliation] coverage: {cov:.1f}% of {total:,} rows")
    print(f"[done] wrote own_affiliation, own_city -> {OUT}")


if __name__ == "__main__":
    main()
