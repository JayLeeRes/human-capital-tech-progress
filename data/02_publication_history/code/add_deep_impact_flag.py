"""Post-process — add `is_deep_impact` column to publication_history.csv.

Matches `publication_name` (normalized: lowercase, alphanumeric only) against
journal names parsed from `deep_impact_journal_list.qmd`.
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PUB_HIST = ROOT / "output" / "publication_history.csv"
JOURNAL_LIST = ROOT / "input" / "deep_impact_journal_list.qmd"
CHUNK_SIZE = 100_000
JOURNAL_ALIASES = {
    "Aea Papers And Proceedings": [
        "American Economic Review Papers and Proceedings",
    ],
    "Annals Of The AAPSS": [
        "Annals of the American Academy of Political and Social Science",
    ],
    "Carnegie-Rochester Conference Series On Public Policy": [
        "Carnegie Rochester Confer Series on Public Policy",
    ],
    "Contemporary Accounting Research/Recherche Comptable Contemporaine": [
        "Contemporary Accounting Research",
    ],
    "Journal Of The Royal Statistical Society: Series A": [
        "Journal of the Royal Statistical Society Series A Statistics in Society",
    ],
    "Oxford Economic Papers N. S.": [
        "Oxford Economic Papers",
    ],
    "Proceedings Of The National Academy Of Sciences": [
        "Proceedings of the National Academy of Sciences of the United States of America",
    ],
}


def normalize(s):
    if not isinstance(s, str):
        return ""
    s = s.lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def parse_journal_list():
    text = JOURNAL_LIST.read_text(encoding="utf-8")
    names = []
    for line in text.splitlines():
        if not line.startswith("|"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 4:
            continue
        try:
            int(parts[0])
        except ValueError:
            continue
        name = parts[1]
        if name:
            names.append(name)
    return list(dict.fromkeys(names))


def normalize_series(series):
    return (
        series.fillna("")
        .astype("string")
        .str.lower()
        .str.strip()
        .str.replace(r"[^a-z0-9]+", " ", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def main():
    names = parse_journal_list()
    if not names:
        raise RuntimeError(f"No journal names parsed from {JOURNAL_LIST}")
    print(f"[load] {len(names)} deep impact journal names from .qmd", flush=True)
    norm_set = {normalize(n) for n in names}
    alias_norms = {
        normalize(canonical): {normalize(alias) for alias in aliases}
        for canonical, aliases in JOURNAL_ALIASES.items()
    }
    norm_set.update(alias for aliases in alias_norms.values() for alias in aliases)

    columns = pd.read_csv(PUB_HIST, nrows=0).columns
    if "publication_name" not in columns:
        raise RuntimeError(f"Missing publication_name column in {PUB_HIST}")

    tmp = PUB_HIST.with_suffix(PUB_HIST.suffix + ".tmp")
    matched = 0
    total = 0
    pub_names_seen = set()
    try:
        for chunk_no, df in enumerate(
            pd.read_csv(PUB_HIST, low_memory=False, chunksize=CHUNK_SIZE)
        ):
            norm_col = normalize_series(df["publication_name"])
            df["is_deep_impact"] = norm_col.isin(norm_set).astype("int8")
            matched += int(df["is_deep_impact"].sum())
            total += len(df)
            pub_names_seen.update(norm_col.dropna().unique())
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
        raise RuntimeError(f"{PUB_HIST} contains no rows")
    tmp.replace(PUB_HIST)
    print(f"[load] {total:,} rows in publication_history.csv", flush=True)
    print(f"[match] is_deep_impact=1: {matched:,d} / {total:,d} "
          f"({matched / total * 100:.1f}%)", flush=True)

    unmatched = [
        n for n in names
        if normalize(n) not in pub_names_seen
        and not (alias_norms.get(normalize(n), set()) & pub_names_seen)
    ]
    print(f"[diagnostic] deep_impact names absent from publication_history: {len(unmatched)}", flush=True)
    for n in unmatched[:15]:
        print(f"  - {n}", flush=True)
    if len(unmatched) > 15:
        print(f"  ... and {len(unmatched) - 15} more", flush=True)

    print(f"[done] {PUB_HIST}", flush=True)


if __name__ == "__main__":
    main()
