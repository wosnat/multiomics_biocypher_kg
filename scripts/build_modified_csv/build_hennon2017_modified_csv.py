#!/usr/bin/env python3
"""Reshape Hennon 2017 S5 (long: gene x category) to wide: 1 row per gene + 9 boolean Y/blank columns.

S5 lists 513 (gene, category) rows over 452 unique AEZ55 genes; ~60 genes
appear in two enrichment categories. Output: 452 rows, one per gene, with:
- DE stats columns (logFC, logCPM, LR, PValue, annot) — same across category-duplicate rows
- 9 boolean columns (`Y` if assigned, blank otherwise), one per category

Output: `ALT gene categories ...moesm31_esm_modified.csv` alongside the original.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Hennon 2017")
SRC_CSV = PAPER_DIR / "ALT gene categories for enrichment 41396_2018_bfismej2017189_moesm31_esm.csv"
OUT_CSV = PAPER_DIR / "ALT gene categories for enrichment 41396_2018_bfismej2017189_moesm31_esm_modified.csv"


def slugify(label: str) -> str:
    s = label.lower().strip()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s


def main() -> None:
    df = pd.read_csv(SRC_CSV, skiprows=2)
    assert "enrichment category" in df.columns, df.columns.tolist()

    categories = sorted(df["enrichment category"].dropna().unique())
    cat_cols = {cat: slugify(cat) for cat in categories}
    print(f"Categories ({len(categories)}): {list(cat_cols.values())}")

    de_cols = ["geneID", "Path", "KO", "chr", "start", "end", "strand",
               "logFC", "logCPM", "LR", "PValue", "annot"]
    # DE stats are identical across a gene's category rows; first() is safe
    base = df.groupby("geneID", as_index=False)[de_cols].first()

    for cat, col in cat_cols.items():
        gene_set = set(df.loc[df["enrichment category"] == cat, "geneID"])
        base[col] = base["geneID"].map(lambda g: "Y" if g in gene_set else "")

    # Sanity
    assert len(base) == df["geneID"].nunique() == 452, len(base)
    n_double = (base[list(cat_cols.values())].eq("Y").sum(axis=1) >= 2).sum()
    print(f"Wrote {len(base)} rows; {n_double} genes carry >=2 category labels")

    base.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
