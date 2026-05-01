#!/usr/bin/env python3
"""Strip percent signs and add log2(MV/Cell) enrichment to Fadeev 2022 strain CSVs.

Source CSVs (one per strain) have:
  - Header row 0: NCBI_PGAP_function, NCBI_PGAP_accession, COG20_CATEGORY_accession,
    COG20_CATEGORY_function, Prop. Abund. MVs, Prop. Abund. Cells
  - The two abundance columns are formatted as percent strings (e.g. "13.3%"); the
    `observations_adapter` calls float() on each cell, which would skip every row.

Output: `<orig stem>_modified.csv` with three numeric columns added:
  - prop_abund_mvs_percent   : Prop. Abund. MVs without the % sign
  - prop_abund_cells_percent : Prop. Abund. Cells without the % sign
  - log2_mv_cell_enrichment  : log2(MVs / Cells); NaN unless both > 0
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/fadeev 2022")
STRAIN_CSVS = [
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - AD45.csv",
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - ATCC27126.csv",
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - BGP6.csv",
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - BS11.csv",
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - HOT1A3.csv",
    "Table_S1-most_abundant_MV_proteins_Supplementary_Data - MIT1002.csv",
]


def _strip_percent(value) -> float:
    if pd.isna(value):
        return float("nan")
    s = str(value).strip()
    if not s:
        return float("nan")
    if s.endswith("%"):
        s = s[:-1].strip()
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _log2_ratio(mv: float, cell: float) -> float:
    if pd.isna(mv) or pd.isna(cell):
        return float("nan")
    if mv <= 0 or cell <= 0:
        return float("nan")
    return math.log2(mv / cell)


def process(csv_name: str) -> None:
    src = PAPER_DIR / csv_name
    df = pd.read_csv(src)

    required = {"Prop. Abund. MVs", "Prop. Abund. Cells", "NCBI_PGAP_accession"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_name}: missing columns {sorted(missing)}")

    df["prop_abund_mvs_percent"] = df["Prop. Abund. MVs"].apply(_strip_percent)
    df["prop_abund_cells_percent"] = df["Prop. Abund. Cells"].apply(_strip_percent)
    df["log2_mv_cell_enrichment"] = [
        _log2_ratio(mv, cell)
        for mv, cell in zip(df["prop_abund_mvs_percent"], df["prop_abund_cells_percent"])
    ]

    df = df.dropna(how="all")
    out = src.with_name(src.stem + "_modified.csv")
    df.to_csv(out, index=False)

    n_rows = len(df)
    n_mv = df["prop_abund_mvs_percent"].notna().sum()
    n_cell = df["prop_abund_cells_percent"].notna().sum()
    n_log2 = df["log2_mv_cell_enrichment"].notna().sum()
    print(
        f"{csv_name}: {n_rows} rows -> {out.name} "
        f"(mv={n_mv}, cell={n_cell}, log2={n_log2})"
    )


def main() -> None:
    for name in STRAIN_CSVS:
        process(name)


if __name__ == "__main__":
    main()
