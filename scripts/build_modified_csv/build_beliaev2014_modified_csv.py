"""Beliaev 2014 Tables S3-S8: re-derive log2 fold change from RPKM with a one-read pseudocount.

Each table has three RPKM columns (axenic control, co-culture lactate, co-culture HCO3-)
and three published fold-change columns that are plain RPKM ratios. A gene with zero
reads is written at a floor (Synechococcus RPKM 0.12, Shewanella 0), so a ratio against
it is meaningless: floor-driven ratios take every one of the top-50 effects in S3, S6, S7
and S8 (up to 498-fold), and 97 genes that drop to zero print as 0.00 and could not be
converted at all.

This script writes `<table>_modified.csv` with

    log2_fold_change = log2((num + r1) / (den + r1))

where a floor value counts as 0 and r1 is the table's one-read RPKM (the smallest RPKM
above the floor: 0.54 for Synechococcus, 0.8 for Shewanella). Rows with both arms at the
floor get a blank value. The tables are the paper's own >= 2-fold calls, so the
paperconfig keeps `prefiltered: true`; only magnitudes change.

Before writing, it checks that the published fold column equals num / den (within RPKM
rounding) on rows off the floor, so a misread column order fails loudly.

Run: uv run python scripts/build_modified_csv/build_beliaev2014_modified_csv.py
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Beliaev 2014")
SKIP_ROWS = 2
# positional RPKM columns: 2 = axenic control, 3 = co-culture lactate, 4 = co-culture HCO3-
# positional fold columns: 5 = lactate / axenic, 6 = HCO3- / axenic, 7 = HCO3- / lactate
AX, LAC, HCO3 = 2, 3, 4
TABLES = {  # table -> (numerator, denominator, published fold column)
    "s3": (LAC, AX, 5), "s4": (LAC, AX, 5),
    "s5": (HCO3, AX, 6), "s6": (HCO3, AX, 6),
    "s7": (HCO3, LAC, 7), "s8": (HCO3, LAC, 7),
}


def build(table: str, num_col: int, den_col: int, fold_col: int) -> None:
    src = PAPER_DIR / f"table {table} 41396_2014_bfismej201469_moesm40_esm.csv"
    dst = src.with_name(src.stem + "_modified.csv")
    title = src.read_text(encoding="utf-8-sig").splitlines()[:SKIP_ROWS]
    df = pd.read_csv(src, skiprows=SKIP_ROWS, encoding="utf-8-sig")

    rpkm = df.iloc[:, [AX, LAC, HCO3]].apply(pd.to_numeric, errors="coerce")
    values = np.unique(rpkm.values[~np.isnan(rpkm.values)])
    floor, one_read = float(values[0]), float(values[1])

    num = pd.to_numeric(df.iloc[:, num_col], errors="coerce")
    den = pd.to_numeric(df.iloc[:, den_col], errors="coerce")
    published = pd.to_numeric(df.iloc[:, fold_col], errors="coerce")

    off_floor = (num > floor) & (den > floor) & published.notna() & (published > 0)
    rel = ((published - num / den).abs() / published)[off_floor]
    if rel.median() > 0.02:
        raise SystemExit(f"{table}: published fold != num/den (median rel err {rel.median():.3f}) — columns misread?")

    n0 = num.where(num > floor, 0.0)
    d0 = den.where(den > floor, 0.0)
    both_zero = (n0 == 0) & (d0 == 0)
    log2fc = [None if (bz or math.isnan(a) or math.isnan(b)) else round(math.log2((a + one_read) / (b + one_read)), 6)
              for a, b, bz in zip(n0, d0, both_zero)]
    df["log2_fold_change"] = log2fc

    with dst.open("w", encoding="utf-8", newline="") as f:
        for line in title:
            f.write(line + "\n")
        df.to_csv(f, index=False)

    lf = pd.Series([v for v in log2fc if v is not None])
    floored = int(((num <= floor) | (den <= floor)).sum())
    pub = np.log2(published[published > 0])
    print(f"{dst.name}: {len(lf)} rows (floor {floor}, one read {one_read}; {floored} rows with an arm at the "
          f"floor, {int(both_zero.sum())} both-floor blanked); log2 {lf.min():.2f} .. {lf.max():.2f} "
          f"(published {pub.min():.2f} .. {pub.max():.2f}); |log2| < 1 now: {int((lf.abs() < 1).sum())}")


def main() -> None:
    for table, (num_col, den_col, fold_col) in TABLES.items():
        build(table, num_col, den_col, fold_col)


if __name__ == "__main__":
    main()
