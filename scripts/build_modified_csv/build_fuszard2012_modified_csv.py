#!/usr/bin/env python3
"""Add `log2_fold_change = log2(Deplete Mean / Control Mean)` to each Fuszard 2012 strain CSV.

The published "Ratio of means" equals Deplete / Control for MIT9312 and SS120, but for
NATL2A it is Deplete Mean / `# peptides` on every row (a spreadsheet error in the
supplement): 62 of 63 NATL2A rows read below 1, where Deplete / Control gives the
paper's 6 up and 1 down past its 1.6 / 0.6 cutoff. So the ratio is always recomputed
from the two mean columns; the script reports, per strain, how many published ratios
agree with the recomputation.

Source CSVs (one per strain) have:
  - row 0: title line (e.g., "SS120 Quantitation data")
  - row 1: actual header (Protein, Gene, ORF, # peptides, Control Mean, Deplete Mean, Deplete SD, Ratio of means)
  - rows 2+: data

Output: `<orig stem>_modified.csv` with proper header as row 0 and a new `log2_fold_change` column.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Fuszard 2012")
CSVS = [
    "table s1 MIT9312 Quantitation data.csv",
    "table s1 NATL2A Quantitation data.csv",
    "table s1 SS120 Quantitation data.csv",
]


def process(csv_name: str) -> None:
    src = PAPER_DIR / csv_name
    # Skip the first row (strain title) so pandas reads the real header
    df = pd.read_csv(src, skiprows=1)

    if "Ratio of means" not in df.columns:
        raise ValueError(f"Column 'Ratio of means' missing in {csv_name}. Columns: {df.columns.tolist()}")
    # Header case differs by strain ("Control MEAN" / "deplete MEAN" in MIT9312)
    control = next(c for c in df.columns if c.strip().lower() == "control mean")
    deplete = next(c for c in df.columns if c.strip().lower() == "deplete mean")
    ratio = pd.to_numeric(df[deplete], errors="coerce") / pd.to_numeric(df[control], errors="coerce")
    published = pd.to_numeric(df["Ratio of means"], errors="coerce")
    has = ratio.notna() & published.notna()
    agree = int(((published - ratio).abs() / ratio < 1e-3)[has].sum())

    def _log2(x) -> float:
        try:
            v = float(x)
        except Exception:
            return float("nan")
        if pd.isna(v) or v <= 0:
            return float("nan")
        return math.log2(v)

    df["log2_fold_change"] = ratio.apply(_log2)

    # Drop completely empty rows (trailing)
    df = df.dropna(how="all")

    out = src.with_name(src.stem + "_modified.csv")
    df.to_csv(out, index=False)
    print(f"{csv_name}: {len(df)} rows -> {out.name} "
          f"(published ratio == Deplete/Control on {agree} of {int(has.sum())})")
    print(f"  log2FC range: {df['log2_fold_change'].min():.3f} to {df['log2_fold_change'].max():.3f}")


def main() -> None:
    for name in CSVS:
        process(name)


if __name__ == "__main__":
    main()
