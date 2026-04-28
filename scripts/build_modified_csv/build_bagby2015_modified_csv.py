#!/usr/bin/env python3
"""Bagby 2015 Table S2 (moesm34) — clean view for `derived_metrics_table` ingest.

The original supp CSV has:
  - 9 lines of title + footnotes before the header row
  - Multi-line quoted column headers with embedded newlines and en-dashes
  - UTF-8 BOM
  - Trailing trailing-empty unnamed columns

The DerivedMetric adapter calls plain `pd.read_csv(...)` and the column names
must be safe to reference from YAML, so this script writes a normalized
sibling `..._modified.csv` with: BOM stripped, preamble dropped, ASCII column
names, and no trailing all-empty columns. Original file is left untouched.

Re-run is idempotent.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/bagby 2015")
SRC = PAPER_DIR / "41396_2015_bfismej201536_moesm34_esm.csv"
DST = PAPER_DIR / "41396_2015_bfismej201536_moesm34_esm_modified.csv"

# 9 preamble rows (title + 5 numbered notes) precede the multi-line header row.
SKIP_ROWS = 9

# Map raw (post-newline-collapse) column names → normalized snake_case names.
RENAME = {
    "Locus tag1": "locus_tag",
    "Alternative locus tag": "alternative_locus_tag",
    "Cyanobase category2": "cyanobase_category",
    "Expression change under  high light shock3": "expression_change_high_light_shock",
    "Expression change under  -CO2 shock4": "expression_change_low_co2_shock",
    "Rapid recovery under  -CO2 shock5": "rapid_recovery_low_co2_shock",
    "Expression change under  -CO2/-O2 shock4": "expression_change_low_co2_no_o2_shock",
    "Rapid recovery under  -CO2/-O2 shock5": "rapid_recovery_low_co2_no_o2_shock",
}


def main() -> None:
    df = pd.read_csv(
        SRC,
        skiprows=SKIP_ROWS,
        header=0,
        engine="python",
        encoding="utf-8-sig",
    )
    # Collapse embedded newlines and replace Unicode en-dash with ASCII '-'.
    df.columns = [c.replace("\n", " ").replace("–", "-").strip() for c in df.columns]
    # Drop trailing empty columns produced by stray commas in the preamble.
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

    missing = [c for c in RENAME if c not in df.columns]
    if missing:
        raise SystemExit(f"Expected columns missing from {SRC.name}: {missing}\nGot: {list(df.columns)}")

    df = df.rename(columns=RENAME)
    df = df[list(RENAME.values())]
    df.to_csv(DST, index=False)
    print(f"{DST.name}: {len(df)} rows, columns={list(df.columns)}")


if __name__ == "__main__":
    main()
