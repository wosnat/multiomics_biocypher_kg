#!/usr/bin/env python3
"""Add log2(FC), normalized p-value, and DEH24 locus-tag columns to Moreno 2023 S3
Alteromonas CSVs (5 cocultures: MED4, SS120, WH7803, WH8102, BL107).

Follows the same convention as ``scripts/build_moreno2023_modified_csv.py``:
reads the original CSVs untouched and writes ``<orig>_modified.csv`` in the
same directory.

Why: the S3 sheets have per-strain column-name quirks (Spanish "Oscuridad" /
English "D" mix, and a ``FC GmM/L`` no-space typo in SS120). This script
normalises them so the paperconfig can reference a single consistent set of
column names across all 5 strains.

Per-CSV source columns (see CLAUDE comment block in the deploy task):

| Strain | t.test darkness cols            | FC darkness cols                           |
|--------|----------------------------------|--------------------------------------------|
| MED4   | t.test G nM/O, t.test G mM/O    | FC G nM/D, FC G mM/D                       |
| SS120  | t.test G nM/O, t.test G mM/O    | FC G nM/D, FC G mM/D, **FC GmM/L (typo)**  |
| WH7803 | t.test G nM/O, t.test G mM/O    | FC G nM/D, FC G mM/D                       |
| WH8102 | t.test G nM/D, t.test G mM/D    | FC G nM/D, FC G mM/D                       |
| BL107  | t.test G nM/O, t.test G mM/O    | FC G nM/D, FC G mM/D                       |

Light columns are consistent across all: ``FC G nM/L``, ``FC G mM/L``,
``t.test G nM/L``, ``t.test G mM/L`` — except SS120's ``FC GmM/L``.

Output columns added to every S3 CSV:
  locus_tag_deh24                     (regex-extracted from Description)
  log2_FC_light_low_glucose           (log2 of FC G nM/L)
  log2_FC_light_high_glucose          (log2 of FC G mM/L  OR FC GmM/L for SS120)
  log2_FC_dark_low_glucose            (log2 of FC G nM/D)
  log2_FC_dark_high_glucose           (log2 of FC G mM/D)
  t_test_light_low_glucose            (verbatim from t.test G nM/L)
  t_test_light_high_glucose           (verbatim from t.test G mM/L)
  t_test_dark_low_glucose             (verbatim from t.test G nM/O  OR t.test G nM/D for WH8102)
  t_test_dark_high_glucose            (verbatim from t.test G mM/O  OR t.test G mM/D for WH8102)
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/moreno 2023")
S3_TEMPLATE = "table s3 Alteromonas in {strain} cultures spectrum.03275-22-s0003.csv"
DEH24_RE = re.compile(r"DEH24_\d+")

# Per-strain column mapping: new_col -> source_col in raw CSV.
# If a source column is missing we emit NaN for that new column and warn.
PER_STRAIN_FC_COLS: dict[str, dict[str, str]] = {
    "MED4": {
        "log2_FC_light_low_glucose": "FC G nM/L",
        "log2_FC_light_high_glucose": "FC G mM/L",
        "log2_FC_dark_low_glucose": "FC G nM/D",
        "log2_FC_dark_high_glucose": "FC G mM/D",
    },
    "SS120": {
        "log2_FC_light_low_glucose": "FC G nM/L",
        "log2_FC_light_high_glucose": "FC GmM/L",  # no-space typo in source
        "log2_FC_dark_low_glucose": "FC G nM/D",
        "log2_FC_dark_high_glucose": "FC G mM/D",
    },
    "WH7803": {
        "log2_FC_light_low_glucose": "FC G nM/L",
        "log2_FC_light_high_glucose": "FC G mM/L",
        "log2_FC_dark_low_glucose": "FC G nM/D",
        "log2_FC_dark_high_glucose": "FC G mM/D",
    },
    "WH8102": {
        "log2_FC_light_low_glucose": "FC G nM/L",
        "log2_FC_light_high_glucose": "FC G mM/L",
        "log2_FC_dark_low_glucose": "FC G nM/D",
        "log2_FC_dark_high_glucose": "FC G mM/D",
    },
    "BL107": {
        "log2_FC_light_low_glucose": "FC G nM/L",
        "log2_FC_light_high_glucose": "FC G mM/L",
        "log2_FC_dark_low_glucose": "FC G nM/D",
        "log2_FC_dark_high_glucose": "FC G mM/D",
    },
}

PER_STRAIN_TTEST_COLS: dict[str, dict[str, str]] = {
    "MED4": {
        "t_test_light_low_glucose": "t.test G nM/L",
        "t_test_light_high_glucose": "t.test G mM/L",
        "t_test_dark_low_glucose": "t.test G nM/O",
        "t_test_dark_high_glucose": "t.test G mM/O",
    },
    "SS120": {
        "t_test_light_low_glucose": "t.test G nM/L",
        "t_test_light_high_glucose": "t.test G mM/L",
        "t_test_dark_low_glucose": "t.test G nM/O",
        "t_test_dark_high_glucose": "t.test G mM/O",
    },
    "WH7803": {
        "t_test_light_low_glucose": "t.test G nM/L",
        "t_test_light_high_glucose": "t.test G mM/L",
        "t_test_dark_low_glucose": "t.test G nM/O",
        "t_test_dark_high_glucose": "t.test G mM/O",
    },
    "WH8102": {
        "t_test_light_low_glucose": "t.test G nM/L",
        "t_test_light_high_glucose": "t.test G mM/L",
        "t_test_dark_low_glucose": "t.test G nM/D",
        "t_test_dark_high_glucose": "t.test G mM/D",
    },
    "BL107": {
        "t_test_light_low_glucose": "t.test G nM/L",
        "t_test_light_high_glucose": "t.test G mM/L",
        "t_test_dark_low_glucose": "t.test G nM/O",
        "t_test_dark_high_glucose": "t.test G mM/O",
    },
}


def _log2(x) -> float:
    try:
        v = float(x)
    except Exception:
        return float("nan")
    if pd.isna(v) or v <= 0:
        return float("nan")
    return math.log2(v)


def _extract_deh24(desc) -> str:
    if not isinstance(desc, str):
        return ""
    m = DEH24_RE.search(desc)
    return m.group(0) if m else ""


def process(strain: str) -> None:
    src = PAPER_DIR / S3_TEMPLATE.format(strain=strain)
    if not src.exists():
        print(f"  MISSING: {src}")
        return
    df = pd.read_csv(src)

    # DEH24 locus tag column
    if "Description" in df.columns:
        df["locus_tag_deh24"] = df["Description"].apply(_extract_deh24)
    else:
        print(f"  WARNING: {strain}: no Description column; locus_tag_deh24 left empty")
        df["locus_tag_deh24"] = ""

    # log2(FC) columns
    fc_map = PER_STRAIN_FC_COLS[strain]
    for new_col, src_col in fc_map.items():
        if src_col in df.columns:
            df[new_col] = df[src_col].apply(_log2)
        else:
            print(f"  WARNING: {strain}: missing source FC column {src_col!r} for {new_col}")
            df[new_col] = float("nan")

    # t-test p-value columns (verbatim)
    ttest_map = PER_STRAIN_TTEST_COLS[strain]
    for new_col, src_col in ttest_map.items():
        if src_col in df.columns:
            df[new_col] = df[src_col]
        else:
            print(f"  WARNING: {strain}: missing source t-test column {src_col!r} for {new_col}")
            df[new_col] = float("nan")

    out = src.with_name(src.stem + "_modified.csv")
    df.to_csv(out, index=False)

    n_deh24 = (df["locus_tag_deh24"] != "").sum()
    print(f"{strain}: {len(df)} rows -> {out.name}")
    print(f"  locus_tag_deh24: {n_deh24}/{len(df)} extracted")
    for new_col in fc_map:
        vals = df[new_col].dropna()
        if len(vals):
            print(f"  {new_col}: {vals.min():.3f} to {vals.max():.3f} (n={len(vals)})")
        else:
            print(f"  {new_col}: no numeric values")


def main() -> None:
    print("=== Alteromonas S3 ===")
    for strain in ["MED4", "SS120", "WH7803", "WH8102", "BL107"]:
        process(strain)


if __name__ == "__main__":
    main()
