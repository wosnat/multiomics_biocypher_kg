#!/usr/bin/env python3
"""Add log2(FC) columns to Moreno 2023 cyano S2 + Marinobacter S4 CSVs.

S3 (Alteromonas) CSVs are intentionally SKIPPED in this batch — the DEH24 → MADE_RS
id bridge is deferred to a follow-up PR.

S2 (cyano) columns:
  FC.L.LGnM, FC.L.LGmM, FC.D.LGnM, FC.D.LGmM
S4 (Marinobacter) columns (linear ratios, named differently):
  FC G nM/L, FC G mM/L, FC G nM/D, FC G mM/D
  (also TTEST columns: t.test G nM/L, t.test G mM/L, t.test G nM/O, t.test G mM/O
   — the 'O' is a typo for 'D' in the sheet author's notation; paired with G nM/D / G mM/D FC)

Output log2 column names (new):
  cyano:   log2FC.L.LGnM, log2FC.L.LGmM, log2FC.D.LGnM, log2FC.D.LGmM
  marino:  log2FC.GnM.L, log2FC.GmM.L, log2FC.GnM.D, log2FC.GmM.D
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/moreno 2023")

CYANO_STRAINS = ["MED4", "SS120", "WH7803", "WH8102", "BL107"]
CYANO_TEMPLATE = "table s2 {strain} spectrum.03275-22-s0002.csv"
CYANO_FC_COLS = {
    "FC.L.LGnM": "log2FC.L.LGnM",
    "FC.L.LGmM": "log2FC.L.LGmM",
    "FC.D.LGnM": "log2FC.D.LGnM",
    "FC.D.LGmM": "log2FC.D.LGmM",
}

MARINO_STRAINS = ["MED4", "SS120", "WH7803", "WH8102", "BL107"]
MARINO_TEMPLATE = "table s4 Marinobacter in {strain} cultures spectrum.03275-22-s0004.csv"
MARINO_FC_COLS = {
    "FC G nM/L": "log2FC.GnM.L",
    "FC G mM/L": "log2FC.GmM.L",
    "FC G nM/D": "log2FC.GnM.D",
    "FC G mM/D": "log2FC.GmM.D",
}
# Per-strain column-name overrides — the SS120 sheet has `FC GmM/L` (no space)
# instead of `FC G mM/L`. Without this override the light_high_glucose analysis
# is silently dropped for SS120 Marinobacter. Confirmed by inspecting the original
# CSV header: `FC G nM/L,FC GmM/L,FC G nM/D,FC G mM/D,...`.
MARINO_FC_OVERRIDES_BY_STRAIN: dict[str, dict[str, str]] = {
    "SS120": {
        "FC GmM/L": "log2FC.GmM.L",  # missing the space — sheet typo
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


def process(src: Path, fc_map: dict[str, str]) -> None:
    if not src.exists():
        print(f"  MISSING: {src}")
        return
    df = pd.read_csv(src)
    missing = [c for c in fc_map if c not in df.columns]
    if missing:
        print(f"  WARNING: missing FC columns in {src.name}: {missing}")
    for orig, new in fc_map.items():
        if orig in df.columns:
            df[new] = df[orig].apply(_log2)
    out = src.with_name(src.stem + "_modified.csv")
    df.to_csv(out, index=False)
    # Report
    print(f"{src.name}: {len(df)} rows -> {out.name}")
    for new in fc_map.values():
        if new in df.columns:
            mn, mx = df[new].min(), df[new].max()
            print(f"  {new}: {mn:.3f} to {mx:.3f}")


def main() -> None:
    print("=== Cyano S2 ===")
    for strain in CYANO_STRAINS:
        src = PAPER_DIR / CYANO_TEMPLATE.format(strain=strain)
        process(src, CYANO_FC_COLS)
    print("=== Marinobacter S4 ===")
    for strain in MARINO_STRAINS:
        src = PAPER_DIR / MARINO_TEMPLATE.format(strain=strain)
        fc_map = dict(MARINO_FC_COLS)
        fc_map.update(MARINO_FC_OVERRIDES_BY_STRAIN.get(strain, {}))
        process(src, fc_map)
    print("=== S3 (Alteromonas) SKIPPED — DEH24 bridge deferred ===")


if __name__ == "__main__":
    main()
