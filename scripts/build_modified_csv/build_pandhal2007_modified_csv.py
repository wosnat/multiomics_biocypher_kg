#!/usr/bin/env python3
"""Aggregate Pandhal 2007 iTRAQ workbook into per-comparison `_modified.csv` files.

Source workbook (`pr060460csi20061210_124449.xls`) holds two iTRAQ runs:
  - iTRAQ1 (148 proteins): tags 114=HL, 115=LL, 116=ML, 117=HL(rep)
      ratios reported as 114:115, 116:115, 117:115
  - iTRAQ2  (92 proteins): tags 114=LL, 115=ML, 116=HL, 117=HL(rep)
      ratios reported as 115:114, 116:114, 117:114

Two pairwise comparisons are measured (LL = reference for both):
  - HL vs LL: iTRAQ1 ratios 114:115 + 117:115, iTRAQ2 ratios 116:114 + 117:114
  - ML vs LL: iTRAQ1 ratio  116:115,           iTRAQ2 ratio  115:114

Per protein x comparison we emit one row with:
  - log2_fold_change = mean of log2(ratio) across contributing iTRAQ ratios
  - max_p_value      = max of P Val across contributing iTRAQ ratios (most conservative)
  - n_measurements   = count of non-null contributing ratios

Outputs (in the paper directory):
  - pandhal2007_HL_vs_LL_modified.csv
  - pandhal2007_ML_vs_LL_modified.csv
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Pandhal 2007")
WORKBOOK = PAPER_DIR / "pr060460csi20061210_124449.xls"

ITRAQ1_SHEET = "iTRAQ1 protein list"
ITRAQ2_SHEET = "iTRAQ2 protein list"

# Per-comparison: list of (sheet, ratio_col, p_val_col)
COMPARISONS: dict[str, list[tuple[str, str, str]]] = {
    "HL_vs_LL": [
        (ITRAQ1_SHEET, "Ratio 114:115", "P Val 114:115"),
        (ITRAQ1_SHEET, "Ratio 117:115", "P Val 117:115"),
        (ITRAQ2_SHEET, "Ratio 116:114", "P Val 116:114"),
        (ITRAQ2_SHEET, "Ratio 117:114", "P Val 117:114"),
    ],
    "ML_vs_LL": [
        (ITRAQ1_SHEET, "Ratio 116:115", "P Val 116:115"),
        (ITRAQ2_SHEET, "Ratio 115:114", "P Val 115:114"),
    ],
}


def _log2(x) -> float:
    try:
        v = float(x)
    except Exception:
        return float("nan")
    if pd.isna(v) or v <= 0:
        return float("nan")
    return math.log2(v)


def _read_sheet(sheet: str) -> pd.DataFrame:
    """Read a protein-list sheet. Row 0 is the header; protein rows start at row 1."""
    df = pd.read_excel(WORKBOOK, sheet_name=sheet, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _build_long(comparison: str) -> pd.DataFrame:
    """Long-form table: one row per (Accession, contributing-ratio).

    Columns: Accession, Protein Name, COG group, Cellular location,
             log2_ratio, p_value, source (which sheet+col supplied this row).
    """
    sources = COMPARISONS[comparison]
    rows: list[dict] = []
    for sheet, ratio_col, p_col in sources:
        df = _read_sheet(sheet)
        for _, r in df.iterrows():
            accession = r.get("Accession")
            if pd.isna(accession) or not str(accession).strip():
                continue
            ratio = r.get(ratio_col)
            if pd.isna(ratio):
                continue
            l2 = _log2(ratio)
            if math.isnan(l2):
                continue
            # ProteinPilot's weighted-error-factor "P Val" is not a true
            # probability and can exceed 1.0 in rare cases (numerical artifact
            # of the estimator). Clamp at 1.0 so the value remains a valid
            # "probability of no effect" -- semantically equivalent to
            # "definitely not significant". Also enforces the KG validity
            # invariant adjusted_p_value in [0, 1].
            raw_p = r.get(p_col)
            if pd.isna(raw_p):
                p_value = float("nan")
            else:
                p_value = min(max(float(raw_p), 0.0), 1.0)
            rows.append(
                {
                    "Accession": str(accession).strip(),
                    "Protein Name": r.get("Protein Name", ""),
                    "COG group": r.get("COG group", ""),
                    "Cellular location": r.get("Cellular location", ""),
                    "log2_ratio": l2,
                    "p_value": p_value,
                    "source": f"{sheet} :: {ratio_col}",
                }
            )
    return pd.DataFrame(rows)


def _aggregate(long_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse long-form to one row per Accession.

    log2_fold_change = mean(log2_ratio across contributing ratios)
    adjusted_p_value = max(p_value across contributing ratios, ignoring NaN)
    n_measurements   = count of contributing ratios
    sources          = '; '-joined source descriptors (deduped, ordered)
    """
    if long_df.empty:
        return long_df
    grp = long_df.groupby("Accession", sort=True)
    agg = grp.agg(
        log2_fold_change=("log2_ratio", "mean"),
        adjusted_p_value=("p_value", "max"),
        n_measurements=("log2_ratio", "count"),
    ).reset_index()

    # First non-null protein/COG/location per Accession (these don't vary across runs)
    def _first(values):
        for v in values:
            if pd.notna(v) and str(v).strip():
                return v
        return ""

    descriptors = grp.agg(
        protein_name=("Protein Name", _first),
        cog_group=("COG group", _first),
        cellular_location=("Cellular location", _first),
        sources=("source", lambda s: " ; ".join(dict.fromkeys(s))),
    ).reset_index()

    return agg.merge(descriptors, on="Accession").rename(
        columns={
            "protein_name": "Protein Name",
            "cog_group": "COG group",
            "cellular_location": "Cellular location",
        }
    )


def main() -> None:
    if not WORKBOOK.exists():
        raise SystemExit(f"Source workbook not found: {WORKBOOK}")

    for comparison in COMPARISONS:
        long_df = _build_long(comparison)
        out_df = _aggregate(long_df)
        # Stable column order
        out_df = out_df[
            [
                "Accession",
                "Protein Name",
                "log2_fold_change",
                "adjusted_p_value",
                "n_measurements",
                "COG group",
                "Cellular location",
                "sources",
            ]
        ]
        out_path = PAPER_DIR / f"pandhal2007_{comparison}_modified.csv"
        out_df.to_csv(out_path, index=False)
        sig = out_df[
            (out_df["log2_fold_change"].abs() >= math.log2(1.25))
            & (out_df["adjusted_p_value"] < 0.05)
        ]
        print(f"{comparison}: {len(out_df)} proteins -> {out_path.name}")
        print(
            f"  log2FC range: {out_df['log2_fold_change'].min():.3f} .. "
            f"{out_df['log2_fold_change'].max():.3f}"
        )
        print(
            f"  candidates passing |log2FC|>=log2(1.25) & p<0.05: {len(sig)}"
        )


if __name__ == "__main__":
    main()
