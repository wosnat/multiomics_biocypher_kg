#!/usr/bin/env python3
"""Convert Oleza 2015 Table S3 (xlsx) to a clean CSV for KG integration.

Source: pmic8102-sup-0003-tables3.xlsx, sheet "Table S3"
  - row 0: long title
  - row 1: super-header (blank | WH7803 | WH8102 | "Prediction for Secretion")
  - row 2: column-name row
  - rows 3+: 72 R. pomeroyi DSS-3 polypeptides with per-condition NSAF triplicates

Output: `table_s3_dss3_exoproteome.csv` — one row per protein, flat columns,
plus a derived `detection_status` column per condition computed from the A/B/C
replicates: detected (3/3 non-zero), sporadic (1-2/3), not_detected (0/3 or all NaN).
The `Average % NSAF` (paper-printed) feeds a numeric DerivedMetric and the
detection_status feeds a categorical DerivedMetric — together they separate
intensity-when-detected from confidence.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Oleza 2015")
SRC = PAPER_DIR / "pmic8102-sup-0003-tables3.xlsx"
OUT = PAPER_DIR / "table_s3_dss3_exoproteome.csv"


def main() -> None:
    df = pd.read_excel(SRC, sheet_name="Table S3", header=[1, 2])

    flat: dict[str, pd.Series] = {}
    for (super_h, sub_h) in df.columns:
        if pd.isna(sub_h) or str(sub_h).startswith("Unnamed"):
            continue
        if str(super_h).startswith("Unnamed"):
            key = str(sub_h).strip().lower().replace(" ", "_").replace("%", "pct")
        else:
            sup = str(super_h).strip().lower()
            sub = str(sub_h).strip().lower()
            sub_norm = (sub.replace("average % nsaf", "avg_nsaf")
                            .replace("% nsaf ", "nsaf_")
                            .replace(" ", "_")
                            .replace("%", "pct"))
            key = f"{sup}_{sub_norm}"
        flat[key] = df[(super_h, sub_h)]

    out = pd.DataFrame(flat)

    rename = {
        "annotated_function_and_designation": "annotated_function",
        "additional_function_search": "additional_function",
    }
    out = out.rename(columns=rename)
    drop = [c for c in out.columns if c.endswith(".1")
            or c == "protein_sequence"
            or c.startswith("prediction for secretion")]
    out = out.drop(columns=drop, errors="ignore")
    out = out.dropna(subset=["ncbi_reference"])

    for cond in ("wh7803", "wh8102"):
        rep_cols = [f"{cond}_nsaf_a", f"{cond}_nsaf_b", f"{cond}_nsaf_c"]
        reps = out[rep_cols]
        all_nan = reps.isna().all(axis=1)
        n_non_zero = (reps.fillna(0) > 0).sum(axis=1)
        status = pd.Series("not_detected", index=out.index)
        status[n_non_zero.between(1, 2)] = "sporadic"
        status[n_non_zero == 3] = "detected"
        status[all_nan] = "not_detected"
        out[f"{cond}_detection_status"] = status

    out.to_csv(OUT, index=False)
    print(f"{SRC.name} -> {OUT.name}")
    print(f"  rows: {len(out)}")
    print(f"  columns: {list(out.columns)}")
    nz_w7803 = out["wh7803_avg_nsaf"].notna().sum()
    nz_w8102 = out["wh8102_avg_nsaf"].notna().sum()
    print(f"  proteins with WH7803 avg_nsaf: {nz_w7803}")
    print(f"  proteins with WH8102 avg_nsaf: {nz_w8102}")


if __name__ == "__main__":
    main()
