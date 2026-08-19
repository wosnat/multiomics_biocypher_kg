#!/usr/bin/env python3
"""Extract Steglich 2010 Additional file 5 (MOESM5) whole-transcriptome
RNA half-life table into a flat CSV.

The legend's file numbering is offset by four (files 1-4 are PPTs), so
"Additional file 1"'s table actually ships as ``13059_2010_2347_MOESM5_ESM.XLS``
(sheet ``Tabelle1``, 2,043 rows). Columns carry two gene annotations
(``old annotation`` = canonical PMM*, ``new annotation`` = PMED4_*), a
log2 expression level at t0, a boolean expression call, a decay-cluster
assignment (1-12; blank for unclustered genes), and the fitted half-life /
decay-rate values with their standard errors.

Output: ``steglich2010_halflife_modified.csv`` alongside the XLS, with
clean snake_case headers, integer-string cluster labels (empty when
unassigned), and no other transformation. The two SE columns are kept in
the CSV for provenance but are wired into the paperconfig only via
``field_description`` (they are precision, not measurement).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Steglich 2010")
SRC_XLS = PAPER_DIR / "13059_2010_2347_MOESM5_ESM.XLS"
OUT_CSV = PAPER_DIR / "steglich2010_halflife_modified.csv"

COLUMN_MAP = {
    "old annotation": "old_annotation",
    "new annotation": "new_annotation",
    "gene name and function": "gene_function",
    "expression at 0 min in log2": "expression_t0_log2",
    "expression": "expressed",
    "cluster": "decay_cluster",
    "half-life time [min]": "half_life_min",
    "standard error of half-life time [min]": "half_life_se_min",
    "decay rate [min]": "decay_rate_per_min",
    "standard error of decay rate [min]": "decay_rate_se_per_min",
    "lower bound of error interval of decay rate [min]": "decay_rate_lower",
    "upper bound of error interval of decay rate [min]": "decay_rate_upper",
}


def main() -> None:
    df = pd.read_excel(SRC_XLS, sheet_name="Tabelle1")
    df = df.rename(columns=COLUMN_MAP)
    assert list(df.columns) == list(COLUMN_MAP.values()), df.columns

    # Cluster: float 1.0-12.0 or NaN → integer string or empty.
    df["decay_cluster"] = df["decay_cluster"].map(
        lambda v: "" if pd.isna(v) else str(int(v))
    )

    df.to_csv(OUT_CSV, index=False)

    clustered = df[df["decay_cluster"] != ""]
    print(f"rows: {len(df)}")
    print(f"non-null half_life_min: {df['half_life_min'].notna().sum()}")
    print(f"clustered rows: {len(clustered)}")
    print("cluster sizes:")
    counts = clustered["decay_cluster"].value_counts()
    for c in sorted(counts.index, key=int):
        print(f"  cluster {c}: {counts[c]}")


if __name__ == "__main__":
    main()
