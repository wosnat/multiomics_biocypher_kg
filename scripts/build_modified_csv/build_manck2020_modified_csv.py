#!/usr/bin/env python3
"""Flatten the Manck 2020 (mSystems 5:e00070-20) supplementary data sets to CSV.

Alteromonas macleodii ATCC 27126, iron- and carbon-limitation RNA-seq (DESeq2).

Inputs (never edited):
  msystems.00070-20-sd001.xlsx  Data Set S1, sheets "Carbon Limitation" /
                                "Iron Limitation": DE ORFs only (|log2FC| >= 1,
                                FDR < 0.05), title row + header on row 2.
  msystems.00070-20-sd002.xlsx  Data Set S2: the 66 putative TonB-dependent
                                transporters (TBDTs) with significant log2FC
                                under Fe / C limitation ("n.s." otherwise).

Outputs (in the paper dir):
  manck2020_sd001_carbon_limitation.csv
  manck2020_sd001_iron_limitation.csv
      snake_case headers; values unchanged.
  manck2020_tbdts.csv
      one row per S2 TBDT: orf_id (MASE_RS*), orf_id_old (MASE_*),
      is_tbdt = "yes", predicted_substrate from main-text Table 1
      (gene-neighborhood prediction for the 26 Fe/C-responsive TBDTs;
      "not_predicted" when Table 1 has a dash; "not_assessed" for the 40 TBDTs
      Table 1 does not list). The S2 fold changes are NOT carried: they
      duplicate the Data Set S1 DE edges.

Data Set S3 (1,804 TBDTs over 31 Alteromonas genomes + network clusters) is
not converted: it needs a multi-organism cluster concept (backlog).

Idempotent: re-running overwrites the outputs with identical content.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Manck 2020")
SD1 = PAPER_DIR / "msystems.00070-20-sd001.xlsx"
SD2 = PAPER_DIR / "msystems.00070-20-sd002.xlsx"

SD1_COLUMNS = {
    "OrfID": "orf_id",
    "OrfID_Old": "orf_id_old",
    "Product": "product",
    "Condition": "condition",
    "Log2FoldChange": "log2_fold_change",
    "Standard Error": "lfc_se",
    "Wald Statistic": "wald_stat",
    "Wald Test p-value": "p_value",
    "FDR adjusted p-value": "fdr",
}

SHEETS = {
    "Carbon Limitation": ("C", "manck2020_sd001_carbon_limitation.csv"),
    "Iron Limitation": ("Fe", "manck2020_sd001_iron_limitation.csv"),
}

# Main-text Table 1: predicted substrate of the 26 Fe- or C-responsive TBDTs
# (gene-neighborhood analysis). Footnote b: MASE_17485 / MASE_17500 share one
# PepSY permease, so which of the two is the siderophore TBDT is ambiguous.
TABLE1_SUBSTRATE = {
    "MASE_00400": "not_predicted",
    "MASE_03260": "exogenous_siderophore",
    "MASE_08200": "exogenous_siderophore",
    "MASE_15245": "exogenous_siderophore",
    "MASE_17080": "exogenous_siderophore",
    "MASE_17485": "exogenous_siderophore",
    "MASE_17500": "exogenous_siderophore",
    "MASE_17930": "endogenous_siderophore",
    "MASE_18610": "not_predicted",
    "MASE_18645": "not_predicted",
    "MASE_18655": "ferric_citrate",
    "MASE_00935": "galactolipids",
    "MASE_03530": "peptide",
    "MASE_04115": "not_predicted",
    "MASE_05185": "not_predicted",
    "MASE_05565": "not_predicted",
    "MASE_06600": "polyhydroxybutyrate",
    "MASE_08095": "not_predicted",
    "MASE_10030": "peptide",
    "MASE_12160": "peptide",
    "MASE_13325": "alpha_glucans",
    "MASE_13870": "n_acetyl_beta_d_glucosaminides",
    "MASE_15330": "not_predicted",
    "MASE_18035": "xyloglucans",
    "MASE_19120": "alpha_arabinosides",
    "MASE_19670": "beta_glucans",
}


def build_sd1() -> None:
    for sheet, (cond, out_name) in SHEETS.items():
        df = pd.read_excel(SD1, sheet_name=sheet, header=1)
        df = df[list(SD1_COLUMNS)].rename(columns=SD1_COLUMNS)
        df = df[df["orf_id"].notna()]
        assert (df["condition"] == cond).all(), sheet
        assert not df["orf_id"].duplicated().any(), sheet
        df.to_csv(PAPER_DIR / out_name, index=False)
        print(f"{out_name}: {len(df)} rows")


def build_sd2() -> None:
    df = pd.read_excel(SD2, header=1)
    df = df.rename(columns={"OrfID": "orf_id", "OrfID_Old": "orf_id_old"})
    df = df[df["orf_id"].notna()][["orf_id", "orf_id_old"]]
    assert len(df) == 66, len(df)
    missing = set(TABLE1_SUBSTRATE) - set(df["orf_id_old"])
    assert not missing, missing
    df["is_tbdt"] = "yes"
    df["predicted_substrate"] = df["orf_id_old"].map(TABLE1_SUBSTRATE).fillna("not_assessed")
    out = PAPER_DIR / "manck2020_tbdts.csv"
    df.to_csv(out, index=False)
    print(f"{out.name}: {len(df)} rows; "
          f"{(df.predicted_substrate != 'not_assessed').sum()} with a Table 1 entry")


def main() -> None:
    build_sd1()
    build_sd2()


if __name__ == "__main__":
    main()
