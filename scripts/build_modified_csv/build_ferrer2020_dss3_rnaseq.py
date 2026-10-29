"""Build the R. pomeroyi DSS-3 RNA-seq DE table for Ferrer 2020 (ISME J).

Source: ``41396_2020_811_moesm4_esm.xlsx`` sheet ``DSS3 DESeq`` (Table S5 /
SI Table S1): DESeq2 results for DSS-3 after 8 h in co-culture with
Thalassiosira pseudonana CCMP1335 vs an axenic DSS-3 glucose control.

The legend says "Log2 fold-differences that are negative are enriched in the
co-culture, differences that are positive are enriched in the glucose control",
i.e. the column is log2(control / coculture). The builder NEGATES it so
``log2fc_coculture_vs_glucose`` > 0 means higher in co-culture.

Skipped (organisms not in the KG): moesm2 (Polaribacter dokdonensis MED152) and
moesm3 (Stenotrophomonas sp. SKA14).

Idempotent; never edits the original xlsx.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/Ferrer 2020"
SRC = PAPER / "41396_2020_811_moesm4_esm.xlsx"
OUT = PAPER / "ferrer2020_dss3_coculture_vs_glucose.csv"


def clean(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip().str.replace("|", ",", regex=False).str.replace("'", "^", regex=False)


def main() -> None:
    df = pd.read_excel(SRC, sheet_name="DSS3 DESeq", skiprows=1)
    df = df[df["Original Locus Tag"].notna()]
    out = pd.DataFrame()
    out["locus_tag"] = clean(df["Original Locus Tag"])
    out["new_locus_tag"] = clean(df["New Locus Tag"])
    out["protein_id"] = clean(df["Protein Tag"])
    out["gene_symbol"] = clean(df["Gene ID"])
    out["product"] = clean(df["Description"])
    out["log2fc_coculture_vs_glucose"] = -pd.to_numeric(df["log2 Fold-difference"], errors="coerce")
    out["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    out = out[out["locus_tag"] != ""]
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(out)} rows")


if __name__ == "__main__":
    main()
