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

The xlsx's "New Locus Tag" / "Protein Tag" columns are not always on the row of
their "Original Locus Tag": in the megaplasmid block they are shifted by one row
(SPOA0353, the aldehyde dehydrogenase, carries SPO_RS21625 instead of NCBI's
SPO_RS21620), and ~370 protein accessions are from an older RefSeq annotation.
Used as gene-unique ids, the shifted tags merged neighbouring genes in the DSS-3
gene-ID mapping (90 Tier-1 conflicts). The builder blanks every new tag or
protein accession that contradicts NCBI's own pairing of the row's original tag
(``old_locus_tag`` in the strain's gene_mapping.csv); the original tag, which
matches the row's product, is kept.

Idempotent; never edits the original xlsx.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/Ferrer 2020"
SRC = PAPER / "41396_2020_811_moesm4_esm.xlsx"
OUT = PAPER / "ferrer2020_dss3_coculture_vs_glucose.csv"
GENE_MAPPING = ROOT / "cache/data/Ruegeria/genomes/DSS-3/gene_mapping.csv"


def clean(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip().str.replace("|", ",", regex=False).str.replace("'", "^", regex=False)


def blank_contradicting_ids(out: pd.DataFrame) -> pd.DataFrame:
    """Blank new tags / protein ids that NCBI pairs with a different original tag."""
    gm = pd.read_csv(GENE_MAPPING)
    old2new, new2pid = {}, dict(zip(gm["locus_tag_ncbi"], gm["protein_id"]))
    for new, olds in zip(gm["locus_tag_ncbi"], gm["old_locus_tags"].fillna("")):
        for old in olds.split(","):
            if old.strip():
                old2new[old.strip()] = new
    expected = out["locus_tag"].map(old2new)
    bad_new = expected.notna() & (out["new_locus_tag"] != "") & (out["new_locus_tag"] != expected)
    expected_pid = expected.map(new2pid)
    bad_pid = expected_pid.notna() & (out["protein_id"] != "") & (out["protein_id"] != expected_pid)
    out.loc[bad_new, "new_locus_tag"] = ""
    out.loc[bad_pid, "protein_id"] = ""
    print(f"blanked {bad_new.sum()} new_locus_tag and {bad_pid.sum()} protein_id cells that contradict NCBI")
    return out


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
    out = blank_contradicting_ids(out)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(out)} rows")


if __name__ == "__main__":
    main()
