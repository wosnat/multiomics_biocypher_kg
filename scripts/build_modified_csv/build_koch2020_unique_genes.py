#!/usr/bin/env python3
"""Koch 2020 Table S7 (strain-unique genes) -> one CSV per KG strain.

Source: 41598_2020_57526_MOESM8_ESM.xlsx, sheet "TableS7_uniques_NCBI": genes present in
one Alteromonas macleodii genome and absent from the other 11 of the paper's pangenome, with
the top NCBI BLAST hit outside A. macleodii and its taxonomy.

Kept strains: the KG strains — MIT1002 (MIT1002_0xxxx), 83-1 (Alt831_xxxxx; the paper's
"83-1/BGP" row, BGP6 is a clonal twin), AD45 (AMBAS45_), ATCC27126 (MASE_), EZ55 (altEZ55_) and
HOT1A3 (Ga0125132_). The last two are JGI IMG locus tags (genomes 2785510739 / 2687453488),
bridged to the NCBI genes by the paperconfig's `id_translation` entries (diamond protein match
over `img_annotation/*.faa` + `.gff`). D7, EC673, Te101 are not KG strains.

Output: `unique_genes/<strain>_unique_genes.csv` with `locus_tag` (paper's tag),
`strain_unique` = "yes", and the closest-relative columns for reference.

Table S3 (MIT1002 24 h / 48 h vs 12 h in coculture with Prochlorococcus NATL2A) is NOT
used: it is a filtered re-analysis of the Biller 2016 data already in the KG.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Koch 2020")
SRC = PAPER_DIR / "41598_2020_57526_MOESM8_ESM.xlsx"
OUT_DIR = PAPER_DIR / "unique_genes"
STRAINS = {
    "MIT1002": "MIT1002",
    "83-1/BGP": "83-1",
    "AD45": "AD45",
    "ATCC27126": "ATCC27126",
    "EZ55": "EZ55",
    "HOT1A3": "HOT1A3",
}


def main() -> None:
    df = pd.read_excel(SRC, sheet_name="TableS7_uniques_NCBI")
    OUT_DIR.mkdir(exist_ok=True)
    for paper_strain, strain in STRAINS.items():
        sub = df[df["Strain"] == paper_strain].copy()
        out = pd.DataFrame({
            "locus_tag": sub["Locus tag"].astype(str).str.strip(),
            "strain_unique": "yes",
            "closest_relative": sub["Closest relative"],
            "closest_relative_phylum": sub["Phylum"],
            "closest_relative_genus": sub["Genus"],
        }).drop_duplicates("locus_tag")
        path = OUT_DIR / f"{strain}_unique_genes.csv"
        out.to_csv(path, index=False)
        print(f"{path}: {len(out)} genes")


if __name__ == "__main__":
    main()
