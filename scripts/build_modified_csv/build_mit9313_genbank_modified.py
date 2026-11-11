#!/usr/bin/env python3
"""Write MIT9313_genbank_modified.tsv: the MIT9313 ID bridge with stale RS locus tags blanked.

In a few rows of MIT9313_genbank.tsv the `locus_tag` (AKG35_RS*) column names a different
gene from the row's `PMTid`: the RS tag belongs to a paralog elsewhere in the genome.
P9313_04931 / PMT_0419 is psbA1 at 468 kb, but the row's RS tag AKG35_RS08155 is psbA2
(PMT1532) at 1.63 Mb. The row's `P9313name` numbering is genome-ordered, and in every such
row it places the gene where the PMTid's gene is, not where the RS tag's gene is, so the
RS tag is the stale cell. UniProt agrees (Q7TTH6: psbA1 = PMT_0419).

A row's `locus_tag` is blanked when the gene it names does not carry the row's PMTid while
another MIT9313 gene does. The other cells are kept.

Run from the repo root after prepare_data step 2 (reads gene_annotations_merged.json).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

RES_DIR = Path("data/Prochlorococcus/papers_and_supp/MIT9313_resources")
SRC = RES_DIR / "MIT9313_genbank.tsv"
OUT = RES_DIR / "MIT9313_genbank_modified.tsv"
MERGED = Path("cache/data/Prochlorococcus/genomes/MIT9313/gene_annotations_merged.json")


def main() -> None:
    genes = json.loads(MERGED.read_text())
    ncbi_to_lt = {g["locus_tag_ncbi"]: lt for lt, g in genes.items() if g.get("locus_tag_ncbi")}
    old_to_lt = {}
    for lt, g in genes.items():
        for old in g.get("old_locus_tags") or []:
            old_to_lt[old] = lt

    df = pd.read_csv(SRC, sep="\t", dtype=str, keep_default_na=False)
    blanked = []
    for i, row in df.iterrows():
        rs_gene = ncbi_to_lt.get(row["locus_tag"])
        pmt_gene = old_to_lt.get(row["PMTid"])
        if rs_gene and pmt_gene and rs_gene != pmt_gene:
            blanked.append((row["P9313name"], row["PMTid"], row["locus_tag"], rs_gene, pmt_gene))
            df.at[i, "locus_tag"] = ""
    df.to_csv(OUT, sep="\t", index=False)

    print(f"{len(blanked)} stale locus_tag cells blanked -> {OUT}")
    for p9313, pmt, rs, rs_gene, pmt_gene in blanked:
        print(f"  {p9313} {pmt}: {rs} is {rs_gene}, PMTid is {pmt_gene}")


if __name__ == "__main__":
    main()
