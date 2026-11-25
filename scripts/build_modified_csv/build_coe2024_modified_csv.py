#!/usr/bin/env python3
"""Coe 2024 Supplemental Table 2: drop RAST rows that repeat an NCBI-model gene.

The table quantifies NATL2A genes against two annotations, NCBI CDS models
(`Gene ID` = ``cds-WP_...``) and RAST models (``fig|1218.94.peg.N``). Most genes
have one row; 7 have both, with identical values where the two models agree and
different values where they don't (PMN2A_0530, _1491, _1754). The two rows share
one expression-edge id, so the graph silently kept whichever came last.

The KG's genes are NCBI models, so a RAST row is dropped when a ``cds-`` row
carries the same ``NCBI ID_3`` locus tag; the 164 RAST rows that are a gene's only
row are kept, and an exact repeat of a row is dropped. Writes ``Supplemental Table 2_modified.csv``; idempotent; never edits
the original.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/coe 2024"
SRC = PAPER / "Supplemental Table 2.csv"
OUT = PAPER / "Supplemental Table 2_modified.csv"
RAST_OUT = PAPER / "coe2024_rast_id_translation.csv"


def drop_rast_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    gene_id = df["Gene ID"].fillna("").astype(str)
    tag = df["NCBI ID_3"].fillna("").astype(str).str.strip()
    ncbi_tags = set(tag[gene_id.str.startswith("cds-") & (tag != "")])
    rast_dupe = gene_id.str.startswith("fig|") & tag.isin(ncbi_tags)
    # nirA (PMN2A_1298) is also listed twice as the same cds- row, values identical.
    return df[~rast_dupe].drop_duplicates()


def main() -> None:
    df = pd.read_csv(SRC, dtype=str, keep_default_na=False)
    out = drop_rast_duplicates(df)
    out.to_csv(OUT, index=False)
    rast = df[df["Gene ID"].str.startswith("fig|") & (df["NCBI ID_3"].str.strip() != "")]
    rast[["Gene ID", "NCBI ID_3"]].drop_duplicates().to_csv(RAST_OUT, index=False)
    print(f"wrote {RAST_OUT.relative_to(ROOT)}: {len(rast)} RAST ids")
    print(f"dropped {len(df) - len(out)} rows (RAST repeats of an NCBI-model gene, exact repeats); "
          f"wrote {OUT.relative_to(ROOT)}: {len(out)} rows")


if __name__ == "__main__":
    main()
