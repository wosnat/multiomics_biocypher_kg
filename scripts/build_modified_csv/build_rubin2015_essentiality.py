#!/usr/bin/env python3
"""Extract Rubin 2015 (PNAS) Dataset S3 — the Synechococcus elongatus PCC 7942
RB-TnSeq gene-essentiality table — into a flat CSV.

Source: ``pnas.1519220112.sd03.xlsx`` (sheet ``gene_essentiality``, 2,723 genes).
Each gene carries an ``essentiality`` call:

* ``essential``      — log2 likelihood ratio < -2 on the insertion-index HMM
                        (lethal when mutated), confirmed by outgrowth fitness;
* ``beneficial``     — insertions reduced fitness below -1.32 (P < 0.01, FDR < 0.1;
                        t test over four standard lab outgrowths, ~6 generations);
* ``non-essential``  — log2 likelihood ratio > 2, no fitness defect;
* ``ambiguous``      — between the two likelihood cut-offs, or essential calls
                        not supported by the outgrowth fitness;
* ``not_analyzed``   — reads could not be mapped uniquely (psbA / psbD / cpc
                        paralogs with near-identical sequence).

Output ``rubin2015_essentiality.csv`` in the paper dir: snake_case headers and the
category normalised to snake_case (``non-essential`` -> ``non_essential``); no
rows dropped. ID columns: ``locus_tag_7942`` (Synpcc7942_NNNN, the GCF_000012525
old locus tag; plasmid genes use the paper's ``pANS_N`` / ``pANL`` style),
``refseq_locus_tag`` (SYNPCC7942_RS*), ``jgi_id`` (IMG gene oid).

Skipped (reported, not wired): sd01/sd02 (raw insertion positions + barcodes),
sd04 (ncRNA essentiality — no gene nodes), sd05 (upstream regulatory-region
essentiality — not genes), sd06 (essential intergenic regions).

Idempotent: re-running rewrites the same CSV.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Rubin 2015")
SRC = PAPER_DIR / "pnas.1519220112.sd03.xlsx"
OUT = PAPER_DIR / "rubin2015_essentiality.csv"

COLUMN_MAP = {
    "7942_ID": "locus_tag_7942",
    "jgi_ID": "jgi_id",
    "RefSeq_ID": "refseq_locus_tag",
    "start": "start",
    "end": "end",
    "strand": "strand",
    "chromosome": "replicon",
    "name": "gene_name",
    "description": "description",
    "essentiality": "essentiality",
    "TIGR_functional_catigory": "tigr_functional_category",
    "cyano_conserved": "cyano_conserved",
    "greencut2_conserved": "greencut2_conserved",
    "function_unknown": "function_unknown",
}

CATEGORY_MAP = {
    "essential": "essential",
    "beneficial": "beneficial",
    "non-essential": "non_essential",
    "ambiguous": "ambiguous",
    "not_analyzed": "not_analyzed",
}


def main() -> None:
    df = pd.read_excel(SRC, sheet_name="gene_essentiality")
    df = df[list(COLUMN_MAP)].rename(columns=COLUMN_MAP)
    df = df[df["locus_tag_7942"].notna()].copy()
    df["essentiality"] = df["essentiality"].str.strip().map(CATEGORY_MAP)
    if df["essentiality"].isna().any():
        raise SystemExit("unmapped essentiality category in sd03")
    df["jgi_id"] = df["jgi_id"].map(lambda v: "" if pd.isna(v) else str(int(v)))
    for col in ("description", "gene_name", "tigr_functional_category"):
        df[col] = df[col].fillna("").astype(str).str.replace("|", ",").str.replace("'", "^")
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(df)} rows)")
    print(df["essentiality"].value_counts().to_string())


if __name__ == "__main__":
    main()
