#!/usr/bin/env python3
"""Collapse Voigt 2014 TSS tables (S4 MED4 / S5 MIT9313) to per-gene scalars.

Tables S4 (moesm59, 4,190 rows) and S5 (moesm60, 8,746 rows) list individual
transcriptional start sites, several rows per gene, keyed to native
``oldLocusTag`` (PMM* / PMT*) with the newer ``locusTag`` (PMED4_* / P9313_*)
alongside. TSS types: gTSS (primary/genic), aTSS (antisense), iTSS
(internal), nTSS (orphan), aArray (array-only antisense).

Reduction to per-gene routing scalars (deliberately NOT a representation of
the paper's TSS/UTR/operon/ncRNA contribution — that stays Tier 3):

- gene universe   = rows of type gTSS/aTSS/aArray/iTSS carrying a gene id
                    (nTSS orphans are dropped: they are not gene TSSs)
- has_primary_tss = any gTSS row for the gene
- antisense_tss_count = count of aTSS + aArray rows
- internal_tss_count  = count of iTSS rows
- minus10_element_score = max "-10 element score" across the gene's gTSS rows
- tss_distance_to_cds   = "TSS distance" of the gene's strongest gTSS
  (max normalized Solexa read value; fallback raw Solexa, then 454, then the
  first gTSS row) — i.e. the 5'UTR length of the dominant primary TSS

Outputs (committed):
    data/Prochlorococcus/papers_and_supp/Voigt 2014/voigt2014_med4_tss_metrics.csv
    data/Prochlorococcus/papers_and_supp/Voigt 2014/voigt2014_mit9313_tss_metrics.csv
plus a row-count report printed to stdout (redirect into the commit message
or the paperconfig comment to keep the numbers checkable).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Voigt 2014")

SOURCES = {
    "med4": PAPER_DIR / "41396_2014_bfismej201457_moesm59_esm.xls",      # Table S4
    "mit9313": PAPER_DIR / "41396_2014_bfismej201457_moesm60_esm.xls",   # Table S5
}

GENE_TSS_TYPES = {"gTSS", "aTSS", "aArray", "iTSS"}


def strongest_gtss(rows: pd.DataFrame) -> pd.Series:
    """Pick the dominant gTSS row: max normalized Solexa reads, with raw
    Solexa then 454 reads as fallbacks; first row if no read values at all."""
    for col in ("read value normalized Solexa", "read value Solexa", "read value 454"):
        vals = pd.to_numeric(rows[col], errors="coerce")
        if vals.notna().any():
            return rows.loc[vals.idxmax()]
    return rows.iloc[0]


def build(strain: str, src: Path) -> pd.DataFrame:
    df = pd.read_excel(src, sheet_name=0)
    df = df[df["type"].isin(GENE_TSS_TYPES)].copy()
    # Gene key: prefer oldLocusTag, but ONLY when it looks like a real native
    # tag (PMM*/PMT*). The column also carries free text ("new ORF determined
    # by 454", IG* intergenic labels) shared across unrelated loci — keying on
    # those would merge distinct genes into one chimeric record. Anything
    # non-native falls back to the row's locusTag (distinct per locus).
    old = df["oldLocusTag"].astype("string")
    native = old.str.match(r"^PM[MT]\d", na=False)
    df["gene_key"] = old.where(native).fillna(df["locusTag"])
    df = df[df["gene_key"].notna()]

    records = []
    for key, rows in df.groupby("gene_key"):
        first = rows.iloc[0]
        g = rows[rows["type"] == "gTSS"]
        rec = {
            "old_locus_tag": first["oldLocusTag"] if pd.notna(first["oldLocusTag"]) else "",
            "locus_tag_new": first["locusTag"] if pd.notna(first["locusTag"]) else "",
            "gene": first.get("gene") if pd.notna(first.get("gene")) else "",
            "product": first.get("product") if pd.notna(first.get("product")) else "",
            "has_primary_tss": "true" if len(g) else "false",
            "antisense_tss_count": int((rows["type"].isin(["aTSS", "aArray"])).sum()),
            "internal_tss_count": int((rows["type"] == "iTSS").sum()),
            "minus10_element_score": np.nan,
            "tss_distance_to_cds": np.nan,
        }
        if len(g):
            scores = pd.to_numeric(g["-10 element score"], errors="coerce")
            if scores.notna().any():
                rec["minus10_element_score"] = scores.max()
            best = strongest_gtss(g)
            dist = pd.to_numeric(pd.Series([best["TSS distance"]]), errors="coerce").iloc[0]
            if pd.notna(dist):
                rec["tss_distance_to_cds"] = dist
        records.append(rec)

    out = pd.DataFrame.from_records(records)
    out_path = PAPER_DIR / f"voigt2014_{strain}_tss_metrics.csv"
    out.to_csv(out_path, index=False)

    print(f"{strain}: {len(out)} genes (from {len(df)} gene-linked TSS rows)")
    print(f"  has_primary_tss true: {(out['has_primary_tss'] == 'true').sum()}")
    print(f"  antisense_tss_count > 0: {(out['antisense_tss_count'] > 0).sum()}")
    print(f"  internal_tss_count > 0: {(out['internal_tss_count'] > 0).sum()}")
    print(f"  minus10_element_score non-null: {out['minus10_element_score'].notna().sum()}")
    print(f"  tss_distance_to_cds non-null: {out['tss_distance_to_cds'].notna().sum()}")
    return out


def main() -> None:
    for strain, src in SOURCES.items():
        build(strain, src)


if __name__ == "__main__":
    main()
