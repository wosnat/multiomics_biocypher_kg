#!/usr/bin/env python3
"""Koch 2019 Table S1 (Alteromonas macleodii 83-1, glucose vs laminarin/alginate/pectin mix) -> CSV.

Source: 41396_2018_252_moesm8_esm.xlsx, sheet "Table S1": 4,227 genes (Alt831_ locus tags, IMG gene
ids) with %NSAF (G1, G2, M1, M2; "nd" = not detected), DESeq2 log2 fold changes and adjusted P for
M1 vs M2, G1 vs G2, G1 vs M1, G2 vs M2, and RPKM. G = glucose, M = polysaccharide mix; phase 1 =
13 h, phase 2 = 21 h (M1 = laminarin-fuelled growth, M2 = alginate + pectin use without net growth).

Sign convention in the paper: "X vs Y" = log2(X / Y) — checked on the laminarin GH16
Alt831_03772 (M1 vs M2 = +7.2) and the pectate lyase Alt831_00890 (G2 vs M2 = -7.8, up in M2).
The KG orients each contrast as treatment vs control (later phase vs earlier, polysaccharide vs
glucose), so every column is negated:
    m2_vs_m1 = -(M1 vs M2), g2_vs_g1 = -(G1 vs G2), m1_vs_g1 = -(G1 vs M1), m2_vs_g2 = -(G2 vs M2).

Output: `koch2019_table_s1.csv` (one row per gene).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Koch 2019")
SRC = PAPER_DIR / "41396_2018_252_moesm8_esm.xlsx"
OUT = PAPER_DIR / "koch2019_table_s1.csv"
COLS = ["locus_tag", "img_gene_id", "product", "nsaf_g1", "nsaf_g2", "nsaf_m1", "nsaf_m2",
        "lfc_m1_vs_m2", "lfc_g1_vs_g2", "lfc_g1_vs_m1", "lfc_g2_vs_m2",
        "padj_m1_vs_m2", "padj_g1_vs_g2", "padj_g1_vs_m1", "padj_g2_vs_m2"]


def main() -> None:
    raw = pd.read_excel(SRC, sheet_name="Table S1", header=None, skiprows=2).iloc[:, :15]
    raw.columns = COLS
    raw = raw[raw["locus_tag"].astype(str).str.startswith("Alt831_")].copy()
    for c in COLS[3:]:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")  # "nd" -> NaN
    out = pd.DataFrame({
        "locus_tag": raw["locus_tag"].str.strip(),
        "img_gene_id": raw["img_gene_id"].astype("Int64").astype(str),
        "product": raw["product"],
        "log2fc_m2_vs_m1": -raw["lfc_m1_vs_m2"], "padj_m2_vs_m1": raw["padj_m1_vs_m2"],
        "log2fc_g2_vs_g1": -raw["lfc_g1_vs_g2"], "padj_g2_vs_g1": raw["padj_g1_vs_g2"],
        "log2fc_m1_vs_g1": -raw["lfc_g1_vs_m1"], "padj_m1_vs_g1": raw["padj_g1_vs_m1"],
        "log2fc_m2_vs_g2": -raw["lfc_g2_vs_m2"], "padj_m2_vs_g2": raw["padj_g2_vs_m2"],
        "nsaf_g1": raw["nsaf_g1"], "nsaf_g2": raw["nsaf_g2"],
        "nsaf_m1": raw["nsaf_m1"], "nsaf_m2": raw["nsaf_m2"],
    })
    out.to_csv(OUT, index=False)
    print(f"{OUT}: {len(out)} genes")


if __name__ == "__main__":
    main()
