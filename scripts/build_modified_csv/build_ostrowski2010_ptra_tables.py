"""Build the Ostrowski 2010 ptrA-mutant DE table (CSV).

Ostrowski et al. 2010 (ISME J 4:908, doi 10.1038/ismej.2010.24) — Synechococcus
sp. WH8102 ptrA (SYNW1019) insertional mutant vs isogenic wild type, both under
P stress (10 uM PO4 Aquil, harvested ~125-150 h after inoculation, just after
the onset of phosphatase induction in the wild type). Whole-genome microarray,
SAM with FDR < 1%. The paper publishes only the genes changed > 2-fold:

  * main-text Table 1 — 68 genes down > 2-fold in the mutant
  * main-text Table 2 — 64 genes up > 2-fold in the mutant
  * Supplementary Table 3 (`..._moesm21_esm.doc`) — 29 further down-regulated
    genes of unknown function (Table 1 footnote; 68 + 29 = the 97 genes the
    text counts)

Each table also carries a "wild type early P stress" log2FC column copied from
Tetu et al. 2009 (an independent experiment, integrated from Tetu 2009's own
supplement). It is NOT emitted here.

Tables 1/2 are read from the article PDF with `pdftotext -layout`; the PDF's
minus glyph comes out as \\x02 and is turned back into "-". Supplementary
Tables 1/2 (moesm19/moesm20) are primer lists — no expression data.

Output (idempotent), in the paper dir:
  ostrowski2010_ptra_vs_wt_gt2fold.csv with columns
    locus_tag, product, sam_score, log2_fold_change, source_table
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Synechococcus/papers_and_supp/Ostrowski 2010"
PDF = PAPER / "41396_2010_article_bfismej201024.pdf"
S3_DOC = PAPER / "41396_2010_bfismej201024_moesm21_esm.doc"
OUT = PAPER / "ostrowski2010_ptra_vs_wt_gt2fold.csv"

NUM = r"(\x02?\d+(?:\.\d+)?)"
ROW_RE = re.compile(
    r"(SYNW\d{4})\s+(\S.*?)\s{2,}" + NUM + r"\s+" + NUM
    + r"\s+(\x02?\d+(?:\.\d+)?|—[a-z]?)\s*$"
)


def _num(s: str) -> float:
    return float(s.replace("\x02", "-"))


def parse_pdf_tables() -> list[dict]:
    text = subprocess.run(
        ["pdftotext", "-layout", str(PDF), "-"],
        check=True, capture_output=True, text=True,
    ).stdout
    rows, sec = [], None
    for line in text.splitlines():
        if "Table 1 Synechococcus sp. WH8102 genes whose" in line:
            sec = "table_1_down"
        elif "Table 2 Synechococcus sp. WH8102 genes whose" in line:
            sec = "table_2_up"
        elif "Comparison of the relative expression" in line:  # Figure 2 caption
            sec = None
        if sec and (m := ROW_RE.search(line)):
            lt, product, sam, lfc, _wt_early = m.groups()
            rows.append({
                "locus_tag": lt,
                "product": product.strip(),
                "sam_score": _num(sam),
                "log2_fold_change": _num(lfc),
                "source_table": sec,
            })
    return rows


def parse_s3() -> list[dict]:
    text = subprocess.run(
        ["antiword", "-w", "0", str(S3_DOC)],
        check=True, capture_output=True, text=True,
    ).stdout
    rows = []
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 4 and re.fullmatch(r"SYNW\d{4}", cells[0]):
            rows.append({
                "locus_tag": cells[0],
                "product": cells[1],
                "sam_score": float(cells[2]),
                "log2_fold_change": float(cells[3]),
                "source_table": "supp_table_3_down",
            })
    return rows


def main() -> None:
    df = pd.DataFrame(parse_pdf_tables() + parse_s3())
    counts = df["source_table"].value_counts().to_dict()
    assert counts == {"table_1_down": 68, "table_2_up": 64, "supp_table_3_down": 29}, counts
    assert df["locus_tag"].is_unique
    down = df["source_table"].str.endswith("_down")
    assert (df.loc[down, "log2_fold_change"] <= -1).all()
    assert (df.loc[~down, "log2_fold_change"] >= 1).all()
    df.to_csv(OUT, index=False)
    print(f"{OUT.name}: {len(df)} rows {counts}")


if __name__ == "__main__":
    main()
