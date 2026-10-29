"""Build Tetu 2009 per-contrast SAM tables (CSV) from the supplementary xls.

Tetu et al. 2009 (ISME J 3:835, doi 10.1038/ismej.2009.31) — Synechococcus sp.
WH8102 whole-genome microarray. Supplementary Table 1 (`..._moesm7_esm.xls`)
has one sheet per contrast, each the SAM output (FDR < 5%) as a "Positive
genes (n)" block followed by a "Negative genes (n)" block:

  S1 (A) Early P stress    WT 5 uM P (alkaline phosphatase first detected) vs WT P-replete 87 uM P mid-log
  S1 (B) Late P stress     WT 5 uM P (growth largely ceased) vs WT P-replete mid-log
  S1 (C) SYNW0947 knockout phoB mutant vs isogenic WT, both P-replete mid-log
  S1 (D) SYNW0948 knockout phoR mutant vs isogenic WT, both P-replete mid-log

Columns: Row, Gene ID ("SYNW1018/or0044"), Gene Name (product), SAM Score(d),
Numerator(r) (= mean log2 ratio, i.e. the log2 fold change; matches the
paper's Table 1), Denominator(s+s0), q-value(%) (percent).

Rows whose Gene ID is blank and Gene Name is "PSEUDOGENE" (2 per sheet at
most) name no gene and are dropped.

Outputs (idempotent), in the paper dir:
  tetu2009_s1<a-d>_<contrast>.csv with columns
    locus_tag, gene_id_raw, orf_id, product, sam_score, log2_fold_change,
    q_value (fraction, = q-value(%) / 100), sam_block (positive|negative)
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Synechococcus/papers_and_supp/Tetu 2009"
SRC = PAPER / "41396_2009_bfismej200931_moesm7_esm.xls"

SHEETS = {
    "S1 (A) Early P stress": "tetu2009_s1a_early_p_stress.csv",
    "S1 (B) Late P stress": "tetu2009_s1b_late_p_stress.csv",
    "S1 (C) SYNW0947 knockout": "tetu2009_s1c_synw0947_phob_knockout.csv",
    "S1 (D) SYNW0948 knockout": "tetu2009_s1d_synw0948_phor_knockout.csv",
}
LT_RE = re.compile(r"(SYNW\d{4})", re.I)
OR_RE = re.compile(r"(or\d+)", re.I)
BLOCK_RE = re.compile(r"(Positive|Negative) genes \((\d+)\)", re.I)


def parse_sheet(raw: pd.DataFrame) -> pd.DataFrame:
    rows, block, expected, pseudo = [], None, {}, {}
    for _, r in raw.iterrows():
        first = r.iloc[0]
        if isinstance(first, str) and (m := BLOCK_RE.search(first)):
            block = m.group(1).lower()
            expected[block] = int(m.group(2))
            continue
        gid = r.iloc[2]
        if block is not None and pd.isna(gid) and r.iloc[3] == "PSEUDOGENE":
            # array probe for an unnamed pseudogene: no gene id, cannot resolve
            pseudo[block] = pseudo.get(block, 0) + 1
            continue
        if block is None or not isinstance(gid, str) or gid.strip() == "Gene ID":
            continue
        lt = LT_RE.search(gid)
        orf = OR_RE.search(gid)
        rows.append({
            "locus_tag": lt.group(1).upper() if lt else "",
            "gene_id_raw": gid.strip(),
            "orf_id": orf.group(1) if orf else "",
            "product": str(r.iloc[3]).strip() if pd.notna(r.iloc[3]) else "",
            "sam_score": float(r.iloc[4]),
            "log2_fold_change": float(r.iloc[5]),
            "q_value": round(float(r.iloc[7]) / 100.0, 8),
            "sam_block": block,
        })
    df = pd.DataFrame(rows)
    for blk, n in expected.items():
        got = int((df["sam_block"] == blk).sum())
        assert got + pseudo.get(blk, 0) == n, f"{blk}: expected {n} rows, parsed {got}"
    # sign check: positive block has log2FC > 0, negative block < 0
    assert (df.loc[df.sam_block == "positive", "log2_fold_change"] > 0).all()
    assert (df.loc[df.sam_block == "negative", "log2_fold_change"] < 0).all()
    df.attrs["pseudogenes_dropped"] = sum(pseudo.values())
    return df


def main() -> None:
    book = pd.read_excel(SRC, sheet_name=None, header=None)
    for sheet, out in SHEETS.items():
        df = parse_sheet(book[sheet])
        df.to_csv(PAPER / out, index=False)
        print(f"{out}: {len(df)} rows, {(df.locus_tag == '').sum()} without SYNW tag, "
              f"{df.attrs['pseudogenes_dropped']} unnamed PSEUDOGENE rows dropped")


if __name__ == "__main__":
    main()
