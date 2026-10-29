#!/usr/bin/env python3
"""Flatten Christie-Oleza et al. 2012 (ISME J) Table S2 into one wide CSV.

Source: ``41396_2012_bfismej201186_moesm53_esm.xls`` — ten sheets ``Panel A`` ..
``Panel J``, each comparing whole-cell proteomes of Ruegeria pomeroyi DSS-3 after
20 h in natural seawater (spectral counts, PatternLab AC-test, duplicates):

  A-D  community present (+) vs absent (-): cape, beach, ferry basin, marina
  E-G  beach / ferry basin / marina (+) vs cape (-), community absent (filtered SW)
  H-J  beach / ferry basin / marina (+) vs cape (-), community present

Quirks handled here:
* The ``Panel A`` sheet holds panel A followed by copies of panels B-J; only the
  block up to the first footnote row is read from every sheet.
* ``Fold change`` is a signed LINEAR fold (footnote b: positive = up, negative =
  down; no value lies in (-1, 1)). Converted to log2: sign(fc) * log2(|fc|).
* ``p-value`` is the AC-test p (not FDR-adjusted). The paper's significance class
  (blue) is p < 0.01 and |fold| > 2, which the cell colours confirm exactly.
* ``Signal +`` / ``Signal -`` look swapped relative to the fold sign (e.g. panel J
  SPOA0370 fold +35 with Signal+ = 1, Signal- = 29); they are kept verbatim as
  ``signal_col7`` / ``signal_col8`` for provenance and are not used.

Output: ``oleza2012_table_s2_wide.csv`` (one row per protein, union over panels;
columns ``log2fc_<X>``, ``fold_change_<X>``, ``p_value_<X>`` per panel).
Idempotent.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Oleza 2012")
SRC = PAPER_DIR / "41396_2012_bfismej201186_moesm53_esm.xls"
OUT = PAPER_DIR / "oleza2012_table_s2_wide.csv"
PANELS = "ABCDEFGHIJ"


def read_panel(letter: str) -> pd.DataFrame:
    raw = pd.read_excel(SRC, sheet_name=f"Panel {letter}", header=None)
    rows = []
    started = False
    for _, r in raw.iterrows():
        ref = r.iloc[1]
        if isinstance(ref, str) and ref.startswith("YP_"):
            started = True
            rows.append(r.iloc[1:8].tolist())
        elif started:
            break  # end of the first block
    df = pd.DataFrame(rows, columns=["ncbi_reference", "locus_tag", "product",
                                     "fold_change", "p_value", "signal_col7", "signal_col8"])
    df["fold_change"] = pd.to_numeric(df["fold_change"])
    if ((df["fold_change"] > -1) & (df["fold_change"] < 1)).any():
        raise SystemExit(f"panel {letter}: fold change inside (-1,1) — not a signed linear fold")
    df["log2fc"] = df["fold_change"].map(lambda v: math.copysign(math.log2(abs(v)), v))
    if df["locus_tag"].duplicated().any():
        raise SystemExit(f"panel {letter}: duplicated locus tag")
    return df


def main() -> None:
    meta = {}
    wide = None
    for letter in PANELS:
        df = read_panel(letter)
        sig = ((df["p_value"] < 0.01) & (df["fold_change"].abs() > 2)).sum()
        print(f"panel {letter}: {len(df)} proteins, {sig} significant")
        for _, r in df.iterrows():
            meta.setdefault(r["locus_tag"], (r["ncbi_reference"], r["product"]))
        part = df.set_index("locus_tag")[["log2fc", "fold_change", "p_value"]]
        part.columns = [f"{c}_{letter}" for c in part.columns]
        wide = part if wide is None else wide.join(part, how="outer")
    wide.insert(0, "ncbi_reference", [meta[t][0] for t in wide.index])
    wide.insert(1, "product", [str(meta[t][1]).replace("|", ",").replace("'", "^") for t in wide.index])
    wide.index.name = "locus_tag"
    wide.sort_index().to_csv(OUT)
    print(f"wrote {OUT} ({len(wide)} proteins)")


if __name__ == "__main__":
    main()
