#!/usr/bin/env python3
"""Export Johnson 2026b (GSE314951) media-2 DESeq2 sheets to flat CSVs.

Every column header in ``media-2.xlsx`` carries an EMBEDDED NEWLINE
(e.g. ``log2FC \\n2H Dark vs Light``), which would have to be quoted
verbatim in the paperconfig; instead each DESeq2 sheet is exported to a
flat CSV with clean snake_case headers. Each sheet holds TWO timepoints
as paired log2FC / p-adjusted columns (4 contrasts total across the two
sheets). The DIMA sheets (per-iModulon activity differences) are NOT
exported — module activity has no schema slot (decision D5; backlog).

Outputs (committed):
    data/.../Johnson 2026b/johnson2026b_dark_vs_light_2_4h.csv
        refseq_locus, log2fc_2h, padj_2h, log2fc_4h, padj_4h
    data/.../Johnson 2026b/johnson2026b_light_vs_dark_14_16h.csv
        refseq_locus, log2fc_14h, padj_14h, log2fc_16h, padj_16h
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Johnson 2026b")
SRC = PAPER_DIR / "media-2.xlsx"

SHEETS = {
    "Dark vs Light 2-4h (DESeq2)": (
        "johnson2026b_dark_vs_light_2_4h.csv",
        ["refseq_locus", "log2fc_2h", "padj_2h", "log2fc_4h", "padj_4h"],
    ),
    "Light vs Dark 14-16h (DESeq2)": (
        "johnson2026b_light_vs_dark_14_16h.csv",
        ["refseq_locus", "log2fc_14h", "padj_14h", "log2fc_16h", "padj_16h"],
    ),
}


def main() -> None:
    for sheet, (out_name, cols) in SHEETS.items():
        df = pd.read_excel(SRC, sheet_name=sheet)
        assert df.shape[1] == len(cols), (sheet, df.columns)
        df.columns = cols
        out = PAPER_DIR / out_name
        df.to_csv(out, index=False)
        print(f"{sheet}: {len(df)} rows -> {out.name}")


if __name__ == "__main__":
    main()
