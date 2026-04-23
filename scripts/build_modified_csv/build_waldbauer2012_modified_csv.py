#!/usr/bin/env python3
"""Extract Waldbauer 2012 Table S2 (diel cycling metrics for 312 genes) from PDF into CSV.

Table S2 is a 9-page PDF with per-gene numeric metrics. Several rows have
descriptions that wrap onto the same baseline as numeric values — pdfplumber's
word grouping interleaves them. Workaround: character-level extraction,
filtering numeric columns to ``[0-9.-]`` chars in the column's x-range.

Output: ``Table_S2_modified.csv`` (9 columns, 312 data rows) alongside the PDF.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import pdfplumber

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Waldbauer  2012")
SRC_PDF = PAPER_DIR / "Table_S2.pdf"
OUT_CSV = PAPER_DIR / "Table_S2_modified.csv"

COL_NAMES = [
    "pmed4_id",
    "gene_name",
    "gene_description",
    "peak_time_transcript_h",
    "peak_time_protein_h",
    "protein_transcript_lag_h",
    "transcript_amplitude_log2",
    "protein_amplitude_log2",
    "transcript_protein_amplitude_ratio",
]
COL_BOUNDS = [60, 128, 170, 395, 460, 524, 580, 640, 690, 760]
NUMERIC_COLS = list(range(3, 9))
NUMERIC_CHARS = set("0123456789.")  # all Table S2 values are non-negative (lag is mod 24)
NUM_RE = re.compile(r"\d+(?:\.\d+)?")
Y_TOL_ANCHOR = 2.5
Y_TOL_DESC_CONT = 12.0


def assign_col(x_center: float) -> int | None:
    for i in range(len(COL_BOUNDS) - 1):
        if COL_BOUNDS[i] <= x_center < COL_BOUNDS[i + 1]:
            return i
    return None


def extract_numeric_cell(chars, anchor_y: float, x_lo: float, x_hi: float) -> str:
    # Description text wraps onto the same baseline as numeric values; stray
    # digits from words like "(A1)" / "(NCBI)" / "P-loop" interleave with the
    # real value. Group numeric chars by x-proximity and pick the longest
    # contiguous group that parses as a number (ties broken by rightmost).
    in_band = [
        c for c in chars
        if abs(c["top"] - anchor_y) <= Y_TOL_ANCHOR
        and x_lo <= (c["x0"] + c["x1"]) / 2 < x_hi
        and c["text"] in NUMERIC_CHARS
    ]
    in_band.sort(key=lambda c: c["x0"])
    groups: list[list[dict]] = []
    GAP = 3.0  # x-distance that separates numeric value from adjacent desc digit
    for c in in_band:
        if groups and c["x0"] - groups[-1][-1]["x1"] <= GAP:
            groups[-1].append(c)
        else:
            groups.append([c])
    best: str | None = None
    best_key: tuple[int, float] = (-1, -1.0)
    for g in groups:
        s = "".join(c["text"] for c in g)
        m = NUM_RE.fullmatch(s)
        if not m:
            continue
        key = (len(s), g[-1]["x1"])
        if key > best_key:
            best_key = key
            best = s
    return best or "".join(c["text"] for c in in_band)


def extract_text_cell(words, anchor_y: float, next_anchor_y: float, x_lo: float, x_hi: float) -> str:
    in_band = [
        w for w in words
        if anchor_y - Y_TOL_ANCHOR <= w["top"] < next_anchor_y - Y_TOL_ANCHOR
        and x_lo <= (w["x0"] + w["x1"]) / 2 < x_hi
    ]
    in_band.sort(key=lambda w: (w["top"], w["x0"]))
    return " ".join(w["text"] for w in in_band).strip()


def main() -> None:
    rows: list[list[str]] = []
    with pdfplumber.open(SRC_PDF) as pdf:
        for page in pdf.pages:
            chars = page.chars
            words = page.extract_words(use_text_flow=False)
            anchors = sorted(
                [w["top"] for w in words if w["text"].startswith("PMED4_")]
            )
            page_bottom = page.height
            for i, anchor_y in enumerate(anchors):
                next_y = anchors[i + 1] if i + 1 < len(anchors) else anchor_y + Y_TOL_DESC_CONT
                row: list[str] = [""] * len(COL_NAMES)
                for ci in range(len(COL_NAMES)):
                    x_lo, x_hi = COL_BOUNDS[ci], COL_BOUNDS[ci + 1]
                    if ci in NUMERIC_COLS:
                        row[ci] = extract_numeric_cell(chars, anchor_y, x_lo, x_hi)
                    else:
                        row[ci] = extract_text_cell(words, anchor_y, next_y, x_lo, x_hi)
                rows.append(row)

    # Sanity checks
    assert len(rows) == 312, f"Expected 312 rows, got {len(rows)}"
    for r in rows:
        assert r[0].startswith("PMED4_"), f"Bad locus_tag: {r[0]!r}"
        for ci in NUMERIC_COLS:
            v = r[ci]
            assert v, f"Empty numeric cell at {r[0]} col {COL_NAMES[ci]}"
            float(v)  # will raise on malformed

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COL_NAMES)
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
