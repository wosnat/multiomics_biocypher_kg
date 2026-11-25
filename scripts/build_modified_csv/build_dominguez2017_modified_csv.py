#!/usr/bin/env python3
"""Merge Dominguez 2017 S3 Upregulated + Downregulated CSVs into a single _modified.csv.

Sign convention:
  - Rows where `Highest mean condition == "Azaserine"` (and Lowest == Control): up in azaserine vs control => positive log2FC
  - Rows where `Highest mean condition == "Control"` (and Lowest == Azaserine): down in azaserine vs control => negative log2FC

Output columns: all original columns plus `log2_fold_change` (signed).
Written to `table s3 Combined_modified.csv` alongside originals.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

GN_RE = re.compile(r"GN=(\S+)")
# UniProt accession pattern (https://www.uniprot.org/help/accession_numbers)
UNIPROT_ACC_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def entry_name(accession_cell) -> str:
    """'1::Q7VEC7_PROMA' -> 'Q7VEC7_PROMA' (the '<n>::' prefix is not part of the id)."""
    if not isinstance(accession_cell, str):
        return ""
    return accession_cell.split("::", 1)[-1].strip()


def accession_from_entry(accession_cell) -> str:
    """'1::Q7VEC7_PROMA' -> 'Q7VEC7'; '' for reviewed entry names ('1::HIS7_PROMA')."""
    stem = entry_name(accession_cell).rsplit("_", 1)[0]
    return stem if UNIPROT_ACC_RE.match(stem) else ""

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Domínguez 2017")
UP_CSV = PAPER_DIR / "table s3 Upregulated prots rel quant sys003172107st8.csv"
DOWN_CSV = PAPER_DIR / "table s3 Downregulated prots rel quant sys003172107st8.csv"
OUT_CSV = PAPER_DIR / "table s3 Combined_modified.csv"


def main() -> None:
    up = pd.read_csv(UP_CSV)
    down = pd.read_csv(DOWN_CSV)

    print(f"Upregulated rows: {len(up)}")
    print(f"Downregulated rows: {len(down)}")

    combined = pd.concat([up, down], ignore_index=True)

    def signed_log2(row: pd.Series) -> float:
        fc = row["Max fold change"]
        if pd.isna(fc) or fc <= 0:
            return float("nan")
        l2 = math.log2(float(fc))
        highest = str(row.get("Highest mean condition", "")).strip()
        # If Control is the higher condition => the comparison Azaserine vs Control is DOWN
        if highest.lower() == "control":
            return -l2
        return l2

    combined["log2_fold_change"] = combined.apply(signed_log2, axis=1)

    # Extract the GN=... token from the free-text Description column so the
    # paperconfig can declare a dedicated id column (gene_name Tier 3) instead
    # of mis-tagging the whole Description string as a locus_tag, which
    # tokenises sentence words into specific_lookup.
    combined["extracted_gn"] = (
        combined["Description"].fillna("").map(
            lambda s: (m := GN_RE.search(s)) and m.group(1) or ""
        )
    )

    # `Accession` is "<n>::<UniProt entry name>" ("1::Q7VEC7_PROMA"); the prefix makes
    # every lookup miss. For unreviewed entries the name is "<accession>_PROMA", so the
    # accession can be recovered. Without it, rows fall back to the shared GN= symbol,
    # and five protein pairs (sufB, hisB, metB, ompR, wecE) collapsed onto one gene each.
    combined["uniprot_entry_name"] = combined["Accession"].map(entry_name)
    combined["uniprot_accession"] = combined["Accession"].map(accession_from_entry)

    combined.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(combined)} rows to {OUT_CSV}")
    print(f"log2FC range: {combined['log2_fold_change'].min():.3f} to {combined['log2_fold_change'].max():.3f}")
    print(f"  positive (up in azaserine): {(combined['log2_fold_change'] > 0).sum()}")
    print(f"  negative (down in azaserine): {(combined['log2_fold_change'] < 0).sum()}")


if __name__ == "__main__":
    main()
