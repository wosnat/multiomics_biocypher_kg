"""Build the R. pomeroyi DSS-3 proteomics DE table for Barber 2026 (Chem Sci).

Source: ``SupplementaryData_1_R_pomeroyi_proteomics.xlsx``, sheet ``Imputed``
(1,838 protein groups; MaxQuant LFQ + Perseus, missing values imputed from a
normal distribution, Student's t-tests with permutation FDR q-values).

Perseus reports ``Difference <A>_<B>`` = mean log2 LFQ(A) - mean log2 LFQ(B).
For the three growth-substrate contrasts the first group is glucose
(``glu_L_cys``, ``glu_R-DHPS``, ``glu_S_DHPS``), so the Difference is
glucose - substrate. The KG wants treatment vs control, so the builder
NEGATES it: ``log2fc_<substrate>_vs_glucose = -Difference``. Check: the DHPS TRAP
SBP HpsK (AAV93906.1) reads Difference glu_R-DHPS = -8.95, i.e. +8.95 log2 on
R-DHPS, which matches the paper (Fig. 2, induced on DHPS).
``Difference R-DHPS_S_DHPS`` is already R - S and is kept as is.

Significance in the sheet ("C: Student's T-test Significant <contrast>" = '+')
equals q < 0.05 for every row (verified), so q-value with threshold 0.05 is the
KG significance call.

Gene id: first accession of ``T: Majority protein IDs`` (GenBank AAV*/AHB*
protein accessions, versioned; REV__ decoys stripped). ``fasta_locus_tag``
holds an SPO locus tag when the FASTA header names one (hypothetical proteins).

Skipped (not R. pomeroyi): Supplementary Data 2 (Bilophila wadsworthia) and
Supplementary Data 3 (TRAP SSN cluster ids).

Idempotent; never edits the original xlsx.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/Barber 2026"
SRC = PAPER / "SupplementaryData_1_R_pomeroyi_proteomics.xlsx"
OUT = PAPER / "barber2026_rpom_proteomics_de.csv"

T = "N: Student's T-test "
CONTRASTS = {
    # out label: (perseus contrast, sign)
    "L_cys_vs_glucose": ("glu_L_cys", -1.0),
    "R_DHPS_vs_glucose": ("glu_R-DHPS", -1.0),
    "S_DHPS_vs_glucose": ("glu_S_DHPS", -1.0),
    "R_DHPS_vs_S_DHPS": ("R-DHPS_S_DHPS", 1.0),
}
SPO_RE = re.compile(r"\b(SPO[A-Z]?\d{4}[a-z]?)\b")


def first_id(cell: str) -> str:
    ids = [i.strip() for i in str(cell).split(";") if i.strip() and not i.startswith("REV__")]
    return ids[0] if ids else ""


def main() -> None:
    df = pd.read_excel(SRC, sheet_name="Imputed")
    out = pd.DataFrame()
    out["protein_id"] = df["T: Majority protein IDs"].map(first_id)
    out["protein_ids_all"] = df["T: Protein IDs"].astype(str).str.replace(";", ",")
    hdr = df["T: Fasta headers"].astype(str).str.split(";").str[0]
    out["product"] = (
        hdr.str.replace(r"^\S+\s+", "", regex=True)
        .str.replace(r"\s*\[Ruegeria pomeroyi DSS-3\]\s*$", "", regex=True)
        .str.replace("|", ",", regex=False)
        .str.replace("'", "^", regex=False)
    )
    out["fasta_locus_tag"] = hdr.map(lambda s: (m := SPO_RE.search(s)) and m.group(1) or "")
    for label, (c, sign) in CONTRASTS.items():
        out[f"log2fc_{label}"] = sign * df[f"{T}Difference {c}"].astype(float)
        out[f"qvalue_{label}"] = df[f"{T}q-value {c}"].astype(float)
        out[f"pvalue_{label}"] = np.power(10.0, -df[f"N: -Log Student's T-test p-value {c}"].astype(float))
    out = out[out["protein_id"] != ""]
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(out)} rows")


if __name__ == "__main__":
    main()
