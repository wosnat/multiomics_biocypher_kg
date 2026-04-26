#!/usr/bin/env python3
"""Extract Biller 2014 vesicle proteomics Tables S2 (MED4), S3 (MIT9313) and the
DNA-seq Table S4 (MED4 top-50 vesicle ORFs) from the supplementary PDF into
clean CSVs.

S2/S3 rows: gi|<digits>|ref|<NP_xxx>| <desc> <mw> <mascot%>  <localization>
S4 rows:    <protein_id> <start> <stop> <description> <avg_reads>
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import pdfplumber

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/biller 2014")
PDF = PAPER_DIR / "1243457.biller.sm.revision1.pdf"
MED4_GFF = Path("cache/data/Prochlorococcus/genomes/MED4/genomic.gff")
GFF_LT_RE = re.compile(r"locus_tag=([^;]+)")

LOCALIZATIONS = (
    "Cytoplasmic Membrane",
    "Outer Membrane",
    "Cytoplasmic",
    "Periplasmic",
    "Extracellular",
    "Unknown",
)
LOC_RE = re.compile(r"\s+(" + "|".join(re.escape(l) for l in LOCALIZATIONS) + r")\s*$")
ROW_RE = re.compile(
    r"^gi\|(?P<gi>\d+)\|ref\|(?P<refseq>NP_[\d\.]+)\|\s+"
    r"(?P<rest>.+?)\s+(?P<mw>[\d,]+\.\d{2})\s+(?P<prob>\d+\.\d+)%$"
)


def parse_table(page_text: str, title_marker: str) -> list[dict]:
    rows: list[dict] = []
    in_table = False
    for raw in page_text.splitlines():
        line = raw.rstrip()
        if title_marker in line:
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        # Skip page numbers (e.g., "25", "26")
        if line.strip().isdigit():
            continue
        loc_m = LOC_RE.search(line)
        if not loc_m:
            continue
        localization = loc_m.group(1)
        head = line[: loc_m.start()].rstrip()
        m = ROW_RE.match(head)
        if not m:
            continue
        rows.append({
            "gi_id": m.group("gi"),
            "refseq_protein_id": m.group("refseq"),
            "description": m.group("rest").strip(),
            "molecular_weight_da": m.group("mw").replace(",", ""),
            "mascot_probability_pct": m.group("prob"),
            "predicted_localization": localization,
            # Use "Y" (matches Biller 2018 convention) — pandas would silently
            # coerce a uniform "true"/"false" column into Python bool, then
            # write it back capitalised, which the strict adapter rejects.
            "vesicle_proteome_member": "Y",
        })
    return rows


S4_ROW_RE = re.compile(
    r"^(?P<pid>\d+)\s+(?P<start>\d+)\s+(?P<stop>\d+)\s+(?P<desc>.+?)\s+(?P<reads>\d+)$"
)


def parse_s4(page_text: str) -> list[dict]:
    rows: list[dict] = []
    in_table = False
    for raw in page_text.splitlines():
        line = raw.rstrip()
        if "Table S4." in line:
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        if line.strip().isdigit():
            continue
        m = S4_ROW_RE.match(line)
        if not m:
            continue
        rows.append({
            "genbank_protein_id": m.group("pid"),
            "start": m.group("start"),
            "stop": m.group("stop"),
            "description": m.group("desc").strip(),
            "avg_reads": m.group("reads"),
        })
    return rows


def main() -> None:
    with pdfplumber.open(PDF) as pdf:
        # Table S2 = page index 26 (1-indexed page 27); S3 = index 27; S4 = index 28
        s2_text = pdf.pages[26].extract_text()
        s3_text = pdf.pages[27].extract_text()
        s4_text = pdf.pages[28].extract_text()

    s2 = parse_table(s2_text, "Table S2.")
    s3 = parse_table(s3_text, "Table S3.")
    s4 = parse_s4(s4_text)

    fields = [
        "gi_id", "refseq_protein_id", "description", "molecular_weight_da",
        "mascot_probability_pct", "predicted_localization", "vesicle_proteome_member",
    ]
    for label, rows in (("table_s2_med4_vesicle_proteome", s2),
                        ("table_s3_mit9313_vesicle_proteome", s3)):
        out = PAPER_DIR / f"{label}.csv"
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        locs = sorted({r["predicted_localization"] for r in rows})
        print(f"{label}: {len(rows)} rows -> {out.name} (localizations: {locs})")

    # Resolve S4 gi numbers to MED4 locus_tags via coordinate match against the
    # current NCBI GFF (Biller 2014 used the original 2003 NCBI annotation;
    # coordinates are identical or shifted by 1-3 bp).
    by_start: dict[int, str] = {}
    by_stop: dict[int, str] = {}
    for line in MED4_GFF.read_text().splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        m = GFF_LT_RE.search(parts[8])
        if not m:
            continue
        s, e, lt = int(parts[3]), int(parts[4]), m.group(1)
        by_start.setdefault(s, lt)
        by_stop.setdefault(e, lt)

    matched = 0
    for r in s4:
        s, e = int(r["start"]), int(r["stop"])
        lt = ""
        for delta in (0, -1, 1, -2, 2, -3, 3):
            if s + delta in by_start:
                lt = by_start[s + delta]; break
        if not lt:
            for delta in (0, -1, 1, -2, 2, -3, 3):
                if e + delta in by_stop:
                    lt = by_stop[e + delta]; break
        r["locus_tag"] = lt
        if lt:
            matched += 1

    s4_fields = ["locus_tag", "genbank_protein_id", "start", "stop", "description", "avg_reads"]
    s4_out = PAPER_DIR / "table_s4_med4_vesicle_dna_top50.csv"
    with s4_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=s4_fields)
        w.writeheader()
        w.writerows(s4)
    print(f"table_s4_med4_vesicle_dna_top50: {len(s4)} rows -> {s4_out.name} (locus_tag matched: {matched}/{len(s4)})")


if __name__ == "__main__":
    main()
