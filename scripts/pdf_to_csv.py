"""
Convert PDF supplementary tables to CSV files for:
- Read 2017 (3h, 12h, 24h post-starvation DE gene tables)
- Thompson 2016 (dark vs light, lysate vs spent medium DE gene tables)
"""

import re
import csv
import subprocess
import sys
from pathlib import Path


def extract_text(pdf_path):
    """Run pdftotext -layout and return text."""
    result = subprocess.run(
        ["pdftotext", "-layout", str(pdf_path), "-"],
        capture_output=True
    )
    return result.stdout.decode("utf-8", errors="replace")


def clean_text(s):
    """Remove PDF encoding artifacts: soft hyphen, tab, CR, non-breaking space.

    Replaces [\t\r\xa0]+ with single space so real multi-space column gaps
    (5+ spaces) are preserved for later splitting.
    """
    s = re.sub(r"[\t\r\xa0]+", " ", s)
    # Remove soft hyphen (U+00AD) and normalize unicode dash/hyphen variants.
    # The PDF encodes a single hyphen as: hyphen + soft-hyphen + figure-hyphen,
    # so after removing U+00AD we get double-hyphen; collapse to single.
    s = s.replace("\u00ad", "")
    s = s.replace("\u2013", "-").replace("\u2010", "-")
    s = re.sub(r"-{2,}", "-", s)
    return s.strip()


def clean_field(s):
    """Collapse all whitespace runs to a single space within a text field."""
    return re.sub(r" {2,}", " ", s).strip()


def parse_read_2017(pdf_path, time_label):
    """
    Parse Read 2017 DE gene tables.
    Columns: Name, Number, log2 Fold Change, Standard Deviation,
             Standard Error, p-value, Category, Definition
    """
    text = extract_text(pdf_path)
    rows = []

    # Regex to match a data row:
    # PMM\d+ followed by numeric fields then free text
    # The minus sign can be regular hyphen or soft-hyphen combinations
    neg = r"[-\u00ad\u2010\u2013]?"
    num_pattern = rf"{neg}\d+\.?\d*"
    # Match: locus_tag  count  log2fc  stddev  stderr  pval  rest
    data_re = re.compile(
        rf"^(PMM\d+)\s+(\d+)\s+({neg}\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(.*)"
    )

    for raw_line in text.split("\n"):
        line = clean_text(raw_line)
        m = data_re.match(line)
        if not m:
            continue

        name = m.group(1)
        number = m.group(2)
        log2fc = m.group(3).replace("\u00ad", "").replace("\u2010", "-")
        stddev = m.group(4)
        stderr = m.group(5)
        pval = m.group(6)
        rest = m.group(7).strip()

        # Category and Definition are separated by a large whitespace gap (5+ spaces).
        # Within each field, the PDF encodes word gaps as tab+CR+space+NBSP which
        # clean_text() converts to 3 spaces; we collapse those after splitting.
        parts = re.split(r" {5,}", rest, maxsplit=1)
        category = clean_field(parts[0]) if parts else ""
        definition = clean_field(parts[1]) if len(parts) > 1 else ""

        rows.append({
            "locus_tag": name,
            "number": number,
            "log2_fold_change": log2fc,
            "standard_deviation": stddev,
            "standard_error": stderr,
            "p_value": pval,
            "category": category,
            "definition": definition,
            "time_point": time_label,
        })

    return rows


def parse_thompson_2016(pdf_path, comparison_type):
    """
    Parse Thompson 2016 DE gene tables.

    Dark vs Light (s008): sections *** Uninfected/Infected: Xh post-inoculation ... ***
    Columns: locus_tag, dark_count, light_count, log2_fold_change, direction,
             gene_function, pathway, infection_status, time_point

    Lysate vs Spent Med (s007): sections *** Light/Dark: Xh post-inoculation ... ***
    Columns: locus_tag, infected_count, uninfected_count, log2_fold_change, direction,
             gene_function, pathway, light_condition, time_point
    """
    text = extract_text(pdf_path)
    rows = []

    current_section = {}
    section_re = re.compile(r"\*{3}\s*(.*?)\s*\*{3}")

    # Pattern for data rows:
    # locus_tag  count1  count2  fold_change  direction  gene_function  pathway
    # locus_tag may be PMM\d+, RNA_\d+, or N/A
    # direction is " (up arrow) or # (down arrow) in the PDF
    neg = r"[-\u2013\u2010\u00ad]?"
    data_re = re.compile(
        rf"^(PMM\d+|RNA_\d+|N/A)\s+(\d+)\s+(\d+)\s+({neg}\d+\.?\d*)\s+([\"#\u201c\u201d\u2019])\s+(.*)"
    )

    for raw_line in text.split("\n"):
        line = clean_text(raw_line)

        # Check for section header
        sec_m = section_re.search(line)
        if sec_m:
            section_text = sec_m.group(1).strip()
            current_section = parse_thompson_section(section_text, comparison_type)
            continue

        # Try to match a data row
        m = data_re.match(line)
        if not m:
            continue

        locus_tag = m.group(1)
        count1 = m.group(2)
        count2 = m.group(3)
        log2fc_raw = m.group(4)
        log2fc = log2fc_raw.replace("\u2013", "-").replace("\u2010", "-").replace("\u00ad", "")
        direction_char = m.group(5)
        direction = "up" if direction_char in ('"', "\u201c", "\u201d") else "down"
        rest = m.group(6).strip()

        # Split gene_function and pathway by 2+ spaces
        parts = re.split(r"  +", rest, maxsplit=1)
        gene_function = parts[0].strip() if parts else rest
        pathway = parts[1].strip() if len(parts) > 1 else ""

        row = {
            "locus_tag": locus_tag,
            "log2_fold_change": log2fc,
            "direction": direction,
            "gene_function": gene_function,
            "pathway": pathway,
        }

        if comparison_type == "dark_vs_light":
            row["dark_count_rpkm"] = count1
            row["light_count_rpkm"] = count2
            row["infection_status"] = current_section.get("infection_status", "")
            row["time_point_h"] = current_section.get("time_point_h", "")
            row["inoculation_type"] = current_section.get("inoculation_type", "")
        else:  # lysate_vs_spent_med
            row["infected_count_rpkm"] = count1
            row["uninfected_count_rpkm"] = count2
            row["light_condition"] = current_section.get("light_condition", "")
            row["time_point_h"] = current_section.get("time_point_h", "")

        rows.append(row)

    return rows


def parse_thompson_section(section_text, comparison_type):
    """Parse section header text into condition metadata."""
    result = {}

    if comparison_type == "dark_vs_light":
        # e.g. "Uninfected: 0.5 h post-inoculation (spent medium)"
        # or   "Infected: 4.5 h post-inoculation (phage)"
        m = re.match(
            r"(Uninfected|Infected):\s*([\d.]+)\s*h\s*post-inoculation\s*\(([^)]+)\)",
            section_text, re.IGNORECASE
        )
        if m:
            result["infection_status"] = m.group(1).lower()
            result["time_point_h"] = m.group(2)
            result["inoculation_type"] = m.group(3).strip()
    else:  # lysate_vs_spent_med
        # e.g. "Light: 0.5 h post-inoculation (phage or spent medium)"
        # or   "Dark: 2.5 h post-inoculation (phage or spent medium)"
        m = re.match(
            r"(Light|Dark):\s*([\d.]+)\s*h\s*post-inoculation",
            section_text, re.IGNORECASE
        )
        if m:
            result["light_condition"] = m.group(1).lower()
            result["time_point_h"] = m.group(2)

    return result


def write_csv(rows, out_path, fieldnames):
    """Write rows to CSV."""
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"  Wrote {len(rows)} rows -> {out_path}")


def main():
    base = Path("/fast_data/Osnat/multiomics_biocypher_kg/data/Prochlorococcus/papers_and_supp")

    # ── Read 2017 ──────────────────────────────────────────────────────────────
    read_dir = base / "Read 2017"
    read_fields = [
        "locus_tag", "number", "log2_fold_change", "standard_deviation",
        "standard_error", "p_value", "category", "definition", "time_point"
    ]

    read_pdfs = [
        (
            read_dir / "DE genes top 50 percent 3 hours post starvation 41396_2017_bfismej201788_moesm47_esm.pdf",
            "3h",
            read_dir / "DE_genes_3h_post_starvation.csv",
        ),
        (
            read_dir / "41396_2017_bfismej201788_moesm48_esm.pdf",
            "12h",
            read_dir / "DE_genes_12h_post_starvation.csv",
        ),
        (
            read_dir / "41396_2017_bfismej201788_moesm49_esm.pdf",
            "24h",
            read_dir / "DE_genes_24h_post_starvation.csv",
        ),
    ]

    print("Processing Read 2017 PDFs...")
    for pdf_path, time_label, out_path in read_pdfs:
        print(f"  Parsing {pdf_path.name} ({time_label})...")
        rows = parse_read_2017(pdf_path, time_label)
        write_csv(rows, out_path, read_fields)

    # ── Thompson 2016 ──────────────────────────────────────────────────────────
    thompson_dir = base / "Thompson 2016"

    print("\nProcessing Thompson 2016 PDFs...")

    # Dark vs Light
    dvl_pdf = thompson_dir / "DE genes and counts dark vs light pone.0165375.s008.pdf"
    dvl_fields = [
        "locus_tag", "dark_count_rpkm", "light_count_rpkm", "log2_fold_change",
        "direction", "infection_status", "time_point_h", "inoculation_type",
        "gene_function", "pathway"
    ]
    print(f"  Parsing {dvl_pdf.name}...")
    rows = parse_thompson_2016(dvl_pdf, "dark_vs_light")
    write_csv(rows, thompson_dir / "DE_genes_dark_vs_light.csv", dvl_fields)

    # Lysate vs Spent Medium
    lsm_pdf = thompson_dir / "DE  genes counts lyzate vs spent med pone.0165375.s007.pdf"
    lsm_fields = [
        "locus_tag", "infected_count_rpkm", "uninfected_count_rpkm", "log2_fold_change",
        "direction", "light_condition", "time_point_h",
        "gene_function", "pathway"
    ]
    print(f"  Parsing {lsm_pdf.name}...")
    rows = parse_thompson_2016(lsm_pdf, "lysate_vs_spent_med")
    write_csv(rows, thompson_dir / "DE_genes_lysate_vs_spent_medium.csv", lsm_fields)

    print("\nDone.")


if __name__ == "__main__":
    main()
