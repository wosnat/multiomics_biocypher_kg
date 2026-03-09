"""Extract supplementary tables from Thompson 2011 PDFs.

Uses pdfplumber table extraction for text, then character-level analysis
for bold/red formatting detection on numeric value cells.
"""

import pdfplumber
import csv
import re
import os
from collections import defaultdict

DIR = os.path.dirname(os.path.abspath(__file__))


def clean_gene_name(gene_text):
    """Clean gene name: remove footnote markers."""
    if not gene_text:
        return ''
    gene_text = gene_text.replace('‡', '').replace('†', '').replace('°', '')
    gene_text = ' '.join(gene_text.split())
    return gene_text.strip()


def clean_description(desc_text):
    """Clean description text."""
    if not desc_text:
        return ''
    # Collapse whitespace from multiline cells
    desc_text = ' '.join(desc_text.split())
    # Fix hyphenated line breaks (e.g., "Gram- negative" -> "Gram-negative")
    desc_text = re.sub(r'(\w)- (\w)', r'\1-\2', desc_text)
    # Fix "PM- 23" -> "PM-23" (number after hyphen)
    desc_text = re.sub(r'- (\d)', r'-\1', desc_text)
    return desc_text.strip()


def detect_value_formatting(page, table_obj, row_idx, col_idx):
    """Detect if a value cell contains bold and/or red numeric text.

    Uses the table's cell bounding boxes and character-level font info.
    """
    cells = table_obj.cells
    if not cells:
        return False, False

    # Group cells by row (y position)
    row_groups = defaultdict(list)
    for cell_bbox in cells:
        x0, y0, x1, y1 = cell_bbox
        row_key = round(y0, 0)
        row_groups[row_key].append(cell_bbox)

    sorted_row_keys = sorted(row_groups.keys())
    if row_idx >= len(sorted_row_keys):
        return False, False

    row_key = sorted_row_keys[row_idx]
    row_cells = sorted(row_groups[row_key], key=lambda c: c[0])

    if col_idx >= len(row_cells):
        return False, False

    x0, y0, x1, y1 = row_cells[col_idx]

    # Find numeric characters strictly within the cell bounding box (no y-margin
    # to avoid picking up characters from adjacent rows)
    x_margin = 3
    cell_chars = [c for c in page.chars
                  if c['x0'] >= x0 - x_margin and c['x0'] <= x1 + x_margin
                  and c['top'] >= y0 and c['top'] <= y1
                  and c['text'].strip() in '-.0123456789']

    if not cell_chars:
        return False, False

    # If chars span multiple y-positions (multiline cell), only use the first line
    # to avoid mixing formatting from the row below
    y_positions = sorted(set(round(c['top'], 0) for c in cell_chars))
    if len(y_positions) > 1:
        first_y = y_positions[0]
        cell_chars = [c for c in cell_chars if round(c['top'], 0) == first_y]

    if not cell_chars:
        return False, False

    # Use majority voting to avoid stray characters from adjacent columns
    bold_count = 0
    red_count = 0
    total = len(cell_chars)
    for c in cell_chars:
        fontname = c.get('fontname', '')
        if 'Bold' in fontname or 'bold' in fontname:
            bold_count += 1
        color = c.get('non_stroking_color', None)
        if color and isinstance(color, (list, tuple)) and len(color) >= 3:
            r, g, b = color[0], color[1], color[2]
            if r > 0.5 and g < 0.3 and b < 0.3:
                red_count += 1

    has_bold = bold_count > total / 2
    has_red = red_count > total / 2

    return has_bold, has_red


def extract_table(pdf_path):
    """Extract table rows with significance markers from a PDF."""
    all_rows = []

    with pdfplumber.open(pdf_path) as pdf:
        for page_num, page in enumerate(pdf.pages):
            settings = {
                "vertical_strategy": "lines",
                "horizontal_strategy": "lines",
                "snap_tolerance": 5,
            }
            found_tables = page.find_tables(settings)

            for table_obj in found_tables:
                table_data = table_obj.extract()
                if not table_data:
                    continue

                for row_idx, row in enumerate(table_data):
                    if not row or len(row) < 8:
                        continue

                    first = (row[0] or '').strip()

                    # Skip header/title/footnote rows
                    if 'Supplementary' in first or first.startswith('Gene'):
                        continue
                    if first.startswith('-') or first.startswith('*') or first.startswith('°'):
                        continue

                    gene = clean_gene_name(row[0])
                    description = clean_description(row[1])
                    cluster = (row[2] or '').strip()

                    # Process 5 value columns (indices 3-7)
                    values = []
                    for col_i in range(3, 8):
                        val_text = (row[col_i] or '').strip()
                        val_text = re.sub(r'\s+', '', val_text)

                        if val_text and re.match(r'^-?\d', val_text):
                            is_bold, is_red = detect_value_formatting(
                                page, table_obj, row_idx, col_i)
                            if is_red and is_bold:
                                val_text += '**'
                            elif is_bold:
                                val_text += '*'

                        values.append(val_text)

                    all_rows.append({
                        'gene': gene,
                        'description': description,
                        'cluster': cluster,
                        'values': values,
                    })

    return all_rows


def write_csv(rows, output_path, time_columns):
    """Write rows to CSV."""
    header = ['Gene', 'Description', 'Cluster'] + time_columns
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            csv_row = [row['gene'], row['description'], row['cluster']] + row['values']
            writer.writerow(csv_row)
    print(f"Wrote {len(rows)} data rows to {output_path}")


def compare_with_old(new_rows, old_csv_path, table_name):
    """Compare new extraction with old CSV and report differences."""
    if not os.path.exists(old_csv_path):
        return

    with open(old_csv_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        old_rows = list(reader)

    print(f"\n=== {table_name} Comparison ===")
    print(f"Old: {len(old_rows)} rows, New: {len(new_rows)} rows")

    # Build gene-indexed maps (strip significance markers from gene names for matching)
    def normalize_gene(g):
        g = g.strip()
        # Remove trailing ** or * from gene name (incorrectly added by old conversion)
        g = re.sub(r'\*+$', '', g)
        return g

    # Compare row by row (positional since gene names may differ)
    max_rows = max(len(old_rows), len(new_rows))
    diffs = 0
    for i in range(max_rows):
        if i >= len(new_rows):
            print(f"\n  Row {i+1}: MISSING in new (old: {old_rows[i][0]})")
            diffs += 1
            continue
        if i >= len(old_rows):
            print(f"\n  Row {i+1}: EXTRA in new ({new_rows[i]['gene']})")
            diffs += 1
            continue

        old_r = old_rows[i]
        new_r = new_rows[i]

        changes = []
        old_gene = old_r[0] if old_r else ''
        new_gene = new_r['gene']
        if normalize_gene(old_gene) != normalize_gene(new_gene):
            changes.append(f"gene: '{old_gene}' -> '{new_gene}'")

        old_desc = old_r[1] if len(old_r) > 1 else ''
        if old_desc != new_r['description']:
            changes.append(f"desc: '{old_desc[:50]}' -> '{new_r['description'][:50]}'")

        old_cluster = old_r[2] if len(old_r) > 2 else ''
        if old_cluster != new_r['cluster']:
            changes.append(f"cluster: '{old_cluster}' -> '{new_r['cluster']}'")

        old_vals = old_r[3:8] if len(old_r) >= 8 else old_r[3:]
        new_vals = new_r['values']
        for j, (ov, nv) in enumerate(zip(old_vals, new_vals)):
            if ov != nv:
                changes.append(f"val[{j}]: '{ov}' -> '{nv}'")

        if changes:
            diffs += 1
            gene_label = new_gene or old_gene
            print(f"\n  Row {i+1} ({gene_label}):")
            for c in changes:
                print(f"    {c}")

    print(f"\n  Total: {diffs} rows with differences out of {max_rows}")


def main():
    pdf1 = os.path.join(DIR, '41396_2011_bfismej201149_moesm54_esm.pdf')
    pdf2 = os.path.join(DIR, '41396_2011_bfismej201149_moesm55_esm.pdf')

    print("Processing Table 1 (MED4)...")
    rows1 = extract_table(pdf1)
    out1 = os.path.join(DIR, 'supp table 1_new.csv')
    write_csv(rows1, out1, ['0', '12', '24', '48', 'R*'])
    compare_with_old(rows1, os.path.join(DIR, 'supp table 1.csv'), 'Table 1 (MED4)')

    print("\n\nProcessing Table 2 (MIT9313)...")
    rows2 = extract_table(pdf2)
    out2 = os.path.join(DIR, 'supp table 2_new.csv')
    write_csv(rows2, out2, ['0', '16', '28', '53', 'R*'])
    compare_with_old(rows2, os.path.join(DIR, 'supp table 2.csv'), 'Table 2 (MIT9313)')


if __name__ == '__main__':
    main()
