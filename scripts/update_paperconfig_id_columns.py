#!/usr/bin/env python3
"""
One-time migration script: add id_columns and product_columns to all existing
paperconfig.yaml csv table entries that are missing them.

Usage:
    uv run python scripts/update_paperconfig_id_columns.py           # dry-run
    uv run python scripts/update_paperconfig_id_columns.py --apply   # write changes
"""

import argparse
import re
from pathlib import Path

# ---------------------------------------------------------------------------
# Mapping: papername (or dir-fragment) -> table_key -> fields to inject
# Each entry has:
#   id_columns:      list of {column, id_type}
#   product_columns: list of {column}   (optional)
# ---------------------------------------------------------------------------

# Helper to build id_col and product_col lists
def idc(*cols):
    return [{"column": col, "id_type": id_type} for col, id_type in cols]

def prodc(*cols):
    return [{"column": col} for col in cols]


# Keyed by the paperconfig directory name (case-sensitive, as on disk)
MAPPINGS = {
    "Aharonovich 2016": {
        # supp_table_1, supp_table_2, supp_table_3 all share same CSV structure
        "_all_csv": {
            "id_columns": idc(("Synonym", "old_locus_tag")),
            "product_columns": prodc("Product"),
        },
    },

    "Al-Hosani 2015": {
        "_all_csv": {
            "id_columns": idc(("locus_tag", "locus_tag")),
        },
    },

    "bagby 2015": {
        "_all_csv": {
            "id_columns": idc(
                ("Locus tag2", "old_locus_tag"),
                ("Alternative locus tag", "alternative_locus_tag"),
                ("Probeset", "probeset"),
            ),
            "product_columns": prodc("Gene description"),
        },
    },

    "barreto 2022": {
        "_all_csv": {
            "id_columns": idc(("symbol", "gene_name")),
            "product_columns": prodc("product"),
        },
    },

    "biller 2016": {
        "supp_table_2": {
            "id_columns": idc(
                ("Original_NCBI_ID", "old_locus_tag"),
                ("Alternate_NCBI_ID", "locus_tag_cyanorak"),
                ("Alternate_NCBI_ID_2", "locus_tag_ncbi"),
                ("Gene name", "gene_name"),
            ),
            "product_columns": prodc("NCBI description", "RAST description"),
        },
        "supp_table_3": {
            "id_columns": idc(
                ("locus_tag", "locus_tag"),
                ("Gene ID", "old_locus_tag"),
                ("Gene name", "gene_name"),
            ),
            "product_columns": prodc("RAST annotation"),
        },
    },

    "Biller 2018": {
        "supp_table_s3": {
            "id_columns": idc(
                ("NCBI ID_3", "old_locus_tag"),
                ("NCBI ID_2", "locus_tag_ncbi"),
                ("NCBI ID", "locus_tag_cyanorak"),
                ("Gene Name", "gene_name"),
            ),
            "product_columns": prodc("Genbank Annotation", "RAST annotation"),
        },
        "supp_table_s6b": {
            "id_columns": idc(("RAST_region_ID", "annotation_specific")),
            "product_columns": prodc("RAST_Annotation"),
        },
    },

    "coe 2024": {
        "supp_table_2": {
            "id_columns": idc(
                ("NCBI ID_3", "old_locus_tag"),
                ("NCBI ID_2", "locus_tag_ncbi"),
                ("NCBI ID", "locus_tag_cyanorak"),
                ("Gene ID", "annotation_specific"),
                ("Gene Name", "gene_name"),
            ),
            "product_columns": prodc("product"),
        },
        "supp_table_4": {
            "id_columns": idc(
                ("NCBI ID", "old_locus_tag"),
                ("Gene ID", "annotation_specific"),
                ("Gene Name", "gene_name"),
            ),
            "product_columns": prodc("product"),
        },
    },

    "Fang 2019": {
        "supp_table_2": {
            "id_columns": idc(
                ("Gene ID", "old_locus_tag"),
                ("Gene", "gene_name"),
            ),
            "product_columns": prodc("Gene annotation"),
        },
    },

    "he 2022": {
        "_all_csv": {
            "id_columns": idc(
                ("locus_tag", "locus_tag"),
                ("Gene Name", "gene_name"),
            ),
            "product_columns": prodc("Definition"),
        },
    },

    "Hennon 2017": {
        "supp_table_2": {
            "id_columns": idc(("locus_tag", "locus_tag")),
            "product_columns": prodc("annotation"),
        },
        "supp_table_3": {
            "id_columns": idc(("locus tag", "locus_tag")),
            "product_columns": prodc("annotation"),
        },
        "supp_table_4": {
            "id_columns": idc(("locus tag", "locus_tag")),
            "product_columns": prodc("annotation"),
        },
    },

    "Labban 2022": {
        "_all_csv": {
            "id_columns": idc(("Gene ID", "locus_tag_ncbi")),
            "product_columns": prodc("Annotation"),
        },
    },

    "Lin 2015": {
        "_all_csv": {
            "id_columns": idc(
                ("ID", "old_locus_tag"),
                ("Gene Name", "gene_name"),
            ),
            "product_columns": prodc("Description"),
        },
    },

    "lindell 2007": {
        "supp_table_3": {
            "id_columns": idc(("ORF ", "old_locus_tag")),
            "product_columns": prodc(" gene name, possible product and function"),
        },
    },

    "martiny 2006": {
        "_all_csv": {
            "id_columns": idc(("ORF", "old_locus_tag")),
            "product_columns": prodc("Description"),
        },
    },

    "steglich 2006": {
        "supp_table_s1": {
            "id_columns": idc(("gene ID", "old_locus_tag")),
            "product_columns": prodc("Description"),
        },
    },

    "tetu 2019": {
        "supp_dataset_3_mit9312": {
            "id_columns": idc(
                ("MIT9312 original locus tag", "old_locus_tag"),
                ("MIT9312 (2017 NCBI)", "locus_tag_ncbi"),
                ("protein_id", "protein_id_refseq"),
            ),
            "product_columns": prodc("product"),
        },
        "supp_dataset_4_natl2a": {
            "id_columns": idc(
                ("NATL2A original locus tag", "old_locus_tag"),
                ("NATL2A locus tag (2017 NCBI)", "locus_tag_ncbi"),
                ("protein_id", "protein_id_refseq"),
            ),
            "product_columns": prodc("product"),
        },
    },

    "thompson 2011": {
        "_all_csv": {
            "id_columns": idc(
                ("locus_tag", "locus_tag"),
                ("Gene", "old_locus_tag"),
            ),
            "product_columns": prodc("Description"),
        },
    },

    "tolonen 2006": {
        "_all_csv": {
            "id_columns": idc(("GENE", "old_locus_tag")),
            "product_columns": prodc("ANNOTATION"),
        },
    },

    "Read 2017": {
        "_all_csv": {
            "id_columns": idc(("locus_tag", "locus_tag")),
            "product_columns": prodc("definition"),
        },
    },

    "Thompson 2016": {
        "_all_csv": {
            "id_columns": idc(("locus_tag", "locus_tag")),
            "product_columns": prodc("gene_function"),
        },
    },
}


def _yaml_scalar(value: str) -> str:
    """Quote a scalar value if it has leading/trailing whitespace."""
    if value != value.strip():
        return f'"{value}"'
    return value


def format_id_columns_yaml(indent: int, id_cols: list) -> str:
    """Format id_columns block as YAML text with given base indent."""
    lines = [f"{' ' * indent}id_columns:"]
    for e in id_cols:
        lines.append(f"{' ' * indent}- column: {_yaml_scalar(e['column'])}")
        lines.append(f"{' ' * (indent + 2)}id_type: {e['id_type']}")
    return "\n".join(lines) + "\n"


def format_product_columns_yaml(indent: int, prod_cols: list) -> str:
    """Format product_columns block as YAML text with given base indent."""
    lines = [f"{' ' * indent}product_columns:"]
    for e in prod_cols:
        lines.append(f"{' ' * indent}- column: {_yaml_scalar(e['column'])}")
    return "\n".join(lines) + "\n"


def find_table_statistical_analyses(text: str, table_key: str):
    """
    Find the position just before `statistical_analyses:` for a given table.
    Returns (insert_pos, indent) or (None, None) if not found.
    """
    import re
    # Find the table key line (e.g. "    supp_table_1:")
    # Allow any leading whitespace, then the exact key, then optional space and colon
    key_pattern = re.compile(r"^( *)(" + re.escape(table_key) + r"):[ \t]*$", re.MULTILINE)
    m = key_pattern.search(text)
    if not m:
        return None, None

    table_indent = len(m.group(1))
    content_indent = table_indent + 2  # standard 2-space yaml indent

    # From that point forward, find `statistical_analyses:` at content_indent
    start = m.end()
    sa_pattern = re.compile(
        r"^( {" + str(content_indent) + r"})statistical_analyses:[ \t]*$",
        re.MULTILINE,
    )
    sa_m = sa_pattern.search(text, start)
    if not sa_m:
        return None, None

    return sa_m.start(), content_indent


def already_has_field(text: str, table_key: str, field: str) -> bool:
    """Check whether a field already exists inside a table entry."""
    import re
    key_pattern = re.compile(r"^( *)(" + re.escape(table_key) + r"):[ \t]*$", re.MULTILINE)
    m = key_pattern.search(text)
    if not m:
        return False
    table_indent = len(m.group(1))
    content_indent = table_indent + 2

    # Look for the field at content_indent between table key and next sibling or EOF
    field_pattern = re.compile(
        r"^( {" + str(content_indent) + r"})" + re.escape(field) + r":[ \t]*$",
        re.MULTILINE,
    )
    # Search only within the table's block (until the next key at table_indent or less)
    next_sibling = re.compile(
        r"^\s{0," + str(table_indent) + r"}\S",
        re.MULTILINE,
    )
    end_m = next_sibling.search(text, m.end() + 1)
    search_end = end_m.start() if end_m else len(text)
    return bool(field_pattern.search(text, m.end(), search_end))


def process_paperconfig(path: Path, apply: bool, _yaml=None) -> bool:
    """Process one paperconfig using text-level insertion to preserve formatting."""

    dir_name = path.parent.name
    paper_mapping = MAPPINGS.get(dir_name)
    if paper_mapping is None:
        return False

    original = path.read_text(encoding="utf-8")
    text = original

    # Collect all tables of type csv and their mappings
    # Find all table keys in the file
    table_key_re = re.compile(r"^( +)(\w[\w_]*):[ \t]*$", re.MULTILINE)
    changed_tables = []

    # We need to apply insertions from bottom to top to preserve positions
    insertions = []  # list of (insert_pos, insert_text, table_key, added_fields)

    # Determine which tables are type csv
    STRUCTURAL_KEYS = {
        "supplementary_materials", "statistical_analyses", "environmental_conditions",
        "publication", "id_columns", "product_columns",
    }
    all_tables = {}
    for m in table_key_re.finditer(text):
        tkey = m.group(2)
        if tkey in STRUCTURAL_KEYS:
            continue
        # Check if there's a `type: csv` shortly after
        snippet = text[m.end():m.end() + 200]
        if re.search(r"type:\s*csv\b", snippet):
            all_tables[tkey] = m.start()

    for table_key in all_tables:
        fields = paper_mapping.get(table_key) or paper_mapping.get("_all_csv")
        if fields is None:
            continue

        fields_to_add = {}
        if "id_columns" in fields and not already_has_field(text, table_key, "id_columns"):
            fields_to_add["id_columns"] = fields["id_columns"]
        if "product_columns" in fields and not already_has_field(text, table_key, "product_columns"):
            fields_to_add["product_columns"] = fields["product_columns"]

        if not fields_to_add:
            continue

        insert_pos, content_indent = find_table_statistical_analyses(text, table_key)
        if insert_pos is None:
            print(f"  WARNING: could not find statistical_analyses for table '{table_key}'")
            continue

        insert_text = ""
        if "id_columns" in fields_to_add:
            insert_text += format_id_columns_yaml(content_indent, fields_to_add["id_columns"])
        if "product_columns" in fields_to_add:
            insert_text += format_product_columns_yaml(content_indent, fields_to_add["product_columns"])

        insertions.append((insert_pos, insert_text, table_key, list(fields_to_add.keys())))

    if not insertions:
        return False

    # Apply insertions in reverse order (bottom to top) to preserve positions
    insertions.sort(key=lambda x: x[0], reverse=True)
    for insert_pos, insert_text, table_key, added_fields in insertions:
        text = text[:insert_pos] + insert_text + text[insert_pos:]
        changed_tables.append((table_key, added_fields))

    print(f"\n{path}")
    for tbl, added_fields in sorted(changed_tables, key=lambda x: original.find(x[0])):
        print(f"  {tbl}: added {', '.join(added_fields)}")

    if apply:
        path.write_text(text, encoding="utf-8")

    return True


def main():
    parser = argparse.ArgumentParser(description="Add id_columns/product_columns to paperconfigs")
    parser.add_argument("--apply", action="store_true", help="Write changes (default: dry-run)")
    parser.add_argument("--paper", help="Process only this paper dir name")
    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    txt_file = project_root / "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
    paths = [project_root / line.strip() for line in txt_file.read_text().splitlines() if line.strip()]

    if args.paper:
        paths = [p for p in paths if args.paper.lower() in p.parent.name.lower()]

    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"=== update_paperconfig_id_columns ({mode}) ===")

    total_changed = 0
    for path in paths:
        if not path.exists():
            print(f"WARNING: {path} not found")
            continue
        if process_paperconfig(path, args.apply):
            total_changed += 1

    print(f"\nDone. {total_changed} paperconfig(s) {'updated' if args.apply else 'would be updated'}.")
    if not args.apply:
        print("Run with --apply to write changes.")


if __name__ == "__main__":
    main()
