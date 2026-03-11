"""
Utility for parsing the Cyanorak functional role hierarchy.

The hierarchy is stored in ``data/cyanorak_roles.csv``, downloaded from
the Cyanorak database.

CSV columns: primary role id, primary role, secondary role id, secondary role,
sub role id, sub role.

The hierarchy has three levels:
- Primary role (letter code like "A", "B", etc.) — parent is None (root)
- Secondary role (code like "A.1", "B.5", etc.) — parent is primary role id
- Sub role (code like "B.5.1", "D.1.3", etc.) — parent is secondary role id
"""

import csv
from pathlib import Path


def parse_cyanorak_role_tree(csv_path: Path) -> dict[str, dict]:
    """
    Parse the Cyanorak roles CSV file into a hierarchy dict.

    Args:
        csv_path: path to ``data/cyanorak_roles.csv``

    Returns:
        ``{code: {"description": str, "parent": str | None}}``

    The returned dict includes all entries in the file: primary, secondary,
    and sub role levels.  The "Unclassified" entry (with ``-`` as primary
    role id) is skipped.
    """
    tree: dict[str, dict] = {}

    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            primary_id = row.get("primary role id", "").strip()
            primary_desc = row.get("primary role", "").strip()
            secondary_id = row.get("secondary role id", "").strip()
            secondary_desc = row.get("secondary role", "").strip()
            sub_id = row.get("sub role id", "").strip()
            sub_desc = row.get("sub role", "").strip()

            # Skip the "Unclassified" entry (primary role id is "-")
            if not primary_id or primary_id == "-":
                continue

            # Add primary role if not already present (root, no parent)
            if primary_id not in tree:
                tree[primary_id] = {
                    "description": primary_desc,
                    "parent": None,
                }

            # Add secondary role if present (parent is primary)
            if secondary_id and secondary_id != "-" and secondary_id not in tree:
                tree[secondary_id] = {
                    "description": secondary_desc,
                    "parent": primary_id,
                }

            # Add sub role if present (parent is secondary)
            if sub_id and sub_id != "-" and sub_id not in tree:
                tree[sub_id] = {
                    "description": sub_desc,
                    "parent": secondary_id,
                }

    return tree
