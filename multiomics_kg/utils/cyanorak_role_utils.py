"""
Utility for parsing the Cyanorak functional role hierarchy.

The hierarchy is stored in ``data/Prochlorococcus/cyanorak_roles.txt``,
copy-pasted from https://cyanorak.sb-roscoff.fr/cyanorak/help.html.

Format of the text file (after the "Cyanorak Roles" header):
- Code line: matches ``^[0-9A-Z]+(\\.\\d+)*$``
- Description: the next non-empty, non-``-`` line after a code line
- Leaf entries (no children in the tree) are followed by one or two ``-`` lines
- The "Unclassified" entry has no code → skipped

Parent derivation: strip the last ``.N`` segment from the code.
  "B.5.1" → parent "B.5"
  "B.5"   → parent "B"
  "B"     → parent None (root)
  "0.2"   → parent "0"
  "0"     → parent None (root)
"""

import re
from pathlib import Path

_CODE_RE = re.compile(r"^[0-9A-Z]+(\.\d+)*$")


def parse_cyanorak_role_tree(txt_path: Path) -> dict[str, dict]:
    """
    Parse the Cyanorak roles text file into a hierarchy dict.

    Args:
        txt_path: path to ``data/Prochlorococcus/cyanorak_roles.txt``

    Returns:
        ``{code: {"description": str, "parent": str | None}}``

    The returned dict includes all entries in the file: both leaf-level codes
    (directly assigned to genes) and intermediate parent codes.
    """
    lines = txt_path.read_text(encoding="utf-8").splitlines()

    tree: dict[str, dict] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1

        # Skip empty lines, dashes, and the header
        if not line or line == "-" or line == "Cyanorak Roles":
            continue

        # Check if this is a code line
        if _CODE_RE.match(line):
            code = line
            # Find the description: next non-empty, non-dash line
            description = ""
            while i < len(lines):
                desc_line = lines[i].strip()
                i += 1
                if desc_line and desc_line != "-":
                    description = desc_line
                    break

            parent = _derive_parent(code)
            tree[code] = {"description": description, "parent": parent}

    return tree


def _derive_parent(code: str) -> str | None:
    """Derive the parent code by stripping the last dot-segment."""
    if "." not in code:
        return None
    return code.rsplit(".", 1)[0]
