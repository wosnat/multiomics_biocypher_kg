"""Named transform functions and registry for annotation build scripts."""

from __future__ import annotations

import re
from typing import Any


def _tx_first_token_space(value: str) -> str:
    """Return first whitespace-separated token: 'dnaN rps3' → 'dnaN'."""
    if not isinstance(value, str) or not value.strip():
        return ""
    return value.strip().split()[0]


def _tx_add_go_prefix(value: str) -> str:
    """Prepend 'GO:' to bare 7-digit numeric IDs: '0009360' → 'GO:0009360'.
    If already 'GO:...' return as-is.
    """
    s = str(value).strip()
    if not s or s == "-":
        return ""
    if s.startswith("GO:"):
        return s
    if re.match(r"^\d{7}$", s):
        return f"GO:{s}"
    return s  # leave non-GO ontology terms unchanged


def _tx_strip_function_prefix(value: str) -> str:
    """Strip 'FUNCTION: ' prefix from UniProt cc_function text."""
    if not isinstance(value, str):
        return ""
    return re.sub(r"^FUNCTION:\s*", "", value.strip(), flags=re.IGNORECASE)


def _tx_strip_prefix_ko(value: str) -> str:
    """Strip 'ko:' prefix from KEGG KO IDs: 'ko:K02710' → 'K02710'."""
    return re.sub(r"^ko:", "", str(value).strip(), flags=re.IGNORECASE)


def _tx_extract_go_from_pipe(value: str) -> str:
    """Extract GO ID from 'term_name|NNNNNNN||evidence' format.

    'DNA replication|0006260||IEA' → 'GO:0006260'
    Falls back to _tx_add_go_prefix if no pipe found.
    """
    s = str(value).strip()
    if not s or s == "-":
        return ""
    if "|" in s:
        parts = s.split("|")
        goid_candidate = parts[1].strip()
        if re.match(r"^\d{7}$", goid_candidate):
            return f"GO:{goid_candidate}"
    # Already GO: prefixed or bare digit?
    return _tx_add_go_prefix(s)


def _tx_extract_pfam_ids(value: str) -> list[str]:
    """Keep only PF* tokens from a comma-separated domain list.

    'TIGR00663,PF00712,IPR022634' → ['PF00712']
    """
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(",")
            if v.strip().startswith("PF")]


# Map from YAML transform name → function
_TRANSFORMS: dict[str, Any] = {
    "first_token_space": _tx_first_token_space,
    "add_go_prefix": _tx_add_go_prefix,
    "strip_function_prefix": _tx_strip_function_prefix,
    "strip_prefix_ko": _tx_strip_prefix_ko,
    "extract_go_from_pipe": _tx_extract_go_from_pipe,
    "extract_pfam_ids": _tx_extract_pfam_ids,
}
