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


def _tx_extract_go_from_brackets(value: str) -> str:
    """Extract GO ID from UniProt bracket notation.

    'DNA polymerase III complex [GO:0009360]' → 'GO:0009360'
    Returns empty string if no GO term found.
    """
    s = str(value).strip()
    if not s or s == "-":
        return ""
    parts = s.split("GO:")
    if len(parts) < 2:
        return ""
    return "GO:" + parts[-1].rstrip("]").strip()


def _tx_extract_pfam_ids(value: str) -> list[str]:
    """Keep only PF* tokens from a comma-separated domain list.

    'TIGR00663,PF00712,IPR022634' → ['PF00712']
    """
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(",")
            if v.strip().startswith("PF")]


_ECO_PATTERN = re.compile(r'\s*\{ECO:[^}]+\}[.,]?\s*')


def _tx_clean_function_description(value: str) -> str:
    """Strip 'FUNCTION:' prefix and inline ECO evidence tags."""
    if not isinstance(value, str):
        return ""
    s = re.sub(r'^FUNCTION:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip('.')
    return s


def _tx_clean_catalytic_activity(value: str) -> str:
    """Strip 'CATALYTIC ACTIVITY:' prefix and ECO tags from one reaction chunk."""
    if not isinstance(value, str):
        return ""
    s = re.sub(r'^CATALYTIC ACTIVITY:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip(';').strip()
    return s


def _tx_extract_cofactor_name(value: str) -> str:
    """'COFACTOR: Name=FMN; Xref=…' → 'FMN'"""
    if not isinstance(value, str):
        return ""
    m = re.match(r'COFACTOR:\s*Name=([^;]+)', value.strip(), re.IGNORECASE)
    return m.group(1).strip() if m else ""


def _tx_extract_pathway_name(value: str) -> str:
    """'PATHWAY: Energy metabolism; oxidative phosphorylation. {ECO:…}.' → clean string"""
    if not isinstance(value, str):
        return ""
    s = re.sub(r'^PATHWAY:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip('.')
    return s


def _tx_extract_tm_range(value: str) -> str:
    """'TRANSMEM 32..50; /note="Helical"; …' → '32..50'"""
    if not isinstance(value, str):
        return ""
    m = re.search(r'TRANSMEM\s+(\d+\.\.\d+)', value)
    return m.group(1) if m else ""


def _tx_extract_signal_range(value: str) -> str:
    """'SIGNAL 1..26; /evidence="…"' → '1..26'"""
    if not isinstance(value, str):
        return ""
    m = re.search(r'SIGNAL\s+(\d+\.\.\d+)', value)
    return m.group(1) if m else ""


# Map from YAML transform name → function
_TRANSFORMS: dict[str, Any] = {
    "first_token_space": _tx_first_token_space,
    "add_go_prefix": _tx_add_go_prefix,
    "strip_function_prefix": _tx_strip_function_prefix,
    "strip_prefix_ko": _tx_strip_prefix_ko,
    "extract_go_from_pipe": _tx_extract_go_from_pipe,
    "extract_pfam_ids": _tx_extract_pfam_ids,
    "extract_go_from_brackets": _tx_extract_go_from_brackets,
    "clean_function_description": _tx_clean_function_description,
    "clean_catalytic_activity": _tx_clean_catalytic_activity,
    "extract_cofactor_name": _tx_extract_cofactor_name,
    "extract_pathway_name": _tx_extract_pathway_name,
    "extract_tm_range": _tx_extract_tm_range,
    "extract_signal_range": _tx_extract_signal_range,
}
