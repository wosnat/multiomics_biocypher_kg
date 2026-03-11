"""Named transform functions and registry for annotation build scripts."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


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


def _tx_split_cog_category(value: str) -> list[str]:
    """Split a multi-letter COG category string into a list of single-char codes.

    'LU' → ['L', 'U']  |  'S' → ['S']  |  '-' → []
    """
    if not isinstance(value, str):
        return []
    s = value.strip()
    if not s or s == "-":
        return []
    return list(s)


def _tx_extract_pfam_ids(value: str) -> list[str]:
    """Keep only PF* tokens from a comma-separated domain list.

    'TIGR00663,PF00712,IPR022634' → ['PF00712']
    """
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(",")
            if v.strip().startswith("PF")]


def _tx_extract_pfam_names(value: str) -> list[str]:
    """Keep only non-PF* tokens from a comma-separated domain list (shortnames).

    'PF00712,DNA_pol3_beta,DNA_pol3_beta_2' → ['DNA_pol3_beta', 'DNA_pol3_beta_2']
    Only relevant for eggnog PFAMs column which mixes PF* IDs and shortnames.
    """
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(",")
            if v.strip() and not v.strip().startswith("PF")]


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


# ── EC number normalization ──────────────────────────────────────────────────

# Lazy-loaded transfer map: {old_ec: [new_ec, ...]}  (empty list = deleted)
_EC_TRANSFER_MAP: dict[str, list[str]] | None = None


def _resolve_ec_chain(
    ec: str, transfer_map: dict[str, list[str]], visited: set[str] | None = None,
) -> list[str]:
    """Follow transfer chain for one EC to its final current successor(s).

    Returns list of current EC numbers, or [] if the chain ends at a deleted entry.
    """
    if visited is None:
        visited = set()
    if ec in visited:
        return [ec]  # cycle guard
    if ec not in transfer_map:
        return [ec]  # current EC — not obsolete
    visited.add(ec)
    successors = transfer_map[ec]
    if not successors:
        return []  # deleted entry
    result = []
    for s in successors:
        result.extend(_resolve_ec_chain(s, transfer_map, visited.copy()))
    return result


def _build_ec_transfer_map() -> dict[str, list[str]]:
    """Build old→new EC mapping from cached Expasy data.

    Returns dict where keys are obsolete EC numbers and values are lists of
    final current successor EC numbers (empty list for deleted entries).
    Chained transfers (A→B→C) are resolved so A maps directly to C.
    """
    cache_path = Path(__file__).resolve().parents[3] / "cache" / "data" / "ec" / "ec_data.json"
    if not cache_path.exists():
        logger.warning(
            "EC cache not found at %s; skipping EC normalization. "
            "Run create_knowledge_graph.py first to populate the cache.",
            cache_path,
        )
        return {}

    with open(cache_path) as f:
        data = json.load(f)

    # Pass 1: build raw transfer map (may contain intermediate entries)
    transfer_map: dict[str, list[str]] = {}
    for ec_num, info in data.get("enzymes", {}).items():
        if not isinstance(info, dict):
            continue
        desc = info.get("de", "")
        if desc.startswith("Transferred entry:"):
            # Parse successors from "Transferred entry: X.Y.Z.W" or
            # "Transferred entry: X.Y.Z.W, A.B.C.D and E.F.G.H"
            rest = desc[len("Transferred entry:"):].strip().rstrip(".")
            # Split on ", " and " and "
            parts = re.split(r",\s*|\s+and\s+", rest)
            successors = [p.strip() for p in parts if re.match(r"^\d+\.[\d\-]+\.[\d\-]+\.[\d\-]+$", p.strip())]
            transfer_map[ec_num] = successors
        elif desc.startswith("Deleted"):
            transfer_map[ec_num] = []

    # Pass 2: resolve chains so every entry points to final current EC(s)
    for ec_num in list(transfer_map.keys()):
        transfer_map[ec_num] = _resolve_ec_chain(ec_num, transfer_map)

    logger.info("EC transfer map: %d transferred, %d deleted",
                sum(1 for v in transfer_map.values() if v),
                sum(1 for v in transfer_map.values() if not v))
    return transfer_map


def _tx_normalize_ec(value: str) -> str | list[str]:
    """Remap obsolete/transferred EC numbers to current equivalents.

    Returns:
      - Original EC string if it's current (not in transfer map)
      - Successor EC string(s) if transferred (list for multi-successor)
      - Empty string if deleted
    """
    global _EC_TRANSFER_MAP
    if _EC_TRANSFER_MAP is None:
        _EC_TRANSFER_MAP = _build_ec_transfer_map()

    s = str(value).strip()
    if not s or s == "-":
        return ""

    if s not in _EC_TRANSFER_MAP:
        return s  # current EC, pass through

    successors = _EC_TRANSFER_MAP[s]
    if not successors:
        return ""  # deleted entry
    if len(successors) == 1:
        return successors[0]
    return successors  # multi-successor: list


# Map from YAML transform name → function
_TRANSFORMS: dict[str, Any] = {
    "first_token_space": _tx_first_token_space,
    "add_go_prefix": _tx_add_go_prefix,
    "strip_function_prefix": _tx_strip_function_prefix,
    "strip_prefix_ko": _tx_strip_prefix_ko,
    "extract_go_from_pipe": _tx_extract_go_from_pipe,
    "extract_pfam_ids": _tx_extract_pfam_ids,
    "extract_pfam_names": _tx_extract_pfam_names,
    "extract_go_from_brackets": _tx_extract_go_from_brackets,
    "clean_function_description": _tx_clean_function_description,
    "clean_catalytic_activity": _tx_clean_catalytic_activity,
    "extract_cofactor_name": _tx_extract_cofactor_name,
    "extract_pathway_name": _tx_extract_pathway_name,
    "extract_tm_range": _tx_extract_tm_range,
    "extract_signal_range": _tx_extract_signal_range,
    "split_cog_category": _tx_split_cog_category,
    "normalize_ec": _tx_normalize_ec,
}
