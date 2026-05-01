"""Pure-Python CAZy classification helpers.

CAZy hierarchy is now derived in-process — no file I/O. The 6 classes are
hardcoded; family + subfamily IDs are parsed from observed eggNOG annotations
in `multiomics_kg/adapters/cazy_adapter.py`.

Public API:
    CAZY_CLASSES — dict of class code → display name (immutable map).
    parse_cazy_id(token) — split a CAZy token into (family, subfamily | None).
    is_valid_cazy(value) — True when value is a recognized class / family / subfamily.
    cazy_ancestors(value) — root-to-parent ancestor chain ([class] or [class, family]).
"""
from __future__ import annotations

import re

CAZY_CLASSES: dict[str, str] = {
    "GH":  "Glycoside Hydrolases",
    "GT":  "GlycosylTransferases",
    "PL":  "Polysaccharide Lyases",
    "CE":  "Carbohydrate Esterases",
    "AA":  "Auxiliary Activities",
    "CBM": "Carbohydrate-Binding Modules",
}

_CAZY_FAMILY_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)(?:_(\d+))?$")


def parse_cazy_id(token: str) -> tuple[str, str | None] | None:
    """Return (family_id, subfamily_id_or_None) or None for malformed tokens.

    Examples:
        'GH13'    → ('GH13', None)
        'GH13_5'  → ('GH13', 'GH13_5')
        'CBM48'   → ('CBM48', None)
        'invalid' → None
    """
    if not token:
        return None
    m = _CAZY_FAMILY_RE.match(token.strip())
    if not m:
        return None
    cls, fam_num, sub_num = m.groups()
    family = f"{cls}{fam_num}"
    subfamily = f"{family}_{sub_num}" if sub_num else None
    return family, subfamily


def is_valid_cazy(value: str) -> bool:
    """True for a recognized class code, family ID, or subfamily ID."""
    if not value:
        return False
    if value in CAZY_CLASSES:
        return True
    return parse_cazy_id(value) is not None


def cazy_ancestors(value: str) -> list[str]:
    """Root-to-parent ancestor list. Empty for class-level or unknown values."""
    if not value or value in CAZY_CLASSES:
        return []
    parsed = parse_cazy_id(value)
    if parsed is None:
        return []
    family, subfamily = parsed
    cls = _CAZY_FAMILY_RE.match(family).group(1)
    if subfamily and value == subfamily:
        return [cls, family]
    return [cls]
