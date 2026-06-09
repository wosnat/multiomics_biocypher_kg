"""Pure-Python PSORTb subcellular-localization vocabulary + helpers.

PSORTb v3.0.3 emits one of a closed 6-value set per protein. Five are real
localizations; ``Unknown`` is the no-signal sentinel (assigned when no class
scores above PSORTb's confidence floor). The KG models the five real classes as
``SubcellularLocalization`` ontology nodes and skips ``Unknown`` — absence of a
``Gene_has_subcellular_localization`` edge encodes "no confident localization".

This module is intentionally file-I/O-free so it is unit-testable in isolation.

Public API:
    LOCALIZATION_VOCAB — dict of raw call → display name (real classes only).
    UNKNOWN_SENTINEL    — the literal string PSORTb uses for no-signal calls.
    is_kept(call)       — True when ``call`` is a real (non-sentinel) localization.
    display_name(call)  — human-readable name for a real localization.
"""
from __future__ import annotations

# The no-signal sentinel — kept verbatim in gene_annotations_merged.json so the
# field is round-trippable, but NOT emitted as a node or edge.
UNKNOWN_SENTINEL = "Unknown"

# Closed vocabulary of the FIVE real PSORTb localizations (Gram-negative model).
# Order = broad → membrane → periplasm → outer → secreted (display only; the
# ontology is flat — level 0 for all).
LOCALIZATION_VOCAB: dict[str, str] = {
    "Cytoplasmic": "Cytoplasmic",
    "CytoplasmicMembrane": "Cytoplasmic membrane",
    "Periplasmic": "Periplasmic",
    "OuterMembrane": "Outer membrane",
    "Extracellular": "Extracellular",
}


def is_kept(call: str | None) -> bool:
    """True when ``call`` is a real localization the KG keeps (not the sentinel)."""
    return bool(call) and call in LOCALIZATION_VOCAB


def display_name(call: str) -> str:
    """Human-readable display name for a real localization call.

    Falls back to the raw call for anything outside the vocabulary (defensive —
    callers should gate on ``is_kept`` first).
    """
    return LOCALIZATION_VOCAB.get(call, call)
