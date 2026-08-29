"""Shared rules for inferring JCVI TIGR roles from NCBIfam TIGR* family hits.

Three consumers must agree on the gate and the naming, so the rules live here:
- step-2 merge (``build_gene_annotations``): ``gene_category`` fill-only rule +
  ``[tigr_role_inferred]`` description lines;
- ``functional_annotation_adapter``: inferred ``Gene_has_tigr_role`` edges +
  two-level ``TigrRole`` nodes;
- ``ncbifam_adapter``: ``Ncbifam_family_has_tigr_role`` bridge (ungated —
  the bridge is ontology→ontology and records every archive link).

Design: docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md
"""

from __future__ import annotations

import json
import re
from pathlib import Path

# Only equivalog families ("same function in every member") transfer a role to
# a gene. subfamily / equivalog_domain / hypoth_equivalog / domain … never do —
# measured 2026-08-29: the gate cuts multi-role genes 375→15 and
# Cyanorak contradictions 7%→5.4% at the cost of ~1/3 of heterotroph reach.
EQUIVALOG_TYPES: frozenset[str] = frozenset({"equivalog"})

_SLUG_RE = re.compile(r"[^a-z0-9]+")
ROLE_SEP = " / "


def load_tigr_roles(cache_root: Path) -> dict | None:
    """Load ``<cache_root>/ncbifam/tigr_roles.json`` (step 9); ``None`` if absent."""
    path = Path(cache_root) / "ncbifam" / "tigr_roles.json"
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def mainrole_slug(mainrole: str) -> str:
    """``"Energy metabolism"`` → ``"energy_metabolism"`` (mainrole node local id)."""
    return _SLUG_RE.sub("_", mainrole.strip().lower()).strip("_")


def role_name(role_id: str, tigr_roles: dict) -> str:
    """Compound display name ``"<mainrole> / <sub1role>"`` (Cyanorak convention)."""
    r = tigr_roles["roles"][role_id]
    return f"{r['mainrole']}{ROLE_SEP}{r['sub1role']}" if r.get("sub1role") else r["mainrole"]


def split_role_name(name: str) -> tuple[str, str | None]:
    """Inverse of :func:`role_name` for compound names carried by Cyanorak."""
    if ROLE_SEP in name:
        main, sub = name.split(ROLE_SEP, 1)
        return main.strip(), sub.strip()
    return name.strip(), None


def inferred_roles(
    ncbifam_ids: list[str] | None,
    tigr_roles: dict,
    ncbifam_ref: dict,
) -> dict[str, list[str]]:
    """Equivalog-gated ``{role_id: [supporting TIGR accessions]}`` for one gene.

    An accession contributes iff ``ncbifam_ref[acc]["family_type"]`` is in
    :data:`EQUIVALOG_TYPES` AND ``tigr_roles["family_role"]`` maps it. Junk
    roles (156/157/…) are returned like any other — the nodes carry
    ``is_uninformative`` downstream; hiding them here would hide that a
    family is "hypothetical".
    """
    out: dict[str, list[str]] = {}
    if not ncbifam_ids or not tigr_roles or not ncbifam_ref:
        return out
    family_role = tigr_roles.get("family_role") or {}
    for acc in ncbifam_ids:
        if not acc or acc not in family_role:
            continue
        if (ncbifam_ref.get(acc) or {}).get("family_type") not in EQUIVALOG_TYPES:
            continue
        for role_id in family_role[acc]:
            out.setdefault(role_id, []).append(acc)
    return {k: sorted(v) for k, v in sorted(out.items())}
