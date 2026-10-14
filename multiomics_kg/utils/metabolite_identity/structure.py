"""Stage G: structure_kind, B2 scrub, generalized-formula fill (spec §3.3)."""
from __future__ import annotations

import re
from collections import defaultdict

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.names import (
    is_generic_formula, is_polymer_formula, parse_elements)

STRUCTURE_KINDS = ("defined", "generic", "polymer")
_POLYMER_NAMES = re.compile(
    r"^(poly|chitosan|chitin|cellulose|starch|glycogen|amylose|peptidoglycan)|\)n$", re.I)


def classify_structure(entry: dict, chebi_term_id: str | None, chebi: ChebiIndex) -> str:
    formula = entry.get("formula")
    if formula is None and chebi_term_id:
        formula = (chebi.term(chebi_term_id) or {}).get("generalized_formula")
    if is_polymer_formula(formula) or _POLYMER_NAMES.search(entry.get("name") or ""):
        return "polymer"
    if is_generic_formula(formula):
        return "generic"
    if chebi_term_id and chebi.is_class(chebi_term_id):
        return "generic"
    return "defined"


def _is_placeholder(formula: str | None) -> bool:
    return is_generic_formula(formula) or is_polymer_formula(formula)


def apply_structure(table, chebi, chebi_of):
    filled, scrubbed = [], []
    for node, e in table.items():
        term_id = chebi_of.get(node)
        gformula = (chebi.term(term_id) or {}).get("generalized_formula") if term_id else None
        kind = classify_structure(e, term_id, chebi)
        e["structure_kind"] = kind
        if kind != "defined":
            concrete = bool(e.get("formula")) and not _is_placeholder(e["formula"])
            if concrete or e.get("inchikey") or e.get("smiles"):
                # B2: drop the representative structure, keep its element set (spec §3.3).
                e["_elements_keep"] = parse_elements(e.get("formula"))
                e["inchikey"] = e["smiles"] = None
                if concrete:
                    e["formula"] = gformula if _is_placeholder(gformula) else None
                scrubbed.append(node)
            e["mass"] = None
        if e.get("formula") is None and gformula:
            e["formula"] = gformula
            filled.append(node)
        if _is_placeholder(e.get("formula")):
            e["mass"] = None
    by_key = defaultdict(list)
    for node, e in table.items():
        if e.get("inchikey"):
            by_key[e["inchikey"]].append(node)
    dup = {k: v for k, v in by_key.items() if len(v) > 1}
    if dup:
        raise ValueError(f"InChIKey on more than one Metabolite after the B2 scrub: {dup}")
    counts = defaultdict(int)
    for e in table.values():
        keep = e.pop("_elements_keep", None)
        e["elements"] = sorted(set(parse_elements(e.get("formula"))) | set(keep or []))
        counts[e["structure_kind"]] += 1
    return {"scrubbed": sorted(scrubbed), "filled": sorted(filled), "counts": dict(counts)}
