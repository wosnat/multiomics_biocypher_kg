"""Stage H: chebi_id pick, id arrays, all_identifiers, identity_basis, name (spec §3.2)."""
from __future__ import annotations

import re

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.names import normalize_hmdb

_KEGG = "kegg.compound:"


def _bare(cid: str) -> str:
    return cid.split(":", 1)[1]


def _entity_rank(chebi: ChebiIndex):
    return lambda c: (chebi.charge(c) != 0, int(_bare(c)))


# Isotopologue names (spec §1.1 B1b): "succinic acid-d4", "UREA C 14", "calcium-43(2+)", "[13C]…", "(15N)…".
# Deliberately narrow: must NOT match "1H-indole", "Vitamin B12", "(2E)-hexenal", "(3S)-…".
_ISOTOPE = re.compile(
    r"(?-i:-d\d+\b)|isotop|deuter|[\[(]\d{2,3}[CNOHSP][\])]|\[[23]H\]|\b(?:C|N|O|H) 1[0-9]\b|[a-z]+-\d{2,3}\(\d?[+-]\)", re.I)


def is_isotopologue_name(name: str | None) -> bool:
    return bool(_ISOTOPE.search(name or ""))


def _usable(p: str | None, node_id: str, chebi: ChebiIndex) -> bool:
    """An entity term for this node: not a class, not an isotopologue, not another KEGG compound's term."""
    if p is None or chebi.is_class(p):
        return False
    t = chebi.term(p) or {}
    if is_isotopologue_name(t.get("name")):
        return False
    xref = t.get("kegg", [])
    return not xref or (node_id.startswith(_KEGG) and node_id[len(_KEGG):] in xref)


def pick_chebi_id(node_id, entry, chebi):
    if node_id.startswith("chebi:"):
        p = chebi.primary(node_id)
        return _bare(p) if p else _bare(node_id)
    cands = []
    if node_id.startswith(_KEGG):
        cands = [c for c in chebi.kegg_terms(node_id[len(_KEGG):]) if not chebi.is_class(c)]
    if not cands:
        aliases = list(entry.get("mnx_chebi_aliases") or []) + ([entry["chebi_id"]] if entry.get("chebi_id") else [])
        cands = sorted({p for p in (chebi.primary(a) for a in aliases) if _usable(p, node_id, chebi)})
    if not cands:
        return None
    return _bare(sorted(cands, key=_entity_rank(chebi))[0])


def build_identity_fields(table, identity: IdentityResult, rekeys, chebi, old_entries):
    rekeyed_into: dict[str, list[str]] = {}
    for old, new in rekeys.items():
        rekeyed_into.setdefault(new, []).append(old)
    index: dict[str, str] = {}
    for node in sorted(table):
        e = table[node]
        if node.startswith("chebi:"):
            t = chebi.term(node)
            if t and t.get("name"):
                e["name"] = t["name"]
        group = identity.members.get(node, [node])
        olds = sorted(set(group) | set(rekeyed_into.get(node, [])))
        kegg = [node[len(_KEGG):]] if node.startswith(_KEGG) else []
        kegg += sorted({m[len(_KEGG):] for m in olds if m.startswith(_KEGG)} - set(kegg))
        e["kegg_compound_id"] = kegg[0] if kegg else None
        e["kegg_compound_ids"] = kegg or None
        chebi_id = pick_chebi_id(node, e, chebi)
        e["chebi_id"] = chebi_id
        chebi_terms = set()
        for m in olds:
            src = old_entries.get(m, {})
            if m.startswith("chebi:"):
                p = chebi.primary(m)
                if p:
                    chebi_terms.add(p)
            if m.startswith(_KEGG):
                chebi_terms |= {c for c in chebi.kegg_terms(m[len(_KEGG):]) if not chebi.is_class(c)}
            # A secondary MNX ChEBI alias maps onto its primary term: every ChEBI id we publish is primary.
            for a in list(src.get("mnx_chebi_aliases") or []) + ([src["chebi_id"]] if src.get("chebi_id") else []):
                p = chebi.primary(a)
                if _usable(p, node, chebi):
                    chebi_terms.add(p)
        if chebi_id is None and chebi_terms:
            chebi_id = _bare(sorted(chebi_terms, key=_entity_rank(chebi))[0])
            e["chebi_id"] = chebi_id
        ids = [chebi_id] if chebi_id else []
        ids += sorted({_bare(c) for c in chebi_terms} - set(ids), key=int)
        e["chebi_ids"] = ids or None
        mnx = [e["mnxm_id"]] if e.get("mnxm_id") else []
        mnx += sorted({old_entries.get(m, {}).get("mnxm_id") for m in olds} - {None} - set(mnx))
        e["mnxm_id"] = mnx[0] if mnx else None
        e["mnxm_ids"] = mnx if len(mnx) > 1 else None
        e["hmdb_id"] = normalize_hmdb(e.get("hmdb_id"))
        member_hmdb = sorted({h for h in (normalize_hmdb(old_entries.get(m, {}).get("hmdb_id")) for m in olds) if h})
        if e["hmdb_id"] is None and member_hmdb:
            e["hmdb_id"] = member_hmdb[0]
        ident = {node} | set(olds)
        ident |= {f"{_KEGG}{k}" for k in kegg}
        for c in chebi_terms | ({f"CHEBI:{chebi_id}"} if chebi_id else set()):
            ident.add(f"chebi:{_bare(c)}")
            ident |= {f"chebi:{_bare(a)}" for a in chebi.secondary_ids(c)}
        ident |= {f"mnx:{m}" for m in mnx}
        hmdbs = {normalize_hmdb(old_entries.get(m, {}).get("hmdb_id")) for m in olds} | {e["hmdb_id"]}
        ident |= {f"hmdb:{h}" for h in hmdbs if h}
        if e.get("inchikey"):
            ident.add(f"inchikey:{e['inchikey'].upper()}")
        e["all_identifiers"] = sorted(ident)
        basis = set(identity.basis.get(node, []))
        if rekeyed_into.get(node):
            basis.add("chebi_kegg_xref")
        e["identity_basis"] = sorted(basis) or None
        for x in e["all_identifiers"]:
            if index.setdefault(x, node) != node:
                raise ValueError(f"identifier {x} on two Metabolite nodes: {index[x]} and {node}")
    return index
