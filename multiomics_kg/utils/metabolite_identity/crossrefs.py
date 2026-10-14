"""Stage H: chebi_id pick, id arrays, all_identifiers, identity_basis, name (spec §3.2)."""
from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Collection

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.canonical import IdentityCollisionError, IdentityResult
from multiomics_kg.utils.metabolite_identity.names import normalize_hmdb

_KEGG = "kegg.compound:"


def _bare(cid: str) -> str:
    return cid.split(":", 1)[1]


def _entity_rank(chebi: ChebiIndex):
    return lambda c: (chebi.charge(c) != 0, int(_bare(c)))


# Isotopologue names (spec §1.1 B1b): "succinic acid-d4", "UREA C 14", "calcium-43(2+)", "[13C]…", "(15N)…".
# Deliberately narrow: must NOT match "1H-indole", "Vitamin B12", "(2E)-hexenal", "(3S)-…",
# "(25S)-…" / "(20S)-…" stereo locants (the bracketed element is upper-case C/N/O/H only, case-
# sensitive) or "Vitamin C 12" (only the labelling isotopes 13C, 14C, 15N, 17O, 18O spelled apart).
_ISOTOPE = re.compile(
    r"(?-i:-d\d+\b)|isotop|deuter|(?-i:[\[(]\d{2,3}[CNOH][\])])|\[[23]H\]"
    r"|\b(?:C 1[34]|N 15|O 1[78])\b|[a-z]+-\d{2,3}\(\d?[+-]\)", re.I)


def is_isotopologue_name(name: str | None) -> bool:
    return bool(_ISOTOPE.search(name or ""))


def shared_chebi_terms(table: dict[str, dict], identity: IdentityResult,
                       chebi: ChebiIndex) -> dict[str, list[str]]:
    """Entity terms that cross-reference >= 2 written KEGG nodes that remain separate.

    A KEGG id counts for the canonical node that carries it (the node itself or a merged
    member), so a term whose KEGG xrefs all landed on one canonical node is not shared.
    A shared term names neither node: the KEGG<->KEGG guard refused to merge them (stereo or
    specific-vs-generic pairs), so publishing it on both would put one identifier on two nodes.
    """
    node_of: dict[str, str] = {}
    for node in sorted(table):
        for m in [node, *identity.members.get(node, [])]:
            if m.startswith(_KEGG):
                node_of.setdefault(m[len(_KEGG):], node)
    by_term: dict[str, set[str]] = defaultdict(set)
    for kid, node in node_of.items():
        for t in chebi.kegg_terms(kid):
            by_term[t].add(node)
    return {t: sorted(nodes) for t, nodes in sorted(by_term.items())
            if len(nodes) > 1 and not chebi.is_class(t)}


def _usable(p: str | None, node_id: str, chebi: ChebiIndex, shared: Collection[str] = ()) -> bool:
    """An entity term for this node: not a class, not an isotopologue, not another KEGG
    compound's term, not a term shared by two separate KEGG nodes."""
    if p is None or p in shared or chebi.is_class(p):
        return False
    t = chebi.term(p) or {}
    if is_isotopologue_name(t.get("name")):
        return False
    xref = t.get("kegg", [])
    return not xref or (node_id.startswith(_KEGG) and node_id[len(_KEGG):] in xref)


def node_olds(node: str, identity: IdentityResult, rekeys: dict[str, str]) -> list[str]:
    """Every pre-canonical id behind ``node``: its merge group plus ids stage D re-keyed onto it."""
    rekeyed_into = [old for old, new in rekeys.items() if new == node]
    return sorted(set(identity.members.get(node, [node])) | set(rekeyed_into))


def chebi_pick(node_id, entry, chebi, *, shared: Collection[str] = (), olds=None, old_entries=None):
    """The one ChEBI pick shared by stage G (structure) and stage H (cross-references).

    Returns ``(chebi_id, entity_pool)``: ``entity_pool`` is every usable entity term of the node
    (not a class, not shared, not an isotopologue) from its own KEGG xrefs, its members' /
    re-keyed ids (their KEGG xrefs and ChEBI terms) and MNX ChEBI aliases (the node's and its
    members'). ``chebi_id`` is the node's own entity term for a ``chebi:`` node; else the best
    entity term, preferring the node's own KEGG xrefs, then its MNX aliases, then the rest of
    the pool. Only when the node has no entity term at all may a class term be ``chebi_id`` —
    the node is then ``generic`` (stage G reads that class term)."""
    olds = list(olds) if olds is not None else [node_id]
    old_entries = old_entries or {}
    rank = _entity_rank(chebi)

    def ok(p):
        if p is None or p in shared:
            return False
        return not is_isotopologue_name((chebi.term(p) or {}).get("name"))

    def entity(p):
        return ok(p) and not chebi.is_class(p)

    own_kegg, alias_terms, member_terms, classes = [], [], set(), []
    if node_id.startswith(_KEGG):
        for c in chebi.kegg_terms(node_id[len(_KEGG):]):
            if entity(c):
                own_kegg.append(c)
            elif ok(c):
                classes.append(c)

    def aliases(e):
        return list(e.get("mnx_chebi_aliases") or []) + ([e["chebi_id"]] if e.get("chebi_id") else [])
    for a in aliases(entry):
        p = chebi.primary(a)
        if _usable(p, node_id, chebi, shared):
            alias_terms.append(p)
    for m in olds:
        if m.startswith("chebi:"):
            p = chebi.primary(m)
            if entity(p):
                member_terms.add(p)
            elif ok(p):
                classes.append(p)
        if m.startswith(_KEGG) and m != node_id:
            for c in chebi.kegg_terms(m[len(_KEGG):]):
                if entity(c):
                    member_terms.add(c)
                elif ok(c):
                    classes.append(c)
        for a in aliases(old_entries.get(m, {})):
            p = chebi.primary(a)
            if _usable(p, node_id, chebi, shared):
                member_terms.add(p)
    pool = set(own_kegg) | set(alias_terms) | member_terms
    if node_id.startswith("chebi:"):
        p = chebi.primary(node_id)
        return (_bare(p) if p else _bare(node_id)), pool
    for tier in (own_kegg, alias_terms, member_terms):
        if tier:
            return _bare(sorted(set(tier), key=rank)[0]), pool
    if classes:
        return _bare(sorted(set(classes), key=lambda c: int(_bare(c)))[0]), pool
    return None, pool


def pick_chebi_id(node_id, entry, chebi, shared: Collection[str] = (), *, olds=None, old_entries=None):
    return chebi_pick(node_id, entry, chebi, shared=shared, olds=olds, old_entries=old_entries)[0]


def build_identity_fields(table, identity: IdentityResult, rekeys, chebi, old_entries,
                          shared: Collection[str] = ()):
    """Stage H. ``shared`` (see ``shared_chebi_terms``) are ChEBI terms published on no node."""
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
        chebi_id, chebi_terms = chebi_pick(node, e, chebi, shared=shared, olds=olds, old_entries=old_entries)
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
                raise IdentityCollisionError(f"identifier {x} on two Metabolite nodes: {index[x]} and {node}")
    return index
