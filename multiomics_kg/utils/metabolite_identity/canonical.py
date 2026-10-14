"""Stages E-F: merge links, union-find with the KEGG guard, winner rule, id rewrite (spec §3.1)."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from multiomics_kg.utils.chebi_utils import ChebiIndex

MERGE_BASES = ("mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual")
_KEGG = "kegg.compound:"


@dataclass(frozen=True, order=True)
class MergeLink:
    a: str
    b: str
    basis: str


def _link(x: str, y: str, basis: str) -> MergeLink:
    a, b = sorted((x, y))
    return MergeLink(a, b, basis)


@dataclass
class IdentityOverrides:
    merges: list[MergeLink]
    keep_separate: set[frozenset[str]]


def load_identity_overrides(path: Path) -> IdentityOverrides:
    if not path.exists():
        return IdentityOverrides([], set())
    doc = yaml.safe_load(path.read_text()) or {}
    merges = [_link(e["a"], e["b"], "manual") for e in (doc.get("merge") or [])]
    keep = {frozenset((e["a"], e["b"])) for e in (doc.get("keep_separate") or [])}
    clash = {frozenset((m.a, m.b)) for m in merges} & keep
    if clash:
        raise ValueError(f"metabolite_identity_overrides.yaml: pairs in both merge and "
                         f"keep_separate: {sorted(sorted(c) for c in clash)}")
    return IdentityOverrides(merges, keep)


def entity_terms(node_id: str, entry: dict, chebi: ChebiIndex) -> set[str]:
    if node_id.startswith(_KEGG):
        return set(chebi.kegg_terms(node_id[len(_KEGG):]))
    if node_id.startswith("chebi:"):
        p = chebi.primary(node_id)
        return {p} if p else set()
    if entry.get("chebi_id"):
        p = chebi.primary(entry["chebi_id"])
        return {p} if p else set()
    return set()


def collect_merge_links(table: dict[str, dict], chebi: ChebiIndex) -> list[MergeLink]:
    links: set[MergeLink] = set()
    by_mnx: dict[str, list[str]] = defaultdict(list)
    terms_of = {n: entity_terms(n, e, chebi) for n, e in table.items()}
    by_term: dict[str, list[str]] = defaultdict(list)
    for n, e in table.items():
        if e.get("mnxm_id"):
            by_mnx[e["mnxm_id"]].append(n)
        for t in terms_of[n]:
            by_term[t].append(n)
    for nodes in by_mnx.values():
        nodes = sorted(nodes)
        links.update(_link(nodes[0], o, "mnx_same_entity") for o in nodes[1:])
    for nodes in by_term.values():
        nodes = sorted(nodes)
        links.update(_link(nodes[0], o, "chebi_kegg_xref") for o in nodes[1:])
    for n, ts in terms_of.items():
        for t in ts:
            tt = chebi.term(t) or {}
            for rel, basis in (("conjugates", "chebi_conjugate"), ("tautomers", "chebi_tautomer")):
                for partner in tt.get(rel, []):
                    p = chebi.primary(partner)
                    for other in by_term.get(p, []):
                        if other != n:
                            links.add(_link(n, other, basis))
    return sorted(links)


@dataclass
class IdentityResult:
    canonical_of: dict[str, str]
    members: dict[str, list[str]]
    basis: dict[str, list[str]]
    member_basis: dict[str, str]
    rejected: list[MergeLink] = field(default_factory=list)
    refused_collisions: list[dict] = field(default_factory=list)
    stale_overrides: list[dict] = field(default_factory=list)


def _choose(members, reaction_refs, pathway_refs, chebi) -> str:
    kegg = [m for m in members if m.startswith(_KEGG)]
    if kegg:
        return sorted(kegg, key=lambda m: (-reaction_refs.get(m, 0), -pathway_refs.get(m, 0), m))[0]
    ch = [m for m in members if m.startswith("chebi:")]
    if ch:
        return sorted(ch, key=lambda m: (chebi.charge(m) != 0, int(m.split(":")[1])))[0]
    return sorted(members)[0]


def canonicalise(table, links, *, overrides, reaction_refs, pathway_refs, chebi, table_pins):
    parent = {n: n for n in table}
    members = {n: {n} for n in table}
    attach: dict[str, str] = {}
    rejected, refused, stale = [], [], []
    pins_by_node: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for tkey, nkey, node in table_pins:
        pins_by_node[node].add((tkey, nkey))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def kegg_ids(root):
        return {m for m in members[root] if m.startswith(_KEGG)}

    def pins(root):
        return set().union(*(pins_by_node.get(m, set()) for m in members[root]))

    for m in overrides.merges:
        if m.a not in table or m.b not in table:
            stale.append({"kind": "merge", "a": m.a, "b": m.b})
    for pair in overrides.keep_separate:
        if not all(x in table for x in pair):
            stale.append({"kind": "keep_separate", "pair": sorted(pair)})
    manual = [m for m in overrides.merges if m.a in table and m.b in table]

    ordered = manual + sorted((l for l in links if l.basis == "mnx_same_entity")) + \
        sorted((l for l in links if l.basis != "mnx_same_entity"))
    for link in ordered:
        ra, rb = find(link.a), find(link.b)
        if ra == rb:
            continue
        ga, gb = members[ra], members[rb]
        if any(frozenset((x, y)) in overrides.keep_separate for x in ga for y in gb):
            continue
        ka, kb = kegg_ids(ra), kegg_ids(rb)
        if link.basis not in ("mnx_same_entity", "manual") and ka and kb and ka != kb:
            rejected.append(link)
            continue
        pa, pb = pins(ra), pins(rb)
        clash = [(t, n1, n2) for t, n1 in pa for t2, n2 in pb if t == t2 and n1 != n2]
        if clash:
            refused.append({"link": [link.a, link.b, link.basis],
                            "tables": sorted({c[0] for c in clash})})
            continue
        parent[rb] = ra
        members[ra] |= members.pop(rb)
        # Each member is attributed the basis of the first link that joined it to anything.
        for m in ga | gb:
            attach.setdefault(m, link.basis)

    canonical_of, groups, basis, member_basis = {}, {}, {}, {}
    for root, ms in members.items():
        canon = _choose(ms, reaction_refs, pathway_refs, chebi)
        for m in ms:
            canonical_of[m] = canon
        if len(ms) > 1:
            groups[canon] = sorted(ms)
            for m in ms:
                if m != canon:
                    member_basis[m] = attach[m]
            basis[canon] = sorted({member_basis[m] for m in ms if m != canon})
    return IdentityResult(canonical_of, groups, basis, member_basis, sorted(rejected), refused, stale)


def _union(a, b):
    return sorted(set(a or []) | set(b or []))


def rewrite_to_canonical(table, reactions, subtree_substrates, canonical_of):
    new: dict[str, dict] = {}
    for node in sorted(table, key=lambda n: (canonical_of[n] != n, n)):
        canon = canonical_of[node]
        src = table[node]
        if canon not in new:
            new[canon] = dict(table[canon]) if canon in table else dict(src)
            new[canon]["pathways"] = sorted(set(new[canon].get("pathways") or []))
            new[canon]["evidence_sources"] = sorted(set(new[canon].get("evidence_sources") or []))
        tgt = new[canon]
        tgt["pathways"] = _union(tgt.get("pathways"), src.get("pathways"))
        tgt["evidence_sources"] = _union(tgt.get("evidence_sources"), src.get("evidence_sources"))
        for k, v in src.items():
            if tgt.get(k) is None and v is not None:
                tgt[k] = v
    bare = {n[len(_KEGG):]: c[len(_KEGG):] for n, c in canonical_of.items()
            if n.startswith(_KEGG) and c.startswith(_KEGG)}
    for r in reactions.values():
        r["compounds"] = sorted({bare.get(c, c) for c in r.get("compounds", [])})
    for tc, ids in subtree_substrates.items():
        subtree_substrates[tc] = sorted({canonical_of.get(i, i) for i in ids})
    return new
