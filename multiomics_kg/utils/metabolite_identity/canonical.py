"""Stages E-F: merge links, union-find with the KEGG guard, winner rule, id rewrite (spec §3.1)."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.stereo import stereo_distinct_kegg

MERGE_BASES = ("mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual")
_KEGG = "kegg.compound:"
# The canonical KEGG id needs this many times the runner-up's reaction references (spec §3.1);
# the same margin as the dominant-stereo-form rule. A closer call needs a `canonical:` pin.
CANONICAL_MARGIN = 2


class IdentityCollisionError(ValueError):
    """An identifier or InChIKey would sit on two Metabolite nodes (stages G/H)."""


@dataclass(frozen=True, order=True)
class MergeLink:
    a: str
    b: str
    basis: str


def _link(x: str, y: str, basis: str) -> MergeLink:
    a, b = sorted((x, y))
    return MergeLink(a, b, basis)


@dataclass(frozen=True)
class CanonicalPin:
    """A ``canonical:`` override: ``id`` is the canonical id of the merge group it ends up in.
    ``group`` is the group as the curator saw it (a pin whose group no longer holds together
    is reported stale, not applied)."""
    id: str
    group: tuple[str, ...]


@dataclass
class IdentityOverrides:
    merges: list[MergeLink]
    keep_separate: set[frozenset[str]]
    canonical: list[CanonicalPin] = field(default_factory=list)


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
    pins = []
    for e in doc.get("canonical") or []:
        group = tuple(sorted(e.get("group") or []))
        if e["id"] not in group:
            raise ValueError(f"metabolite_identity_overrides.yaml: canonical pin {e['id']} is not "
                             f"in its own group {list(group)}")
        pins.append(CanonicalPin(e["id"], group))
    return IdentityOverrides(merges, keep, pins)


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
    # MNX links refused because the two sides are stereo forms of one structure.
    stereo_refused: list[dict] = field(default_factory=list)
    # Groups the margin rule could not decide (stage P fails on any).
    canonical_unresolved: list[dict] = field(default_factory=list)
    canonical_pinned: list[str] = field(default_factory=list)


def _margin_winner(kegg, refs):
    """The KEGG id with >= CANONICAL_MARGIN times the runner-up's references, else None."""
    ranked = sorted(kegg, key=lambda m: (-refs.get(m, 0), m))
    top, second = refs.get(ranked[0], 0), refs.get(ranked[1], 0)
    return ranked[0] if top > 0 and top >= CANONICAL_MARGIN * second else None


def _choose(members, reaction_refs, all_reaction_refs, chebi, pin=None):
    """``(canonical id, None)``, or ``(None, problem)`` when a KEGG group needs a pin (spec §3.1):
    a pin, else most gene-reachable reaction references with a >= 2x margin, else (only when
    that measure is 0 vs 0) most KEGG reaction references with the same margin."""
    kegg = [m for m in members if m.startswith(_KEGG)]
    if pin:
        return pin, None
    if len(kegg) == 1:
        return kegg[0], None
    if kegg:
        win = _margin_winner(kegg, reaction_refs)
        if win is None and not any(reaction_refs.get(m, 0) for m in kegg):
            win = _margin_winner(kegg, all_reaction_refs)
        if win:
            return win, None
        return None, {"members": sorted(kegg),
                      "gene_reachable_reactions": {m: reaction_refs.get(m, 0) for m in sorted(kegg)},
                      "all_kegg_reactions": {m: all_reaction_refs.get(m, 0) for m in sorted(kegg)}}
    ch = [m for m in members if m.startswith("chebi:")]
    if ch:
        return sorted(ch, key=lambda m: (chebi.charge(m) != 0, int(m.split(":")[1])))[0], None
    return sorted(members)[0], None


def _stereo_conflict(ka, kb, chebi):
    for x in sorted(ka):
        for y in sorted(kb):
            why = stereo_distinct_kegg(x[len(_KEGG):], y[len(_KEGG):], chebi)
            if why:
                return x, y, why
    return None


def canonicalise(table, links, *, overrides, reaction_refs, chebi, table_pins, all_reaction_refs=None):
    all_reaction_refs = all_reaction_refs or {}
    parent = {n: n for n in table}
    members = {n: {n} for n in table}
    attach: dict[str, str] = {}
    rejected, refused, stale, stereo_refused = [], [], [], []
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
        # MNX sometimes calls a stereo-unspecified compound, or a sugar and its anomer, one
        # entity; §1 keeps stereo forms apart (linked by stage I, never merged).
        if link.basis == "mnx_same_entity" and ka and kb and ka != kb:
            conflict = _stereo_conflict(ka, kb, chebi)
            if conflict:
                rejected.append(link)
                stereo_refused.append({"link": [link.a, link.b], "kegg": [conflict[0], conflict[1]],
                                       "chebi_names": [list(p) for p in conflict[2]]})
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

    pin_of: dict[str, str] = {}
    for pin in overrides.canonical:
        if pin.id not in table:
            stale.append({"kind": "canonical", "id": pin.id, "group": list(pin.group)})
            continue
        root = find(pin.id)
        if {find(m) for m in pin.group if m in table} != {root} or len(members[root]) < 2:
            stale.append({"kind": "canonical", "id": pin.id, "group": list(pin.group)})
            continue
        if root in pin_of:
            raise ValueError(f"metabolite_identity_overrides.yaml: two canonical pins for one group: "
                             f"{pin_of[root]} and {pin.id}")
        pin_of[root] = pin.id

    canonical_of, groups, basis, member_basis = {}, {}, {}, {}
    unresolved, pinned = [], []
    for root, ms in members.items():
        canon, problem = _choose(ms, reaction_refs, all_reaction_refs, chebi, pin_of.get(root))
        if problem:
            unresolved.append(problem)
            canon = problem["members"][0]       # placeholder; stage P refuses to write
        if root in pin_of:
            pinned.append(canon)
        for m in ms:
            canonical_of[m] = canon
        if len(ms) > 1:
            groups[canon] = sorted(ms)
            for m in ms:
                if m != canon:
                    member_basis[m] = attach[m]
            basis[canon] = sorted({member_basis[m] for m in ms if m != canon})
    return IdentityResult(canonical_of, groups, basis, member_basis, sorted(rejected), refused, stale,
                          stereo_refused=sorted(stereo_refused, key=lambda r: r["link"]),
                          canonical_unresolved=sorted(unresolved, key=lambda u: u["members"]),
                          canonical_pinned=sorted(pinned))


# MNX-derived fields: when MNX says two nodes are one entity but they stay apart (stereo guard,
# keep_separate), they belong to one node only, else stage H puts one mnx: / inchikey: / hmdb:
# identifier on two nodes (chebi_id too: for a KEGG entry it is MNX's ChEBI xref). formula / mass
# stay (identical across stereo forms).
_MNX_FIELDS = ("mnxm_id", "smiles", "hmdb_id", "chebi_id", "mnx_chebi_aliases")


def _structure_terms(nodes, table, chebi):
    terms = set()
    for n in nodes:
        terms |= entity_terms(n, table[n], chebi)
    return [t for t in sorted(terms) if not chebi.is_class(t) and (chebi.term(t) or {}).get("inchikey")]


def split_shared_mnx(table, canonical_of, chebi, reaction_refs):
    """Which pre-canonical nodes lose their MNX fields: ``{node: replacement inchikey | None}``,
    plus a report row per split MNX id (explorer ask 2026-10-15, answer 1).

    The MNX fields stay with the canonical group whose own ChEBI structure matches MNX's InChIKey
    (protonation ignored); with no single match, with the group holding more gene-reachable
    reactions (then the lowest id). Every other group drops ``_MNX_FIELDS`` and takes its InChIKey
    from its own ChEBI term (one sharing MNX's connectivity block first, then the lowest id)."""
    by_mnx: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for n, e in table.items():
        if e.get("mnxm_id"):
            by_mnx[e["mnxm_id"]][canonical_of.get(n, n)].append(n)
    strip, report = {}, []
    for mnx, groups in sorted(by_mnx.items()):
        if len(groups) < 2:
            continue
        ik = next((table[n]["inchikey"] for ns in groups.values() for n in sorted(ns)
                   if table[n].get("inchikey")), None)
        terms = {c: _structure_terms(ns, table, chebi) for c, ns in groups.items()}
        match = [c for c in groups if ik and any(
            chebi.term(t)["inchikey"][:23] == ik[:23] for t in terms[c])]
        pool = match if len(match) == 1 else (match or list(groups))
        keeper = sorted(pool, key=lambda c: (-reaction_refs.get(c, 0), c))[0]
        for c, ns in sorted(groups.items()):
            if c == keeper:
                continue
            own = sorted(terms[c], key=lambda t: (not (ik and chebi.term(t)["inchikey"][:14] == ik[:14]),
                                                  int(t.split(":")[1])))
            new_ik = chebi.term(own[0])["inchikey"] if own else None
            for n in ns:
                strip[n] = new_ik
        report.append({"mnxm_id": mnx, "kept_on": keeper,
                       "stripped_from": sorted(c for c in groups if c != keeper),
                       "basis": "structure_match" if len(match) == 1 else "reaction_refs"})
    return strip, report


def apply_mnx_split(entries, strip):
    for n, new_ik in strip.items():
        e = entries.get(n)
        if e is None:
            continue
        for k in _MNX_FIELDS:
            e[k] = [] if k == "mnx_chebi_aliases" else None
        e["inchikey"] = new_ik


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
