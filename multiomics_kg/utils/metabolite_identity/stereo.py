"""Stages I and L: stereo families (linked, not merged) and the anomer guard (spec §3.5, §3.6)."""
from __future__ import annotations

from collections import defaultdict

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.names import strip_anomer_tokens, strip_geometric_tokens


def _block(entry: dict) -> str | None:
    ik = entry.get("inchikey")
    return ik.split("-")[0] if ik else None


def _names(n, table, terms_of, chebi) -> set[str]:
    return {x for x in {table[n].get("name") or ""} |
            {(chebi.term(t) or {}).get("name") or "" for t in terms_of.get(n, ())} if x}


def is_geometric_pair(a, b, table, terms_of, chebi) -> bool:
    """Names (KEGG or ChEBI) that differ only by cis/trans / E/Z descriptors (spec §3.6)."""
    na, nb = _names(a, table, terms_of, chebi), _names(b, table, terms_of, chebi)
    return any(x.casefold() != y.casefold() and strip_geometric_tokens(x) == strip_geometric_tokens(y)
               for x in na for y in nb)


def _kegg_structures(kegg_id: str, chebi: ChebiIndex) -> list[tuple[str, str, str]]:
    """(connectivity block, stereo layer, name) of each non-class ChEBI term KEGG ``kegg_id``
    cross-references. The node's own ``inchikey`` can't be used: it comes from MNX, so every
    member of an MNX group shares it."""
    out = []
    for t in chebi.kegg_terms(kegg_id):
        tt = chebi.term(t) or {}
        ik = tt.get("inchikey")
        if ik and not chebi.is_class(t):
            out.append((ik[:14], ik[15:23], tt.get("name") or ""))
    return out


def stereo_distinct_kegg(a: str, b: str, chebi: ChebiIndex) -> list[tuple[str, str]]:
    """Why two KEGG ids (bare C numbers) are different stereo forms of one structure, as
    ``[(name in one, name in the other), …]``; empty when they aren't (explorer ask 2026-10-15).

    A term of one id shares a connectivity block with the other id's terms but carries a stereo
    layer none of them has (stereo-unspecified vs specific, D vs L, generic sugar vs anomer).
    Names that differ only by cis/trans / E/Z are exempt: geometric isomers are not stereo forms
    (spec §3.6), and KEGG lists generic + trans caffeate as two ids of one compound."""
    sa, sb = _kegg_structures(a, chebi), _kegg_structures(b, chebi)
    out = []
    for xs, ys in ((sa, sb), (sb, sa)):
        for blk, layer, name in xs:
            same_block = [y for y in ys if y[0] == blk]
            if not same_block or layer in {y[1] for y in same_block}:
                continue
            if all(name.casefold() != y[2].casefold()
                   and strip_geometric_tokens(name) == strip_geometric_tokens(y[2]) for y in same_block):
                continue
            out.extend((name, y[2]) for y in same_block)
    return sorted(set(out))


# A stereo-form path climbs ChEBI is_a at most this many steps (ring form > D-/L- form > generic).
STEREO_PATH_MAX = 3


def _stereo_parents(t: str, chebi: ChebiIndex):
    """ChEBI terms reachable from ``t`` by <= STEREO_PATH_MAX ``is_a`` steps without changing the
    molecular formula or the charge (explorer ask 8, 2026-10-16). ChEBI never links a specific
    form to its stereo-unspecified parent directly (D-glucopyranose > D-glucose > glucose), and
    gives that parent no structure, only a formula; formula + charge is the chemical test that
    survives that. A term ChEBI gives no formula may be passed through (L-rhamnopyranose >
    L-rhamnose > rhamnose) but is never a result; a term with a different formula or charge ends
    the path. ``t`` itself needs a formula."""
    tt = chebi.term(t) or {}
    formula, charge = tt.get("generalized_formula"), tt.get("charge")
    if not formula:
        return set()
    out, frontier, seen = set(), [t], {t}
    for _ in range(STEREO_PATH_MAX):
        nxt = []
        for x in frontier:
            for p in (chebi.term(x) or {}).get("is_a", []):
                pt = chebi.term(p) or {}
                pf = pt.get("generalized_formula")
                if p in seen or (pf is not None and (pf != formula or pt.get("charge") != charge)):
                    continue
                seen.add(p)
                nxt.append(p)
                if pf is not None:
                    out.add(p)
        frontier = nxt
    return out


def _transitive_reduction(edges):
    """Drop (a, b) when b is still reachable from a without it (beta-D-glucose -> D-glucose ->
    glucose keeps two edges, not three)."""
    out: dict[str, set[str]] = defaultdict(set)
    for a, b in edges:
        out[a].add(b)

    def reachable(a, b):
        stack, seen = [m for m in out[a] if m != b], set()
        while stack:
            x = stack.pop()
            if x == b:
                return True
            if x not in seen:
                seen.add(x)
                stack.extend(out[x])
        return False
    return {(a, b) for a, b in edges if not reachable(a, b)}


def stereo_links(table, terms_of, chebi: ChebiIndex, *, curated=(), blocked=frozenset()):
    """Stage I: ``(stereo_pairs, enantiomer_pairs, report)``.

    ``terms_of``: node -> the ChEBI terms the node *publishes* (``chebi_ids``, after the owner
    rule), so a term ChEBI cross-references sloppily to a second KEGG id doesn't link it.

    ``stereo_form_of`` (specific -> generic): some published term of the specific node reaches a
    published term of the generic node through ``_stereo_parents``; the generic's term is not a
    ChEBI class (E: class membership is not a stereo relation); both nodes are
    ``structure_kind = defined``; not a geometric pair (spec §3.6); neither node in ``blocked``;
    plus ``curated`` (specific, generic) pairs for generics ChEBI gives no formula. Pairs that
    point both ways are dropped (reported); the result is transitively reduced.

    ``enantiomer_of``: ChEBI ``is_enantiomer_of`` between published terms of two nodes whose
    InChIKeys share a connectivity block (enantiomers always have structures)."""
    by_term: dict[str, list[str]] = defaultdict(list)
    for n, ts in terms_of.items():
        if n in table:
            for t in ts:
                by_term[t].append(n)
    defined = {n for n, e in table.items() if e.get("structure_kind", "defined") == "defined"}
    stereo, enan = set(), set()
    for n, ts in terms_of.items():
        if n not in table or n in blocked:
            continue
        for t in ts:
            if n in defined:
                for p in _stereo_parents(t, chebi):
                    if chebi.is_class(p):
                        continue
                    for o in by_term.get(p, []):
                        if o == n or o in blocked or o not in defined:
                            continue
                        if is_geometric_pair(n, o, table, terms_of, chebi):
                            continue                       # geometric isomers are never stereo-linked
                        stereo.add((n, o))
            if not _block(table[n]):
                continue
            for partner in (chebi.term(t) or {}).get("enantiomers", []):
                for o in by_term.get(chebi.primary(partner) or "", []):
                    if o != n and _block(table[o]) == _block(table[n]) \
                            and not is_geometric_pair(n, o, table, terms_of, chebi):
                        enan.add(tuple(sorted((n, o))))
    stale = [list(c) for c in curated if c[0] not in table or c[1] not in table]
    cur = {tuple(c) for c in curated if c[0] in table and c[1] in table}
    stereo = (stereo - {(b, a) for a, b in cur}) | cur          # a curated direction wins
    mutual = sorted({tuple(sorted(p)) for p in stereo if (p[1], p[0]) in stereo})
    stereo = {p for p in stereo if tuple(sorted(p)) not in mutual}
    stereo -= {p for p in stereo if tuple(sorted(p)) in enan}     # enantiomers are not stereo forms
    stereo = _transitive_reduction(stereo)
    return stereo, enan, {"mutual_dropped": [list(m) for m in mutual], "stale_curated": stale}


def stereo_families(nodes, stereo_pairs, enantiomer_pairs):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in list(stereo_pairs) + list(enantiomer_pairs):
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    comps: dict[str, set[str]] = defaultdict(set)
    for n in nodes:
        comps[find(n)].add(n)
    return {n: frozenset(comps[find(n)]) for n in nodes}


def unspecified_members(stereo_pairs):
    return {b for _, b in stereo_pairs}


def is_anomer_pair(a, b, table, terms_of, chebi):
    na, nb = _names(a, table, terms_of, chebi), _names(b, table, terms_of, chebi)
    same = any(x.casefold() != y.casefold() and strip_anomer_tokens(x) == strip_anomer_tokens(y)
               for x in na for y in nb)
    sugar = any(chebi.is_carbohydrate(t) for t in terms_of.get(a, set()) | terms_of.get(b, set()))
    return same and sugar


def finalize_stereo(table, stereo_pairs, enantiomer_pairs, survivors):
    for e in table.values():
        e.pop("stereo_form_of", None)
        e.pop("enantiomer_of", None)
    n_s = n_e = 0
    for a, b in sorted(stereo_pairs):
        if a in survivors and b in survivors:
            table[a].setdefault("stereo_form_of", []).append(b)
            n_s += 1
    for a, b in sorted(enantiomer_pairs):
        if a in survivors and b in survivors:
            table[a].setdefault("enantiomer_of", []).append(b)
            n_e += 1
    for e in table.values():
        for k in ("stereo_form_of", "enantiomer_of"):
            if k in e:
                e[k] = sorted(set(e[k]))
    return {"stereo_form_of": n_s, "enantiomer_of": n_e}
