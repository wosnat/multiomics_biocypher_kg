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


def stereo_links(table, terms_of, chebi: ChebiIndex):
    by_term: dict[str, list[str]] = defaultdict(list)
    for n, ts in terms_of.items():
        if n in table:
            for t in ts:
                by_term[t].append(n)
    stereo, enan = set(), set()
    for n, ts in terms_of.items():
        if n not in table or not _block(table[n]):
            continue
        for t in ts:
            tt = chebi.term(t) or {}
            for rel, out in (("is_a", stereo), ("enantiomers", enan)):
                for partner in tt.get(rel, []):
                    for o in by_term.get(chebi.primary(partner) or "", []):
                        if o == n or _block(table[o]) != _block(table[n]):
                            continue
                        if is_geometric_pair(n, o, table, terms_of, chebi):
                            continue                       # geometric isomers are never stereo-linked
                        out.add((n, o) if rel == "is_a" else tuple(sorted((n, o))))
    return stereo, enan


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
