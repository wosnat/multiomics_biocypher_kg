"""Stage C: role-typed TCDB substrates → ChemicalRole hierarchy (spec §3.4)."""
from __future__ import annotations

from collections import deque

from multiomics_kg.utils.chebi_utils import ROLE_BRANCHES, ROLE_ROOT, ChebiIndex


def role_node_id(primary: str) -> str:
    return f"chebi.role:{primary.split(':', 1)[1]}"


def _primary_of_string(s: str, chebi: ChebiIndex) -> str | None:
    head = s.partition(";")[0]
    prefix, _, local = head.partition(":")
    if prefix.strip().upper() != "CHEBI" or not local.strip():
        return None
    return chebi.primary(local.strip())


def partition_role_strings(strings, chebi):
    molecules, roles = set(), {}
    for s in sorted(strings):
        p = _primary_of_string(s, chebi)
        if p is not None and chebi.is_role(p):
            roles[s] = role_node_id(p)
        else:
            molecules.add(s)
    return molecules, roles


def _level(p: str, chebi: ChebiIndex) -> int:
    q, seen = deque([(p, 0)]), {p}
    while q:
        x, d = q.popleft()
        if x in ROLE_BRANCHES:
            return d
        for y in (chebi.term(x) or {}).get("is_a", []):
            if y not in seen:
                seen.add(y)
                q.append((y, d + 1))
    raise ValueError(f"role {p} has no path to a ChEBI role branch")


def build_role_hierarchy(role_ids, chebi):
    primaries = {f"CHEBI:{r.split(':', 1)[1]}" for r in role_ids}
    keep = set()
    for p in primaries:
        keep.add(p)
        keep |= {a for a in chebi.ancestors(p) if a != ROLE_ROOT and chebi.is_role(a)}
    out = {}
    for p in sorted(keep, key=lambda c: int(c.split(":")[1])):
        rid = role_node_id(p)
        t = chebi.term(p)
        out[rid] = {
            "name": t["name"],
            "parents": sorted(role_node_id(x) for x in t.get("is_a", []) if x in keep),
            "level": _level(p, chebi),
            "level_kind": "chebi_role_branch" if p in ROLE_BRANCHES else "chebi_role",
            "all_identifiers": sorted([rid, f"chebi:{p.split(':')[1]}"]
                                      + [f"chebi:{a.split(':')[1]}" for a in chebi.secondary_ids(p)]),
        }
    return out


def build_subtree_roles(subtree_strings, role_map):
    out = {}
    for tc, strings in sorted(subtree_strings.items()):
        ids = sorted({role_map[s] for s in strings if s in role_map})
        if ids:
            out[tc] = ids
    return out
