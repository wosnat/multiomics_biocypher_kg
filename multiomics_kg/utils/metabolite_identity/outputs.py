"""Stage P: validation, atomic multi-file write, resolved tables, mapping v2, report rows."""
from __future__ import annotations

import io
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable

import pandas as pd

from multiomics_kg.utils.paperconfig_utils import PROJECT_ROOT

REASON_ORDER = ("alias_override", "dominant_stereo_form", "case_fold_convergence",
                "direct_id_canonicalised", "placeholder_recurated", "remeasured_elsewhere", "merged",
                "name_re_resolved")


def validate_build(*, table, roles, subtree_substrates, subtree_roles, class_trees, resolution,
                   chebi_missing_needed, go_substrates=None, go_roles=None):
    errs = []
    errs += [f"ambiguous paper name (needs an override): {a}" for a in resolution.ambiguous]
    errs += [f"same-table collision: {c}" for c in resolution.collisions]
    errs += list(resolution.errors)
    owner: dict[str, str] = {}
    for kind, items in (("Metabolite", table), ("ChemicalRole", roles)):
        for node, e in items.items():
            for x in e.get("all_identifiers", []):
                prev = owner.setdefault(x, f"{kind} {node}")
                if prev != f"{kind} {node}":
                    errs.append(f"identifier {x} on two nodes: {prev} and {kind} {node}")
    iks = Counter(e["inchikey"] for e in table.values() if e.get("inchikey"))
    errs += [f"InChIKey {k} on {n} Metabolite nodes" for k, n in iks.items() if n > 1]
    for node, e in table.items():
        for k in ("stereo_form_of", "enantiomer_of"):
            errs += [f"{node}.{k} → {t} is not a written compound" for t in e.get(k, []) if t not in table]
    for tc, ids in subtree_substrates.items():
        errs += [f"subtree_substrates[{tc}] → {i} is not a written compound" for i in ids if i not in table]
    for tc, ids in subtree_roles.items():
        errs += [f"subtree_roles[{tc}] → {i} is not a ChemicalRole" for i in ids if i not in roles]
    for go_id, targets in sorted((go_substrates or {}).items()):
        errs += [f"go_transport_substrates[{go_id}] → {t} is not a written compound"
                 for t in sorted(targets) if t not in table]
    for go_id, targets in sorted((go_roles or {}).items()):
        errs += [f"go_transport_substrates[{go_id}] → {t} is not a ChemicalRole"
                 for t in sorted(targets) if t not in roles]
    for r, e in roles.items():
        errs += [f"ChemicalRole {r} parent {p} missing" for p in e.get("parents", []) if p not in roles]
    for code, tree in class_trees.items():
        ids = {n["id"] for n in tree["nodes"]}
        errs += [f"{code}: parent {n['parent']} missing" for n in tree["nodes"] if n["parent"] and n["parent"] not in ids]
        errs += [f"{code}: member {m} dangles" for m in tree["members"]
                 if m["metabolite_id"] not in table or m["class_id"] not in ids]
    for key, d in resolution.decisions.items():
        if d.metabolite_id and d.metabolite_id not in table:
            errs.append(f"row {key} resolved to {d.metabolite_id}, which step 6 did not write")
        if d.metabolite_id and not d.stereo_assignment:
            errs.append(f"row {key} resolved without a stereo_assignment")
    errs += [f"ChEBI id {c} needed by a canonical node is missing from the committed ChEBI cache"
             for c in sorted(chebi_missing_needed)]
    return errs


def write_all_atomic(files):
    tmps = []
    try:
        for path, text in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(path.name + ".tmp")
            tmps.append((tmp, path))  # registered before writing so a partial .tmp is cleaned
            tmp.write_text(text)
    except BaseException:
        for tmp, _ in tmps:
            tmp.unlink(missing_ok=True)
        raise
    # A failure inside this loop is not rolled back; all validation happens before.
    for tmp, path in tmps:
        os.replace(tmp, path)


def _repo_relative(path: Path) -> str:
    """``path`` relative to the project root (machine-independent reports); the raw path when
    it lies outside the repository."""
    try:
        return Path(path).resolve().relative_to(Path(PROJECT_ROOT).resolve()).as_posix()
    except ValueError:
        return str(path)


def write_resolved_tables(tables, decisions):
    out = {}
    for t in tables:
        df = pd.read_csv(t.csv_path, dtype=str, keep_default_na=False, skiprows=t.skip_rows or None)
        mid, meth, sa = [""] * len(df), ["unresolved"] * len(df), [""] * len(df)
        for row in t.rows:
            d = decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                mid[row.row_index], meth[row.row_index], sa[row.row_index] = (
                    d.metabolite_id, d.resolution_method, d.stereo_assignment)
        df["metabolite_id"], df["resolution_method"], df["stereo_assignment"] = mid, meth, sa
        resolved_path = t.csv_path.with_name(t.csv_path.stem + "_resolved.csv")
        buf = io.StringIO()
        df.to_csv(buf, index=False)
        out[resolved_path] = buf.getvalue()
        n_res = sum(1 for m in mid if m)
        report = {"source_csv": _repo_relative(t.csv_path), "resolved_csv": _repo_relative(resolved_path),
                  "total_rows": len(df), "resolved": n_res, "unresolved": len(df) - n_res,
                  "match_rate": round(n_res / len(df), 4) if len(df) else 0.0,
                  "method_counts": dict(sorted(Counter(meth).items()))}
        out[t.csv_path.with_name(t.csv_path.stem + "_resolution_report.json")] = json.dumps(report, indent=2, sort_keys=True) + "\n"
    return out


def build_mapping_v2(*, id_aliases, resolution, overrides, tables, compounds_summary, per_paper_stats):
    decisions = defaultdict(lambda: defaultdict(dict))
    for t in tables:
        for row in t.rows:
            d = resolution.decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                key = row.direct_id if row.direct_id and d.resolution_method.endswith("_direct") else row.name_key
                decisions[row.paper][row.table][key] = {"metabolite_id": d.metabolite_id,
                                                        "resolution_method": d.resolution_method,
                                                        "stereo_assignment": d.stereo_assignment}
    unresolved = defaultdict(lambda: defaultdict(list))
    for u in resolution.unresolved:
        unresolved[u["paper"]][u["table"]].append(u.get("name"))
    return {
        "schema_version": 2,
        "id_aliases": dict(sorted(id_aliases.items())),
        "name_lookup": dict(sorted(resolution.name_lookup.items())),
        "paper_overrides": {p: {k: {"id": o.id, "stereo_assignment": o.stereo_assignment, "reason": o.reason}
                                for k, o in sorted(ovs.items())} for p, ovs in sorted(overrides.items())},
        "decisions": json.loads(json.dumps(decisions, sort_keys=True)),
        "unresolved": json.loads(json.dumps(unresolved, sort_keys=True)),
        "ambiguous": {},
        "compounds": compounds_summary,
        "conflicts": {},
        "per_paper_stats": per_paper_stats,
    }


def repoint_reason(*, old_id, new, row, identity, rekeys, dropped, case_fold_keys, written=()):
    """Why a row's metabolite_id changed. ``old_id`` must already be normalised (as the
    caller passes it) so it keys ``dropped`` directly. ``merged`` only when identity links the
    old id to the new node (stage E merge or stage D re-key); otherwise the row's name simply
    resolved somewhere else (``name_re_resolved``) — including when the old target is still
    a written node (``written``, kept for the caller's record)."""
    if not old_id or old_id == new.metabolite_id:
        return None
    if new.resolution_method == "alias_override":
        return "alias_override"
    if new.resolution_method == "dominant_stereo_form":
        return "dominant_stereo_form"
    if row.name_key in case_fold_keys and new.resolution_method == "name_match":
        return "case_fold_convergence"
    if new.resolution_method.endswith("_direct"):
        return "direct_id_canonicalised"
    if old_id.startswith("mnx:"):
        return "placeholder_recurated"
    if old_id in dropped:
        return "remeasured_elsewhere"
    if old_id not in written:
        canon = identity.canonical_of.get(old_id)
        if canon is None and old_id in rekeys:
            canon = identity.canonical_of.get(rekeys[old_id], rekeys[old_id])
        if canon == new.metabolite_id:
            return "merged"
    return "name_re_resolved"


def _remap_fallback(old, norm, *, table, id_index, previous_entries, is_isotopologue, is_class):
    """Bases for an old id none of the direct rules explain (spec §6.10 remap):
    (i) the old id (raw or normalised) is an identifier of a written node;
    (ii) the previous entry's MNX id is an identifier of a written node (the substrate string
         now resolves there; the basis says why the old term no longer stands on its own);
    (iii) otherwise the old node is no longer reachable from any source."""
    for x in (old, norm):
        node = id_index.get(x)
        if node in table:
            # The basis label is approximate: besides genuine ChEBI-KEGG re-keys and primary-id
            # changes, ids relocated as a conjugate acid/base or zwitterion of another node's
            # term also land here (they are that node's identifiers).
            basis = ("rekey_chebi_kegg_xref" if node.startswith("kegg.compound:") and old.startswith("chebi:")
                     else "chebi_primary_id")
            return node, basis
    prev = previous_entries.get(old) or previous_entries.get(norm) or {}
    mnxm = prev.get("mnxm_id")
    node = id_index.get(f"mnx:{mnxm}") if mnxm else None
    if node in table:
        term = old if old.startswith("chebi:") else (f"chebi:{prev['chebi_id']}" if prev.get("chebi_id") else None)
        if is_isotopologue(prev.get("name")):
            return node, "isotopologue_rejected"
        if term and is_class(term):
            return node, "class_term_rejected"
        return node, "mnx_same_entity"
    return "", "no_longer_reachable"


def _substrate_targets(old_string_to_primary, string_to_primary, id_index, table):
    """Old substrate primary → the written node its TCDB substrate string(s) now resolve to
    (through ``id_index``); the node most of its strings reach wins, ties by id."""
    votes: dict[str, Counter] = defaultdict(Counter)
    for s, old in (old_string_to_primary or {}).items():
        new = (string_to_primary or {}).get(s)
        if new is None:
            continue
        node = id_index.get(new, new)
        if node in table:
            votes[old][node] += 1
    return {old: sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] for old, c in votes.items()}


def _substrate_basis(old, node, *, table, chebi, previous_entries, is_isotopologue):
    """Why an old substrate id now lives on ``node``, read from its ChEBI term: an isotope
    variant was rejected, a class term was rejected, it is a conjugate acid/base of the node's
    term, else the node carries the same entity under another (primary) ChEBI id."""
    prev = previous_entries.get(old) or {}
    term = old if old.startswith("chebi:") else (f"chebi:{prev['chebi_id']}" if prev.get("chebi_id") else None)
    p = chebi.primary(term) if (chebi is not None and term) else None
    t = (chebi.term(p) or {}) if p else {}
    if is_isotopologue(prev.get("name")) or is_isotopologue(t.get("name")) or (p and chebi.is_isotope_variant(p)):
        return "isotopologue_rejected"
    if p and chebi.is_class(p):
        return "class_term_rejected"
    e = table[node]
    node_terms = {f"CHEBI:{c}" for c in [e.get("chebi_id"), *(e.get("chebi_ids") or [])] if c}
    if p and (set(t.get("conjugates", [])) & node_terms
              or any(p in (chebi.term(n) or {}).get("conjugates", []) for n in sorted(node_terms))):
        return "chebi_conjugate"
    return "chebi_primary_id"


def remap_rows(*, old_ids, table, identity, rekeys, role_of, dropped, normalize,
               id_index=None, previous_entries=None, is_isotopologue=None, is_class=None,
               old_string_to_primary=None, string_to_primary=None, chebi=None):
    """``id_index`` given → the fallback chain (``_remap_fallback``) accounts for every old id;
    without it, ids the direct rules miss are returned as ``unaccounted``. An old id the
    fallback finds unreachable but whose TCDB substrate string (``old_string_to_primary``, the
    previous build's string → primary map) now resolves to a written node (``string_to_primary``
    taken through ``id_index``) remaps there, with a basis from its ChEBI term
    (``_substrate_basis``)."""
    rows, unaccounted = [], []
    by_string = _substrate_targets(old_string_to_primary, string_to_primary, id_index or {}, table)
    # An old node id that names a more specific stereo form than any KG metabolite (narrowed D):
    # the node it was taken off is the stereo-unspecified form, recorded as such, not as a merge.
    removed_on = {x: n for n, e in table.items() for x in (e.get("removed_stereo_xrefs") or [])}
    for old in sorted(old_ids):
        if old in table:
            continue
        norm = normalize(old)
        if old in role_of:
            rows.append({"old_id": old, "new_id": "", "basis": "removed_role", "role_id": role_of[old]})
        elif norm in dropped:
            basis = "removed_placeholder" if old.startswith("mnx:") else "remeasured_elsewhere"
            rows.append({"old_id": old, "new_id": dropped[norm], "basis": basis, "role_id": ""})
        elif norm in rekeys:
            new = identity.canonical_of.get(rekeys[norm], rekeys[norm])
            rows.append({"old_id": old, "new_id": new, "basis": "rekey_chebi_kegg_xref", "role_id": ""})
        elif norm in identity.member_basis:
            rows.append({"old_id": old, "new_id": identity.canonical_of[norm],
                         "basis": identity.member_basis[norm], "role_id": ""})
        elif old in removed_on or norm in removed_on:
            rows.append({"old_id": old, "new_id": removed_on.get(old) or removed_on[norm],
                         "basis": "stereo_form_removed", "role_id": ""})
        elif norm in table and old.startswith("chebi:"):
            rows.append({"old_id": old, "new_id": norm, "basis": "chebi_primary_id", "role_id": ""})
        elif id_index is not None:
            new, basis = _remap_fallback(old, norm, table=table, id_index=id_index,
                                         previous_entries=previous_entries or {},
                                         is_isotopologue=is_isotopologue or (lambda _n: False),
                                         is_class=is_class or (lambda _c: False))
            node = by_string.get(old) or by_string.get(norm)
            if basis == "no_longer_reachable" and node:
                new, basis = node, _substrate_basis(
                    old, node, table=table, chebi=chebi, previous_entries=previous_entries or {},
                    is_isotopologue=is_isotopologue or (lambda _n: False))
            rows.append({"old_id": old, "new_id": new, "basis": basis, "role_id": ""})
        else:
            unaccounted.append(old)
    return rows, unaccounted


def build_identity_report(**sections):
    return json.loads(json.dumps(sections, sort_keys=True, default=sorted))
