"""Stages E-M, O and the stage-P validation of prepare-data step 6 as one pure function.

Spec 2026-10-13 §5.2. ``run_identity_stages`` takes everything stages A-D produced (in memory)
and returns every file's content in memory; ``render_outputs`` turns that into
``{path: text}`` for ``outputs.write_all_atomic``. Nothing here writes to disk, so a failed run
(``Step6ValidationError``) writes nothing. The only disk read is the paper CSVs, inside
``write_resolved_tables``.

Determinism: every list is sorted, every JSON is ``indent=2, sort_keys=True``; ``today`` reaches
only newly frozen lock entries and a newly seeded lumping pin.
"""
from __future__ import annotations

import copy
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from multiomics_kg.utils.chebi_utils import ChebiIndex, build_chebi_cache
from multiomics_kg.utils.metabolite_identity.canonical import (
    IdentityCollisionError, IdentityOverrides, canonicalise, collect_merge_links, entity_terms, rewrite_to_canonical)
from multiomics_kg.utils.metabolite_identity.classes import (
    COMPOUND_BRITE_TREES, parse_compound_brite, prune_compound_brite)
from multiomics_kg.utils.metabolite_identity.crossrefs import (
    build_identity_fields, is_isotopologue_name, pick_chebi_id, shared_chebi_terms)
from multiomics_kg.utils.metabolite_identity.lumping import update_lumping_pin
from multiomics_kg.utils.metabolite_identity.outputs import (
    build_identity_report, build_mapping_v2, remap_rows, repoint_reason, validate_build,
    write_resolved_tables)
from multiomics_kg.utils.metabolite_identity.paper import (
    Override, PaperTable, StereoContext, resolve_paper_rows)
from multiomics_kg.utils.metabolite_identity.stereo import (
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)
from multiomics_kg.utils.metabolite_identity.structure import apply_structure

_KEGG = "kegg.compound:"
_TRANSIENT = ("mnx_chebi_aliases",)


class Step6ValidationError(RuntimeError):
    def __init__(self, errors: list[str]):
        super().__init__(f"step 6 refused to write ({len(errors)} problem(s)):\n  " + "\n  ".join(errors))
        self.errors = errors


@dataclass
class Step6Inputs:
    table: dict[str, dict]                      # stage D: full node id → entry (+ transient mnx_chebi_aliases)
    kegg_rest: dict                             # kos, pathways, subcategories, categories, reactions
    subtree_substrates: dict[str, list[str]]
    chemical_roles: dict[str, dict]
    subtree_roles: dict[str, list[str]]
    rekeys: dict[str, str]                      # "chebi:<primary n>" → "kegg.compound:C…"
    role_of_old: dict[str, str]                 # old substrate id → chebi.role: id
    tables: list[PaperTable]
    paper_overrides: dict[str, dict[str, Override]]
    name_candidates: dict[str, set[str]]
    lock: dict
    lumping_pin: dict | None
    identity_overrides: IdentityOverrides
    chebi: ChebiIndex
    chebi_from_raw: bool
    rxn_genes: dict[str, set[str]]
    compound_to_pathways: dict[str, list[str]]
    base_pathways: set[str]                     # catalysis ∪ transport-reachable pathways
    brite_raw: dict[str, dict]
    brite_committed: dict[str, dict]
    tcdb_kept: set[str]
    tcdb_levels: dict[str, int]
    tcdb_static: dict                           # kept_tcdb_ids, pfam_bridge, go_bridge, seed_aliases
    previous_ids: set[str]
    old_resolved: dict[tuple[str, str, int], str]
    today: str
    # Previous kegg_data entries by node id (remap fallback: the old entry's mnx id / name).
    previous_entries: dict[str, dict] = field(default_factory=dict)


@dataclass
class Step6Outputs:
    kegg_data: dict
    tcdb_pruned: dict
    class_trees: dict[str, dict]
    chebi_cache: dict | None
    mapping: dict
    lock: dict
    report: dict
    resolved_files: dict[Path, str]


@dataclass
class Step6Paths:
    kegg_data: Path
    tcdb_pruned: Path
    brite_dir: Path
    chebi_cache: Path
    mapping: Path
    lock: Path
    report: Path


def _catalyst_gene_counts(reactions, rxn_genes):
    genes: dict[str, set[str]] = defaultdict(set)
    for rid, r in reactions.items():
        for c in r.get("compounds", []):
            genes[f"{_KEGG}{c}"] |= rxn_genes.get(rid, set())
    return {n: len(g) for n, g in genes.items()}


def _normalizer(chebi: ChebiIndex):
    """Old / paper ids → current form: a secondary ChEBI id becomes its primary."""
    def norm(x: str) -> str:
        if x.startswith("chebi:"):
            p = chebi.primary(x)
            return f"chebi:{p.split(':')[1]}" if p else x
        return x
    return norm


def _measured_set(table, res, old_resolved, norm, id_index):
    """Stage K (spec §3.6a): the metabolomics tag only on final targets; drop candidates with
    no catalysis / transport / measurement reason to exist. A drop whose old paper rows moved
    goes to the target's ``recurated_from_ids`` (never ``all_identifiers``). Old row ids are
    taken through ``id_index``, so an old id that stage D re-keyed or stage E merged onto the
    dropped node still links the drop."""
    targets = {d.metabolite_id for d in res.decisions.values() if d.metabolite_id}
    missing = sorted(t for t in targets if t not in table)
    if missing:
        raise Step6ValidationError([f"paper row target {t} is not a compound step 6 holds" for t in missing])
    for n, e in table.items():
        ev = set(e.get("evidence_sources") or []) - {"metabolomics"}
        if n in targets:
            ev.add("metabolomics")
        e["evidence_sources"] = sorted(ev)
    moved: dict[str, Counter] = defaultdict(Counter)
    for key, old in old_resolved.items():
        d = res.decisions.get(key)
        if old and d and d.metabolite_id:
            n = norm(old)
            moved[id_index.get(n, n)][d.metabolite_id] += 1
    dropped, unmoved = {}, []
    for n in sorted(n for n, e in table.items() if not e["evidence_sources"]):
        if moved.get(n):
            tgt = sorted(moved[n].items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            dropped[n] = tgt
            table[tgt]["recurated_from_ids"] = sorted(set(table[tgt].get("recurated_from_ids") or []) | {n})
        else:
            unmoved.append(n)
        del table[n]
    return targets, dropped, unmoved


def _reprune_pathways(kegg_rest, table, targets, inp):
    """Keep base (catalysis ∪ transport) pathways plus those of the final paper targets, and
    filter every compound / KO / reaction ``pathways`` list to what is kept."""
    allowed = set(inp.base_pathways)
    for t in targets:
        if t.startswith(_KEGG):
            allowed |= set(inp.compound_to_pathways.get(t[len(_KEGG):], []))
    pathways = {p: v for p, v in kegg_rest["pathways"].items() if p in allowed}
    for entries in (table, kegg_rest.get("kos", {}), kegg_rest.get("reactions", {})):
        for e in entries.values():
            if e.get("pathways"):
                e["pathways"] = sorted(p for p in e["pathways"] if p in pathways)
    subs = {v.get("subcategory") for v in pathways.values() if v.get("subcategory")}
    subcategories = {s: v for s, v in kegg_rest["subcategories"].items() if s in subs}
    cats = {v.get("category") for v in subcategories.values() if v.get("category")}
    categories = {c: v for c, v in kegg_rest["categories"].items() if c in cats}
    return pathways, subcategories, categories


def _class_trees(inp, id_index, survivors):
    """Stage M. A tree whose KEGG download failed falls back to the committed pruned file,
    its members re-mapped through the current id aliases (spec §6.2)."""
    out, errors = {}, []
    for code in sorted(COMPOUND_BRITE_TREES):
        raw = inp.brite_raw.get(code)
        if raw:
            nodes, members = parse_compound_brite(code, raw)
        elif code in inp.brite_committed:
            committed = inp.brite_committed[code]
            nodes = {n["id"]: {"parent": n["parent"], "name": n["name"], "level": n["level"]}
                     for n in committed["nodes"]}
            members = []
            for m in committed["members"]:
                mid = id_index.get(m["metabolite_id"], m["metabolite_id"])
                if mid.startswith(_KEGG):           # prune_compound_brite takes bare C numbers
                    members.append((mid[len(_KEGG):], m["class_id"]))
        else:
            errors.append(f"no KEGG BRITE {code}: download failed and no committed file")
            continue
        out[code] = prune_compound_brite(code, nodes, members, canonical_of=id_index, written=survivors)
    if errors:
        raise Step6ValidationError(errors)
    return out


def _groups_by_basis(identity):
    by: dict[str, set[str]] = defaultdict(set)
    for canon, bases in identity.basis.items():
        for b in bases:
            by[b].add(canon)
    return {b: sorted(v) for b, v in sorted(by.items())}


def _exact_name_pairs(table):
    """Spec §1.2 'exact-name match not covered above': report only, never applied."""
    by_name = defaultdict(list)
    for n, e in table.items():
        if e.get("name"):
            by_name[e["name"].casefold()].append(n)
    return [sorted(v) for _, v in sorted(by_name.items()) if len(v) > 1]


def run_identity_stages(inp: Step6Inputs) -> Step6Outputs:
    chebi = inp.chebi
    norm = _normalizer(chebi)
    table = copy.deepcopy(inp.table)
    old_entries = copy.deepcopy(table)
    kegg_rest = copy.deepcopy(inp.kegg_rest)
    reactions = kegg_rest["reactions"]
    subs = {k: list(v) for k, v in inp.subtree_substrates.items()}

    # E — canonicalise. Rows pinned by an override or a direct id refuse merges that would put
    # two of one table's names on one node (assay collision).
    reaction_refs = Counter(f"{_KEGG}{c}" for r in reactions.values() for c in r.get("compounds", []))
    pathway_refs = {n: len(e.get("pathways") or []) for n, e in table.items()}
    pins = []
    for t in inp.tables:
        for row in t.rows:
            ov = inp.paper_overrides.get(t.paper, {}).get(row.name_key)
            fixed = (ov.id if ov else None) or row.direct_id
            fixed = norm(fixed) if fixed else None
            if fixed and fixed in table:
                pins.append((t.key, row.name_key, fixed))
    identity = canonicalise(table, collect_merge_links(table, chebi), overrides=inp.identity_overrides,
                            reaction_refs=reaction_refs, pathway_refs=pathway_refs, chebi=chebi,
                            table_pins=sorted(set(pins)))
    # F — rewrite every id to canonical (reactions and subtree_substrates in place)
    table = rewrite_to_canonical(table, reactions, subs, identity.canonical_of)
    # ChEBI terms that xref two KEGG nodes the guard kept apart are published on neither.
    shared = shared_chebi_terms(table, identity, chebi)
    try:
        # G — structure (B2 scrub before any identifier is published)
        chebi_of = {}
        for n, e in table.items():
            cid = pick_chebi_id(n, e, chebi, shared)
            chebi_of[n] = f"CHEBI:{cid}" if cid else None
        structure = apply_structure(table, chebi, chebi_of)
        # H — cross-reference fields
        identifier_index = build_identity_fields(table, identity, inp.rekeys, chebi, old_entries,
                                                 shared=shared)
    except IdentityCollisionError as exc:        # identifier / InChIKey on two nodes
        raise Step6ValidationError([str(exc)]) from exc
    id_index = dict(identifier_index)
    for old, canon in identity.canonical_of.items():
        id_index.setdefault(old, canon)
    for old, new in inp.rekeys.items():
        id_index.setdefault(old, identity.canonical_of.get(new, new))
    # I — stereo families (in memory only)
    terms_of = {n: entity_terms(n, e, chebi) for n, e in table.items()}
    stereo_pairs, enan_pairs = stereo_links(table, terms_of, chebi)
    families = stereo_families(set(table), stereo_pairs, enan_pairs)
    # J — resolve every paper row
    ctx = StereoContext(families, unspecified_members(stereo_pairs),
                        lambda a, b: is_anomer_pair(a, b, table, terms_of, chebi),
                        _catalyst_gene_counts(reactions, inp.rxn_genes))
    name_candidates = {k: {id_index.get(norm(c), norm(c)) for c in v}
                       for k, v in sorted(inp.name_candidates.items())}
    res = resolve_paper_rows(inp.tables, overrides=inp.paper_overrides, name_candidates=name_candidates,
                             id_index=id_index, ctx=ctx, lock=inp.lock, today=inp.today)
    # K — measured set; drop measurement-only losers
    targets, dropped, dropped_unmoved = _measured_set(table, res, inp.old_resolved, norm, id_index)
    pathways, subcategories, categories = _reprune_pathways(kegg_rest, table, targets, inp)
    survivors = set(table)
    structure["counts"] = dict(sorted(Counter(e["structure_kind"] for e in table.values()).items()))
    structure["scrubbed"] = [n for n in structure["scrubbed"] if n in survivors]
    structure["filled"] = [n for n in structure["filled"] if n in survivors]
    # L — stereo lists on survivors only
    stereo_counts = finalize_stereo(table, stereo_pairs, enan_pairs, survivors)
    # M — compound classes
    class_trees = _class_trees(inp, id_index, survivors)
    # O — lumping pin, on canonical subtree_substrates
    pin, lumping = update_lumping_pin(inp.lumping_pin, inp.tcdb_kept, inp.tcdb_levels,
                                      {tc: len(ids) for tc, ids in subs.items()}, inp.today)
    # P — validate (writing is the caller's job, only after this returns)
    missing = set()
    if not inp.chebi_from_raw:
        missing = {n for n in table if n.startswith("chebi:") and chebi.primary(n) is None}
    errors = validate_build(table=table, roles=inp.chemical_roles, subtree_substrates=subs,
                            subtree_roles=inp.subtree_roles, class_trees=class_trees, resolution=res,
                            chebi_missing_needed=missing)

    # An old id whose canonical node was dropped at stage K remaps like the dropped node.
    def old_norm(x: str) -> str:
        n = norm(x)
        c = id_index.get(n, n)
        return c if c in dropped else n
    remap, unaccounted = remap_rows(old_ids=inp.previous_ids, table=table, identity=identity,
                                    rekeys=inp.rekeys, role_of=inp.role_of_old, dropped=dropped,
                                    normalize=old_norm, id_index=id_index,
                                    previous_entries=inp.previous_entries,
                                    is_isotopologue=is_isotopologue_name, is_class=chebi.is_class)
    errors += [f"remap {r['old_id']} → {r['new_id']} ({r['basis']}): not a compound step 6 writes"
               for r in remap if r["new_id"] and r["new_id"] not in table]
    if errors:
        raise Step6ValidationError(errors)

    compounds, additional = {}, {}
    for n in sorted(table):
        e = {k: v for k, v in table[n].items() if k not in _TRANSIENT and v not in (None, [], "", {})}
        if n.startswith(_KEGG):
            compounds[n[len(_KEGG):]] = e
        else:
            additional[n] = e
    kegg_data = dict(kegg_rest, pathways=pathways, subcategories=subcategories, categories=categories,
                     compounds=compounds, additional_compounds=additional)
    tcdb_pruned = dict(inp.tcdb_static, subtree_substrates=dict(sorted(subs.items())),
                       chemical_roles=inp.chemical_roles, subtree_roles=inp.subtree_roles, lumping_pin=pin)

    # Repoints (release TSV) and the old-id remap.
    raw_names: dict[str, set[str]] = defaultdict(set)
    for t in inp.tables:
        for row in t.rows:
            raw_names[row.name_key].add(row.name)
    case_fold_keys = {k for k, v in raw_names.items() if len(v) > 1}
    repoints = []
    for t in inp.tables:
        for row in t.rows:
            key = (row.paper, row.table, row.row_index)
            d, old = res.decisions.get(key), inp.old_resolved.get(key)
            if d is None:
                continue
            reason = repoint_reason(old_id=old and old_norm(old), new=d, row=row, identity=identity,
                                    rekeys=inp.rekeys, dropped=dropped, case_fold_keys=case_fold_keys)
            if reason:
                repoints.append({"paper": row.paper, "table": row.table, "paper_name": row.name,
                                 "old_id": old, "new_id": d.metabolite_id,
                                 "stereo_assignment": d.stereo_assignment, "reason": reason})

    # Mapping v2 (a record of the paper decisions; no step reads it).
    summary: dict[str, dict[str, set]] = defaultdict(lambda: {"aliases": set(), "methods": set()})
    per_paper: dict[str, dict[str, int]] = defaultdict(lambda: {"resolved": 0, "unresolved": 0, "total": 0})
    for t in inp.tables:
        for row in t.rows:
            d = res.decisions.get((row.paper, row.table, row.row_index))
            s = per_paper[row.paper]
            s["total"] += 1
            if d and d.metabolite_id:
                summary[d.metabolite_id]["aliases"].add(row.name)
                summary[d.metabolite_id]["methods"].add(d.resolution_method)
                s["resolved"] += 1
            else:
                s["unresolved"] += 1
    compounds_summary = {k: {"aliases": sorted(v["aliases"]), "methods": sorted(v["methods"])}
                         for k, v in sorted(summary.items())}
    id_aliases = {k: v for k, v in sorted(id_index.items()) if k != v and v in table}
    mapping = build_mapping_v2(id_aliases=id_aliases,
                               resolution=res, overrides=inp.paper_overrides, tables=inp.tables,
                               compounds_summary=compounds_summary,
                               per_paper_stats=dict(sorted(per_paper.items())))

    # The pin's frozen_on already records a seeding run; a per-run "seeded" flag would make a
    # re-run with the committed pin differ byte-wise, so it is left out of the report.
    lumping_report = {k: v for k, v in lumping.items() if k != "seeded"}
    report = build_identity_report(
        chebi_release=chebi.release,
        merge_groups={c: {"members": m, "basis": identity.basis.get(c, [])}
                      for c, m in sorted(identity.members.items())},
        merge_groups_by_basis=_groups_by_basis(identity),
        id_aliases=id_aliases,
        rekeys=inp.rekeys, shared_chebi_terms=shared, rejected_links=[[l.a, l.b, l.basis] for l in identity.rejected],
        refused_assay_collisions=identity.refused_collisions,
        stale_identity_override=identity.stale_overrides, stale_alias_override=res.stale_overrides,
        stale_lock_entry=res.stale_lock, lock_target_remapped=res.lock_remapped,
        roles_removed=inp.role_of_old, dropped=dropped, dropped_without_rows=dropped_unmoved,
        dropped_same_entity_review=[{"dropped": d, "target": t,
                                     "dropped_name": old_entries.get(d, {}).get("name"),
                                     "target_name": table[t].get("name")}
                                    for d, t in sorted(dropped.items())],
        heuristic_name_candidates=_exact_name_pairs(table),
        ambiguous_paper_names=res.ambiguous, unresolved_paper_names=res.unresolved,
        dominant_form_drift=res.drift, lumping=lumping_report, structure=structure,
        stereo_counts=stereo_counts, remap=remap, unaccounted_old_ids=unaccounted, repoints=repoints,
        counts={"metabolites": len(table), "kegg": len(compounds), "additional": len(additional),
                "merged_away": sum(len(m) - 1 for m in identity.members.values()),
                "rekeyed": len(inp.rekeys), "roles": len(inp.chemical_roles), "dropped": len(dropped),
                "dropped_without_rows": len(dropped_unmoved),
                "class_nodes": sum(len(t["nodes"]) for t in class_trees.values()),
                "class_members": sum(len(t["members"]) for t in class_trees.values())})
    return Step6Outputs(kegg_data, tcdb_pruned, class_trees,
                        build_chebi_cache(chebi) if inp.chebi_from_raw else None,
                        mapping, res.lock, report, write_resolved_tables(inp.tables, res.decisions))


def _dump(obj) -> str:
    return json.dumps(obj, indent=2, sort_keys=True) + "\n"


def render_outputs(out: Step6Outputs, paths: Step6Paths) -> dict[Path, str]:
    files = {paths.kegg_data: _dump(out.kegg_data), paths.tcdb_pruned: _dump(out.tcdb_pruned),
             paths.mapping: _dump(out.mapping), paths.lock: _dump(out.lock), paths.report: _dump(out.report)}
    for code, tree in sorted(out.class_trees.items()):
        files[paths.brite_dir / f"brite_{code}.json"] = _dump(tree)
    if out.chebi_cache is not None:
        files[paths.chebi_cache] = _dump(out.chebi_cache)
    files.update(out.resolved_files)
    return files
