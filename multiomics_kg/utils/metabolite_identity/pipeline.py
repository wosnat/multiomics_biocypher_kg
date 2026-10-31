"""Stages E-M and the stage-P validation of prepare-data step 6 as one pure function.

Spec 2026-10-13 §5.2. ``run_identity_stages`` takes everything stages A-D produced (in memory)
and returns every file's content in memory; ``render_outputs`` turns that into
``{path: text}`` for ``outputs.write_all_atomic``. Nothing here writes to disk, so a failed run
(``Step6ValidationError``) writes nothing. The only disk read is the paper CSVs, inside
``write_resolved_tables``.

Determinism: every list is sorted, every JSON is ``indent=2, sort_keys=True``; ``today`` reaches
only newly frozen lock entries.
"""
from __future__ import annotations

import copy
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from multiomics_kg.utils.chebi_utils import ChebiIndex, build_chebi_cache
from multiomics_kg.utils.metabolite_identity.canonical import (
    IdentityCollisionError, IdentityOverrides, apply_mnx_split, canonicalise, collect_merge_links,
    entity_terms, rewrite_to_canonical, split_shared_mnx)
from multiomics_kg.utils.metabolite_identity.classes import (
    COMPOUND_BRITE_TREES, parse_compound_brite, prune_compound_brite)
from multiomics_kg.utils.metabolite_identity.go_transport import rewrite_go_targets, target_distribution
from multiomics_kg.utils.metabolite_identity.crossrefs import (
    _own_family, apply_xref_overrides, build_identity_fields, withhold_contested_terms, chebi_node_structure, is_isotopologue_name, node_olds, pick_chebi_id, shared_chebi_terms)
from multiomics_kg.utils.metabolite_identity.outputs import (
    build_identity_report, build_mapping_v2, remap_rows, repoint_reason, validate_build,
    write_resolved_tables)
from multiomics_kg.utils.metabolite_identity.paper import (
    Override, PaperTable, StereoContext, resolve_paper_rows)
from multiomics_kg.utils.metabolite_identity.stereo import (
    chebi_stereo_forms, finalize_stereo, is_anomer_pair, kegg_stereo_forms, stereo_distinct_kegg,
    stereo_families, stereo_links, unspecified_members)
from multiomics_kg.utils.metabolite_identity.structure import apply_structure

_KEGG = "kegg.compound:"
_TRANSIENT = ("mnx_chebi_aliases",)
_REPORT_BASELINE = ("previous kegg_data.json (the last step-6 run, not a release; the "
                    "release-to-release record is docs/kg-changes/metabolite-identity/"
                    "metabolite_id_remap.tsv)")


class Step6ValidationError(RuntimeError):
    def __init__(self, errors: list[str]):
        super().__init__(f"step 6 refused to write ({len(errors)} problem(s)):\n  " + "\n  ".join(errors))
        self.errors = errors


@dataclass
class Step6Inputs:
    table: dict[str, dict]                      # stage D: full node id → entry (+ transient mnx_chebi_aliases)
    kegg_rest: dict                             # kos, pathways, subcategories, categories, reactions
    system_substrates: dict[str, list[str]]     # {tc_system: [raw-resolved primary ids]}
    chemical_roles: dict[str, dict]
    system_roles: dict[str, list[str]]
    rekeys: dict[str, str]                      # "chebi:<primary n>" → "kegg.compound:C…"
    role_of_old: dict[str, str]                 # old substrate id → chebi.role: id
    tables: list[PaperTable]
    paper_overrides: dict[str, dict[str, Override]]
    name_candidates: dict[str, set[str]]
    lock: dict
    identity_overrides: IdentityOverrides
    chebi: ChebiIndex
    chebi_from_raw: bool
    rxn_genes: dict[str, set[str]]
    compound_to_pathways: dict[str, list[str]]
    base_pathways: set[str]                     # catalysis ∪ transport-reachable pathways
    brite_raw: dict[str, dict]
    brite_committed: dict[str, dict]
    tcdb_static: dict                           # kept_tcdb_ids, pfam_bridge, go_bridge, seed_aliases
    previous_ids: set[str]
    old_resolved: dict[tuple[str, str, int], str]
    today: str
    # Previous kegg_data entries by node id (remap fallback: the old entry's mnx id / name).
    previous_entries: dict[str, dict] = field(default_factory=dict)
    # TCDB substrate string → primary: the previous build's map (re-resolved with chebi=None)
    # and this build's. An old id nothing else explains remaps to where its string now goes.
    old_string_to_primary: dict[str, str] = field(default_factory=dict)
    string_to_primary: dict[str, str] = field(default_factory=dict)
    # KEGG compound → every KEGG reaction listing it (not only gene-reachable ones): the
    # canonical-id fallback measure when no member has a gene-reachable reaction (spec §3.1).
    all_reaction_refs: dict[str, int] = field(default_factory=dict)
    # GO transport targets (spec 2026-10-18 §5.4), from stage D: {go_id: {primary: [go_chebi_ids]}}
    # for metabolite targets (pre-canonical ids) and {go_id: {chebi.role:n: [ids]}} for roles.
    go_substrates: dict = field(default_factory=dict)
    go_roles: dict = field(default_factory=dict)
    # The target filter's ChEBI ids resolved to primaries (plus each entry's own chebi:<n>);
    # matched again after canonicalisation. go_filter_primaries: {"chebi:<n>": primary | None}.
    go_excluded_primaries: set = field(default_factory=set)
    go_filter_primaries: dict = field(default_factory=dict)
    go_release: str | None = None
    # {strain: [[GO terms of one gene that carry a GO transport target], ...]} (distribution)
    go_distribution_input: dict = field(default_factory=dict)
    go_drops: dict = field(default_factory=dict)          # stage-D drop counts
    go_axioms_in_corpus: int = 0                          # KG GO terms with a transport axiom
    go_unresolved_ids: list = field(default_factory=list)  # GO ChEBI ids no node could take
    # TCDB substrate target filter (config/tcdb_substrate_target_filter.yaml): applied at stage D
    # on raw primaries, matched again after canonicalisation (a kept target merged onto a
    # filtered node is dropped). tcdb_filter_primaries: {"chebi:<n>": primary | None}.
    tcdb_excluded_primaries: set = field(default_factory=set)
    tcdb_filter_primaries: dict = field(default_factory=dict)
    tcdb_filter_drops: dict = field(default_factory=dict)  # stage-D drops per raw primary


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


def _measured_set(table, res, old_resolved, norm, id_index, previous_entries=None):
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
    # Carry forward earlier builds' recurated_from_ids (spec §3.6a): a rerun over its own
    # outputs sees no moved rows, so the provenance must come from the previous entry of each
    # surviving node (previous node ids taken through id_index). Ids that are live nodes again
    # are not recurations any more.
    carried: dict[str, set[str]] = defaultdict(set)
    for old_node, prev in sorted((previous_entries or {}).items()):
        ids = prev.get("recurated_from_ids") or []
        if not ids:
            continue
        n = norm(old_node)
        tgt = id_index.get(n, n)
        if tgt in table:
            carried[tgt] |= {i for i in ids if i not in table}
    for tgt, ids in sorted(carried.items()):
        merged = set(table[tgt].get("recurated_from_ids") or []) | ids
        if merged:
            table[tgt]["recurated_from_ids"] = sorted(merged)
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


def _identifier_coverage_errors(table, identity, rekeys, dropped, dropped_unmoved=()):
    """Spec §2 goal 12: every re-keyed id and every merged-away member id resolves through
    exactly one written node's all_identifiers. Ids whose canonical node was dropped at stage
    K go with it: a drop with moved rows is remapped (remeasured_elsewhere); a measurement-only
    candidate no row kept was never a written node."""
    owners: dict[str, list[str]] = defaultdict(list)
    for n, e in table.items():
        for x in e.get("all_identifiers") or []:
            owners[x].append(n)
    wanted = {old: identity.canonical_of.get(new, new) for old, new in rekeys.items()}
    for canon, members in identity.members.items():
        for m in members:
            if m != canon:
                wanted.setdefault(m, canon)
    errs = []
    for x, canon in sorted(wanted.items()):
        if canon in dropped or canon in dropped_unmoved:
            continue
        n = len(owners.get(x, []))
        if n != 1:
            errs.append(f"identifier {x} (re-keyed / merged into {canon}) resolves through {n} "
                        f"written nodes, expected exactly 1")
    return errs


def _canonical_margin_errors(unresolved, table):
    """Stage P (spec §3.1): a merge group whose canonical id the margin rule can't decide.
    The message is meant to be pasted from: members, names, both counts, a ready pin."""
    errs = []
    for u in unresolved:
        rows = "; ".join(f"{m} {table.get(m, {}).get('name')!r}: {u['gene_reachable_reactions'][m]} "
                         f"gene-reachable / {u['all_kegg_reactions'][m]} all KEGG reactions"
                         for m in u["members"])
        group = ", ".join(f'"{m}"' for m in u["members"])
        errs.append(
            f"canonical_margin_unresolved: no member of merge group [{rows}] has >= 2x the "
            f"runner-up's reactions. Pick the canonical id and add to "
            f"config/metabolite_identity_overrides.yaml:\n"
            f"canonical:\n  - id: \"<one of the members>\"\n    group: [{group}]\n"
            f"    reason: \"<why this member>\"")
    return errs


def _stereo_group_errors(identity, chebi):
    """Stage P: no written node holds two KEGG ids that are stereo forms of one structure
    (explorer ask 2026-10-15, ask 4). A curated ``manual`` merge is the one exception."""
    errs = []
    for canon, members in sorted(identity.members.items()):
        if "manual" in identity.basis.get(canon, []):
            continue
        kegg = sorted(m[len(_KEGG):] for m in members if m.startswith(_KEGG))
        for i, a in enumerate(kegg):
            for b in kegg[i + 1:]:
                why = stereo_distinct_kegg(a, b, chebi)
                if why:
                    errs.append(f"merge group {canon} holds stereo forms {a} / {b} "
                                f"({why[0][0]!r} vs {why[0][1]!r}); the stereo guard should have refused it")
    return errs


def _anomer_xref_review(table, stereo_pairs, terms_of, chebi):
    """Explorer ask 5: anomer nodes whose generic sibling is also written, with the HMDB / ChEBI
    xrefs each side carries. A review list; nothing is moved automatically."""
    rows = []
    for a, b in sorted(stereo_pairs):
        if a in table and b in table and is_anomer_pair(a, b, table, terms_of, chebi):
            rows.append({"anomer": a, "anomer_name": table[a].get("name"),
                         "anomer_hmdb_id": table[a].get("hmdb_id"),
                         "anomer_chebi_ids": table[a].get("chebi_ids") or [],
                         "generic": b, "generic_name": table[b].get("name"),
                         "generic_hmdb_id": table[b].get("hmdb_id"),
                         "generic_chebi_ids": table[b].get("chebi_ids") or []})
    return rows


def _exact_name_pairs(table, keep_separate=frozenset()):
    """Spec §1.2 'exact-name match not covered above': report only, never applied. Pairs a
    human already decided to keep separate are not repeated."""
    by_name = defaultdict(list)
    for n, e in table.items():
        if e.get("name"):
            by_name[e["name"].casefold()].append(n)
    out = []
    for _, v in sorted(by_name.items()):
        v = sorted(v)
        if len(v) == 2 and frozenset(v) in keep_separate:
            continue
        if len(v) > 1:
            out.append(v)
    return out


def run_identity_stages(inp: Step6Inputs) -> Step6Outputs:
    chebi = inp.chebi
    norm = _normalizer(chebi)
    table = copy.deepcopy(inp.table)
    old_entries = copy.deepcopy(table)
    kegg_rest = copy.deepcopy(inp.kegg_rest)
    reactions = kegg_rest["reactions"]
    subs = {k: list(v) for k, v in inp.system_substrates.items()}

    # E — canonicalise. Rows pinned by an override or a direct id refuse merges that would put
    # two of one table's names on one node (assay collision).
    reaction_refs = Counter(f"{_KEGG}{c}" for r in reactions.values() for c in r.get("compounds", []))
    all_reaction_refs = {f"{_KEGG}{c}": n for c, n in inp.all_reaction_refs.items()}
    pins = []
    for t in inp.tables:
        for row in t.rows:
            ov = inp.paper_overrides.get(t.paper, {}).get(row.name_key)
            fixed = (ov.id if ov else None) or row.direct_id
            fixed = norm(fixed) if fixed else None
            if fixed and fixed in table:
                pins.append((t.key, row.name_key, fixed))
    identity = canonicalise(table, collect_merge_links(table, chebi), overrides=inp.identity_overrides,
                            reaction_refs=reaction_refs, all_reaction_refs=all_reaction_refs, chebi=chebi,
                            table_pins=sorted(set(pins)))
    stage_e_errors = (_canonical_margin_errors(identity.canonical_unresolved, table)
                      + _stereo_group_errors(identity, chebi))
    if stage_e_errors:
        raise Step6ValidationError(stage_e_errors)
    # MNX fields of nodes MNX joins but stage E kept apart go to one node (old_entries too:
    # stage H reads a node's own pre-merge entry back).
    mnx_strip, mnx_split = split_shared_mnx(table, identity.canonical_of, chebi, reaction_refs)
    apply_mnx_split(table, mnx_strip)
    apply_mnx_split(old_entries, mnx_strip)
    # F — rewrite every id to canonical (reactions and system_substrates in place)
    table = rewrite_to_canonical(table, reactions, subs, identity.canonical_of)
    # GO transport targets follow the same rewrite; the target filter is matched again on
    # canonical ids (a kept target merged onto a filtered node is dropped). Roles pass through.
    tcdb_excluded = {identity.canonical_of.get(p, p) for p in inp.tcdb_excluded_primaries}
    tcdb_filtered_after = 0
    for tc in list(subs):
        kept_subs = [p for p in subs[tc] if p not in tcdb_excluded]
        tcdb_filtered_after += len(subs[tc]) - len(kept_subs)
        if kept_subs:
            subs[tc] = kept_subs
        else:
            del subs[tc]
    go_excluded = {identity.canonical_of.get(p, p) for p in inp.go_excluded_primaries}
    go_subs, go_filtered_after = rewrite_go_targets(inp.go_substrates, identity.canonical_of, go_excluded)
    go_roles = {g: {r: sorted(ids) for r, ids in m.items()} for g, m in inp.go_roles.items() if m}
    # ChEBI terms that xref two KEGG nodes the guard kept apart are published on neither.
    shared = shared_chebi_terms(table, identity, chebi)
    shared = withhold_contested_terms(table, identity, inp.rekeys, chebi, old_entries, shared)
    try:
        # G — structure (B2 scrub before any identifier is published)
        # The same pick as stage H (members and re-keyed ids included), so a node whose only
        # ChEBI term is a class is classified generic here and gets that class as chebi_id there.
        chebi_of = {}
        for n, e in table.items():
            cid = pick_chebi_id(n, e, chebi, shared, olds=node_olds(n, identity, inp.rekeys),
                                old_entries=old_entries)
            chebi_of[n] = f"CHEBI:{cid}" if cid else None
        structure = apply_structure(table, chebi, chebi_of)
        # H — cross-reference fields
        identifier_index = build_identity_fields(table, identity, inp.rekeys, chebi, old_entries,
                                                 shared=shared)
        xref_applied, xref_stale = apply_xref_overrides(
            table, identifier_index, inp.identity_overrides.xrefs if inp.identity_overrides else [], chebi)
        # narrowed D: a chebi: node's structure / ChEBI ids are its own term's (after stage G)
        chebi_structure = chebi_node_structure(table, identifier_index, chebi)
    except IdentityCollisionError as exc:        # identifier / InChIKey on two nodes
        raise Step6ValidationError([str(exc)]) from exc
    id_index = dict(identifier_index)
    for old, canon in identity.canonical_of.items():
        id_index.setdefault(old, canon)
    for old, new in inp.rekeys.items():
        id_index.setdefault(old, identity.canonical_of.get(new, new))
    # I — stereo families (in memory only)
    terms_of = {n: entity_terms(n, e, chebi) for n, e in table.items()}
    # Stereo links read each node's own ChEBI identity: the KEGG xrefs (or own term) of the node and
    # its merged members, minus terms the owner rule withholds from it, plus their ChEBI conjugate
    # acids / bases. Not MNX's ChEBI aliases,
    # which drag in classes ("N-acetyl-D-hexosamine" on GlcNAc, "cresol" on benzyl alcohol).
    published = {}
    for n in table:
        own = set()
        for m in node_olds(n, identity, inp.rekeys):
            own |= entity_terms(m, old_entries.get(m, {}), chebi)
        own = {t for t in own if n not in (shared.get(t) or ())}
        # conjugate acids / bases are the same compound (explorer ask 8, condition 1)
        published[n] = set().union(*(_own_family(t, chebi) for t in own)) if own else set()
    ov = inp.identity_overrides
    stereo_pairs, enan_pairs, stereo_report = stereo_links(
        table, published, chebi, curated=ov.stereo_links if ov else (),
        blocked=ov.no_stereo_link if ov else frozenset())
    families = stereo_families(set(table), stereo_pairs, enan_pairs)
    # J — resolve every paper row
    ctx = StereoContext(families, unspecified_members(stereo_pairs),
                        lambda a, b: is_anomer_pair(a, b, table, terms_of, chebi),
                        _catalyst_gene_counts(reactions, inp.rxn_genes),
                        chebi_stereo=chebi_stereo_forms(table, published, chebi)
                        | kegg_stereo_forms(table, chebi))
    name_candidates = {k: {id_index.get(norm(c), norm(c)) for c in v}
                       for k, v in sorted(inp.name_candidates.items())}
    res = resolve_paper_rows(inp.tables, overrides=inp.paper_overrides, name_candidates=name_candidates,
                             id_index=id_index, ctx=ctx, lock=inp.lock, today=inp.today)
    # K — measured set; drop measurement-only losers
    with chebi.untracked():          # old row ids are baseline-only lookups
        targets, dropped, dropped_unmoved = _measured_set(table, res, inp.old_resolved, norm, id_index,
                                                     inp.previous_entries)
    pathways, subcategories, categories = _reprune_pathways(kegg_rest, table, targets, inp)
    survivors = set(table)
    structure["counts"] = dict(sorted(Counter(e["structure_kind"] for e in table.values()).items()))
    structure["scrubbed"] = [n for n in structure["scrubbed"] if n in survivors]
    structure["filled"] = [n for n in structure["filled"] if n in survivors]
    anomer_xrefs = _anomer_xref_review(table, stereo_pairs, terms_of, chebi)
    # L — stereo lists on survivors only
    stereo_counts = finalize_stereo(table, stereo_pairs, enan_pairs, survivors)
    # M — compound classes
    class_trees = _class_trees(inp, id_index, survivors)
    # P — validate (writing is the caller's job, only after this returns)
    missing = set()
    if not inp.chebi_from_raw:
        missing = {n for n in table if n.startswith("chebi:") and chebi.primary(n) is None}
    errors = validate_build(table=table, roles=inp.chemical_roles, system_substrates=subs,
                            system_roles=inp.system_roles, class_trees=class_trees, resolution=res,
                            chebi_missing_needed=missing, go_substrates=go_subs, go_roles=go_roles)

    # An old id whose canonical node was dropped at stage K remaps like the dropped node.
    def old_norm(x: str) -> str:
        n = norm(x)
        c = id_index.get(n, n)
        return c if c in dropped else n
    with chebi.untracked():          # previous ids are baseline-only lookups
        remap, unaccounted = remap_rows(old_ids=inp.previous_ids, table=table, identity=identity,
                                        rekeys=inp.rekeys, role_of=inp.role_of_old, dropped=dropped,
                                        normalize=old_norm, id_index=id_index,
                                        previous_entries=inp.previous_entries,
                                        is_isotopologue=is_isotopologue_name, is_class=chebi.is_class,
                                        old_string_to_primary=inp.old_string_to_primary,
                                        string_to_primary=inp.string_to_primary, chebi=chebi)
    errors += [f"remap {r['old_id']} → {r['new_id']} ({r['basis']}): not a compound step 6 writes"
               for r in remap if r["new_id"] and r["new_id"] not in table]
    errors += _identifier_coverage_errors(table, identity, inp.rekeys, dropped, set(dropped_unmoved))
    errors += [f"xref override {x['id']} → {x['to']}: target dropped at stage K"
               for x in xref_applied if x["to"] not in table]
    if errors:
        raise Step6ValidationError(errors)

    compounds, additional = {}, {}
    for n in sorted(table):
        e = {k: v for k, v in table[n].items() if k not in _TRANSIENT and v not in (None, [], "", {})}
        if n.startswith(_KEGG):
            compounds[n[len(_KEGG):]] = e
        else:
            additional[n] = e
    go_targets: dict[str, dict[str, list[str]]] = {}
    for g in sorted(set(go_subs) | set(go_roles)):
        go_targets[g] = dict(sorted({**go_subs.get(g, {}), **go_roles.get(g, {})}.items()))
    go_excluded_targets = sorted({
        (c if c in table else cid)
        for cid, p in inp.go_filter_primaries.items()
        for c in [identity.canonical_of.get(p, p) if p else None]})
    tcdb_excluded_targets = sorted({
        (c if c in table else cid)
        for cid, p in inp.tcdb_filter_primaries.items()
        for c in [identity.canonical_of.get(p, p) if p else None]})
    kegg_data = dict(kegg_rest, pathways=pathways, subcategories=subcategories, categories=categories,
                     compounds=compounds, additional_compounds=additional,
                     go_transport_substrates=go_targets,
                     go_transport_excluded_targets=go_excluded_targets,
                     go_transport_release=inp.go_release,
                     tcdb_substrate_excluded_targets=tcdb_excluded_targets)
    go_metabolite_targets = {t for m in go_subs.values() for t in m}
    go_report = {
        "release": inp.go_release,
        "axioms_in_corpus": inp.go_axioms_in_corpus,
        "edges": sum(len(m) for m in go_targets.values()),
        "metabolite_edges": sum(len(m) for m in go_subs.values()),
        "role_edges": sum(len(m) for m in go_roles.values()),
        "terms": len(go_targets),
        "targets": len({t for m in go_targets.values() for t in m}),
        "metabolite_targets": len(go_metabolite_targets),
        "role_targets": len({t for m in go_roles.values() for t in m}),
        # written only because of GO: no catalysis, TCDB or measurement reason
        "new_metabolites": sum(1 for n in go_metabolite_targets
                               if not set(table[n].get("evidence_sources") or [])
                               & {"metabolism", "transport", "metabolomics"}),
        "drops": dict(sorted({**inp.go_drops, "filtered_after_canonical": go_filtered_after}.items())),
        "unresolved_chebi_ids": sorted(inp.go_unresolved_ids),
        "excluded_targets": go_excluded_targets,
        "distribution": target_distribution(inp.go_distribution_input, go_targets),
    }
    tcdb_pruned = dict(inp.tcdb_static, system_substrates=dict(sorted(subs.items())),
                       chemical_roles=inp.chemical_roles, system_roles=inp.system_roles,
                       substrate_filter={"excluded_targets": tcdb_excluded_targets,
                                         "dropped_pairs": inp.tcdb_filter_drops,
                                         "dropped_after_canonical": tcdb_filtered_after})

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
            with chebi.untracked():
                old_n = old and old_norm(old)
            reason = repoint_reason(old_id=old_n, new=d, row=row, identity=identity,
                                    rekeys=inp.rekeys, dropped=dropped, case_fold_keys=case_fold_keys,
                                    written=table)
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
    # Like id_aliases: only names whose node survived stage K.
    res.name_lookup = {k: v for k, v in res.name_lookup.items() if v in table}
    mapping = build_mapping_v2(id_aliases=id_aliases,
                               resolution=res, overrides=inp.paper_overrides, tables=inp.tables,
                               compounds_summary=compounds_summary,
                               per_paper_stats=dict(sorted(per_paper.items())))

    report = build_identity_report(
        # remap / repoints / dropped are diagnostics against this baseline; a no-change rerun
        # over the committed outputs may leave them empty. The baseline string says so, so a
        # reader of the committed report does not take an emptied remap for lost data.
        baseline=_REPORT_BASELINE,
        chebi_release=chebi.release,
        merge_groups={c: {"members": m, "basis": identity.basis.get(c, [])}
                      for c, m in sorted(identity.members.items())},
        merge_groups_by_basis=_groups_by_basis(identity),
        id_aliases=id_aliases,
        rekeys=inp.rekeys, shared_chebi_terms=shared, rejected_links=[[l.a, l.b, l.basis] for l in identity.rejected],
        stereo_refused_links=identity.stereo_refused, mnx_split=mnx_split,
        canonical_pinned=identity.canonical_pinned, anomer_xref_review=anomer_xrefs,
        stereo_link_review=stereo_report, chebi_node_structure=chebi_structure,
        refused_assay_collisions=identity.refused_collisions,
        stale_identity_override=identity.stale_overrides + xref_stale, xref_overrides_applied=xref_applied, stale_alias_override=res.stale_overrides,
        alias_stereo_review=res.alias_stereo_review,
        stale_lock_entry=res.stale_lock, lock_target_remapped=res.lock_remapped,
        roles_removed=inp.role_of_old, dropped=dropped, dropped_without_rows=dropped_unmoved,
        dropped_same_entity_review=[{"dropped": d, "target": t,
                                     "dropped_name": old_entries.get(d, {}).get("name"),
                                     "target_name": table[t].get("name")}
                                    for d, t in sorted(dropped.items())],
        heuristic_name_candidates=_exact_name_pairs(
            table, inp.identity_overrides.keep_separate if inp.identity_overrides else frozenset()),
        ambiguous_paper_names=res.ambiguous, unresolved_paper_names=res.unresolved,
        dominant_form_drift=res.drift, structure=structure,
        go_transport=go_report,
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
