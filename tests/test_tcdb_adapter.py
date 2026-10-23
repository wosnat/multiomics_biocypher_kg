"""Tests for the TCDB ontology adapter (per-strain + multi-strain).

TCDB cleanup (spec 2026-10-22, §5): one Gene_has_tcdb_family edge per
(gene, TC node, tcdb_call) stamped with tcdb_gate; substrate / role edges leave
tc_system nodes only; nodes carry direction instead of substrate_breadth.
"""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.tcdb_adapter import (
    MultiTcdbAnnotationAdapter,
    TcdbAnnotationAdapter,
)

_LEVEL_KINDS = ["tc_class", "tc_subclass", "tc_family", "tc_subfamily", "tc_system"]


def _node(tc: str, name: str = "", **extra) -> dict:
    level = tc.count(".")
    return {"name": name, "level": level, "level_kind": _LEVEL_KINDS[level],
            "parent": tc.rsplit(".", 1)[0] if "." in tc else None, **extra}


# A small hierarchy covering an ABC uptake branch, a named efflux family, and a
# range-only system (kept, no gene points at it, unnamed).
_HIER = {
    "1": _node("1", "Channels and Pores"),
    "1.A": _node("1.A"),
    "1.A.1": _node("1.A.1", "VIC family"),
    "1.A.1.5": _node("1.A.1.5"),
    "1.A.1.5.2": _node("1.A.1.5.2", superfamily="VIC Superfamily",
                       description="Epsilon leaf curated text with 'quote' and |pipe|."),
    "2": _node("2", "Electrochemical Potential-driven Transporters"),
    "2.A": _node("2.A", "Porters"),
    "2.A.6": _node("2.A.6", "The Resistance-Nodulation-Cell Division (RND) Efflux Family"),
    "2.A.6.1": _node("2.A.6.1"),
    "2.A.6.1.1": _node("2.A.6.1.1", "Multidrug pump AcrB"),
    "2.A.6.1.9": _node("2.A.6.1.9"),   # range-only, unnamed
    "3": _node("3", "Primary Active Transporters"),
    "3.A": _node("3.A"),
    "3.A.1": _node("3.A.1", "The ATP-binding Cassette (ABC) Superfamily"),
    "3.A.1.4": _node("3.A.1.4", "The HAAT Family"),
    "3.A.1.4.4": _node("3.A.1.4.4", "The high-affinity urea porter"),
    "3.A.1.4.5": _node("3.A.1.4.5", "Another urea porter"),
    "9": _node("9", "Incompletely characterized"),   # in the full tree, NOT kept
}
_KEPT = sorted(t for t in _HIER if t != "9")

_ROLES = {
    "chebi.role:24432": {"name": "biological role", "parents": [], "level": 0,
                         "level_kind": "chebi_role_branch",
                         "all_identifiers": ["chebi.role:24432", "chebi:24432"]},
    "chebi.role:26672": {"name": "siderophore", "parents": ["chebi.role:24432"], "level": 1,
                         "level_kind": "chebi_role",
                         "all_identifiers": ["chebi.role:26672", "chebi:26672"]},
}


def _hit(tcid, identity, bitscore, qcov=95.0, scov=90.0, evalue=1e-50):
    return {"tcid": tcid, "identity": identity, "qcov": qcov, "scov": scov,
            "evalue": evalue, "bitscore": bitscore}


def _call(tcid, tier, identity, qcov=90.0, scov=85.0, evalue=1e-30):
    return {"tcid": tcid, "tier": tier, "identity": identity, "qcov": qcov, "scov": scov,
            "evalue": evalue, "confidence_score": 0.5, "consensus_n": 3, "length": 300,
            "level_kind": _LEVEL_KINDS[tcid.count(".")]}


def _best(*hits):
    return {"rule": "best_hit_v1", "tcids": [h["tcid"] for h in hits],
            "identity": hits[0]["identity"], "hits": list(hits)}


def _setup(tmp_path, genes: dict, calls: dict | None, *, aliases=None, hier=None,
           kept=None, system_substrates=None, system_roles=None, extra=None,
           **adapter_kwargs):
    """Stage strain + cache, return a downloaded MultiTcdbAnnotationAdapter."""
    d = tmp_path / "MED4"
    d.mkdir(parents=True, exist_ok=True)
    (d / "gene_annotations_merged.json").write_text(json.dumps(genes))
    if calls is not None:
        (d / "tcdb").mkdir(exist_ok=True)
        (d / "tcdb" / "MED4.tcdb.calls.json").write_text(json.dumps(calls))
    cache = tmp_path / "cache"
    (cache / "tcdb").mkdir(parents=True, exist_ok=True)
    (cache / "tcdb" / "tcdb_hierarchy.json").write_text(json.dumps(hier or _HIER))
    pruned = {
        "kept_tcdb_ids": kept if kept is not None else _KEPT,
        "system_substrates": system_substrates if system_substrates is not None else {
            "1.A.1.5.2": ["chebi:9999", "kegg.compound:C00208"],
            "3.A.1.4.4": ["kegg.compound:C00086"],
            "2.A.6.1.9": ["kegg.compound:C00001"],
        },
        "system_roles": system_roles if system_roles is not None else {
            "1.A.1.5.2": ["chebi.role:26672"]},
        "chemical_roles": _ROLES,
        "seed_aliases": aliases or {},
        "pfam_bridge": {},
        "go_bridge": {},
    }
    pruned.update(extra or {})
    (cache / "tcdb" / "tcdb_pruned.json").write_text(json.dumps(pruned))
    csv_path = cache / "genomes.csv"
    csv_path.write_text(f"strain_name,data_dir\nMED4,{d}\n")
    a = MultiTcdbAnnotationAdapter(genome_config_file=str(csv_path), cache_root=cache,
                                   **adapter_kwargs)
    a.download_data()
    return a


def _gene_edges(adapter):
    return [e for e in adapter.get_edges() if e[3] == "gene_has_tcdb_family"]


def _by_call(adapter, locus_tag="PMM_0001"):
    """{(target, tcdb_call): props} for one gene."""
    out = {}
    for _id, src, tgt, _lbl, props in _gene_edges(adapter):
        if src == f"ncbigene:{locus_tag}":
            out[(tgt, props["tcdb_call"])] = props
    return out


def _gene(lt="PMM_0001", wp="WP_1", eggnog=(), **extra):
    return {lt: {"locus_tag": lt, "protein_id": wp, "tcdb_eggnog_ids": list(eggnog), **extra}}


# ============================================================================
# Gene edges: one per tcdb_call
# ============================================================================


def test_one_edge_per_tcdb_call(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.4.4"]), {
        "WP_1": {"calls": [_call("3.A.1.4.4", 1, 88.0)],
                 "best_hit": _best(_hit("3.A.1.4.4", 88.0, 400.0))}})
    edges = _gene_edges(a)
    ids = sorted(e[0] for e in edges)
    assert ids == ["PMM_0001-has_tcdb-3.A.1.4.4-diamond_best_hit",
                   "PMM_0001-has_tcdb-3.A.1.4.4-diamond_tier",
                   "PMM_0001-has_tcdb-3.A.1.4.4-eggnog"]
    assert {e[2] for e in edges} == {"tcdb:3.A.1.4.4"}
    by = {e[4]["tcdb_call"]: e[4] for e in edges}
    assert by["eggnog"]["sources"] == ["eggnog"]
    assert by["eggnog"]["evidence"] == "family_inferred"
    for kind in ("diamond_tier", "diamond_best_hit"):
        assert by[kind]["sources"] == ["tcdb_diamond"]
        assert by[kind]["evidence"] == "homology"


def test_edge_properties_by_call_kind(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.4"]), {
        "WP_1": {"calls": [_call("3.A.1.4", 2, 55.0)],
                 "best_hit": _best(_hit("3.A.1.4.4", 58.0, 300.0),
                                   _hit("3.A.1.4.5", 52.0, 280.0))}})
    by = _by_call(a)
    egg = by[("tcdb:3.A.1.4", "eggnog")]
    assert set(egg) == {"sources", "evidence", "tcdb_call", "tcdb_gate"}
    tier = by[("tcdb:3.A.1.4", "diamond_tier")]
    assert set(tier) == {"sources", "evidence", "tcdb_call", "tcdb_gate",
                         "tier", "identity", "qcov", "scov", "evalue"}
    assert tier["tier"] == 2 and isinstance(tier["tier"], int)
    assert tier["identity"] == 55.0 and tier["scov"] == 85.0
    bh = by[("tcdb:3.A.1.4.5", "diamond_best_hit")]
    assert set(bh) == {"sources", "evidence", "tcdb_call", "tcdb_gate", "rank", "tier",
                       "identity", "qcov", "scov", "evalue", "bitscore", "bitscore_ratio"}
    assert bh["rank"] == 2 and bh["tier"] == 2 and bh["identity"] == 52.0
    assert bh["bitscore"] == 280.0
    assert bh["bitscore_ratio"] == pytest.approx(280.0 / 300.0, abs=1e-4)
    for props in by.values():
        for k in props:
            assert k not in ("confidence_score", "consensus_n") and not k.startswith("best_hit")


@pytest.mark.parametrize("kind, tc, tier_or_rank, identity, gate", [
    ("diamond_best_hit", "3.A.1.4.4", 1, 85.0, "direct"),
    ("diamond_best_hit", "3.A.1.4.4", 1, 39.0, "below_gate"),
    ("diamond_best_hit", "3.A.1.4.5", 2, 85.0, "below_gate"),
    ("eggnog", "3.A.1.4.4", None, None, "direct"),
    ("eggnog", "3.A.1.4", None, None, "annotated"),
    ("diamond_tier", "3.A.1.4.4", 1, 80.0, "direct"),
    ("diamond_tier", "3.A.1.4", 2, 50.0, "annotated"),
    ("diamond_tier", "3.A.1", 3, 30.0, "below_gate"),
])
def test_tcdb_gate_rules(tmp_path, kind, tc, tier_or_rank, identity, gate):
    calls: dict = {"WP_1": {"calls": []}}
    eggnog: list = []
    if kind == "eggnog":
        eggnog = [tc]
    elif kind == "diamond_tier":
        calls["WP_1"]["calls"] = [_call(tc, tier_or_rank, identity)]
    else:
        hits = [_hit(tc, identity, 400.0)]
        if tier_or_rank == 2:
            hits = [_hit("3.A.1.4.4", 30.0, 410.0), _hit(tc, identity, 400.0)]
        calls["WP_1"]["best_hit"] = _best(*hits)
    by = _by_call(_setup(tmp_path, _gene(eggnog=eggnog), calls))
    props = by[(f"tcdb:{tc}", kind)]
    assert props["tcdb_gate"] == gate


@pytest.mark.parametrize("kind, gate", [
    ("diamond_tier", "annotated"),       # tier-1 call re-anchored onto a subfamily
    ("diamond_best_hit", "annotated"),   # rank-1 tier-1 best hit re-anchored
    ("eggnog", "annotated"),
])
def test_direct_requires_a_tc_system_target(tmp_path, kind, gate):
    tc = "3.A.1.4.99"   # not kept, aliased to the subfamily 3.A.1.4
    calls: dict = {"WP_1": {"calls": []}}
    eggnog: list = []
    if kind == "eggnog":
        eggnog = [tc]
    elif kind == "diamond_tier":
        calls["WP_1"]["calls"] = [_call(tc, 1, 85.0)]
    else:
        calls["WP_1"]["best_hit"] = _best(_hit(tc, 85.0, 400.0))
    a = _setup(tmp_path, _gene(eggnog=eggnog), calls, aliases={tc: "3.A.1.4"})
    props = _by_call(a)[("tcdb:3.A.1.4", kind)]
    assert props["tcdb_gate"] == gate


def test_eggnog_alias_to_ancestor_is_annotated(tmp_path):
    """A retired eggNOG 5-part id re-anchored to its subfamily is NOT direct."""
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.4.99"]), None,
               aliases={"3.A.1.4.99": "3.A.1.4"})
    edges = _gene_edges(a)
    assert len(edges) == 1
    _id, _s, tgt, _l, props = edges[0]
    assert _id == "PMM_0001-has_tcdb-3.A.1.4-eggnog"
    assert tgt == "tcdb:3.A.1.4"
    assert props["tcdb_gate"] == "annotated"


def test_eggnog_alias_dropped_emits_nothing(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["7.Z.1.1.1", "8.Q.2"]), None,
               aliases={"7.Z.1.1.1": ""})
    assert _gene_edges(a) == []


def test_two_eggnog_ids_aliasing_to_one_node_give_one_edge(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.35", "3.A.1.36"]), None,
               aliases={"3.A.1.35": "3.A.1", "3.A.1.36": "3.A.1"})
    edges = _gene_edges(a)
    assert [e[0] for e in edges] == ["PMM_0001-has_tcdb-3.A.1-eggnog"]


def test_within_kind_duplicates_keep_higher_identity(tmp_path):
    """Two tier calls landing on one node (here via an alias) collapse to one edge
    carrying the higher-identity call's statistics."""
    a = _setup(tmp_path, _gene(), {"WP_1": {"calls": [
        _call("3.A.1.4", 2, 45.0, qcov=70.0), _call("3.A.1.4.77", 1, 72.0, qcov=99.0)]}},
        aliases={"3.A.1.4.77": "3.A.1.4"})
    edges = [e for e in _gene_edges(a) if e[4]["tcdb_call"] == "diamond_tier"]
    assert len(edges) == 1
    assert edges[0][4]["identity"] == 72.0 and edges[0][4]["qcov"] == 99.0


def test_within_kind_best_hit_duplicates_keep_rank_one(tmp_path):
    """Two best-hit systems aliasing to one node collapse to the rank-1 hit, not the
    higher-identity one; the gate follows rank 1 (rank 2 would read below_gate)."""
    a = _setup(tmp_path, _gene(), {"WP_1": {"calls": [], "best_hit": _best(
        _hit("3.A.1.4.77", 50.0, 400.0), _hit("3.A.1.4.78", 85.0, 390.0))}},
        aliases={"3.A.1.4.77": "3.A.1.4.4", "3.A.1.4.78": "3.A.1.4.4"})
    edges = [e for e in _gene_edges(a) if e[4]["tcdb_call"] == "diamond_best_hit"]
    assert len(edges) == 1
    props = edges[0][4]
    assert edges[0][2] == "tcdb:3.A.1.4.4"
    assert props["rank"] == 1 and props["identity"] == 50.0
    assert props["tcdb_gate"] == "direct"


def test_best_hit_without_call_emits_best_hit_edge(tmp_path):
    """Proteins with a best hit but no tier call (top HSP < 50 aa) get an edge."""
    a = _setup(tmp_path, _gene(), {"WP_1": {"calls": [],
                                            "best_hit": _best(_hit("3.A.1.4.4", 75.0, 90.0))}})
    edges = _gene_edges(a)
    assert len(edges) == 1
    props = edges[0][4]
    assert props["tcdb_call"] == "diamond_best_hit"
    assert props["rank"] == 1 and props["tier"] == 1 and props["tcdb_gate"] == "direct"


def test_missing_calls_json_yields_eggnog_only(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.4.4"],
                               tcdb_diamond_ids=["3.A.1.4"],
                               tcdb_best_hit_ids=["3.A.1.4.4"]), None)
    edges = _gene_edges(a)
    assert [e[4]["tcdb_call"] for e in edges] == ["eggnog"]


def test_bitscore_ratio_is_one_at_rank_one(tmp_path):
    a = _setup(tmp_path, _gene(), {"WP_1": {"calls": [], "best_hit": _best(
        _hit("3.A.1.4.4", 60.0, 512.5), _hit("3.A.1.4.5", 58.0, 470.0))}})
    by = _by_call(a)
    assert by[("tcdb:3.A.1.4.4", "diamond_best_hit")]["bitscore_ratio"] == 1.0
    assert by[("tcdb:3.A.1.4.5", "diamond_best_hit")]["bitscore_ratio"] >= 0.9


def test_unkept_targets_dropped(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["9"]), {"WP_1": {"calls": [_call("9.B.1", 3, 20.0)]}})
    assert _gene_edges(a) == []


def test_shared_protein_id_gives_both_genes_edges(tmp_path):
    genes = {**_gene("PMM_0002", "WP_2"), **_gene("PMM_0003", "WP_2")}
    a = _setup(tmp_path, genes, {"WP_2": {"calls": [],
                                          "best_hit": _best(_hit("3.A.1.4.4", 90.0, 99.0))}})
    assert {e[1] for e in _gene_edges(a)} == {"ncbigene:PMM_0002", "ncbigene:PMM_0003"}


def test_gene_edge_ids_unique(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["3.A.1.4.4", "3.A.1.4"]), {"WP_1": {
        "calls": [_call("3.A.1.4.4", 1, 80.0), _call("3.A.1.4", 2, 50.0)],
        "best_hit": _best(_hit("3.A.1.4.4", 80.0, 400.0), _hit("3.A.1.4.5", 79.0, 390.0))}})
    ids = [e[0] for e in _gene_edges(a)]
    assert len(ids) == len(set(ids)) == 6


def test_per_strain_yields_raw_claims(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps(_gene(eggnog=["3.A.1.4.4"])))
    (d / "tcdb").mkdir()
    (d / "tcdb" / "S.tcdb.calls.json").write_text(json.dumps({"WP_1": {
        "calls": [_call("3.A.1.4", 2, 50.0)],
        "best_hit": _best(_hit("3.A.1.4.4", 80.0, 400.0))}}))
    claims = list(TcdbAnnotationAdapter(genome_dir=d).get_edges())
    assert [(lt, kind, tc) for lt, kind, tc, _p in claims] == [
        ("PMM_0001", "eggnog", "3.A.1.4.4"),
        ("PMM_0001", "diamond_tier", "3.A.1.4"),
        ("PMM_0001", "diamond_best_hit", "3.A.1.4.4"),
    ]


# ============================================================================
# Nodes
# ============================================================================


def test_node_props(tmp_path):
    a = _setup(tmp_path, {}, None)
    nodes = {n[0]: n[2] for n in a.get_nodes() if n[1] == "tcdb family"}
    assert set(nodes) == {f"tcdb:{t}" for t in _KEPT}
    for props in nodes.values():
        assert props["level_kind"] in _LEVEL_KINDS
        assert "substrate_breadth" not in props
        assert props["direction"] in {"uptake", "export", "unknown"}
        assert ("direction_basis" in props) == (props["direction"] != "unknown")
    assert nodes["tcdb:1"]["name"] == "Channels and Pores"
    leaf = nodes["tcdb:1.A.1.5.2"]
    assert leaf["level"] == 4 and leaf["level_kind"] == "tc_system"
    assert leaf["name"] == "1.A.1.5.2"
    assert leaf["superfamily"] == "VIC Superfamily"
    assert leaf["description"] == "Epsilon leaf curated text with ^quote^ and pipe."
    assert "description" not in nodes["tcdb:1.A.1.5"]
    # ABC numbering + named family
    assert (nodes["tcdb:3.A.1.4.4"]["direction"], nodes["tcdb:3.A.1.4.4"]["direction_basis"]) \
        == ("uptake", "abc_numbering")
    assert nodes["tcdb:3.A.1"]["direction"] == "unknown"
    assert (nodes["tcdb:2.A.6"]["direction"], nodes["tcdb:2.A.6"]["direction_basis"]) \
        == ("export", "tcdb_name")


def test_unnamed_system_falls_back_and_inherits_direction(tmp_path):
    a = _setup(tmp_path, {}, None)
    nodes = {n[0]: n[2] for n in a.get_nodes() if n[1] == "tcdb family"}
    p = nodes["tcdb:2.A.6.1.9"]
    assert p["name"] == "2.A.6.1.9"
    assert p["direction"] == "export" and p["direction_basis"] == "tcdb_name"


def test_orchestrator_emits_hierarchy_edges(tmp_path):
    a = _setup(tmp_path, {}, None)
    parents = {(e[1], e[2]) for e in a.get_edges() if e[3] == "tcdb_family_is_a_tcdb_family"}
    assert ("tcdb:1.A.1.5.2", "tcdb:1.A.1.5") in parents
    assert ("tcdb:2.A.6.1.9", "tcdb:2.A.6.1") in parents
    assert len(parents) == len(_KEPT) - 3   # three kept roots


# ============================================================================
# Substrate / role edges: tc_system only, no properties
# ============================================================================


def test_substrate_edges_only_from_systems(tmp_path):
    a = _setup(tmp_path, {}, None)
    edges = list(a.get_edges())
    sub = [e for e in edges if e[3] in ("tcdb_family_transports_metabolite",
                                        "tcdb_family_transports_chemical_role")]
    assert sub
    for _id, src, _tgt, _lbl, props in sub:
        assert src.removeprefix("tcdb:").count(".") == 4
        assert props == {}
    mets = {(e[0], e[1], e[2]) for e in sub if e[3] == "tcdb_family_transports_metabolite"}
    assert ("1.A.1.5.2-transports-chebi:9999", "tcdb:1.A.1.5.2", "chebi:9999") in mets
    assert ("2.A.6.1.9-transports-kegg.compound:C00001", "tcdb:2.A.6.1.9",
            "kegg.compound:C00001") in mets
    roles = [(e[0], e[1], e[2]) for e in sub if e[3] == "tcdb_family_transports_chemical_role"]
    assert roles == [("t2r:1.A.1.5.2:chebi.role:26672", "tcdb:1.A.1.5.2", "chebi.role:26672")]


def test_substrate_edges_skip_unkept_and_non_system_keys(tmp_path):
    a = _setup(tmp_path, {}, None, system_substrates={
        "9.A.1.1.1": ["kegg.compound:C1"], "3.A.1.4": ["kegg.compound:C2"]},
        system_roles={})
    assert not [e for e in a.get_edges() if e[3] == "tcdb_family_transports_metabolite"]


def test_chemical_role_nodes_and_isa(tmp_path):
    a = _setup(tmp_path, {}, None)
    nodes = {n[0]: n for n in a.get_nodes()}
    r = nodes["chebi.role:26672"]
    assert r[1] == "chemical role" and r[2]["level_kind"] == "chebi_role" and r[2]["level"] == 1
    assert ("chebi.role:26672-isa-chebi.role:24432", "chebi.role:26672", "chebi.role:24432",
            "chemical_role_is_a_chemical_role", {}) in list(a.get_edges())


def test_legacy_pruned_keys_ignored(tmp_path):
    """A legacy file's roll-up keys never reach the graph."""
    a = _setup(tmp_path, {}, None, extra={
        "subtree_substrates": {"3.A.1": ["kegg.compound:C9"]},
        "lumping_pin": {"lumping_tcdb_ids": ["3.A.1"]}})
    assert not any(e[1] == "tcdb:3.A.1" and e[3] == "tcdb_family_transports_metabolite"
                   for e in a.get_edges())


# ============================================================================
# Cross-ontology bridge edges (TcdbFamily -> Pfam / GO), on all kept nodes
# ============================================================================

_GO_TERMS = {
    "GO:0005215": "molecular function",
    "GO:0006810": "biological process",
    "GO:0016020": "cellular component",
}

_BRIDGES = {
    "pfam_bridge": {
        "1.A.1": {"PF00001": ["1.A.1.9.1"], "PF99999": ["1.A.1.9.2"]},
        "2.A.6.1.9": {"PF00873": ["2.A.6.1.9"]},   # range-only system
        "9.Z.9": {"PF00001": ["9.Z.9.9.9"]},       # not a kept node
    },
    "go_bridge": {
        "1.A.1": {
            "GO:0005215": ["1.A.1.9.1"],
            "GO:0006810": ["1.A.1.9.1"],
            "GO:0016020": ["1.A.1.9.1"],
            "GO:0000001": ["1.A.1.9.1"],   # absent from the KG -> skipped
        },
    },
}

_NOT_BRIDGE = ("tcdb_family_is_a_tcdb_family", "tcdb_family_transports_metabolite",
               "tcdb_family_transports_chemical_role")


def _bridge_edges(adapter):
    return [e for e in adapter.get_edges()
            if e[3].startswith("tcdb_family_") and e[3] not in _NOT_BRIDGE]


def test_bridges_absent_when_node_sets_not_injected(tmp_path):
    a = _setup(tmp_path, {}, None, extra=_BRIDGES)
    assert _bridge_edges(a) == []


def test_pfam_bridge_pruned_to_existing_pfam_nodes(tmp_path):
    a = _setup(tmp_path, {}, None, extra=_BRIDGES, pfam_node_ids={"PF00001"})
    edges = [e for e in _bridge_edges(a) if e[3] == "tcdb_family_has_pfam_domain"]
    assert len(edges) == 1
    _id, source, target, _label, props = edges[0]
    assert (source, target) == ("tcdb:1.A.1", "pfam:PF00001")
    assert props["curated_tcids"] == ["1.A.1.9.1"]


def test_pfam_bridge_on_range_only_system(tmp_path):
    """Bridges cover range-only kept nodes too (spec §5.5)."""
    a = _setup(tmp_path, {}, None, extra=_BRIDGES, pfam_node_ids={"PF00001", "PF00873"})
    pairs = {(e[1], e[2]) for e in _bridge_edges(a)}
    assert ("tcdb:2.A.6.1.9", "pfam:PF00873") in pairs
    assert not any(s == "tcdb:9.Z.9" for s, _t in pairs)


def test_go_bridge_routes_each_namespace_to_its_own_edge_label(tmp_path):
    a = _setup(tmp_path, {}, None, extra=_BRIDGES, go_terms=_GO_TERMS)
    by_label = {e[3]: e for e in _bridge_edges(a)}
    assert by_label["tcdb_family_enables_molecular_function"][2] == "go:0005215"
    assert by_label["tcdb_family_involved_in_biological_process"][2] == "go:0006810"
    assert by_label["tcdb_family_located_in_cellular_component"][2] == "go:0016020"
    assert not any(e[2] == "go:0000001" for e in _bridge_edges(a))


def test_bridges_do_not_touch_genes(tmp_path):
    a = _setup(tmp_path, _gene(eggnog=["1.A.1"]), None, extra=_BRIDGES,
               pfam_node_ids={"PF00001"}, go_terms=_GO_TERMS)
    for _id, source, target, _label, _props in _bridge_edges(a):
        assert source.startswith("tcdb:")
        assert not target.startswith("ncbigene:")


# ============================================================================
# House rules
# ============================================================================


def test_no_schema_property_is_declared_bool():
    """R5: native bool is forbidden graph-wide."""
    import pathlib
    import re
    text = pathlib.Path("config/schema_config.yaml").read_text()
    offenders = re.findall(r"^\s+([a-z_]+): bool\s*(?:#.*)?$", text, re.M)
    assert offenders == [], f"native bool properties remain: {offenders}"


def test_is_promiscuous_is_gone_from_post_import():
    import pathlib
    for p in ("scripts/post-import.cypher", "scripts/post-import.sh"):
        text = pathlib.Path(p).read_text()
        assert "is_promiscuous" not in text, f"{p} still sets or reads is_promiscuous"


def test_string_props_are_sanitized(tmp_path):
    a = _setup(tmp_path, {}, None)
    for _nid, label, props in a.get_nodes():
        for v in props.values():
            vals = v if isinstance(v, list) else [v]
            for s in vals:
                if isinstance(s, str):
                    assert "'" not in s and "|" not in s, (label, s)


# ============================================================================
# Post-import (TCDB cleanup, spec §5.1 / §5.3 / §7)
# ============================================================================

_POST_IMPORT_FILES = ("scripts/post-import.sh", "scripts/post-import.cypher")


def _cypher_code(text: str) -> str:
    """Strip `//` comment lines so absence checks ignore history comments."""
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("//")
    )


def _tcdb_statements(text: str) -> list[str]:
    """Every Cypher statement that touches TCDB, comments stripped and
    whitespace-normalised, in file order."""
    import re
    if "<<'CYPHER'" in text:  # post-import.sh: keep only the heredoc bodies
        text = "\n".join(re.findall(r"<<'CYPHER'\n(.*?)\nCYPHER\n", text, re.S))
    code = _cypher_code(text)
    stmts = [re.sub(r"\s+", " ", s).strip() for s in code.split(";")]
    return [s for s in stmts if re.search(r"tcdb|Tcdb", s)]


@pytest.mark.parametrize("path", _POST_IMPORT_FILES)
def test_post_import_uses_tcdb_gate(path):
    import pathlib
    text = pathlib.Path(path).read_text()
    code = _cypher_code(text)
    # the three annotation gates (annotation_quality, annotation_types,
    # informative_annotation_types) + the TcdbFamily / Gene family counts
    assert code.count("tcdb_gate IN ['direct','annotated']") >= 5, path
    assert "tcdb_gate: 'direct'" in code, path
    assert "'eggnog' IN rt.sources OR rt.tier <= 2" not in code, path
    for gone in ("substrate_breadth", "attachment_depth", "substrate_depth",
                 "tcdb_evidence_score_max", "transport_substrate_resolution",
                 "source_agreement", "go_support", "t.metabolite_count",
                 "r.tcdb_family_count", "r.direct_tcdb_family_count"):
        assert gone not in code, f"{path} still references {gone!r}"
    for present in ("t.tcdb_substrate_count", "r.tcdb_system_count",
                    "r.direct_tcdb_system_count", "r.direct_transporter_count",
                    "g.tcdb_family_count", "g.transported_metabolite_count",
                    "m.transporter_count", "m.transporter_gene_count"):
        assert present in code, f"{path} missing {present!r}"


def test_post_import_sh_and_cypher_tcdb_blocks_match():
    import pathlib
    sh = _tcdb_statements(pathlib.Path("scripts/post-import.sh").read_text())
    cy = _tcdb_statements(pathlib.Path("scripts/post-import.cypher").read_text())
    assert len(sh) > 10
    assert sh == cy
