"""Tests for TCDB ontology adapter (per-strain + multi-strain)."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.tcdb_adapter import (
    TcdbAnnotationAdapter,
    MultiTcdbAnnotationAdapter,
)


@pytest.fixture
def strain_dir(tmp_path):
    """Stage a minimal gene_annotations_merged.json for one strain."""
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001",
                     "transporter_classification": ["1.A.1.5.2"]},
        "PMM_0002": {"locus_tag": "PMM_0002",
                     "transporter_classification": ["3.A.1"]},
        "PMM_0003": {"locus_tag": "PMM_0003"},  # no TCDB
    }))
    return d


@pytest.fixture
def cache_root(tmp_path):
    """Stage a minimal cache_root with tcdb_hierarchy.json + tcdb_pruned.json."""
    tcdb_dir = tmp_path / "tcdb"
    tcdb_dir.mkdir(parents=True)
    (tcdb_dir / "tcdb_hierarchy.json").write_text(json.dumps({
        "1": {"name": "Channels and Pores", "level": 0,
              "level_kind": "tc_class", "parent": None},
        "1.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"name": "VIC family", "level": 2,
                  "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"name": "", "level": 3, "level_kind": "tc_subfamily",
                    "parent": "1.A.1"},
        "1.A.1.5.2": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "superfamily": "VIC Superfamily"},
    }))
    (tcdb_dir / "tcdb_pruned.json").write_text(json.dumps({
        "kept_tcdb_ids": ["1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2"],
        # Pre-rolled-up substrates: every ancestor of the leaf carries the
        # leaf's primaries. Step 6 emits this map directly; adapter just
        # iterates and yields edges.
        "subtree_substrates": {
            "1": ["chebi:9999", "kegg.compound:C00208"],
            "1.A": ["chebi:9999", "kegg.compound:C00208"],
            "1.A.1": ["chebi:9999", "kegg.compound:C00208"],
            "1.A.1.5": ["chebi:9999", "kegg.compound:C00208"],
            "1.A.1.5.2": ["chebi:9999", "kegg.compound:C00208"],
        },
    }))
    return tmp_path


def test_per_strain_get_all_tcdb_ids(strain_dir):
    a = TcdbAnnotationAdapter(genome_dir=strain_dir)
    assert a.get_all_tcdb_ids() == {"1.A.1.5.2", "3.A.1"}


def test_per_strain_get_edges(strain_dir):
    a = TcdbAnnotationAdapter(genome_dir=strain_dir)
    edges = list(a.get_edges())
    edge_targets = {(e[1], e[2], e[3]) for e in edges}
    assert ("ncbigene:PMM_0001", "tcdb:1.A.1.5.2", "gene_has_tcdb_family") in edge_targets
    assert ("ncbigene:PMM_0002", "tcdb:3.A.1", "gene_has_tcdb_family") in edge_targets
    assert len(edges) == 2  # PMM_0003 has no TCDB → no edge


def _make_orchestrator(cache_root, strain_dir):
    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv),
        cache_root=cache_root,
        test_mode=False,
        cache=True,
    )
    a.download_data(cache=True)
    return a


def test_orchestrator_emits_one_node_per_kept_id(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    node_ids = {n[0] for n in a.get_nodes()}
    assert node_ids == {"tcdb:1", "tcdb:1.A", "tcdb:1.A.1", "tcdb:1.A.1.5", "tcdb:1.A.1.5.2"}


def test_orchestrator_emits_hierarchy_edges(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    edges = list(a.get_edges())
    parent_edges = {(e[1], e[2]) for e in edges if e[3] == "tcdb_family_is_a_tcdb_family"}
    assert ("tcdb:1.A", "tcdb:1") in parent_edges
    assert ("tcdb:1.A.1", "tcdb:1.A") in parent_edges
    assert ("tcdb:1.A.1.5", "tcdb:1.A.1") in parent_edges
    assert ("tcdb:1.A.1.5.2", "tcdb:1.A.1.5") in parent_edges
    assert len(parent_edges) == 4


def test_orchestrator_rolls_up_substrate_edges_to_ancestors(cache_root, strain_dir):
    """Substrate edges are emitted from every ancestor of a leaf with substrates,
    not just the leaf. Recover leaf-only semantics by filtering on source level_kind."""
    a = _make_orchestrator(cache_root, strain_dir)
    edges = list(a.get_edges())
    sub_edges = {(e[1], e[2]) for e in edges if e[3] == "tcdb_family_transports_metabolite"}
    expected = set()
    for ancestor in ("tcdb:1", "tcdb:1.A", "tcdb:1.A.1", "tcdb:1.A.1.5", "tcdb:1.A.1.5.2"):
        for met in ("kegg.compound:C00208", "chebi:9999"):
            expected.add((ancestor, met))
    assert sub_edges == expected


def test_orchestrator_node_props(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    nodes = {n[0]: n for n in a.get_nodes()}

    # Class node: name from _TC_CLASS_NAMES, level 0
    nid, _label, props = nodes["tcdb:1"]
    assert props["level"] == 0
    assert props["level_kind"] == "tc_class"
    assert props["name"] == "Channels and Pores"
    assert props["tcdb_id"] == "1"

    # Leaf node: name falls back to tcdb_id when source name is empty;
    # superfamily is set
    nid, _label, props = nodes["tcdb:1.A.1.5.2"]
    assert props["level"] == 4
    assert props["level_kind"] == "tc_specificity"
    assert props["tcdb_id"] == "1.A.1.5.2"
    assert props["name"] == "1.A.1.5.2"  # fallback
    assert props["superfamily"] == "VIC Superfamily"


def test_orchestrator_drops_gene_edges_for_unpruned_ids(cache_root, strain_dir):
    """Gene→TCDB edges to TCIDs absent from BOTH kept_tcdb_ids and seed_aliases
    are dropped (e.g. retired TCDB IDs with no curated ancestor).
    """
    a = _make_orchestrator(cache_root, strain_dir)
    edges = list(a.get_edges())
    gene_edges = [(e[1], e[2]) for e in edges if e[3] == "gene_has_tcdb_family"]
    targets = {t for _, t in gene_edges}
    assert "tcdb:3.A.1" not in targets
    assert ("ncbigene:PMM_0001", "tcdb:1.A.1.5.2") in gene_edges


def test_orchestrator_remaps_gene_edges_via_seed_aliases(tmp_path):
    """When seed_aliases is present, gene edges to retired TCIDs are re-anchored
    onto the nearest curated ancestor instead of being dropped."""
    # Strain has a gene annotated with retired ID 1.A.1.99.X
    strain_dir = tmp_path / "MED4"
    strain_dir.mkdir()
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001",
                     "transporter_classification": ["1.A.1.99.X"]},
    }))

    cache_root = tmp_path / "cache"
    tcdb_dir = cache_root / "tcdb"
    tcdb_dir.mkdir(parents=True)
    (tcdb_dir / "tcdb_hierarchy.json").write_text(json.dumps({
        "1": {"name": "Channels", "level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"name": "VIC", "level": 2, "level_kind": "tc_family", "parent": "1.A"},
    }))
    (tcdb_dir / "tcdb_pruned.json").write_text(json.dumps({
        "kept_tcdb_ids": ["1", "1.A", "1.A.1"],
        "subtree_substrates": {},
        "seed_aliases": {"1.A.1.99.X": "1.A.1"},
    }))

    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv), cache_root=cache_root,
        test_mode=False, cache=True,
    )
    a.download_data(cache=True)

    edges = list(a.get_edges())
    gene_edges = [(e[1], e[2]) for e in edges if e[3] == "gene_has_tcdb_family"]
    # Edge re-anchored to the family-level ancestor, not dropped.
    assert ("ncbigene:PMM_0001", "tcdb:1.A.1") in gene_edges


# ============================================================================
# Cross-ontology bridge edges (TcdbFamily -> Pfam / GO)
# ============================================================================


def _genome_csv(cache_root, strain_dir) -> str:
    """Write the minimal genomes.csv the orchestrator reads. Returns its path."""
    csv_path = cache_root / "genomes.csv"
    csv_path.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    return str(csv_path)


@pytest.fixture
def bridge_cache(cache_root):
    """Add pfam_bridge / go_bridge blocks to the staged tcdb_pruned.json."""
    p = cache_root / "tcdb" / "tcdb_pruned.json"
    data = json.loads(p.read_text())
    data["pfam_bridge"] = {
        "1.A.1": {"PF00001": ["1.A.1.9.1"], "PF99999": ["1.A.1.9.2"]},
    }
    data["go_bridge"] = {
        "1.A.1": {
            "GO:0005215": ["1.A.1.9.1"],   # molecular function
            "GO:0006810": ["1.A.1.9.1"],   # biological process
            "GO:0016020": ["1.A.1.9.1"],   # cellular component
            "GO:0000001": ["1.A.1.9.1"],   # absent from the KG -> must be skipped
        },
    }
    p.write_text(json.dumps(data))
    return cache_root


_GO_TERMS = {
    "GO:0005215": "molecular function",
    "GO:0006810": "biological process",
    "GO:0016020": "cellular component",
}


def _bridge_edges(adapter):
    return [e for e in adapter.get_edges() if e[3].startswith("tcdb_family_")
            and e[3] not in ("tcdb_family_is_a_tcdb_family",
                             "tcdb_family_transports_metabolite")]


def test_bridges_absent_when_node_sets_not_injected(bridge_cache, strain_dir):
    """No injected endpoint set -> no bridge edges (they could dangle)."""
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
    )
    a.download_data()
    assert _bridge_edges(a) == []


def test_pfam_bridge_pruned_to_existing_pfam_nodes(bridge_cache, strain_dir):
    """PF99999 has no Pfam node, so its edge must not be emitted."""
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        pfam_node_ids={"PF00001"},
    )
    a.download_data()
    edges = [e for e in _bridge_edges(a) if e[3] == "tcdb_family_has_pfam_domain"]
    assert len(edges) == 1
    _id, source, target, label, props = edges[0]
    assert source == "tcdb:1.A.1"
    assert target == "pfam:PF00001"
    assert props["curated_tcids"] == ["1.A.1.9.1"]


def test_pfam_bridge_direction_is_family_to_pfam(bridge_cache, strain_dir):
    """Direction is semantic: the family CONTAINS the domain.

    Reading these backwards (has domain -> is a transporter) is only ~31%
    precise on real data, so TcdbFamily must be the source.
    """
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        pfam_node_ids={"PF00001"},
    )
    a.download_data()
    e = next(e for e in _bridge_edges(a) if e[3] == "tcdb_family_has_pfam_domain")
    assert e[1].startswith("tcdb:") and e[2].startswith("pfam:")


def test_go_bridge_routes_each_namespace_to_its_own_edge_label(bridge_cache, strain_dir):
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        go_terms=_GO_TERMS,
    )
    a.download_data()
    by_label = {e[3]: e for e in _bridge_edges(a)}
    assert by_label["tcdb_family_enables_molecular_function"][2] == "go:0005215"
    assert by_label["tcdb_family_involved_in_biological_process"][2] == "go:0006810"
    assert by_label["tcdb_family_located_in_cellular_component"][2] == "go:0016020"


def test_go_bridge_skips_terms_with_no_node(bridge_cache, strain_dir):
    """GO:0000001 is not in the injected term map -> no edge, no dangle."""
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        go_terms=_GO_TERMS,
    )
    a.download_data()
    assert not any(e[2] == "go:0000001" for e in _bridge_edges(a))


def test_bridge_targets_are_always_kept_nodes(bridge_cache, strain_dir):
    """Bridge sources must exist in kept_tcdb_ids (dangling-proof on both ends)."""
    p = bridge_cache / "tcdb" / "tcdb_pruned.json"
    data = json.loads(p.read_text())
    data["pfam_bridge"]["9.Z.9"] = {"PF00001": ["9.Z.9.9.9"]}  # not a kept node
    p.write_text(json.dumps(data))
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        pfam_node_ids={"PF00001"},
    )
    a.download_data()
    assert not any(e[1] == "tcdb:9.Z.9" for e in _bridge_edges(a))


def test_bridges_do_not_touch_genes(bridge_cache, strain_dir):
    """Design invariant: bridges are ontology->ontology ONLY.

    Gene-level transporter identity must keep coming from Gene_has_tcdb_family
    alone; no bridge edge may introduce a gene endpoint.
    """
    a = MultiTcdbAnnotationAdapter(
        genome_config_file=_genome_csv(bridge_cache, strain_dir),
        cache_root=bridge_cache,
        pfam_node_ids={"PF00001"},
        go_terms=_GO_TERMS,
    )
    a.download_data()
    for _id, source, target, _label, _props in _bridge_edges(a):
        assert not source.startswith("ncbigene:")
        assert not target.startswith("ncbigene:")


# ============================================================================
# Two-source provenance on Gene_has_tcdb_family
# ============================================================================


@pytest.fixture
def dual_source_strain(tmp_path):
    """A strain whose genes carry eggNOG-only, diamond-only and BOTH TC calls."""
    d = tmp_path / "MED4"
    d.mkdir(parents=True, exist_ok=True)
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {  # both sources agree on 1.A.1.5.2
            "locus_tag": "PMM_0001", "protein_id": "WP_001.1",
            "transporter_classification": ["1.A.1.5.2"],
            "tcdb_eggnog_ids": ["1.A.1.5.2"],
            "tcdb_diamond_ids": ["1.A.1.5.2"],
        },
        "PMM_0002": {  # eggNOG only
            "locus_tag": "PMM_0002", "protein_id": "WP_002.1",
            "transporter_classification": ["1.A.1.5.2"],
            "tcdb_eggnog_ids": ["1.A.1.5.2"],
        },
        "PMM_0003": {  # diamond only, tier 3
            "locus_tag": "PMM_0003", "protein_id": "WP_003.1",
            "transporter_classification": ["1.A.1"],
            "tcdb_diamond_ids": ["1.A.1"],
        },
    }))
    tcdb_dir = d / "tcdb"
    tcdb_dir.mkdir()
    (tcdb_dir / "MED4.tcdb.calls.json").write_text(json.dumps({
        "WP_001.1": {"calls": [{
            "tcid": "1.A.1.5.2", "tier": 1, "confidence_score": 0.85,
            "identity": 87.4, "qcov": 92.1, "scov": 89.7, "evalue": 1e-180,
            "length": 412, "consensus_n": 3, "consensus_agreement": "5_part",
            "incompletely_characterized": False,
        }]},
        "WP_003.1": {"calls": [{
            "tcid": "1.A.1", "tier": 3, "confidence_score": 0.135,
            "identity": 30.0, "qcov": 45.0, "scov": 30.0, "evalue": 1e-10,
            "length": 150, "consensus_n": 1, "consensus_agreement": "3_part",
            "incompletely_characterized": False,
        }]},
    }))
    return d


def _gene_edges(strain_dir):
    a = TcdbAnnotationAdapter(genome_dir=strain_dir)
    return {(e[1], e[2]): e[4] for e in a.get_edges()}


def test_both_sources_produce_ONE_edge_with_both_in_sources(dual_source_strain):
    """A TC called by eggNOG AND diamond is a single edge — that agreement is the
    corroboration signal, and two parallel edges would destroy it."""
    edges = _gene_edges(dual_source_strain)
    props = edges[("ncbigene:PMM_0001", "tcdb:1.A.1.5.2")]
    assert sorted(props["sources"]) == ["diamond", "eggnog"]
    assert len([k for k in edges if k[0] == "ncbigene:PMM_0001"]) == 1


def test_eggnog_only_edge_carries_no_diamond_evidence(dual_source_strain):
    props = _gene_edges(dual_source_strain)[("ncbigene:PMM_0002", "tcdb:1.A.1.5.2")]
    assert props["sources"] == ["eggnog"]
    for k in ("tier", "confidence_score", "identity", "qcov", "evalue", "consensus_n"):
        assert k not in props, f"{k} must be sparse on an eggNOG-only edge"


def test_diamond_edge_carries_evidence_from_calls_json(dual_source_strain):
    props = _gene_edges(dual_source_strain)[("ncbigene:PMM_0001", "tcdb:1.A.1.5.2")]
    assert props["tier"] == 1
    assert props["confidence_score"] == 0.85
    assert props["identity"] == 87.4
    assert props["qcov"] == 92.1
    assert props["consensus_n"] == 3
    assert isinstance(props["tier"], int)
    assert isinstance(props["identity"], float)


def test_tier3_diamond_only_edge_is_still_emitted(dual_source_strain):
    """Tier 3 is conservative remote homology, not noise — it becomes an edge.
    Suppression happens downstream in post-import's annotation_types gate."""
    props = _gene_edges(dual_source_strain)[("ncbigene:PMM_0003", "tcdb:1.A.1")]
    assert props["sources"] == ["diamond"]
    assert props["tier"] == 3


def test_missing_calls_json_degrades_to_sources_only(tmp_path):
    """A strain not yet tcdb-diamond-run still emits its eggNOG edges."""
    d = tmp_path / "XX"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "X_0001": {"locus_tag": "X_0001", "protein_id": "WP_9.1",
                   "transporter_classification": ["1.A.1"],
                   "tcdb_eggnog_ids": ["1.A.1"]},
    }))
    props = _gene_edges(d)[("ncbigene:X_0001", "tcdb:1.A.1")]
    assert props["sources"] == ["eggnog"]
    assert "tier" not in props
