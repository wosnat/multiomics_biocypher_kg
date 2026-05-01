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
        "leaf_substrates": {
            "1.A.1.5.2": ["kegg.compound:C00208", "chebi:9999"],
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


def test_orchestrator_emits_substrate_edges_only_on_leaves(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    edges = list(a.get_edges())
    sub_edges = {(e[1], e[2]) for e in edges if e[3] == "tcdb_family_transports_metabolite"}
    assert sub_edges == {
        ("tcdb:1.A.1.5.2", "kegg.compound:C00208"),
        ("tcdb:1.A.1.5.2", "chebi:9999"),
    }


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
    """Gene→TCDB edges must only target TCDB IDs in kept_tcdb_ids.

    Per-strain adapter yields an edge to '3.A.1' (PMM_0002) but the orchestrator's
    pruned set has only the 1.* family. The orchestrator must drop the 3.A.1 edge.
    """
    a = _make_orchestrator(cache_root, strain_dir)
    edges = list(a.get_edges())
    gene_edges = [(e[1], e[2]) for e in edges if e[3] == "gene_has_tcdb_family"]
    targets = {t for _, t in gene_edges}
    assert "tcdb:3.A.1" not in targets
    assert ("ncbigene:PMM_0001", "tcdb:1.A.1.5.2") in gene_edges
