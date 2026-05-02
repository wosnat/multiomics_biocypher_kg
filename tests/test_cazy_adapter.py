"""Tests for CAZy ontology adapter (per-strain + multi-strain)."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.cazy_adapter import (
    CazyAnnotationAdapter,
    MultiCazyAnnotationAdapter,
)


@pytest.fixture
def strain_dir(tmp_path):
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001", "cazy_ids": ["GH13"]},
        "PMM_0002": {"locus_tag": "PMM_0002", "cazy_ids": ["GH13_5", "CBM48"]},
        "PMM_0003": {"locus_tag": "PMM_0003"},  # no CAZy
        "PMM_0004": {"locus_tag": "PMM_0004", "cazy_ids": ["GH13_5", "garbage"]},
    }))
    return d


def test_per_strain_get_all_cazy_ids_filters_garbage(strain_dir):
    a = CazyAnnotationAdapter(genome_dir=strain_dir)
    ids = a.get_all_cazy_ids()
    assert ids == {"GH13", "GH13_5", "CBM48"}
    assert "garbage" not in ids


def test_per_strain_get_edges_attaches_at_most_specific_level(strain_dir):
    a = CazyAnnotationAdapter(genome_dir=strain_dir)
    edges = list(a.get_edges())
    edge_targets = {(e[1], e[2]) for e in edges}
    assert ("ncbigene:PMM_0001", "cazy:GH13") in edge_targets
    assert ("ncbigene:PMM_0002", "cazy:GH13_5") in edge_targets
    assert ("ncbigene:PMM_0002", "cazy:CBM48") in edge_targets
    assert ("ncbigene:PMM_0004", "cazy:GH13_5") in edge_targets
    assert all("garbage" not in t for _, t in edge_targets)


def _make_orchestrator(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    return MultiCazyAnnotationAdapter(
        genome_config_file=str(config_csv),
        test_mode=False,
    )


def test_orchestrator_emits_class_family_subfamily(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    nodes = list(a.get_nodes())
    node_ids = {n[0] for n in nodes}
    assert node_ids == {"cazy:GH", "cazy:GH13", "cazy:GH13_5", "cazy:CBM", "cazy:CBM48"}


def test_orchestrator_class_node_has_human_readable_name(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    nodes = {n[0]: n for n in a.get_nodes()}
    _, _label, props = nodes["cazy:GH"]
    assert props["name"] == "Glycoside Hydrolases"
    assert props["level"] == 0
    assert props["level_kind"] == "cazy_class"
    # Family fallback: name == cazy_id when source has no human name
    _, _label, props = nodes["cazy:GH13"]
    assert props["name"] == "GH13"
    assert props["level"] == 1
    assert props["level_kind"] == "cazy_family"
    # Subfamily
    _, _label, props = nodes["cazy:GH13_5"]
    assert props["name"] == "GH13_5"
    assert props["level"] == 2
    assert props["level_kind"] == "cazy_subfamily"


def test_orchestrator_emits_parent_edges(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    edges = list(a.get_edges())
    parents = {(e[1], e[2]) for e in edges if e[3] == "cazy_family_is_a_cazy_family"}
    assert ("cazy:GH13", "cazy:GH") in parents
    assert ("cazy:GH13_5", "cazy:GH13") in parents
    assert ("cazy:CBM48", "cazy:CBM") in parents
    # No self-edges, no class→class edges
    assert ("cazy:GH", "cazy:GH") not in parents


def test_orchestrator_emits_gene_edges(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    edges = list(a.get_edges())
    gene_edges = {(e[1], e[2]) for e in edges if e[3] == "gene_has_cazy_family"}
    assert ("ncbigene:PMM_0001", "cazy:GH13") in gene_edges
    assert ("ncbigene:PMM_0002", "cazy:GH13_5") in gene_edges
    assert ("ncbigene:PMM_0002", "cazy:CBM48") in gene_edges
    assert ("ncbigene:PMM_0004", "cazy:GH13_5") in gene_edges


def test_orchestrator_does_not_emit_unobserved_classes(strain_dir, tmp_path):
    """Spec: 'observed-only — every emitted node has at least one gene-annotation source.'

    The fixture only annotates GH and CBM. The orchestrator must NOT emit GT, PL,
    CE, or AA class nodes (they have no observed members)."""
    a = _make_orchestrator(strain_dir, tmp_path)
    nodes = list(a.get_nodes())
    node_ids = {n[0] for n in nodes}
    assert "cazy:GT" not in node_ids
    assert "cazy:PL" not in node_ids
    assert "cazy:CE" not in node_ids
    assert "cazy:AA" not in node_ids
