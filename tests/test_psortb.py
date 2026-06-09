"""Tests for the PSORTb subcellular-localization integration.

Covers the pure vocabulary helpers (multiomics_kg/utils/psortb.py) and the
ontology adapter (per-strain + multi-strain orchestrator). All pure / fixture-
backed — no live Neo4j, no network.
"""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.psortb_adapter import (
    SubcellularLocalizationAdapter,
    MultiSubcellularLocalizationAdapter,
    _node_id,
)
from multiomics_kg.utils.psortb import (
    LOCALIZATION_VOCAB,
    UNKNOWN_SENTINEL,
    display_name,
    is_kept,
)


# ── pure vocab / helpers ──────────────────────────────────────────────────────

def test_vocab_has_five_real_classes_and_excludes_unknown():
    assert set(LOCALIZATION_VOCAB) == {
        "Cytoplasmic", "CytoplasmicMembrane", "Periplasmic",
        "OuterMembrane", "Extracellular",
    }
    assert UNKNOWN_SENTINEL == "Unknown"
    assert UNKNOWN_SENTINEL not in LOCALIZATION_VOCAB


def test_is_kept_skips_sentinel_and_nulls():
    assert is_kept("OuterMembrane") is True
    assert is_kept("Cytoplasmic") is True
    assert is_kept("Unknown") is False
    assert is_kept(None) is False
    assert is_kept("") is False
    assert is_kept("Nonsense") is False


def test_display_name_maps_membrane_classes():
    assert display_name("OuterMembrane") == "Outer membrane"
    assert display_name("CytoplasmicMembrane") == "Cytoplasmic membrane"
    assert display_name("Cytoplasmic") == "Cytoplasmic"
    # defensive fallback
    assert display_name("Nonsense") == "Nonsense"


# ── adapter fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def strain_dir(tmp_path):
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001",
                     "psortb_localization": "OuterMembrane", "psortb_score": 9.97},
        "PMM_0002": {"locus_tag": "PMM_0002",
                     "psortb_localization": "Cytoplasmic", "psortb_score": 8.96},
        "PMM_0003": {"locus_tag": "PMM_0003",
                     "psortb_localization": "Unknown", "psortb_score": 2.0},  # sentinel
        "PMM_0004": {"locus_tag": "PMM_0004"},  # no psortb call at all
        "PMM_0005": {"locus_tag": "PMM_0005",
                     "psortb_localization": "Extracellular", "psortb_score": 10.0},
    }))
    return d


def _make_orchestrator(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    return MultiSubcellularLocalizationAdapter(
        genome_config_file=str(config_csv), test_mode=False,
    )


# ── per-strain adapter ────────────────────────────────────────────────────────

def test_per_strain_get_all_localizations_excludes_unknown(strain_dir):
    a = SubcellularLocalizationAdapter(genome_dir=strain_dir)
    assert a.get_all_localizations() == {"OuterMembrane", "Cytoplasmic", "Extracellular"}


def test_per_strain_edges_skip_unknown_and_missing(strain_dir):
    a = SubcellularLocalizationAdapter(genome_dir=strain_dir)
    edges = list(a.get_edges())
    sources = {e[1] for e in edges}
    # PMM_0003 (Unknown sentinel) and PMM_0004 (no call) yield no edge.
    # Gene node ids normalize via bioregistry → "ncbigene:" (colon); only the
    # psortb target falls back to "psortb_" (underscore, non-bioregistry prefix).
    assert sources == {"ncbigene:PMM_0001", "ncbigene:PMM_0002", "ncbigene:PMM_0005"}


def test_per_strain_edge_carries_score_and_label(strain_dir):
    a = SubcellularLocalizationAdapter(genome_dir=strain_dir)
    edges = {e[1]: e for e in a.get_edges()}
    eid, src, tgt, label, props = edges["ncbigene:PMM_0001"]
    assert label == "gene_has_subcellular_localization"
    assert tgt == _node_id("OuterMembrane") == "psortb_OuterMembrane"
    assert props == {"score": 9.97}
    assert isinstance(props["score"], float)


# ── orchestrator ──────────────────────────────────────────────────────────────

def test_orchestrator_emits_only_observed_real_classes(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    node_ids = {n[0] for n in a.get_nodes()}
    assert node_ids == {
        "psortb_OuterMembrane", "psortb_Cytoplasmic", "psortb_Extracellular",
    }
    # No Unknown node, no unobserved-class nodes.
    assert "psortb_Unknown" not in node_ids
    assert "psortb_Periplasmic" not in node_ids


def test_orchestrator_nodes_are_flat_level_zero(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    nodes = {n[0]: n for n in a.get_nodes()}
    _id, label, props = nodes["psortb_OuterMembrane"]
    assert label == "subcellular localization"
    assert props["level"] == 0
    assert "level_kind" not in props        # flat ontology
    assert props["psortb_id"] == "OuterMembrane"
    assert props["name"] == "Outer membrane"


def test_orchestrator_emits_no_hierarchy_edges(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    labels = {e[3] for e in a.get_edges()}
    # Flat ontology: only gene edges, no <x>_is_a_<x>.
    assert labels == {"gene_has_subcellular_localization"}


def test_orchestrator_gene_edges_delegate_to_strain(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    gene_edges = {(e[1], e[2]) for e in a.get_edges()}
    assert ("ncbigene:PMM_0001", "psortb_OuterMembrane") in gene_edges
    assert ("ncbigene:PMM_0005", "psortb_Extracellular") in gene_edges
    # Unknown / missing produce no edge.
    assert all(src != "ncbigene:PMM_0003" for src, _ in gene_edges)
    assert all(src != "ncbigene:PMM_0004" for src, _ in gene_edges)
