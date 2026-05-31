"""Tests for the SignalP signal-peptide integration.

Covers the pure vocabulary + parser helpers (multiomics_kg/utils/signalp.py)
and the ontology adapter (per-strain + multi-strain orchestrator). All pure /
fixture-backed — no live Neo4j, no network, no Docker.
"""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.signalp_adapter import (
    SignalPeptideAdapter,
    MultiSignalPeptideAdapter,
    _node_id,
)
from multiomics_kg.utils.signalp import (
    OTHER_SENTINEL,
    PROB_COLUMNS,
    SIGNALP_VOCAB,
    display_name,
    first_token,
    is_kept,
    parse_cs_pos,
    parse_prediction_results,
)


# ── pure vocab / helpers ──────────────────────────────────────────────────────

def test_vocab_has_five_real_types_and_excludes_other():
    assert set(SIGNALP_VOCAB) == {"SP", "LIPO", "TAT", "TATLIPO", "PILIN"}
    assert OTHER_SENTINEL == "OTHER"
    assert OTHER_SENTINEL not in SIGNALP_VOCAB


def test_prob_columns_fixed_signalp6_order():
    assert PROB_COLUMNS == ["OTHER", "SP", "LIPO", "TAT", "TATLIPO", "PILIN"]


def test_is_kept_skips_sentinel_and_nulls():
    assert is_kept("SP") is True
    assert is_kept("LIPO") is True
    assert is_kept("TATLIPO") is True
    assert is_kept("OTHER") is False
    assert is_kept(None) is False
    assert is_kept("") is False
    assert is_kept("Nonsense") is False


def test_display_name_maps_types():
    assert display_name("SP") == "Signal peptide (Sec/SPI)"
    assert display_name("LIPO") == "Lipoprotein signal peptide (Sec/SPII)"
    assert display_name("PILIN") == "Pilin-like signal peptide (Sec/SPIII)"
    # defensive fallback
    assert display_name("Nonsense") == "Nonsense"


def test_first_token_extracts_wp_accession():
    assert first_token("WP_011131644.1 MULTISPECIES: foo [Prochlorococcus]") == "WP_011131644.1"
    assert first_token("WP_002805124.1") == "WP_002805124.1"
    assert first_token("") == ""
    assert first_token("   ") == ""


# ── cleavage-site parsing ──────────────────────────────────────────────────────

def test_parse_cs_pos_tsv_form():
    assert parse_cs_pos("CS pos: 21-22. Pr: 0.3982") == (21, 0.3982)
    assert parse_cs_pos("CS pos: 23-24. Pr: 0.9964") == (23, 0.9964)


def test_parse_cs_pos_verbose_form():
    # output.json style — the verbose envelope SignalP also writes.
    assert parse_cs_pos("Cleavage site between pos. 26 and 27. Probability 0.862107") == (
        26, 0.862107,
    )


def test_parse_cs_pos_empty_returns_nones():
    assert parse_cs_pos("") == (None, None)
    assert parse_cs_pos("   ") == (None, None)
    assert parse_cs_pos(None) == (None, None)


# ── full prediction_results.txt parsing ────────────────────────────────────────

# Mirrors the real SignalP 6.0 prediction_results.txt layout (verified against
# cache/data/Prochlorococcus/genomes/MED4/signalp/prediction_results.txt).
_RAW = (
    "# SignalP-6.0\tOrganism: Other\tTimestamp: 20260513041522\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tLIPO(Sec/SPII)\tTAT(Tat/SPI)\t"
    "TATLIPO(Tat/SPII)\tPILIN(Sec/SPIII)\tCS Position\n"
    "WP_001.1 MULTISPECIES: photosystem II [Prochlorococcus]\tOTHER\t"
    "0.999987\t0.000002\t0.000002\t0.000000\t0.000000\t0.000007\t\n"
    "WP_002.1 tetratricopeptide repeat protein [Prochlorococcus]\tSP\t"
    "0.452442\t0.542014\t0.004536\t0.000303\t0.000252\t0.000422\tCS pos: 21-22. Pr: 0.3982\n"
    "WP_003.1 Ycf48 [Prochlorococcus]\tLIPO\t"
    "0.000000\t0.000001\t1.000000\t0.000000\t0.000000\t0.000000\tCS pos: 23-24. Pr: 0.9964\n"
)


def test_parse_prediction_results_keys_on_wp():
    recs = parse_prediction_results(_RAW)
    assert set(recs) == {"WP_001.1", "WP_002.1", "WP_003.1"}


def test_parse_prediction_results_other_record():
    recs = parse_prediction_results(_RAW)
    r = recs["WP_001.1"]
    assert r["signalp_type"] == "OTHER"
    assert r["probability"] == pytest.approx(0.999987)
    assert r["cleavage_site"] is None
    assert r["cleavage_probability"] is None


def test_parse_prediction_results_sp_record_uses_winning_class_prob():
    recs = parse_prediction_results(_RAW)
    r = recs["WP_002.1"]
    assert r["signalp_type"] == "SP"
    # probability = the SP column (0.542014), NOT the OTHER column (0.452442).
    assert r["probability"] == pytest.approx(0.542014)
    assert r["cleavage_site"] == 21
    assert r["cleavage_probability"] == pytest.approx(0.3982)


def test_parse_prediction_results_lipo_record():
    recs = parse_prediction_results(_RAW)
    r = recs["WP_003.1"]
    assert r["signalp_type"] == "LIPO"
    assert r["probability"] == pytest.approx(1.0)
    assert r["cleavage_site"] == 23


def test_parse_prediction_results_skips_comments_and_blanks():
    recs = parse_prediction_results("# header\n\n# ID\tPrediction\n   \n")
    assert recs == {}


# ── adapter fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def strain_dir(tmp_path):
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001", "signalp_type": "SP",
                     "signalp_probability": 0.95,
                     "signalp_cleavage_site": 21, "signalp_cleavage_probability": 0.88},
        "PMM_0002": {"locus_tag": "PMM_0002", "signalp_type": "LIPO",
                     "signalp_probability": 0.99,
                     "signalp_cleavage_site": 19, "signalp_cleavage_probability": 0.71},
        "PMM_0003": {"locus_tag": "PMM_0003", "signalp_type": "OTHER",
                     "signalp_probability": 1.0},                          # sentinel
        "PMM_0004": {"locus_tag": "PMM_0004"},                             # no signalp call
        "PMM_0005": {"locus_tag": "PMM_0005", "signalp_type": "TAT",
                     "signalp_probability": 0.80,
                     "signalp_cleavage_site": 30, "signalp_cleavage_probability": 0.62},
    }))
    return d


def _make_orchestrator(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    return MultiSignalPeptideAdapter(genome_config_file=str(config_csv), test_mode=False)


# ── per-strain adapter ────────────────────────────────────────────────────────

def test_per_strain_get_all_types_excludes_other(strain_dir):
    a = SignalPeptideAdapter(genome_dir=strain_dir)
    assert a.get_all_types() == {"SP", "LIPO", "TAT"}


def test_per_strain_edges_skip_other_and_missing(strain_dir):
    a = SignalPeptideAdapter(genome_dir=strain_dir)
    sources = {e[1] for e in a.get_edges()}
    # PMM_0003 (OTHER sentinel) and PMM_0004 (no call) yield no edge.
    assert sources == {"ncbigene:PMM_0001", "ncbigene:PMM_0002", "ncbigene:PMM_0005"}


def test_per_strain_edge_carries_probability_and_cleavage(strain_dir):
    a = SignalPeptideAdapter(genome_dir=strain_dir)
    edges = {e[1]: e for e in a.get_edges()}
    eid, src, tgt, label, props = edges["ncbigene:PMM_0001"]
    assert label == "gene_has_signal_peptide_type"
    assert tgt == _node_id("SP") == "signalp_SP"
    assert props["probability"] == pytest.approx(0.95)
    assert props["cleavage_site"] == 21
    assert props["cleavage_probability"] == pytest.approx(0.88)
    assert isinstance(props["probability"], float)


def test_per_strain_edge_omits_null_cleavage(strain_dir, tmp_path):
    d = tmp_path / "S2"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "G1": {"locus_tag": "G1", "signalp_type": "SP", "signalp_probability": 0.7},
    }))
    a = SignalPeptideAdapter(genome_dir=d)
    _eid, _src, _tgt, _label, props = next(iter(a.get_edges()))
    assert "cleavage_site" not in props
    assert "cleavage_probability" not in props
    assert props["probability"] == pytest.approx(0.7)


# ── orchestrator ──────────────────────────────────────────────────────────────

def test_orchestrator_emits_only_observed_real_types(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    node_ids = {n[0] for n in a.get_nodes()}
    assert node_ids == {"signalp_SP", "signalp_LIPO", "signalp_TAT"}
    assert "signalp_OTHER" not in node_ids
    assert "signalp_TATLIPO" not in node_ids


def test_orchestrator_nodes_are_flat_level_zero(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    nodes = {n[0]: n for n in a.get_nodes()}
    _id, label, props = nodes["signalp_SP"]
    assert label == "signal peptide type"
    assert props["level"] == 0
    assert "level_kind" not in props        # flat ontology
    assert props["signalp_id"] == "SP"
    assert props["name"] == "Signal peptide (Sec/SPI)"


def test_orchestrator_emits_no_hierarchy_edges(strain_dir, tmp_path):
    a = _make_orchestrator(strain_dir, tmp_path)
    labels = {e[3] for e in a.get_edges()}
    assert labels == {"gene_has_signal_peptide_type"}
