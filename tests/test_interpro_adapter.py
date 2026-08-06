"""Unit tests for the InterPro adapter (pure evidence aggregation + ancestor walk)
plus a fast real-cache sanity check on MED4."""

from __future__ import annotations

from pathlib import Path

import pytest

from multiomics_kg.adapters.interpro_adapter import (
    InterproAnnotationAdapter,
    aggregate_match_evidence,
    kept_ids,
)

MED4_DIR = Path("cache/data/Prochlorococcus/genomes/MED4")


# ── pure: evidence aggregation ────────────────────────────────────────────────

def test_aggregate_envelope_and_best_scores():
    matches = [
        {"start": 10, "end": 100, "evalue": 1e-20, "score": 50.0, "library": "PFAM"},
        {"start": 5, "end": 80, "evalue": 1e-30, "score": 70.0, "library": "NCBIFAM"},
    ]
    props = aggregate_match_evidence(matches)
    assert props["start"] == 5           # min
    assert props["end"] == 100           # max
    assert props["evalue"] == 1e-30      # best (min)
    assert props["score"] == 70.0        # best (max)
    assert props["libraries"] == ["NCBIFAM", "PFAM"]  # sorted distinct
    assert props["match_count"] == 2


def test_aggregate_null_evalue_score_omitted():
    """Profile/pattern hits have no e-value/score → keys omitted, not zeroed."""
    matches = [{"start": 1, "end": 30, "evalue": None, "score": None, "library": "HAMAP"}]
    props = aggregate_match_evidence(matches)
    assert "evalue" not in props
    assert "score" not in props
    assert props["start"] == 1 and props["end"] == 30
    assert props["libraries"] == ["HAMAP"]
    assert props["match_count"] == 1


def test_aggregate_mixed_null_evalue_uses_present_ones():
    matches = [
        {"start": 1, "end": 10, "evalue": None, "score": 12.0, "library": "PANTHER"},
        {"start": 1, "end": 10, "evalue": 3e-9, "score": None, "library": "PFAM"},
    ]
    props = aggregate_match_evidence(matches)
    assert props["evalue"] == 3e-9
    assert props["score"] == 12.0


# ── pure: ancestor pruning ────────────────────────────────────────────────────

REF = {
    "IPR000004": {"parent": "IPR000002", "level": 2},
    "IPR000002": {"parent": "IPR000001", "level": 1},
    "IPR000001": {"parent": None, "level": 0},
    "IPR000685": {"parent": None, "level": 0},
}


def test_kept_ids_pulls_full_ancestor_chain():
    kept = kept_ids({"IPR000004"}, REF)
    assert kept == {"IPR000004", "IPR000002", "IPR000001"}  # + ancestors


def test_kept_ids_entry_absent_from_reference_kept_as_is():
    kept = kept_ids({"IPR999999"}, REF)
    assert kept == {"IPR999999"}  # no ancestors, but itself kept (node still emitted)


def test_kept_ids_cycle_safe():
    cyclic = {"A": {"parent": "B"}, "B": {"parent": "A"}}
    assert kept_ids({"A"}, cyclic) == {"A", "B"}  # terminates


# ── fast real-cache sanity (MED4) ─────────────────────────────────────────────

@pytest.mark.skipif(not MED4_DIR.exists(), reason="MED4 cache not present")
def test_med4_edges_wellformed():
    a = InterproAnnotationAdapter(MED4_DIR)
    edges = list(a.get_edges())
    assert len(edges) > 1000
    for eid, src, tgt, label, props in edges[:500]:
        assert label == "gene_has_interpro_entry"
        assert src.startswith("ncbigene:")
        assert tgt.startswith("interpro:IPR")
        assert props["match_count"] >= 1
        assert props["libraries"]                       # non-empty
        if "start" in props and "end" in props:
            assert props["start"] <= props["end"]        # valid envelope
        if "evalue" in props:
            assert props["evalue"] >= 0
