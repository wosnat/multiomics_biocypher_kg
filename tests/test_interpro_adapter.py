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
    # R1: house-rule lowercase snake_case, even though InterProScan ships
    # member-DB names UPPERCASE.
    assert props["libraries"] == ["ncbifam", "pfam"]  # sorted distinct
    assert props["match_count"] == 2


def test_aggregate_null_evalue_score_omitted():
    """Profile/pattern hits have no e-value/score → keys omitted, not zeroed."""
    matches = [{"start": 1, "end": 30, "evalue": None, "score": None, "library": "HAMAP"}]
    props = aggregate_match_evidence(matches)
    assert "evalue" not in props
    assert "score" not in props
    assert props["start"] == 1 and props["end"] == 30
    assert props["libraries"] == ["hamap"]
    assert props["match_count"] == 1


def test_aggregate_mixed_null_evalue_uses_present_ones():
    matches = [
        {"start": 1, "end": 10, "evalue": None, "score": 12.0, "library": "PANTHER"},
        {"start": 1, "end": 10, "evalue": 3e-9, "score": None, "library": "PFAM"},
    ]
    props = aggregate_match_evidence(matches)
    assert props["evalue"] == 3e-9
    assert props["score"] == 12.0


def test_aggregate_libraries_match_the_declared_vocabulary():
    """Every library aggregate_match_evidence can emit must be a declared value
    of Gene_has_interpro_entry.libraries (config/controlled_vocabularies.yaml)."""
    from multiomics_kg.utils.controlled_vocab import VOCAB

    raw_libraries = ["CDD", "GENE3D", "HAMAP", "NCBIFAM", "PANTHER", "PFAM",
                      "PIRSF", "PRINTS", "PROSITE_PATTERNS", "PROSITE_PROFILES",
                      "SFLD", "SMART", "SUPERFAMILY"]
    matches = [{"start": 1, "end": 1, "evalue": None, "score": None, "library": lib}
               for lib in raw_libraries]
    props = aggregate_match_evidence(matches)
    for value in props["libraries"]:
        VOCAB.check("Gene_has_interpro_entry", "libraries", value)


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


# ── Layer A: InterproEntry → EC/CAZy router edges (design 2026-08-10) ──────────

def test_observed_ec_node_ids_wellformed():
    a = InterproAnnotationAdapter(MED4_DIR)
    obs = a.observed_ec_node_ids()
    assert obs, "MED4 has EC-annotated genes"
    assert all(nid.startswith("ec:") for nid in obs)


def test_layer_a_related_ec_pruned_and_marked(tmp_path):
    from multiomics_kg.adapters.interpro_adapter import (
        MultiInterproAnnotationAdapter, _ec_node_id,
    )
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")
    m = MultiInterproAnnotationAdapter(
        str(cfg), pfam_node_ids=set(), ec_node_ids=_all_med4_ec_node_ids(),
    )

    strain = m._strain_adapters[0]
    observed = strain.observed_ec_node_ids()
    # pick a real MED4 EC token
    some_ec = next(
        ec for g in strain._genes.values() for ec in (g.get("ec_numbers") or [])
        if _ec_node_id(ec) in observed
    )
    # FAMILY single-EC → not ambiguous; DOMAIN multi-EC (2nd EC unobserved) → ambiguous + pruned
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": [some_ec]},
        "IPR900002": {"type": "DOMAIN", "name": "dom", "ec_numbers": [some_ec, "9.9.9.9"]},
    }
    m._observed_ids = lambda: {"IPR900001", "IPR900002"}

    ec_edges = [e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"]
    by_src = {e[1]: e[4] for e in ec_edges}
    fam_id = _interpro_id("IPR900001")
    dom_id = _interpro_id("IPR900002")
    assert by_src[fam_id]["ambiguous"] is False          # FAMILY + single EC
    assert by_src[fam_id]["source_db"] == "interpro.xml"
    assert by_src[dom_id]["ambiguous"] is True           # DOMAIN / multi-EC
    # the unobserved 9.9.9.9 was pruned (dangling-proof): only one edge per entry
    assert sum(1 for e in ec_edges if e[1] == dom_id) == 1


def _interpro_id(acc):
    from multiomics_kg.adapters.interpro_adapter import _interpro_node_id
    return _interpro_node_id(acc)


def _all_med4_ec_node_ids() -> set[str]:
    """Every EC node id MED4's genes reference — stands in for the real Expasy
    node set in tests that only need 'the endpoint exists'."""
    from multiomics_kg.adapters.interpro_adapter import _ec_node_id
    return InterproAnnotationAdapter(MED4_DIR).observed_ec_node_ids() | {
        _ec_node_id("9.9.9.9")
    }


def test_layer_a_related_ec_pruned_to_existing_ec_nodes(tmp_path):
    """An EC a gene carries but Expasy has NO node for (obsolete InterPro xref,
    e.g. 1.2.8.1) must not produce a Layer-A router edge — it would dangle."""
    from multiomics_kg.adapters.interpro_adapter import (
        MultiInterproAnnotationAdapter, _ec_node_id,
    )
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")

    strain = InterproAnnotationAdapter(MED4_DIR)
    observed = strain.observed_ec_node_ids()
    some_ec = next(
        ec for g in strain._genes.values() for ec in (g.get("ec_numbers") or [])
        if _ec_node_id(ec) in observed
    )
    # The gene set carries `some_ec` AND the obsolete one, but only `some_ec`
    # has an EcNumber node.
    obsolete = "1.2.8.1"
    m = MultiInterproAnnotationAdapter(
        str(cfg), pfam_node_ids=set(), ec_node_ids={_ec_node_id(some_ec)},
    )
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": [some_ec]},
        "IPR900003": {"type": "FAMILY", "name": "obs", "ec_numbers": [obsolete]},
    }
    m._observed_ids = lambda: {"IPR900001", "IPR900003"}
    # Pretend a gene carries the obsolete EC, so `observed_ec` alone would keep it.
    for a in m._strain_adapters:
        a.observed_ec_node_ids = lambda: observed | {_ec_node_id(obsolete)}

    ec_edges = [e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"]
    targets = {e[2] for e in ec_edges}
    assert _ec_node_id(some_ec) in targets
    assert _ec_node_id(obsolete) not in targets


def test_layer_a_related_ec_suppressed_without_ec_node_ids(tmp_path):
    """No EC node set injected → no Layer-A EC edges (mirrors pfam_node_ids)."""
    from multiomics_kg.adapters.interpro_adapter import MultiInterproAnnotationAdapter
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": ["2.7.7.7"]},
    }
    m._observed_ids = lambda: {"IPR900001"}
    assert not [
        e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"
    ]
