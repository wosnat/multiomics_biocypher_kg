# tests/test_tcdb_diamond.py
"""Unit tests for multiomics_kg.utils.tcdb_diamond."""
import json

from multiomics_kg.utils import tcdb_diamond  # noqa: F401
from multiomics_kg.utils.tcdb_diamond import truncate_tcid


def test_truncate_tcid_keeps_first_n_parts():
    assert truncate_tcid("1.A.11.1.5", 3) == "1.A.11"
    assert truncate_tcid("1.A.11.1.5", 4) == "1.A.11.1"
    assert truncate_tcid("1.A.11.1.5", 5) == "1.A.11.1.5"


def test_truncate_tcid_passthrough_when_already_short():
    assert truncate_tcid("1.A.11", 5) == "1.A.11"
    assert truncate_tcid("1.A", 3) == "1.A"


def test_truncate_tcid_empty_input_returns_empty():
    assert truncate_tcid("", 3) == ""


from multiomics_kg.utils.tcdb_diamond import classify_hit


def _hit(identity, qcov, scov, length=200, evalue=1e-10):
    """Build a hit dict with sensible defaults — only override what the test cares about."""
    return {
        "identity": identity, "qcov": qcov, "scov": scov,
        "length": length, "evalue": evalue,
    }


def test_classify_tier_1_high_identity_high_qcov():
    # 80%/85% identity over 80%/75% qcov — tier 1
    assert classify_hit(_hit(80.0, 85.0, 50.0)) == 1
    assert classify_hit(_hit(70.0, 70.0, 50.0)) == 1


def test_classify_tier_2_mid_identity_mid_qcov():
    # 50% identity / 65% qcov — passes tier 2 floor, fails tier 1
    assert classify_hit(_hit(50.0, 65.0, 30.0)) == 2
    assert classify_hit(_hit(40.0, 60.0, 30.0)) == 2


def test_classify_tier_3_gblast3_floor():
    # 35% identity (no tier 1/2), qcov >= 40 — tier 3
    assert classify_hit(_hit(35.0, 45.0, 30.0)) == 3
    # qcov < 40 BUT scov >= 40 — still tier 3 (gblast3 OR rule)
    assert classify_hit(_hit(35.0, 25.0, 50.0)) == 3
    # 25% identity is fine for tier 3 (no identity floor)
    assert classify_hit(_hit(25.0, 45.0, 30.0)) == 3


def test_classify_drops_hit_below_floor():
    # qcov AND scov both < 40 — fails gblast3 OR rule
    assert classify_hit(_hit(35.0, 30.0, 30.0)) is None
    # length < 50 — fails HSP-length floor
    assert classify_hit(_hit(80.0, 90.0, 90.0, length=49)) is None
    # e-value > 0.001 — fails e-value gate
    assert classify_hit(_hit(80.0, 90.0, 90.0, evalue=0.01)) is None


def test_classify_handles_boundary_values_inclusively():
    # All thresholds are >= (not >) — boundary inputs should pass
    assert classify_hit(_hit(70.0, 70.0, 0.0)) == 1
    assert classify_hit(_hit(40.0, 60.0, 0.0)) == 2
    assert classify_hit(_hit(0.0, 40.0, 0.0)) == 3


from multiomics_kg.utils.tcdb_diamond import consensus_collapse


def test_consensus_all_agree_at_5_part():
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.5"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1.5", "agreement": "5_part", "n": 3}


def test_consensus_demote_to_4_part():
    # 4-part agreement, disagreement at position 5
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.7"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1", "agreement": "4_part", "n": 2}


def test_consensus_demote_to_3_part():
    # 3-part agreement only
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.2.3"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11", "agreement": "3_part", "n": 2}


def test_consensus_reject_below_3_part_agreement():
    # Different families
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "2.A.7.4.3"}]
    result = consensus_collapse(hits)
    assert result is None


def test_consensus_single_hit_keeps_5_part():
    hits = [{"tcid": "1.A.11.1.5"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1.5", "agreement": "5_part", "n": 1}


def test_consensus_empty_input_returns_none():
    assert consensus_collapse([]) is None


from multiomics_kg.utils.tcdb_diamond import is_class_9


def test_is_class_9_matches_top_class():
    assert is_class_9("9.B.82.1.5") is True
    assert is_class_9("9.A.1") is True
    assert is_class_9("9") is True


def test_is_class_9_excludes_other_classes():
    assert is_class_9("1.A.11.1.5") is False
    assert is_class_9("8.A.1") is False  # auxiliary, not incompletely characterized
    assert is_class_9("19.A.1") is False  # not real but tests prefix matching


def test_is_class_9_handles_empty():
    assert is_class_9("") is False


from multiomics_kg.utils.tcdb_diamond import parse_diamond_row


def test_parse_diamond_row_extracts_typed_fields():
    line = "WP_011131900.1\tlcl|Q9I3F6-1.A.11.1.5\t87.4\t92.1\t89.7\t412\t1.2e-180\t650.5"
    row = parse_diamond_row(line)
    assert row["query_id"] == "WP_011131900.1"
    assert row["subject_id"] == "lcl|Q9I3F6-1.A.11.1.5"
    assert row["identity"] == 87.4
    assert row["qcov"] == 92.1
    assert row["scov"] == 89.7
    assert row["length"] == 412
    assert row["evalue"] == 1.2e-180
    assert row["bitscore"] == 650.5


def test_parse_diamond_row_returns_none_for_short_line():
    assert parse_diamond_row("only\ttwo\tcolumns") is None
    assert parse_diamond_row("") is None


def test_parse_diamond_row_returns_none_for_invalid_numeric():
    line = "WP_111.1\tlcl|X-1.A\tNOT_A_NUMBER\t90\t90\t100\t1e-10\t300"
    assert parse_diamond_row(line) is None


from multiomics_kg.utils.tcdb_diamond import parse_tcdb_subject_id


def test_parse_tcdb_subject_id_lcl_prefix():
    assert parse_tcdb_subject_id("lcl|Q9I3F6-1.A.11.1.5") == ("Q9I3F6", "1.A.11.1.5")


def test_parse_tcdb_subject_id_no_lcl_prefix():
    # If header was written without lcl|, accept it
    assert parse_tcdb_subject_id("Q9I3F6-1.A.11.1.5") == ("Q9I3F6", "1.A.11.1.5")


def test_parse_tcdb_subject_id_handles_dashes_in_accession():
    # UniProt accessions can contain dashes for isoforms — split on the LAST dash
    # before a dotted TCID
    result = parse_tcdb_subject_id("lcl|P12345-2-1.A.11.1.5")
    assert result == ("P12345-2", "1.A.11.1.5")


def test_parse_tcdb_subject_id_returns_none_when_no_tcid():
    assert parse_tcdb_subject_id("lcl|Q9I3F6") is None
    assert parse_tcdb_subject_id("") is None


def test_parse_tcdb_subject_id_validates_tcid_shape():
    # TCID must be at least 3 dotted parts to be plausible
    assert parse_tcdb_subject_id("lcl|Q9I3F6-1.A") is None


from multiomics_kg.utils.tcdb_diamond import build_strain_calls


def test_build_strain_calls_full_pipeline(tmp_path):
    tsv_content = (
        # 5 strong identical hits -> tier 1 5-part call
        "WP_AAA.1\tlcl|Q1-1.A.11.1.5\t87.4\t92.1\t89.7\t412\t1e-180\t650\n"
        "WP_AAA.1\tlcl|Q2-1.A.11.1.5\t86.0\t91.0\t88.0\t410\t1e-179\t640\n"
        "WP_AAA.1\tlcl|Q3-1.A.11.1.5\t85.0\t90.0\t87.0\t408\t1e-178\t630\n"
        # Hits scattered across families -> reject (consensus < 3-part)
        "WP_BBB.1\tlcl|Q4-1.A.11.1.5\t75.0\t75.0\t60.0\t300\t1e-100\t450\n"
        "WP_BBB.1\tlcl|Q5-2.A.7.4.3\t72.0\t73.0\t60.0\t300\t1e-95\t440\n"
        # Single hit, mid identity -> tier 2 with 4-part TCID
        "WP_CCC.1\tlcl|Q6-1.A.11.1.5\t50.0\t65.0\t40.0\t250\t1e-50\t300\n"
        # Single hit, low identity passing tier 3 floor -> tier 3 with 3-part
        "WP_DDD.1\tlcl|Q7-1.A.11.1.5\t30.0\t45.0\t30.0\t150\t1e-10\t120\n"
        # Class 9 hit -> tagged
        "WP_EEE.1\tlcl|Q8-9.B.82.1.5\t80.0\t85.0\t75.0\t300\t1e-150\t500\n"
        # Hit failing the floor -> dropped
        "WP_FFF.1\tlcl|Q9-1.A.11.1.5\t30.0\t30.0\t30.0\t150\t1e-10\t100\n"
    )
    tsv = tmp_path / "test.tcdb.tsv"
    tsv.write_text(tsv_content)

    calls, summary = build_strain_calls(tsv)

    def call0(pid):
        """Return the first (highest-confidence) candidate for a protein."""
        return calls[pid]["calls"][0]

    # WP_AAA.1: 3 strong consensus hits at 5-part -> tier 1 specificity
    c = call0("WP_AAA.1")
    assert c["tcid"] == "1.A.11.1.5"
    assert c["tier"] == 1
    assert c["level_kind"] == "tc_specificity"
    assert c["consensus_agreement"] == "5_part"
    assert c["consensus_n"] == 3
    assert c["incompletely_characterized"] is False
    # Single-family protein -> one candidate
    assert len(calls["WP_AAA.1"]["calls"]) == 1

    # WP_BBB.1: scattered across 2 families -> NOW recovers BOTH as candidates
    # (was rejected by global consensus_collapse in the old single-call design)
    assert "WP_BBB.1" in calls
    bbb_calls = calls["WP_BBB.1"]["calls"]
    assert len(bbb_calls) == 2
    bbb_tcids = {c["tcid"] for c in bbb_calls}
    assert bbb_tcids == {"1.A.11.1.5", "2.A.7.4.3"}
    # First candidate is the higher-confidence one (75/75 -> 0.5625 > 72/73 -> 0.5256)
    assert bbb_calls[0]["confidence_score"] > bbb_calls[1]["confidence_score"]

    # WP_CCC.1: single tier-2 hit -> 4-part subfamily
    c = call0("WP_CCC.1")
    assert c["tcid"] == "1.A.11.1"
    assert c["tier"] == 2
    assert c["level_kind"] == "tc_subfamily"

    # WP_DDD.1: tier 3 -> 3-part family
    c = call0("WP_DDD.1")
    assert c["tcid"] == "1.A.11"
    assert c["tier"] == 3
    assert c["level_kind"] == "tc_family"

    # WP_EEE.1: class 9 -> tagged; identity is high (80%) so still tier 1 — no demotion
    # despite class 9 (spec §6.4-C: "No demotion — let merge / downstream consumers decide")
    c = call0("WP_EEE.1")
    assert c["tcid"] == "9.B.82.1.5"
    assert c["incompletely_characterized"] is True
    assert c["tier"] == 1
    assert c["level_kind"] == "tc_specificity"

    # WP_FFF.1: floor failure -> not in calls
    assert "WP_FFF.1" not in calls

    # confidence_score on the WINNING candidate per protein
    assert abs(call0("WP_AAA.1")["confidence_score"] - 0.8050) < 1e-3
    assert abs(call0("WP_CCC.1")["confidence_score"] - 0.3250) < 1e-3
    assert abs(call0("WP_DDD.1")["confidence_score"] - 0.1350) < 1e-3
    assert abs(call0("WP_EEE.1")["confidence_score"] - 0.6800) < 1e-3

    # Summary — counts CANDIDATES, not proteins
    assert summary["raw_hit_lines"] == 9
    assert summary["proteins_with_call"] == 5      # was 4 (WP_BBB.1 recovered)
    assert summary["total_candidates"] == 6        # 4 single-family + 2 from WP_BBB.1
    assert "proteins_rejected_by_consensus" not in summary  # field removed
    # WP_BBB.1 contributes 2 tier-1 candidates (75/75 + 72/73 both clear T1 floor)
    assert summary["tier_distribution"] == {"1": 4, "2": 1, "3": 1}
    # 1-candidate proteins: AAA, CCC, DDD, EEE = 4; 2-candidate: BBB = 1
    assert summary["candidates_per_protein_distribution"] == {"1": 4, "2": 1}


def test_build_strain_calls_mixed_tier_hits_use_best_not_worst(tmp_path):
    """When top-N hits agree at consensus depth but vary in identity, the
    BEST hit's tier drives effective_tier (not the worst). Coupled with the
    consensus-depth floor via max(...), this gives the more informative call
    when strong evidence exists for the family but the top-N includes a weak
    homolog at the floor.
    """
    # 3 hits all agree at 4-part 1.A.11.1 (disagree at 5-part):
    #   - hit 1: identity 75%, qcov 75% -> tier 1 (>=70/>=70)
    #   - hit 2: identity 50%, qcov 65% -> tier 2
    #   - hit 3: identity 30%, qcov 45% -> tier 3 (floor only)
    # Old rule (worst_tier): max(3, 2) -> tier 3 -> 3-part 1.A.11
    # New rule (best_tier):  max(1, 2) -> tier 2 -> 4-part 1.A.11.1
    tsv_content = (
        "WP_GGG.1\tlcl|Q10-1.A.11.1.5\t75.0\t75.0\t60.0\t300\t1e-100\t450\n"
        "WP_GGG.1\tlcl|Q11-1.A.11.1.7\t50.0\t65.0\t60.0\t300\t1e-50\t300\n"
        "WP_GGG.1\tlcl|Q12-1.A.11.1.9\t30.0\t45.0\t30.0\t150\t1e-10\t120\n"
    )
    tsv = tmp_path / "test.tcdb.tsv"
    tsv.write_text(tsv_content)
    calls, _ = build_strain_calls(tsv)
    c = calls["WP_GGG.1"]["calls"][0]
    assert c["tier"] == 2
    assert c["tcid"] == "1.A.11.1"
    assert c["level_kind"] == "tc_subfamily"
    # Best hit (75%/75%) drives metadata + score
    assert c["identity"] == 75.0
    assert c["qcov"] == 75.0
    # confidence: 0.75 * 0.75 * 0.85 (4_part) = 0.4781
    assert abs(c["confidence_score"] - 0.4781) < 1e-3


def test_build_strain_calls_emits_multi_family_candidates(tmp_path):
    """The headline multi-call recovery: a protein with diamond hits in two
    distinct 3-part families (e.g. RND + MFP partner-protein confusion) gets
    BOTH as candidates instead of being rejected by global consensus.
    """
    tsv_content = (
        # 2 strong hits in family 2.A.6 (RND)
        "WP_RND.1\tlcl|Q1-2.A.6.1.4\t100.0\t100.0\t100.0\t400\t1e-200\t800\n"
        "WP_RND.1\tlcl|Q2-2.A.6.1.5\t98.0\t99.0\t99.0\t398\t1e-198\t790\n"
        # 1 weaker hit in family 8.A.1 (MFP)
        "WP_RND.1\tlcl|Q3-8.A.1.2.1\t60.0\t80.0\t80.0\t300\t1e-100\t450\n"
    )
    tsv = tmp_path / "x.tcdb.tsv"
    tsv.write_text(tsv_content)

    calls, summary = build_strain_calls(tsv)
    rec = calls["WP_RND.1"]
    assert len(rec["calls"]) == 2
    # First candidate is the higher-confidence one (RND, 100/100)
    assert rec["calls"][0]["tcid"].startswith("2.A.6")
    assert rec["calls"][0]["confidence_score"] > rec["calls"][1]["confidence_score"]
    # Both families recovered
    assert {c["tcid"][:5] for c in rec["calls"]} == {"2.A.6", "8.A.1"}
    # Summary reflects per-candidate counts
    assert summary["total_candidates"] == 2
    assert summary["proteins_with_call"] == 1
    assert summary["candidates_per_protein_distribution"] == {"2": 1}


# ============================================================================
# calls.json consumption helpers
# ============================================================================

from multiomics_kg.utils.tcdb_diamond import (
    load_calls_json,
    iter_candidates,
    best_call,
    call_tc_families,
    call_tcids,
)


def _cand(tcid, score, **extra):
    return {"tcid": tcid, "confidence_score": score, **extra}


def _rec(*candidates):
    """Per-protein record — the on-disk shape is just {"calls": [...]}."""
    return {"calls": list(candidates)}


def test_load_calls_json_reads_on_disk_shape(tmp_path):
    p = tmp_path / "x.tcdb.calls.json"
    payload = {"WP_X.1": _rec(_cand("1.A.11.1.5", 0.8))}
    p.write_text(json.dumps(payload))
    assert load_calls_json(p) == payload


def test_load_calls_json_missing_file_raises():
    import pytest
    with pytest.raises(FileNotFoundError):
        load_calls_json("/nonexistent/path/calls.json")


def test_iter_candidates_yields_every_candidate():
    """Nothing is pre-filtered — the tier policy already gated quality at build
    time, so consumers see the full candidate set and threshold explicitly.
    """
    calls = {
        "WP_A.1": _rec(_cand("1.A.11", 0.5), _cand("2.A.6", 0.1)),
        "WP_B.1": _rec(_cand("8.A.1", 0.05)),
    }
    assert len(list(iter_candidates(calls))) == 3


def test_iter_candidates_preserves_per_protein_order():
    rec = _rec(_cand("1.A.11.1.5", 0.9), _cand("2.A.6.1", 0.4), _cand("8.A.1", 0.2))
    seen = [c["tcid"] for _pid, c in iter_candidates({"WP_X.1": rec})]
    assert seen == ["1.A.11.1.5", "2.A.6.1", "8.A.1"]


def test_best_call_returns_highest_confidence():
    rec = _rec(_cand("2.A.6.1", 0.85), _cand("8.A.1.2", 0.55))
    assert best_call(rec)["tcid"] == "2.A.6.1"


def test_best_call_none_on_empty():
    assert best_call({"calls": []}) is None
    assert best_call({}) is None  # missing calls key


def test_call_tc_families_collapses_to_3_part_and_dedupes():
    rec = _rec(_cand("1.A.11.1.5", 0.9), _cand("1.A.11.2", 0.5), _cand("2.A.6.1", 0.4))
    assert call_tc_families(rec) == ["1.A.11", "2.A.6"]


def test_call_tc_families_empty_record():
    assert call_tc_families({"calls": []}) == []


def test_call_tcids_preserves_depth_and_order():
    rec = _rec(_cand("1.A.11.1.5", 0.9), _cand("2.A.6.1", 0.4), _cand("8.A.1", 0.2))
    assert call_tcids(rec) == ["1.A.11.1.5", "2.A.6.1", "8.A.1"]


# ============================================================================
# Regression: the artifact carries no downstream-derived state
# ============================================================================


def test_calls_carry_only_sequence_derived_fields(tmp_path):
    """The runner must not reach into eggNOG or gene_annotations_merged.json.

    Guards the cycle fixed in the 2026-08-06 spec (§3.2): the old version read
    step 2's OUTPUT to compute agreement tags, while step 2 is where these
    calls are destined to be merged.
    """
    tsv = tmp_path / "x.tcdb.tsv"
    tsv.write_text("WP_A.1\tlcl|Q1-2.A.6.1.4\t100.0\t100.0\t100.0\t400\t1e-200\t800\n")

    calls, summary = build_strain_calls(tsv)

    rec = calls["WP_A.1"]
    # Protein-level: nothing but `calls`
    assert set(rec) == {"calls"}
    # Candidate-level: no cross-source or filter state
    cand = rec["calls"][0]
    forbidden = {"egn_agreement", "egn_tcids", "pfam_agreement", "pfam_ids",
                 "pfam_tc_families", "filter_action"}
    assert forbidden.isdisjoint(cand)
    assert set(cand) == {
        "tcid", "level_kind", "tier", "confidence_score", "identity", "qcov",
        "scov", "evalue", "length", "consensus_n", "consensus_agreement",
        "incompletely_characterized",
    }
    assert forbidden.isdisjoint(summary)


def test_build_strain_calls_takes_only_the_tsv():
    """Signature guard: extra positional args would reintroduce the cycle."""
    import inspect
    params = list(inspect.signature(build_strain_calls).parameters)
    assert params == ["tsv_path"]


def test_verdict_is_independent_of_unrelated_siblings(tmp_path):
    """A candidate's fields must not change because the protein happens to have
    an unrelated second-domain hit.

    The deleted filter chain violated this: 1,159 candidates were kept purely
    for lacking a sibling, and 1,452 were dropped by a sibling that was itself
    dropped (spec §3.1).
    """
    hit_a = "WP_X.1\tlcl|Q1-2.A.6.1.4\t100.0\t100.0\t100.0\t400\t1e-200\t800\n"
    hit_b = "WP_X.1\tlcl|Q3-8.A.1.2.1\t60.0\t80.0\t80.0\t300\t1e-100\t450\n"

    alone = tmp_path / "alone.tsv"
    alone.write_text(hit_a)
    with_sibling = tmp_path / "sibling.tsv"
    with_sibling.write_text(hit_a + hit_b)

    solo, _ = build_strain_calls(alone)
    multi, _ = build_strain_calls(with_sibling)

    solo_cand = solo["WP_X.1"]["calls"][0]
    multi_cand = next(c for c in multi["WP_X.1"]["calls"] if c["tcid"].startswith("2.A.6"))
    assert solo_cand == multi_cand


from multiomics_kg.utils.tcdb_diamond import confidence_score


# ── consensus depth can never exceed the shallowest hit ──────────────────────
#
# `parse_tcdb_subject_id` accepts 3-5 part TCIDs, and list slicing does not pad,
# so a group of 4-part hits used to match at depth 5 and report "5_part" —
# inflating the agreement weight to 1.0 and labelling a subfamily-depth call as
# tc_specificity. TCDB ships only 5-part headers today (40,520/40,520 candidates
# across 42 strains), so this guards an invariant rather than a live bug.


def test_consensus_never_claims_more_depth_than_the_hits_have():
    hits = [{"tcid": "1.A.11.1"}, {"tcid": "1.A.11.1"}]
    out = consensus_collapse(hits)
    assert out["agreement"] == "4_part"
    assert out["tcid"] == "1.A.11.1"


def test_consensus_depth_capped_by_the_shallowest_hit():
    """A 3-part hit alongside 5-part ones caps the consensus at 3 parts."""
    hits = [{"tcid": "2.A.1.2.9"}, {"tcid": "2.A.1"}]
    out = consensus_collapse(hits)
    assert out["agreement"] == "3_part"
    assert out["tcid"] == "2.A.1"


def test_five_part_hits_still_reach_full_depth():
    """The guard must not cost depth in the normal all-5-part case."""
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.5"}]
    out = consensus_collapse(hits)
    assert out["agreement"] == "5_part"
    assert out["tcid"] == "1.A.11.1.5"


def test_agreement_weight_is_defined_for_every_returned_depth():
    """confidence_score() indexes _AGREEMENT_WEIGHT by the agreement string — a
    depth the table lacks would raise KeyError at build time."""
    for tcids, expected in (
        (["1.A.11.1.5", "1.A.11.1.5"], "5_part"),
        (["1.A.11.1", "1.A.11.1"], "4_part"),
        (["1.A.11.1.5", "1.A.11.2.7"], "3_part"),
    ):
        out = consensus_collapse([{"tcid": t} for t in tcids])
        assert out["agreement"] == expected
        assert confidence_score(80.0, 80.0, out["agreement"]) > 0
