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


from multiomics_kg.utils.tcdb_diamond import compute_egn_agreement


def test_egn_agreement_confirms_identical():
    assert compute_egn_agreement("1.A.11", ["1.A.11"]) == "confirms"
    assert compute_egn_agreement("1.A.11.1.5", ["1.A.11.1.5"]) == "confirms"


def test_egn_agreement_confirms_diamond_descendant():
    # eggNOG family-level, diamond a strict descendant — still confirms (consistent)
    assert compute_egn_agreement("1.A.11.1.5", ["1.A.11"]) == "refines"
    assert compute_egn_agreement("1.A.11.1", ["1.A.11"]) == "refines"


def test_egn_agreement_extends_when_eggnog_missing():
    assert compute_egn_agreement("1.A.11.1.5", None) == "extends"
    assert compute_egn_agreement("1.A.11.1.5", []) == "extends"
    assert compute_egn_agreement("1.A.11.1.5", "") == "extends"


def test_egn_agreement_conflicts_different_family():
    assert compute_egn_agreement("1.A.11.1.5", ["2.A.7"]) == "conflicts"
    assert compute_egn_agreement("1.A.11", ["1.B.5"]) == "conflicts"


def test_egn_agreement_eggnog_descendant_of_diamond_is_conflict():
    # Diamond at 3-part, eggNOG at 5-part below it. Same family — not a conflict;
    # this case is rare (eggNOG rarely emits 5-part) but treat as confirms.
    assert compute_egn_agreement("1.A.11", ["1.A.11.1.5"]) == "confirms"


def test_egn_agreement_accepts_legacy_string_form():
    # Backward compat: callers passing a single str (not list) still work
    assert compute_egn_agreement("1.A.11", "1.A.11") == "confirms"
    assert compute_egn_agreement("1.A.11", "2.A.7") == "conflicts"


def test_egn_agreement_multi_value_picks_best_match():
    # MreB/MreBCD case: eggNOG carries both 1.A.33.1 (legacy Hsp70-channel call)
    # and 9.B.157.1 (correct MreBCD-family call) for rod-shape-determining proteins.
    # Diamond's 9.B.157.1 must be tagged `confirms`, not `conflicts`.
    assert compute_egn_agreement("9.B.157.1", ["1.A.33.1", "9.B.157.1"]) == "confirms"
    # Order shouldn't matter
    assert compute_egn_agreement("9.B.157.1", ["9.B.157.1", "1.A.33.1"]) == "confirms"


def test_egn_agreement_multi_value_refines_when_one_is_ancestor():
    # eggNOG has both an unrelated TC and a family-level ancestor of diamond's call —
    # `refines` wins over `conflicts` because it's the stronger match.
    assert compute_egn_agreement("1.A.11.1.5", ["2.A.7", "1.A.11"]) == "refines"


def test_egn_agreement_multi_value_all_conflict():
    # Every eggNOG value disagrees at family level — still conflicts
    assert compute_egn_agreement("1.A.11.1.5", ["2.A.7", "3.A.1"]) == "conflicts"


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


import textwrap
from multiomics_kg.utils.tcdb_diamond import load_eggnog_kegg_tc


def test_load_eggnog_kegg_tc_extracts_per_protein_value(tmp_path):
    # Real eggNOG-mapper v2.1 column order, 21 columns
    content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_011131852.1\t59919.PMM0213\t2.76e-223\t617.0\tCOG3329@1|root\t1117\tS\tdesc\tsbtA\t-\t-\tko:K07086\t-\t-\t-\t-\tko00000\t1.A.11\t-\t-\tSbt_1
        WP_011131900.1\t59919.PMM0263\t0.0\t934.0\tCOG0004@1|root\t1117\tP\tAmm\tamtB\t-\t-\tko:K03320\t-\t-\t-\t-\tko00000,ko02000\t1.A.11\t-\t-\tAmmonium_transp
        WP_NOTC.1\t59919.PMM0001\t0.0\t100.0\tCOGZZZ\t1117\tS\tdesc\tx\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-
        ## end of file
        """)
    path = tmp_path / "MED4.emapper.annotations"
    path.write_text(content)

    result = load_eggnog_kegg_tc(path)
    assert result == {
        "WP_011131852.1": ["1.A.11"],
        "WP_011131900.1": ["1.A.11"],
        # WP_NOTC.1 has KEGG_TC = "-" -> not in dict
    }


def test_load_eggnog_kegg_tc_handles_missing_file(tmp_path):
    assert load_eggnog_kegg_tc(tmp_path / "no_such_file") == {}


def test_load_eggnog_kegg_tc_skips_empty_or_dash_values(tmp_path):
    content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_A.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t-\t-\t-\t-
        WP_B.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t\t-\t-\t-
        WP_C.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.11,3.A.1.27\t-\t-\t-
        """)
    path = tmp_path / "x.emapper.annotations"
    path.write_text(content)
    result = load_eggnog_kegg_tc(path)
    # Multi-value KEGG_TC: keep ALL values; compute_egn_agreement picks the best match
    assert result == {"WP_C.1": ["1.A.11", "3.A.1.27"]}


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

    egn_content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_AAA.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.11\t-\t-\t-
        WP_DDD.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t9.A.1\t-\t-\t-
        """)
    egn = tmp_path / "test.emapper.annotations"
    egn.write_text(egn_content)

    calls, summary = build_strain_calls(tsv, egn)

    def call0(pid):
        """Return the first (highest-confidence) candidate for a protein."""
        return calls[pid]["calls"][0]

    # WP_AAA.1: 3 strong consensus hits at 5-part -> tier 1 specificity, refines eggNOG
    c = call0("WP_AAA.1")
    assert c["tcid"] == "1.A.11.1.5"
    assert c["tier"] == 1
    assert c["level_kind"] == "tc_specificity"
    assert c["consensus_agreement"] == "5_part"
    assert c["consensus_n"] == 3
    assert c["egn_agreement"] == "refines"
    assert calls["WP_AAA.1"]["egn_tcids"] == ["1.A.11"]
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
    assert c["egn_agreement"] == "extends"

    # WP_DDD.1: tier 3 -> 3-part family; eggNOG conflict (different family)
    c = call0("WP_DDD.1")
    assert c["tcid"] == "1.A.11"
    assert c["tier"] == 3
    assert c["level_kind"] == "tc_family"
    assert c["egn_agreement"] == "conflicts"
    assert calls["WP_DDD.1"]["egn_tcids"] == ["9.A.1"]

    # WP_EEE.1: class 9 -> tagged; identity is high (80%) so still tier 1 — no demotion
    # despite class 9 (spec §6.4-C: "No demotion — let merge / downstream consumers decide")
    c = call0("WP_EEE.1")
    assert c["tcid"] == "9.B.82.1.5"
    assert c["incompletely_characterized"] is True
    assert c["egn_agreement"] == "extends"
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
    assert summary["agreement_distribution"]["refines"] == 1
    assert summary["agreement_distribution"]["extends"] == 4   # CCC, EEE, + both BBB candidates (no eggNOG)
    assert summary["agreement_distribution"]["conflicts"] == 1
    assert summary["agreement_distribution"]["confirms"] == 0
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
    egn = tmp_path / "missing.emapper.annotations"  # no eggNOG file

    calls, _ = build_strain_calls(tsv, egn)
    c = calls["WP_GGG.1"]["calls"][0]
    assert c["tier"] == 2
    assert c["tcid"] == "1.A.11.1"
    assert c["level_kind"] == "tc_subfamily"
    # Best hit (75%/75%) drives metadata + score
    assert c["identity"] == 75.0
    assert c["qcov"] == 75.0
    # confidence: 0.75 * 0.75 * 0.85 (4_part) = 0.4781
    assert abs(c["confidence_score"] - 0.4781) < 1e-3


from multiomics_kg.utils.tcdb_diamond import (
    load_pfam_to_tc_map,
    load_gene_pfams,
    compute_pfam_agreement,
    gene_pfam_implied_families,
)


def test_load_pfam_to_tc_map_parses_3col_tsv(tmp_path):
    p = tmp_path / "pfam_map.tsv"
    p.write_text(
        "PF07885\t1.A.1.1.1\tThe Voltage-gated Ion Channel (VIC) Superfamily\n"
        "PF07885\t1.A.1.29.1\tThe Voltage-gated Ion Channel (VIC) Superfamily\n"
        "PF00520\t1.A.1.2.2\tThe Voltage-gated Ion Channel (VIC) Superfamily\n"
        "PF02462\t1.B.6.2.7\tThe OmpA-OmpF Porin (OOP) Family\n"
    )
    m = load_pfam_to_tc_map(p)
    # Both PF07885 5-part hits collapse to the same 3-part family
    assert m["PF07885"] == {"1.A.1"}
    assert m["PF00520"] == {"1.A.1"}
    assert m["PF02462"] == {"1.B.6"}


def test_load_pfam_to_tc_map_skips_garbage_rows(tmp_path):
    p = tmp_path / "pfam_map.tsv"
    p.write_text(
        "PF07885\t1.A.1.1.1\tValid row\n"
        "header line\tnot tabbed correctly\n"
        "TesT\t3.A.1.1.1\tnon-PF row that the source page sometimes has\n"
        "PF12345\t1.A\ttoo-short TC id\n"
    )
    m = load_pfam_to_tc_map(p)
    assert m == {"PF07885": {"1.A.1"}}


def test_load_pfam_to_tc_map_missing_file_returns_empty(tmp_path):
    assert load_pfam_to_tc_map(tmp_path / "nope.tsv") == {}


def test_load_gene_pfams_extracts_protein_id_to_pfams(tmp_path):
    p = tmp_path / "gene_annotations_merged.json"
    p.write_text(json.dumps({
        "MIT1002_00220": {"protein_id": "WP_014947909.1", "pfam_ids": ["PF06723"]},
        "MIT1002_00500": {"protein_id_refseq": "WP_222.1", "pfam_ids": ["PF00001", "PF00002"]},
        "MIT1002_NOPID": {"locus_tag": "x"},  # no protein_id -> skipped
        "MIT1002_NOPFAM": {"protein_id": "WP_333.1", "pfam_ids": []},
    }))
    result = load_gene_pfams(p)
    assert result == {
        "WP_014947909.1": ["PF06723"],
        "WP_222.1": ["PF00001", "PF00002"],
    }


def test_gene_pfam_implied_families_unions_and_sorts():
    pfam_map = {"PF1": {"2.A.6"}, "PF2": {"8.A.1", "1.B.17"}, "PF3": {"2.A.6"}}
    assert gene_pfam_implied_families(["PF1", "PF2", "PF3"], pfam_map) == ["1.B.17", "2.A.6", "8.A.1"]


def test_gene_pfam_implied_families_empty_when_no_pfams_or_no_map_hits():
    assert gene_pfam_implied_families([], {"PF1": {"2.A.6"}}) == []
    assert gene_pfam_implied_families(["PF999"], {"PF1": {"2.A.6"}}) == []


def test_compute_pfam_agreement_neutral_when_implied_empty():
    assert compute_pfam_agreement("1.A.11", ["2.A.7"], []) == "neutral"


def test_compute_pfam_agreement_confirms_diamond():
    # Pfam-implied family matches diamond's, but no eggNOG family
    assert compute_pfam_agreement("9.B.157.1", ["1.A.33.1"], ["9.B.157"]) == "confirms_diamond"


def test_compute_pfam_agreement_confirms_eggnog():
    # Pfam supports an eggNOG family, not diamond's
    assert compute_pfam_agreement("2.A.6.1", ["8.A.1.2.1"], ["8.A.1"]) == "confirms_eggnog"


def test_compute_pfam_agreement_confirms_both():
    # Pfam covers BOTH diamond's family AND an eggNOG family
    assert compute_pfam_agreement("2.A.6.1", ["8.A.1.2.1"], ["2.A.6", "8.A.1"]) == "confirms_both"


def test_compute_pfam_agreement_contradicts_both():
    # Pfam implies a TC family, but neither diamond nor eggNOG hit it
    assert compute_pfam_agreement("1.A.11", ["2.A.7"], ["5.A.1"]) == "contradicts_both"


def test_build_strain_calls_pfam_aware_pipeline(tmp_path):
    """End-to-end Pfam path: a single MreB protein with PF06723 (MreB_Mbl)
    should get pfam_agreement = 'confirms_diamond' when diamond's 9.B.157 call
    matches the Pfam mapping. Validates the new schema fields and the new
    pfam_agreement_distribution summary key.
    """
    tsv = tmp_path / "mreb.tsv"
    tsv.write_text(
        "WP_MREB.1\tlcl|Q1-9.B.157.1\t89.9\t100.0\t95.0\t340\t1e-180\t650\n"
        "WP_MREB.1\tlcl|Q2-9.B.157.1\t85.0\t99.0\t94.0\t338\t1e-170\t620\n"
    )
    egn = tmp_path / "mreb.emapper.annotations"
    import textwrap
    egn.write_text(textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_MREB.1\tx\tx\tx\tx\tx\tD\tx\tmreB\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.33.1,9.B.157.1\t-\t-\tMreB_Mbl
        """))
    gene_ann = tmp_path / "gene_annotations_merged.json"
    gene_ann.write_text(json.dumps({
        "MIT1002_00220": {"protein_id": "WP_MREB.1", "pfam_ids": ["PF06723"]},
    }))
    pfam_map = tmp_path / "tcdb_pfam_map.tsv"
    pfam_map.write_text(
        "PF06723\t9.B.157.1.1\tThe Cell Shape-determining MreBCD (MreBCD) Family\n"
    )

    calls, summary = build_strain_calls(tsv, egn, gene_ann, pfam_map)
    rec = calls["WP_MREB.1"]
    # Protein-level fields
    assert rec["pfam_ids"] == ["PF06723"]
    assert rec["pfam_tc_families"] == ["9.B.157"]
    # Single 3-part family -> 1 candidate
    assert len(rec["calls"]) == 1
    c = rec["calls"][0]
    # eggNOG already carries 9.B.157.1 (the multi-valued case), so Pfam's
    # 9.B.157 corroborates BOTH diamond AND one of the eggNOG values.
    assert c["pfam_agreement"] == "confirms_both"
    # eggNOG agreement on the headline call
    assert c["egn_agreement"] == "confirms"
    pad = summary["pfam_agreement_distribution"]
    assert pad["confirms_both"] == 1
    assert pad["confirms_diamond"] == 0
    assert pad["confirms_eggnog"] == 0
    assert pad["contradicts_both"] == 0
    assert pad["neutral"] == 0


def test_build_strain_calls_pfam_paths_optional(tmp_path):
    """When pfam_map / gene_annotations paths are not supplied, calls still
    get pfam_* fields with neutral defaults so the schema is stable."""
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        "WP_X.1\tlcl|Q-1.A.11.1.5\t87.4\t92.1\t89.7\t412\t1e-180\t650\n"
    )
    egn = tmp_path / "missing.emapper.annotations"  # no eggnog file
    calls, summary = build_strain_calls(tsv, egn)  # no gene_ann / pfam_map
    rec = calls["WP_X.1"]
    assert rec["pfam_ids"] == []
    assert rec["pfam_tc_families"] == []
    assert rec["calls"][0]["pfam_agreement"] == "neutral"
    assert summary["pfam_agreement_distribution"]["neutral"] == 1


def test_build_strain_calls_multi_value_eggnog_kegg_tc(tmp_path):
    """Regression: eggNOG's KEGG_TC may carry multiple comma-separated values.
    Diamond's call matching ANY of them must produce confirms/refines, not
    conflicts. This was the MreB false-positive (~14% of phase-1 conflicts):
    eggNOG's `1.A.33.1,9.B.157.1` for MreBCD-family proteins was being
    misclassified as conflicts because only the first value was inspected.
    """
    tsv_content = (
        # Strong consensus on 9.B.157.1 (MreBCD)
        "WP_MREB.1\tlcl|Q1-9.B.157.1\t89.9\t100.0\t95.0\t340\t1e-180\t650\n"
        "WP_MREB.1\tlcl|Q2-9.B.157.1\t85.0\t99.0\t94.0\t338\t1e-170\t620\n"
    )
    tsv = tmp_path / "mreb.tcdb.tsv"
    tsv.write_text(tsv_content)

    import textwrap
    egn = tmp_path / "mreb.emapper.annotations"
    egn.write_text(textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_MREB.1\tx\tx\tx\tx\tx\tD\tx\tmreB\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.33.1,9.B.157.1\t-\t-\tMreB_Mbl
        """))

    calls, summary = build_strain_calls(tsv, egn)
    rec = calls["WP_MREB.1"]
    c = rec["calls"][0]
    assert c["tcid"] == "9.B.157.1"
    assert c["egn_agreement"] == "confirms"
    assert rec["egn_tcids"] == ["1.A.33.1", "9.B.157.1"]
    # Summary should reflect 1 confirm, 0 conflicts (the bug-era output had 1 conflict)
    assert summary["agreement_distribution"]["confirms"] == 1
    assert summary["agreement_distribution"]["conflicts"] == 0


from multiomics_kg.utils.tcdb_diamond import annotate_candidate_filters


def _mk_cand(tcid, score, egn="extends", pfam="neutral", consensus_n=3):
    """Test factory — builds a minimal candidate dict for filter tests.
    Default consensus_n=3 keeps tests focused on the rules under test (not
    on the singleton-low-score rule); set consensus_n=1 to exercise that path.
    """
    return {
        "tcid": tcid, "confidence_score": score,
        "egn_agreement": egn, "pfam_agreement": pfam,
        "consensus_n": consensus_n,
    }


def test_filter_single_candidate_kept_when_pfam_contradicts():
    # Pfam-contradicts and egn-conflicts rules require multi-candidate context;
    # a single-candidate protein can't be dropped by those even if its tags say so.
    # The candidate has consensus_n=3 so the singleton rule doesn't fire either.
    calls = {"WP_X.1": {"calls": [_mk_cand("9.X.1", 0.05, pfam="contradicts_both")]}}
    stats = annotate_candidate_filters(calls)
    assert calls["WP_X.1"]["calls"][0]["filter_action"] == "keep"
    assert stats == {"keep": 1}


def test_filter_singleton_low_score_drops_even_single_candidate():
    # consensus_n=1 AND score < 0.20 -> drop_singleton_low_score, regardless
    # of whether there are siblings. Targets the weakest-evidence case.
    calls = {"WP_S.1": {"calls": [_mk_cand("1.A.1", 0.15, consensus_n=1)]}}
    stats = annotate_candidate_filters(calls)
    assert calls["WP_S.1"]["calls"][0]["filter_action"] == "drop_singleton_low_score"
    assert stats == {"drop_singleton_low_score": 1}


def test_filter_singleton_above_threshold_kept():
    # consensus_n=1 BUT score >= threshold -> keep
    calls = {"WP_S.1": {"calls": [_mk_cand("1.A.1", 0.25, consensus_n=1)]}}
    annotate_candidate_filters(calls)
    assert calls["WP_S.1"]["calls"][0]["filter_action"] == "keep"


def test_filter_multi_hit_below_threshold_kept():
    # consensus_n>1 with low score: singleton rule doesn't fire (it requires
    # consensus_n==1). Without siblings to compare, relative-confidence rule
    # can't fire either. -> keep.
    calls = {"WP_M.1": {"calls": [_mk_cand("1.A.1", 0.15, consensus_n=3)]}}
    annotate_candidate_filters(calls)
    assert calls["WP_M.1"]["calls"][0]["filter_action"] == "keep"


def test_filter_drops_pfam_contradicting_candidate():
    # When one candidate has Pfam support and another has Pfam contradicting,
    # drop the contradicting one.
    calls = {"WP_Y.1": {"calls": [
        _mk_cand("2.A.6.1", 0.80, pfam="confirms_diamond"),
        _mk_cand("9.B.99",  0.50, pfam="contradicts_both"),
    ]}}
    annotate_candidate_filters(calls)
    actions = [c["filter_action"] for c in calls["WP_Y.1"]["calls"]]
    assert actions == ["keep", "drop_pfam_contradicts"]


def test_filter_drops_egn_conflicting_candidate():
    # When one candidate confirms eggNOG and another conflicts, drop the conflict.
    calls = {"WP_Z.1": {"calls": [
        _mk_cand("2.A.6.1", 0.80, egn="conflicts"),
        _mk_cand("8.A.1.2", 0.55, egn="confirms"),
    ]}}
    annotate_candidate_filters(calls)
    actions = [c["filter_action"] for c in calls["WP_Z.1"]["calls"]]
    assert actions == ["drop_egn_conflicts", "keep"]


def test_filter_drops_low_relative_confidence():
    # Without Pfam/eggNOG signals, the low-confidence candidate is dropped
    # if its score is < 0.25 × the best sibling's score.
    calls = {"WP_W.1": {"calls": [
        _mk_cand("2.A.6.1", 0.80),
        _mk_cand("3.A.1",   0.10),  # 0.10 < 0.25*0.80 = 0.20
    ]}}
    annotate_candidate_filters(calls)
    actions = [c["filter_action"] for c in calls["WP_W.1"]["calls"]]
    assert actions == ["keep", "drop_low_confidence"]


def test_filter_does_not_drop_low_confidence_when_alone():
    # Single candidate with score 0.05 still kept (no sibling to dominate it).
    calls = {"WP_V.1": {"calls": [_mk_cand("9.X.1", 0.05)]}}
    annotate_candidate_filters(calls)
    assert calls["WP_V.1"]["calls"][0]["filter_action"] == "keep"


def test_filter_priority_pfam_over_egn_over_confidence():
    # A candidate that fails multiple rules gets tagged by the highest-priority one.
    calls = {"WP_M.1": {"calls": [
        _mk_cand("1.A.1", 0.80, egn="confirms", pfam="confirms_diamond"),
        # Candidate below fails all 3 rules; Pfam priority wins
        _mk_cand("9.B.1", 0.05, egn="conflicts", pfam="contradicts_both"),
    ]}}
    annotate_candidate_filters(calls)
    assert calls["WP_M.1"]["calls"][1]["filter_action"] == "drop_pfam_contradicts"


def test_filter_keeps_all_when_no_confirming_sibling():
    # All candidates are 'contradicts_both' — no Pfam-supported sibling to compare,
    # so the pfam rule doesn't fire. Same for egn 'conflicts' when no eggNOG confirm.
    calls = {"WP_N.1": {"calls": [
        _mk_cand("9.A.1", 0.50, egn="conflicts", pfam="contradicts_both"),
        _mk_cand("9.B.1", 0.45, egn="conflicts", pfam="contradicts_both"),
    ]}}
    annotate_candidate_filters(calls)
    # Neither pfam_contradicts nor egn_conflicts rule fires; low_confidence
    # requires score < 0.25 × max (0.50). 0.45 > 0.125 so it stays.
    actions = [c["filter_action"] for c in calls["WP_N.1"]["calls"]]
    assert actions == ["keep", "keep"]


def test_filter_action_distribution_in_summary(tmp_path):
    """end-to-end: build_strain_calls populates filter_action and includes
    filter_action_distribution in the summary."""
    # 2-candidate protein where the second should be Pfam-contradicted
    tsv = tmp_path / "x.tcdb.tsv"
    tsv.write_text(
        # Strong hits in family 2.A.6 (Pfam-supported)
        "WP_F.1\tlcl|Q1-2.A.6.1.4\t100.0\t100.0\t100.0\t400\t1e-200\t800\n"
        # Weak hit in family 9.B.99 (Pfam contradicts)
        "WP_F.1\tlcl|Q2-9.B.99.1.1\t30.0\t40.0\t40.0\t150\t1e-10\t100\n"
    )
    egn = tmp_path / "missing.emapper.annotations"
    gene_ann = tmp_path / "ga.json"
    gene_ann.write_text(json.dumps({
        "LT01": {"protein_id": "WP_F.1", "pfam_ids": ["PFA"]},
    }))
    pfam_map = tmp_path / "pmap.tsv"
    pfam_map.write_text("PFA\t2.A.6.1.1\tRND-superfamily\n")
    calls, summary = build_strain_calls(tsv, egn, gene_ann, pfam_map)
    actions = [c["filter_action"] for c in calls["WP_F.1"]["calls"]]
    # First (2.A.6) is confirms_diamond, kept. Second (9.B.99) is
    # contradicts_both with a confirming sibling -> dropped.
    assert "keep" in actions
    assert "drop_pfam_contradicts" in actions
    assert summary["filter_action_distribution"]["drop_pfam_contradicts"] == 1


def test_build_strain_calls_emits_multi_family_candidates(tmp_path):
    """The headline multi-call recovery: a protein with diamond hits in two
    distinct 3-part families (e.g. RND + MFP partner-protein confusion) gets
    BOTH as candidates instead of being rejected by global consensus. Each
    candidate carries its own egn/pfam agreement.
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
    import textwrap
    egn = tmp_path / "x.emapper.annotations"
    egn.write_text(textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_RND.1\tx\tx\tx\tx\tx\tM\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t8.A.1.2.1\t-\t-\tHlyD_D23
        """))

    calls, summary = build_strain_calls(tsv, egn)
    rec = calls["WP_RND.1"]
    assert len(rec["calls"]) == 2
    # First candidate is the higher-confidence one (RND, 100/100)
    assert rec["calls"][0]["tcid"].startswith("2.A.6")
    assert rec["calls"][0]["confidence_score"] > rec["calls"][1]["confidence_score"]
    # Per-candidate egn agreement: RND vs eggNOG's MFP -> conflicts;
    # MFP vs eggNOG's MFP -> confirms
    by_tcid = {c["tcid"][:5]: c for c in rec["calls"]}
    assert by_tcid["2.A.6"]["egn_agreement"] == "conflicts"
    assert by_tcid["8.A.1"]["egn_agreement"] == "confirms"
    # Summary reflects per-candidate counts
    assert summary["total_candidates"] == 2
    assert summary["proteins_with_call"] == 1
    assert summary["candidates_per_protein_distribution"] == {"2": 1}


# ============================================================================
# Phase 2 — calls.json consumption helpers
# ============================================================================

from multiomics_kg.utils.tcdb_diamond import (
    KEEP_ACTION,
    FILTER_ACTIONS,
    load_calls_json,
    iter_kept_candidates,
    best_kept_call,
    kept_tc_families,
    kept_call_tcids,
    summarize_filter_actions,
)


def _mk_rec(*candidates):
    """Test factory — build a per-protein record from positional candidates.
    Each candidate is a dict (use _mk_kept_cand / _mk_drop_cand helpers).
    """
    return {
        "egn_tcids": [], "pfam_ids": [], "pfam_tc_families": [],
        "calls": list(candidates),
    }


def _kept(tcid, score):
    return {"tcid": tcid, "confidence_score": score, "filter_action": "keep"}


def _dropped(tcid, score, reason="drop_low_confidence"):
    return {"tcid": tcid, "confidence_score": score, "filter_action": reason}


def test_load_calls_json_reads_on_disk_shape(tmp_path):
    p = tmp_path / "x.tcdb.calls.json"
    payload = {"WP_X.1": _mk_rec(_kept("1.A.11.1.5", 0.8))}
    p.write_text(json.dumps(payload))
    loaded = load_calls_json(p)
    assert loaded == payload


def test_load_calls_json_missing_file_raises():
    import pytest
    with pytest.raises(FileNotFoundError):
        load_calls_json("/nonexistent/path/calls.json")


def test_iter_kept_candidates_yields_only_keeps():
    calls = {
        "WP_A.1": _mk_rec(_kept("1.A.1", 0.8), _dropped("2.A.7", 0.1)),
        "WP_B.1": _mk_rec(_dropped("9.X.1", 0.05, "drop_singleton_low_score")),
        "WP_C.1": _mk_rec(_kept("3.A.5", 0.6)),
    }
    seen = list(iter_kept_candidates(calls))
    assert [(pid, c["tcid"]) for pid, c in seen] == [("WP_A.1", "1.A.1"), ("WP_C.1", "3.A.5")]


def test_iter_kept_candidates_preserves_per_protein_order():
    # calls[] is already sorted by confidence_score desc; iter should follow
    rec = _mk_rec(
        _kept("2.A.6.1", 0.9),
        _kept("8.A.1", 0.5),
        _kept("9.B.99", 0.3),
    )
    seen = [c["tcid"] for _pid, c in iter_kept_candidates({"WP_X.1": rec})]
    assert seen == ["2.A.6.1", "8.A.1", "9.B.99"]


def test_best_kept_call_returns_top_keep():
    rec = _mk_rec(
        _dropped("1.A.1", 0.95),  # higher score but dropped — should NOT win
        _kept("2.A.6.1", 0.8),
        _kept("8.A.1", 0.5),
    )
    assert best_kept_call(rec)["tcid"] == "2.A.6.1"


def test_best_kept_call_none_when_all_dropped():
    rec = _mk_rec(_dropped("9.X.1", 0.05, "drop_singleton_low_score"))
    assert best_kept_call(rec) is None


def test_best_kept_call_empty_calls_list():
    assert best_kept_call({"calls": []}) is None
    assert best_kept_call({}) is None  # missing calls key


def test_kept_tc_families_collapses_to_3_part_and_dedupes():
    rec = _mk_rec(
        _kept("1.A.11.1.5", 0.8),     # -> 1.A.11
        _kept("1.A.11.2.1", 0.5),     # -> 1.A.11 (dedup)
        _kept("2.A.6.1", 0.4),        # -> 2.A.6
        _dropped("9.B.99.1", 0.05, "drop_singleton_low_score"),  # excluded
    )
    assert kept_tc_families(rec) == ["1.A.11", "2.A.6"]


def test_kept_tc_families_empty_when_none_kept():
    rec = _mk_rec(_dropped("1.A.1", 0.05))
    assert kept_tc_families(rec) == []


def test_kept_call_tcids_preserves_depth_and_order():
    rec = _mk_rec(
        _kept("1.A.11.1.5", 0.8),    # 5-part
        _kept("2.A.6.1", 0.5),       # 4-part
        _kept("8.A.1", 0.3),         # 3-part
        _dropped("9.B", 0.1),
    )
    assert kept_call_tcids(rec) == ["1.A.11.1.5", "2.A.6.1", "8.A.1"]


def test_summarize_filter_actions_recounts():
    calls = {
        "WP_A.1": _mk_rec(_kept("1.A", 0.5), _dropped("2.A", 0.1, "drop_low_confidence")),
        "WP_B.1": _mk_rec(_dropped("9.X", 0.05, "drop_singleton_low_score")),
    }
    counts = summarize_filter_actions(calls)
    assert counts["keep"] == 1
    assert counts["drop_low_confidence"] == 1
    assert counts["drop_singleton_low_score"] == 1
    assert counts["drop_pfam_contradicts"] == 0  # zero-init
    # All keys present
    assert set(counts) == set(FILTER_ACTIONS)


def test_phase2_vocabulary_constants_exposed():
    # Phase 2 callers can import these instead of hardcoding strings
    assert KEEP_ACTION == "keep"
    assert "keep" in FILTER_ACTIONS
    assert "drop_singleton_low_score" in FILTER_ACTIONS
    assert len(FILTER_ACTIONS) == 5
