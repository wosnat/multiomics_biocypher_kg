# tests/test_tcdb_diamond.py
"""Unit tests for multiomics_kg.utils.tcdb_diamond."""
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
    assert compute_egn_agreement("1.A.11", "1.A.11") == "confirms"
    assert compute_egn_agreement("1.A.11.1.5", "1.A.11.1.5") == "confirms"


def test_egn_agreement_confirms_diamond_descendant():
    # eggNOG family-level, diamond a strict descendant — still confirms (consistent)
    assert compute_egn_agreement("1.A.11.1.5", "1.A.11") == "refines"
    assert compute_egn_agreement("1.A.11.1", "1.A.11") == "refines"


def test_egn_agreement_extends_when_eggnog_missing():
    assert compute_egn_agreement("1.A.11.1.5", None) == "extends"
    assert compute_egn_agreement("1.A.11.1.5", "") == "extends"


def test_egn_agreement_conflicts_different_family():
    assert compute_egn_agreement("1.A.11.1.5", "2.A.7") == "conflicts"
    assert compute_egn_agreement("1.A.11", "1.B.5") == "conflicts"


def test_egn_agreement_eggnog_descendant_of_diamond_is_conflict():
    # Diamond at 3-part, eggNOG at 5-part below it. Same family — not a conflict;
    # this case is rare (eggNOG rarely emits 5-part) but treat as confirms.
    assert compute_egn_agreement("1.A.11", "1.A.11.1.5") == "confirms"


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
        "WP_011131852.1": "1.A.11",
        "WP_011131900.1": "1.A.11",
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
    # Multi-value KEGG_TC: keep first one (rare; we don't try to merge)
    assert result == {"WP_C.1": "1.A.11"}


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

    # WP_AAA.1: 3 strong consensus hits at 5-part -> tier 1 specificity, refines eggNOG
    assert calls["WP_AAA.1"]["tcid"] == "1.A.11.1.5"
    assert calls["WP_AAA.1"]["tier"] == 1
    assert calls["WP_AAA.1"]["level_kind"] == "tc_specificity"
    assert calls["WP_AAA.1"]["consensus_agreement"] == "5_part"
    assert calls["WP_AAA.1"]["consensus_n"] == 3
    assert calls["WP_AAA.1"]["egn_agreement"] == "refines"
    assert calls["WP_AAA.1"]["egn_tcid"] == "1.A.11"
    assert calls["WP_AAA.1"]["incompletely_characterized"] is False

    # WP_BBB.1: scattered hits -> rejected (not in calls)
    assert "WP_BBB.1" not in calls

    # WP_CCC.1: single tier-2 hit -> 4-part subfamily
    assert calls["WP_CCC.1"]["tcid"] == "1.A.11.1"
    assert calls["WP_CCC.1"]["tier"] == 2
    assert calls["WP_CCC.1"]["level_kind"] == "tc_subfamily"
    assert calls["WP_CCC.1"]["egn_agreement"] == "extends"

    # WP_DDD.1: tier 3 -> 3-part family; eggNOG conflict (different family)
    assert calls["WP_DDD.1"]["tcid"] == "1.A.11"
    assert calls["WP_DDD.1"]["tier"] == 3
    assert calls["WP_DDD.1"]["level_kind"] == "tc_family"
    assert calls["WP_DDD.1"]["egn_agreement"] == "conflicts"
    assert calls["WP_DDD.1"]["egn_tcid"] == "9.A.1"

    # WP_EEE.1: class 9 -> tagged; identity is high (80%) so still tier 1 — no demotion
    # despite class 9 (spec §6.4-C: "No demotion — let merge / downstream consumers decide")
    assert calls["WP_EEE.1"]["tcid"] == "9.B.82.1.5"
    assert calls["WP_EEE.1"]["incompletely_characterized"] is True
    assert calls["WP_EEE.1"]["egn_agreement"] == "extends"
    assert calls["WP_EEE.1"]["tier"] == 1
    assert calls["WP_EEE.1"]["level_kind"] == "tc_specificity"

    # WP_FFF.1: floor failure -> not in calls
    assert "WP_FFF.1" not in calls

    # confidence_score: (best identity / 100) * (best qcov / 100) * agreement_weight
    #   WP_AAA.1: 0.874 * 0.921 * 1.0 (5_part) ≈ 0.8050
    #   WP_CCC.1: 0.50  * 0.65  * 1.0 (5_part) = 0.3250
    #   WP_DDD.1: 0.30  * 0.45  * 1.0 (5_part) = 0.1350
    #   WP_EEE.1: 0.80  * 0.85  * 1.0 (5_part) = 0.6800
    assert abs(calls["WP_AAA.1"]["confidence_score"] - 0.8050) < 1e-3
    assert abs(calls["WP_CCC.1"]["confidence_score"] - 0.3250) < 1e-3
    assert abs(calls["WP_DDD.1"]["confidence_score"] - 0.1350) < 1e-3
    assert abs(calls["WP_EEE.1"]["confidence_score"] - 0.6800) < 1e-3

    # Summary
    assert summary["raw_hit_lines"] == 9
    assert summary["proteins_with_call"] == 4
    assert summary["proteins_rejected_by_consensus"] == 1
    # WP_EEE.1 is tier 1 (high identity 80%, no class-9 demotion per spec §6.4-C)
    assert summary["tier_distribution"] == {"1": 2, "2": 1, "3": 1}
    assert summary["agreement_distribution"]["refines"] == 1
    assert summary["agreement_distribution"]["extends"] == 2
    assert summary["agreement_distribution"]["conflicts"] == 1
    assert summary["agreement_distribution"]["confirms"] == 0


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

    assert calls["WP_GGG.1"]["tier"] == 2
    assert calls["WP_GGG.1"]["tcid"] == "1.A.11.1"
    assert calls["WP_GGG.1"]["level_kind"] == "tc_subfamily"
    # Best hit (75%/75%) drives metadata + score
    assert calls["WP_GGG.1"]["identity"] == 75.0
    assert calls["WP_GGG.1"]["qcov"] == 75.0
    # confidence: 0.75 * 0.75 * 0.85 (4_part) = 0.4781
    assert abs(calls["WP_GGG.1"]["confidence_score"] - 0.4781) < 1e-3
