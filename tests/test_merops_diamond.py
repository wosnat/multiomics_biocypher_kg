# tests/test_merops_diamond.py
"""Unit tests for multiomics_kg.utils.merops_diamond."""
import json

from multiomics_kg.utils.merops_diamond import (
    best_call,
    build_strain_calls,
    call_codes,
    call_families,
    classify_hit,
    confidence_score,
    family_group_agreement,
    id_family,
    id_kind,
    iter_candidates,
    load_calls_json,
    parse_diamond_row,
    parse_family_txt,
    parse_merops_subject_id,
    parse_scan_lib_header,
    _truncate_call,
)


# ----------------------------------------------------------------------------
# Header / identifier parsing
# ----------------------------------------------------------------------------

HEADER_WITH_SUBFAM = (
    ">MER0000002 - chymotrypsin A (cattle-type) (Bos taurus) [S01.001]#S01A#"
    "{peptidase unit: 34-263}~source XP_608091~\n"
)
HEADER_NO_SUBFAM = (
    ">MER0026205 - prolyl aminopeptidase (Thermoplasma-type) (Thermoplasma "
    "acidophilum) [S33.008]#S33#{peptidase unit: 1-280}~source Q9HIA0~\n"
)


def test_parse_scan_lib_header_with_subfamily():
    rec = parse_scan_lib_header(HEADER_WITH_SUBFAM)
    assert rec == {
        "mernum": "MER0000002",
        "name": "chymotrypsin A (cattle-type) (Bos taurus)",
        "merops_id": "S01.001",
        "subfamily": "S01A",
    }


def test_parse_scan_lib_header_family_only_token():
    rec = parse_scan_lib_header(HEADER_NO_SUBFAM)
    assert rec["merops_id"] == "S33.008"
    assert rec["subfamily"] == "S33"


def test_parse_scan_lib_header_rejects_garbage():
    assert parse_scan_lib_header(">not a merops header\n") is None
    assert parse_scan_lib_header("ACDEFGHIKLMNPQRSTVWY\n") is None


def test_id_family():
    assert id_family("S01.001") == "S01"
    assert id_family("M16.A03") == "M16"
    assert id_family("I04.953") == "I04"
    assert id_family("garbage") is None
    assert id_family("") is None


def test_id_kind_holotype_putative_homolog():
    assert id_kind("S01.001") == "holotype"
    assert id_kind("S01.152") == "holotype"
    assert id_kind("S01.A33") == "putative"
    assert id_kind("S01.B01") == "putative"
    assert id_kind("S01.971") == "nonpeptidase_homolog"
    assert id_kind("C19.900") == "nonpeptidase_homolog"
    assert id_kind("nonsense") is None


def test_parse_merops_subject_id_with_real_subfamily():
    rec = parse_merops_subject_id("MER0000002|S01.001|S01A")
    assert rec == {
        "mernum": "MER0000002",
        "merops_id": "S01.001",
        "family": "S01",
        "subfamily": "S01A",
    }


def test_parse_merops_subject_id_family_token_maps_to_none():
    rec = parse_merops_subject_id("MER0026205|S33.008|S33")
    assert rec["subfamily"] is None
    assert rec["family"] == "S33"


def test_parse_merops_subject_id_rejects_malformed():
    assert parse_merops_subject_id("") is None
    assert parse_merops_subject_id("MER123-S01.001") is None
    assert parse_merops_subject_id("MER123|garbage|S01A") is None
    assert parse_merops_subject_id("MER123|S01.001") is None


# ----------------------------------------------------------------------------
# family.txt parsing
# ----------------------------------------------------------------------------

FAMILY_TXT = (
    "A01\tAA\t\tpepsin\t\tBiochem.J. 290:205-218 (1993)\tpeptidase\tYes\t\n"
    "A05\tA-\t\tthermopsin\t\tProteolysis in Cell Function (1997)\tpeptidase\tYes\t\n"
    "I04\tID\t\tserpin\t\tsome ref\tinhibitor\tYes\t\n"
    "short\tline\n"
)


def test_parse_family_txt_basic_row():
    fams = parse_family_txt(FAMILY_TXT)
    assert fams["A01"] == {"clan": "AA", "name": "pepsin", "entry_type": "peptidase"}
    assert fams["I04"]["entry_type"] == "inhibitor"


def test_parse_family_txt_unassigned_clan_sentinel_is_none():
    fams = parse_family_txt(FAMILY_TXT)
    assert fams["A05"]["clan"] is None


def test_parse_family_txt_skips_short_rows():
    fams = parse_family_txt(FAMILY_TXT)
    assert "short" not in fams


# ----------------------------------------------------------------------------
# Tier policy (thresholds are deliberate tcdb-diamond parity)
# ----------------------------------------------------------------------------

def _hit(identity, qcov, scov, length=200, evalue=1e-10):
    return {
        "identity": identity, "qcov": qcov, "scov": scov,
        "length": length, "evalue": evalue,
    }


def test_classify_tier_1():
    assert classify_hit(_hit(80.0, 85.0, 50.0)) == 1
    assert classify_hit(_hit(70.0, 70.0, 50.0)) == 1


def test_classify_tier_2():
    assert classify_hit(_hit(50.0, 65.0, 30.0)) == 2
    assert classify_hit(_hit(40.0, 60.0, 30.0)) == 2


def test_classify_tier_3_floor():
    assert classify_hit(_hit(35.0, 45.0, 30.0)) == 3
    assert classify_hit(_hit(35.0, 25.0, 50.0)) == 3  # scov OR rule
    assert classify_hit(_hit(25.0, 45.0, 30.0)) == 3  # no identity floor


def test_classify_drops_below_floor():
    assert classify_hit(_hit(35.0, 30.0, 30.0)) is None      # both cov < 40
    assert classify_hit(_hit(80.0, 90.0, 90.0, length=49)) is None
    assert classify_hit(_hit(80.0, 90.0, 90.0, evalue=0.01)) is None


# ----------------------------------------------------------------------------
# Consensus + truncation
# ----------------------------------------------------------------------------

def _subject(merops_id, subfamily, mernum="MER1", tier=1):
    return {
        "mernum": mernum, "merops_id": merops_id,
        "family": id_family(merops_id), "subfamily": subfamily,
        "tier": tier, "identity": 80.0, "qcov": 80.0, "scov": 80.0,
        "evalue": 1e-50, "length": 200,
    }


def test_agreement_id_level():
    hits = [_subject("S08.001", "S08A"), _subject("S08.001", "S08A")]
    assert family_group_agreement(hits) == "id"


def test_agreement_subfamily_level():
    hits = [_subject("S08.001", "S08A"), _subject("S08.036", "S08A")]
    assert family_group_agreement(hits) == "subfamily"


def test_agreement_family_level_mixed_subfamilies():
    hits = [_subject("S08.001", "S08A"), _subject("S08.100", "S08B")]
    assert family_group_agreement(hits) == "family"


def test_agreement_family_level_when_no_subfamilies_and_ids_differ():
    # Family without subfamilies (token == family → parsed subfamily None):
    # id disagreement can only collapse to family, never "subfamily".
    hits = [_subject("S33.001", None), _subject("S33.008", None)]
    assert family_group_agreement(hits) == "family"


def test_truncate_tier1_claims_identifier():
    assert _truncate_call("S08", "S08A", "S08.036", 1) == ("S08.036", "merops_id", "S08A")


def test_truncate_tier2_claims_subfamily():
    assert _truncate_call("S08", "S08A", "S08.036", 2) == ("S08A", "merops_subfamily", "S08A")


def test_truncate_tier2_ragged_family_without_subfamilies():
    # The deliberate TCDB deviation: tier 2 confidence, family-level claim.
    assert _truncate_call("S33", None, "S33.008", 2) == ("S33", "merops_family", None)


def test_truncate_tier3_claims_family():
    assert _truncate_call("S08", "S08A", "S08.036", 3) == ("S08", "merops_family", None)


def test_confidence_score_formula():
    assert confidence_score(100.0, 100.0, "id") == 1.0
    assert confidence_score(100.0, 100.0, "subfamily") == 0.85
    assert abs(confidence_score(50.0, 80.0, "family") - 0.28) < 1e-9


# ----------------------------------------------------------------------------
# build_strain_calls end-to-end
# ----------------------------------------------------------------------------

def _tsv_row(query, subject, identity, qcov, scov, length=200, evalue=1e-50, bits=500.0):
    return f"{query}\t{subject}\t{identity}\t{qcov}\t{scov}\t{length}\t{evalue}\t{bits}\n"


FAMILY_META = {
    "S08": {"clan": "SB", "name": "subtilisin", "entry_type": "peptidase"},
    "S33": {"clan": "SC", "name": "prolyl aminopeptidase", "entry_type": "peptidase"},
    "I04": {"clan": "ID", "name": "serpin", "entry_type": "inhibitor"},
}


def test_build_strain_calls_id_consensus(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        _tsv_row("WP_1", "MER1|S08.036|S08A", 85.0, 90.0, 90.0)
        + _tsv_row("WP_1", "MER2|S08.036|S08A", 80.0, 85.0, 85.0)
    )
    calls, summary = build_strain_calls(tsv, FAMILY_META)
    cand = calls["WP_1"]["calls"][0]
    assert cand["code"] == "S08.036"
    assert cand["level_kind"] == "merops_id"
    assert cand["tier"] == 1
    assert cand["clan"] == "SB"
    assert cand["catalytic_type"] == "S"
    assert cand["entry_type"] == "peptidase"
    assert cand["consensus_n"] == 2
    assert cand["consensus_agreement"] == "id"
    assert cand["best_hit_mernum"] == "MER1"  # highest identity
    assert cand["best_hit_kind"] == "holotype"
    assert cand["homolog_hit_fraction"] == 0.0
    assert summary["proteins_with_call"] == 1
    assert summary["tier_distribution"] == {"1": 1}


def test_build_strain_calls_subfamily_disagreement_truncates(tmp_path):
    # Two different ids in the same subfamily, tier-1 evidence → depth caps at
    # subfamily → effective tier 2, subfamily-level claim.
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        _tsv_row("WP_1", "MER1|S08.036|S08A", 85.0, 90.0, 90.0)
        + _tsv_row("WP_1", "MER2|S08.001|S08A", 80.0, 85.0, 85.0)
    )
    calls, _ = build_strain_calls(tsv, FAMILY_META)
    cand = calls["WP_1"]["calls"][0]
    assert cand["code"] == "S08A"
    assert cand["level_kind"] == "merops_subfamily"
    assert cand["tier"] == 2
    assert cand["consensus_agreement"] == "subfamily"


def test_build_strain_calls_multi_family_emits_multiple_candidates(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        _tsv_row("WP_1", "MER1|S08.036|S08A", 85.0, 90.0, 90.0)
        + _tsv_row("WP_1", "MER2|I04.001|I04A", 45.0, 65.0, 65.0)
    )
    calls, summary = build_strain_calls(tsv, FAMILY_META)
    cands = calls["WP_1"]["calls"]
    assert len(cands) == 2
    assert cands[0]["family"] == "S08"  # higher confidence first
    inhibitor = cands[1]
    assert inhibitor["entry_type"] == "inhibitor"
    assert inhibitor["catalytic_type"] is None
    assert summary["entry_type_distribution"] == {"peptidase": 1, "inhibitor": 1}
    assert summary["catalytic_type_distribution"] == {"S": 1, "inhibitor": 1}


def test_build_strain_calls_homolog_flagging(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        _tsv_row("WP_1", "MER1|S08.971|S08A", 85.0, 90.0, 90.0)
        + _tsv_row("WP_1", "MER2|S08.972|S08A", 60.0, 80.0, 80.0)
    )
    calls, summary = build_strain_calls(tsv, FAMILY_META)
    cand = calls["WP_1"]["calls"][0]
    assert cand["best_hit_kind"] == "nonpeptidase_homolog"
    assert cand["homolog_hit_fraction"] == 1.0
    assert summary["best_hit_kind_distribution"] == {"nonpeptidase_homolog": 1}


def test_build_strain_calls_unknown_family_gets_null_clan(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(_tsv_row("WP_1", "MER1|M99.001|M99", 85.0, 90.0, 90.0))
    calls, _ = build_strain_calls(tsv, {})
    cand = calls["WP_1"]["calls"][0]
    assert cand["clan"] is None
    assert cand["family"] == "M99"


def test_build_strain_calls_counts_parse_failures(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(
        "not\tenough\tcolumns\n"
        + _tsv_row("WP_1", "unparseable-subject", 85.0, 90.0, 90.0)
        + _tsv_row("WP_2", "MER1|S08.036|S08A", 85.0, 90.0, 90.0)
    )
    calls, summary = build_strain_calls(tsv, FAMILY_META)
    assert summary["parse_failures"] == 2
    assert summary["raw_hit_lines"] == 3
    assert list(calls) == ["WP_2"]


def test_build_strain_calls_drops_below_floor_rows(tmp_path):
    tsv = tmp_path / "x.tsv"
    tsv.write_text(_tsv_row("WP_1", "MER1|S08.036|S08A", 85.0, 30.0, 30.0))
    calls, summary = build_strain_calls(tsv, FAMILY_META)
    assert calls == {}
    assert summary["proteins_with_hits"] == 0
    assert summary["parse_failures"] == 0


# ----------------------------------------------------------------------------
# Consumption helpers
# ----------------------------------------------------------------------------

def test_consumption_helpers(tmp_path):
    rec = {
        "calls": [
            {"code": "S08.036", "family": "S08", "confidence_score": 0.9},
            {"code": "I04", "family": "I04", "confidence_score": 0.3},
        ]
    }
    path = tmp_path / "X.merops.calls.json"
    path.write_text(json.dumps({"WP_1": rec}))

    calls = load_calls_json(path)
    assert [pid for pid, _ in iter_candidates(calls)] == ["WP_1", "WP_1"]
    assert best_call(rec)["code"] == "S08.036"
    assert call_families(rec) == ["I04", "S08"]
    assert call_codes(rec) == ["S08.036", "I04"]
    assert best_call({"calls": []}) is None


def test_parse_diamond_row_malformed():
    assert parse_diamond_row("a\tb\tc\n") is None
    assert parse_diamond_row("q\ts\tNaNo\t1\t2\t3\t4\t5\n") is None
    row = parse_diamond_row("q\ts\t85.5\t90\t80\t200\t1e-50\t500\n")
    assert row["identity"] == 85.5 and row["length"] == 200


# ============================================================================
# Phase-2: Pfam bridge + cleavage specificity parsers
# ============================================================================

from multiomics_kg.utils.merops_diamond import (
    STANDARD_AA_3,
    aggregate_cleavages,
    cleavage_properties,
    parse_interpro_txt_stream,
)


def test_parse_interpro_txt_stream_counts_distinct_identifiers():
    lines = [
        '"A01A","A01.001","pepsin A","63-388","P0DJD8","PF00026/PF14543","IPR001461"',
        '"A01A","A01.001","pepsin A","63-388","B7Z719","PF00026/PF14543","IPR001461"',  # same id, dup accession
        '"A01A","A01.002","pepsin B","60-380","P27821","PF00026","IPR001461"',
        '"S01A","S01.001","chymotrypsin A","34-263","P00766","PF00089",""',
        '"M10B","M10.051","serralysin","1-470","P23694","",""',  # no Pfam -> ignored
    ]
    bridge = parse_interpro_txt_stream(lines)
    assert bridge["A01"]["PF00026"] == 2      # A01.001 + A01.002, dedup on identifier
    assert bridge["A01"]["PF14543"] == 1      # only A01.001
    assert bridge["S01"]["PF00089"] == 1
    assert "M10" not in bridge


def test_aggregate_cleavages_filters_nonstandard_p1():
    q = chr(39)  # single quote — the file quotes fields with it
    def row(ident, p1, kind):
        f = [q + "CLE1" + q, q + ident + q, "s", "s", "-", "-", "-",
             q + p1 + q, "-", "-", "-", "-", "ref", "NULL", "NULL", "NULL",
             "e", "NULL", "NULL", "NULL", "NULL", "NULL", q + kind + q]
        return "\t".join(f)
    agg = aggregate_cleavages([
        row("S01.001", "Lys", "physiological"),
        row("S01.001", "Arg", "synthetic"),
        row("S01.151", "TyI", "non-physiological"),   # modified residue: counted in total, not in p1
        row("A01.001", "Phe", "physiological"),
    ])
    assert agg["S01"]["total"] == 3
    assert agg["S01"]["physiological"] == 1
    assert agg["S01"]["p1"]["Lys"] == 1 and agg["S01"]["p1"]["Arg"] == 1
    assert "TyI" not in agg["S01"]["p1"]
    assert agg["A01"]["total"] == 1


def test_cleavage_properties_shapes():
    from collections import Counter
    props = cleavage_properties({
        "p1": Counter({"Lys": 36, "Arg": 34, "Glu": 11, "Ala": 9}),
        "physiological": 25, "total": 100,
    })
    assert props["cleavage_p1_residues"] == ["Lys", "Arg", "Glu"]  # >=10% share, max 3 (Ala at exactly 10% loses to top-3)
    assert props["known_cleavage_count"] == 100
    s = props["cleavage_summary"]
    # shares over the standard-P1 subtotal (90): 36/90=40%, 34/90=38%, 11/90=12%
    assert s == "cleaves after Lys (40%) / Arg (38%) / Glu (12%) - 100 known cleavages (25% physiological)"
    assert "'" not in s and "|" not in s
    # no standard-P1 data -> count-only summary, no residues key
    props2 = cleavage_properties({"p1": Counter(), "physiological": 3, "total": 10})
    assert "cleavage_p1_residues" not in props2
    assert props2["cleavage_summary"] == "10 known cleavages (30% physiological)"
    assert props2["known_cleavage_count"] == 10
    # no data at all -> empty (sparse discipline)
    assert cleavage_properties({"p1": Counter(), "physiological": 0, "total": 0}) == {}
