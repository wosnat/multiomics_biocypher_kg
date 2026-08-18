"""Tests for the faceted InterProScan parser (multi-ontology redesign)."""
import pytest
from multiomics_kg.utils.interproscan import parse_interproscan_json, summarize


def _loc(start, end, evalue=None, score=None):
    return {"start": start, "end": end, "evalue": evalue, "score": score}


def _match(library, acc, desc, entry=None, locations=None):
    sig = {"accession": acc, "description": desc,
           "signatureLibraryRelease": {"library": library}}
    if entry:
        sig["entry"] = entry
    return {"signature": sig, "locations": locations or [_loc(1, 50)]}


ENTRY_FAM = {
    "accession": "IPR003686", "description": "Photosystem II PsbI", "type": "FAMILY",
    "goXRefs": [{"id": "GO:0015979"}, {"id": "GO:0009523"}],
    "pathwayXRefs": [{"databaseName": "MetaCyc", "id": "PWY-101"}],  # must be DROPPED
}
ENTRY_SF = {
    "accession": "IPR037271", "description": "PsbI superfamily",
    "type": "HOMOLOGOUS_SUPERFAMILY", "goXRefs": [{"id": "GO:0015979"}],
    "pathwayXRefs": [],
}

RAW = {"results": [{
    "md5": "abc", "xref": [{"id": "WP_000001.1"}],
    "matches": [
        _match("PFAM", "PF02532.18", "PSII PsbI", ENTRY_FAM,
               [_loc(1, 36, evalue=4.1e-18, score=76.3)]),
        _match("HAMAP", "MF_01316", "PSII reaction center I", ENTRY_FAM,
               [_loc(1, 36, score=17.4)]),
        _match("NCBIFAM", "NF002735.2", "photosystem II protein I", None,
               [_loc(1, 38, evalue=3.3e-23, score=92.7)]),
        _match("SUPERFAMILY", "SSF161041", "PsbI", ENTRY_SF, [_loc(1, 35)]),
    ],
}, {
    "md5": "def", "xref": [{"id": "WP_000002.1"}], "matches": [],
}]}


@pytest.fixture()
def calls():
    return parse_interproscan_json(RAW)


def test_libraries_facet_sparse_and_version_stripped(calls):
    rec = calls["WP_000001.1"]
    assert set(rec["libraries"]) == {"PFAM", "HAMAP", "NCBIFAM", "SUPERFAMILY"}
    pf = rec["libraries"]["PFAM"][0]
    assert pf["accession"] == "PF02532"          # version stripped
    assert pf["ipr"] == "IPR003686"
    assert pf["evalue"] == 4.1e-18 and pf["score"] == 76.3
    nf = rec["libraries"]["NCBIFAM"][0]
    assert nf["accession"] == "NF002735" and nf["ipr"] is None


def test_interpro_rollup_no_score_evalue_attributed(calls):
    ent = calls["WP_000001.1"]["interpro_entries"]["IPR003686"]
    assert ent["type"] == "FAMILY"
    assert ent["libraries"] == ["HAMAP", "PFAM"]
    assert ent["match_count"] == 2
    assert ent["evalue"] == 4.1e-18 and ent["evalue_library"] == "PFAM"
    assert "score" not in ent                     # count-don't-combine
    sf = calls["WP_000001.1"]["interpro_entries"]["IPR037271"]
    assert sf["evalue"] is None and sf["evalue_library"] is None


def test_go_terms_carry_entry_attribution(calls):
    go = calls["WP_000001.1"]["go_terms"]
    assert go["GO:0015979"] == ["IPR003686", "IPR037271"]
    assert go["GO:0009523"] == ["IPR003686"]


def test_no_pathways_anywhere(calls):
    rec = calls["WP_000001.1"]
    assert "pathways" not in rec
    assert all("pathways" not in e for e in rec["interpro_entries"].values())


def test_zero_match_protein_kept(calls):
    rec = calls["WP_000002.1"]
    assert rec["match_count"] == 0 and rec["libraries"] == {} \
        and rec["interpro_entries"] == {} and rec["go_terms"] == {}


def test_summarize_qc():
    calls = parse_interproscan_json(RAW)
    s = summarize(calls, strain="X", input_proteins=2,
                  tool_version="5.78-109.0", applications="ALL_DEFAULT")
    assert s["calls_made"] == 2 and s["proteins_no_match"] == 1
    assert s["total_matches"] == 4
    assert s["interpro_integrated_matches"] == 3
    assert s["distribution"] == {"HAMAP": 1, "NCBIFAM": 1, "PFAM": 1, "SUPERFAMILY": 1}
    assert s["proteins_with_go_terms"] == 1 and s["distinct_go_terms"] == 2
    assert "pathway_databases" not in s and "distinct_pathways" not in s
