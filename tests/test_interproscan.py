"""Unit tests for multiomics_kg.utils.interproscan (pure JSON parsing/summary)."""

from multiomics_kg.utils.interproscan import (
    normalize_pathway_xref,
    parse_entry_xrefs,
    parse_interproscan_json,
    summarize,
)


# A trimmed but structurally faithful InterProScan --formats JSON document:
# two proteins with matches (one integrated + one un-integrated + one
# pattern-only signature with no evalue), and one processed protein with no
# matches (the sentinel case).
SAMPLE = {
    "results": [
        {
            "md5": "aaa",
            "xref": [{"id": "WP_002805854.1", "name": "WP_002805854.1 rbcL"}],
            "matches": [
                {
                    "signature": {
                        "accession": "PF00016",
                        "name": "RuBisCO_large",
                        "description": "Ribulose bisphosphate carboxylase large chain",
                        "signatureLibraryRelease": {"library": "PFAM", "version": "37.0"},
                        "entry": {
                            "accession": "IPR000685",
                            "description": "RuBisCO large subunit, C-terminal",
                            "type": "DOMAIN",
                            "goXRefs": [{"id": "GO:0016984", "databaseName": "GO"}],
                            "pathwayXRefs": [{"id": "00710", "databaseName": "KEGG"}],
                        },
                    },
                    "locations": [{"start": 20, "end": 460, "evalue": 1e-80, "score": 300.0}],
                },
                {
                    # Un-integrated member-DB hit (entry is null).
                    "signature": {
                        "accession": "G3DSA:3.20.20.110",
                        "name": None,
                        "description": None,
                        "signatureLibraryRelease": {"library": "GENE3D", "version": "4.3"},
                        "entry": None,
                    },
                    "locations": [{"start": 15, "end": 470, "evalue": 2e-60, "score": 210.0}],
                },
            ],
        },
        {
            "md5": "bbb",
            "xref": [{"id": "WP_002805169.1", "name": "WP_002805169.1 atpH"}],
            "matches": [
                {
                    # Pattern hit — no evalue/score.
                    "signature": {
                        "accession": "PS00000",
                        "name": "ATP_C",
                        "description": "ATP synthase c subunit signature",
                        "signatureLibraryRelease": {"library": "PROSITE_PATTERNS", "version": "2024"},
                        "entry": {
                            "accession": "IPR000454",
                            "description": "ATP synthase F0, subunit c",
                            "type": "FAMILY",
                            "goXRefs": [],
                            "pathwayXRefs": [],
                        },
                    },
                    "locations": [{"start": 5, "end": 40}],
                },
            ],
        },
        {
            "md5": "ccc",
            "xref": [{"id": "WP_999999999.1", "name": "WP_999999999.1 hypothetical"}],
            "matches": [],
        },
    ]
}


def test_parse_keys_by_wp_accession():
    calls = parse_interproscan_json(SAMPLE)
    assert set(calls) == {"WP_002805854.1", "WP_002805169.1", "WP_999999999.1"}


def test_parse_flattens_matches_and_aggregates():
    calls = parse_interproscan_json(SAMPLE)
    rbcl = calls["WP_002805854.1"]
    assert rbcl["md5"] == "aaa"
    assert rbcl["match_count"] == 2
    # sorted by (start, evalue, accession): GENE3D (start 15) before PFAM (start 20)
    assert [m["library"] for m in rbcl["matches"]] == ["GENE3D", "PFAM"]
    assert rbcl["interpro_entries"] == ["IPR000685"]  # GENE3D hit is un-integrated
    assert rbcl["go_terms"] == ["GO:0016984"]
    assert rbcl["pathways"] == ["KEGG:00710"]
    assert rbcl["libraries"] == ["GENE3D", "PFAM"]


def test_unintegrated_match_has_null_interpro_fields():
    calls = parse_interproscan_json(SAMPLE)
    g3d = next(m for m in calls["WP_002805854.1"]["matches"] if m["library"] == "GENE3D")
    assert g3d["interpro_accession"] is None
    assert g3d["interpro_description"] is None
    assert g3d["interpro_type"] is None
    # xrefs are not stored per match — they hang off the entry (see entry_xrefs)
    assert "go_terms" not in g3d
    assert "pathways" not in g3d


def test_pattern_hit_has_null_evalue_and_score():
    calls = parse_interproscan_json(SAMPLE)
    m = calls["WP_002805169.1"]["matches"][0]
    assert m["evalue"] is None
    assert m["score"] is None
    assert m["interpro_accession"] == "IPR000454"


def test_zero_match_protein_is_kept_as_sentinel():
    calls = parse_interproscan_json(SAMPLE)
    sentinel = calls["WP_999999999.1"]
    assert sentinel["match_count"] == 0
    assert sentinel["matches"] == []
    assert sentinel["interpro_entries"] == []


def test_multiple_xrefs_fan_out_to_each_accession():
    doc = {"results": [{
        "md5": "d", "xref": [{"id": "WP_A"}, {"id": "WP_B"}],
        "matches": [{
            "signature": {"accession": "PF1", "description": "x",
                          "signatureLibraryRelease": {"library": "PFAM"}, "entry": None},
            "locations": [{"start": 1, "end": 10, "evalue": 1e-5, "score": 20.0}],
        }],
    }]}
    calls = parse_interproscan_json(doc)
    assert calls["WP_A"]["match_count"] == 1
    assert calls["WP_B"]["match_count"] == 1


def test_empty_document():
    assert parse_interproscan_json({}) == {}
    assert parse_interproscan_json({"results": []}) == {}


def test_entry_xrefs_side_table_is_lossless_join():
    """Per-entry GO/pathway detail is recoverable from `interpro_accession`."""
    xrefs = parse_entry_xrefs(SAMPLE)
    # Only entries carrying at least one xref appear (IPR000454 has none).
    assert set(xrefs) == {"IPR000685"}
    assert xrefs["IPR000685"] == {"go_terms": ["GO:0016984"], "pathways": ["KEGG:00710"]}
    # Join back: rbcL's PFAM match points at the entry that holds the xrefs.
    calls = parse_interproscan_json(SAMPLE)
    pfam = next(m for m in calls["WP_002805854.1"]["matches"] if m["library"] == "PFAM")
    assert xrefs[pfam["interpro_accession"]]["go_terms"] == calls["WP_002805854.1"]["go_terms"]


def test_reactome_species_projections_collapse_to_stable_id():
    """InterPro lists every species projection of one curated Reactome event;
    for marine bacteria the species is noise, so only the stable id is kept."""
    assert normalize_pathway_xref("Reactome", "R-HSA-73817") == "Reactome:73817"
    assert normalize_pathway_xref("Reactome", "R-DME-73817") == "Reactome:73817"
    assert normalize_pathway_xref("Reactome", "R-MTU-2408557") == "Reactome:2408557"
    # Non-Reactome databases pass through untouched.
    assert normalize_pathway_xref("MetaCyc", "PWY-6349") == "MetaCyc:PWY-6349"
    assert normalize_pathway_xref("KEGG", "00710") == "KEGG:00710"
    # Anything not matching the species pattern is left alone rather than mangled.
    assert normalize_pathway_xref("Reactome", "73817") == "Reactome:73817"


def test_reactome_collapse_dedups_protein_level_pathways():
    doc = {
        "results": [{
            "md5": "a", "xref": [{"id": "WP_1.1", "name": "n"}],
            "matches": [{
                "signature": {
                    "accession": "PF1", "name": "n", "description": "d",
                    "signatureLibraryRelease": {"library": "PFAM", "version": "1"},
                    "entry": {
                        "accession": "IPR1", "description": "d", "type": "FAMILY",
                        "goXRefs": [],
                        "pathwayXRefs": [
                            {"id": f"R-{sp}-73817", "databaseName": "Reactome"}
                            for sp in ("HSA", "MMU", "DME", "CEL", "BTA")
                        ],
                    },
                },
                "locations": [{"start": 1, "end": 9, "evalue": 1e-9, "score": 1.0}],
            }],
        }]
    }
    calls = parse_interproscan_json(doc)
    assert calls["WP_1.1"]["pathways"] == ["Reactome:73817"]  # 5 species → 1 id


def test_summarize_qc_fields():
    calls = parse_interproscan_json(SAMPLE)
    s = summarize(
        calls, strain="MED4", input_proteins=3, tool_version="5.78-109.0",
        applications="ALL_DEFAULT", image_digest="sha256:x", wallclock_s=12.34,
    )
    assert s["strain"] == "MED4"
    assert s["calls_made"] == 3
    assert s["proteins_no_match"] == 1
    assert s["total_matches"] == 3
    assert s["interpro_integrated_matches"] == 2  # rbcL PFAM + atpH pattern
    assert s["distribution"] == {"GENE3D": 1, "PFAM": 1, "PROSITE_PATTERNS": 1}
    assert s["sentinel_rate"] == round(1 / 3, 4)
    assert s["parse_failures"] == 0
    assert s["image_digest"] == "sha256:x"
    assert s["wallclock_s"] == 12.3


def test_summarize_omits_optional_fields_when_absent():
    s = summarize({}, strain="X", input_proteins=0, tool_version="v",
                  applications="ALL_DEFAULT")
    assert "image_digest" not in s
    assert "wallclock_s" not in s
    assert "xrefs_requested" not in s
    assert s["sentinel_rate"] == 0.0


def test_summarize_xref_coverage():
    """GO/pathway coverage counters — these are what prove --goterms/--pathways
    took effect; a run without them yields entries with empty xref arrays and
    therefore zeros across all four counters."""
    calls = parse_interproscan_json(SAMPLE)
    s = summarize(calls, strain="MED4", input_proteins=3, tool_version="v",
                  applications="ALL_DEFAULT", xrefs_requested=True)
    # Only rbcL's IPR000685 carries xrefs in SAMPLE; atpH's entry has empty ones.
    assert s["proteins_with_go_terms"] == 1
    assert s["distinct_go_terms"] == 1
    assert s["proteins_with_pathways"] == 1
    assert s["distinct_pathways"] == 1
    assert s["pathway_databases"] == {"KEGG": 1}
    assert s["xrefs_requested"] is True


def test_summarize_xref_counters_zero_without_goterms():
    """A --no-xrefs run: every entry has empty goXRefs/pathwayXRefs."""
    stripped = {
        "results": [
            {
                "md5": "aaa",
                "xref": [{"id": "WP_1.1", "name": "WP_1.1"}],
                "matches": [{
                    "signature": {
                        "accession": "PF00016", "name": "n", "description": "d",
                        "signatureLibraryRelease": {"library": "PFAM", "version": "37.0"},
                        "entry": {"accession": "IPR000685", "description": "d",
                                  "type": "DOMAIN", "goXRefs": [], "pathwayXRefs": []},
                    },
                    "locations": [{"start": 1, "end": 9, "evalue": 1e-9, "score": 1.0}],
                }],
            }
        ]
    }
    s = summarize(parse_interproscan_json(stripped), strain="X", input_proteins=1,
                  tool_version="v", applications="ALL_DEFAULT", xrefs_requested=False)
    assert s["proteins_with_go_terms"] == 0
    assert s["distinct_go_terms"] == 0
    assert s["proteins_with_pathways"] == 0
    assert s["distinct_pathways"] == 0
    assert s["pathway_databases"] == {}
    assert s["xrefs_requested"] is False
