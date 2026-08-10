"""Unit tests for edge-level provenance/confidence (design 2026-08-10 §5.3)."""

from multiomics_kg.utils.annotation_provenance import annotation_edge_props


def _gene():
    return {
        "go_terms": ["GO:1", "GO:2", "GO:3"],
        "go_terms_source": {
            "GO:1": ["eggnog", "uniprot", "interpro"],  # multi-source, curated
            "GO:2": ["interpro"],                        # interpro-only
            "GO:3": ["cyanorak"],                        # curated, single
        },
        "go_terms_evidence": {"GO:2": "domain_inferred"},  # sparse: only interpro-touched
        "pfam_ids": ["PF1"],
        "pfam_ids_source": {"PF1": ["eggnog", "interpro"]},  # dependent pair
        "pfam_ids_evidence": {"PF1": "curated"},
    }


def test_curated_multisource_high_score():
    p = annotation_edge_props(_gene(), "go_terms", "GO:1")
    assert p["sources"] == ["eggnog", "uniprot", "interpro"]
    assert p["evidence"] == "curated"
    assert p["evidence_score"] == 3          # >=2 indep + curated + not domain


def test_domain_inferred_low_score():
    p = annotation_edge_props(_gene(), "go_terms", "GO:2")
    assert p["evidence"] == "domain_inferred"
    assert p["evidence_score"] == 0          # single source, not curated, domain


def test_curated_single_source_defaults_evidence():
    p = annotation_edge_props(_gene(), "go_terms", "GO:3")
    assert p["evidence"] == "curated"        # no evidence entry → default curated
    assert p["evidence_score"] == 2          # curated + not domain, but 1 source


def test_pfam_eggnog_interpro_not_independent():
    """eggNOG-Pfam and InterPro-Pfam are the same signal — no +1 for corroboration."""
    p = annotation_edge_props(_gene(), "pfam_ids", "PF1")
    # curated (+1) + not domain (+1) = 2; the dependent pair does NOT add the corroboration +1
    assert p["evidence_score"] == 2


def test_missing_token_is_curated_no_sources():
    p = annotation_edge_props({"go_terms": ["GO:9"]}, "go_terms", "GO:9")
    assert p["evidence"] == "curated"
    assert "sources" not in p
    assert p["evidence_score"] == 2
