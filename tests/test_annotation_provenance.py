"""Unit tests for edge-level provenance/confidence (design 2026-08-10 §5.3)."""

import pytest

from multiomics_kg.utils.annotation_provenance import annotation_edge_props


def _gene():
    return {
        "go_terms": ["GO:1", "GO:2", "GO:3"],
        "go_terms_source": {
            "GO:1": ["eggnog", "uniprot", "interproscan"],  # multi-source, curated
            "GO:2": ["interproscan"],                      # interproscan-only
            "GO:3": ["cyanorak"],                        # curated, single
        },
        "go_terms_evidence": {"GO:2": "domain_inferred"},  # sparse: only interproscan-touched
        "pfam_ids": ["PF1"],
        "pfam_ids_source": {"PF1": ["eggnog", "interproscan"]},  # dependent pair
        "pfam_ids_evidence": {"PF1": "curated"},
    }


def test_curated_multisource_high_score():
    p = annotation_edge_props(_gene(), "go_terms", "GO:1")
    assert p["sources"] == ["eggnog", "uniprot", "interproscan"]
    assert p["evidence"] == "curated"
    assert p["evidence_score"] == 1.0          # >=2 indep + curated + not domain; 3/3


def test_domain_inferred_low_score():
    p = annotation_edge_props(_gene(), "go_terms", "GO:2")
    assert p["evidence"] == "domain_inferred"
    assert p["evidence_score"] == 0.0          # single source, not curated, domain; 0/3


def test_curated_single_source_defaults_evidence():
    p = annotation_edge_props(_gene(), "go_terms", "GO:3")
    assert p["evidence"] == "curated"        # no evidence entry → default curated
    assert p["evidence_score"] == 0.667          # curated + not domain, but 1 source; 2/3


def test_pfam_eggnog_interpro_not_independent():
    """eggNOG-Pfam and InterPro-Pfam are the same signal — no +1 for corroboration."""
    p = annotation_edge_props(_gene(), "pfam_ids", "PF1")
    # curated (+1) + not domain (+1) = 2/3; the dependent pair does NOT add the corroboration +1
    assert p["evidence_score"] == 0.667


def test_missing_token_is_curated_no_sources():
    p = annotation_edge_props({"go_terms": ["GO:9"]}, "go_terms", "GO:9")
    assert p["evidence"] == "curated"
    assert "sources" not in p
    assert p["evidence_score"] == 0.667


def test_source_values_are_all_data_source_ids():
    """R2: every declared sources value must be an id in gene_annotations_config."""
    import yaml
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    cfg = yaml.safe_load(open("config/gene_annotations_config.yaml"))
    ds_ids = {ls["id"]
              for src in cfg["sources"].values()
              for ls in src["logical_sources"]}
    for entry in load_vocabularies().values():
        if entry.property == "sources":
            undeclared = set(entry.values) - ds_ids
            assert not undeclared, (
                f"{entry.id} declares source value(s) {sorted(undeclared)} with "
                f"no matching DataSource id. Known ids: {sorted(ds_ids)}")


def test_interpro_source_label_is_interproscan():
    gene = {"go_terms_source": {"GO:1": ["interproscan"]},
            "go_terms_evidence": {"GO:1": "domain_inferred"}}
    props = annotation_edge_props(gene, "go_terms", "GO:1")
    assert props["sources"] == ["interproscan"]


@pytest.mark.parametrize("sources,evidence,expected", [
    (["uniprot", "eggnog"], "curated",         1.0),    # 3/3
    (["eggnog"],            "curated",         0.667),  # 2/3
    (["eggnog"],            "family_inferred", 0.333),  # 1/3
    (["interproscan"],      "domain_inferred", 0.0),    # 0/3
])
def test_evidence_score_is_normalized_to_unit_interval(sources, evidence, expected):
    gene = {"go_terms_source": {"GO:1": sources},
            "go_terms_evidence": {"GO:1": evidence}}
    assert annotation_edge_props(gene, "go_terms", "GO:1")["evidence_score"] == expected


def test_score_is_rounded_to_three_decimals():
    gene = {"go_terms_source": {"GO:1": ["eggnog"]},
            "go_terms_evidence": {"GO:1": "family_inferred"}}
    score = annotation_edge_props(gene, "go_terms", "GO:1")["evidence_score"]
    assert isinstance(score, float)
    assert score == round(score, 3)


def test_round_recovers_the_raw_signal_count():
    """KG-IPT-012: round, never truncate — 0.333 * 3 = 0.999.

    Exercises the actual recovery path against a real emitted score, not
    Python's `round` in isolation: single source + family_inferred fires
    exactly one of the three signals ("not domain_inferred"), so
    evidence_score = 1/3 = 0.333, and round(score * signal_count) must
    recover that raw count of 1.
    """
    gene = {"go_terms_source": {"GO:1": ["eggnog"]},
            "go_terms_evidence": {"GO:1": "family_inferred"}}
    score = annotation_edge_props(gene, "go_terms", "GO:1")["evidence_score"]
    assert score == 0.333
    assert round(score * 3) == 1
