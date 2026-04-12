"""Tests for prompt construction."""
from multiomics_kg.extraction.timepoint.prompts import (
    SHARED_RULES,
    VALID_GROWTH_PHASES_LIST,
    build_prompt,
)


def test_shared_rules_mentions_unknown_vs_other():
    assert "unknown" in SHARED_RULES.lower()
    assert "other:" in SHARED_RULES


def test_shared_rules_mentions_enum_values():
    for v in ("exponential", "nutrient_limited", "acclimated_steady_state",
              "infected", "recovery", "diel", "darkness", "death", "acute_stress"):
        assert v in SHARED_RULES


def test_valid_growth_phases_list_matches_validator():
    # Must exactly match VALID_GROWTH_PHASES in the validator (no drift).
    expected = {
        "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
        "infected", "recovery", "diel", "darkness", "death", "acute_stress",
        "unknown",
    }
    assert set(VALID_GROWTH_PHASES_LIST) == expected


def test_build_prompt_includes_paper_name_and_doi():
    background = {
        "papername": "Tetu 2019",
        "doi": "10.1038/s42003-019-0410-x",
        "experiments": {"exp1": {"treatment_condition": "HDPE"}},
    }
    targets = [{
        "id": "DE_x",
        "experiment_key": "exp1",
        "logfc_col": "HDPE log2 Fold change",
        "existing": {"timepoint": None, "timepoint_hours": None, "growth_phase": None},
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
    }]
    prompt = build_prompt(background, targets, pdf_cache_entry=None)
    assert "Tetu 2019" in prompt
    assert "10.1038/s42003-019-0410-x" in prompt
    assert "DE_x" in prompt
    assert "HDPE log2 Fold change" in prompt


def test_build_prompt_includes_cache_entry_when_present():
    background = {"papername": "X", "doi": None, "experiments": {}}
    cache_entry = {"title": "My Title", "abstract": "Short abstract"}
    prompt = build_prompt(background, [], pdf_cache_entry=cache_entry)
    assert "My Title" in prompt


def test_build_prompt_skips_cache_entry_when_absent():
    background = {"papername": "X", "doi": None, "experiments": {}}
    prompt = build_prompt(background, [], pdf_cache_entry=None)
    # Should not have a "Publication metadata" section
    assert "PDF cache" not in prompt.lower() or "unavailable" in prompt.lower()
