"""Unit tests for multiomics_kg/vocab/non_de_evidence.py.

Covers:
- COMPARTMENTS: membership of the 5 canonical values
- EXTENDED_OMICS_TYPES: includes existing 5 + PAIRED_RNASEQ_PROTEOME
- VALUE_KINDS: the 3-value enum the adapter dispatches on
- KNOWN_METRIC_TYPES: dict[str, str] mapping each seeded metric_type to its
  value_kind. The registry's ONLY job is the stable metric_type → value_kind
  contract; everything else (unit, rankable, has_p_value, allowed_categories)
  lives on paperconfig entries.
- DEFAULT_SKIP_TOKENS, VALID_BLANK_POLICIES, BUCKET_THRESHOLD_* constants
"""
import pytest

from multiomics_kg.vocab.non_de_evidence import (
    COMPARTMENTS,
    EXTENDED_OMICS_TYPES,
    VALUE_KINDS,
    KNOWN_METRIC_TYPES,
    DEFAULT_SKIP_TOKENS,
    VALID_BLANK_POLICIES,
    BUCKET_THRESHOLD_TOP_DECILE,
    BUCKET_THRESHOLD_TOP_QUARTILE,
    BUCKET_THRESHOLD_MID,
)


def test_compartments_has_expected_values():
    # `extracellular` added in Phase 2 metabolomics (66bca38) — general
    # outside-the-cell pool, distinct from the proteomics-specific `exoproteome`.
    assert COMPARTMENTS == {
        "whole_cell", "vesicle", "exoproteome", "extracellular",
        "spent_medium", "lysate",
    }


def test_extended_omics_types_extends_existing_set():
    assert {"RNASEQ", "MICROARRAY", "PROTEOMICS", "EXOPROTEOMICS", "METABOLOMICS"} <= EXTENDED_OMICS_TYPES
    assert "PAIRED_RNASEQ_PROTEOME" in EXTENDED_OMICS_TYPES
    assert {"DNASEQ", "VESICLE_PROTEOMICS", "VESICLE_DNASEQ"} <= EXTENDED_OMICS_TYPES


def test_value_kinds_enum():
    assert VALUE_KINDS == {"numeric", "boolean", "categorical"}


def test_known_metric_types_biller_2018_entries():
    assert KNOWN_METRIC_TYPES["periodic_in_axenic_LD"] == "boolean"
    assert KNOWN_METRIC_TYPES["periodic_in_coculture_LD"] == "boolean"
    assert KNOWN_METRIC_TYPES["periodic_in_axenic_extended_darkness"] == "boolean"
    assert KNOWN_METRIC_TYPES["periodic_in_coculture_extended_darkness"] == "boolean"
    assert KNOWN_METRIC_TYPES["darkness_survival_class"] == "categorical"


def test_known_metric_types_numeric_backlog_entries():
    """Numeric metric_type names seeded for upcoming zinser 2009 / Waldbauer 2012 retrofits.
    Registered so a future paper using them can't accidentally clash with a boolean/categorical
    variant, but all other metadata (unit, rankable, has_p_value, p_value_threshold) is declared
    inline on those future paperconfigs when they land."""
    for name in ("fourier_score", "peak_time_h", "peak_fit_r_squared",
                 "protein_transcript_lag_h", "damping_ratio", "diel_amplitude"):
        assert KNOWN_METRIC_TYPES[name] == "numeric", name


def test_known_metric_types_value_kinds_are_valid():
    """Every value in the registry must be a member of VALUE_KINDS."""
    for name, kind in KNOWN_METRIC_TYPES.items():
        assert kind in VALUE_KINDS, f"{name} → {kind!r}"


def test_known_metric_types_is_a_plain_dict_no_dataclass():
    """Registry is deliberately a flat dict, not a dataclass container.
    Per-metric metadata lives on paperconfig entries."""
    assert isinstance(KNOWN_METRIC_TYPES, dict)
    for v in KNOWN_METRIC_TYPES.values():
        assert isinstance(v, str)


def test_token_defaults_and_bucket_thresholds():
    assert DEFAULT_SKIP_TOKENS == ("NA", "N/A", "n/a", "#N/A")
    assert VALID_BLANK_POLICIES == ("skip", "true", "false")
    assert BUCKET_THRESHOLD_TOP_DECILE == 90
    assert BUCKET_THRESHOLD_TOP_QUARTILE == 75
    assert BUCKET_THRESHOLD_MID == 25
