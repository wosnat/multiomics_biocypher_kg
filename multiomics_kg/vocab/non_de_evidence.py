# multiomics_kg/vocab/non_de_evidence.py
"""Controlled vocabularies for non-DE-evidence schema (DerivedMetric today;
AbundanceAnalysis in a future slice).

Single source of truth imported by:
- scripts/validate_paperconfig.py
- multiomics_kg/adapters/observations_adapter.py   (Plan 2)
- multiomics_kg/download/build_gene_id_mapping.py
- multiomics_kg/download/resolve_paper_ids.py
- tests/

Design principle (option C): the vocab is deliberately narrow. The ONLY
per-metric fact enforced centrally is value_kind (which drives adapter
edge-type dispatch: numeric → quantifies_gene, boolean → flags_gene,
categorical → classifies_gene). Every other per-metric datum — unit,
rankable, has_p_value, p_value_threshold, allowed_categories,
field_description — is declared inline on paperconfig entries. This keeps
the central registry stable as new papers land and per-paper quirks don't
leak into the vocabulary.
"""
from __future__ import annotations


# ─── Compartment vocabulary (Experiment.compartment) ──────────────────────────

COMPARTMENTS: frozenset[str] = frozenset({
    "whole_cell",       # intracellular (default)
    "vesicle",          # extracellular vesicle fraction
    "exoproteome",      # secreted proteins in medium
    "spent_medium",     # culture supernatant
    "lysate",           # cell lysate
})


# ─── Omics-type vocabulary (extends existing VALID_TYPES) ──────────────────────

EXTENDED_OMICS_TYPES: frozenset[str] = frozenset({
    "RNASEQ",
    "MICROARRAY",
    "PROTEOMICS",
    "EXOPROTEOMICS",
    "METABOLOMICS",
    "PAIRED_RNASEQ_PROTEOME",  # Waldbauer 2012 et al.
})


# ─── value_kind enum — adapter edge-type discriminator ─────────────────────────

VALUE_KINDS: frozenset[str] = frozenset({"numeric", "boolean", "categorical"})


# ─── Token-parsing defaults for boolean derived_metrics_table entries ──────────

# Literal CSV cell values that mean "not tested" (no edge emitted).
DEFAULT_SKIP_TOKENS: tuple[str, ...] = ("NA", "N/A", "n/a", "#N/A")

# Allowed values of the `blank_policy` paperconfig field.
VALID_BLANK_POLICIES: tuple[str, ...] = ("skip", "true", "false")


# ─── Percentile cutoffs pinned by parent spec ──────────────────────────────────

BUCKET_THRESHOLD_TOP_DECILE: int = 90    # percentile >= 90 → "top_decile"
BUCKET_THRESHOLD_TOP_QUARTILE: int = 75  # 75 <= percentile < 90 → "top_quartile"
BUCKET_THRESHOLD_MID: int = 25           # 25 <= percentile < 75 → "mid", else "low"


# ─── KNOWN_METRIC_TYPES registry (filled in Task 3) ────────────────────────────
# Maps metric_type → value_kind. Nothing else. A paperconfig that declares
# a metric_type in this registry must use the matching value_kind; a
# metric_type absent from the registry is accepted with a validator warning
# (authors may introduce new names; the registry grows slowly and only
# records the one thing future papers must agree on).

KNOWN_METRIC_TYPES: dict[str, str] = {
    # ── Numeric (backlog papers: zinser 2009, Waldbauer 2012) ──
    # Registered so a future paper using one of these names can't silently
    # re-declare it as boolean/categorical. All other metadata (unit, rankable,
    # has_p_value, p_value_threshold) is declared inline on those paperconfigs.
    "fourier_score":            "numeric",
    "peak_time_h":               "numeric",
    "peak_fit_r_squared":        "numeric",
    "protein_transcript_lag_h":  "numeric",
    "damping_ratio":             "numeric",
    "diel_amplitude":            "numeric",
    # Paired-modality variants (Waldbauer 2012 reports peak-time and amplitude
    # separately for transcript and protein from one time course; the generics
    # above stay for single-modality papers like zinser 2009).
    "peak_time_transcript_h":       "numeric",
    "peak_time_protein_h":          "numeric",
    "diel_amplitude_transcript_log2": "numeric",
    "diel_amplitude_protein_log2":    "numeric",

    # ── Boolean (Biller 2018 S4A + S4B) ──
    "periodic_in_axenic_LD":                  "boolean",
    "periodic_in_coculture_LD":                "boolean",
    "periodic_in_axenic_extended_darkness":    "boolean",
    "periodic_in_coculture_extended_darkness": "boolean",

    # ── Categorical (Biller 2018 S5) ──
    # The paperconfig entry that uses this metric_type declares its
    # `allowed_categories` inline (Task 11). The registry only locks the
    # value_kind — class vocabularies are per-paper.
    "darkness_survival_class": "categorical",
}
