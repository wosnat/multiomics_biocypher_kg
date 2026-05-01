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
    "PROTEOMICS",              # whole-cell default
    "EXOPROTEOMICS",           # compartment baked in: exoproteome (Kaur 2018, Oleza 2017)
    "VESICLE_PROTEOMICS",      # compartment baked in: vesicle (Biller 2014, Biller 2022)
    "METABOLOMICS",
    "PAIRED_RNASEQ_PROTEOME",  # Waldbauer 2012 et al.
    "DNASEQ",                  # whole-cell DNA-seq default
    "VESICLE_DNASEQ",          # compartment baked in: vesicle (Biller 2014 S4)
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

    # ── Numeric, vesicle compartment proteomics (Biller 2022 S3 + S4) ──
    "cell_abundance_biovolume_normalized":    "numeric",
    "vesicle_abundance_biovolume_normalized": "numeric",
    "log2_vesicle_cell_enrichment":           "numeric",

    # ── Vesicle compartment, presence + quality + localization (Biller 2014 S2 + S3) ──
    "vesicle_proteome_member":              "boolean",
    "mascot_identification_probability":    "numeric",
    "predicted_subcellular_localization":   "categorical",

    # ── Numeric, vesicle DNA-seq read coverage (Biller 2014 S4) ──
    "vesicle_dna_avg_read_coverage":        "numeric",

    # ── Boolean, KEGG/topic enrichment set membership (Hennon 2017 S5) ──
    "enrichment_bacterial_chemotaxis":       "boolean",
    "enrichment_biosynthesis_of_amino_acids": "boolean",
    "enrichment_carbon_metabolism":           "boolean",
    "enrichment_fatty_acid_metabolism":       "boolean",
    "enrichment_flagellar_assembly":          "boolean",
    "enrichment_redox_genes":                 "boolean",
    "enrichment_ribosome":                    "boolean",
    "enrichment_tonb_associated_genes":       "boolean",

    # ── Numeric, replicate-count proteomics detection per compartment (Lu 2025 S1) ──
    # 0-N integer = number of biological replicates in which the protein was
    # detected by LC-MS/MS in the named compartment. Not rankable (mass ties on
    # a small ordinal) and no p-value. Compartment is encoded by the parent
    # Experiment (exoproteome / whole_cell), so distinct metric_types per
    # compartment let queries filter by detection compartment directly.
    "exoproteome_detection_replicates":  "numeric",
    "whole_cell_detection_replicates":   "numeric",
}
