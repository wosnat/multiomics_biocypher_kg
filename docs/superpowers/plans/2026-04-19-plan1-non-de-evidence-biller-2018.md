# Plan 1 — Non-DE Evidence: Vocab, Schema, Paperconfig Preprocessing, Biller 2018 Authoring

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the foundation for the Biller 2018 non-DE-evidence slice — a vocabulary module, schema additions (`DerivedMetric` node + 6 new edge types + `Experiment.compartment`), `derived_metrics_table` support in the paperconfig validator / gene-ID preprocessing / paper-ID pre-resolution, and the Biller 2018 paperconfig authoring for Tables S4A (split across axenic + coculture entries), S4B, and S5 — **4** `derived_metrics_table` entries total. The new paperconfig entries are no-ops in the KG until Plan 2 wires an adapter; the full KG build must still succeed unchanged vs the Step 0 baseline.

**Architecture:** A new `multiomics_kg/vocab/` package owns controlled vocabularies — `COMPARTMENTS`, `EXTENDED_OMICS_TYPES`, `VALUE_KINDS` (the 3-value `{numeric, boolean, categorical}` enum the adapter dispatches on), and `KNOWN_METRIC_TYPES: dict[str, str]` (metric_type → value_kind — stable type-to-edge-shape contract only; a deliberately narrow registry, **not** a per-paper metadata catalogue). All other per-metric metadata (`rankable`, `has_p_value`, `p_value_threshold`, `unit`, `allowed_categories`, `field_description`) is declared inline **on the paperconfig entry** per metric, so one paper's idiosyncratic `periodic_in_axenic_LD` flag doesn't force global registry churn. Validator cross-checks paperconfig `value_kind` against `KNOWN_METRIC_TYPES[metric_type]` when the metric_type is known (hard error on mismatch) and warns on novel metric_types (sanity nudge, not a block). The existing paperconfig pipeline — `paperconfig_utils.py` → `validate_paperconfig.py` → `build_gene_id_mapping.py` → `resolve_paper_ids.py` — is extended to recognize a fourth supplementary-materials type `derived_metrics_table` alongside `csv` / `id_translation` / `annotation_gff` / `gene_clusters`. The validator's monolithic `validate()` is refactored to expose a pure `validate_paperconfig_content(config, path) -> (errors, warnings)` function for unit testing, while the CLI entry keeps its current prints + exit codes. Biller 2018's paperconfig is retrofitted with 3 new entries; the 3 source CSVs (Tables S4A, S4B, S5) are read as-is — no content edits, no `_modified.csv`.

**Tech Stack:** Python 3.12 (uv), pandas, PyYAML, pytest, BioCypher schema YAML.

---

## Scope boundary

**In scope:** vocab module, schema additions, paperconfig_utils iterator, validator refactor + per-`value_kind` dispatch for `derived_metrics_table`, `build_gene_id_mapping` extractor, `resolve_paper_ids` resolver, Biller 2018 paperconfig authoring for S4A/S4B/S5, preprocessing + full-paperconfig regression, test-mode KG build smoke check.

**Out of scope (Plan 2):** the `observations_adapter.py` module, its wiring into `create_knowledge_graph.py`, the synthetic numeric-DM test fixture, emission of `DerivedMetric` nodes / `derived_metric_{flags,classifies,quantifies}_gene` edges.

**Out of scope (Plan 3):** post-import Cypher (rank/percentile/bucket/significance/rollups), `snapshot_data.json` regeneration, KG validity assertions, the `docs/kg-changes/non-de-evidence-extension.md` deliverable, skill-doc updates.

**Never in this slice:** `AbundanceAnalysis`, `abundance_table` / `timeseries_abundance_table`, `VALUE_TYPES`, abundance post-import passes, abundance rollups on Experiment/Publication/OrganismTaxon.

---

## File Structure

**New files**

- `multiomics_kg/vocab/__init__.py` — package marker; re-exports `non_de_evidence` symbols at package level for short imports.
- `multiomics_kg/vocab/non_de_evidence.py` — controlled vocabularies: `COMPARTMENTS`, `EXTENDED_OMICS_TYPES`, `VALUE_KINDS`, `KNOWN_METRIC_TYPES: dict[str, str]`, `DEFAULT_SKIP_TOKENS`, `VALID_BLANK_POLICIES`, `BUCKET_THRESHOLD_TOP_DECILE`, `BUCKET_THRESHOLD_TOP_QUARTILE`, `BUCKET_THRESHOLD_MID`. **No `MetricSpec` dataclass.** Per-metric metadata (`rankable`, `has_p_value`, `p_value_threshold`, `unit`, `allowed_categories`, `field_description`) lives on paperconfig entries, not in central vocab.
- `tests/test_vocab_non_de_evidence.py` — unit tests for the vocab module.

**Modified files (high-level responsibilities)**

- `config/schema_config.yaml` — add `DerivedMetric` node block, 3 binding edges (`publication_has_derived_metric`, `experiment_has_derived_metric`, `derived_metric_belongs_to_organism`), 3 measurement edges (`derived_metric_quantifies_gene`, `derived_metric_flags_gene`, `derived_metric_classifies_gene`), and `compartment: str` property on the `experiment:` block.
- `multiomics_kg/utils/paperconfig_utils.py` — add `iter_derived_metrics_tables(config)` iterator.
- `tests/test_paperconfig_utils.py` — add tests for `iter_derived_metrics_tables`.
- `scripts/validate_paperconfig.py` — refactor: extract `validate_paperconfig_content(config, path) -> (errors, warnings)`; add vocabulary checks (`compartment`, `PAIRED_RNASEQ_PROTEOME`); relax `control_condition` / `test_type` to warnings for experiments whose only supplementary entries are `derived_metrics_table` or `gene_clusters`; dispatch on `type: derived_metrics_table` with a per-`value_kind` branch (numeric / boolean / categorical), including CSV dry-runs for boolean (hard error) and categorical (warning).
- `tests/test_paperconfig_validation.py` — add unit tests for each `value_kind` path, vocabulary enforcement, DE-field relaxation, and the refactored function entry.
- `multiomics_kg/download/build_gene_id_mapping.py` — add `extract_rows_from_derived_metrics_table(entry, paper_name, table_key)` alongside the other extractors; add a dispatcher case in `process_strain()`.
- `multiomics_kg/download/resolve_paper_ids.py` — add `resolve_derived_metrics_entry(pub_name, table_key, table_config, force=False)` and dispatch in `main()` before the generic `resolve_table`.
- `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml` — append 3 `derived_metrics_table` entries (S4A NATL2A periodicity — 4 boolean metrics; S4B MIT1002 periodicity — 2 boolean metrics; S5 NATL2A darkness survival — 1 categorical metric).

**Not touched in this plan** (reserved for Plan 2 / Plan 3)

- `multiomics_kg/adapters/` (no adapter yet)
- `create_knowledge_graph.py` (no new adapter to wire)
- `scripts/post-import.sh`, `scripts/post-import.cypher`
- `tests/kg_validity/` (no DerivedMetric edges in the graph yet)
- `docs/kg-changes/`
- Skill docs

---

## Canonical references baked into this plan

### 1. `KNOWN_METRIC_TYPES` seed entries (drives Task 3)

The central registry is **only** a `metric_type → value_kind` map. All other per-metric metadata is declared inline on paperconfig entries (Task 11). The seed set covers what Plan 1 + upcoming Plan 2 fixture + backlog papers reference; new names can be added in later slices without widening the data structure.

| `metric_type` | `value_kind` | Used by (current or near-term) |
|---|---|---|
| `fourier_score` | `numeric` | zinser 2009 (future retrofit) |
| `peak_time_h` | `numeric` | zinser 2009 |
| `peak_fit_r_squared` | `numeric` | zinser 2009 |
| `protein_transcript_lag_h` | `numeric` | Waldbauer 2012 |
| `damping_ratio` | `numeric` | Waldbauer 2012 |
| `diel_amplitude` | `numeric` | Waldbauer 2012 |
| `periodic_in_axenic_LD` | `boolean` | Biller 2018 S4A |
| `periodic_in_coculture_LD` | `boolean` | Biller 2018 S4A, S4B |
| `periodic_in_axenic_extended_darkness` | `boolean` | Biller 2018 S4A |
| `periodic_in_coculture_extended_darkness` | `boolean` | Biller 2018 S4A, S4B |
| `darkness_survival_class` | `categorical` | Biller 2018 S5 |

### 2. Per-value_kind paperconfig field rules (drives Task 8)

Every metric entry under `metrics:` in a `derived_metrics_table` declares the following, and the validator enforces all three (existence, gating, vocabulary):

| Field | `numeric` | `boolean` | `categorical` |
|---|---|---|---|
| `metric_type` | required | required | required |
| `value_kind` | required; `"numeric"` | required; `"boolean"` | required; `"categorical"` |
| `value_col` | required; must be a CSV column | required | required |
| `rankable` | required; `"true"` or `"false"` | **forbidden** (adapter sets `"false"` at ingest) | **forbidden** (adapter sets `"false"` at ingest) |
| `has_p_value` | required; `"true"` or `"false"` | **forbidden** (adapter sets `"false"` at ingest) | **forbidden** (adapter sets `"false"` at ingest) |
| `p_value_threshold` | required iff `has_p_value == "true"`; float ∈ (0, 1] | forbidden | forbidden |
| `p_value_col` | at least one of `p_value_col` / `adjusted_p_value_col` required iff `has_p_value == "true"`; forbidden otherwise | forbidden | forbidden |
| `adjusted_p_value_col` | as above | forbidden | forbidden |
| `unit` | optional (free string or omitted; omit for dimensionless) | forbidden | forbidden |
| `true_tokens` | forbidden | **required; non-empty list of literal cell strings** | forbidden |
| `false_tokens` | forbidden | optional; default `[]` | forbidden |
| `skip_tokens` | forbidden | optional; default `DEFAULT_SKIP_TOKENS` | forbidden |
| `blank_policy` | forbidden | optional; default `"skip"`; one of `VALID_BLANK_POLICIES` | forbidden |
| `allowed_categories` | forbidden | forbidden | **required; non-empty list** |
| `field_description` | **required** (non-empty string) | **required** (non-empty string) | **required** (non-empty string) |

**Cross-check against registry:**
- If `metric_type ∈ KNOWN_METRIC_TYPES` and the paperconfig's declared `value_kind != KNOWN_METRIC_TYPES[metric_type]` → **hard error** (a metric_type's edge type can't change).
- If `metric_type ∉ KNOWN_METRIC_TYPES` → **warning** ("novel metric_type; consider adding to multiomics_kg/vocab/non_de_evidence.py if it will be reused") — still accepted.

**Why `field_description` is required on every metric:** it's the human-readable explanation of what the metric measures (e.g., `"RAIN periodicity FDR<0.05 in NATL2A axenic L:D (Table S4A)"`). Free text — no vocabulary constraint, no length limit beyond non-empty. Downstream MCP tools surface it verbatim so users understand what a numeric / boolean / categorical value represents; a missing field turns the `DerivedMetric` node into an opaque ID. Required for all three value_kinds.

**Why `allowed_categories` is declared on the paperconfig, not in vocab:** the vocabulary of classes for a categorical metric is an artifact of what the paper actually wrote in its table. For Biller 2018 S5, the 3 distinct strings are exactly what's in the CSV column `darkness_cluster` (verified: 100 / 90 / 79 rows). Another paper's `darkness_survival_class` might use entirely different labels; forcing them through one global set would either restrict or destabilize the registry.

### 3. Biller 2018 post-Step-0 paperconfig state

`data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml` already contains (as of commit 25f0cd8):
- 3 experiments: `darkness_extended_darkness_natl2a_rnaseq_axenic`, `darkness_extended_darkness_natl2a_rnaseq_coculture`, `darkness_extended_darkness_mit1002_rnaseq`
- 3 supplementary-material entries other than DE: `annotation_gff_mit1002_old_draft`, `id_translation_mit1002_rast_diamond`, `id_translation_mit1002_rast_names`
- 2 DE csv entries: `supp_table_s3`, `supp_table_s6b`
- **no** `gene_clusters` entries (removed in Step 0)

### 4. Source CSV column names (Task 11 depends on these exact strings)

- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4.csv`
  - Columns (post-header): `RAST annotation gene number`, `NCBI ID`, `NCBI ID_2`, `NCBI ID_3`, `Genbank Annotation`, `RAST annotation`, `Periodic in axenic, L:D cultures`, `Periodic in co-cultured, L:D cultures`, `Periodic in axenic, extended darkness cultures`, `Periodic in co-cultured, extended darkness cultures`, `periodicity_cluster`
  - Gene ID column for NATL2A locus tags: `NCBI ID_2` (id_type `locus_tag_ncbi`, PMN2A_RS* format)
  - Each Y/blank column: `{"Y": N, blank: M}` — no other tokens
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4.csv`
  - Columns: `RAST annotation gene number`, `Locus ID`, `RAST RegionID (Chromosome_start_stop)`, `RAST Annotation`, `Periodic in co-cultured, L:D cultures`, `Periodic in co-cultured, extended darkness cultures`, `periodicity_cluster`
  - Gene ID column for MIT1002: `Locus ID` (id_type `locus_tag`, MIT1002_NNNN with trailing whitespace — strip on read)
  - Each Y/blank column: Y / blank
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5.csv`
  - Columns: `NCBI ID`, `NCBI ID_2`, `NCBI ID_3`, `Gene Name`, `Annotation`, `Consistently in extended darkness axenic cultures, 72-144 hours`, `Consistently in extended darkness co-cultures, 72-144 hours`, `Unique to axenic cultures`, `Unique to co-cultures`, `darkness_cluster`
  - Gene ID column: `NCBI ID_2` (id_type `locus_tag_ncbi`, PMN2A_RS*)
  - `darkness_cluster`: exactly 3 values seen (listed in table above)

---

## Tasks

### Task 1: Create `multiomics_kg/vocab/` package skeleton + failing test

**Files:**
- Create: `multiomics_kg/vocab/__init__.py`
- Create: `multiomics_kg/vocab/non_de_evidence.py`
- Create: `tests/test_vocab_non_de_evidence.py`

- [ ] **Step 1: Write the failing test file**

```python
# tests/test_vocab_non_de_evidence.py
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
    assert COMPARTMENTS == {"whole_cell", "vesicle", "exoproteome", "spent_medium", "lysate"}


def test_extended_omics_types_extends_existing_set():
    assert {"RNASEQ", "MICROARRAY", "PROTEOMICS", "EXOPROTEOMICS", "METABOLOMICS"} <= EXTENDED_OMICS_TYPES
    assert "PAIRED_RNASEQ_PROTEOME" in EXTENDED_OMICS_TYPES


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
```

- [ ] **Step 2: Run test — expect ModuleNotFoundError**

Run: `uv run pytest tests/test_vocab_non_de_evidence.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'multiomics_kg.vocab'`.

- [ ] **Step 3: Create the package marker**

```python
# multiomics_kg/vocab/__init__.py
"""Controlled vocabularies shared across the knowledge-graph pipeline.

Single source of truth for string enumerations so the validator, the
adapter (Plan 2), preprocessing scripts, and tests can't drift.
"""
```

- [ ] **Step 4: Create a stub `non_de_evidence.py` that makes import succeed**

```python
# multiomics_kg/vocab/non_de_evidence.py
"""Controlled vocabularies for non-DE-evidence schema (DerivedMetric today;
AbundanceAnalysis in a future slice).

Design note: the registry is deliberately narrow. Only value_kind (which
drives adapter edge-type dispatch) is enforced centrally. Per-metric
metadata — rankable, has_p_value, p_value_threshold, unit,
allowed_categories, field_description — is declared inline on paperconfig
entries, so one paper's idiosyncratic metric doesn't force global vocab churn.
"""
from __future__ import annotations

COMPARTMENTS: frozenset[str] = frozenset()
EXTENDED_OMICS_TYPES: frozenset[str] = frozenset()
VALUE_KINDS: frozenset[str] = frozenset()
KNOWN_METRIC_TYPES: dict[str, str] = {}

DEFAULT_SKIP_TOKENS: tuple[str, ...] = ()
VALID_BLANK_POLICIES: tuple[str, ...] = ()

BUCKET_THRESHOLD_TOP_DECILE: int = 0
BUCKET_THRESHOLD_TOP_QUARTILE: int = 0
BUCKET_THRESHOLD_MID: int = 0
```

- [ ] **Step 5: Verify test now fails on value assertions, not import**

Run: `uv run pytest tests/test_vocab_non_de_evidence.py -v`
Expected: test collection succeeds; test bodies fail with `AssertionError` (empty vocabs).

- [ ] **Step 6: Commit (red tests are OK at this checkpoint)**

```bash
git add multiomics_kg/vocab/__init__.py multiomics_kg/vocab/non_de_evidence.py tests/test_vocab_non_de_evidence.py
git commit -m "vocab: scaffold non_de_evidence module with KNOWN_METRIC_TYPES + failing tests"
```

---

### Task 2: Fill in `COMPARTMENTS`, `EXTENDED_OMICS_TYPES`, token/bucket constants

**Files:**
- Modify: `multiomics_kg/vocab/non_de_evidence.py`

- [ ] **Step 1: Populate the flat vocabularies + bucket thresholds + VALUE_KINDS**

Replace the body of `multiomics_kg/vocab/non_de_evidence.py` with:

```python
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

KNOWN_METRIC_TYPES: dict[str, str] = {}   # filled in Task 3
```

- [ ] **Step 2: Re-run tests**

Run: `uv run pytest tests/test_vocab_non_de_evidence.py -v`
Expected: `test_compartments_has_expected_values`, `test_extended_omics_types_extends_existing_set`, `test_value_kinds_enum`, `test_known_metric_types_is_a_plain_dict_no_dataclass`, `test_token_defaults_and_bucket_thresholds` PASS; `test_known_metric_types_biller_2018_entries`, `test_known_metric_types_numeric_backlog_entries`, `test_known_metric_types_value_kinds_are_valid` still FAIL (empty registry).

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/vocab/non_de_evidence.py
git commit -m "vocab: populate COMPARTMENTS, EXTENDED_OMICS_TYPES, VALUE_KINDS, token + bucket constants"
```

---

### Task 3: Seed `KNOWN_METRIC_TYPES` registry

**Files:**
- Modify: `multiomics_kg/vocab/non_de_evidence.py`

- [ ] **Step 1: Append the registry**

Replace the line `KNOWN_METRIC_TYPES: dict[str, str] = {}   # filled in Task 3` in `multiomics_kg/vocab/non_de_evidence.py` with:

```python
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
```

- [ ] **Step 2: Run tests — all vocab tests should pass**

Run: `uv run pytest tests/test_vocab_non_de_evidence.py -v`
Expected: all 7 test functions PASS.

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/vocab/non_de_evidence.py
git commit -m "vocab: seed KNOWN_METRIC_TYPES with Biller 2018 + backlog numeric names"
```

---

### Task 4: Add `iter_derived_metrics_tables` to `paperconfig_utils.py`

**Files:**
- Modify: `multiomics_kg/utils/paperconfig_utils.py`
- Modify: `tests/test_paperconfig_utils.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_paperconfig_utils.py`:

```python
# ─── iter_derived_metrics_tables ────────────────────────────────────────


def test_iter_derived_metrics_tables():
    """iter_derived_metrics_tables yields (key, config) for derived_metrics_table entries only."""
    from multiomics_kg.utils.paperconfig_utils import iter_derived_metrics_tables

    config = {
        "publication": {
            "supplementary_materials": {
                "supp_table_1": {"type": "csv", "filename": "data.csv"},
                "cluster_table_1": {"type": "gene_clusters", "filename": "c.csv"},
                "dm_s4a": {
                    "type": "derived_metrics_table",
                    "filename": "s4a.csv",
                    "organism": "Prochlorococcus NATL2A",
                    "experiment": "darkness_extended_darkness_natl2a_rnaseq_axenic",
                    "name_col": "NCBI ID_2",
                    "metrics": [
                        {"metric_type": "periodic_in_axenic_LD",
                         "value_kind": "boolean",
                         "value_col": "Periodic in axenic, L:D cultures",
                         "true_tokens": ["Y"]},
                    ],
                },
                "id_trans": {"type": "id_translation", "filename": "ids.csv"},
            }
        }
    }

    results = list(iter_derived_metrics_tables(config))
    assert len(results) == 1
    key, table = results[0]
    assert key == "dm_s4a"
    assert table["type"] == "derived_metrics_table"


def test_iter_derived_metrics_tables_empty():
    from multiomics_kg.utils.paperconfig_utils import iter_derived_metrics_tables

    assert list(iter_derived_metrics_tables({})) == []
    assert list(iter_derived_metrics_tables({"publication": {"supplementary_materials": {}}})) == []
    # Purely DE paperconfig → no derived_metrics_table entries
    config = {"publication": {"supplementary_materials": {
        "supp_1": {"type": "csv", "filename": "d.csv"},
    }}}
    assert list(iter_derived_metrics_tables(config)) == []
```

- [ ] **Step 2: Run test — expect ImportError**

Run: `uv run pytest tests/test_paperconfig_utils.py::test_iter_derived_metrics_tables tests/test_paperconfig_utils.py::test_iter_derived_metrics_tables_empty -v`
Expected: FAIL with `ImportError: cannot import name 'iter_derived_metrics_tables'`.

- [ ] **Step 3: Add the iterator**

In `multiomics_kg/utils/paperconfig_utils.py`, add the following function directly after `iter_cluster_tables` (around line 95):

```python
def iter_derived_metrics_tables(config: dict):
    """Yield (table_key, table_config) for derived_metrics_table-type supplementary tables."""
    supp = get_supplementary_materials(config)
    for key, table in supp.items():
        if not isinstance(table, dict):
            continue
        if table.get("type") == "derived_metrics_table":
            yield key, table
```

- [ ] **Step 4: Run test — expect PASS**

Run: `uv run pytest tests/test_paperconfig_utils.py::test_iter_derived_metrics_tables tests/test_paperconfig_utils.py::test_iter_derived_metrics_tables_empty -v`
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/paperconfig_utils.py tests/test_paperconfig_utils.py
git commit -m "paperconfig_utils: add iter_derived_metrics_tables iterator"
```

---

### Task 5: Schema — add `DerivedMetric` node, 6 new edges, `Experiment.compartment`

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add `compartment` to the Experiment node properties**

Open `config/schema_config.yaml`. Under the `experiment:` block (starts around line 33), find the `properties:` map and add `compartment: str` immediately after `light_intensity: str`. The resulting snippet in place of what's currently at `config/schema_config.yaml:53`:

```yaml
    light_intensity: str
    compartment: str                 # adapter-emitted; default "whole_cell"; vocab in multiomics_kg/vocab/non_de_evidence.py
    table_scope: str
    table_scope_detail: str
```

- [ ] **Step 2: Add the `DerivedMetric` node definition**

Insert the following block **after** the existing `clustering analysis:` node block (which ends near `config/schema_config.yaml:104`) and **before** `# Publication → Experiment`:

```yaml
# DerivedMetric node: one per (Experiment × metric_type). Carries column-level
# metric semantics. Edge type emitted to Gene depends on value_kind.
# Adapter emits in Plan 2; post-import computes total_gene_count + growth_phases in Plan 3.
derived metric:
  is_a: information content entity
  represented_as: node
  preferred_id: derived_metric_id
  label_in_input: derived_metric
  properties:
    name: str
    experiment_id: str
    organism_name: str
    publication_doi: str
    compartment: str
    omics_type: str
    treatment_type: str[]
    background_factors: str[]
    treatment: str
    light_condition: str
    experimental_context: str
    metric_type: str
    value_kind: str                  # "numeric" | "boolean" | "categorical"
    unit: str
    rankable: str                    # "true" | "false" (string enum)
    has_p_value: str                 # "true" | "false" (string enum)
    p_value_threshold: float
    allowed_categories: str[]
    field_description: str
    total_gene_count: int            # post-import (Plan 3)
    growth_phases: str[]             # post-import (Plan 3)
```

- [ ] **Step 3: Add the 3 binding edges**

Insert **after** the `experiment has clustering analysis:` block and **before** `# GeneCluster → Gene membership edges`:

```yaml
# Publication → DerivedMetric (mirrors Publication_has_clustering_analysis)
publication has derived metric:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: publication_has_derived_metric
  source: publication
  target: derived metric

# Experiment → DerivedMetric (mirrors Experiment_has_clustering_analysis)
experiment has derived metric:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: experiment_has_derived_metric
  source: experiment
  target: derived metric

# DerivedMetric → OrganismTaxon (mirrors Clusteringanalysis_belongs_to_organism)
derived metric belongs to organism:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: derived_metric_belongs_to_organism
  source: derived metric
  target: organism taxon
```

- [ ] **Step 4: Add the 3 measurement edges (DerivedMetric → Gene)**

Insert directly after the binding edges from Step 3:

```yaml
# DerivedMetric → Gene (numeric metric_type; value_kind="numeric")
derived metric quantifies gene:
  is_a: Association
  represented_as: edge
  label_as_edge: derived_metric_quantifies_gene
  label_in_input: derived_metric_quantifies_gene
  source: derived metric
  target: gene
  properties:
    metric_type: str
    value: float
    p_value: float
    adjusted_p_value: float
    rank_by_metric: int              # post-import (Plan 3; rankable="true" only)
    metric_percentile: float         # post-import (Plan 3; rankable="true" only)
    metric_bucket: str               # post-import (Plan 3; rankable="true" only)
    significant: str                 # post-import (Plan 3; has_p_value="true" only)

# DerivedMetric → Gene (boolean metric_type; value_kind="boolean")
derived metric flags gene:
  is_a: Association
  represented_as: edge
  label_as_edge: derived_metric_flags_gene
  label_in_input: derived_metric_flags_gene
  source: derived metric
  target: gene
  properties:
    metric_type: str
    value_flag: str                  # "true" | "false" (string enum)

# DerivedMetric → Gene (categorical metric_type; value_kind="categorical")
derived metric classifies gene:
  is_a: Association
  represented_as: edge
  label_as_edge: derived_metric_classifies_gene
  label_in_input: derived_metric_classifies_gene
  source: derived metric
  target: gene
  properties:
    metric_type: str
    value_text: str                  # must be in parent DerivedMetric.allowed_categories
```

- [ ] **Step 5: Verify schema loads (no adapter emits to it yet)**

Run: `uv run pytest tests/test_schema_config.py -v`
Expected: PASS (this suite validates YAML well-formedness and existing entries; it does not require adapter emissions).

- [ ] **Step 6: Verify test build succeeds unchanged**

Run: `uv run python create_knowledge_graph.py --test 2>&1 | tail -30`
Expected: completes cleanly, no errors about unknown node labels. No new CSVs appear under the output directory for `DerivedMetric` because no adapter emits yet.

- [ ] **Step 7: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add DerivedMetric node, 6 binding+measurement edges, Experiment.compartment"
```

---

### Task 6: Refactor validator — extract `validate_paperconfig_content(config, path)`

**Files:**
- Modify: `scripts/validate_paperconfig.py`

Goal: The current `validate(config_path)` in `scripts/validate_paperconfig.py` is a ~400-line CLI entry that both reads the YAML and emits `print()` side effects. Extract the content-validation logic into a pure `validate_paperconfig_content(config: dict, path: str) -> tuple[list[str], list[str]]` function that returns `(errors, warnings)`. Keep `validate(config_path)` as a thin wrapper that loads YAML, calls `validate_paperconfig_content`, and handles the CLI side effects.

- [ ] **Step 1: Write failing test**

Append to `tests/test_paperconfig_validation.py`:

```python
# ---------------------------------------------------------------------------
# validate_paperconfig_content — pure function entry
# ---------------------------------------------------------------------------

def test_validate_paperconfig_content_is_importable():
    """The refactored pure function must be exported from validate_paperconfig."""
    from validate_paperconfig import validate_paperconfig_content

    errors, warnings = validate_paperconfig_content({}, "/tmp/fake.yaml")
    assert isinstance(errors, list)
    assert isinstance(warnings, list)


def test_validate_paperconfig_content_accepts_minimal_valid_config(tmp_path):
    """A minimal well-formed config returns empty errors (warnings OK)."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "paper.pdf"
    pdf.write_bytes(b"%PDF-1.4\n")  # placeholder; validator only checks existence
    config = {
        "publication": {
            "papername": "Test 2026",
            "papermainpdf": str(pdf),
            "experiments": {
                "exp1": {
                    "name": "MED4 N-starvation",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                    "treatment_condition": "N-limited",
                    "control_condition": "N-replete",
                },
            },
            "supplementary_materials": {
                "t1": {
                    "type": "csv",
                    "filename": str(csv_file),
                    "statistical_analyses": [
                        {
                            "id": "an1",
                            "experiment": "exp1",
                            "name_col": "Gene",
                            "logfc_col": "log2FC",
                            "adjusted_p_value_col": "padj",
                            "timepoint": "unknown",
                            "timepoint_hours": None,
                            "growth_phase": "exponential",
                        },
                    ],
                },
            },
        },
    }
    cfg_path = _write_config(tmp_path, config)
    errors, _warnings = validate_paperconfig_content(config, str(cfg_path))
    assert errors == [], f"Expected no errors; got: {errors}"
```

- [ ] **Step 2: Run test — expect ImportError**

Run: `uv run pytest tests/test_paperconfig_validation.py::test_validate_paperconfig_content_is_importable -v`
Expected: FAIL.

- [ ] **Step 3: Refactor `scripts/validate_paperconfig.py`**

At the end of the imports section (around `scripts/validate_paperconfig.py:88`), leave all constants in place. Replace the current `def validate(config_path: str) -> bool:` signature with two functions. The body that currently sits inside `validate()` after the `load_paperconfig` call is moved verbatim into `validate_paperconfig_content`, except that:
- every `errors.append(...)` / `warnings.append(...)` already uses the two local lists (keep as-is)
- the final `_print_results(errors, warnings)` call is moved OUT of the content function back into `validate()` (so prints stay in the CLI wrapper)
- `validate_paperconfig_content` returns `(errors, warnings)`

Concretely:

```python
def validate_paperconfig_content(
    config: dict, config_path: str,
) -> tuple[list[str], list[str]]:
    """Validate a parsed paperconfig dict; returns (errors, warnings).

    Pure function — no stdout prints about results, no SystemExit. Callers
    render the diagnostics they want.

    NOTE: the current implementation still prints progress lines via
    print(); those are debug output, not results, and remain for parity
    with the previous CLI behavior when invoked via validate().
    """
    errors: list[str] = []
    warnings: list[str] = []

    abs_config = os.path.abspath(config_path)
    _project_root = os.path.dirname(abs_config)
    for _ in range(10):
        candidate = os.path.join(
            _project_root, "data", "Prochlorococcus", "treatment_organisms.csv"
        )
        if os.path.exists(candidate):
            break
        _project_root = os.path.dirname(_project_root)
    else:
        _project_root = None

    treatment_organisms = _load_treatment_organisms(_project_root or "")
    all_canonical_organisms = CANONICAL_GENOMIC_ORGANISMS | treatment_organisms
    genome_accessions = _load_genome_accessions(_project_root or "")

    # ---- everything the old validate() did between the try/except YAML parse
    # and the final _print_results call goes here VERBATIM, using `config`
    # instead of reloading YAML ----
    # (paste existing body from after "config = load_paperconfig(Path(config_path))"
    # through the final "_validate_organism_consistency(...)" call, minus the
    # _print_results line.)

    return errors, warnings


def validate(config_path: str) -> bool:
    """CLI entry: parse YAML, validate, print results, return success bool."""
    try:
        config = load_paperconfig(Path(config_path))
    except Exception as e:
        print(f"FAIL: Cannot parse YAML: {e}")
        return False

    errors, warnings = validate_paperconfig_content(config, config_path)
    _print_results(errors, warnings)
    return len(errors) == 0
```

- [ ] **Step 4: Manually copy the body of the old `validate()` — precise range + early-exit conversion**

Copy from the current `validate()` function into `validate_paperconfig_content` with these exact boundaries:

**Do NOT copy (already re-computed at the top of `validate_paperconfig_content` by Step 3):**
- The `abs_config = os.path.abspath(config_path)` block through the four `_project_root` / `treatment_organisms` / `all_canonical_organisms` / `genome_accessions` assignments (current file lines ~789–806). These are now at the top of `validate_paperconfig_content`.
- The YAML parse `try: config = load_paperconfig(Path(config_path))` / `except: ...` block (lines ~809–813). The CLI `validate()` wrapper owns YAML parsing.
- The final `_print_results(errors, warnings)` call and the `return len(errors) == 0` line. The CLI wrapper owns these.

**DO copy** (around lines 816 through 1178 of the current file) — starting from:
```python
is_resource_config = bool(config) and "publication" not in config
```
through (inclusive):
```python
    # --- Organism consistency check ---
    if experiments and supp:
        _validate_organism_consistency(
            experiments, supp, config_path, errors, warnings,
        )
```

**Convert early-exit `return False` paths to append-and-return-tuple.** The moved body has 2 `return False` paths that need conversion so `validate_paperconfig_content` returns `(errors, warnings)` consistently:

1. Current line ~824: `if not config or "publication" not in config: print("FAIL: Missing top-level 'publication' key"); return False`
   → Replace with: `if not config or "publication" not in config: errors.append(f"{config_path} | Missing top-level 'publication' key"); return errors, warnings`

2. Current line ~896: inside `if not supp: errors.append("Missing 'supplementary_materials'"); _print_results(...); return len(errors) == 0`
   → Replace with: `if not supp: errors.append(f"{config_path} | Missing 'supplementary_materials'"); return errors, warnings` (drop the `_print_results` call; the CLI wrapper prints).

Every other `errors.append(...)` / `warnings.append(...)` in the moved body stays as-is — they already use the two local lists and don't early-return.

At the very end of `validate_paperconfig_content`, after the copied body, add:
```python
    return errors, warnings
```

- [ ] **Step 5: Run the full existing validation test suite to verify no regressions**

Run: `uv run pytest tests/test_paperconfig_validation.py -v`
Expected: all previously green tests still green; the two new tests (Step 1) now PASS.

- [ ] **Step 6: Run the real-paper regression**

Run: `uv run python scripts/validate_paperconfig.py --all 2>&1 | tail -10`
Expected: `Results: N/N passed, 0/N failed` (same N as before the refactor).

- [ ] **Step 7: Commit**

```bash
git add scripts/validate_paperconfig.py tests/test_paperconfig_validation.py
git commit -m "validate_paperconfig: extract pure validate_paperconfig_content(config, path)"
```

---

### Task 7: Validator — add `compartment` vocab check + `PAIRED_RNASEQ_PROTEOME` + DE-field relaxation

**Files:**
- Modify: `scripts/validate_paperconfig.py`
- Modify: `tests/test_paperconfig_validation.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_paperconfig_validation.py`:

```python
# ---------------------------------------------------------------------------
# compartment vocabulary check (Task 7)
# ---------------------------------------------------------------------------

def test_validate_rejects_unknown_compartment(tmp_path):
    """Experiment.compartment must be in COMPARTMENTS."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X",
            "papermainpdf": str(pdf),
            "experiments": {
                "exp1": {
                    "name": "e", "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ", "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                    "treatment_condition": "A", "control_condition": "B",
                    "compartment": "nucleus",   # ← not in vocabulary
                },
            },
            "supplementary_materials": {
                "t1": {"type": "csv", "filename": str(csv_file),
                       "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                                 "name_col": "Gene", "logfc_col": "log2FC",
                                                 "adjusted_p_value_col": "padj",
                                                 "timepoint": "unknown", "timepoint_hours": None,
                                                 "growth_phase": "exponential"}]},
            },
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("compartment" in e and "nucleus" in e for e in errors), errors


def test_validate_accepts_whole_cell_compartment(tmp_path):
    """'whole_cell' is the default; explicit 'whole_cell' must be accepted."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp1": {
                "name": "e", "organism": "Prochlorococcus MED4",
                "omics_type": "RNASEQ", "test_type": "DESeq2",
                "treatment_type": ["nitrogen"],
                "treatment_condition": "A", "control_condition": "B",
                "compartment": "whole_cell",
            }},
            "supplementary_materials": {"t1": {
                "type": "csv", "filename": str(csv_file),
                "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                          "name_col": "Gene", "logfc_col": "log2FC",
                                          "adjusted_p_value_col": "padj",
                                          "timepoint": "unknown", "timepoint_hours": None,
                                          "growth_phase": "exponential"}],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert not any("compartment" in e for e in errors), errors


def test_validate_accepts_paired_rnaseq_proteome_omics_type(tmp_path):
    """PAIRED_RNASEQ_PROTEOME is a canonical omics_type after this slice lands."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp1": {
                "name": "e", "organism": "Prochlorococcus MED4",
                "omics_type": "PAIRED_RNASEQ_PROTEOME",
                "test_type": "DESeq2", "treatment_type": ["diel"],
                "treatment_condition": "A", "control_condition": "B",
            }},
            "supplementary_materials": {"t1": {
                "type": "csv", "filename": str(csv_file),
                "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                          "name_col": "Gene", "logfc_col": "log2FC",
                                          "adjusted_p_value_col": "padj",
                                          "timepoint": "unknown", "timepoint_hours": None,
                                          "growth_phase": "exponential"}],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert not any("omics_type" in e for e in errors), errors


def test_validate_relaxes_de_fields_for_derived_metrics_only_experiment(tmp_path):
    """Experiments whose only supplementary entries are derived_metrics_table may
    omit control_condition and test_type (warn, don't error)."""
    from validate_paperconfig import validate_paperconfig_content

    # A minimal CSV with a Y/blank boolean column, referenced by the DM entry.
    csv_file = tmp_path / "dm.csv"
    csv_file.write_text("locus_tag,flag_col\nPMN2A_RS00015,Y\nPMN2A_RS00020,\n")
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp_dm_only": {
                "name": "e", "organism": "Prochlorococcus NATL2A",
                "omics_type": "RNASEQ",
                "treatment_type": ["darkness"],
                "treatment_condition": "A",
                # ← control_condition + test_type deliberately omitted
            }},
            "supplementary_materials": {"dm1": {
                "type": "derived_metrics_table",
                "filename": str(csv_file),
                "organism": "Prochlorococcus NATL2A",
                "experiment": "exp_dm_only",
                "name_col": "locus_tag",
                "id_columns": [{"column": "locus_tag", "id_type": "locus_tag_ncbi"}],
                "metrics": [{
                    "metric_type": "periodic_in_axenic_LD",
                    "value_kind": "boolean",
                    "value_col": "flag_col",
                    "true_tokens": ["Y"],
                }],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert not any("control_condition" in e for e in errors), errors
    assert not any("test_type" in e for e in errors), errors
```

- [ ] **Step 2: Run tests — expect 4 failures**

Run: `uv run pytest tests/test_paperconfig_validation.py -k "compartment or paired_rnaseq or relaxes_de_fields" -v`
Expected: the 4 new tests FAIL (vocab unknown; omics unknown; relaxation not wired).

- [ ] **Step 3: Add vocab imports at top of `scripts/validate_paperconfig.py`**

Near the existing imports (around line 45), add:

```python
from multiomics_kg.vocab.non_de_evidence import (
    COMPARTMENTS,
    EXTENDED_OMICS_TYPES,
    VALUE_KINDS,
    KNOWN_METRIC_TYPES,
    DEFAULT_SKIP_TOKENS,
    VALID_BLANK_POLICIES,
)
```

- [ ] **Step 4: Replace `VALID_TYPES` with the extended set**

Find the line `VALID_TYPES = {"RNASEQ", "MICROARRAY", "PROTEOMICS", "EXOPROTEOMICS", "METABOLOMICS"}` (around line 81) and replace with:

```python
# Accepted omics_type values on experiment entries.
# Extended set includes PAIRED_RNASEQ_PROTEOME (Waldbauer 2012 et al.).
VALID_TYPES = set(EXTENDED_OMICS_TYPES)
```

- [ ] **Step 5: Add compartment vocab check in `_validate_experiments`**

Inside `_validate_experiments` (after the `# Canonical treatment_organism (optional...)` block, around line 484), add:

```python
        # Canonical compartment (optional; default "whole_cell" at ingest time)
        compartment = exp.get("compartment", "")
        if compartment and compartment not in COMPARTMENTS:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "compartment", compartment, COMPARTMENTS,
                )
            )
        elif compartment:
            print(f"    compartment '{compartment}': OK")
```

- [ ] **Step 6: Extend DE-field relaxation so derived_metrics_table experiments are DE-field-optional**

In the `_validate_experiments` signature and body, `has_analyses` is currently computed as `referenced_exp_keys and exp_key in referenced_exp_keys` (around line 385). The referenced set is built from `statistical_analyses` only. This already correctly treats `derived_metrics_table`-only experiments as DE-field-optional: they have no `statistical_analyses`, so `exp_key` won't enter `referenced_exp_keys`, and `DE_REQUIRED_EXPERIMENT_FIELDS` become warnings.

The test from Step 1 (`test_validate_relaxes_de_fields_for_derived_metrics_only_experiment`) asserts this. Run it now:

Run: `uv run pytest tests/test_paperconfig_validation.py::test_validate_relaxes_de_fields_for_derived_metrics_only_experiment -v`
Expected: PASS without further edits (relaxation is automatic via the existing `has_analyses` check).

If it fails because some other part of the body hard-errors on `derived_metrics_table` (unknown type), defer the fix to Task 8 — that task adds the dispatcher branch that makes the supplementary_materials loop not reject the new type.

- [ ] **Step 7: Re-run the 4 new tests**

Run: `uv run pytest tests/test_paperconfig_validation.py -k "compartment or paired_rnaseq or relaxes_de_fields" -v`
Expected: `test_validate_rejects_unknown_compartment`, `test_validate_accepts_whole_cell_compartment`, `test_validate_accepts_paired_rnaseq_proteome_omics_type` PASS. `test_validate_relaxes_de_fields_for_derived_metrics_only_experiment` may still fail here — it depends on Task 8.

- [ ] **Step 8: Run the real-paper regression**

Run: `uv run python scripts/validate_paperconfig.py --all 2>&1 | tail -5`
Expected: `Results: N/N passed, 0/N failed` — no existing paperconfig is affected (none currently uses `compartment` or `PAIRED_RNASEQ_PROTEOME`).

- [ ] **Step 9: Commit**

```bash
git add scripts/validate_paperconfig.py tests/test_paperconfig_validation.py
git commit -m "validate_paperconfig: add compartment vocab + PAIRED_RNASEQ_PROTEOME omics_type"
```

---

### Task 8: Validator — dispatch + per-`value_kind` validation for `derived_metrics_table`

**Files:**
- Modify: `scripts/validate_paperconfig.py`
- Modify: `tests/test_paperconfig_validation.py`

- [ ] **Step 1: Write failing tests (one per dispatch path + required-field + registry-sanity)**

Append to `tests/test_paperconfig_validation.py`:

```python
# ---------------------------------------------------------------------------
# derived_metrics_table dispatch (Task 8) — option C: per-metric metadata
# declared inline on paperconfig; KNOWN_METRIC_TYPES only pins value_kind.
# ---------------------------------------------------------------------------

def _dm_boolean_csv(tmp_path: Path) -> Path:
    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,Periodic in axenic L D cultures\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,\n"
        "PMN2A_RS00025,NA\n"
    )
    return csv


def _dm_wrapper_config(tmp_path: Path, dm_entry: dict) -> dict:
    csv_file = _write_minimal_csv(tmp_path)  # unused but keeps structure
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    return {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp_dm_only": {
                "name": "e", "organism": "Prochlorococcus NATL2A",
                "omics_type": "RNASEQ",
                "treatment_type": ["darkness"],
                "treatment_condition": "ED",
            }},
            "supplementary_materials": {"dm1": dm_entry},
        },
    }


def test_validate_boolean_dm_entry_accepts_shape(tmp_path):
    """Minimal well-formed boolean DM entry — no rankable/has_p_value declared
    (adapter sets them to 'false' at ingest since boolean is definitionally
    non-rankable and has no p-value). field_description is required."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table",
        "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "blank_policy": "skip",
            "field_description": "RAIN periodicity FDR<0.05 in axenic L:D",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert errors == [], errors


def test_validate_dm_requires_field_description(tmp_path):
    """field_description is required on every metric (all value_kinds) —
    free-text human-readable explanation surfaced by downstream MCP tools."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            # ← no field_description
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("field_description" in e for e in errors), errors


def test_validate_dm_rejects_empty_field_description(tmp_path):
    """field_description must be non-empty (whitespace-only is also rejected)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "field_description": "   ",   # whitespace-only
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("field_description" in e for e in errors), errors


def test_validate_boolean_dm_rejects_unknown_token(tmp_path):
    """CSV dry-run finds an unclassified token → hard error."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,flag\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,MAYBE\n"   # ← unknown token
    )
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean", "value_col": "flag",
            "true_tokens": ["Y"], "blank_policy": "skip",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("MAYBE" in e for e in errors), errors


def test_validate_boolean_dm_rejects_missing_true_tokens(tmp_path):
    """true_tokens is required (non-empty list) for value_kind=boolean."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            # ← no true_tokens
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("true_tokens" in e for e in errors), errors


def test_validate_boolean_dm_forbids_rankable_and_has_p_value(tmp_path):
    """rankable and has_p_value are forbidden on boolean entries — adapter sets
    them to "false" at ingest since they're definitionally "false" for boolean.
    Declaring them on the paperconfig (even as "false") is redundant and
    rejected to prevent copy-paste noise."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "rankable": "false",       # ← forbidden even when "false"
            "has_p_value": "false",    # ← forbidden even when "false"
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("rankable" in e for e in errors), errors
    assert any("has_p_value" in e for e in errors), errors


def test_validate_categorical_dm_forbids_rankable_and_has_p_value(tmp_path):
    """Same forbid rule as boolean — adapter sets at ingest."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text("NCBI ID_2,darkness_cluster\nPMN2A_RS00015,a\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            "allowed_categories": ["a"],
            "rankable": "false",       # ← forbidden
            "has_p_value": "false",    # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("rankable" in e for e in errors), errors
    assert any("has_p_value" in e for e in errors), errors


def test_validate_boolean_dm_rejects_forbidden_numeric_fields(tmp_path):
    """unit / p_value_col / p_value_threshold are forbidden on boolean entries."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "unit": "h",                  # ← forbidden
            "p_value_threshold": 0.05,    # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("unit" in e for e in errors) and any("p_value_threshold" in e for e in errors), errors


def test_validate_categorical_dm_requires_allowed_categories(tmp_path):
    """value_kind=categorical MUST declare allowed_categories inline (non-empty list)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text("NCBI ID_2,darkness_cluster\nPMN2A_RS00015,darkness_axenic+darkness_coculture\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            # ← no allowed_categories
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("allowed_categories" in e for e in errors), errors


def test_validate_categorical_dm_warns_on_out_of_set_values(tmp_path):
    """CSV dry-run warns if a value_col cell lies outside the paperconfig's
    declared allowed_categories. (Hard error at ingest is the adapter's job
    in Plan 2.)"""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text(
        "NCBI ID_2,darkness_cluster\n"
        "PMN2A_RS00015,darkness_axenic+darkness_coculture\n"
        "PMN2A_RS00020,totally_new_category\n"   # ← not in allowed_categories
    )
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            "allowed_categories": [
                "darkness_axenic+darkness_coculture",
                "darkness_coculture+unique_coculture",
                "darkness_axenic+unique_axenic",
            ],
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert any("totally_new_category" in w for w in warnings), warnings


def test_validate_numeric_dm_requires_value_col_and_rankable(tmp_path):
    """value_kind=numeric requires value_col + rankable + has_p_value."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser_s1.csv"
    csv.write_text("locus_tag,Fourier,Fourier FDR\nPMM0001,0.82,0.01\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "fourier_score",
            "value_kind": "numeric",
            # ← value_col, rankable, has_p_value all missing
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    joined = " | ".join(errors)
    assert "value_col" in joined and "rankable" in joined and "has_p_value" in joined, errors


def test_validate_numeric_dm_has_p_value_true_requires_threshold(tmp_path):
    """When has_p_value='true', p_value_threshold is required and at least one
    of p_value_col / adjusted_p_value_col must be present."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser.csv"
    csv.write_text("locus_tag,Fourier\nPMM0001,0.82\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "fourier_score",
            "value_kind": "numeric",
            "value_col": "Fourier",
            "rankable": "true",
            "has_p_value": "true",
            # ← p_value_threshold and p_value_col/adjusted_p_value_col missing
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("p_value_threshold" in e for e in errors), errors
    assert any("p_value_col" in e or "adjusted_p_value_col" in e for e in errors), errors


def test_validate_numeric_dm_has_p_value_false_forbids_threshold(tmp_path):
    """When has_p_value='false', p_value_threshold is forbidden."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser.csv"
    csv.write_text("locus_tag,Peak\nPMM0001,6.0\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "peak_time_h",
            "value_kind": "numeric",
            "value_col": "Peak",
            "unit": "h",
            "rankable": "false",
            "has_p_value": "false",
            "p_value_threshold": 0.05,  # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("p_value_threshold" in e for e in errors), errors


def test_validate_rejects_value_kind_mismatch_against_registry(tmp_path):
    """KNOWN_METRIC_TYPES pins metric_type → value_kind. A paperconfig that
    re-declares a known metric_type under a different value_kind is a hard
    error (prevents silent edge-type changes across papers)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",   # registry: boolean
            "value_kind": "numeric",                    # ← mismatch
            "value_col": "Periodic in axenic L D cultures",
            "rankable": "true",
            "has_p_value": "false",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("periodic_in_axenic_LD" in e and "value_kind" in e for e in errors), errors


def test_validate_novel_metric_type_warns_not_errors(tmp_path):
    """A metric_type absent from KNOWN_METRIC_TYPES is accepted with a warning —
    authors may introduce new names; the registry grows slowly."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "a_brand_new_per_paper_flag",   # not in registry
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert not any("a_brand_new_per_paper_flag" in e for e in errors), errors
    assert any("a_brand_new_per_paper_flag" in w and ("novel" in w.lower() or "not in" in w.lower())
               for w in warnings), warnings


def test_validate_rejects_unknown_value_kind(tmp_path):
    """value_kind must be a member of VALUE_KINDS (numeric|boolean|categorical)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "ordinal",  # ← not a VALUE_KINDS member
            "value_col": "Periodic in axenic L D cultures",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("value_kind" in e and "ordinal" in e for e in errors), errors


def test_validate_rejects_missing_top_level_fields(tmp_path):
    """Top-level required fields: filename, organism, experiment, name_col, metrics."""
    from validate_paperconfig import validate_paperconfig_content

    dm_entry = {
        "type": "derived_metrics_table",
        # everything else missing
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    joined = " | ".join(errors)
    for req in ("filename", "organism", "experiment", "name_col", "metrics"):
        assert req in joined, f"'{req}' should be flagged as missing; errors={errors}"
```

- [ ] **Step 2: Run tests — expect failures**

Run: `uv run pytest tests/test_paperconfig_validation.py -k "dm_ or derived_metrics" -v`
Expected: ~16 new dm_ tests FAIL (`derived_metrics_table` type is unrecognized; dispatcher not written).

- [ ] **Step 3: Add the `derived_metrics_table` branch to the supplementary_materials loop — BEFORE the filename short-circuits**

In `scripts/validate_paperconfig.py` inside `validate_paperconfig_content` (moved there by Task 6), find the supplementary_materials loop. The current structure is:

```python
for table_key, table in supp.items():
    table_type = table.get("type", "csv")
    print(f"\n  [{table_key}] (type: {table_type})")

    fn = table.get("filename", "")
    if not fn:
        errors.append(f"{table_key}: Missing 'filename'")
        continue
    if not os.path.exists(fn):
        errors.append(f"{table_key}: File not found: {fn}")
        continue

    print(f"    file: {fn}")
    # ... annotation_gff / id_translation / gene_clusters / csv branches ...
```

**Insert the `derived_metrics_table` dispatcher BEFORE the `fn = table.get("filename", "")` line** — right after the `print(f"\n  [{table_key}] (type: {table_type})")` line. This lets `_validate_derived_metrics_entry` produce its full required-field error set (filename / organism / experiment / name_col / metrics all missing) even when `filename` is absent, rather than being short-circuited by the generic "Missing 'filename'" check that would produce only a single error:

```python
    for table_key, table in supp.items():
        table_type = table.get("type", "csv")
        print(f"\n  [{table_key}] (type: {table_type})")

        # ── derived_metrics_table (Plan 1 / non-DE-evidence slice) ──────────
        # Dispatched BEFORE the generic filename short-circuits so the helper's
        # detailed required-field validation runs for every DM entry, including
        # entries with missing/bad filename.
        if table_type == "derived_metrics_table":
            _validate_derived_metrics_entry(
                table_key, table, config_path,
                all_canonical_organisms, errors, warnings,
            )
            continue

        fn = table.get("filename", "")
        if not fn:
            errors.append(f"{table_key}: Missing 'filename'")
            continue
        # ... rest unchanged ...
```

The `_validate_derived_metrics_entry` helper itself handles filename-missing and file-not-found via its own required-field check + CSV-load try/except (see Step 4).

- [ ] **Step 4: Add the `_validate_derived_metrics_entry` helper**

Insert this function at module scope in `scripts/validate_paperconfig.py`, above `validate_paperconfig_content` (after `_validate_gene_clusters_entry`):

```python
def _validate_derived_metrics_entry(
    key: str, table: dict, config_path: str,
    all_canonical_organisms: set,
    errors: list, warnings: list,
) -> None:
    """Validate a type: derived_metrics_table supplementary entry.

    Design (option C): per-metric metadata (rankable, has_p_value,
    p_value_threshold, unit, allowed_categories) lives inline on each metric
    entry — it's paper-specific. KNOWN_METRIC_TYPES only pins metric_type →
    value_kind to catch edge-type drift across papers.

    Validation steps per metric:
      1. required top-level fields present (filename, organism, experiment,
         name_col, metrics)
      2. organism is in canonical list
      3. CSV file readable; name_col + value_col exist as columns
      4. metric_type non-empty; value_kind ∈ VALUE_KINDS
      5. if metric_type ∈ KNOWN_METRIC_TYPES: declared value_kind must match;
         if not: warn (novel name)
      6. field_description: required non-empty string (free text) — surfaced
         by downstream MCP tools
      7. per-value_kind required / forbidden field gating:
           numeric     — require value_col, rankable, has_p_value;
                         rankable/has_p_value ∈ {"true","false"};
                         if has_p_value=="true": require p_value_threshold ∈ (0,1]
                           and at least one of p_value_col/adjusted_p_value_col;
                         if has_p_value=="false": those three forbidden;
                         unit optional (free string or omitted);
                         forbidden: true_tokens/false_tokens/skip_tokens/
                                   blank_policy/allowed_categories
           boolean     — require value_col, true_tokens (non-empty list);
                         optional false_tokens / skip_tokens / blank_policy
                         (∈ VALID_BLANK_POLICIES, default "skip");
                         forbidden: unit, rankable, has_p_value, p_value_*,
                                   allowed_categories
                                   (rankable/has_p_value are definitionally
                                   "false" for boolean — adapter sets them
                                   at ingest, paperconfig must not declare);
                         CSV dry-run: hard-error on any unclassified cell value
           categorical — require value_col, allowed_categories (non-empty list);
                         forbidden: unit, rankable, has_p_value,
                                   true_tokens/false_tokens/skip_tokens/
                                   blank_policy, p_value_*
                                   (rankable/has_p_value are definitionally
                                   "false" for categorical — adapter sets them
                                   at ingest, paperconfig must not declare);
                         CSV dry-run: warn on cells outside allowed_categories
                         (adapter hard-errors at ingest in Plan 2)
    """
    # ── Required top-level fields ──
    for req in ("filename", "organism", "experiment", "name_col", "metrics"):
        if req not in table or table.get(req) in (None, "", []):
            errors.append(
                f"{config_path} | {key} | derived_metrics_table missing required "
                f"field '{req}'"
            )

    # ── organism vocab ──
    organism = table.get("organism", "")
    if organism and organism not in all_canonical_organisms:
        errors.append(
            _canonical_organism_error(config_path, key, organism, all_canonical_organisms)
        )

    # ── experiment reference type-check (cross-reference happens elsewhere) ──
    exp_ref = table.get("experiment", "")
    if exp_ref and not isinstance(exp_ref, str):
        errors.append(f"{config_path} | {key} | 'experiment' must be a string key")

    # ── CSV load for dry-runs ──
    filename = table.get("filename", "")
    sep = table.get("sep", ",")
    skip = int(table.get("skip_rows", 0) or 0)
    df = None
    cols: list[str] = []
    if filename and os.path.exists(filename):
        try:
            df = pd.read_csv(filename, sep=sep, skiprows=skip if skip else None, dtype=str)
            cols = list(df.columns)
            print(f"    columns: {cols}")
        except Exception as e:
            errors.append(f"{config_path} | {key} | cannot read CSV: {e}")
    elif filename:
        errors.append(f"{config_path} | {key} | file not found: {filename}")

    # ── name_col presence ──
    name_col = table.get("name_col", "")
    if name_col and cols and name_col not in cols:
        errors.append(f"{config_path} | {key} | name_col '{name_col}' not in CSV columns")

    # ── id_columns validation (reuse existing helper) ──
    id_columns = table.get("id_columns", []) or []
    if cols:
        _validate_id_columns(id_columns, cols, key, errors, warnings)

    # ── metrics ──
    metrics = table.get("metrics") or []
    if not isinstance(metrics, list) or not metrics:
        errors.append(f"{config_path} | {key} | 'metrics' must be a non-empty list")
        return

    for i, m in enumerate(metrics):
        ctx = f"{key}.metrics[{i}]"
        if not isinstance(m, dict):
            errors.append(f"{config_path} | {ctx} | must be a mapping")
            continue

        metric_type = m.get("metric_type", "")
        value_kind = m.get("value_kind", "")
        value_col = m.get("value_col", "")

        # ── metric_type presence ──
        if not metric_type:
            errors.append(f"{config_path} | {ctx} | missing metric_type")
            continue

        # ── value_kind presence + vocabulary ──
        if not value_kind:
            errors.append(f"{config_path} | {ctx} | missing value_kind")
            continue
        if value_kind not in VALUE_KINDS:
            errors.append(
                f"{config_path} | {ctx} | value_kind '{value_kind}' not in "
                f"{sorted(VALUE_KINDS)}"
            )
            continue

        # ── registry cross-check (known metric_type must match value_kind) ──
        if metric_type in KNOWN_METRIC_TYPES:
            expected_kind = KNOWN_METRIC_TYPES[metric_type]
            if value_kind != expected_kind:
                errors.append(
                    f"{config_path} | {ctx} | metric_type '{metric_type}' is "
                    f"registered as value_kind='{expected_kind}' in KNOWN_METRIC_TYPES; "
                    f"paperconfig declares value_kind='{value_kind}'. Use a different "
                    f"metric_type name if this paper measures a different kind."
                )
                continue
        else:
            warnings.append(
                f"{config_path} | {ctx} | metric_type '{metric_type}' is novel "
                f"(not in KNOWN_METRIC_TYPES). Consider adding it to "
                f"multiomics_kg/vocab/non_de_evidence.py if future papers will reuse it."
            )

        # ── value_col presence + CSV column match ──
        if not value_col:
            errors.append(f"{config_path} | {ctx} | missing value_col")
        elif cols and value_col not in cols:
            errors.append(
                f"{config_path} | {ctx} | value_col '{value_col}' not in CSV columns"
            )

        # ── field_description: required non-empty string (free text) ──
        # Human-readable explanation surfaced by downstream MCP tools;
        # without it a DerivedMetric node is opaque.
        fd = m.get("field_description")
        if not isinstance(fd, str) or not fd.strip():
            errors.append(
                f"{config_path} | {ctx} | 'field_description' is required and must be "
                f"a non-empty string (free text; no vocabulary constraint)"
            )

        # ── per-value_kind gating ──
        if value_kind == "numeric":
            # rankable + has_p_value required and must be "true"|"false"
            for flag_field in ("rankable", "has_p_value"):
                if flag_field not in m:
                    errors.append(
                        f"{config_path} | {ctx} | '{flag_field}' is required for "
                        f"value_kind=numeric (\"true\" or \"false\")"
                    )
                elif m[flag_field] not in ("true", "false"):
                    errors.append(
                        f"{config_path} | {ctx} | '{flag_field}' must be \"true\" "
                        f"or \"false\" (string enum), got {m[flag_field]!r}"
                    )
            # p-value field gating
            has_p = m.get("has_p_value")
            threshold = m.get("p_value_threshold")
            pval_col = m.get("p_value_col")
            adj_col = m.get("adjusted_p_value_col")
            if has_p == "true":
                if threshold is None:
                    errors.append(
                        f"{config_path} | {ctx} | p_value_threshold is required "
                        f"when has_p_value='true'"
                    )
                elif not isinstance(threshold, (int, float)) or not (0 < threshold <= 1):
                    errors.append(
                        f"{config_path} | {ctx} | p_value_threshold must be a "
                        f"number in (0, 1], got {threshold!r}"
                    )
                if not pval_col and not adj_col:
                    errors.append(
                        f"{config_path} | {ctx} | at least one of p_value_col / "
                        f"adjusted_p_value_col is required when has_p_value='true'"
                    )
                for col_field in ("p_value_col", "adjusted_p_value_col"):
                    declared = m.get(col_field)
                    if declared and cols and declared not in cols:
                        errors.append(
                            f"{config_path} | {ctx} | {col_field} '{declared}' "
                            f"not in CSV columns"
                        )
            elif has_p == "false":
                for forbid in ("p_value_threshold", "p_value_col", "adjusted_p_value_col"):
                    if forbid in m:
                        errors.append(
                            f"{config_path} | {ctx} | '{forbid}' is forbidden when "
                            f"has_p_value='false'"
                        )
            # Forbidden-field gating
            for forbid in ("true_tokens", "false_tokens", "skip_tokens",
                           "blank_policy", "allowed_categories"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=numeric"
                    )

        elif value_kind == "boolean":
            # true_tokens required
            true_tokens = m.get("true_tokens")
            if not isinstance(true_tokens, list) or not true_tokens:
                errors.append(
                    f"{config_path} | {ctx} | 'true_tokens' is required and must be "
                    f"a non-empty list for value_kind=boolean"
                )
            # Forbidden-field gating. rankable/has_p_value are definitionally
            # "false" for boolean — adapter sets them at ingest, paperconfig
            # must not declare them (forbids redundant / copy-paste noise).
            for forbid in ("unit", "rankable", "has_p_value",
                           "p_value_col", "adjusted_p_value_col",
                           "p_value_threshold", "allowed_categories"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=boolean"
                    )
            # blank_policy vocab
            blank_policy = m.get("blank_policy", "skip")
            if blank_policy not in VALID_BLANK_POLICIES:
                errors.append(
                    f"{config_path} | {ctx} | blank_policy '{blank_policy}' not in "
                    f"{list(VALID_BLANK_POLICIES)}"
                )
            # CSV dry-run: every non-blank cell must be classified
            if df is not None and value_col and value_col in df.columns:
                true_set = set(m.get("true_tokens") or [])
                false_set = set(m.get("false_tokens") or [])
                skip_set = set(m.get("skip_tokens") or list(DEFAULT_SKIP_TOKENS))
                classified = true_set | false_set | skip_set
                distinct = set()
                for v in df[value_col].fillna("").astype(str).str.strip().tolist():
                    if v == "":
                        continue
                    distinct.add(v)
                unknown = distinct - classified
                if unknown:
                    errors.append(
                        f"{config_path} | {ctx} | unclassified token(s) in value_col "
                        f"'{value_col}': {sorted(unknown)}. Add them to true_tokens, "
                        f"false_tokens, or skip_tokens."
                    )

        elif value_kind == "categorical":
            # allowed_categories required
            allowed_cats = m.get("allowed_categories")
            if not isinstance(allowed_cats, list) or not allowed_cats:
                errors.append(
                    f"{config_path} | {ctx} | 'allowed_categories' is required and "
                    f"must be a non-empty list for value_kind=categorical"
                )
            # Forbidden-field gating. rankable/has_p_value are definitionally
            # "false" for categorical — adapter sets them at ingest, paperconfig
            # must not declare them (forbids redundant / copy-paste noise).
            for forbid in ("unit", "rankable", "has_p_value",
                           "true_tokens", "false_tokens", "skip_tokens",
                           "blank_policy", "p_value_col", "adjusted_p_value_col",
                           "p_value_threshold"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=categorical"
                    )
            # CSV dry-run: warn on cells outside allowed_categories
            if (df is not None and value_col and value_col in df.columns
                    and isinstance(allowed_cats, list) and allowed_cats):
                allowed = set(allowed_cats)
                distinct = set()
                for v in df[value_col].fillna("").astype(str).str.strip().tolist():
                    if v:
                        distinct.add(v)
                unknown = distinct - allowed
                if unknown:
                    warnings.append(
                        f"{config_path} | {ctx} | value_col '{value_col}' contains "
                        f"{len(unknown)} value(s) outside declared allowed_categories: "
                        f"{sorted(unknown)}. Adapter will hard-error on these rows at "
                        f"ingest (Plan 2). Either expand allowed_categories or fix the CSV."
                    )
```

- [ ] **Step 5: Run DM tests**

Run: `uv run pytest tests/test_paperconfig_validation.py -k "dm_ or derived_metrics or relaxes_de_fields" -v`
Expected: all ~16 new DM tests PASS; `test_validate_relaxes_de_fields_for_derived_metrics_only_experiment` also PASSes now.

- [ ] **Step 6: Run real-paper regression**

Run: `uv run python scripts/validate_paperconfig.py --all 2>&1 | tail -5`
Expected: `Results: N/N passed` — existing paperconfigs unchanged.

- [ ] **Step 7: Commit**

```bash
git add scripts/validate_paperconfig.py tests/test_paperconfig_validation.py
git commit -m "validate_paperconfig: dispatch derived_metrics_table per value_kind with CSV dry-runs"
```

---

### Task 9: `build_gene_id_mapping.py` — add `extract_rows_from_derived_metrics_table`

**Files:**
- Modify: `multiomics_kg/download/build_gene_id_mapping.py`
- Create: `tests/test_build_gene_id_mapping_derived_metrics.py`

- [ ] **Step 1: Write failing test**

Create `tests/test_build_gene_id_mapping_derived_metrics.py`:

```python
"""Unit tests for extract_rows_from_derived_metrics_table in build_gene_id_mapping.

A derived_metrics_table entry uses id_columns identically to a csv entry,
but without the statistical_analyses indirection. The extractor harvests
(id_val, id_type) pairs per row for GeneIdGraph seeding.
"""
from pathlib import Path

import pandas as pd

from multiomics_kg.download.build_gene_id_mapping import (
    extract_rows_from_derived_metrics_table,
)


def test_extract_rows_basic(tmp_path, monkeypatch):
    """name_col + id_columns yields one row per non-empty CSV row."""
    from multiomics_kg.download import build_gene_id_mapping as bgm

    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,Gene Name,flag\n"
        "PMN2A_RS00015,dnaN,Y\n"
        "PMN2A_RS00020,,Y\n"
    )
    # Patch PROJECT_ROOT so the extractor's relative-path resolution lands in tmp_path.
    monkeypatch.setattr(bgm, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": csv.name,
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp1",
        "name_col": "NCBI ID_2",
        "id_columns": [
            {"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"},
            {"column": "Gene Name", "id_type": "gene_name"},
        ],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "flag",
            "true_tokens": ["Y"],
        }],
    }
    rows = extract_rows_from_derived_metrics_table(entry, "Biller 2018", "dm_s4a")
    # Row 1: (PMN2A_RS00015, locus_tag_ncbi), (dnaN, gene_name) — gene_name is second-order
    # Row 2: (PMN2A_RS00020, locus_tag_ncbi) — gene_name blank → not emitted
    assert len(rows) == 2
    pairs_row_1 = {p for p in rows[0][0]}
    assert ("PMN2A_RS00015", "locus_tag_ncbi") in pairs_row_1
    assert ("dnaN", "gene_name") in pairs_row_1
    pairs_row_2 = {p for p in rows[1][0]}
    assert ("PMN2A_RS00020", "locus_tag_ncbi") in pairs_row_2
    assert not any(pt == "gene_name" for _, pt in pairs_row_2)


def test_extract_rows_missing_id_columns_returns_empty(tmp_path, monkeypatch):
    from multiomics_kg.download import build_gene_id_mapping as bgm

    csv = tmp_path / "s4a.csv"
    csv.write_text("NCBI ID_2,flag\nPMN2A_RS00015,Y\n")
    monkeypatch.setattr(bgm, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": csv.name,
        "name_col": "NCBI ID_2",
        # no id_columns → extractor returns no useful pairs
    }
    rows = extract_rows_from_derived_metrics_table(entry, "X", "dm_x")
    assert rows == []
```

- [ ] **Step 2: Run test — expect ImportError**

Run: `uv run pytest tests/test_build_gene_id_mapping_derived_metrics.py -v`
Expected: FAIL — `cannot import name 'extract_rows_from_derived_metrics_table'`.

- [ ] **Step 3: Add the extractor**

In `multiomics_kg/download/build_gene_id_mapping.py`, add immediately after `extract_rows_from_gene_clusters_table` (around line 413):

```python
def extract_rows_from_derived_metrics_table(
    entry: dict,
    paper_name: str,
    table_key: str,
) -> list[tuple[list[tuple[str, str]], str]]:
    """Extract rows from a derived_metrics_table entry (id_columns only).

    Structure mirrors extract_rows_from_csv_table, minus the statistical_analyses
    indirection. name_col is at the table level; id_columns live at the table
    level. The `metrics` list is irrelevant to gene-ID mapping.
    """
    id_columns: list[dict] = entry.get("id_columns") or []
    name_col = entry.get("name_col", "")

    if not id_columns and not name_col:
        return []

    filename = entry.get("filename", "")
    sep = entry.get("sep", ",")
    skip_rows = entry.get("skip_rows", 0)

    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] derived_metrics file not found: {path}", file=sys.stderr)
        return []

    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip_rows, dtype=str,
                         encoding="utf-8-sig")
    except Exception as e:
        print(f"    [warn] could not read {path}: {e}", file=sys.stderr)
        return []

    source_label = f"{paper_name}/{table_key}"
    result: list[tuple[list[tuple[str, str]], str]] = []

    # name_col value appears with its declared id_type (defaulting to locus_tag),
    # matching the csv-entry convention.
    name_col_id_type = next(
        (c.get("id_type", "other") for c in id_columns if c.get("column") == name_col),
        "locus_tag",
    )
    extra_cols = [
        cs for cs in id_columns
        if cs.get("column", "") in df.columns and cs.get("column", "") != name_col
    ]

    for _, row in df.iterrows():
        row_pairs: list[tuple[str, str]] = []

        if name_col and name_col in df.columns:
            val = _safe_str(row.get(name_col, ""))
            if val:
                row_pairs.append((val, name_col_id_type))

        for col_spec in extra_cols:
            col = col_spec["column"]
            id_type = col_spec.get("id_type", "other")
            val = _safe_str(row.get(col, ""))
            if val:
                row_pairs.append((val, id_type))

        if row_pairs:
            result.append((row_pairs, source_label))

    return result
```

- [ ] **Step 4: Add dispatcher case in `process_strain()`**

In `process_strain()` in `build_gene_id_mapping.py` (around line 657), find the `elif entry_type == "gene_clusters":` branch. Append immediately after its `print(f"    collected {len(rows)} rows")` line:

```python
            elif entry_type == "derived_metrics_table":
                rows = extract_rows_from_derived_metrics_table(entry_config, paper_name, table_key)
                all_rows.extend(rows)
                print(f"    collected {len(rows)} rows")
```

- [ ] **Step 5: Run tests**

Run: `uv run pytest tests/test_build_gene_id_mapping_derived_metrics.py -v`
Expected: both tests PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_gene_id_mapping.py tests/test_build_gene_id_mapping_derived_metrics.py
git commit -m "build_gene_id_mapping: extract IDs from derived_metrics_table entries"
```

---

### Task 10: `resolve_paper_ids.py` — add `resolve_derived_metrics_entry`

**Files:**
- Modify: `multiomics_kg/download/resolve_paper_ids.py`
- Create: `tests/test_resolve_paper_ids_derived_metrics.py`

- [ ] **Step 1: Write failing test**

Create `tests/test_resolve_paper_ids_derived_metrics.py`:

```python
"""Unit tests for resolve_derived_metrics_entry in resolve_paper_ids.

The resolver writes <stem>_resolved.csv next to the source, adding
a resolved_locus_tag column and a resolution_method column.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from multiomics_kg.download import resolve_paper_ids as rpi
from multiomics_kg.utils.gene_id_utils import MappingData


def _seed_mapping(genome_dir: Path) -> None:
    genome_dir.mkdir(parents=True, exist_ok=True)
    gene_id_mapping = {
        "version": 2,
        "strain": "NATL2A",
        "specific_lookup": {
            "PMN2A_RS00015": "PMN2A_RS00015",
            "PMN2A_RS00020": "PMN2A_RS00020",
        },
        "multi_lookup": {},
        "conflicts": {},
        "genes": {
            "PMN2A_RS00015": {"tier1_ids": ["PMN2A_RS00015"], "tier2_ids": [], "tier3_ids": []},
            "PMN2A_RS00020": {"tier1_ids": ["PMN2A_RS00020"], "tier2_ids": [], "tier3_ids": []},
        },
        "stats": {"n_genes": 2, "n_specific": 2, "n_multi": 0, "n_conflicts": 0, "passes": 1},
    }
    (genome_dir / "gene_id_mapping.json").write_text(json.dumps(gene_id_mapping))


def test_resolve_derived_metrics_entry_writes_resolved_csv(tmp_path, monkeypatch):
    src = tmp_path / "s4a.csv"
    src.write_text(
        "NCBI ID_2,flag\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,\n"
        "UNKNOWN_ID,Y\n"
    )

    genome_dir = tmp_path / "cache/data/Prochlorococcus/genomes/NATL2A"
    _seed_mapping(genome_dir)

    monkeypatch.setattr(
        "multiomics_kg.utils.gene_id_utils.get_genome_dir",
        lambda organism, project_root: str(genome_dir),
    )
    monkeypatch.setattr(rpi, "get_genome_dir",
                        lambda organism, project_root=None: str(genome_dir))
    monkeypatch.setattr(rpi, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": str(src.relative_to(tmp_path)),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp1",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{"metric_type": "periodic_in_axenic_LD", "value_kind": "boolean",
                     "value_col": "flag", "true_tokens": ["Y"]}],
    }

    result = rpi.resolve_derived_metrics_entry("Biller 2018", "dm_s4a", entry, force=True)
    assert result is not None
    assert not result.get("skipped")

    resolved = tmp_path / "s4a_resolved.csv"
    assert resolved.exists()

    out = pd.read_csv(resolved)
    assert "resolved_locus_tag" in out.columns
    assert "resolution_method" in out.columns
    assert out.loc[out["NCBI ID_2"] == "PMN2A_RS00015", "resolved_locus_tag"].iloc[0] == "PMN2A_RS00015"
    assert pd.isna(out.loc[out["NCBI ID_2"] == "UNKNOWN_ID", "resolved_locus_tag"].iloc[0])
```

- [ ] **Step 2: Run test — expect ImportError / AttributeError**

Run: `uv run pytest tests/test_resolve_paper_ids_derived_metrics.py -v`
Expected: FAIL — `resolve_derived_metrics_entry` does not exist.

- [ ] **Step 3: Add the resolver**

In `multiomics_kg/download/resolve_paper_ids.py`, add immediately after `resolve_gene_clusters_table` (around line 243):

```python
def resolve_derived_metrics_entry(
    pub_name: str,
    table_key: str,
    table_config: dict,
    force: bool = False,
) -> dict | None:
    """Resolve gene IDs in one derived_metrics_table supplementary entry.

    Same shape as resolve_table but:
      - name_col is at the entry level, not under statistical_analyses
      - organism is at the entry level (or derived via get_organism_for_entry)
    """
    if table_config.get("type") != "derived_metrics_table":
        return None

    filename = table_config.get("filename")
    if not filename:
        return None

    src = PROJECT_ROOT / filename
    if not src.exists():
        print(f"  [warn] derived_metrics CSV not found: {src}", file=sys.stderr)
        return None

    name_col = table_config.get("name_col")
    if not name_col:
        return {"skipped": True, "reason": "no name_col defined"}

    if name_col == "locus_tag":
        return {"skipped": True, "reason": "name_col is already locus_tag"}

    organism = table_config.get("organism", "")
    if not organism:
        print(
            f"  [warn] No organism for derived_metrics {pub_name}/{table_key} — skipping",
            file=sys.stderr,
        )
        return None

    resolved_path = get_resolved_path(src)
    sep = table_config.get("sep", ",")
    skip_rows = int(table_config.get("skip_rows", 0) or 0)

    if not force and resolved_path.exists():
        if resolved_path.stat().st_mtime >= src.stat().st_mtime:
            return {
                "skipped": True, "reason": "up to date",
                "resolved_path": str(resolved_path.relative_to(PROJECT_ROOT)),
            }

    genome_dir = get_genome_dir(organism, str(PROJECT_ROOT))
    if not genome_dir:
        print(
            f"  [warn] No genome dir for organism '{organism}' in {pub_name}/{table_key} — skipping",
            file=sys.stderr,
        )
        return None

    mapping_data = load_mapping_v2(genome_dir)
    if mapping_data is None:
        from multiomics_kg.utils.gene_id_utils import build_id_lookup
        lookup, locus_tags, _supp_keys = build_id_lookup(genome_dir)
        if lookup is None:
            print(
                f"  [warn] No gene annotation data for '{organism}' in {pub_name}/{table_key} — skipping",
                file=sys.stderr,
            )
            return None
        mapping_data = MappingData(
            specific_lookup=lookup, multi_lookup={}, conflicts={},
            locus_tags=locus_tags if locus_tags else set(), version=0,
        )

    id_columns: list[dict] = table_config.get("id_columns") or []

    try:
        df = pd.read_csv(src, sep=sep, skiprows=skip_rows if skip_rows else None)
    except Exception as e:
        print(f"  [error] Could not read {src}: {e}", file=sys.stderr)
        return None

    if name_col not in df.columns:
        print(
            f"  [warn] name_col '{name_col}' not in {src.name} — skipping {pub_name}/{table_key}",
            file=sys.stderr,
        )
        return None

    locus_tag_col: list[str | None] = []
    method_col: list[str] = []
    method_counts: Counter = Counter()
    unresolved_rows: list[dict] = []
    total = 0
    n_resolved = 0

    for row_idx, row in df.iterrows():
        name_val = row.get(name_col)
        if pd.isna(name_val) or str(name_val).strip() == "":
            locus_tag_col.append(None)
            method_col.append("empty")
            method_counts["empty"] += 1
            continue

        total += 1
        lt, method = resolve_row(row.to_dict(), name_col, id_columns, mapping_data)
        locus_tag_col.append(lt)
        method_col.append(method)
        method_counts[method] += 1

        if lt:
            n_resolved += 1
        else:
            unresolved_rows.append({
                "row": row_idx,
                "raw_id": str(name_val).strip(),
                "method": method,
                "id_col_vals": {
                    c.get("column", ""): str(row.get(c.get("column", ""), ""))
                    for c in id_columns
                    if c.get("column", "") in df.columns
                },
            })

    RESOLVED_COL = "resolved_locus_tag"
    df_out = df.copy()
    if RESOLVED_COL in df_out.columns:
        print(
            f"  [error] Column '{RESOLVED_COL}' already exists in {src.name} — "
            f"cannot write resolved output for {pub_name}/{table_key}",
            file=sys.stderr,
        )
        return None
    df_out[RESOLVED_COL] = locus_tag_col
    df_out["resolution_method"] = method_col
    try:
        df_out.to_csv(resolved_path, index=False)
    except Exception as e:
        print(f"  [error] Could not write {resolved_path}: {e}", file=sys.stderr)
        return None

    _write_report(
        pub_name, table_key, name_col, id_columns, src, resolved_path,
        total, n_resolved, method_counts, unresolved_rows,
    )

    return {
        "skipped": False,
        "pub_name": pub_name,
        "table_key": table_key,
        "organism": organism,
        "name_col": name_col,
        "src": str(src.relative_to(PROJECT_ROOT)),
        "resolved_path": str(resolved_path.relative_to(PROJECT_ROOT)),
        "total": total,
        "resolved": n_resolved,
        "method_counts": dict(method_counts),
        "unresolved_rows": unresolved_rows,
    }
```

- [ ] **Step 4: Dispatch the new type in `main()`**

In `main()` of `resolve_paper_ids.py` (around line 594), find the `if table_type == "gene_clusters":` branch. Extend the if/else so the new type dispatches correctly:

```python
            if table_type == "gene_clusters":
                result = resolve_gene_clusters_table(
                    pub_name, table_key, table_config, force=args.force
                )
            elif table_type == "derived_metrics_table":
                result = resolve_derived_metrics_entry(
                    pub_name, table_key, table_config, force=args.force
                )
            else:
                result = resolve_table(
                    pub_name, table_key, table_config, force=args.force, config=config
                )
```

- [ ] **Step 5: Run tests**

Run: `uv run pytest tests/test_resolve_paper_ids_derived_metrics.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/resolve_paper_ids.py tests/test_resolve_paper_ids_derived_metrics.py
git commit -m "resolve_paper_ids: resolve gene IDs in derived_metrics_table entries"
```

---

### Task 11: Author Biller 2018 `derived_metrics_table` entries (S4A, S4B, S5)

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`

- [ ] **Step 1: Append 3 supplementary-materials entries**

Open `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`. Under the existing `publication.supplementary_materials:` block, **after** the `supp_table_s6b:` entry (i.e., at the end of the file, preserving existing 2-space indent), append:

```yaml
    # S4A is split into two derived_metrics_table entries — one per parent
    # Experiment — because the spec pins one DerivedMetric per (Experiment ×
    # metric_type) and S4A reports periodicity signals from BOTH the axenic
    # NATL2A experiment (cols "in axenic, L:D" + "in axenic, extended darkness")
    # and the coculture NATL2A experiment (cols "in co-cultured, L:D" + "in
    # co-cultured, extended darkness"). Both entries share the same source CSV
    # file and id_columns but select disjoint value_col sets. The deleted
    # gene_clusters entry preserved this dual attribution via a 2-element
    # `experiments:` list; derived_metrics_table has no such list — the split
    # restores the attribution.
    natl2a_periodicity_axenic_derived_metrics:
      type: derived_metrics_table
      filename: data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4.csv
      organism: Prochlorococcus NATL2A
      experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
      name_col: NCBI ID_2
      id_columns:
      - column: NCBI ID_2
        id_type: locus_tag_ncbi
      - column: NCBI ID_3
        id_type: old_locus_tag
      - column: NCBI ID
        id_type: locus_tag_cyanorak
      metrics:
      - metric_type: periodic_in_axenic_LD
        value_kind: boolean
        value_col: "Periodic in axenic, L:D cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in NATL2A axenic L:D (Table S4A)"
      - metric_type: periodic_in_axenic_extended_darkness
        value_kind: boolean
        value_col: "Periodic in axenic, extended darkness cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in NATL2A axenic extended darkness (Table S4A)"

    natl2a_periodicity_coculture_derived_metrics:
      type: derived_metrics_table
      filename: data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4.csv
      organism: Prochlorococcus NATL2A
      experiment: darkness_extended_darkness_natl2a_rnaseq_coculture
      name_col: NCBI ID_2
      id_columns:
      - column: NCBI ID_2
        id_type: locus_tag_ncbi
      - column: NCBI ID_3
        id_type: old_locus_tag
      - column: NCBI ID
        id_type: locus_tag_cyanorak
      metrics:
      - metric_type: periodic_in_coculture_LD
        value_kind: boolean
        value_col: "Periodic in co-cultured, L:D cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in NATL2A coculture L:D (Table S4A)"
      - metric_type: periodic_in_coculture_extended_darkness
        value_kind: boolean
        value_col: "Periodic in co-cultured, extended darkness cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in NATL2A coculture extended darkness (Table S4A)"

    mit1002_periodicity_derived_metrics:
      type: derived_metrics_table
      filename: data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4.csv
      organism: Alteromonas macleodii MIT1002
      experiment: darkness_extended_darkness_mit1002_rnaseq
      name_col: Locus ID
      id_columns:
      # S4B's `Locus ID` column contains 4-digit MIT1002_0002 values. Canonical
      # NCBI MIT1002 locus tags are 5-digit (MIT1002_00013); the 4-digit form
      # matches the `genbank_id` column of MIT1002_systematicnames_conversiontable.csv,
      # which the existing id_translation_mit1002_rast_names entry declares as
      # `alternative_locus_tag`. Use the same id_type here so build_gene_id_mapping
      # bridges 4-digit → canonical via the existing transitive closure.
      - column: Locus ID
        id_type: alternative_locus_tag
      metrics:
      - metric_type: periodic_in_coculture_LD
        value_kind: boolean
        value_col: "Periodic in co-cultured, L:D cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in MIT1002 coculture L:D (Table S4B)"
      - metric_type: periodic_in_coculture_extended_darkness
        value_kind: boolean
        value_col: "Periodic in co-cultured, extended darkness cultures"
        true_tokens: ["Y"]
        false_tokens: []
        skip_tokens: ["NA", "N/A"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in MIT1002 coculture extended darkness (Table S4B)"

    # S5's darkness_cluster is a SINGLE composite column mixing axenic and
    # coculture signals into one categorical label per gene; it cannot be split
    # column-wise like S4A. The deleted gene_clusters entry linked S5 to both
    # axenic and coculture NATL2A experiments via its `experiments:` list;
    # derived_metrics_table has a single `experiment:` field, so we attribute
    # to the axenic experiment as the canonical parent. The composite label
    # itself encodes the cross-condition information, so the cross-experiment
    # attribution is not lost — it lives in the category names. All 4 id_columns
    # from the deleted entry are preserved, including `Gene Name` → gene_name.
    natl2a_darkness_survival_derived_metrics:
      type: derived_metrics_table
      filename: data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5.csv
      organism: Prochlorococcus NATL2A
      experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
      name_col: NCBI ID_2
      id_columns:
      - column: NCBI ID_2
        id_type: locus_tag_ncbi
      - column: NCBI ID_3
        id_type: old_locus_tag
      - column: NCBI ID
        id_type: locus_tag_cyanorak
      - column: Gene Name
        id_type: gene_name
      metrics:
      - metric_type: darkness_survival_class
        value_kind: categorical
        value_col: darkness_cluster
        # allowed_categories are the 3 distinct composite strings actually
        # present in the CSV's darkness_cluster column (verified: 100 / 90 /
        # 79 rows). Under option C, this list lives on the paperconfig, not
        # the central vocab — categorical class vocabularies are per-paper.
        allowed_categories:
        - darkness_axenic+darkness_coculture
        - darkness_coculture+unique_coculture
        - darkness_axenic+unique_axenic
        field_description: "Transcript presence at 72-144h extended darkness, axenic vs coculture (Biller 2018 Table S5)"
```

**Notes:**
- S4B's `Locus ID` column contains values with a trailing space (`"MIT1002_0002 "`). The extractor / resolver both call `pd.read_csv(..., dtype=str)` and strip via `_safe_str()` / `str(name_val).strip()`, so the trailing space is absorbed. No CSV edit required.
- S4A is authored as **two** `derived_metrics_table` entries (axenic + coculture) that share the same source CSV but select disjoint `value_col` sets — restoring the dual-experiment attribution the deleted `natl2a_periodicity` cluster entry carried via its 2-element `experiments:` list.
- S5 stays as a single entry because its `darkness_cluster` is one composite column mixing axenic + coculture signals; it's linked to the axenic NATL2A experiment as the canonical parent. The cross-condition information is preserved inside the category names themselves.
- All id_columns from the deleted cluster entries are carried over verbatim — including `Gene Name` → `gene_name` on S5 (which was absent from an earlier draft of this plan).
- No `rankable` / `has_p_value` on bool/categorical entries. Option-C paperconfig field rules forbid these on bool/categorical (see Task 8 dispatch table + `_validate_derived_metrics_entry`). The Plan-2 adapter sets both to `"false"` on the emitted `DerivedMetric` node when `value_kind != "numeric"` — they are definitionally false, not a paper-specific claim.
- No `_modified.csv` is created; the source CSVs remain untouched (per the "never edit a paper's source CSV in place" convention).

- [ ] **Step 2: Validate the retrofitted paperconfig**

Run: `uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml" 2>&1 | tail -20`
Expected: `VALIDATION PASSED`. If warnings appear about unknown `id_type` on `NCBI ID` → `locus_tag_cyanorak`, or about the composite-label semantics, read them but do not treat as failures.

- [ ] **Step 3: Run the full regression**

Run: `uv run python scripts/validate_paperconfig.py --all 2>&1 | tail -5`
Expected: `Results: N/N passed, 0/N failed`.

- [ ] **Step 4: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
git commit -m "biller 2018: add derived_metrics_table entries for S4A/S4B/S5"
```

---

### Task 12: Preprocessing regeneration — refresh NATL2A + MIT1002

**Files:**
- Modify: `cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json` (regenerated)
- Modify: `cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json` (regenerated)
- Create: `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv`
- Create: `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv`
- Create: `data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv`
- Create: corresponding `_resolved_report.txt` alongside each

- [ ] **Step 1: Regenerate gene_id_mapping.json for both strains**

Run: `bash scripts/prepare_data.sh --steps 3 --strains NATL2A MIT1002 --force`
Expected: two strains processed; logs at `logs/prepare_data_step3.log`. Each strain prints the number of collected rows from the new `derived_metrics_table` entries.

- [ ] **Step 2: Spot-check mapping counts**

Run: `uv run python -c "
import json
for strain in ['NATL2A', 'MIT1002']:
    for org in ['Prochlorococcus', 'Alteromonas']:
        import pathlib
        p = pathlib.Path(f'cache/data/{org}/genomes/{strain}/gene_id_mapping.json')
        if p.exists():
            d = json.loads(p.read_text())
            print(f'{strain}: n_genes={d[\"stats\"][\"n_genes\"]} n_specific={d[\"stats\"][\"n_specific\"]} n_multi={d[\"stats\"][\"n_multi\"]} n_conflicts={d[\"stats\"][\"n_conflicts\"]}')
"`
Expected: NATL2A and MIT1002 stats print; `n_conflicts` unchanged from Step-0 baseline. If `n_conflicts` increases, the new entries' `id_columns` contain a Tier-1 ID that disagrees with an existing lock — re-check the table to fix (usually means the wrong `id_type` was declared).

- [ ] **Step 3: Pre-resolve the 3 new source CSVs via 4 paperconfig entries**

Run: `uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Biller 2018" --force 2>&1 | tail -40`
Expected: **4 resolver invocations** (one per new paperconfig entry — `natl2a_periodicity_axenic_derived_metrics`, `natl2a_periodicity_coculture_derived_metrics`, `mit1002_periodicity_derived_metrics`, `natl2a_darkness_survival_derived_metrics`) plus re-resolution of the 2 existing DE csv entries (unchanged). Because the two S4A entries share one source CSV, they produce **3 unique `_resolved.csv` files** — the second S4A invocation overwrites the first with identical content (all rows resolve against the same `name_col: NCBI ID_2` with the same `id_columns`; value_col selection differs, but resolve_paper_ids only cares about gene-ID columns, not per-metric columns). Reported resolution rates: S4A ≥98%, S4B ≥98%, S5 ≥98%. Typical unresolved rows: 0–1 per table (row-1 header artifact if any; any NCBI IDs withdrawn between 2018 and today).

- [ ] **Step 4: Run `/check-gene-ids` (skill) for the 3 new tables**

Run: `uv run python -c "
import pandas as pd
for fn, name_col in [
    ('data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv', 'NCBI ID_2'),
    ('data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv', 'Locus ID'),
    ('data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv',  'NCBI ID_2'),
]:
    df = pd.read_csv(fn)
    total = len(df)
    resolved = df['resolved_locus_tag'].notna().sum()
    print(f'{fn.split(chr(47))[-1]}: resolved {resolved}/{total} ({100*resolved/max(total,1):.1f}%)')
"`
Expected: each table ≥98% resolved.

- [ ] **Step 5: Commit regenerated mapping + resolved CSVs**

```bash
git add cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json \
        cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping_report.json \
        cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json \
        cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping_report.json \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv" \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved_report.txt" \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv" \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved_report.txt" \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv" \
        "data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved_report.txt"
git commit -m "biller 2018: regenerate gene_id_mapping + pre-resolved DM CSVs (NATL2A, MIT1002)"
```

---

### Task 13: Smoke build + definition-of-done check

**Files:**
- No code modifications — validation only.

- [ ] **Step 1: Run the full unit-test suite**

Run: `uv run pytest -m "not slow and not kg" -v 2>&1 | tail -20`
Expected: all tests pass. Count should be ≥ the pre-Plan-1 count (existing tests + Plan 1's additions).

- [ ] **Step 2: Run test-mode KG build**

Run: `uv run python create_knowledge_graph.py --test 2>&1 | tee /tmp/plan1_build.log | tail -30`
Expected: pipeline completes. No `DerivedMetric` CSV appears under `biocypher-log/example_knowledge_graph/` (no adapter emits yet). The `MultiOMICSAdapter` and `MultiClusterAdapter` should process Biller 2018 identically to the Step-0 baseline — the new `derived_metrics_table` entries are skipped by both because they filter on `type == "csv"` / `type == "gene_clusters"` respectively.

- [ ] **Step 3: Verify CSV output unchanged vs Step-0 baseline**

Run: `uv run python -c "
import pathlib
out = pathlib.Path('biocypher-log/example_knowledge_graph')
by_label = {}
for p in out.glob('**/*.csv'):
    label = p.name.split('-')[0]
    by_label.setdefault(label, 0)
    by_label[label] += sum(1 for _ in p.open()) - 1
for k, v in sorted(by_label.items()):
    print(f'{k}: {v}')
"`
Expected: output matches pre-Plan-1 structure — no new `DerivedMetric*` CSV files; `Changes_expression_of` row count is the same number emitted by Step 0's test-mode build (full-mode rebuild deferred to Plan 2).

- [ ] **Step 4: Run the real-paper validator regression**

Run: `uv run python scripts/validate_paperconfig.py --all 2>&1 | tail -5`
Expected: `Results: N/N passed, 0/N failed` (same N as Step-0 baseline).

- [ ] **Step 5: Commit nothing (this task is pure verification)**

If Steps 1–4 pass, Plan 1 is complete per its Definition of Done. Open Plan 2 next.

---

## Definition of Done (end-to-end checklist for Plan 1)

- [ ] `multiomics_kg/vocab/non_de_evidence.py` exports `COMPARTMENTS`, `EXTENDED_OMICS_TYPES`, `VALUE_KINDS`, `KNOWN_METRIC_TYPES: dict[str, str]`, token + bucket constants. No `MetricSpec` dataclass.
- [ ] `tests/test_vocab_non_de_evidence.py` green.
- [ ] `config/schema_config.yaml` declares `DerivedMetric` node, 3 binding edges, 3 measurement edges, `Experiment.compartment`.
- [ ] `test build`: `uv run python create_knowledge_graph.py --test` succeeds; no `DerivedMetric` CSVs emitted (no adapter yet).
- [ ] `paperconfig_utils.iter_derived_metrics_tables` exported + tested.
- [ ] `validate_paperconfig_content(config, path) -> (errors, warnings)` is a pure function importable from `scripts/validate_paperconfig`.
- [ ] Validator enforces: `Experiment.compartment` vocab; `PAIRED_RNASEQ_PROTEOME` omics_type accepted; DE fields (`control_condition`, `test_type`) optional for experiments whose only supplementary entries are `derived_metrics_table` / `gene_clusters`.
- [ ] Validator dispatches `derived_metrics_table` entries per `value_kind`: all top-level required fields present (filename, organism, experiment, name_col, metrics); `value_kind` ∈ `VALUE_KINDS`; `KNOWN_METRIC_TYPES` cross-check (hard error on value_kind mismatch; warning on novel metric_type); numeric path requires `value_col` + `rankable` + `has_p_value` and gates p-value fields on `has_p_value=="true"`; boolean path requires `true_tokens` (non-empty) and CSV-dry-runs for unclassified tokens (hard error); categorical path requires `allowed_categories` (non-empty list) inline and CSV-dry-runs for out-of-set values (warning). Forbidden-field gating enforced per value_kind.
- [ ] `build_gene_id_mapping.extract_rows_from_derived_metrics_table` ingests `id_columns` from `derived_metrics_table` entries (tested).
- [ ] `resolve_paper_ids.resolve_derived_metrics_entry` writes `<stem>_resolved.csv` next to each of Biller 2018's 3 new CSVs (tested).
- [ ] Biller 2018 paperconfig validates with `uv run python scripts/validate_paperconfig.py` against its new entries.
- [ ] `prepare_data.sh --steps 3 --strains NATL2A MIT1002 --force` regenerates `gene_id_mapping.json` for both strains; no new conflicts vs the Step-0 baseline.
- [ ] `resolve_paper_ids --papers "Biller 2018" --force` produces 3 new `_resolved.csv` files (from 4 paperconfig invocations — S4A-axenic + S4A-coculture share a source CSV) with ≥98% resolution each.
- [ ] `uv run python scripts/validate_paperconfig.py --all` — all paperconfigs pass.
- [ ] `uv run pytest -m "not slow and not kg" -v` — all tests pass.
- [ ] `uv run python create_knowledge_graph.py --test` succeeds; CSV output shape matches Step-0 baseline.

---

## Risks and open questions

- **S5's `darkness_cluster` values vs parent-spec example categories.** The parent spec lists `{present_in_axenic, present_in_coculture, unique_to_axenic, unique_to_coculture}` as illustrative ("e.g., …"). Real CSV values are composite labels (`darkness_axenic+darkness_coculture` etc.). Under option C, `allowed_categories` lives on the paperconfig entry, so the vocabulary is per-paper by construction — no central-vocab collision risk. A future paper that wants to re-encode "darkness survival" with cleaner labels can either reuse `metric_type: darkness_survival_class` with different `allowed_categories`, or introduce a new metric_type name (e.g., `darkness_survival_binary`). Both are fine; the registry only pins value_kind.
- **Option C vs option A (central-registry) design trade-off.** Central registry shrinks to `metric_type → value_kind` only. Benefit: no vocabulary bloat as new papers ship. Cost: authors must declare `rankable` / `has_p_value` / `allowed_categories` on every paperconfig entry — the validator catches forgetting, but there's no per-metric-default inheritance. Acceptable: each per-metric declaration is a specific claim about what *this* paper did, which is the right place to localize it.
- **S4A is split across two `derived_metrics_table` entries** (one per parent Experiment) so the axenic-vs-coculture attribution the deleted cluster entry preserved via its 2-element `experiments:` list is not lost. Plan 2's adapter processes the two entries independently; each emits 2 `DerivedMetric` nodes (axenic LD + axenic ED, or coculture LD + coculture ED) bound to the appropriate Experiment. Shared source CSV, disjoint `value_col` columns.
- **S5 stays as a single categorical entry** bound to the axenic NATL2A experiment. Its `darkness_cluster` column mixes signals from both axenic and coculture conditions into one composite category string per gene — the column cannot be split, so the axenic-vs-coculture distinction is encoded in the category names themselves (`darkness_axenic+darkness_coculture` etc.).
- **MIT1002 `Locus ID` trailing-whitespace handling.** Relies on `pd.read_csv(dtype=str)` + `_safe_str` / `str.strip()` absorption — which both the extractor (Task 9) and the resolver (Task 10) already do. Verified by inspection.
- **MIT1002 4-digit `Locus ID` bridging.** S4B's `Locus ID` column uses `MIT1002_0002` 4-digit — the same form as `genbank_id` in `MIT1002_systematicnames_conversiontable.csv`. The existing `id_translation_mit1002_rast_names` entry declares that column as `alternative_locus_tag`, so Task 11 uses the same id_type here to let `build_gene_id_mapping`'s transitive closure bridge 4-digit → canonical 5-digit (`MIT1002_NNNNN`) via rast_fig_id / diamond. If Task 12 Step 3 reports <98% resolution for S4B, re-check the chain by inspecting the resolved report — the fix is typically to declare an extra id_column bridge, not to change the metric vocabulary.
- **Validator CSV dry-run cost.** The boolean dry-run reads the full `value_col`. Biller 2018's largest CSV is 2,198 rows (S4A) — sub-second. The parent spec flags this as a future concern if CSVs exceed 100K rows; not a concern at this slice size.
- **`validate_paperconfig_content` refactor scope.** Task 6 is a large mechanical move. Risk: a missed `_print_results(errors, warnings)` call or a subtle variable-scope change breaks existing real-paper regressions. Mitigation: Task 6 Step 6 runs the full `--all` regression immediately after the refactor, before any new behavior lands.
- **No KG-validity coverage in this plan.** KG validity tests for DerivedMetric land in Plan 3, by which point the Plan 2 adapter has produced real edges for those tests to inspect.

## Dependencies on other plans

- **Upstream:** Step 0 (commit 25f0cd8) already dropped the 3 `gene_clusters` entries from Biller 2018 and captured the baseline. This plan builds directly on that state.
- **Downstream (Plan 2):** will create `multiomics_kg/adapters/observations_adapter.py`, consume the Biller 2018 `derived_metrics_table` entries authored here, wire into `create_knowledge_graph.py`, and commit the synthetic numeric-DM paperconfig fixture under `tests/fixtures/non_de/`. Plan 2 does NOT need to modify the schema or vocab — both land here.
- **Downstream (Plan 3):** post-import Cypher, KG validity assertions, `snapshot_data.json` regeneration, `docs/kg-changes/non-de-evidence-extension.md`, skill docs. Plan 3 relies on Plan 2's fixture being ingested.
