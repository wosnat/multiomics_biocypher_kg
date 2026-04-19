# Non-DE Evidence — Biller 2018 Slice

**Date:** 2026-04-19
**Status:** Approved scope — awaiting implementation plans
**Parent spec:** `docs/superpowers/specs/2026-04-17-non-de-evidence-extension-design.md`

## Context

The parent spec covers the full non-DE-evidence surface: `AbundanceAnalysis` + `DerivedMetric` nodes with 10 new edge types, post-import rollups, KG validity tests. This doc narrows the **first implementation slice** to what a single real paper — **Biller 2018** — drives: full `DerivedMetric` infrastructure (numeric, boolean, categorical), with abundance deferred.

Slicing rationale: a real paper as driver surfaces gene-ID resolution edge cases, retrofit regression, and denormalization invariants against actual data rather than synthetic fixtures. Biller 2018 exercises boolean + categorical DerivedMetric through its Tables S4A / S4B / S5; numeric DM infrastructure lands in full and is covered by a synthetic fixture (~100 rows) sized to exercise bucket thresholds meaningfully.

Abundance is deferred — no paper in this slice requires it, and the future `biller 2022` / `Oleza 2015` slices are its natural drivers.

## In scope

### Schema (`config/schema_config.yaml`)

- `DerivedMetric` node with full properties (numeric metadata included even though the real paper doesn't use it, so the zinser 2009 follow-up is a paperconfig-authoring task only)
- Binding edges: `publication_has_derived_metric`, `experiment_has_derived_metric`, `derived_metric_belongs_to_organism`
- Measurement edges: `derived_metric_quantifies_gene` (numeric), `derived_metric_flags_gene` (boolean), `derived_metric_classifies_gene` (categorical)
- `Experiment.compartment` (adapter-emitted, default `"whole_cell"`)
- All post-import-computed properties declared in schema (filled in Plan 3)

### Vocabulary module

Create `multiomics_kg/vocab/non_de_evidence.py` exporting:

- `COMPARTMENTS` (5 values per parent spec §3)
- `METRIC_TYPES` (full set per parent spec §2: 6 numeric + 5 boolean + 1 categorical), as a `dict[str, MetricSpec]` dataclass
- `EXTENDED_OMICS_TYPES` (existing 5 + `PAIRED_RNASEQ_PROTEOME`)
- `DEFAULT_SKIP_TOKENS`, `VALID_BLANK_POLICIES`, `BUCKET_THRESHOLD_*` constants

Not exported in this slice: `VALUE_TYPES` (abundance-only).

### Paperconfig layer

- New supplementary_materials type `derived_metrics_table` recognized by validator, `build_gene_id_mapping`, `resolve_paper_ids`, and adapter
- Optional `compartment` on the `experiments:` block (vocab-checked against `COMPARTMENTS`)
- DE-only field relaxation extended: experiments whose only supplementary entries are `derived_metrics_table` or `gene_clusters` may omit `control_condition` and `test_type`
- Validator dispatch:
  - `numeric` branch: vocab-match on `value_kind`, `unit`, `has_p_value`; `p_value_col` / `adjusted_p_value_col` / `p_value_threshold` gated by vocab `has_p_value`
  - `boolean` branch: `true_tokens` required + CSV dry-run (every non-blank cell must be in `true_tokens ∪ false_tokens ∪ skip_tokens`; unexpected tokens are hard errors)
  - `categorical` branch: CSV dry-run warning on out-of-`allowed_categories` values (adapter hard-errors at ingest)
- Extracted testable `validate_paperconfig_content(config, path) → (errors, warnings)` function lifted from the existing CLI-only `validate()` entry

Deferred (no abundance entries in this slice):
- Cross-entry "no redundant per-timepoint abundance" invariant

### Preprocessing pipeline

- `multiomics_kg/download/build_gene_id_mapping.py` — one new extractor + dispatcher case for `derived_metrics_table` (shape identical to `extract_rows_from_csv_table` minus the `statistical_analyses` detour)
- `multiomics_kg/download/resolve_paper_ids.py` — `resolve_derived_metrics_entry` + main-loop integration
- `paperconfig_utils.py` — `iter_derived_metrics_tables` iterator

### Adapter

Create `multiomics_kg/adapters/observations_adapter.py` as a sibling to `cluster_adapter.py` and `omics_adapter.py`:

- `ObservationsAdapter` — consumes one paperconfig; for each `derived_metrics_table` entry, emits one `DerivedMetric` node per `metric_type` + 3 binding edges + one of three measurement edges per metric based on vocab `value_kind`
- `MultiObservationsAdapter` — wraps `paperconfig_files.txt`, builds `_organism_lookup` from `cyanobacteria_genomes.csv`, delegates per-paperconfig
- **ID conventions match cluster_adapter** (verified at `multiomics_kg/adapters/cluster_adapter.py`):
  - Publication = `doi:{doi}` (with `doi:` prefix)
  - Experiment = `{doi}_{exp_key}` (no prefix on doi)
  - OrganismTaxon target = `_organism_lookup[preferred_name]` → `insdc.gcf:{accession}`
  - DerivedMetric node ID = `derived_metric:{doi_short}:{entry_key}:{metric_type}`
- Boolean token parser with hard-error on unexpected tokens
- Categorical hard-error on out-of-`allowed_categories` values
- Denormalized fields on DerivedMetric computed from parent Experiment, never re-read from paperconfig — drift impossible
- String sanitization via local `_clean_str()` helper per CLAUDE.md convention (`"'" → "^"`, `"|" → ""`)
- Wired into `create_knowledge_graph.py` after `MultiClusterAdapter`

### Post-import (`scripts/post-import.sh` + `scripts/post-import.cypher`, byte-identical per CLAUDE.md)

- **Rank / percentile / bucket pass** for `derived_metric_quantifies_gene` edges grouped by `(DerivedMetric)`, only when parent `rankable="true"`. Thresholds pinned: `top_decile` (≥90), `top_quartile` (75–<90), `mid` (25–<75), `low` (<25)
- **Significance derivation** on numeric edges only when parent `has_p_value="true"` AND `p_value_threshold IS NOT NULL` AND the edge's `adjusted_p_value IS NOT NULL`: `significant = "true"` iff `adjusted_p_value < p_value_threshold`
- **DerivedMetric rollups**: `total_gene_count`, `growth_phases` (from parent Experiment)
- **Experiment rollups**: `reports_fold_change`, `reports_derived_metric_types`, `derived_metric_count`, `derived_metric_value_kinds`, `derived_metric_gene_count`
- **Publication rollups** (same shape, DerivedMetric-scoped)
- **OrganismTaxon rollups** (same shape, DerivedMetric-scoped)
- **Gene routing counts**: `numeric_metric_count`, `classifier_flag_count`, `classifier_label_count`, `numeric_metric_types_observed`, `classifier_flag_types_observed`, `classifier_label_types_observed`, `compartments_observed` (in this slice `compartments_observed` is populated only from DerivedMetric children — future abundance slice extends to include AbundanceAnalysis)
- **New scalar indexes**: `derived_metric_metric_type_idx`, `derived_metric_value_kind_idx`, `derived_metric_compartment_idx`, `derived_metric_omics_type_idx`, `derived_metric_treatment_type_idx`, `derived_metric_organism_idx`, `derived_metric_experiment_idx`, `experiment_compartment_idx`
- **New full-text index**: `derivedMetricFullText` on `DerivedMetric(name, field_description)`
- **Empty-state defaults** on all new int and `str[]` rollups (parent spec §Success criteria)

Deferred post-import passes: abundance rank/percentile/bucket, `Experiment.abundance_time_point_*`, abundance binding rollups on Publication/OrganismTaxon.

### Tests

- Adapter unit tests per value_kind path (numeric / boolean / categorical emission), token parser (unknown tokens → `ValueError`), denormalization correctness
- Validator unit tests per `derived_metrics_table` shape (field requirements, per-value_kind gating, CSV dry-runs)
- `build_gene_id_mapping` tests for the new extractor
- `resolve_paper_ids` tests for the new resolver
- KG validity tests — the parent-spec §Success-criteria assertions filtered to DerivedMetric scope:
  - Edge target is `:Gene`
  - One edge type per DerivedMetric matching `value_kind`
  - Denormalized fields on DerivedMetric equal parent Experiment's
  - `value_text` ∈ `allowed_categories` on every classifies edge
  - `value_flag` ∈ `{"true","false"}` on every flags edge
  - Bucket ↔ percentile consistency on quantifies edges
  - `significant` gating consistency
  - Rollup values match direct query-time counts (Experiment / Publication / OrganismTaxon / Gene)
  - Empty-state defaults (int → 0, str[] → [])
- **Synthetic numeric-DM paperconfig fixture** (~100 rows, committed under `tests/fixtures/non_de/`) exercising numeric DM rank/percentile/bucket/significance against live Neo4j

### Docs and skills

- `docs/kg-changes/non-de-evidence-extension.md` — scoped to DerivedMetric shape; explicitly flags AbundanceAnalysis as a follow-up
- `.claude/skills/paperconfig/SKILL.md` — new `derived_metrics_table` section + `compartment` field + `PAIRED_RNASEQ_PROTEOME` omics_type
- `.claude/skills/cypher-queries/SKILL.md` — query templates for flags / classifies / quantifies edges
- `.claude/skills/omics-edge-snapshot/SKILL.md` — extend to report DerivedMetric edge counts alongside existing `Changes_expression_of`
- `CLAUDE.md` — key-facts update (DerivedMetric node ID format, the 3 new edge types, Biller 2018 retrofit notes)

## Out of scope — deferred to future slices

- `AbundanceAnalysis` node + 3 binding edges + `abundance_analysis_measures_gene`
- `abundance_table` + `timeseries_abundance_table` paperconfig types
- `VALUE_TYPES` vocab
- Abundance post-import passes (`rank_in_sample`, `abundance_percentile`, `abundance_bucket`)
- `Experiment.abundance_time_point_*` parallel arrays
- `Gene.abundance_analysis_count`
- Real integration of `PAIRED_RNASEQ_PROTEOME` papers (Waldbauer 2012) — infra lands, no paper tests it
- Cross-entry "no redundant per-timepoint abundance" validator invariant (relevant only when abundance entries exist)

Future slice drivers:

- biller 2022 vesicles (S2 + S3) → `abundance_table` + `compartment: vesicle`
- Oleza 2015 / Oleza 2017 → `abundance_table` + `compartment: exoproteome`
- Waldbauer 2012 → `PAIRED_RNASEQ_PROTEOME` + real numeric DM
- zinser 2009 → real numeric DM retrofit (pure paperconfig-authoring task once this slice ships)

## Execution shape

### Step 0 — Retrofit removal + baseline capture

Single small PR, not a full plan.

**Scope**
- Edit `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`: remove the 3 `gene_clusters` entries (`natl2a_periodicity`, `mit1002_periodicity`, `natl2a_darkness_survival`)
- Delete the stale `<stem>_resolved.csv` files generated for the 3 dropped CSVs (Tables S4A / S4B / S5) — prevents ghost preprocessing state
- `bash scripts/prepare_data.sh --steps 3 4 --strains NATL2A MIT1002 --force` — refresh preprocessing after the removal
- `docker compose down && docker compose up -d` — rebuild KG end-to-end
- `pytest -m "not slow and not kg"` + `pytest -m kg` — fix any snapshot mismatch by regenerating `tests/kg_validity/snapshot_data.json`

**Captured artifacts** (paths finalized at implementation time; proposed under `docs/kg-changes/biller-2018-retrofit-baseline/`):
- `/omics-edge-snapshot` JSON output
- Canned Biller 2018 Cypher dump — Publication + Experiments + remaining ClusteringAnalysis (if any) + Changes_expression_of counts + Gene_in_gene_cluster counts
- Regenerated `snapshot_data.json` committed alongside

**Definition of done**
- `pytest -m "not slow and not kg"` green
- `pytest -m kg` green
- `/omics-edge-snapshot` shows expected drop in `Gene_in_gene_cluster` edges + zero unrelated regressions
- `Changes_expression_of` count unchanged
- Baseline artifacts committed

**Transient state flag**: between step 0 and Plan 2, the 3 retrofitted tables' evidence is absent from the KG. MCP queries against those cluster IDs will break in the interim. `docs/kg-changes/non-de-evidence-extension.md` (Plan 3) documents this.

### Plan 1 — Vocab + schema + paperconfig preprocessing + Biller 2018 paperconfig authoring

**Produces**

- `multiomics_kg/vocab/non_de_evidence.py` with unit tests
- `DerivedMetric` schema block + 3 binding edges + 3 measurement edges + `Experiment.compartment` in `config/schema_config.yaml`
- `iter_derived_metrics_tables` in `paperconfig_utils.py`
- Testable `validate_paperconfig_content(config, path)` extracted from CLI-only `validate()`
- Per-value_kind `derived_metrics_table` validation + `compartment` vocab check + `PAIRED_RNASEQ_PROTEOME` accepted + DE-only field relaxation
- `extract_rows_from_derived_metrics_table` in `build_gene_id_mapping.py` + dispatcher case
- `resolve_derived_metrics_entry` in `resolve_paper_ids.py` + main-loop integration
- **Biller 2018 paperconfig authoring** — add 3 `derived_metrics_table` entries (S4A: 4 boolean metrics; S4B: 2 boolean metrics; S5: 1 categorical metric). These sit in the paperconfig as no-ops until Plan 2's adapter consumes them (the existing `MultiClusterAdapter` only reads `gene_clusters`, so it ignores them; the new `MultiObservationsAdapter` doesn't exist yet)

**Tested against**

- Synthetic inline YAML strings in validator unit tests
- **Retrofitted Biller 2018 paperconfig** — the file produced by step 0 + this plan's authoring. Validator accepts the new entries and `build_gene_id_mapping` + `resolve_paper_ids` process them cleanly (harvest `id_columns`, write `_resolved.csv`)
- Existing paperconfig regression via `uv run python scripts/validate_paperconfig.py --all`

**Definition of done**

- All unit tests pass
- Retrofitted Biller 2018 paperconfig passes validation
- `bash scripts/prepare_data.sh --steps 3 4 --strains NATL2A MIT1002` regenerates resolved CSVs for S4A / S4B / S5
- `uv run python create_knowledge_graph.py --test` emits the exact same CSVs as step-0 baseline (schema loaded but no adapter emits to it yet; new entries in Biller 2018 paperconfig are no-ops)

### Plan 2 — Adapter + end-to-end build

**Produces**

- `multiomics_kg/adapters/observations_adapter.py` with `ObservationsAdapter` + `MultiObservationsAdapter`
- `create_knowledge_graph.py` wiring after `MultiClusterAdapter`
- **Synthetic numeric-DM paperconfig fixture** (~100 rows, exercises all 3 `value_kinds` including numeric) under `tests/fixtures/non_de/`; wired into an optional test-only `paperconfig_files.txt` entry

(Biller 2018 paperconfig authoring already landed in Plan 1; Plan 2's work is the code that consumes it.)

**Tested against**

- Biller 2018 real paperconfig (authored in Plan 1) — expected node + edge counts against step-0 baseline + S4A / S4B / S5 known row counts
- Synthetic numeric-DM fixture — full numeric DM path end-to-end

**Definition of done**

- `uv run python create_knowledge_graph.py --test` succeeds
- `docker compose up -d` reaches `import` stage clean (zero skipped relationships per `import.report`)
- Edge deltas vs step-0 baseline match expectations: `derived_metric_flags_gene` + `derived_metric_classifies_gene` + `derived_metric_quantifies_gene` counts > 0
- `Changes_expression_of` count unchanged vs step-0 baseline (DE path untouched)
- Adapter unit tests green

### Plan 3 — Post-import + KG validity + docs + skills

**Produces**

- Post-import Cypher additions (see §Post-import above) with byte-identical DE-only regression validated via `scripts/post-import-validate.sh`
- KG validity test additions per §Tests above
- Regenerated `snapshot_data.json` with retrofitted Biller 2018 evidence in new shape
- `/omics-edge-snapshot` skill extension (DerivedMetric counts)
- `docs/kg-changes/non-de-evidence-extension.md`
- Skill doc updates (paperconfig, cypher-queries, omics-edge-snapshot)
- `CLAUDE.md` key-facts update

**Definition of done**

- `pytest -m kg` green
- `pytest -m "not slow and not kg"` green
- `scripts/post-import-validate.sh` diff empty on DE-only subset vs step-0 baseline
- `/omics-edge-snapshot` reports expected DerivedMetric counts
- All docs / skill files committed

## Dependencies

```
Step 0 ──▶ Plan 1 ──▶ Plan 2 ──▶ Plan 3
                           │
                 (synthetic numeric-DM fixture file
                  created here; consumed by Plan 3
                  post-import tests against Neo4j)
```

- **Plan 1** uses inline YAML strings in its validator unit tests. No committed fixture file.
- **Plan 2** creates the synthetic numeric-DM paperconfig fixture file under `tests/fixtures/non_de/` and threads it through the adapter into CSV output.
- **Plan 3** relies on the Plan-2 fixture being in Neo4j so post-import Cypher (rank / percentile / bucket / significance) has real edges to compute on.

## Risks and open questions

- **Synthetic numeric-DM fixture design** — ~100 rows, organism / metric distribution / p-value distribution finalized at Plan 2 authoring. Fixture must be stable enough that Plan 3's bucket / percentile / significance assertions have deterministic expected values.
- **Biller 2018 post-retrofit gene-ID resolution** — the 3 retrofitted tables use `locus_tag` / `locus_tag_ncbi` columns that should already resolve (gene_id_mapping.json covers NATL2A and MIT1002). Verified before Plan 2 KG rebuild via `/check-gene-ids`.
- **Snapshot regeneration cadence** — step 0 and Plan 3 both regenerate. If the sample in `snapshot_data.json` happens to touch DerivedMetric nodes that only exist post-Plan-2, a second regeneration mid-Plan-3 is needed. Minor.
- **Parent spec's "one omics_type per Experiment" invariant** — Biller 2018's retrofitted Experiments are already correctly split (darkness_extended_darkness_natl2a_* carries only RNASEQ). No further splits needed.

## Handoff

This doc is execution scoping; implementation detail (task-by-task) lands in per-plan documents produced via the `superpowers:writing-plans` skill, one plan at a time, starting with Step 0.
