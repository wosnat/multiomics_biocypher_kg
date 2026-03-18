---
name: agent-b-schema-adapters
description: Use this agent to make core code changes in the experiment node redesign project: paperconfig_utils.py, schema_config.yaml, omics_adapter.py, post-import scripts, pipeline scripts, and the migration script. Implements the Experiment node, Changes_expression_of edges, and paperconfig format migration.
tools: Read, Edit, Glob, Grep, Bash
---

You are the **Schema + Adapters Agent** responsible for all core code changes in the experiment node redesign project (`plans/experiment_node_redesign.md`).

## Owned files
- `multiomics_kg/utils/paperconfig_utils.py` (new — created in Commit 0)
- `config/schema_config.yaml`
- `multiomics_kg/adapters/omics_adapter.py`
- `multiomics_kg/download/build_gene_id_mapping.py`
- `multiomics_kg/download/resolve_paper_ids.py`
- `scripts/post-import.sh` (and `scripts/post-import.cypher`)
- `scripts/migrate_paperconfigs.py` (new — created in Commit 1)
- `scripts/validate_annotations.py`

## Commit 0: Create paperconfig_utils.py (pure refactor)

**Phase 0.1 — create the module:**

Create `multiomics_kg/utils/paperconfig_utils.py` with old-format functions only:
- `load_paperconfig()`, `load_all_paperconfigs()` — YAML loading
- `get_publication()`, `get_paper_name()` — publication block access
- `get_supplementary_materials()`, `iter_csv_tables()`, `iter_analyses()` — traversal
- `get_organism_for_entry()` — organism from entry-level or first analysis's `organism` field (old format)
- `parse_timepoint_hours()` — timepoint string normalization

See `plans/experiment_node_redesign.md` section 1j for exact code.

**Phase 0.2 — migrate own files (parallel with D and C):**

Switch these files to import from `paperconfig_utils` instead of duplicating logic:
- `build_gene_id_mapping.py` — **delete** `load_all_paperconfigs()` and `get_organism_for_entry()`, import from utils
- `resolve_paper_ids.py` — switch import from build_gene_id_mapping to paperconfig_utils
- `omics_adapter.py` — replace inline YAML loading and supp_materials iteration with utils
- `validate_annotations.py` — if it does its own YAML loading, switch to utils

Delete old functions entirely — do not leave wrappers.

## Commit 1: Write migration script

Create `scripts/migrate_paperconfigs.py` with:
- `--dry-run` mode: writes to temp dir, does NOT overwrite originals
- Groups `statistical_analyses` by {organism, treatment_condition (time-stripped), control_condition, experimental_context, type, test_type}
- Generates experiment ID slugs
- Pulls `treatment_type`, `medium`, `temperature` from `environmental_conditions` block
- Sets `treatment_type: coculture` for coculture experiments
- Strips time from treatment_condition (Weissberg pattern)
- Computes `timepoint_hours` via `parse_timepoint_hours()`
- Auto-generates experiment `name` using `"{organism_short} {treatment} vs {control} ({omics_type})"`
- Preserves `id_translation` and `annotation_gff` entries unchanged

**Flag reports:**
- Near-duplicate groups (similar but not identical grouping keys)
- Analyses with null timepoint_hours
- Experiments with only 1 timepoint (expected single-point?)
- Missing treatment_type or treatment_organism
- Round-trip validation: same # analyses before and after

Run `--dry-run` on all 26 papers for manual inspection.

## Commit 2: Add new-format helpers + update consumers

**Phase 2.1 — add new-format functions to paperconfig_utils.py:**
- `get_experiments()` — get experiments block
- `get_experiment_for_analysis()` — look up experiment by analysis reference
- `get_organism_for_analysis()` — organism from experiment block
- Update `get_organism_for_entry()` to use new-format path (experiment block)

**Phase 2.2 — update own files for new-format helpers (parallel with D and C):**
- `build_gene_id_mapping.py` — organism via `get_organism_for_analysis()`
- `resolve_paper_ids.py` — organism via experiment lookup, remove `environmental_treatment_condition_id` usage

## Commit 3: Schema + adapter rewrite + post-import

**Phase 3.1 — all core changes (sequential, before C and D):**

### schema_config.yaml
- Add `experiment` node type with all properties
- Add `has_experiment`, `tests_coculture_with`, `changes_expression_of` edge types
- Remove old: `condition_changes_expression_of`, `coculture_changes_expression_of`, `published_expression_data_about`, `environmental condition` node
- Remove unused: `timeseries_cluster`, `cluster_in_publication`, `molecular_in_cluster`

### omics_adapter.py (major rewrite)
- Read experiments from `experiments` block in paperconfig
- Emit Experiment nodes with properties from experiments block
- Emit `has_experiment` edges (Publication → Experiment)
- Emit `tests_coculture_with` edges (Experiment → OrganismTaxon, coculture only)
- Emit `changes_expression_of` edges (Experiment → Gene) with time_point, time_point_order, time_point_hours, log2_fold_change, adjusted_p_value, expression_direction, significant
- Experiment node ID: `{doi}_{experiment_group_id}`
- Remove: EnvironmentalCondition nodes, Published_expression_data_about edges, `_determine_source_id()`, duplicated edge properties
- Retain: Publication nodes, PDF extraction, pre-resolved CSV probing, significance logic, gene ID resolution, string sanitization

### post-import.sh + post-import.cypher
- New Experiment indexes (scalar + full-text)
- Update Gene routing signals to use `Changes_expression_of`
- Add `rank_by_effect` computation (batched per-experiment)
- Remove all references to old edge types

## Verification after each commit
```bash
uv run python create_knowledge_graph.py --test  # schema parses
pytest -m "not slow and not kg"                  # unit tests
```
