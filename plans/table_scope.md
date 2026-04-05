# Table Scope Enhancement for Experiment Nodes

## Motivation

Differential expression tables from publications vary widely in completeness:
- Some include **all detected genes** (~100% of genome)
- Some include only **genes significant in any timepoint** (time-course papers)
- Some include only **significant genes** at the author's threshold
- Some include a **top N** or other filtered subset

Without tracking this, cross-experiment queries conflate "gene not affected" with "gene not reported." A gene absent from a 46-row table (Lindell 2007) is very different from a gene absent from a 1,715-row table (Aharonovich 2016).

## New Properties

### On Experiment nodes (schema + adapter)

| Property | Type | Required | Description |
|---|---|---|---|
| `table_scope` | str | yes | What genes the source table contains |
| `table_scope_detail` | str | no | Free-text clarification for edge cases |

### `table_scope` values

| Value | Meaning |
|---|---|
| `all_detected_genes` | All genes detected in the assay |
| `significant_any_timepoint` | Genes significant in at least one timepoint, all timepoints reported |
| `significant_only` | Only genes passing the author's significance threshold |
| `top_n` | Top N genes ranked by effect size or significance |
| `filtered_subset` | Other author-defined subset |

### Existing fields (no change)

`pvalue_threshold`, `logfc_threshold`, and `prefiltered` remain per-analysis fields on the paperconfig statistical_analyses entries. They control significance determination during edge creation and are not duplicated on Experiment nodes.

## Implementation Steps

### 1. Schema (`config/schema_config.yaml`)
- Add `table_scope: str` and `table_scope_detail: str` to the `experiment` node properties

### 2. Adapter (`multiomics_kg/adapters/omics_adapter.py`)
- Read `table_scope` and `table_scope_detail` from the experiment block in paperconfig
- Add to `exp_props` dict in `get_nodes()`
- Default `table_scope` to empty string if not set (backward compat during rollout)

### 3. Paperconfig validation (`.claude/skills/paperconfig/validate_paperconfig.py`)
- Add `table_scope` to `RECOMMENDED_EXPERIMENT_FIELDS`
- Validate `table_scope` value is one of the allowed enum values
- Warn if `table_scope` is missing

### 4. Paperconfig skill (`.claude/skills/paperconfig/SKILL.md`)
- Add `table_scope` to the experiment block template
- Document allowed values

### 5. Tests
- Update unit tests for omics adapter to include `table_scope`
- Add validation test for enum values

### 6. Documentation
- Create `docs/kg-changes/table-scope.md` (what's changed doc)
- Update CLAUDE.md Experiment node description

## Out of Scope (Task 2)

Auditing all 24 paperconfigs and assigning `table_scope` values is a separate task requiring per-paper judgment.
