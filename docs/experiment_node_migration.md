# Experiment Node Migration -- What Changed

## Summary

The experiment node redesign (2026-03-20) replaces `EnvironmentalCondition` nodes and the split expression edge types (`Condition_changes_expression_of`, `Coculture_changes_expression_of`) with a unified model: `Experiment` nodes linked to `Gene` nodes via a single `Changes_expression_of` edge type. Publication-to-source linkage (`Published_expression_data_about`) is replaced by `Has_experiment` edges, and coculture partner information is captured by `Tests_coculture_with` edges. Experiment-level metadata (organism, conditions, omics type) that was previously duplicated on every expression edge now lives on the Experiment node.

## Before vs After

| Aspect | Before | After |
|---|---|---|
| Expression source node | `EnvironmentalCondition` or `OrganismTaxon` | `Experiment` |
| Expression edge types | `Condition_changes_expression_of` + `Coculture_changes_expression_of` | `Changes_expression_of` |
| Publication link | `Published_expression_data_about` -> source nodes | `Has_experiment` -> `Experiment` |
| Coculture partner | Edge source was `OrganismTaxon` | `Tests_coculture_with` edge from `Experiment` to `OrganismTaxon` |
| Metadata on edges | `publications`, `omics_type`, `organism_strain`, `treatment_condition`, `control_condition`, `experimental_context`, `statistical_test`, `analysis_name` | Moved to `Experiment` node properties |
| Edge properties | fold change, p-value, direction, significance + metadata fields | fold change, p-value, direction, significance + `time_point`, `time_point_order`, `time_point_hours` |
| Condition nodes | ~40 `EnvironmentalCondition` nodes with `condition_category`, `medium`, `temperature` | Removed (absorbed into Experiment) |
| Paperconfig format | `environmental_conditions` block + per-analysis metadata | `experiments` block; analyses reference experiment via key |

## New graph structure

```
Publication --Has_experiment--> Experiment --Changes_expression_of--> Gene
                                    |
                                    +--Tests_coculture_with--> OrganismTaxon
                                       (coculture/viral experiments only)
```

## Experiment node properties

| Property | Type | Description |
|---|---|---|
| `name` | str | Human-readable experiment description |
| `organism_strain` | str | Target organism (e.g., "Prochlorococcus MED4") |
| `treatment_type` | str | Category for filtering (e.g., "nitrogen", "coculture", "light") |
| `treatment` | str | Treatment condition description |
| `control` | str | Control condition description |
| `experimental_context` | str | Additional context (medium, light, etc.) |
| `coculture_partner` | str | Name of coculture partner organism (empty for non-coculture) |
| `omics_type` | str | RNASEQ, PROTEOMICS, METABOLOMICS, or MICROARRAY |
| `statistical_test` | str | Statistical method (e.g., "DESeq2") |
| `is_time_course` | str | "true" or "false" |
| `medium` | str | Growth medium (e.g., "Pro99") |
| `temperature` | str | Growth temperature (e.g., "24C") |
| `light_condition` | str | Light regime (e.g., "continuous light", "diel cycle") |
| `light_intensity` | str | PAR value if reported |

Node ID format: `{doi}_{experiment_key}` (e.g., `10.1038/s41586-024-07246-9_n_limitation_med4_rnaseq`)

## Changes_expression_of edge properties

| Property | Type | Description |
|---|---|---|
| `time_point` | str | Time point label (e.g., "24h", "day 18") |
| `time_point_order` | int | Ordinal position within the experiment time series (0 for single-point) |
| `time_point_hours` | float | Numeric hours for the time point (enables arithmetic comparisons) |
| `log2_fold_change` | float | Log2 fold change value |
| `adjusted_p_value` | float | Adjusted p-value (may be null if not reported) |
| `expression_direction` | str | "up" or "down" |
| `significant` | str | "significant" or "not_significant" |

Post-import computed property:
- `rank_by_effect` (int) -- rank within each (experiment, time_point_order) group by descending |log2_fold_change| (1 = strongest effect)

## Cypher migration examples

### 1. Find expression changes for a gene

**Before:**
```cypher
MATCH (ec:EnvironmentalCondition)-[r:Condition_changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN ec.name, r.log2_fold_change, r.adjusted_p_value
UNION
MATCH (o:OrganismTaxon)-[r:Coculture_changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN o.preferred_name, r.log2_fold_change, r.adjusted_p_value
```

**After:**
```cypher
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN e.name, r.log2_fold_change, r.adjusted_p_value, e.treatment_type
```

### 2. Find coculture experiments

**Before:**
```cypher
MATCH (o:OrganismTaxon)-[r:Coculture_changes_expression_of]->(g:Gene)
RETURN o.preferred_name, g.locus_tag, r.log2_fold_change
```

**After:**
```cypher
MATCH (e:Experiment)-[:Tests_coculture_with]->(o:OrganismTaxon)
MATCH (e)-[r:Changes_expression_of]->(g:Gene)
RETURN e.name, o.preferred_name, g.locus_tag, r.log2_fold_change
```

### 3. Get publication for an expression edge

**Before:**
```cypher
MATCH (p:Publication)-[:Published_expression_data_about]->(ec:EnvironmentalCondition)
MATCH (ec)-[r:Condition_changes_expression_of]->(g:Gene)
RETURN p.title, ec.name, g.locus_tag
```

**After:**
```cypher
MATCH (p:Publication)-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(g:Gene)
RETURN p.title, e.name, g.locus_tag
```

### 4. Filter by treatment type

**Before:**
```cypher
MATCH (ec:EnvironmentalCondition {condition_category: 'nitrogen'})-[r:Condition_changes_expression_of]->(g:Gene)
RETURN ec.name, g.locus_tag, r.log2_fold_change
```

**After:**
```cypher
MATCH (e:Experiment {treatment_type: 'nitrogen'})-[r:Changes_expression_of]->(g:Gene)
RETURN e.name, g.locus_tag, r.log2_fold_change
```

### 5. Time series query

**Before:** (no structured time series support)
```cypher
MATCH (ec:EnvironmentalCondition)-[r:Condition_changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE r.time_point IS NOT NULL
RETURN r.time_point, r.log2_fold_change
```

**After:**
```cypher
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE e.is_time_course = 'true'
RETURN e.name, r.time_point, r.time_point_order, r.time_point_hours, r.log2_fold_change
ORDER BY r.time_point_order
```

## Paperconfig format changes

The `environmental_conditions` block in paperconfig.yaml is replaced by an `experiments` block. Each experiment defines the full experimental context. Individual `statistical_analyses` entries reference their parent experiment via the `experiment` key and carry only per-CSV/per-timepoint fields (`id`, `experiment`, `timepoint`, `timepoint_hours`, `name_col`, `logfc_col`, `adjusted_p_value_col`, `pvalue_threshold`).

Fields that moved from analyses to experiments: `type`/`omics_type`, `test_type`, `organism`, `control_condition`, `treatment_condition`, `experimental_context`, `name`.

## What was removed

- **`EnvironmentalCondition` nodes** (~40 nodes) -- absorbed into Experiment nodes
- **`Condition_changes_expression_of` edges** (~170K) -- replaced by `Changes_expression_of`
- **`Coculture_changes_expression_of` edges** (~17K) -- replaced by `Changes_expression_of`
- **`Published_expression_data_about` edges** -- replaced by `Has_experiment`
- **Expression edge properties** (moved to Experiment node): `publications`, `omics_type`, `organism_strain`, `treatment_condition`, `control_condition`, `experimental_context`, `statistical_test`, `analysis_name`
- **Paperconfig `environmental_conditions` block** -- replaced by `experiments` block
- **Paperconfig analysis fields**: `environmental_control_condition_id`, `environmental_treatment_condition_id` -- replaced by `experiment` reference key