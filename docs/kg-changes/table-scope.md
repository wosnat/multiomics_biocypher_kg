# Table Scope on Experiment Nodes

## What Changed

Two new properties added to `Experiment` nodes:

| Property | Type | Description |
|---|---|---|
| `table_scope` | str | What genes the source DE table contains |
| `table_scope_detail` | str | Free-text clarification (optional) |

## Why

Differential expression tables from publications vary widely in completeness. Some papers report all ~1,700 detected genes; others report only 30-50 significant genes. Without tracking this, downstream queries and LLM agents cannot distinguish "gene not affected" from "gene not reported."

## `table_scope` Values

| Value | When to use |
|---|---|
| `all_detected_genes` | Table includes all genes detected in the assay, regardless of significance |
| `significant_any_timepoint` | Time-course table includes genes significant in at least one timepoint, with FC values reported across all timepoints |
| `significant_only` | Table includes only genes passing the author's significance threshold |
| `top_n` | Table includes top N genes ranked by effect size or significance |
| `filtered_subset` | Other author-defined subset (use `table_scope_detail` to explain) |

## How to Set

In `paperconfig.yaml`, add `table_scope` to each experiment block:

```yaml
experiments:
  my_experiment:
    name: "MED4 nitrogen limitation RNA-seq"
    organism: "Prochlorococcus MED4"
    table_scope: all_detected_genes
    # ... other fields
```

For ambiguous cases, add `table_scope_detail`:

```yaml
    table_scope: filtered_subset
    table_scope_detail: "Top 50% of genes by expression level"
```

## Impact on Queries

- `list_experiments` MCP tool can filter by `table_scope`
- Gene detail views can annotate "no expression data" differently based on whether the experiment could have reported the gene
- Cross-experiment comparisons can restrict to `all_detected_genes` for fair comparison

## Relationship to Existing Fields

`pvalue_threshold`, `logfc_threshold`, and `prefiltered` remain per-analysis fields. They control how individual rows are classified as significant during edge creation. `table_scope` is a higher-level property describing the completeness of the source data.
