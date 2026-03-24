# KG Change: Experiment Node Stats + coculture_partner Fix

Date: 2026-03-22

## Summary

Precomputed expression stats added to Experiment nodes. Empty-string `coculture_partner` fixed to null on non-coculture experiments.

## New Experiment Node Properties (post-import computed)

| Property | Type | Description |
|---|---|---|
| `gene_count` | int | Total `Changes_expression_of` edges from this experiment |
| `significant_count` | int | Edges where `significant = 'significant'` |
| `time_point_count` | int | Number of distinct time points |
| `time_point_labels` | str[] | Ordered time point labels (e.g. `["2h", "4h", "24h"]`); `""` = no label (non-time-course) |
| `time_point_orders` | int[] | Parallel array of sort orders (1-indexed) |
| `time_point_hours` | float[] | Parallel array of hours values; `-1.0` = unknown conversion |
| `time_point_totals` | int[] | Parallel array of per-timepoint gene counts |
| `time_point_significants` | int[] | Parallel array of per-timepoint significant counts |

All time_point arrays are parallel (same length = `time_point_count`), ordered by `time_point_order`.

### Edge cases

- **Neo4j constraint**: arrays cannot contain nulls, so sentinel values are used:
  - `time_point_labels`: `""` (empty string) = no label (non-time-course experiments)
  - `time_point_hours`: `-1.0` = unknown hours conversion
- **Non-time-course experiments** (single implicit time point): arrays have one entry with `time_point_labels = [""]`, `time_point_hours = [-1.0]`
- **Experiments with 0 expression edges** (2 exist): `gene_count = 0`, `significant_count = 0`, `time_point_count = 0`, all arrays are empty `[]`

## Property Semantics Change

| Property | Before | After |
|---|---|---|
| `coculture_partner` on experiments without treatment organism | `""` (empty string) | **absent/null** |
| `coculture_partner` on coculture experiments | real value (unchanged) | real value (unchanged) |
| `coculture_partner` on viral experiments | real value (unchanged) | real value (unchanged) |

Both `coculture` and `viral` treatment types have a `coculture_partner` (the treatment organism — e.g., Alteromonas strain or Phage). All other treatment types have null `coculture_partner`.

Explorer code that checks `e.coculture_partner == ""` or uses string truthiness should switch to `IS NULL` / `IS NOT NULL`.

## No Changes To

- Node labels or types
- Edge labels or types
- Indexes
- Other node/edge properties

## KG Repo Files Changed

| File | Change |
|---|---|
| `multiomics_kg/adapters/omics_adapter.py` | `coculture_partner` only set when non-empty |
| `config/schema_config.yaml` | 8 new Experiment properties |
| `scripts/post-import.sh` | Experiment stats Cypher (2-pass with COALESCE sentinels) |
| `scripts/post-import.cypher` | Same (kept in sync) |
| `tests/kg_validity/test_post_import.py` | 7 new KG validity tests |
