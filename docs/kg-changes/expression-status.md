# KG Change: expression_status + Directional Significant Counts

Date: 2026-03-24

## Summary

Add `expression_status` as a derived property on `Changes_expression_of` edges. Split all precomputed significant counts into directional `significant_up_count` / `significant_down_count` on Experiment and Gene nodes.

## Motivation

The MCP `differential_expression_by_gene` tool needs a single enum field to filter and aggregate expression calls without combining two predicates (`significant` + `expression_direction`). Directional counts on Experiment and Gene nodes give immediate up/down signal without query-time aggregation.

## New Edge Property

| Edge | Property | Type | Derivation |
|---|---|---|---|
| `Changes_expression_of` (all ~188K edges) | `expression_status` | str | See rule below |

### Derivation rule

```
expression_status =
  "significant_up"   when significant = "significant" AND expression_direction = "up"
  "significant_down" when significant = "significant" AND expression_direction = "down"
  "not_significant"  otherwise
```

Computed by post-import Cypher. The source properties `significant` and `expression_direction` are retained (not removed).

## Property Renames on Experiment Nodes

| Before | After | Notes |
|---|---|---|
| `significant_count` (int) | `significant_up_count` (int) + `significant_down_count` (int) | Sum = old `significant_count` |
| `time_point_significants` (int[]) | `time_point_significant_up` (int[]) + `time_point_significant_down` (int[]) | Parallel arrays, same length as other tp arrays |

Unchanged: `gene_count`, `time_point_count`, `time_point_labels`, `time_point_orders`, `time_point_hours`, `time_point_totals`.

## Property Renames on Gene Nodes

| Before | After | Notes |
|---|---|---|
| `significant_expression_count` (int) | `significant_up_count` (int) + `significant_down_count` (int) | Sum = old `significant_expression_count` |

Unchanged: `expression_edge_count`.

## Relationships

- `not_significant_count` is not stored; derive as `gene_count - significant_up_count - significant_down_count`
- All time_point arrays remain parallel (same length = `time_point_count`)

## No Changes To

- Node labels or types
- Edge labels or types
- Indexes
- `significant` or `expression_direction` edge properties (retained)
- `rank_by_effect` edge property

## KG Repo Files Changed

| File | Change |
|---|---|
| `config/schema_config.yaml` | `expression_status` on edge; renamed properties on Experiment + Gene |
| `scripts/post-import.sh` | `expression_status` derivation step; directional counts in experiment stats + gene routing signals |
| `scripts/post-import.cypher` | Same (kept in sync) |
| `tests/kg_validity/test_expression.py` | 4 new tests for `expression_status` |
| `tests/kg_validity/test_post_import.py` | Updated: `significant_count` → `significant_up/down_count`, `time_point_significants` → `time_point_significant_up/down`, parallel array alignment check |

## MCP Impact

Tools that reference the old property names will need updating:
- `significant_count` → `significant_up_count` + `significant_down_count`
- `significant_expression_count` → `significant_up_count` + `significant_down_count`
- `time_point_significants` → `time_point_significant_up` + `time_point_significant_down`
- New `expression_status` available for filtering/aggregation
