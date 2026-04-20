# Biller 2018 Retrofit Removal — Step 0 Baseline

**Date:** 2026-04-20
**Plan:** `docs/superpowers/plans/2026-04-19-step0-biller-2018-retrofit-removal.md`
**Spec:** `docs/superpowers/specs/2026-04-19-non-de-evidence-biller-2018-slice.md`

## What changed

Removed 3 `gene_clusters` entries from Biller 2018's paperconfig and rebuilt the KG:

- `natl2a_periodicity` — sourced from Table S4A
- `mit1002_periodicity` — sourced from Table S4B
- `natl2a_darkness_survival` — sourced from Table S5

Plan 2 re-expresses this evidence through `derived_metrics_table` entries emitting `derived_metric_flags_gene` (periodicity) + `derived_metric_classifies_gene` (darkness survival) edges. Between step 0 and Plan 2, this evidence is absent from the KG — Cypher queries referencing the retrofitted cluster node IDs (`clustering_analysis:mSystems.00040-18:natl2a_periodicity` etc.) return 0 rows.

## Biller 2018 per-paper state

DOI: `10.1128/mSystems.00040-18`; Publication node id: `doi:10.1128/mSystems.00040-18`

| Metric | Before | After | Delta |
|---|---:|---:|---:|
| `Publication` node exists | 1 | 1 | 0 |
| `Has_experiment` children | 3 | 3 | 0 |
| `Changes_expression_of` edges | 12,074 | 12,074 | **0 (unchanged)** |
| `PublicationHasClusteringAnalysis` children | 3 | 0 | -3 |
| Retrofitted `ClusteringAnalysis` nodes (graph-wide) | 3 | 0 | -3 |
| Retrofitted `GeneCluster` nodes | 18 | 0 | -18 |
| Retrofitted `Gene_in_gene_cluster` edges | 2,852 | 0 | -2,852 |
| Retrofitted `ExperimentHasClusteringAnalysis` edges | 5 | 0 | -5 |

## `/omics-edge-snapshot` delta (cross-paper)

Biller 2018's `Changes_expression_of` count (12,074) is identical before and after. All other publications' counts are also unchanged. Project-wide DE path is untouched.

```
========================================================================
OMICS EDGE SNAPSHOT COMPARISON
  Before: 'before_biller2018_retrofit_removal'  (2026-04-19T17:49:43)
  After:  'after_biller2018_retrofit_removal'  (2026-04-20T04:40:10)
========================================================================

Total expression edges: 227,361 -> 227,361  (+0)

By target organism (all 24 organisms +0):
  Alteromonas (MarRef v6)                 652 →     652  (+0)
  Alteromonas macleodii EZ55            1,746 →   1,746  (+0)
  Alteromonas macleodii HOT1A3         51,132 →  51,132  (+0)
  Alteromonas macleodii MIT1002        28,822 →  28,822  (+0)
  Marinobacter (MarRef v6)              5,045 →   5,045  (+0)
  Prochlorococcus AS9601                1,887 →   1,887  (+0)
  Prochlorococcus MED4                 46,536 →  46,536  (+0)
  Prochlorococcus MIT9301                 234 →     234  (+0)
  Prochlorococcus MIT9303                  25 →      25  (+0)
  Prochlorococcus MIT9312               1,157 →   1,157  (+0)
  Prochlorococcus MIT9313              27,487 →  27,487  (+0)
  Prochlorococcus NATL1A                   81 →      81  (+0)
  Prochlorococcus NATL2A               29,034 →  29,034  (+0)
  Prochlorococcus SS120                 4,140 →   4,140  (+0)
  Pseudomonas putida KT2440             5,331 →   5,331  (+0)
  Ruegeria pomeroyi DSS-3               4,524 →   4,524  (+0)
  Shewanella sp. W3-18-1                2,960 →   2,960  (+0)
  Synechococcus CC9311                     27 →      27  (+0)
  Synechococcus PCC 7002                2,486 →   2,486  (+0)
  Synechococcus WH7803                  5,906 →   5,906  (+0)
  Synechococcus WH8102                  4,165 →   4,165  (+0)
  Synechococcus elongatus PCC 7942      3,154 →   3,154  (+0)
  Synechococcus elongatus UTEX 2973       163 →     163  (+0)
  Synechococcus sp. BL107                 667 →     667  (+0)

✓  No regressions — no publication lost edges
Unchanged: 31 publication(s)

  ✓  No genes lost across any publication

EXIT 0: graph looks healthy.
```

## Test-suite side effect

Removing 3 `ClusteringAnalysis` nodes + 18 `GeneCluster` nodes dropped the graph's cluster counts below the pre-existing `>= 15` / `>= 100` floors in `tests/kg_validity/test_clustering_analysis.py`. Step 0 lowered these 5 thresholds:

- `>= 15` → `>= 12` for `ClusteringAnalysis` node count
- `>= 100` → `>= 90` for `GeneCluster` node count
- `>= 15` → `>= 12` for `PublicationHasClusteringAnalysis` edge count
- `>= 100` → `>= 90` for `ClusteringAnalysisHasGeneCluster` edge count
- `>= 15` → `>= 12` for `ClusteringanalysisBelongsToOrganism` edge count

New floors are the current counts (no headroom) — tests now detect any further cluster loss. Plan 2 adds DerivedMetric (not ClusteringAnalysis) nodes, so these floors stay valid across subsequent slices.

## Artifacts committed alongside this document

- `before.json` — Biller-2018-specific Cypher dump before removal
- `after.json` — Biller-2018-specific Cypher dump after removal
- `capture.cypher` — the Cypher script used for both dumps (reusable; uses live-graph PascalCase labels for ClusteringAnalysis binding edges)
- `.claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json` — cross-paper DE snapshot
- `.claude/skills/omics-edge-snapshot/snapshots/after_biller2018_retrofit_removal.json` — cross-paper DE snapshot

## Follow-up

Plan 1 (vocab + schema + paperconfig preprocessing + Biller 2018 paperconfig authoring) adds the 3 `derived_metrics_table` entries that re-cover this evidence in the new shape; Plan 2 wires the adapter so the KG emits `derived_metric_flags_gene` + `derived_metric_classifies_gene` edges from those entries; Plan 3 adds post-import rollups and KG validity assertions.
