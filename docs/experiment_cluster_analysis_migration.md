# What Changed: ClusteringAnalysis Nodes (April 2026)

## Summary

Gene expression clusters now have a **two-level hierarchy**: `ClusteringAnalysis` (the analysis that produced clusters) sits between `Publication` and `GeneCluster`. This replaces the flat structure where clusters were directly attached to publications.

## Why

Papers often contain multiple independent clustering analyses. Tolonen 2006, for example, ran K-means separately on MED4 and MIT9313 — producing two sets of clusters from the same publication. Previously, these were only distinguishable by properties on individual cluster nodes. Now each analysis is an explicit node, grouping its clusters and linking to the organism and (optionally) the experiment that produced the data.

## New Node: ClusteringAnalysis

| Property | Example | Description |
|---|---|---|
| `name` | "MED4 K-means N-starvation clusters" | Human-readable label |
| `organism_name` | "Prochlorococcus MED4" | Same vocabulary as Gene/Experiment |
| `cluster_method` | "K-means (K=9)" | Algorithm used |
| `cluster_type` | "response_pattern" | `diel_cycle` \| `time_series_dynamics` \| `response_pattern` |
| `cluster_count` | 9 | Number of child GeneCluster nodes |
| `total_gene_count` | 410 | Total genes in CSV (may exceed edge count if some IDs don't resolve) |
| `omics_type` | "MICROARRAY" | Platform |
| `treatment_type` | ["nitrogen_stress"] | Array of treatment categories |
| `background_factors` | [] | Array of background experimental factors |
| `treatment` | "N-starvation time course (0, 3, 6, 12, 24, 48h)" | Free-text condition |
| `light_condition` | "continuous light" | Light regime |
| `experimental_context` | "Custom Affymetrix microarray..." | Setup details |

**Current data:** 2 ClusteringAnalysis nodes (Tolonen 2006 MED4 and MIT9313).

## Changed Node: GeneCluster

**Added properties:**
- `id` — node identifier string (e.g., `cluster:msb4100087:med4_kmeans_nstarvation:1`)

**Removed properties:**
- `source_paper` — now reachable via ClusteringAnalysis → Publication path

**Denormalized from ClusteringAnalysis** (for direct query convenience):
`organism_name`, `cluster_method`, `cluster_type`, `treatment_type`, `background_factors`, `treatment`, `omics_type`, `light_condition`

**Descriptions:** `functional_description` and `behavioral_description` come from LLM extraction JSON files. Currently 3/16 clusters have descriptions populated (only clusters with extraction verdict = "pass"). The remaining clusters have empty description fields.

**Current data:** 16 GeneCluster nodes (9 MED4, 7 MIT9313), 928 gene membership edges.

## New Graph Structure

```
Publication
  └── PublicationHasClusteringAnalysis → ClusteringAnalysis
                                            ├── ClusteringAnalysisHasGeneCluster → GeneCluster
                                            │                                        └── Gene_in_gene_cluster → Gene
                                            └── ClusteringanalysisBelongsToOrganism → OrganismTaxon
```

## New Edges

| Edge | From → To | Count | Description |
|---|---|---|---|
| `PublicationHasClusteringAnalysis` | Publication → ClusteringAnalysis | 2 | One per analysis |
| `ClusteringAnalysisHasGeneCluster` | ClusteringAnalysis → GeneCluster | 16 | One per cluster |
| `ClusteringanalysisBelongsToOrganism` | ClusteringAnalysis → OrganismTaxon | 2 | One per analysis |
| `Gene_in_gene_cluster` | GeneCluster → Gene | 928 | Unchanged from before |

## Removed Edges

| Edge | Replacement |
|---|---|
| `Publication_has_gene_cluster` | `Publication → ClusteringAnalysis → GeneCluster` (two-hop) |
| `Genecluster_belongs_to_organism` | `ClusteringAnalysis → OrganismTaxon` (organism edge on analysis, not cluster) |

## Example Queries

### Find all clustering analyses for a publication
```cypher
MATCH (p:Publication)-[:PublicationHasClusteringAnalysis]->(ca:ClusteringAnalysis)
WHERE p.title CONTAINS 'Prochlorococcus'
RETURN ca.name, ca.organism_name, ca.cluster_count, ca.cluster_method
```

### Find clusters and their genes for a specific analysis
```cypher
MATCH (ca:ClusteringAnalysis {name: 'MED4 K-means N-starvation clusters'})
      -[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
      -[:Gene_in_gene_cluster]->(g:Gene)
RETURN gc.name, gc.member_count, gc.functional_description, count(g) as genes
ORDER BY gc.member_count DESC
```

### Find which organism an analysis belongs to
```cypher
MATCH (ca:ClusteringAnalysis)-[:ClusteringanalysisBelongsToOrganism]->(o:OrganismTaxon)
RETURN ca.name, o.preferred_name
```

### Find clusters by treatment type (using denormalized property)
```cypher
MATCH (gc:GeneCluster)
WHERE 'nitrogen_stress' IN gc.treatment_type
RETURN gc.name, gc.organism_name, gc.member_count, gc.functional_description
```

## Data Pipeline

Cluster data comes from three sources:

1. **Paperconfig YAML** — analysis-level metadata (organism, method, conditions)
2. **CSV files** — cluster membership (which genes belong to which cluster)
3. **Extraction JSON** — per-cluster descriptions (id, name, functional/behavioral descriptions) produced by the LLM extraction pipeline

The adapter reads all three and merges them. Per-cluster descriptions are only used when the extraction validation verdict is "pass".
