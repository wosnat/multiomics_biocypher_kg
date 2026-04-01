# ClusteringAnalysis Node Design Spec

**Date:** 2026-03-31
**Status:** Approved

## Problem

The current GeneCluster framework has flat structure: clusters sit directly under Publications with no intermediate grouping. Papers like Tolonen 2006 contain multiple distinct clustering analyses (MED4 K-means, MIT9313 K-means) that are only distinguishable by properties on individual cluster nodes. This loses the structural relationship between clusters that belong to the same analysis and prevents linking clustering analyses to the Experiment nodes that produced the underlying data.

## Solution

Add a `ClusteringAnalysis` intermediate node between Publication and GeneCluster. Each `type: gene_clusters` entry in a paperconfig becomes one ClusteringAnalysis node, grouping the individual clusters underneath it.

## New Node: ClusteringAnalysis

**BioCypher type:** `clustering analysis` (renders as `ClusteringAnalysis`)

**Node ID format:** `clustering_analysis:{doi_short}:{entry_key}`
- Example: `clustering_analysis:msb4100087:cluster_table_med4`

**Properties:**

| Property | Type | Source | Description |
|---|---|---|---|
| `name` | str | paperconfig (new required field) | Human-readable label (e.g., "MED4 K-means N-starvation clusters") |
| `organism_name` | str | paperconfig `organism` | Matches Gene/Experiment vocabulary |
| `cluster_method` | str | paperconfig | Algorithm (e.g., "K-means (K=9)", "Mfuzz") |
| `cluster_count` | int | computed at build time | Number of GeneCluster children |
| `total_gene_count` | int | computed at build time | Total genes across all clusters |
| `omics_type` | str | paperconfig | MICROARRAY, RNASEQ, etc. |
| `treatment_type` | str[] | paperconfig | Category array for filtering |
| `treatment` | str | paperconfig | Free-text condition description |
| `light_condition` | str | paperconfig | Light regime |
| `experimental_context` | str | paperconfig | Setup details |

## GeneCluster Node Changes

No properties removed. Analysis-level properties (`organism_name`, `cluster_method`, `omics_type`, `treatment_type`, `treatment`, `light_condition`, `experimental_context`) remain as denormalized copies on GeneCluster nodes for query convenience. This matches the existing pattern of denormalizing `organism_name` onto Gene nodes.

## Edge Changes

### New edges (4)

| Edge type | Source → Target | Cardinality | Properties |
|---|---|---|---|
| `Publication_has_clustering_analysis` | Publication → ClusteringAnalysis | 1:N | none |
| `Clustering_analysis_has_gene_cluster` | ClusteringAnalysis → GeneCluster | 1:N | none |
| `Clusteringanalysis_belongs_to_organism` | ClusteringAnalysis → OrganismTaxon | N:1 | none |
| `Experiment_has_clustering_analysis` | Experiment → ClusteringAnalysis | N:M (optional) | none |

### Removed edges (2)

| Edge type | Replacement |
|---|---|
| `Publication_has_gene_cluster` | Publication → ClusteringAnalysis → GeneCluster path |
| `Genecluster_belongs_to_organism` | `Clusteringanalysis_belongs_to_organism` (organism relationship belongs at analysis level) |

### Unchanged edges

| Edge type | Notes |
|---|---|
| `Gene_in_gene_cluster` | GeneCluster → Gene, stays exactly as-is |

### Experiment linkage

The `Experiment_has_clustering_analysis` edge is optional. It is emitted when the paperconfig `type: gene_clusters` entry includes an `experiments` list referencing experiment keys from the same publication. A clustering analysis may link to zero, one, or multiple experiments (N:M). The adapter resolves experiment keys to Experiment node IDs (`{doi}_{experiment_key}`).

## Paperconfig Changes

Each `type: gene_clusters` entry already maps 1:1 to a ClusteringAnalysis node. Two fields are added:

### New required field: `name`

```yaml
cluster_table_med4:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"  # human-readable name for ClusteringAnalysis node
  ...
```

### New optional field: `experiments`

```yaml
cluster_table_med4:
  type: gene_clusters
  experiments: [n_starvation_med4]  # references experiment keys in same paperconfig
  ...
```

### Full example

```yaml
cluster_table_med4:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"
  filename: "data/.../med4_kmeans_clusters.csv"
  organism: "Prochlorococcus MED4"
  cluster_method: "K-means (K=9)"
  omics_type: MICROARRAY
  treatment_type: ["nitrogen_stress"]
  treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
  light_condition: "continuous light"
  experimental_context: "Custom Affymetrix microarray..."
  experiments: [n_starvation_med4]
  gene_id_col: "gene_id"
  cluster_col: "cluster"
  clusters:
    "1":
      id: "up_n_transport"
      name: "MED4 cluster 1 (up, N transport)"
      cluster_type: "stress_response"
      functional_description: "N transport genes..."
      behavioral_description: "Rapid upregulation within 3h..."
```

The entry key (`cluster_table_med4`) is used as the analysis identifier for node ID construction.

## Adapter Changes

All changes are confined to `multiomics_kg/adapters/cluster_adapter.py`.

### ClusterAdapter.get_nodes()

Yields two node types:
1. **ClusteringAnalysis nodes** — one per `type: gene_clusters` entry. Properties built from entry-level fields plus computed `cluster_count` and `total_gene_count`.
2. **GeneCluster nodes** — same as today, retaining denormalized copies of analysis-level properties.

### ClusterAdapter.get_edges()

Yields five edge types:
1. `Publication_has_clustering_analysis` — one per analysis
2. `Clustering_analysis_has_gene_cluster` — one per cluster
3. `Clusteringanalysis_belongs_to_organism` — one per analysis
4. `Experiment_has_clustering_analysis` — one per referenced experiment (when `experiments` field is present)
5. `Gene_in_gene_cluster` — unchanged

### MultiClusterAdapter

No structural changes. Passes through the new node/edge types from each ClusterAdapter.

## Schema Config Changes

In `config/schema_config.yaml`:

1. Add `clustering analysis` node type with all properties listed above
2. Add four new edge types: `publication has clustering analysis`, `clustering analysis has gene cluster`, `clusteringanalysis belongs to organism`, `experiment has clustering analysis`
3. Remove two edge types: `publication has gene cluster`, `genecluster belongs to organism`

## Post-Import and Indexes

### New indexes

- Scalar: `clustering_analysis_organism_idx` on `ClusteringAnalysis(organism_name)`
- Scalar: `clustering_analysis_method_idx` on `ClusteringAnalysis(cluster_method)`
- Full-text: `clusteringAnalysisFullText` on ClusteringAnalysis `name`, `treatment`, `experimental_context`

### Removed indexes

- `gene_cluster_organism_idx` — organism queries go through ClusteringAnalysis; the denormalized property on GeneCluster still allows direct filtering if needed

### Unchanged indexes

- `geneClusterFullText` on GeneCluster `name`, `functional_description`, `behavioral_description`

### No post-import computed properties

`cluster_count` and `total_gene_count` are computed by the adapter at build time, not by post-import Cypher.

## Graph Structure Summary

```
Publication
  ├── Has_experiment → Experiment
  │                      └── Experiment_has_clustering_analysis → ClusteringAnalysis (optional, N:M)
  └── Publication_has_clustering_analysis → ClusteringAnalysis
                                              ├── Clustering_analysis_has_gene_cluster → GeneCluster
                                              │                                           └── Gene_in_gene_cluster → Gene
                                              └── Clusteringanalysis_belongs_to_organism → OrganismTaxon
```

## Scope

### In scope
- New ClusteringAnalysis node type and 4 new edge types
- Remove 2 old edge types
- Update cluster_adapter.py
- Update schema_config.yaml
- Update existing paperconfigs (add `name` field to `type: gene_clusters` entries)
- Update post-import.sh indexes
- Update tests

### Out of scope
- Changes to omics_adapter.py or any other adapter
- MCP tool changes (future work)
- Cluster description extraction pipeline (separate effort)
- Post-import computed properties on ClusteringAnalysis
