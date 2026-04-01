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
- Example: `clustering_analysis:msb4100087:med4_kmeans_nstarvation`
- The entry key in the paperconfig must be short, meaningful, and unique within the paper (e.g., `med4_kmeans_nstarvation` not `cluster_table_med4`). This is a paperconfig creation requirement.

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
| `cluster_type` | str | paperconfig | Category enum: `"diel_cycle"`, `"time_series_dynamics"`, `"response_pattern"` |
| `experimental_context` | str | paperconfig | Setup details |

## GeneCluster Node Changes

**Node ID format (updated):** `cluster:{doi_short}:{analysis_key}:{csv_cluster_key}`
- Example: `cluster:msb4100087:med4_kmeans_nstarvation:1`
- The `csv_cluster_key` is the raw value from the CSV `cluster_col` (e.g., `"1"`, `"A"`)
- The extracted `id` from JSON (e.g., `up_n_transport`) is a property, not part of the node ID

**GeneCluster-specific properties** (from extraction JSON + CSV):
- `id` — short snake_case identifier from extraction JSON (e.g., `up_n_transport`), property not part of node ID
- `name` — from extraction JSON
- `member_count` — computed from CSV
- `functional_description` — from extraction JSON
- `behavioral_description` — from extraction JSON
- `peak_time_hours` — from extraction JSON (diel only, nullable)
- `period_hours` — from extraction JSON (diel only, nullable)

**Removed:** `source_paper` (reachable via ClusteringAnalysis → Publication), `cluster_type` (moved to analysis level).

**Denormalized copies from ClusteringAnalysis** (for query convenience): `organism_name`, `cluster_method`, `cluster_type`, `omics_type`, `treatment_type`, `treatment`, `light_condition`, `experimental_context`. This matches the existing pattern of denormalizing `organism_name` onto Gene nodes.

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

Each `type: gene_clusters` entry maps 1:1 to a ClusteringAnalysis node. Changes from current format:

### New required field: `name`

```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"  # human-readable name for ClusteringAnalysis node
  ...
```

### New optional field: `experiments`

```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  experiments: [n_starvation_med4]  # references experiment keys in same paperconfig
  ...
```

### Full example

```yaml
med4_kmeans_nstarvation:
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
  cluster_type: "response_pattern"
  # No clusters: block — cluster keys derived from unique values in CSV cluster_col.
  # Extraction pipeline reads cluster keys from CSV, not paperconfig.
  # Per-cluster data (id, name, descriptions, peak_time_hours, period_hours)
  # comes entirely from extraction JSON output.
```

The entry key (`med4_kmeans_nstarvation`) is used as the analysis identifier for node ID construction.

## Extraction JSON Integration

The LLM-based extraction pipeline (`multiomics_kg/extraction/cluster/`) produces per-analysis JSON files alongside each paperconfig. The adapter reads these to populate cluster descriptions, IDs, and other extracted fields — making the extraction output a persistent data source alongside the paperconfig.

### File convention

For a paperconfig entry with key `med4_kmeans_nstarvation`, the extraction output is:
```
data/.../cluster_extraction_med4_kmeans_nstarvation.json
```
Located in the same directory as the paperconfig. The filename is `cluster_extraction_{entry_key}.json`.

### Extraction JSON structure (read by adapter)

```json
{
  "metadata": {
    "paper": "Tolonen 2006",
    "organism": "Prochlorococcus MED4",
    "table_key": "med4_kmeans_nstarvation"
  },
  "stage2_results": {
    "1": {
      "id": "med4_up_transport_amino_acid_biosynthesis",
      "name": "MED4 cluster 1 (up, N transport)",
      "functional_description": "...",
      "behavioral_description": "...",
      "peak_time_hours": null,
      "period_hours": null
    }
  },
  "stage3_validation": {
    "1": {
      "verdict": "pass",
      "explanation": "..."
    }
  }
}
```

### Data sources by level

| Level | Source | Fields |
|---|---|---|
| **ClusteringAnalysis** | paperconfig | `name`, `organism`, `cluster_method`, `cluster_type`, `omics_type`, `treatment_type`, `treatment`, `light_condition`, `experimental_context`, `experiments` |
| **ClusteringAnalysis** | computed from CSV | `cluster_count`, `total_gene_count` |
| **GeneCluster** | CSV `cluster_col` | cluster keys (unique values), `member_count` |
| **GeneCluster** | extraction JSON | `id`, `name`, `functional_description`, `behavioral_description`, `peak_time_hours`, `period_hours` |
| **GeneCluster** | denormalized from analysis | `organism_name`, `cluster_method`, `cluster_type`, `omics_type`, `treatment_type`, `treatment`, `light_condition`, `experimental_context` |
| **Gene_in_gene_cluster** | CSV | `membership_score` (from `score_col`, optional), `p_value` (optional) |

The adapter only uses extraction results where `stage3_validation` verdict is `"pass"` (skips failed/flagged clusters, leaving them for manual review). Clusters without a passing extraction result get empty description fields.

**Extraction pipeline changes required:** The synthesis prompt (`prompts.py` `SYNTHESIS_PROMPT`) currently produces `id`, `functional_description`, `behavioral_description`. Add three more output fields:
- `name` — short human-readable cluster label (e.g., "MED4 cluster 1 (up, N transport)")
- `peak_time_hours` — peak expression time (float, diel clusters only; null otherwise)
- `period_hours` — oscillation period (float, diel clusters only; null otherwise)

Update `synthesis.py` and `validation.py` accordingly.

### Adapter lookup logic

For each `type: gene_clusters` entry:
1. Read unique cluster keys from CSV `cluster_col`
2. Resolve the extraction JSON path from the paperconfig directory + entry key
3. If the JSON file exists, load `stage2_results` and `stage3_validation`; populate per-cluster properties from passing results
4. If the JSON file doesn't exist, cluster nodes get empty description fields

## Adapter Changes

All changes are confined to `multiomics_kg/adapters/cluster_adapter.py`.

### ClusterAdapter.get_nodes()

Yields two node types:
1. **ClusteringAnalysis nodes** — one per `type: gene_clusters` entry. Properties built from entry-level fields plus computed `cluster_count` and `total_gene_count`.
2. **GeneCluster nodes** — one per unique value in CSV `cluster_col`. Per-cluster properties from extraction JSON, denormalized analysis-level properties from paperconfig.

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

## Paperconfig Skill and Validator Updates

The paperconfig skill (`.claude/skills/paperconfig/SKILL.md`) and validator (`validate_paperconfig.py`) must be updated first, before any adapter or schema changes. They are the entry point for creating new paperconfigs.

### Validator changes (`validate_paperconfig.py`)

1. **Updated required fields:**
   ```python
   REQUIRED_CLUSTER_TABLE_FIELDS = [
       "name", "filename", "organism", "gene_id_col", "cluster_col", "cluster_type",
   ]
   ```
   - `name` added (human-readable label for the ClusteringAnalysis node)
   - `cluster_type` added (moved from per-cluster to analysis level; enum: `diel_cycle`, `time_series_dynamics`, `response_pattern`)
   - `clusters` removed (no longer in paperconfig; cluster keys derived from CSV)

2. **Entry key validation:** Warn if entry key looks like a generic name (e.g., starts with `cluster_table_`). Entry keys should be short, meaningful IDs (e.g., `med4_kmeans_nstarvation`).

3. **Optional field validation:**
   - `experiments` — if present, must be a list of strings. Validate that each referenced experiment key exists in the paperconfig's `experiments:` block (cross-reference check).

4. **Extraction JSON validation:** When a `cluster_extraction_*.json` file exists alongside the paperconfig, validate:
   - `metadata.table_key` matches the entry key
   - Every cluster key in `stage2_results` corresponds to a value in the CSV `cluster_col`
   - Stage 3 `verdict` is reported (warn on `"fail"` clusters)

5. **Removed validation:** Per-cluster `cluster_type`, `name`, `functional_description`, `behavioral_description` checks (these no longer live in the paperconfig).

### Skill changes (`SKILL.md`)

1. Update `type: gene_clusters` documentation:
   - `name` is required (analysis-level, for ClusteringAnalysis node)
   - `cluster_type` is required at analysis level (enum: `diel_cycle`, `time_series_dynamics`, `response_pattern`)
   - `experiments` is optional (list of experiment keys from same paperconfig)
   - Entry key must be short and meaningful (used in node ID)
2. Remove `clusters:` block documentation — per-cluster data comes from extraction JSON
3. Add note about extraction JSON as the source for per-cluster descriptions, IDs, and diel parameters
4. Document the data flow: paperconfig (analysis metadata) + CSV (membership) + extraction JSON (per-cluster properties)

## Scope

### In scope
- Update paperconfig skill and validator (first step)
- New ClusteringAnalysis node type and 4 new edge types
- Remove 2 old edge types
- Update cluster_adapter.py (including extraction JSON merge logic)
- Update schema_config.yaml
- Migrate existing paperconfigs (rename entry keys, add `name` + `cluster_type`, remove `clusters:` block)
- Update post-import.sh indexes
- Update tests

### Out of scope
- Changes to omics_adapter.py or any other adapter
- MCP tool changes (future work)
- The extraction pipeline itself (already built separately), except: adding `name`/`peak_time_hours`/`period_hours` to synthesis output, and reading cluster keys from CSV instead of paperconfig
- Post-import computed properties on ClusteringAnalysis

## Implementation Order

Tolonen 2006 is the validation case throughout — each step is verified against it before moving on.

1. **Paperconfig skill + validator** — update to recognize new format (`name`, `cluster_type` at analysis level, `experiments`, no `clusters:` block, extraction JSON validation)
2. **Tolonen 2006 paperconfig** — migrate to new format (meaningful entry keys, `name`, `cluster_type`, remove `clusters:` block). Update extraction JSON `table_key` to match new entry keys.
3. **Extraction pipeline** — add `name`, `peak_time_hours`, `period_hours` to synthesis output. Update to read cluster keys from CSV instead of paperconfig `clusters:` block. Re-run on Tolonen 2006 to produce updated JSON.
4. **Schema config** — add ClusteringAnalysis node, new edges, remove old edges
5. **Adapter** — restructure cluster_adapter.py for new node/edge types + extraction JSON read
6. **Post-import** — update indexes
7. **Tests** — update/add unit tests, verify with Tolonen 2006 end-to-end

### Separate follow-up task
- Create paperconfigs + CSVs for additional papers (Zinser 2009, Alonso 2023, Bagby 2015, Wang 2014)
- Run extraction pipeline on each
- Review and commit results
