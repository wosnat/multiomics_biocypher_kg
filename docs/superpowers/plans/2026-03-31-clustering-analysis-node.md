# ClusteringAnalysis Node Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a ClusteringAnalysis intermediate node between Publication and GeneCluster, removing the `clusters:` block from paperconfigs, and teaching the adapter to read per-cluster data from extraction JSON files.

**Architecture:** Each `type: gene_clusters` paperconfig entry becomes a ClusteringAnalysis node. GeneCluster nodes are keyed by CSV cluster column values. Per-cluster properties (id, name, descriptions, diel params) come from extraction JSON. Analysis-level properties come from the paperconfig and are denormalized onto GeneCluster nodes.

**Tech Stack:** Python, YAML, Neo4j Cypher, pytest

**Spec:** `docs/superpowers/specs/2026-03-31-clustering-analysis-node-design.md`

---

## File Map

| File | Action | Purpose |
|---|---|---|
| `.claude/skills/paperconfig/validate_paperconfig.py` | Modify | Update required fields, cluster_type enum, remove clusters: validation |
| `.claude/skills/paperconfig/SKILL.md` | Modify | Update gene_clusters docs |
| `data/.../tolonen 2006/paperconfig.yaml` | Modify | Migrate to new format |
| `data/.../tolonen 2006/cluster_extraction_med4.json` | Modify | Update table_key |
| `data/.../tolonen 2006/cluster_extraction_mit9313.json` | Modify | Update table_key |
| `multiomics_kg/extraction/cluster/prompts.py` | Modify | Add name/peak_time_hours/period_hours to synthesis prompt |
| `multiomics_kg/extraction/cluster/pipeline.py` | Modify | Output filename uses entry_key instead of organism_short |
| `config/schema_config.yaml` | Modify | Add ClusteringAnalysis node, new edges, remove old edges |
| `multiomics_kg/adapters/cluster_adapter.py` | Modify | Major restructure: two node types, five edge types, extraction JSON read |
| `scripts/post-import.sh` | Modify | Add ClusteringAnalysis indexes |
| `scripts/post-import.cypher` | Modify | Same (kept in sync) |
| `tests/test_cluster_adapter.py` | Modify | Update all tests for new structure |

---

## Task 1: Update Paperconfig Validator

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py:134-151` (constants), `:631-700` (validation logic)

- [ ] **Step 1: Update constants**

In `validate_paperconfig.py`, replace the cluster constants (around lines 134-151):

```python
VALID_CLUSTER_TYPES = {
    "diel_cycle",
    "time_series_dynamics",
    "response_pattern",
}

REQUIRED_CLUSTER_TABLE_FIELDS = [
    "name", "filename", "organism", "gene_id_col", "cluster_col", "cluster_type",
]

# Removed: REQUIRED_CLUSTER_FIELDS, RECOMMENDED_CLUSTER_FIELDS
# Per-cluster data comes from extraction JSON, not paperconfig
```

- [ ] **Step 2: Rewrite gene_clusters validation function**

Replace the gene_clusters validation block (starts around line 631) to:
- Validate `cluster_type` against `VALID_CLUSTER_TYPES`
- Validate `experiments` field if present: must be a list, each entry must exist in the paperconfig `experiments:` block
- Warn if entry key starts with `cluster_table_` (should be meaningful)
- Remove all `clusters:` block validation (no more per-cluster checks)
- Add extraction JSON validation: if `cluster_extraction_{entry_key}.json` exists alongside the paperconfig, check `metadata.table_key` matches

```python
def _validate_gene_clusters_entry(key, table, config, paperconfig_dir):
    """Validate a type: gene_clusters entry."""
    errors = []
    warnings = []

    # Required fields
    for field in REQUIRED_CLUSTER_TABLE_FIELDS:
        if field not in table:
            errors.append(f"  [{key}] Missing required field: {field}")

    # Organism check
    organism = table.get("organism", "")
    if organism:
        print(f"  Organism: {organism}")
        if organism not in CANONICAL_ORGANISMS:
            warnings.append(f"  [{key}] Organism '{organism}' not in canonical list")

    # cluster_type enum
    ct = table.get("cluster_type", "")
    if ct and ct not in VALID_CLUSTER_TYPES:
        warnings.append(f"  [{key}] cluster_type '{ct}' not in {VALID_CLUSTER_TYPES}")

    # omics_type
    ot = table.get("omics_type", "")
    if ot and ot not in VALID_TYPES:
        warnings.append(f"  [{key}] omics_type '{ot}' not in {VALID_TYPES}")

    # treatment_type must be a list
    tt = table.get("treatment_type")
    if tt is not None:
        if not isinstance(tt, list):
            errors.append(f"  [{key}] treatment_type must be a list, got {type(tt).__name__}")
        else:
            for t in tt:
                if t not in CANONICAL_TREATMENT_TYPES:
                    warnings.append(f"  [{key}] treatment_type '{t}' not canonical")

    # Entry key naming
    if key.startswith("cluster_table_"):
        warnings.append(f"  [{key}] Entry key should be a meaningful ID (e.g., 'med4_kmeans_nstarvation'), not generic")

    # experiments cross-reference
    experiments_ref = table.get("experiments")
    if experiments_ref is not None:
        if not isinstance(experiments_ref, list):
            errors.append(f"  [{key}] experiments must be a list")
        else:
            pub_experiments = config.get("publication", {}).get("experiments", {})
            for exp_key in experiments_ref:
                if exp_key not in pub_experiments:
                    errors.append(f"  [{key}] experiments references '{exp_key}' not found in publication.experiments")

    # CSV validation
    filename = table.get("filename", "")
    if filename:
        csv_path = Path(filename)
        if not csv_path.is_absolute():
            csv_path = _find_project_root() / csv_path
        if csv_path.exists():
            try:
                sep = table.get("separator", ",")
                skip = table.get("skip_rows", 0)
                df = pd.read_csv(csv_path, sep=sep, skiprows=skip, nrows=5)
                cols = set(df.columns)
                for col_field in ("gene_id_col", "cluster_col", "score_col"):
                    col_name = table.get(col_field)
                    if col_name and col_name not in cols:
                        errors.append(f"  [{key}] {col_field}='{col_name}' not found in CSV columns: {sorted(cols)}")
            except Exception as e:
                warnings.append(f"  [{key}] Could not read CSV: {e}")
        else:
            warnings.append(f"  [{key}] CSV file not found: {csv_path}")

    # Extraction JSON validation
    if paperconfig_dir:
        json_path = paperconfig_dir / f"cluster_extraction_{key}.json"
        if json_path.exists():
            try:
                with open(json_path) as f:
                    extraction = json.load(f)
                meta_key = extraction.get("metadata", {}).get("table_key", "")
                if meta_key != key:
                    warnings.append(f"  [{key}] Extraction JSON table_key='{meta_key}' does not match entry key")
                stage3 = extraction.get("stage3_validation", {})
                for ck, cv in stage3.items():
                    if isinstance(cv, dict) and cv.get("verdict") == "fail":
                        warnings.append(f"  [{key}] Extraction cluster '{ck}' has verdict=fail")
            except Exception as e:
                warnings.append(f"  [{key}] Could not read extraction JSON: {e}")

    # Warn if clusters: block still present (deprecated)
    if "clusters" in table:
        warnings.append(f"  [{key}] 'clusters' block is deprecated — per-cluster data comes from extraction JSON")

    return errors, warnings
```

- [ ] **Step 3: Run the validator on tolonen to verify it works with old format (expect warnings)**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

Expected: warnings about missing `name`, `cluster_type` at table level, deprecated `clusters` block, and generic entry key names.

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py
git commit -m "feat: update paperconfig validator for ClusteringAnalysis format

- Required fields: name, cluster_type (at analysis level), remove clusters
- New cluster_type enum: diel_cycle, time_series_dynamics, response_pattern
- Validate experiments cross-reference, extraction JSON consistency
- Warn on deprecated clusters: block and generic entry keys"
```

---

## Task 2: Update Paperconfig Skill Documentation

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md:264-297`

- [ ] **Step 1: Replace the gene_clusters section in SKILL.md**

Find the gene_clusters documentation (around lines 264-297) and replace with:

````markdown
### `type: gene_clusters`

Co-expression cluster membership tables. Processed by `cluster_adapter.py` (NOT by `omics_adapter`).

**Entry key** must be short, meaningful, and unique within the paper (used in graph node IDs):
- Good: `med4_kmeans_nstarvation`, `mit9313_mfuzz_diel`
- Bad: `cluster_table_1`, `supp_table_clusters`

**Required fields (analysis level):**

| Field | Type | Description |
|---|---|---|
| `name` | str | Human-readable label for the clustering analysis (e.g., "MED4 K-means N-starvation clusters") |
| `filename` | str | Path to cluster membership CSV |
| `organism` | str | Target organism (canonical name) |
| `gene_id_col` | str | CSV column with gene identifiers |
| `cluster_col` | str | CSV column with cluster assignment |
| `cluster_type` | str | `diel_cycle` \| `time_series_dynamics` \| `response_pattern` |

**Optional fields (analysis level):**

| Field | Type | Description |
|---|---|---|
| `score_col` | str | CSV column with fuzzy membership score |
| `cluster_method` | str | Algorithm (e.g., "K-means (K=9)", "Mfuzz") |
| `omics_type` | str | MICROARRAY, RNASEQ, PROTEOMICS, METABOLOMICS |
| `light_condition` | str | Light regime (e.g., "14:10 L:D") |
| `treatment_type` | str[] | Array of canonical treatment types |
| `treatment` | str | Experiment description |
| `experimental_context` | str | Additional context |
| `experiments` | str[] | List of experiment keys from the same paperconfig (links ClusteringAnalysis → Experiment) |

**Per-cluster data** comes from extraction JSON (`cluster_extraction_{entry_key}.json`), not from the paperconfig. The extraction pipeline reads cluster keys from the CSV `cluster_col` and produces: `id`, `name`, `functional_description`, `behavioral_description`, `peak_time_hours`, `period_hours`.

**Data flow:** paperconfig (analysis metadata) + CSV (cluster membership) → extraction pipeline → extraction JSON (per-cluster descriptions) → adapter reads all three.

**Example:**

```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"
  filename: "data/.../med4_kmeans_clusters.csv"
  organism: "Prochlorococcus MED4"
  gene_id_col: "gene_id"
  cluster_col: "cluster"
  cluster_type: "response_pattern"
  cluster_method: "K-means (K=9)"
  omics_type: MICROARRAY
  treatment_type: ["nitrogen_stress"]
  treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
  light_condition: "continuous light"
  experimental_context: "Custom Affymetrix microarray..."
  experiments: [n_starvation_med4]
```

**Notes:**
- Gene IDs go through the same step 4 resolution pipeline as DE tables
- For papers with separate clusters per organism, use separate `type: gene_clusters` entries (one per organism/analysis)
- `treatment_type` must be an array (same enum as experiments)
````

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md
git commit -m "docs: update paperconfig skill for ClusteringAnalysis format

- Entry key naming convention, new required/optional fields
- Remove clusters: block docs, add extraction JSON data flow
- New cluster_type enum values"
```

---

## Task 3: Migrate Tolonen 2006 Paperconfig

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml:238-361`
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_med4.json`
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_mit9313.json`

- [ ] **Step 1: Rename entry keys and update paperconfig**

In `paperconfig.yaml`, replace `cluster_table_med4:` (line 238) and its content with:

```yaml
  med4_kmeans_nstarvation:
    type: gene_clusters
    name: "MED4 K-means N-starvation clusters"
    filename: "data/Prochlorococcus/papers_and_supp/tolonen 2006/med4_kmeans_clusters.csv"
    organism: "Prochlorococcus MED4"
    gene_id_col: "gene_id"
    cluster_col: "cluster"
    cluster_type: "response_pattern"
    cluster_method: "K-means (K=9)"
    omics_type: MICROARRAY
    light_condition: "continuous light"
    treatment_type: ["nitrogen_stress"]
    treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
    experimental_context: "Custom Affymetrix microarray MD4-9313; MED4 cells in Pro99 medium under continuous light (27 µmol photons m-2 s-1, 24°C); sampled at 0, 3, 6, 12, 24, 48h post N-deprivation"
```

Replace `cluster_table_mit9313:` (line 306) and its content with:

```yaml
  mit9313_kmeans_nstarvation:
    type: gene_clusters
    name: "MIT9313 K-means N-starvation clusters"
    filename: "data/Prochlorococcus/papers_and_supp/tolonen 2006/mit9313_kmeans_clusters.csv"
    organism: "Prochlorococcus MIT9313"
    gene_id_col: "gene_id"
    cluster_col: "cluster"
    cluster_type: "response_pattern"
    cluster_method: "K-means (K=7)"
    omics_type: MICROARRAY
    light_condition: "continuous light"
    treatment_type: ["nitrogen_stress"]
    treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
    experimental_context: "Custom Affymetrix microarray MD4-9313; MIT9313 cells in Pro99 medium under continuous light (27 µmol photons m-2 s-1, 24°C); sampled at 0, 3, 6, 12, 24, 48h post N-deprivation"
```

Remove the entire `clusters:` block from both entries (all per-cluster definitions).

- [ ] **Step 2: Rename extraction JSON files**

```bash
mv "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_med4.json" \
   "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_med4_kmeans_nstarvation.json"
mv "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_med4.md" \
   "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_med4_kmeans_nstarvation.md"
mv "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_mit9313.json" \
   "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_mit9313_kmeans_nstarvation.json"
mv "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_mit9313.md" \
   "data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_extraction_mit9313_kmeans_nstarvation.md"
```

- [ ] **Step 3: Update table_key in extraction JSONs**

In `cluster_extraction_med4_kmeans_nstarvation.json`, change `metadata.table_key` from `"cluster_table_med4"` to `"med4_kmeans_nstarvation"`.

In `cluster_extraction_mit9313_kmeans_nstarvation.json`, change `metadata.table_key` from `"cluster_table_mit9313"` to `"mit9313_kmeans_nstarvation"`.

- [ ] **Step 4: Run validator on migrated paperconfig**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

Expected: passes with no errors. May warn about deprecated `clusters` block if removal was incomplete.

- [ ] **Step 5: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/tolonen 2006/"
git commit -m "data: migrate tolonen 2006 paperconfig to ClusteringAnalysis format

- Rename entry keys: cluster_table_med4 → med4_kmeans_nstarvation,
  cluster_table_mit9313 → mit9313_kmeans_nstarvation
- Add name, cluster_type at analysis level
- Remove clusters: block (per-cluster data in extraction JSON)
- Rename and update extraction JSON files to match"
```

---

## Task 4: Update Extraction Pipeline

**Files:**
- Modify: `multiomics_kg/extraction/cluster/prompts.py:112-159`
- Modify: `multiomics_kg/extraction/cluster/pipeline.py:155-159`

- [ ] **Step 1: Update SYNTHESIS_PROMPT in prompts.py**

In `prompts.py`, find the `SYNTHESIS_PROMPT` (line 112). In the "Write three outputs per cluster:" section (around line 120), replace with four outputs:

```python
SYNTHESIS_PROMPT = """\
You are writing descriptions for gene expression clusters from a scientific paper.

For each cluster below, you are given extracted data from multiple sources \
(table=structured data files, visual=PDF figures, semantic=paper text). \
Each field is a list of extractions tagged with source and confidence level \
(very_high > high > medium > low).

Write these outputs per cluster:

1. **id** — Short snake_case identifier prefixed with organism. \
Format: {{organism}}_{{direction}}_{{theme}}. \
Examples: med4_up_n_transport, mit9313_down_translation, med4_diel_dawn_psi. \
Must be unique across all clusters in this paper.

2. **name** — Short human-readable cluster label. \
Format: "{{Organism}} cluster {{key}} ({{direction}}, {{theme}})". \
Examples: "MED4 cluster 1 (up, N transport)", "MIT9313 cluster 3 (down, translation)". \
Keep under 60 characters.

3. **functional_description** — What types of genes are in this cluster. \
Include enrichment category, p-value, and specific gene names. Prefer \
higher-confidence sources. Be precise.

4. **behavioral_description** — The temporal/response pattern. Include \
direction (up/down), timing, and magnitude if available.

5. **peak_time_hours** — Peak expression time in hours (float). \
Only for diel/periodic clusters. null for non-periodic clusters.

6. **period_hours** — Oscillation period in hours (float). \
Only for periodic clusters. null for non-periodic clusters.

RULES:
- Never leave fields empty. Use these sentinel values:
  - "not described in paper" — paper doesn't discuss this aspect
  - "insufficient data" — sources found but too ambiguous
  - "conflicting sources: [source A says X, source B says Y]" — paths disagree
- Partial descriptions are fine; incorrect descriptions are not.
- Prefer higher-confidence sources when they agree on meaning.
- If uncertain, use the sentinel value — do NOT guess.
- IDs must be unique within this paper (organism prefix helps).
- peak_time_hours and period_hours must be null (not "null") for non-periodic clusters.

Paper: {paper_name}
Organism: {organism}
Treatment: {treatment}
Cluster method: {cluster_method}

{cluster_blocks}

Return valid JSON using the original cluster keys:
{{
  "<cluster_key>": {{
    "id": "organism_direction_theme",
    "name": "Organism cluster N (direction, theme)",
    "functional_description": "...",
    "behavioral_description": "...",
    "peak_time_hours": null,
    "period_hours": null
  }}
}}
"""
```

- [ ] **Step 2: Update pipeline output filename**

In `pipeline.py` line 159, change:

```python
        out_path = paper_dir / f"cluster_extraction_{organism_short}.json"
```

to:

```python
        out_path = paper_dir / f"cluster_extraction_{tkey}.json"
```

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/cluster/prompts.py multiomics_kg/extraction/cluster/pipeline.py
git commit -m "feat: add name/peak_time_hours/period_hours to extraction synthesis

- Synthesis prompt now produces 6 fields per cluster
- Pipeline output filename uses entry_key instead of organism_short"
```

---

## Task 5: Update Schema Config

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add clustering analysis node type**

After the `gene cluster` node definition (around line 86), add:

```yaml
  clustering analysis:
    is_a: information content entity
    represented_as: node
    preferred_id: clustering_analysis_id
    label_in_input: clustering_analysis
    properties:
      name: str
      organism_name: str
      cluster_method: str
      cluster_count: int
      total_gene_count: int
      omics_type: str
      treatment_type: str[]
      treatment: str
      light_condition: str
      cluster_type: str
      experimental_context: str
```

- [ ] **Step 2: Update gene cluster node properties**

In the existing `gene cluster` node (lines 64-86), update properties to:

```yaml
  gene cluster:
    is_a: GroupingClass
    represented_as: node
    preferred_id: cluster_id
    label_in_input: gene_cluster
    properties:
      id: str
      name: str
      organism_name: str
      cluster_method: str
      cluster_type: str
      treatment_type: str[]
      treatment: str
      omics_type: str
      light_condition: str
      member_count: int
      functional_description: str
      behavioral_description: str
      peak_time_hours: float
      period_hours: float
      experimental_context: str
```

Changes: added `id` property, removed `source_paper`.

- [ ] **Step 3: Add new edge types**

Add these edge types:

```yaml
  publication has clustering analysis:
    is_a: related to at instance level
    represented_as: edge
    label_in_input: publication_has_clustering_analysis
    source: publication
    target: clustering analysis

  clustering analysis has gene cluster:
    is_a: related to at instance level
    represented_as: edge
    label_in_input: clustering_analysis_has_gene_cluster
    source: clustering analysis
    target: gene cluster

  clusteringanalysis belongs to organism:
    is_a: related to at instance level
    represented_as: edge
    label_in_input: clusteringanalysis_belongs_to_organism
    source: clustering analysis
    target: organism taxon

  experiment has clustering analysis:
    is_a: related to at instance level
    represented_as: edge
    label_in_input: experiment_has_clustering_analysis
    source: experiment
    target: clustering analysis
```

- [ ] **Step 4: Remove old edge types**

Remove these edge definitions:
- `publication has gene cluster` (lines 105-112)
- `genecluster belongs to organism` (lines 126-133)

- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add ClusteringAnalysis node, restructure cluster edges

- New node: clustering analysis with 11 properties
- New edges: publication_has_clustering_analysis,
  clustering_analysis_has_gene_cluster,
  clusteringanalysis_belongs_to_organism,
  experiment_has_clustering_analysis
- Removed edges: publication_has_gene_cluster,
  genecluster_belongs_to_organism
- GeneCluster: added id property, removed source_paper"
```

---

## Task 6: Restructure Cluster Adapter — Tests First

**Files:**
- Modify: `tests/test_cluster_adapter.py`

- [ ] **Step 1: Write test for ClusteringAnalysis node emission**

Add a new test that verifies the adapter yields ClusteringAnalysis nodes:

```python
def test_get_nodes_emits_clustering_analysis_node(tmp_path, monkeypatch):
    """ClusterAdapter should yield one ClusteringAnalysis node per gene_clusters entry."""
    # Create a minimal CSV
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,1\nPMM0003,2\n")

    # Create paperconfig
    paperconfig = {
        "publication": {
            "papername": "Test 2006",
        },
        "supplementary_materials": {
            "med4_kmeans_test": {
                "type": "gene_clusters",
                "name": "MED4 K-means test clusters",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
                "cluster_method": "K-means (K=2)",
                "omics_type": "MICROARRAY",
                "treatment_type": ["nitrogen_stress"],
                "treatment": "N-starvation",
                "light_condition": "continuous light",
                "experimental_context": "test context",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    import yaml
    config_path.write_text(yaml.dump(paperconfig))

    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter
    adapter = ClusterAdapter(str(config_path))

    nodes = list(adapter.get_nodes())

    # Should have 1 ClusteringAnalysis + 2 GeneCluster = 3 nodes
    analysis_nodes = [(nid, label, props) for nid, label, props in nodes if label == "clustering_analysis"]
    cluster_nodes = [(nid, label, props) for nid, label, props in nodes if label == "gene_cluster"]

    assert len(analysis_nodes) == 1
    assert len(cluster_nodes) == 2

    # Verify analysis node
    a_id, a_label, a_props = analysis_nodes[0]
    assert "med4_kmeans_test" in a_id
    assert a_props["name"] == "MED4 K-means test clusters"
    assert a_props["organism_name"] == "Prochlorococcus MED4"
    assert a_props["cluster_count"] == 2
    assert a_props["total_gene_count"] == 3
    assert a_props["cluster_type"] == "response_pattern"
    assert a_props["cluster_method"] == "K-means (K=2)"
```

- [ ] **Step 2: Write test for GeneCluster node ID format**

```python
def test_gene_cluster_node_id_includes_analysis_key(tmp_path):
    """GeneCluster node ID should be cluster:{doi}:{analysis_key}:{csv_key}."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,A\nPMM0002,B\n")

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "med4_test_analysis": {
                "type": "gene_clusters",
                "name": "Test analysis",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    import yaml
    config_path.write_text(yaml.dump(paperconfig))

    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter
    adapter = ClusterAdapter(str(config_path))

    nodes = list(adapter.get_nodes())
    cluster_nodes = [(nid, label, props) for nid, label, props in nodes if label == "gene_cluster"]

    cluster_ids = sorted([nid for nid, _, _ in cluster_nodes])
    # IDs should contain the analysis key and CSV cluster key
    assert any("med4_test_analysis" in cid and ":A" in cid for cid in cluster_ids)
    assert any("med4_test_analysis" in cid and ":B" in cid for cid in cluster_ids)
```

- [ ] **Step 3: Write test for extraction JSON reading**

```python
def test_get_nodes_reads_extraction_json(tmp_path):
    """Adapter should populate descriptions from extraction JSON."""
    import json
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    extraction = {
        "metadata": {"table_key": "test_analysis"},
        "stage2_results": {
            "1": {
                "id": "med4_up_transport",
                "name": "MED4 cluster 1 (up, transport)",
                "functional_description": "Transport genes",
                "behavioral_description": "Rapid upregulation",
                "peak_time_hours": None,
                "period_hours": None,
            },
            "2": {
                "id": "med4_down_photosynthesis",
                "name": "MED4 cluster 2 (down, photosynthesis)",
                "functional_description": "PSI genes",
                "behavioral_description": "Downregulated at 6h",
                "peak_time_hours": None,
                "period_hours": None,
            },
        },
        "stage3_validation": {
            "1": {"verdict": "pass"},
            "2": {"verdict": "pass"},
        },
    }
    json_path = tmp_path / "cluster_extraction_test_analysis.json"
    json_path.write_text(json.dumps(extraction))

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters",
                "name": "Test analysis",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    import yaml
    config_path.write_text(yaml.dump(paperconfig))

    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter
    adapter = ClusterAdapter(str(config_path))

    nodes = list(adapter.get_nodes())
    cluster_nodes = {props.get("id"): props
                     for _, label, props in nodes if label == "gene_cluster"}

    assert cluster_nodes["med4_up_transport"]["functional_description"] == "Transport genes"
    assert cluster_nodes["med4_up_transport"]["name"] == "MED4 cluster 1 (up, transport)"
    assert cluster_nodes["med4_down_photosynthesis"]["behavioral_description"] == "Downregulated at 6h"
```

- [ ] **Step 4: Write test for extraction JSON with failed validation**

```python
def test_get_nodes_skips_failed_extraction(tmp_path):
    """Adapter should skip clusters with verdict=fail in extraction JSON."""
    import json
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    extraction = {
        "metadata": {"table_key": "test_analysis"},
        "stage2_results": {
            "1": {"id": "good", "name": "Good", "functional_description": "OK",
                  "behavioral_description": "OK", "peak_time_hours": None, "period_hours": None},
            "2": {"id": "bad", "name": "Bad", "functional_description": "Wrong",
                  "behavioral_description": "Wrong", "peak_time_hours": None, "period_hours": None},
        },
        "stage3_validation": {
            "1": {"verdict": "pass"},
            "2": {"verdict": "fail"},
        },
    }
    json_path = tmp_path / "cluster_extraction_test_analysis.json"
    json_path.write_text(json.dumps(extraction))

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters",
                "name": "Test",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    import yaml
    config_path.write_text(yaml.dump(paperconfig))

    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter
    adapter = ClusterAdapter(str(config_path))

    nodes = list(adapter.get_nodes())
    cluster_nodes = {props.get("id", ""): props
                     for _, label, props in nodes if label == "gene_cluster"}

    # Cluster 1 should have descriptions from extraction
    assert cluster_nodes["good"]["functional_description"] == "OK"
    # Cluster 2 should have empty descriptions (verdict=fail)
    failed = [p for p in cluster_nodes.values() if p.get("id") in ("bad", "")]
    assert len(failed) == 1
    assert failed[0]["functional_description"] == ""
```

- [ ] **Step 5: Write test for new edge types**

```python
def test_get_edges_new_structure(tmp_path, monkeypatch):
    """Adapter should emit 5 edge types including ClusteringAnalysis edges."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    paperconfig = {
        "publication": {
            "papername": "Test 2006",
            "experiments": {
                "exp1": {"name": "test experiment"},
            }
        },
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters",
                "name": "Test",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
                "experiments": ["exp1"],
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    import yaml
    config_path.write_text(yaml.dump(paperconfig))

    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter
    adapter = ClusterAdapter(str(config_path))
    # Need organism_lookup for belongs_to_organism
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}

    edges = list(adapter.get_edges())
    edge_types = set(e[3] for e in edges)

    assert "gene_in_gene_cluster" in edge_types
    assert "publication_has_clustering_analysis" in edge_types
    assert "clustering_analysis_has_gene_cluster" in edge_types
    assert "clusteringanalysis_belongs_to_organism" in edge_types
    assert "experiment_has_clustering_analysis" in edge_types

    # Should NOT have old edge types
    assert "publication_has_gene_cluster" not in edge_types
    assert "genecluster_belongs_to_organism" not in edge_types

    # Count: 2 membership + 1 pub→analysis + 2 analysis→cluster + 1 analysis→org + 1 exp→analysis = 7
    assert len(edges) == 7
```

- [ ] **Step 6: Run tests to verify they fail**

```bash
pytest tests/test_cluster_adapter.py -v -k "clustering_analysis or new_structure or extraction_json or failed_extraction or analysis_key"
```

Expected: all new tests FAIL (adapter not yet updated).

- [ ] **Step 7: Commit test file**

```bash
git add tests/test_cluster_adapter.py
git commit -m "test: add tests for ClusteringAnalysis node and extraction JSON read"
```

---

## Task 7: Restructure Cluster Adapter — Implementation

**Files:**
- Modify: `multiomics_kg/adapters/cluster_adapter.py`

- [ ] **Step 1: Add extraction JSON loading helper**

Add near the top of the file (after `_clean_str`):

```python
def _load_extraction_json(paperconfig_dir: Path, entry_key: str) -> dict:
    """Load extraction JSON for a gene_clusters entry. Returns empty dict if missing."""
    json_path = paperconfig_dir / f"cluster_extraction_{entry_key}.json"
    if not json_path.exists():
        return {}
    try:
        with open(json_path) as f:
            data = json.load(f)
        return data
    except Exception:
        logger.warning("Failed to load extraction JSON: %s", json_path)
        return {}


def _get_extraction_cluster_data(extraction: dict, cluster_key: str) -> dict:
    """Get per-cluster data from extraction JSON, respecting validation verdict."""
    stage2 = extraction.get("stage2_results", {})
    stage3 = extraction.get("stage3_validation", {})

    cluster_data = stage2.get(str(cluster_key), {})
    verdict = stage3.get(str(cluster_key), {}).get("verdict", "")

    if verdict != "pass":
        return {}  # Skip failed/missing validation
    return cluster_data
```

Add `import json` at the top if not already present.

- [ ] **Step 2: Update `_make_cluster_id` for new format**

Replace the existing `_make_cluster_id` function (lines 34-47) with:

```python
def _make_cluster_id(doi: str, paper_name: str, entry_key: str, cluster_key: str) -> str:
    """Build cluster node ID: cluster:{doi_short}:{entry_key}:{cluster_key}."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"cluster:{doi_short}:{entry_key}:{cluster_key}"


def _make_analysis_id(doi: str, paper_name: str, entry_key: str) -> str:
    """Build clustering analysis node ID."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"clustering_analysis:{doi_short}:{entry_key}"
```

- [ ] **Step 3: Rewrite `get_nodes()` to yield both node types**

Replace the `get_nodes()` method of `ClusterAdapter` with:

```python
def get_nodes(self):
    """Yield ClusteringAnalysis nodes and GeneCluster nodes."""
    import pandas as pd

    paper_name = self.config.get("publication", {}).get("papername", "")
    paperconfig_dir = Path(self.config_file).parent

    for entry_key, table in iter_cluster_tables(self.config):
        organism = table.get("organism", "")
        csv_path = self._resolve_csv_path(table)
        if csv_path is None:
            continue

        sep = table.get("separator", ",")
        skip_rows = table.get("skip_rows", 0)
        try:
            df = pd.read_csv(csv_path, sep=sep, skiprows=skip_rows)
        except Exception as e:
            logger.warning("Cannot read CSV %s: %s", csv_path, e)
            continue

        cluster_col = table.get("cluster_col", "cluster")
        cluster_keys = [str(k) for k in sorted(
            df[cluster_col].unique(),
            key=lambda x: (not str(x).isdigit(), str(x))
        )]

        # Load extraction JSON
        extraction = _load_extraction_json(paperconfig_dir, entry_key)

        # Analysis-level properties (from paperconfig)
        analysis_props = {
            "name": _clean_str(table.get("name", "")),
            "organism_name": _clean_str(organism),
            "cluster_method": _clean_str(table.get("cluster_method", "")),
            "cluster_type": _clean_str(table.get("cluster_type", "")),
            "omics_type": _clean_str(table.get("omics_type", "")),
            "treatment_type": table.get("treatment_type", []),
            "treatment": _clean_str(table.get("treatment", "")),
            "light_condition": _clean_str(table.get("light_condition", "")),
            "experimental_context": _clean_str(table.get("experimental_context", "")),
            "cluster_count": len(cluster_keys),
            "total_gene_count": 0,  # updated below
        }

        # Yield ClusteringAnalysis node
        analysis_id = _make_analysis_id(self.doi, paper_name, entry_key)

        total_genes = 0
        cluster_nodes = []
        for ck in cluster_keys:
            members = df[df[cluster_col].astype(str) == ck]
            member_count = len(members)
            total_genes += member_count

            ext_data = _get_extraction_cluster_data(extraction, ck)

            cluster_id = _make_cluster_id(self.doi, paper_name, entry_key, ck)
            props = {
                "id": _clean_str(ext_data.get("id", "")),
                "name": _clean_str(ext_data.get("name", "")),
                "member_count": member_count,
                "functional_description": _clean_str(ext_data.get("functional_description", "")),
                "behavioral_description": _clean_str(ext_data.get("behavioral_description", "")),
                "peak_time_hours": ext_data.get("peak_time_hours"),
                "period_hours": ext_data.get("period_hours"),
                # Denormalized from analysis
                "organism_name": analysis_props["organism_name"],
                "cluster_method": analysis_props["cluster_method"],
                "cluster_type": analysis_props["cluster_type"],
                "omics_type": analysis_props["omics_type"],
                "treatment_type": analysis_props["treatment_type"],
                "treatment": analysis_props["treatment"],
                "light_condition": analysis_props["light_condition"],
                "experimental_context": analysis_props["experimental_context"],
            }
            cluster_nodes.append((cluster_id, "gene_cluster", props))

        analysis_props["total_gene_count"] = total_genes
        yield analysis_id, "clustering_analysis", analysis_props

        for node in cluster_nodes:
            yield node
```

- [ ] **Step 4: Rewrite `get_edges()` for five edge types**

Replace the `get_edges()` method of `ClusterAdapter` with:

```python
def get_edges(self):
    """Yield edges for ClusteringAnalysis → GeneCluster → Gene structure."""
    import pandas as pd

    paper_name = self.config.get("publication", {}).get("papername", "")
    doi = self.doi or ""

    for entry_key, table in iter_cluster_tables(self.config):
        organism = table.get("organism", "")
        csv_path = self._resolve_csv_path(table)
        if csv_path is None:
            continue

        sep = table.get("separator", ",")
        skip_rows = table.get("skip_rows", 0)
        try:
            df = pd.read_csv(csv_path, sep=sep, skiprows=skip_rows)
        except Exception:
            continue

        cluster_col = table.get("cluster_col", "cluster")
        gene_col = "locus_tag" if "locus_tag" in df.columns else table.get("gene_id_col", "gene_id")
        score_col = table.get("score_col")
        cluster_keys = [str(k) for k in sorted(
            df[cluster_col].unique(),
            key=lambda x: (not str(x).isdigit(), str(x))
        )]

        analysis_id = _make_analysis_id(doi, paper_name, entry_key)

        # 1. Publication_has_clustering_analysis
        if doi:
            yield (
                f"{analysis_id}__pub",
                doi,
                analysis_id,
                "publication_has_clustering_analysis",
                {},
            )

        # 2. Clusteringanalysis_belongs_to_organism
        org_id = self._organism_lookup.get(organism)
        if org_id:
            yield (
                f"{analysis_id}__org",
                analysis_id,
                org_id,
                "clusteringanalysis_belongs_to_organism",
                {},
            )

        # 3. Experiment_has_clustering_analysis
        experiments_ref = table.get("experiments", [])
        for exp_key in experiments_ref:
            exp_id = f"{doi}_{exp_key}" if doi else exp_key
            yield (
                f"{analysis_id}__exp__{exp_key}",
                exp_id,
                analysis_id,
                "experiment_has_clustering_analysis",
                {},
            )

        # Per-cluster edges
        for ck in cluster_keys:
            cluster_id = _make_cluster_id(doi, paper_name, entry_key, ck)

            # 4. Clustering_analysis_has_gene_cluster
            yield (
                f"{analysis_id}__{ck}",
                analysis_id,
                cluster_id,
                "clustering_analysis_has_gene_cluster",
                {},
            )

            # 5. Gene_in_gene_cluster
            members = df[df[cluster_col].astype(str) == ck]
            if self.test_mode:
                members = members.head(100)

            for _, row in members.iterrows():
                gene_locus = str(row.get(gene_col, ""))
                if not gene_locus or gene_locus == "nan":
                    continue

                if not gene_locus.startswith("ncbigene:"):
                    gene_locus = f"ncbigene:{gene_locus}"

                edge_props = {}
                if score_col and score_col in df.columns:
                    score_val = row.get(score_col)
                    if pd.notna(score_val):
                        edge_props["membership_score"] = float(score_val)

                yield (
                    f"{cluster_id}__{gene_locus}",
                    cluster_id,
                    gene_locus,
                    "gene_in_gene_cluster",
                    edge_props,
                )
```

- [ ] **Step 5: Run the new tests**

```bash
pytest tests/test_cluster_adapter.py -v -k "clustering_analysis or new_structure or extraction_json or failed_extraction or analysis_key"
```

Expected: all new tests PASS.

- [ ] **Step 6: Run all cluster adapter tests**

```bash
pytest tests/test_cluster_adapter.py -v
```

Expected: old tests may fail due to changed API. Fix any that need updating (old tests reference `publication_has_gene_cluster`, `genecluster_belongs_to_organism`, etc. — update them to new edge types and node structure).

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/adapters/cluster_adapter.py tests/test_cluster_adapter.py
git commit -m "feat: restructure cluster adapter for ClusteringAnalysis nodes

- Yields ClusteringAnalysis + GeneCluster nodes (two types)
- Five edge types: pub→analysis, analysis→cluster, analysis→org,
  exp→analysis (optional), cluster→gene
- Reads extraction JSON for per-cluster descriptions
- Cluster ID format: cluster:{doi}:{analysis_key}:{csv_key}
- Skips extraction results with verdict != pass"
```

---

## Task 8: Update Post-Import Scripts

**Files:**
- Modify: `scripts/post-import.sh:76-91`
- Modify: `scripts/post-import.cypher:61-78`

- [ ] **Step 1: Update post-import.sh**

Replace the GeneCluster index section (lines 76-91) with:

```bash
# ── ClusteringAnalysis indexes ──
CALL {
  CREATE INDEX clustering_analysis_organism_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.organism_name);
  CREATE INDEX clustering_analysis_method_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.cluster_method);
  CREATE INDEX clustering_analysis_type_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.cluster_type);

  CREATE FULLTEXT INDEX clusteringAnalysisFullText IF NOT EXISTS
    FOR (ca:ClusteringAnalysis) ON EACH [ca.name, ca.treatment, ca.experimental_context];
} IN TRANSACTIONS;

# ── GeneCluster indexes (updated — removed organism idx, kept fulltext) ──
CALL {
  CREATE INDEX gene_cluster_treatment_type_idx IF NOT EXISTS FOR (gc:GeneCluster) ON (gc.treatment_type);
  CREATE INDEX gene_cluster_type_idx IF NOT EXISTS FOR (gc:GeneCluster) ON (gc.cluster_type);

  CREATE FULLTEXT INDEX geneClusterFullText IF NOT EXISTS
    FOR (gc:GeneCluster) ON EACH [gc.name, gc.functional_description, gc.behavioral_description, gc.experimental_context];
} IN TRANSACTIONS;

# ── GeneCluster member_count verification ──
CALL {
  MATCH (gc:GeneCluster)
  OPTIONAL MATCH (gc)-[r:Gene_in_gene_cluster]->()
  WITH gc, count(r) AS actual_count
  SET gc.member_count = actual_count;
} IN TRANSACTIONS;
```

- [ ] **Step 2: Update post-import.cypher with identical changes**

Apply the same index changes to `scripts/post-import.cypher` (lines 61-78). Keep both files in sync.

- [ ] **Step 3: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "post-import: add ClusteringAnalysis indexes, remove gene_cluster_organism_idx"
```

---

## Task 9: Run Full Unit Tests

- [ ] **Step 1: Run all unit tests**

```bash
pytest -m "not slow and not kg" -v
```

Expected: all pass. Fix any failures.

- [ ] **Step 2: Commit any fixes**

If fixes were needed:
```bash
git add -u
git commit -m "fix: address test failures from ClusteringAnalysis restructure"
```

---

## Task 10: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update key graph facts**

Update the relevant sections of CLAUDE.md:
- Add ClusteringAnalysis to the node list
- Add new edge types to the relationship list
- Remove `Publication_has_gene_cluster` and `Genecluster_belongs_to_organism`
- Update GeneCluster property list (add `id`, remove `source_paper`, note `cluster_type` is denormalized)
- Add ClusteringAnalysis to post-import indexes section

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for ClusteringAnalysis node structure"
```
