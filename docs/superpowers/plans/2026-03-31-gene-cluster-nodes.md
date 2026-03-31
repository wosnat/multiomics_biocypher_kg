# GeneCluster Nodes Implementation Plan (Phase 1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add GeneCluster nodes and membership edges to the KG, integrating published co-expression clusters from Tolonen 2006 (N-stress, MED4 + MIT9313) as the first data source.

**Architecture:** New `cluster_adapter.py` reads `type: gene_clusters` entries from paperconfigs. Follows the MultiOMICSAdapter pattern (multi-wrapper iterates paperconfig_files.txt, single adapter processes one paperconfig). Gene IDs go through the existing step 4 resolution pipeline. Post-import creates indexes and full-text search.

**Tech Stack:** Python, BioCypher, PyYAML, pandas, Neo4j Cypher (post-import)

**Spec:** `docs/superpowers/specs/2026-03-31-gene-cluster-nodes-design.md`

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `multiomics_kg/adapters/cluster_adapter.py` | ClusterAdapter + MultiClusterAdapter |
| Create | `tests/test_cluster_adapter.py` | Unit tests for cluster adapter |
| Modify | `config/schema_config.yaml` | Add gene_cluster node + 3 edge types |
| Modify | `create_knowledge_graph.py` | Wire up MultiClusterAdapter |
| Modify | `scripts/post-import.sh` | Indexes + member_count verification |
| Modify | `scripts/post-import.cypher` | Same (reference copy) |
| Modify | `.claude/skills/paperconfig/validate_paperconfig.py` | Validate `type: gene_clusters` entries |
| Modify | `.claude/skills/paperconfig/SKILL.md` | Document gene_clusters format |
| Modify | `multiomics_kg/download/resolve_paper_ids.py` | Handle `gene_clusters` table type |
| Modify | `multiomics_kg/utils/paperconfig_utils.py` | Add `iter_cluster_tables()` helper |
| Create | `data/Prochlorococcus/papers_and_supp/tolonen 2006/cluster_data/` | Tolonen cluster CSV + paperconfig updates |

---

### Task 1: Schema — Add GeneCluster node and edges to schema_config.yaml

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add gene cluster node type**

Add after the `experiment` node block in `config/schema_config.yaml`:

```yaml
gene cluster:
  is_a: gene cluster
  represented_as: node
  preferred_id: cluster_id
  label_in_input: gene_cluster
  properties:
    name: str
    source_paper: str
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

- [ ] **Step 2: Add publication_has_gene_cluster edge type**

```yaml
publication has gene cluster:
  is_a: Association
  represented_as: edge
  label_as_edge: publication_has_gene_cluster
  label_in_input: publication_has_gene_cluster
  source: publication
  target: gene cluster
```

- [ ] **Step 3: Add gene_in_gene_cluster edge type**

```yaml
gene in gene cluster:
  is_a: Association
  represented_as: edge
  label_as_edge: gene_in_gene_cluster
  label_in_input: gene_in_gene_cluster
  source: gene cluster
  target: gene
  properties:
    membership_score: float
    p_value: float
```

- [ ] **Step 4: Add genecluster_belongs_to_organism edge type**

```yaml
gene cluster belongs to organism:
  is_a: Association
  represented_as: edge
  label_as_edge: genecluster_belongs_to_organism
  label_in_input: genecluster_belongs_to_organism
  source: gene cluster
  target: organism taxon
```

- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add GeneCluster node and 3 edge types"
```

---

### Task 2: Paperconfig utilities — Add iter_cluster_tables helper

**Files:**
- Modify: `multiomics_kg/utils/paperconfig_utils.py`
- Test: `tests/test_paperconfig_utils.py` (if exists, otherwise inline verification)

- [ ] **Step 1: Write test for iter_cluster_tables**

Create or append to `tests/test_paperconfig_utils.py`:

```python
def test_iter_cluster_tables():
    """iter_cluster_tables yields (key, config) for gene_clusters entries only."""
    from multiomics_kg.utils.paperconfig_utils import iter_cluster_tables

    config = {
        "publication": {
            "supplementary_materials": {
                "supp_table_1": {"type": "csv", "filename": "data.csv"},
                "cluster_table_1": {
                    "type": "gene_clusters",
                    "filename": "clusters.csv",
                    "organism": "Prochlorococcus MED4",
                    "gene_id_col": "ORF",
                    "cluster_col": "cluster",
                    "clusters": {
                        "c1": {"name": "Cluster 1", "cluster_type": "stress_response"},
                    },
                },
                "id_trans": {"type": "id_translation", "filename": "ids.csv"},
            }
        }
    }

    results = list(iter_cluster_tables(config))
    assert len(results) == 1
    key, table = results[0]
    assert key == "cluster_table_1"
    assert table["type"] == "gene_clusters"
    assert table["organism"] == "Prochlorococcus MED4"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_paperconfig_utils.py::test_iter_cluster_tables -v
```
Expected: FAIL — `iter_cluster_tables` not defined.

- [ ] **Step 3: Implement iter_cluster_tables**

Add to `multiomics_kg/utils/paperconfig_utils.py`, near `iter_csv_tables`:

```python
def iter_cluster_tables(config: dict):
    """Yield (table_key, table_config) for gene_clusters-type supplementary tables."""
    supp = get_supplementary_materials(config)
    for key, table in supp.items():
        if not isinstance(table, dict):
            continue
        if table.get("type") == "gene_clusters":
            yield key, table
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/test_paperconfig_utils.py::test_iter_cluster_tables -v
```
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/paperconfig_utils.py tests/test_paperconfig_utils.py
git commit -m "feat: add iter_cluster_tables to paperconfig_utils"
```

---

### Task 3: Cluster adapter — ClusterAdapter and MultiClusterAdapter

**Files:**
- Create: `multiomics_kg/adapters/cluster_adapter.py`
- Create: `tests/test_cluster_adapter.py`

- [ ] **Step 1: Write failing test for ClusterAdapter.get_nodes**

Create `tests/test_cluster_adapter.py`:

```python
"""Tests for ClusterAdapter — GeneCluster nodes and membership edges."""
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest
import yaml


@pytest.fixture
def temp_cluster_dir():
    """Create temp dir with cluster CSV and paperconfig."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create cluster CSV: gene_id, cluster, membership_score
        df = pd.DataFrame({
            "ORF": ["PMM0001", "PMM0002", "PMM0003", "PMM0004", "PMM0005"],
            "cluster": [1, 1, 2, 2, 2],
            "membership": [0.95, 0.87, 0.92, 0.78, 0.65],
        })
        csv_path = os.path.join(tmpdir, "clusters.csv")
        df.to_csv(csv_path, index=False)

        # Also create a resolved version (what step 4 would produce)
        df_resolved = df.copy()
        df_resolved["locus_tag"] = df_resolved["ORF"]  # already locus tags
        df_resolved["resolution_method"] = "tier1:locus_tag"
        resolved_path = os.path.join(tmpdir, "clusters_resolved.csv")
        df_resolved.to_csv(resolved_path, index=False)

        # Create paperconfig
        config = {
            "publication": {
                "papername": "Tolonen 2006",
                "papermainpdf": os.path.join(tmpdir, "paper.pdf"),
                "experiments": {},
                "supplementary_materials": {
                    "cluster_table_med4": {
                        "type": "gene_clusters",
                        "filename": csv_path,
                        "organism": "Prochlorococcus MED4",
                        "gene_id_col": "ORF",
                        "cluster_col": "cluster",
                        "score_col": "membership",
                        "omics_type": "MICROARRAY",
                        "light_condition": "continuous light",
                        "treatment_type": ["nitrogen_stress"],
                        "treatment": "N-starvation time course",
                        "experimental_context": "0, 3, 6, 12, 24, 48h N-deprivation",
                        "clusters": {
                            "1": {
                                "name": "Cluster 1",
                                "cluster_type": "stress_response",
                                "functional_description": "Photosynthesis genes",
                                "behavioral_description": "Early upregulation, sustained",
                            },
                            "2": {
                                "name": "Cluster 2",
                                "cluster_type": "stress_response",
                                "functional_description": "Ribosomal proteins",
                                "behavioral_description": "Gradual downregulation",
                            },
                        },
                    }
                },
            }
        }
        config_path = os.path.join(tmpdir, "paperconfig.yaml")
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        # Create dummy PDF (adapter needs it to exist for publication extraction)
        Path(os.path.join(tmpdir, "paper.pdf")).touch()

        yield tmpdir, config_path, config


def test_cluster_adapter_get_nodes(temp_cluster_dir):
    """ClusterAdapter emits GeneCluster nodes with correct properties."""
    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter

    tmpdir, config_path, config = temp_cluster_dir
    adapter = ClusterAdapter(config_file=config_path)

    nodes = list(adapter.get_nodes())
    # Should have 2 GeneCluster nodes (cluster 1 and 2)
    cluster_nodes = [n for n in nodes if n[1] == "gene_cluster"]
    assert len(cluster_nodes) == 2

    # Check first cluster properties
    node_ids = {n[0] for n in cluster_nodes}
    # ID format: cluster:{doi_short}:{cluster_id} — but no DOI in test, so fallback
    first = cluster_nodes[0]
    props = first[2]
    assert props["name"] in ("Cluster 1", "Cluster 2")
    assert props["organism_name"] == "Prochlorococcus MED4"
    assert props["cluster_method"] == ""  # not specified in this config
    assert props["treatment_type"] == ["nitrogen_stress"]
    assert props["omics_type"] == "MICROARRAY"
    assert props["light_condition"] == "continuous light"
    assert isinstance(props["member_count"], int)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_cluster_adapter.py::test_cluster_adapter_get_nodes -v
```
Expected: FAIL — `cluster_adapter` module not found.

- [ ] **Step 3: Implement ClusterAdapter**

Create `multiomics_kg/adapters/cluster_adapter.py`:

```python
"""
Gene Cluster Adapter

Reads gene_clusters entries from paperconfig.yaml files and emits
GeneCluster nodes with Gene_in_gene_cluster membership edges.
"""
import logging
from pathlib import Path
from typing import Union

import pandas as pd

from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig,
    load_all_paperconfigs,
    get_paper_name,
    get_supplementary_materials,
    iter_cluster_tables,
)

logger = logging.getLogger(__name__)


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV output."""
    if not isinstance(value, str):
        return str(value) if value is not None else ""
    return value.replace("'", "^").replace("|", ",")


def _make_cluster_id(doi: str, paper_name: str, table_key: str, cluster_key: str) -> str:
    """Build cluster node ID: cluster:{doi_short}:{cluster_key}.

    Falls back to paper_name if DOI is not available.
    """
    if doi:
        # Shorten DOI: "10.1038/msb4100087" → "msb4100087"
        doi_short = doi.rsplit("/", 1)[-1] if "/" in doi else doi
    else:
        doi_short = paper_name.lower().replace(" ", "_")
    return f"cluster:{doi_short}:{cluster_key}"


def _resolve_csv_path(csv_path: str) -> Path:
    """Find the resolved CSV if it exists, else return original path."""
    p = Path(csv_path)
    resolved = p.parent / f"{p.stem}_resolved{p.suffix}"
    if resolved.exists():
        return resolved
    return p


class ClusterAdapter:
    """Adapter for one paperconfig's gene_clusters entries."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.doi = self._extract_doi()
        self._cluster_tables = list(iter_cluster_tables(self.config))

    def _extract_doi(self) -> str:
        """Try to get DOI from publication block."""
        pub = self.config.get("publication", {})
        return pub.get("doi", "")

    def get_nodes(self) -> list[tuple]:
        """Emit GeneCluster nodes."""
        nodes = []
        for table_key, table in self._cluster_tables:
            csv_path = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                logger.warning(f"Cluster CSV not found: {csv_path}")
                continue

            df = pd.read_csv(csv_path)
            cluster_col = table["cluster_col"]
            gene_id_col = "locus_tag" if "locus_tag" in df.columns else table["gene_id_col"]

            clusters_meta = table.get("clusters", {})

            for cluster_key, cluster_info in clusters_meta.items():
                cluster_id = _make_cluster_id(self.doi, self.paper_name, table_key, cluster_key)

                # Count members for this cluster
                mask = df[cluster_col].astype(str) == str(cluster_key)
                member_count = int(mask.sum())

                props = {
                    "name": _clean_str(cluster_info.get("name", f"Cluster {cluster_key}")),
                    "source_paper": _clean_str(self.paper_name),
                    "organism_name": _clean_str(table.get("organism", "")),
                    "cluster_method": _clean_str(table.get("cluster_method", "")),
                    "cluster_type": _clean_str(cluster_info.get("cluster_type", "")),
                    "treatment_type": table.get("treatment_type", []),
                    "treatment": _clean_str(table.get("treatment", "")),
                    "omics_type": _clean_str(table.get("omics_type", "")),
                    "light_condition": _clean_str(table.get("light_condition", "")),
                    "member_count": member_count,
                    "functional_description": _clean_str(cluster_info.get("functional_description", "")),
                    "behavioral_description": _clean_str(cluster_info.get("behavioral_description", "")),
                    "peak_time_hours": cluster_info.get("peak_time_hours"),
                    "period_hours": cluster_info.get("period_hours"),
                    "experimental_context": _clean_str(table.get("experimental_context", "")),
                }

                nodes.append((cluster_id, "gene_cluster", props))

        return nodes

    def get_edges(self) -> list[tuple]:
        """Emit Gene_in_gene_cluster, Publication_has_gene_cluster, and
        Genecluster_belongs_to_organism edges."""
        edges = []

        for table_key, table in self._cluster_tables:
            csv_path = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                continue

            df = pd.read_csv(csv_path)
            cluster_col = table["cluster_col"]
            gene_id_col = "locus_tag" if "locus_tag" in df.columns else table["gene_id_col"]
            score_col = table.get("score_col")
            organism = table.get("organism", "")

            clusters_meta = table.get("clusters", {})

            for cluster_key in clusters_meta:
                cluster_id = _make_cluster_id(self.doi, self.paper_name, table_key, cluster_key)

                # Gene_in_gene_cluster edges
                mask = df[cluster_col].astype(str) == str(cluster_key)
                for _, row in df[mask].iterrows():
                    gene_locus = row.get("locus_tag", row.get(gene_id_col, ""))
                    if pd.isna(gene_locus) or not gene_locus:
                        continue
                    gene_id = f"ncbigene:{gene_locus}"

                    edge_props = {}
                    if score_col and score_col in row.index:
                        val = row[score_col]
                        if pd.notna(val):
                            edge_props["membership_score"] = float(val)

                    edge_id = f"{cluster_id}__{gene_locus}"
                    edges.append((edge_id, cluster_id, gene_id, "gene_in_gene_cluster", edge_props))

                # Publication_has_gene_cluster edge
                if self.doi:
                    pub_id = self.doi
                    pub_edge_id = f"pub_cluster__{cluster_id}"
                    edges.append((pub_edge_id, pub_id, cluster_id, "publication_has_gene_cluster", {}))

                # Genecluster_belongs_to_organism edge
                if organism:
                    # Look up organism ID from canonical name
                    org_id = self._organism_id(organism)
                    if org_id:
                        org_edge_id = f"cluster_org__{cluster_id}"
                        edges.append((org_edge_id, cluster_id, org_id, "genecluster_belongs_to_organism", {}))

        return edges

    def _organism_id(self, organism_name: str) -> str:
        """Map organism name to node ID. Organism IDs are insdc.gcf: accessions.

        This is a simplified lookup — the full mapping comes from
        cyanobacteria_genomes.csv. For now, use a static lookup of known organisms.
        """
        # This will be populated from cyanobacteria_genomes.csv at init time
        # in MultiClusterAdapter. For single adapter, return empty to skip.
        return getattr(self, "_organism_lookup", {}).get(organism_name, "")

    def download_data(self, **kwargs):
        """No download needed — cluster data is in paperconfig CSVs."""
        pass


class MultiClusterAdapter:
    """Wrapper that reads paperconfig_files.txt and delegates to ClusterAdapter instances."""

    def __init__(self, config_list_file: str, genome_config_file: str = None, **kwargs):
        """
        Args:
            config_list_file: Path to paperconfig_files.txt (one path per line).
            genome_config_file: Path to cyanobacteria_genomes.csv (for organism ID lookup).
            **kwargs: Passed to each ClusterAdapter (e.g., test_mode).
        """
        self._organism_lookup = {}
        if genome_config_file:
            self._organism_lookup = self._build_organism_lookup(genome_config_file)

        self.adapters = []
        paperconfigs = load_all_paperconfigs(Path(config_list_file))
        for pc_path, config in paperconfigs:
            # Only create adapter if this config has gene_clusters entries
            supp = config.get("publication", {}).get("supplementary_materials", {})
            has_clusters = any(
                isinstance(v, dict) and v.get("type") == "gene_clusters"
                for v in supp.values()
            )
            if has_clusters:
                adapter = ClusterAdapter(config_file=str(pc_path), **kwargs)
                adapter._organism_lookup = self._organism_lookup
                self.adapters.append(adapter)

    def _build_organism_lookup(self, genome_config_file: str) -> dict:
        """Build organism name → node ID mapping from cyanobacteria_genomes.csv."""
        lookup = {}
        try:
            df = pd.read_csv(genome_config_file)
            for _, row in df.iterrows():
                name = row.get("preferred_name", "")
                accession = row.get("assembly_accession", "")
                if name and accession:
                    lookup[name] = f"insdc.gcf:{accession}"
        except Exception as e:
            logger.warning(f"Could not load genome config: {e}")
        return lookup

    def download_data(self, **kwargs):
        pass

    def get_nodes(self) -> list[tuple]:
        nodes = []
        for adapter in self.adapters:
            nodes.extend(adapter.get_nodes())
        return nodes

    def get_edges(self) -> list[tuple]:
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges())
        return edges
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/test_cluster_adapter.py::test_cluster_adapter_get_nodes -v
```
Expected: PASS

- [ ] **Step 5: Write test for get_edges**

Add to `tests/test_cluster_adapter.py`:

```python
def test_cluster_adapter_get_edges(temp_cluster_dir):
    """ClusterAdapter emits membership edges with correct gene IDs and scores."""
    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter

    tmpdir, config_path, config = temp_cluster_dir
    adapter = ClusterAdapter(config_file=config_path)

    edges = list(adapter.get_edges())

    # Gene_in_gene_cluster edges: 5 genes total (2 in cluster 1, 3 in cluster 2)
    membership_edges = [e for e in edges if e[3] == "gene_in_gene_cluster"]
    assert len(membership_edges) == 5

    # Check gene IDs are prefixed
    gene_targets = {e[2] for e in membership_edges}
    assert all(g.startswith("ncbigene:") for g in gene_targets)
    assert "ncbigene:PMM0001" in gene_targets

    # Check membership scores are present
    scores = [e[4].get("membership_score") for e in membership_edges if "membership_score" in e[4]]
    assert len(scores) == 5
    assert all(isinstance(s, float) for s in scores)


def test_cluster_adapter_skips_unresolved_genes(temp_cluster_dir):
    """Genes with NaN locus_tag are skipped (no dangling edges)."""
    from multiomics_kg.adapters.cluster_adapter import ClusterAdapter

    tmpdir, config_path, config = temp_cluster_dir

    # Overwrite resolved CSV with one NaN locus_tag
    df = pd.DataFrame({
        "ORF": ["PMM0001", "UNKNOWN", "PMM0003"],
        "cluster": [1, 1, 2],
        "membership": [0.95, 0.87, 0.92],
        "locus_tag": ["PMM0001", float("nan"), "PMM0003"],
        "resolution_method": ["tier1:locus_tag", "unresolved", "tier1:locus_tag"],
    })
    resolved_path = os.path.join(tmpdir, "clusters_resolved.csv")
    df.to_csv(resolved_path, index=False)

    adapter = ClusterAdapter(config_file=config_path)
    edges = list(adapter.get_edges())
    membership_edges = [e for e in edges if e[3] == "gene_in_gene_cluster"]

    # Only 2 genes resolved, UNKNOWN skipped
    assert len(membership_edges) == 2
```

- [ ] **Step 6: Run edge tests**

```bash
pytest tests/test_cluster_adapter.py -v
```
Expected: ALL PASS

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/adapters/cluster_adapter.py tests/test_cluster_adapter.py
git commit -m "feat: add ClusterAdapter and MultiClusterAdapter"
```

---

### Task 4: Validation — Extend validate_paperconfig.py for gene_clusters

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py`

- [ ] **Step 1: Add vocabulary constants**

Add to the constants section near `CANONICAL_CONDITION_TYPES`:

```python
# Add new treatment_type values for cluster data
CANONICAL_CONDITION_TYPES.update({"diel", "oxygen_stress"})

# Valid cluster_type values
VALID_CLUSTER_TYPES = {
    "diel_periodicity",
    "stress_response",
    "expression_level",
}

# Required fields on gene_clusters supplementary entries
REQUIRED_CLUSTER_TABLE_FIELDS = [
    "filename", "organism", "gene_id_col", "cluster_col", "clusters",
]

# Required fields per cluster definition
REQUIRED_CLUSTER_FIELDS = ["name", "cluster_type"]

# Recommended fields per cluster definition (warn if missing)
RECOMMENDED_CLUSTER_FIELDS = [
    "functional_description", "behavioral_description",
]
```

- [ ] **Step 2: Add gene_clusters validation block**

In the main validation loop over `supplementary_materials`, add a block after the `id_translation` handler (around line 555):

```python
# gene_clusters type
if table_type == "gene_clusters":
    # Required table-level fields
    for field in REQUIRED_CLUSTER_TABLE_FIELDS:
        if field not in table:
            errors.append(f"{table_key}: gene_clusters requires '{field}'")

    # Validate organism
    organism = table.get("organism", "")
    if organism and organism not in CANONICAL_GENOMIC_ORGANISMS:
        warnings.append(
            f"{table_key}: organism '{organism}' not in canonical list"
        )

    # Validate omics_type
    omics = table.get("omics_type", "")
    if omics and omics not in VALID_TYPES:
        warnings.append(
            f"{table_key}: omics_type '{omics}' not in {VALID_TYPES}"
        )

    # Validate treatment_type array
    tt = table.get("treatment_type", [])
    if isinstance(tt, list):
        for val in tt:
            if val not in CANONICAL_CONDITION_TYPES:
                warnings.append(
                    f"{table_key}: treatment_type '{val}' not in canonical list"
                )
    elif isinstance(tt, str):
        warnings.append(f"{table_key}: treatment_type should be a list, got string")

    # Validate CSV and columns
    fn = table.get("filename", "")
    if fn:
        df, cols = _read_csv_safe(fn, ",", 0, errors, table_key)
        if cols is not None:
            for col_field in ["gene_id_col", "cluster_col", "score_col"]:
                col_name = table.get(col_field)
                if col_name and col_name not in cols:
                    errors.append(
                        f"{table_key}: {col_field} '{col_name}' not found in CSV headers {list(cols)}"
                    )

    # Validate clusters block
    clusters = table.get("clusters", {})
    if not clusters:
        warnings.append(f"{table_key}: 'clusters' block is empty")
    else:
        for ck, cv in clusters.items():
            if not isinstance(cv, dict):
                errors.append(f"{table_key}.clusters.{ck}: expected dict")
                continue
            for rf in REQUIRED_CLUSTER_FIELDS:
                if rf not in cv:
                    errors.append(f"{table_key}.clusters.{ck}: missing required field '{rf}'")
            ct = cv.get("cluster_type", "")
            if ct and ct not in VALID_CLUSTER_TYPES:
                warnings.append(
                    f"{table_key}.clusters.{ck}: cluster_type '{ct}' not in {VALID_CLUSTER_TYPES}"
                )
            for rec in RECOMMENDED_CLUSTER_FIELDS:
                if not cv.get(rec):
                    warnings.append(
                        f"{table_key}.clusters.{ck}: missing recommended field '{rec}'"
                    )
    continue
```

- [ ] **Step 3: Test validation manually**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py data/Prochlorococcus/papers_and_supp/tolonen\ 2006/paperconfig.yaml
```

Expected: passes existing validation (no gene_clusters entries yet — that comes in Task 6).

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py
git commit -m "feat: validate gene_clusters entries in paperconfig"
```

---

### Task 5: Gene ID resolution — Extend resolve_paper_ids.py for gene_clusters

**Files:**
- Modify: `multiomics_kg/download/resolve_paper_ids.py`

- [ ] **Step 1: Identify where table types are filtered**

In `resolve_paper_ids.py`, the `resolve_table()` function checks `table_config.get("type")` and skips non-csv types. Find this filter and extend it.

- [ ] **Step 2: Add gene_clusters handling**

In `resolve_table()`, add handling for `type: gene_clusters`. The gene_clusters type uses `gene_id_col` instead of `name_col`, but the resolution logic is the same — resolve each row's gene ID to a locus_tag.

```python
# In resolve_table(), after checking table_type:
table_type = table_config.get("type", "csv")

if table_type == "gene_clusters":
    # gene_clusters uses gene_id_col instead of name_col
    name_col = table_config.get("gene_id_col")
    if not name_col:
        logger.warning(f"  {table_key}: gene_clusters missing gene_id_col, skipping")
        return {"skipped": True}
    # id_columns from the table config (same as csv type)
    id_columns = table_config.get("id_columns", [])
    # Get organism from table level (not from analysis/experiment)
    organism = table_config.get("organism", "")
elif table_type == "csv":
    # existing csv handling...
```

The rest of the resolution logic (loading mapping, resolving rows, writing `_resolved.csv`) stays the same.

- [ ] **Step 3: Test resolution on a dummy cluster CSV**

Create a small test cluster CSV, add a gene_clusters entry to a paperconfig, and run:

```bash
uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Tolonen 2006" --force
```

Expected: produces `clusters_resolved.csv` alongside the original.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/download/resolve_paper_ids.py
git commit -m "feat: extend resolve_paper_ids to handle gene_clusters tables"
```

---

### Task 6: Data — Prepare Tolonen 2006 cluster data and paperconfig

**Files:**
- Create: cluster CSV(s) in `data/Prochlorococcus/papers_and_supp/tolonen 2006/`
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`

- [ ] **Step 1: Extract cluster membership data from Tolonen 2006 supplementary**

The Tolonen 2006 supplementary material contains K-means cluster assignments. Check what files already exist in the directory:

```bash
ls -la "data/Prochlorococcus/papers_and_supp/tolonen 2006/"
ls -la "data/Prochlorococcus/papers_and_supp/tolonen 2006/44320_2006_BFMSB4100087_MOESM1_ESM/"
```

The supplementary archive has a `Kmeans/` directory with `clustering.pdf` — this likely contains the cluster assignments. The actual gene-to-cluster mapping may need to be extracted from GEO or the supplementary data files. Inspect available data and extract to CSV format:

```
gene_id,cluster
PMM0001,1
PMM0002,1
...
```

- [ ] **Step 2: Add gene_clusters entry to Tolonen 2006 paperconfig**

Add under `supplementary_materials` in the existing `paperconfig.yaml`:

```yaml
  cluster_table_med4:
    type: gene_clusters
    filename: "data/Prochlorococcus/papers_and_supp/tolonen 2006/med4_kmeans_clusters.csv"
    organism: "Prochlorococcus MED4"
    gene_id_col: "gene_id"
    cluster_col: "cluster"
    cluster_method: "K-means (9 clusters)"
    omics_type: MICROARRAY
    light_condition: "continuous light"
    treatment_type: ["nitrogen_stress"]
    treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
    experimental_context: "Custom Affymetrix microarray, MED4 grown in Pro99, 21C, 14 umol photons m-2 s-1"
    clusters:
      # Fill via LLM extraction or manual reading of paper
      # Each cluster needs: name, cluster_type, functional_description, behavioral_description
```

(Per-cluster descriptions to be filled by reading the paper or using the LLM extraction tool.)

- [ ] **Step 3: Run resolve_paper_ids on the new cluster CSV**

```bash
uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Tolonen 2006" --force
```

- [ ] **Step 4: Validate paperconfig**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

Expected: passes with possible warnings about missing cluster descriptions.

- [ ] **Step 5: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/tolonen 2006/"
git commit -m "data: add Tolonen 2006 MED4 K-means cluster data"
```

---

### Task 7: Pipeline integration — Wire up MultiClusterAdapter in create_knowledge_graph.py

**Files:**
- Modify: `create_knowledge_graph.py`

- [ ] **Step 1: Add MultiClusterAdapter import and instantiation**

Add after the omics adapter section in `create_knowledge_graph.py`:

```python
from multiomics_kg.adapters.cluster_adapter import MultiClusterAdapter

# ... after omics_adapter.get_edges() ...

# Gene cluster adapter
cluster_adapter = MultiClusterAdapter(
    config_list_file='data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
    genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
    test_mode=TEST_MODE,
)
cluster_adapter.download_data(cache=CACHE)
bc.write_nodes(cluster_adapter.get_nodes())
bc.write_edges(cluster_adapter.get_edges())
```

- [ ] **Step 2: Test in test mode**

```bash
uv run python create_knowledge_graph.py --test
```

Expected: runs without errors. Check output CSV for `GeneCluster` nodes.

- [ ] **Step 3: Commit**

```bash
git add create_knowledge_graph.py
git commit -m "feat: wire MultiClusterAdapter into KG build pipeline"
```

---

### Task 8: Post-import — Add indexes and member_count verification

**Files:**
- Modify: `scripts/post-import.sh`
- Modify: `scripts/post-import.cypher`

- [ ] **Step 1: Add scalar indexes**

Add to both `post-import.sh` and `post-import.cypher`:

```cypher
// ── GeneCluster indexes ─────────────────────────────────────────────────
CREATE INDEX gene_cluster_organism_idx IF NOT EXISTS FOR (gc:GeneCluster) ON (gc.organism_name);
CREATE INDEX gene_cluster_treatment_type_idx IF NOT EXISTS FOR (gc:GeneCluster) ON (gc.treatment_type);
CREATE INDEX gene_cluster_type_idx IF NOT EXISTS FOR (gc:GeneCluster) ON (gc.cluster_type);
```

- [ ] **Step 2: Add full-text index**

```cypher
CREATE FULLTEXT INDEX geneClusterFullText IF NOT EXISTS
FOR (gc:GeneCluster)
ON EACH [gc.name, gc.functional_description, gc.behavioral_description, gc.experimental_context];
```

- [ ] **Step 3: Add member_count verification**

```cypher
// ── GeneCluster member_count verification ──────────────────────────────
// Update member_count from actual edge count (authoritative)
MATCH (gc:GeneCluster)
OPTIONAL MATCH (gc)-[r:Gene_in_gene_cluster]->()
WITH gc, count(r) AS actual_count
SET gc.member_count = actual_count;
```

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "post-import: add GeneCluster indexes and member_count verification"
```

---

### Task 9: Skill docs — Update paperconfig SKILL.md

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md`

- [ ] **Step 1: Add gene_clusters documentation**

Add a new section after the existing `id_translation` documentation:

```markdown
### Gene Cluster Tables (`type: gene_clusters`)

For papers that report co-expression clusters, diel periodicity groups, or
expression-level classifications with gene membership lists.

```yaml
  cluster_table_1:
    type: gene_clusters
    filename: "data/.../clusters.csv"
    organism: "Prochlorococcus MED4"
    gene_id_col: "ORF"              # column with gene identifiers
    cluster_col: "cluster"           # column with cluster assignment
    score_col: "membership"          # optional: fuzzy membership score
    cluster_method: "Mfuzz soft clustering"
    omics_type: MICROARRAY
    light_condition: "14:10 L:D"
    treatment_type: ["diel"]         # array — same vocabulary as experiments
    treatment: "Diel transcriptome, 2h sampling"
    experimental_context: "Custom Affymetrix array, 14:10 L:D cycle, 2 days"
    clusters:
      cluster_1:
        name: "Cluster 1"
        cluster_type: "diel_periodicity"  # diel_periodicity | stress_response | expression_level
        functional_description: "PSI and PSII genes (FDR 1.5e-9)"
        behavioral_description: "Peaks at dawn, drops through day"
        peak_time_hours: 2.0          # optional: for diel clusters
        period_hours: 24.0            # optional: for periodic clusters
```

**Required fields:** `filename`, `organism`, `gene_id_col`, `cluster_col`, `clusters`

**Per-cluster required:** `name`, `cluster_type`

**Per-cluster recommended:** `functional_description`, `behavioral_description`

**treatment_type values:** Same enum as experiments — `nitrogen_stress`, `carbon_stress`,
`diel`, `oxygen_stress`, `light_stress`, `temperature_stress`, etc. Use an array
because some clusters span multiple conditions.

Gene IDs go through the same step 4 resolution pipeline as DE tables.
```

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md
git commit -m "docs: add gene_clusters format to paperconfig skill"
```

---

### Task 10: Run full test suite and verify

**Files:** (no changes — verification only)

- [ ] **Step 1: Run unit tests**

```bash
pytest -m "not slow and not kg" -v
```

Expected: ALL PASS, including new cluster adapter tests.

- [ ] **Step 2: Run test mode build**

```bash
uv run python create_knowledge_graph.py --test
```

Expected: completes without errors. GeneCluster CSV files appear in output.

- [ ] **Step 3: Inspect output**

```bash
# Check that GeneCluster nodes were written
ls -la biocypher-log/example_knowledge_graph/ | grep -i cluster
head -5 biocypher-log/example_knowledge_graph/GeneCluster-part*.csv 2>/dev/null
head -5 biocypher-log/example_knowledge_graph/Gene_in_gene_cluster-part*.csv 2>/dev/null
```

- [ ] **Step 4: Commit any fixes**

If any issues found, fix and commit before proceeding.
