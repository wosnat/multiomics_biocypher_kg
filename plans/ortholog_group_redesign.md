# Plan: OrthologGroup Cluster Nodes — Unified Orthology Model

Replace pairwise `Gene_is_homolog_of_gene` edges, materialized
`*_expression_of_ortholog` edges, **and** `Cyanorak_cluster` nodes /
`Gene_in_cyanorak_cluster` edges with a single `OrthologGroup` node type and
`Gene_in_ortholog_group` edge type, distinguished by properties.

Expression propagation moves to query-time in the explorer.

Reference: `/home/osnat/github/multiomics_explorer/plans/redefine_mcp_tools/kg_homolog_redesign.md`

---

## Current State

### What exists today

| Element | Count | Source |
|---|---|---|
| `Cyanorak_cluster` nodes | 5,619 | CyanorakNcbi adapter |
| `Gene_in_cyanorak_cluster` edges | 20,657 | CyanorakNcbi adapter |
| `Gene_is_homolog_of_gene` edges | 365,840 | post-import.sh |
| `Condition_changes_expression_of_ortholog` | 2,005,209 | post-import.sh |
| `Coculture_changes_expression_of_ortholog` | 82,782 | post-import.sh |

### What replaces them

| New element | Estimated count | Source |
|---|---|---|
| `OrthologGroup` nodes (cyanorak) | ~5,600 | adapter |
| `OrthologGroup` nodes (eggnog lowest-level) | ~3,000 | adapter |
| `OrthologGroup` nodes (eggnog bacteria-level) | ~2,700 | adapter |
| **Total OrthologGroup nodes** | **~11.3K** | |
| `Gene_in_ortholog_group` edges (cyanorak) | ~20,600 | adapter |
| `Gene_in_ortholog_group` edges (eggnog lowest-level) | ~28,000 | adapter |
| `Gene_in_ortholog_group` edges (eggnog bacteria-level) | ~22,000 | adapter |
| **Total membership edges** | **~70K** | |

**Net: 2.47M edges removed, ~81K nodes+edges added.**

---

## Data Sources

### 1. Cyanorak clusters (Pro/Syn only)

From `gene_annotations_merged.json` → `cluster_number` field.
- 5,619 unique clusters, 20,657 gene memberships across 10 Pro/Syn strains
- Alteromonas strains have no Cyanorak clusters
- Cluster IDs: `CK_00000001`, `CK_00000003`, etc.

### 2. EggNOG OGs (all strains)

From `gene_annotations_merged.json` → `eggnog_ogs` list field.
At-format entries: `OG_ID@taxon_id|level_name` (pipe preserved in JSON; only
converted to `,` by `clean_text()` when written to Neo4j gene node properties).

**Two levels extracted per gene:**

| Level | How to find | Example | Coverage |
|---|---|---|---|
| Bacteria (COG) | Entry with `@2\|Bacteria` in `eggnog_ogs` | `COG0592@2\|Bacteria` | ~70% of genes |
| Lowest-level | Entry matching the organism group's **target taxon_id** (whitelist lookup) | `1MKTR@1212\|Prochloraceae` | ~81–89% of genes |

**⚠ Lowest-level selection uses a whitelist, NOT max(taxon_id).** NCBI taxon IDs
are not monotonically ordered by depth, and cross-lineage OGs appear in the data
(e.g., a MED4 gene mapping to Pleurocapsales@52604 instead of Prochloraceae@1212).

Lowest-level per organism group (whitelist):

| Organism group | data_dir contains | Target level | Taxon ID | Fallback level | Fallback Taxon ID |
|---|---|---|---|---|---|
| Prochlorococcus (8 strains) | `Prochlorococcus` | Prochloraceae | 1212 | Cyanobacteria | 1117 |
| Synechococcus/Parasynechococcus (2 strains) | `Synechococcus` | Synechococcus | 1129 | Cyanobacteria | 1117 |
| Alteromonas (3 strains) | `Alteromonas` | Alteromonadaceae | 72275 | Gammaproteobacteria | 1236 |

Observed coverage (target level / fallback / bacteria-only / no-eggnog):

| Strain | Target | Fallback | Bacteria-only | No eggnog |
|---|---|---|---|---|
| MED4 (n=1976) | 88.6% | 2.4% | 2.3% | 6.7% |
| CC9311 (n=3052) | 81.1% | 0.9% | 1.2% | 16.7% |
| MIT1002 (n=4028) | 85.0% | 8.6% | 2.9% | 3.4% |

---

## Unified OrthologGroup Schema

### Node: `OrthologGroup`

| Property | Type | Values |
|---|---|---|
| `id` | str | Unique ID (see scheme below) |
| `name` | str | Raw OG identifier without prefix (e.g., `"CK_00000364"`, `"COG0592@2"`, `"1MKTR@1212"`) |
| `source` | str | `"cyanorak"` or `"eggnog"` |
| `taxonomic_level` | str | `"curated"`, `"Prochloraceae"`, `"Synechococcus"`, `"Alteromonadaceae"`, `"Bacteria"`, etc. |
| `taxon_id` | int | NCBI taxon ID of the level (0 for cyanorak curated) |

### Edge: `Gene_in_ortholog_group`

```
(Gene)-[:Gene_in_ortholog_group]->(OrthologGroup)
```

A gene may have 1–3 memberships:
- Cyanorak cluster (if Pro/Syn with cluster_number)
- Lowest-level eggNOG OG
- Bacteria-level COG (if different from lowest-level)

### ID Scheme

| Source | ID format | Example |
|---|---|---|
| Cyanorak cluster | `cyanorak:<cluster_number>` | `cyanorak:CK_00000001` |
| eggNOG bacteria-level | `eggnog:<COG>@2` | `eggnog:COG0592@2` |
| eggNOG lowest-level | `eggnog:<OG>@<taxon_id>` | `eggnog:1MKTR@1212` |

---

## Phased Implementation

Each phase stops for user approval before proceeding to the next.

### Agent Assignments

| Agent | Role | Phases |
|---|---|---|
| **Agent B** (schema-adapters) | Annotation pipeline, schema, adapter, pipeline, post-import code changes | 1, 2, 3 |
| **Agent C** (tests) | Unit tests, KG validity tests, snapshot | 1, 2, 3 |
| **Agent D** (skills) | Update cypher-queries and paperconfig skills | 4 |
| **Agent E** (docs) | Update CLAUDE.md, memory | 4 |
| **Agent F** (validation-runner) | Run validation after KG rebuilds | 2, 3 |
| **Agent G** (code-review) | Review changes at each phase gate | 1, 2, 3 |
| **Agent H** (plan-manager) | Track phase progress, record gate results | 1, 2, 3, 4 |

### Execution Flow

```
Phase 1:  Agent B (extraction utils + annotation pipeline) → Agent C (unit tests)
          → rebuild annotations → verify → Agent G (code review)
          ⏸ STOP

Phase 2:  Agent B (schema + adapter + pipeline + CyanorakNcbi cleanup)
          → Agent C (unit tests) → rebuild KG → Agent F (validation)
          → Agent G (code review)
          ⏸ STOP

Phase 3:  Agent B (post-import) → rebuild → Agent C (KG tests + snapshot)
          → Agent F (validation) → Agent G (code review)
          ⏸ STOP

Phase 4:  Agent D (skills) ∥ Agent E (docs)  [parallel]
          ⏸ DONE
```

---

### Phase 1: Build Annotations

Add `ortholog_groups` field to `gene_annotations_merged.json` for all strains.
No adapter/schema/graph changes yet — this phase is purely data pipeline.

**1.1 Config** — `config/gene_annotations_config.yaml` `[Agent B]`

Add ortholog_groups as a computed field (documentary entry — the actual logic
lives in the extraction module, called by `build_gene_annotations.py` post-merge):

```yaml
# ── Ortholog group memberships ──────────────────────────────────────────────
# Computed post-merge by extract_ortholog_groups(). Each gene gets a list of
# {og_id, source, taxonomic_level, taxon_id} dicts. Not a merge-rule field.
ortholog_groups:
  type: computed
  description: "List of ortholog group memberships, extracted from eggnog_ogs and cluster_number"
```

**1.2 Extraction logic** — `multiomics_kg/download/utils/ortholog_group_utils.py` (new file) `[Agent B]`

Standalone module with `extract_ortholog_groups()`. Called by
`build_gene_annotations.py` post-merge to write `ortholog_groups` into each
gene's entry in `gene_annotations_merged.json`:

```python
# Whitelist: organism group → (target_taxon_id, fallback_taxon_id)
# Determined from actual eggnog_ogs data. The "lowest level" is the most specific
# OG within the organism's OWN lineage — NOT max(taxon_id), which would pick
# cross-lineage OGs (e.g., Pleurocapsales for a Prochlorococcus gene).
ORGANISM_GROUP_LEVELS = {
    "Prochlorococcus": (1212, 1117),   # Prochloraceae, fallback Cyanobacteria
    "Synechococcus":   (1129, 1117),   # Synechococcus, fallback Cyanobacteria
    "Alteromonas":     (72275, 1236),  # Alteromonadaceae, fallback Gammaproteobacteria
}


def organism_group_from_path(data_dir: str) -> str:
    """Derive organism group from data_dir path.

    E.g. 'cache/data/Prochlorococcus/genomes/MED4/' → 'Prochlorococcus'
    Note: WH8102 (Parasynechococcus) is under Synechococcus/ in data_dir.
    """
    for group in ORGANISM_GROUP_LEVELS:
        if group in data_dir:
            return group
    return "unknown"


def _parse_eggnog_ogs(eggnog_ogs: list[str]) -> dict[int, tuple[str, str]]:
    """Parse eggnog_ogs entries into {taxon_id: (og_id, level_name)}.

    Format in gene_annotations_merged.json is "OG_ID@taxon_id|level_name"
    (pipe-separated — clean_text is NOT applied to annotation JSON).
    Entries without '@' are legacy short names and are skipped.
    """
    parsed = {}
    for entry in eggnog_ogs:
        if "@" not in entry:
            continue
        og_part, rest = entry.split("@", 1)
        parts = rest.split("|", 1)
        if len(parts) != 2:
            continue
        try:
            taxon_id = int(parts[0])
        except ValueError:
            continue
        parsed[taxon_id] = (og_part, parts[1])
    return parsed


def extract_ortholog_groups(gene: dict, organism_group: str) -> list[dict]:
    """Extract ortholog group memberships from a gene's annotations.

    Args:
        gene: Gene dict from gene_annotations_merged.json.
        organism_group: One of "Prochlorococcus", "Synechococcus", "Alteromonas"
            (derived from data_dir path). Used to select the correct target
            taxonomic level for lowest-level OG.

    Returns:
        List of {og_id, source, taxonomic_level, taxon_id} dicts. Deduplicated
        by og_id (a gene cannot belong to the same OG twice).
    """
    groups = []
    seen_ids = set()

    # 1. Cyanorak cluster (Pro/Syn only — Alt genes won't have cluster_number)
    cluster = gene.get("cluster_number")
    if cluster:
        og_id = f"cyanorak:{cluster}"
        groups.append({
            "og_id": og_id,
            "source": "cyanorak",
            "taxonomic_level": "curated",
            "taxon_id": 0,
        })
        seen_ids.add(og_id)

    # 2. Parse eggnog_ogs for bacteria-level and lowest-level
    eggnog_ogs = gene.get("eggnog_ogs") or []
    parsed = _parse_eggnog_ogs(eggnog_ogs)

    # 2a. Bacteria-level COG
    if 2 in parsed:
        og_part, level_name = parsed[2]
        og_id = f"eggnog:{og_part}@2"
        if og_id not in seen_ids:
            groups.append({
                "og_id": og_id,
                "source": "eggnog",
                "taxonomic_level": level_name,
                "taxon_id": 2,
            })
            seen_ids.add(og_id)

    # 2b. Lowest-level OG — whitelist lookup with fallback
    target_tid, fallback_tid = ORGANISM_GROUP_LEVELS.get(
        organism_group, (None, None)
    )
    lowest = None
    if target_tid and target_tid in parsed:
        lowest = (target_tid, *parsed[target_tid])
    elif fallback_tid and fallback_tid in parsed:
        lowest = (fallback_tid, *parsed[fallback_tid])

    if lowest:
        tid, og_part, level_name = lowest
        og_id = f"eggnog:{og_part}@{tid}"
        if og_id not in seen_ids:
            groups.append({
                "og_id": og_id,
                "source": "eggnog",
                "taxonomic_level": level_name,
                "taxon_id": tid,
            })
            seen_ids.add(og_id)

    return groups
```

Wire into `build_gene_annotations.py` — after main merge loop, per strain:

```python
from multiomics_kg.download.utils.ortholog_group_utils import (
    extract_ortholog_groups,
    organism_group_from_path,
)

organism_group = organism_group_from_path(str(genome_dir))
for locus_tag, gene in merged.items():
    gene["ortholog_groups"] = extract_ortholog_groups(gene, organism_group)
```

**Expected output** in `gene_annotations_merged.json`:

```json
{
  "PMM0001": {
    "ortholog_groups": [
      {"og_id": "cyanorak:CK_00000001", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0},
      {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2},
      {"og_id": "eggnog:1MKTR@1212", "source": "eggnog", "taxonomic_level": "Prochloraceae", "taxon_id": 1212}
    ],
    ...
  }
}
```

**1.3 Unit tests** — `tests/test_ortholog_group_extraction.py` (new file) `[Agent C]`

- `test_extract_ortholog_groups_pro_gene` — Pro gene with cluster_number + eggnog_ogs → 3 groups
- `test_extract_ortholog_groups_alt_gene` — Alt gene without cluster → 2 groups (bacteria + lowest)
- `test_extract_ortholog_groups_no_eggnog` — gene with only cluster_number → 1 group
- `test_extract_ortholog_groups_empty` — gene with no annotations → empty list
- `test_extract_ortholog_groups_legacy_format_skipped` — entries without `@` ignored
- `test_lowest_level_uses_whitelist_not_max_taxon_id` — gene with cross-lineage OG (e.g., Pleurocapsales@52604) → picks Prochloraceae@1212, not Pleurocapsales
- `test_lowest_level_falls_back_to_cyanobacteria` — gene missing target level but has Cyanobacteria@1117 → uses fallback
- `test_dedup_og_ids` — gene with duplicate OG entries → no duplicate in output

**1.4 Rebuild annotations** — `bash scripts/prepare_data.sh --steps 2` `[Agent B]`

Run for all 13 strains. Verify `ortholog_groups` field present in output:

```bash
uv run python -c "
import json
with open('cache/data/Prochlorococcus/genomes/MED4/gene_annotations_merged.json') as f:
    data = json.load(f)
sample = data['PMM0001']
print(json.dumps(sample.get('ortholog_groups', []), indent=2))
# Expect 3 entries: cyanorak + bacteria COG + Prochloraceae-level
"
```

**1.5 Code review** `[Agent G]` — verify:
- `extract_ortholog_groups()` handles all edge cases (legacy format, missing fields, cross-lineage OGs)
- Whitelist covers all 3 organism groups with correct taxon IDs
- Dedup logic prevents duplicate OG entries
- Config YAML entry is syntactically valid
- `ortholog_groups` field present and correctly structured in rebuilt JSON

**Phase 1 acceptance criteria:**
- [ ] `extract_ortholog_groups()` unit tests pass (including whitelist, fallback, dedup tests)
- [ ] Unit tests verify Pro genes → up to 3 groups, Alt genes → up to 2 groups
- [ ] All 13 strains have `ortholog_groups` in `gene_annotations_merged.json`
- [ ] Pro genes have up to 3 groups (cyanorak + 2 eggnog levels)
- [ ] Alt genes have up to 2 groups (no cyanorak)
- [ ] No existing fields or tests broken: `pytest -m "not slow and not kg"` all green
- [ ] Code review passed

**⏸ STOP — wait for user approval before Phase 2**

---

### Phase 2: Schema, Adapter, Pipeline

Add OrthologGroup nodes and Gene_in_ortholog_group edges to the KG via a new
adapter. Remove old Cyanorak_cluster emission. Verify with KG rebuild.

**2.1 Schema** — `config/schema_config.yaml` `[Agent B]`

**Add:**

```yaml
ortholog group:
  is_a: GroupingClass
  represented_as: node
  preferred_id: id
  label_in_input: ortholog_group
  properties:
    source: str            # "cyanorak" | "eggnog"
    taxonomic_level: str   # "curated" | "Prochloraceae" | "Alteromonadaceae" | "Bacteria" | etc.
    taxon_id: int          # NCBI taxon ID of the level (0 for curated)

gene in ortholog group:
  is_a: Association
  represented_as: edge
  label_as_edge: gene_in_ortholog_group
  source: gene
  target: ortholog group
  label_in_input: gene_in_ortholog_group
```

**Remove these schema entries:**

- `cyanorak_cluster` (node) — lines 140–147
- `gene_in_cyanorak_cluster` (edge) — lines 148–154
- `gene to gene homology association` (edge) — lines 758–768
- `coculture to ortholog gene expression association` (edge) — lines 77–103
- `condition to ortholog gene expression association` (edge) — lines 104–129

**Remove these gene properties** (now covered by OrthologGroup edges):

- `cluster_number` — line 215
- `alteromonadaceae_og` — line 236
- `bacteria_cog_og` — line 237

**2.2 New adapter** — `multiomics_kg/adapters/ortholog_group_adapter.py` `[Agent B]`

Follow the `functional_annotation_adapter.py` pattern:

```python
class OrthologGroupAdapter:
    """Per-strain: reads pre-computed ortholog_groups from gene_annotations_merged.json."""

    def __init__(self, genome_dir: Path, test_mode: bool = False):
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes = {}
        self._load()

    def _load(self):
        path = self.genome_dir / "gene_annotations_merged.json"
        if path.exists():
            with open(path) as fh:
                self._genes = json.load(fh)

    def get_og_memberships(self) -> list[tuple[str, dict]]:
        """Return (locus_tag, og_dict) pairs for all genes.

        Reads pre-computed ortholog_groups field from JSON
        (written by build_gene_annotations.py in Phase 1).
        """
        results = []
        for lt, gene in self._genes.items():
            for og in gene.get("ortholog_groups", []):
                results.append((lt, og))
            if self.test_mode and len(results) >= 100:
                break
        return results


class MultiOrthologGroupAdapter:
    """Multi-strain: yields OrthologGroup nodes + Gene_in_ortholog_group edges."""

    def __init__(self, genome_config_file: str, test_mode: bool = False):
        self.adapters: list[OrthologGroupAdapter] = []
        self._build_adapters(genome_config_file, test_mode)
        self.test_mode = test_mode

    def _build_adapters(self, genome_config_file: str, test_mode: bool):
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            self.adapters.append(
                OrthologGroupAdapter(genome_dir=Path(data_dir), test_mode=test_mode)
            )

    def get_nodes(self):
        """Yield unique OrthologGroup nodes across all strains."""
        seen = set()
        for adapter in self.adapters:
            for lt, og in adapter.get_og_memberships():
                og_id = og["og_id"]
                if og_id not in seen:
                    seen.add(og_id)
                    yield (
                        og_id,
                        "ortholog_group",
                        {
                            "source": og["source"],
                            "taxonomic_level": og["taxonomic_level"],
                            "taxon_id": og["taxon_id"],
                        },
                    )

    def get_edges(self):
        """Yield Gene_in_ortholog_group edges."""
        for adapter in self.adapters:
            for lt, og in adapter.get_og_memberships():
                gene_id = f"ncbigene:{lt}"
                yield (
                    f"{lt}-og-{og['og_id']}",
                    gene_id,
                    og["og_id"],
                    "gene_in_ortholog_group",
                    {},
                )
```

**2.3 Pipeline** — `create_knowledge_graph.py` `[Agent B]`

Add after CyanorakNcbi write_nodes/write_edges (line ~56):

```python
from multiomics_kg.adapters.ortholog_group_adapter import MultiOrthologGroupAdapter

og_adapter = MultiOrthologGroupAdapter(
    genome_config_file=genome_config_file,
    test_mode=test_mode,
)
bc.write_nodes(og_adapter.get_nodes())
bc.write_edges(og_adapter.get_edges())
```

**2.4 Stop old emission + remove redundant properties** — `cyanorak_ncbi_adapter.py` `[Agent B]`

- Remove `_get_cluster_nodes()` call from `get_nodes()` (line 427)
- Remove `_get_cluster_nodes()` method entirely (lines 347–383)
- Skip `GENE_IN_CLUSTER` edge type in `get_edges()` (lines 440–470)
- Remove `ClusterNodeField` enum (lines 180–181)
- Remove `cluster_node_fields` from `CyanorakNcbiModel` (line 200) and `__init__` (line 213)
- Remove `set_cluster_node_fields()` method (lines 547–551)
- Remove `CLUSTER_NUMBER`, `ALTEROMONADACEAE_OG`, `BACTERIA_COG_OG` from `GeneField` enum (lines 117, 137–138) — these properties will no longer be written to gene nodes

**Do NOT change** `config/gene_annotations_config.yaml` — the `alteromonadaceae_og`
and `bacteria_cog_og` field definitions stay. They're still built into
`gene_annotations_merged.json` (used as inputs by `extract_ortholog_groups()`).
Only the graph output changes (these fields are no longer written to Gene nodes).

**2.5 Unit tests** `[Agent C]`

**`tests/test_ortholog_group_adapter.py`** (new file):
- `test_adapter_yields_nodes_and_edges` — mock gene_annotations_merged.json
- `test_multi_adapter_deduplicates_nodes` — same OG from two strains → one node
- `test_node_properties` — verify source, taxonomic_level, taxon_id
- `test_edge_ids_unique` — no duplicate edge IDs

**`tests/test_cyanorak_ncbi_adapter.py`** (modify):
- Delete/update `cluster_nodes` assertion (line 500) — no longer emits Cyanorak_cluster nodes
- Delete/update `gene_cluster` edge assertion (line 578) — no longer emits Gene_in_cyanorak_cluster edges
- Remove assertions for `cluster_number`, `bacteria_cog_og`, `alteromonadaceae_og` gene properties

**`tests/test_build_gene_annotations.py`** — no changes needed. The
`test_alteromonadaceae_og_*` and `test_bacteria_cog_og_*` tests validate the
annotation pipeline, which is unchanged. These fields still exist in the JSON.

**2.6 Rebuild KG** — `docker compose up --build` `[Agent F]`

At this point the KG will have:
- OrthologGroup nodes + Gene_in_ortholog_group edges (new)
- No Cyanorak_cluster nodes or Gene_in_cyanorak_cluster edges
- Still has Gene_is_homolog_of_gene + *_expression_of_ortholog (from post-import, not yet removed)

**2.7 Code review** `[Agent G]` — verify:
- String literals consistent across schema, adapter, tests
- No old `cyanorak_cluster` or `gene_in_cyanorak_cluster` vocabulary in adapter output
- OrthologGroup node counts match expectations (~11.3K nodes, ~70K edges)
- Existing expression edges unaffected

**Phase 2 acceptance criteria:**
- [ ] New adapter unit tests pass
- [ ] Updated CyanorakNcbi adapter unit tests pass
- [ ] `pytest -m "not slow and not kg"` all green
- [ ] KG rebuild succeeds
- [ ] OrthologGroup nodes visible in Neo4j: `MATCH (og:OrthologGroup) RETURN og.source, og.taxonomic_level, count(og)`
- [ ] Per-source node count validation (approximate thresholds):
  - `source: "cyanorak"` — ~5,600 nodes (±200)
  - `source: "eggnog", taxonomic_level: "Bacteria"` — ~2,700 nodes (±300)
  - `source: "eggnog", taxonomic_level` in Prochloraceae/Synechococcus/Alteromonadaceae — ~3,000 nodes (±300)
  - Zero nodes from any source means a bug in the extraction or adapter
- [ ] Per-source edge count validation:
  - cyanorak membership edges — ~20,600 (±500)
  - eggnog lowest-level edges — ~28,000 (±2,000)
  - eggnog bacteria-level edges — ~22,000 (±2,000)
- [ ] No Cyanorak_cluster nodes: `MATCH (c:Cyanorak_cluster) RETURN count(c)` = 0

**⏸ STOP — wait for user approval before Phase 3**

---

### Phase 3: Post-Import Cleanup, KG Rebuild, KG Tests

Remove all pairwise homolog edges and expression propagation from post-import.
Final KG rebuild. Rewrite KG validity tests.

**3.1 Post-import** — `scripts/post-import.sh` `[Agent B]`

**Remove** lines 32–148 entirely (7 blocks):
1. Cyanorak homolog edges (lines 32–53)
2. Alteromonadaceae eggNOG homolog edges (lines 55–70)
3. Cross-phylum COG homolog edges (lines 72–84)
4. Homolog edge count reporting (lines 86–87)
5. Condition expression propagation (lines 89–115)
6. Coculture expression propagation (lines 117–144)
7. Final edge count reporting (lines 146–148)

**Add** OrthologGroup indexes (after existing full-text indexes):

```cypher
CREATE INDEX ortholog_group_id_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.id);
CREATE INDEX ortholog_group_level_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.taxonomic_level);
```

Post-import.sh becomes: scalar indexes + full-text indexes + OG indexes only.

**3.2 Sync** — `scripts/post-import.cypher` `[Agent B]`

Update reference copy to match post-import.sh.

**3.3 Rebuild KG** — `docker compose up --build` `[Agent F]`

This is the definitive rebuild with no legacy homolog/ortholog edges.

**3.4 KG validity tests** — rewrite `[Agent C]`

**`tests/kg_validity/test_post_import.py`** — delete ALL existing tests
(lines 86–452), replace with:

| New test | What it checks |
|---|---|
| `test_ortholog_group_nodes_exist` | `MATCH (og:OrthologGroup) RETURN count(og)` > 5000 |
| `test_ortholog_group_sources` | `og.source` ∈ {cyanorak, eggnog} only |
| `test_ortholog_group_per_source_counts` | cyanorak >5000, eggnog bacteria >2000, eggnog lowest >2500 — catches silent zero from any source |
| `test_ortholog_group_has_taxonomic_level` | All nodes have non-null `taxonomic_level` |
| `test_gene_in_ortholog_group_edges_exist` | Membership edge count > 50000 |
| `test_gene_in_ortholog_group_per_source_counts` | cyanorak edges >18K, eggnog lowest >24K, eggnog bacteria >18K |
| `test_cyanorak_og_pro_genes_share_group` | Two Pro genes with same cluster_number → shared OG |
| `test_eggnog_lowest_alt_genes_share_group` | Two Alt genes with same alteromonadaceae_og → shared OG |
| `test_cross_phylum_cog_bridging` | Pro gene + Alt gene with same bacteria_cog_og → shared OG |
| `test_no_old_homolog_edges` | `Gene_is_homolog_of_gene` count = 0 |
| `test_no_old_ortholog_expression_edges` | `*_expression_of_ortholog` count = 0 |
| `test_no_cyanorak_cluster_nodes` | `Cyanorak_cluster` count = 0 |
| `test_no_gene_in_cyanorak_cluster_edges` | `Gene_in_cyanorak_cluster` count = 0 |

**`tests/kg_validity/test_structure.py`**:
- Replace `Cyanorak_cluster` with `OrthologGroup` in `EXPECTED_NODE_TYPES` (line 27)
- Delete `test_cyanorak_cluster_count` (lines 84–89)
- Delete `test_prochlorococcus_genes_in_cyanorak_clusters` (lines 208–228)
- Add `test_ortholog_group_count` — >5000 OrthologGroup nodes
- Add `test_prochlorococcus_genes_in_ortholog_groups` — ≥80% of Pro genes have ≥1 OG membership

**`tests/kg_validity/generate_snapshot.py`**:
- Remove `Cyanorak_cluster` from `NODE_PROPERTIES` (line 69)
- Remove `Gene_in_cyanorak_cluster`, `Gene_is_homolog_of_gene`, `*_expression_of_ortholog` from `EDGE_PROPERTIES` (lines 89–103)
- Add `OrthologGroup` to `NODE_PROPERTIES` with `["source", "taxonomic_level", "taxon_id"]`
- Add `Gene_in_ortholog_group` to `EDGE_PROPERTIES` with `[]`
- Remove special homolog-per-source sampling logic (lines 200–220)

**`tests/kg_validity/snapshot_data.json`** — regenerate after rebuild

**3.5 Run all tests** `[Agent F]`

```bash
pytest -m kg -v                    # KG validity
pytest -m "not slow and not kg" -v # unit tests still pass
```

**3.6 Code review** `[Agent G]` — verify:
- No `Gene_is_homolog_of_gene` or `*_expression_of_ortholog` edges in graph
- No `Cyanorak_cluster` nodes in graph
- Direct expression edges unaffected (same counts as before)
- OrthologGroup counts match Phase 2

**Phase 3 acceptance criteria:**
- [ ] Post-import.sh is indexes-only (no Cypher data creation)
- [ ] KG rebuild succeeds
- [ ] `Gene_is_homolog_of_gene` count = 0
- [ ] `*_expression_of_ortholog` count = 0
- [ ] `Cyanorak_cluster` count = 0
- [ ] OrthologGroup nodes and edges present and correct
- [ ] All KG validity tests pass
- [ ] All unit tests pass
- [ ] Snapshot regenerated and committed

**⏸ STOP — wait for user approval before Phase 4**

---

### Phase 4: Documentation

Update all docs, skills, and memory to reflect the new orthology model.

**4.1 CLAUDE.md** `[Agent E]`

- **Line 321** (`test_post_import.py` description): Replace with OrthologGroup test description
- **Line 334** (Node labels): Replace `Cyanorak_cluster` → `OrthologGroup`
- **Line 335** (Relationship labels): Remove `Gene_is_homolog_of_gene`,
  `Gene_in_cyanorak_cluster`, `Condition_changes_expression_of_ortholog`,
  `Coculture_changes_expression_of_ortholog`. Add `Gene_in_ortholog_group`.
- **Lines 339–341** (Key graph facts): Remove ortholog edge descriptions. Add
  OrthologGroup description (sources, levels, `name` property for raw OG identifier,
  membership edges, query-time propagation).
- **Architecture section**: Document `ORGANISM_GROUP_LEVELS` whitelist in
  `multiomics_kg/download/utils/ortholog_group_utils.py` — maps organism group
  (Prochlorococcus/Synechococcus/Alteromonas) to target + fallback taxon IDs
  for lowest-level eggNOG OG selection. Must be updated when adding new organism groups.
- **Line 348** (Post-import description): Update to reflect indexes-only.

**4.2 `.claude/skills/cypher-queries/SKILL.md`** `[Agent D]`

- **Schema tables** (lines 47, 99–106): Remove `Cyanorak_cluster`, `Gene_is_homolog_of_gene`,
  `*_expression_of_ortholog`, `Gene_in_cyanorak_cluster`. Add `OrthologGroup` and
  `Gene_in_ortholog_group`. Remove `bacteria_cog_og`, `alteromonadaceae_og` from
  Gene property list (line 42).
- **Homolog edge sources** (lines 131–137): Replace with OrthologGroup sources description.
- **Homolog query templates** (lines 382–417): Replace with 2-hop OrthologGroup patterns.
- **Ortholog expression templates** (lines 423–437): Replace with 3-hop query-time join.
- **Cyanorak cluster queries** (lines 552–562): Replace with OrthologGroup pattern.
- **Line 692** (ortholog-inferred edges note): Remove/replace.

**4.3 `.claude/skills/paperconfig/SKILL.md`** `[Agent D]`

- **Lines 255–256**: Remove `*_expression_of_ortholog` edge references.
  Replace with note about query-time propagation via OrthologGroup.

**4.4 Memory** `[Agent E]`

- Update `memory/MEMORY.md` edge counts if changed.

**4.5 Explorer migration guide** — `docs/explorer_ortholog_migration.md` `[Agent E]`

Write a standalone document for the `multiomics_explorer` project detailing the KG
schema changes and their impact on MCP tools / Cypher queries. Include:

- **Summary of removed elements**: `Gene_is_homolog_of_gene` edges (365K),
  `Condition_changes_expression_of_ortholog` edges (2M),
  `Coculture_changes_expression_of_ortholog` edges (83K),
  `Cyanorak_cluster` nodes, `Gene_in_cyanorak_cluster` edges,
  gene properties `cluster_number`, `bacteria_cog_og`, `alteromonadaceae_og`.

- **Summary of added elements**: `OrthologGroup` nodes (~11.3K),
  `Gene_in_ortholog_group` edges (~70K), with `source` and `taxonomic_level`
  properties.

- **Cypher migration table** — old query → new query for each pattern:

  | Use case | Old Cypher | New Cypher |
  |---|---|---|
  | Find homologs of a gene | `MATCH (g:Gene)-[:Gene_is_homolog_of_gene]-(h:Gene) WHERE g.locus_tag = $lt RETURN h` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene) WHERE g.locus_tag = $lt AND g <> h RETURN h, og.source, og.taxonomic_level` |
  | Find homologs at a specific level | `MATCH (g)-[:Gene_is_homolog_of_gene {source: "cyanorak"}]-(h)` | `MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "cyanorak"})<-[:Gene_in_ortholog_group]-(h) WHERE g <> h RETURN h` |
  | Ortholog expression propagation | `MATCH (c)-[e:Condition_changes_expression_of_ortholog]->(g:Gene) WHERE g.locus_tag = $lt RETURN e` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)<-[e:Condition_changes_expression_of]-(c) WHERE g.locus_tag = $lt AND g <> h RETURN e, h.locus_tag, og.source, og.taxonomic_level` |
  | Coculture ortholog propagation | `MATCH (o)-[e:Coculture_changes_expression_of_ortholog]->(g:Gene) WHERE g.locus_tag = $lt RETURN e` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)<-[e:Coculture_changes_expression_of]-(o) WHERE g.locus_tag = $lt AND g <> h RETURN e, h.locus_tag, og.source, og.taxonomic_level` |
  | Cyanorak cluster members | `MATCH (g:Gene)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster) WHERE c.cluster_number = $ck RETURN g` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "cyanorak"}) WHERE og.id = "cyanorak:" + $ck RETURN g` |
  | Cross-phylum COG bridging | `MATCH (g)-[:Gene_is_homolog_of_gene {source: "bacteria_cog"}]-(h)` | `MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "eggnog", taxonomic_level: "Bacteria"})<-[:Gene_in_ortholog_group]-(h) WHERE g <> h RETURN h` |

- **Performance notes**:
  - 2-hop homolog queries (via OrthologGroup) replace 1-hop direct edges — expect
    slightly higher latency but dramatically fewer stored edges
  - 3-hop ortholog expression (gene→OG←gene←expression) replaces 1-hop materialized
    edges — add `og.taxonomic_level` filter to control result breadth
  - Index on `OrthologGroup.id` and `OrthologGroup.taxonomic_level` available

- **Breaking changes checklist** for explorer MCP tools — list every tool/query
  that currently references the removed edge types or properties, with the
  replacement pattern

**Phase 4 acceptance criteria:**
- [ ] CLAUDE.md reflects new node/edge types, no old vocabulary
- [ ] cypher-queries skill has working OrthologGroup query templates
- [ ] paperconfig skill has no ortholog edge references
- [ ] Explorer migration guide written with complete Cypher migration table
- [ ] No remaining references to `Cyanorak_cluster`, `Gene_is_homolog_of_gene`,
      `*_expression_of_ortholog`, or `Gene_in_cyanorak_cluster` in docs/skills

---

## What Gets Removed (full list)

### Gene properties (redundant with OrthologGroup edges)

| Property | Where to remove | Notes |
|---|---|---|
| `cluster_number` | `schema_config.yaml` gene properties (line 215), `schema_config.yaml` cyanorak_cluster properties (line 146) | Was used for Cyanorak_cluster nodes + display; OG node ID encodes it |
| `bacteria_cog_og` | `schema_config.yaml` gene properties (line 237) | Covered by OrthologGroup `{source: "eggnog", taxonomic_level: "Bacteria"}` |
| `alteromonadaceae_og` | `schema_config.yaml` gene properties (line 236) | Covered by OrthologGroup `{source: "eggnog", taxonomic_level: "Alteromonadaceae"}` |

**Files affected by gene property removal:**

| File | What to change |
|---|---|
| `config/schema_config.yaml` | Remove `cluster_number`, `bacteria_cog_og`, `alteromonadaceae_og` from gene properties |
| `config/gene_annotations_config.yaml` | **No changes** — field definitions stay (annotation pipeline unchanged) |
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | Remove `CLUSTER_NUMBER`, `ALTEROMONADACEAE_OG`, `BACTERIA_COG_OG` from `GeneField` enum (lines 117, 137–138); remove `ClusterNodeField` enum entirely (lines 180–181); remove `cluster_node_fields` from `CyanorakNcbiModel` (line 200) and `__init__` (line 213); remove `set_cluster_node_fields()` (line 547–551) |
| `tests/test_build_gene_annotations.py` | **No changes** — annotation tests stay (pipeline unchanged) |
| `tests/test_cyanorak_ncbi_adapter.py` | Remove cluster_number property assertions if any; remove cluster node/edge assertions (lines 500, 578) |
| `tests/kg_validity/test_post_import.py` | Remove tests referencing these properties |
| `.claude/skills/cypher-queries/SKILL.md` | Remove `bacteria_cog_og`, `alteromonadaceae_og` from Gene property list (line 42); remove `cluster_number` from Cyanorak_cluster (line 47) |

**Note**: `extract_ortholog_groups()` runs during `build_gene_annotations.py`
(step 2), reading `cluster_number` and `eggnog_ogs` from the merged gene dict
and writing the `ortholog_groups` list back into `gene_annotations_merged.json`.
The `cluster_number`, `bacteria_cog_og`, and `alteromonadaceae_og` properties are
removed from the **graph schema** (Gene nodes in Neo4j), not from the annotation
JSON. The JSON fields continue to exist as build-time inputs.

### Schema entries

| Element | File |
|---|---|
| `cyanorak_cluster` node schema | `config/schema_config.yaml` lines 140–147 |
| `gene_in_cyanorak_cluster` edge schema | `config/schema_config.yaml` lines 148–154 |
| `gene to gene homology association` edge schema | `config/schema_config.yaml` lines 758–768 |
| `coculture to ortholog gene expression` edge schema | `config/schema_config.yaml` lines 77–103 |
| `condition to ortholog gene expression` edge schema | `config/schema_config.yaml` lines 104–129 |

### Adapter code

| Element | File |
|---|---|
| `_get_cluster_nodes()` call in `get_nodes()` | `cyanorak_ncbi_adapter.py` line 427 |
| `GENE_IN_CLUSTER` edge emission in `get_edges()` | `cyanorak_ncbi_adapter.py` lines 440–470 |
| `ClusterNodeField` enum | `cyanorak_ncbi_adapter.py` lines 180–181 |
| `cluster_node_fields` in model/init/setter | `cyanorak_ncbi_adapter.py` lines 200, 213, 547–551 |

### Post-import

| Element | File |
|---|---|
| Homolog edge creation (3 Cypher blocks) | `scripts/post-import.sh` lines 32–84 |
| Homolog count reporting | `scripts/post-import.sh` lines 86–87 |
| Expression propagation (2 Cypher blocks) | `scripts/post-import.sh` lines 89–144 |
| Final count reporting | `scripts/post-import.sh` lines 146–148 |
| Same Cypher blocks | `scripts/post-import.cypher` (reference copy) |

### Tests

| Element | File |
|---|---|
| All homolog/ortholog tests | `tests/kg_validity/test_post_import.py` lines 86–452 |
| `test_cyanorak_cluster_count` | `tests/kg_validity/test_structure.py` lines 84–89 |
| `test_prochlorococcus_genes_in_cyanorak_clusters` | `tests/kg_validity/test_structure.py` lines 208–228 |
| Cyanorak_cluster in `EXPECTED_NODE_TYPES` | `tests/kg_validity/test_structure.py` line 27 |
| Cyanorak_cluster/homolog snapshot entries | `tests/kg_validity/snapshot_data.json` |
| Cyanorak_cluster/homolog in snapshot generator | `tests/kg_validity/generate_snapshot.py` lines 69–220 |
| Cluster node/edge assertions | `tests/test_cyanorak_ncbi_adapter.py` lines 500, 578 |
| ~~`test_alteromonadaceae_og_*` tests~~ | `tests/test_build_gene_annotations.py` — **kept** (annotation pipeline unchanged) |
| ~~`test_bacteria_cog_og_*` tests~~ | `tests/test_build_gene_annotations.py` — **kept** (annotation pipeline unchanged) |

## What Stays

| Element | Why |
|---|---|
| `eggnog_ogs` gene property | Source data, contains ALL levels (not just 2 extracted); used by other analyses |
| `cluster_number` in `gene_annotations_merged.json` | Build-time input for `extract_ortholog_groups()`; just not written to graph |
| `bacteria_cog_og` in `gene_annotations_merged.json` | Stays in JSON as build artifact; just not written to graph |
| `alteromonadaceae_og` in `gene_annotations_merged.json` | Stays in JSON as build artifact; just not written to graph |
| Direct expression edges | Real experimental data, unchanged |
