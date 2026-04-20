# Non-DE Evidence — Biller 2018 Slice, Plan 3 (Post-import Cypher + KG validity + Docs + Skills)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add post-import Cypher for DerivedMetric analysis-node rollups, Experiment/Publication/OrganismTaxon/Gene DM rollups, numeric-DM rank/percentile/bucket/significance, and new indexes; add KG validity tests covering the §Success-criteria assertions filtered to DerivedMetric scope; regenerate `snapshot_data.json` with DerivedMetric anchors; extend `/omics-edge-snapshot` skill with DM counts; write the downstream comms doc and skill/CLAUDE.md updates.

**Architecture:** Extends `scripts/post-import.sh` and byte-identical `scripts/post-import.cypher` (per CLAUDE.md). New Cypher blocks follow the existing Group 1/2/3 organization — Group 1 (indexes), Group 2 (small-table aggregations), Group 3 (heavy writes via `CALL IN TRANSACTIONS`). Mid-plan Docker rebuild wires a synthetic numeric-DM paperconfig fixture (from Plan 2, Task 13) into `create_knowledge_graph.py` so numeric-DM Cypher has data to rank/bucket/score. End-of-plan Docker rebuild unwires the fixture so the regenerated `snapshot_data.json` reflects production shape.

**Tech Stack:** Neo4j 5.x + APOC core, BioCypher, Python 3.13, pytest with `@pytest.mark.kg` marker, `docker exec deploy cypher-shell` for live-graph sanity checks.

---

## File Structure

| Path | Action | Purpose |
|---|---|---|
| `config/schema_config.yaml` | Modify | Declare new post-import-computed properties on Experiment / Publication / OrganismTaxon / Gene |
| `scripts/post-import.sh` | Modify | Add 8 Cypher blocks (indexes, DM analysis-node rollups, Experiment/Publication/OrganismTaxon DM rollups, numeric rank/percentile/bucket, significance, Gene routing counts) |
| `scripts/post-import.cypher` | Modify | Byte-identical Cypher mirror of .sh changes (CLAUDE.md invariant) |
| `scripts/post-import-validate.sh` | Modify | Extend deterministic-dump tool with DM-scoped sections (used for DE-only byte-identical regression diff) |
| `create_knowledge_graph.py` | Modify (temp) | Wire `tests/fixtures/non_de/paperconfig_files.txt` into MultiObservationsAdapter mid-plan; unwire before final rebuild |
| `tests/kg_validity/test_derived_metric.py` | Create | Parent-spec §Success-criteria KG validity assertions filtered to DerivedMetric scope (boolean/categorical; real Biller 2018 data) |
| `tests/kg_validity/test_numeric_derived_metric.py` | Create | Numeric-DM rank/percentile/bucket/significance assertions against the 100-row synthetic fixture |
| `tests/kg_validity/test_post_import.py` | Modify | Add 9 new index names to `EXPECTED_INDEXES`; add empty-state-default assertions for DM rollups on nodes without DM children |
| `tests/kg_validity/generate_snapshot.py` | Modify | Add DerivedMetric + 5 new edge types to snapshot scope |
| `tests/kg_validity/snapshot_data.json` | Regenerate | Post-rebuild with retrofitted Biller 2018 evidence in new shape |
| `.claude/skills/omics-edge-snapshot/omics_edge_snapshot.py` | Modify | Capture DM edge counts per publication alongside `Changes_expression_of` |
| `.claude/skills/omics-edge-snapshot/SKILL.md` | Modify | Document DM-count extension |
| `.claude/skills/paperconfig/SKILL.md` | Modify | Add `derived_metrics_table` section + `compartment` field doc |
| `.claude/skills/cypher-queries/SKILL.md` | Modify | Add query templates for flags / classifies / quantifies edges |
| `docs/kg-changes/non-de-evidence-extension.md` | Create | Downstream comms doc — scoped to DerivedMetric shape, flags AbundanceAnalysis as follow-up |
| `CLAUDE.md` | Modify | Key-facts block: DerivedMetric node ID format, 3 new edge types, Biller 2018 retrofit, new indexes |

---

## Conventions reused verbatim from existing post-import.sh

| Convention | Example | Source |
|---|---|---|
| 3-group structure: indexes → small aggregations → heavy writes | Group 1/2/3 comment banners | `scripts/post-import.sh:25-27,107-111,264-268` |
| Group 3 heavy writes use `CALL { … } IN TRANSACTIONS OF N ROWS` | `annotation_types`, `rank_by_effect` | `scripts/post-import.sh:272-286, 302-313` |
| Empty-state defaults pass BEFORE compute pass (compute overrides) | Experiment stats defaults → computation | `scripts/post-import.sh:152-198` |
| `apoc.coll.sort(apoc.coll.toSet(...))` for list-valued rollups | `p.treatment_types = apoc.coll.sort(tts)` | `scripts/post-import.sh:137-142` |
| `[timing]` wrapper via `time cypher-shell <<'CYPHER' … CYPHER` | One per group | `scripts/post-import.sh:29, 113, 270` |
| Boolean-semantic properties as string enums `"true"` / `"false"` | `expression_status`, `rankable`, `has_p_value` | Parent spec §673 (Ingest invariants) |
| KG validity tests marked `@pytest.mark.kg`, use `run_query` fixture | `pytestmark = pytest.mark.kg` | `tests/kg_validity/conftest.py`, `test_post_import.py` |

## BioCypher-output edge labels (confirmed via live graph after Plan 2 rebuild)

| Schema `label_in_input` | Neo4j label (used in Cypher) | Reason |
|---|---|---|
| `publication_has_derived_metric` | `PublicationHasDerivedMetric` | No `label_as_edge` → BioCypher CamelCases |
| `experiment_has_derived_metric` | `ExperimentHasDerivedMetric` | Same |
| `derived_metric_belongs_to_organism` | `DerivedMetricBelongsToOrganism` | Same |
| `derived_metric_quantifies_gene` | `Derived_metric_quantifies_gene` | `label_as_edge` set → first letter capitalized |
| `derived_metric_flags_gene` | `Derived_metric_flags_gene` | Same |
| `derived_metric_classifies_gene` | `Derived_metric_classifies_gene` | Same |

Verified via `docker exec deploy cypher-shell "CALL db.relationshipTypes()"` on 2026-04-20 after Plan 2 rebuild.

## Live-graph pre-Plan-3 state (HEAD `c4ac43b` on `dev`, rebuilt 2026-04-20)

| Metric | Count |
|---|---|
| `Changes_expression_of` edges | 227,361 (unchanged from step-0 baseline) |
| `DerivedMetric` nodes | 7 (6 boolean + 1 categorical; all Biller 2018) |
| `Derived_metric_flags_gene` edges | 4,160 |
| `Derived_metric_classifies_gene` edges | 258 |
| `Derived_metric_quantifies_gene` edges | 0 (needs synthetic fixture) |
| `PublicationHasDerivedMetric` / `ExperimentHasDerivedMetric` / `DerivedMetricBelongsToOrganism` | 7 each |
| `Experiment.compartment` coverage | 166/166 (all experiments; Plan 2 Task 1 emitted `"whole_cell"` default) |

Plan 3 blocks written against Biller 2018 boolean+categorical can be live-sanity-checked immediately; numeric blocks no-op against current data until the fixture is wired in mid-plan.

---

# Tasks

### Task 0: Capture DE-only regression baseline

**Purpose:** Gives us a byte-identical diff reference. The existing `post-import-validate.sh` dumps only DE-scoped fields. Adding new DM-scoped Cypher must not change any existing field's value — running validate.sh after all Plan 3 changes must produce an identical DE-only section vs. this baseline.

**Files:**
- Output: `/tmp/before_plan3.txt` (ephemeral; not committed)

- [ ] **Step 1: Confirm Docker deploy container is running**

Run: `docker compose ps deploy`
Expected: status `Up`, port `7687/tcp` mapped

- [ ] **Step 2: Capture DE-only baseline dump**

Run: `scripts/post-import-validate.sh > /tmp/before_plan3.txt`
Expected: exit 0, file size > 50 KB (includes full Experiment / Publication / OrganismTaxon / ClusteringAnalysis / GeneCluster / BriteCategory dumps + Gene aggregates + rank subsamples)

- [ ] **Step 3: Sanity-check the baseline**

Run: `wc -l /tmp/before_plan3.txt && head -20 /tmp/before_plan3.txt`
Expected: ~5K+ lines starting with `======== INDEXES ========`

No commit — baseline is ephemeral. Plan 3's final byte-identical diff gate (Task 20) is the consumer.

---

### Task 1: Declare new post-import-computed properties in schema

**Purpose:** BioCypher doesn't enforce schema shape at import (neo4j-admin import accepts any properties on any label), but the schema file is the source-of-truth for what downstream consumers (MCP, snapshot tests, post-import-validate.sh) expect. Declaring the new computed properties keeps the schema, post-import script, and tests aligned.

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add Experiment DM rollup properties**

Add these lines inside the `experiment:` `properties:` block, right after `growth_phases` (currently the last line):

```yaml
    reports_fold_change: str         # post-import (Plan 3): "true" | "false"
    reports_derived_metric_types: str[]  # post-import (Plan 3): distinct metric_type from child DMs
    derived_metric_count: int        # post-import (Plan 3): count of child DerivedMetric nodes
    derived_metric_value_kinds: str[]  # post-import (Plan 3): distinct value_kind from child DMs
    derived_metric_gene_count: int   # post-import (Plan 3): distinct genes reached via any child DM edge
```

- [ ] **Step 2: Add Publication DM rollup properties**

Add inside the `publication:` `properties:` block, right after `growth_phases`:

```yaml
    derived_metric_count: int        # post-import (Plan 3)
    derived_metric_gene_count: int   # post-import (Plan 3): distinct genes via child DMs
    compartments: str[]              # post-import (Plan 3): distinct compartment from child Experiments
    derived_metric_types: str[]      # post-import (Plan 3): distinct metric_type from child DMs
    derived_metric_value_kinds: str[]  # post-import (Plan 3)
```

- [ ] **Step 3: Add OrganismTaxon DM rollup properties**

Find `organism taxon:` (~line 412) `properties:` block. Add after the last existing post-import-computed property:

```yaml
    derived_metric_count: int        # post-import (Plan 3)
    derived_metric_gene_count: int   # post-import (Plan 3)
    compartments: str[]              # post-import (Plan 3)
    derived_metric_types: str[]      # post-import (Plan 3)
    derived_metric_value_kinds: str[]  # post-import (Plan 3)
```

- [ ] **Step 4: Add Gene DM routing properties**

Find `gene:` `properties:` block. Add after the last existing routing property (`cluster_types`):

```yaml
    numeric_metric_count: int                # post-import (Plan 3): distinct DMs reaching via quantifies_gene
    classifier_flag_count: int               # post-import (Plan 3): distinct DMs reaching via flags_gene
    classifier_label_count: int              # post-import (Plan 3): distinct DMs reaching via classifies_gene
    numeric_metric_types_observed: str[]     # post-import (Plan 3): sorted distinct metric_type over quantifies
    classifier_flag_types_observed: str[]    # post-import (Plan 3): sorted distinct metric_type over flags
    classifier_label_types_observed: str[]   # post-import (Plan 3): sorted distinct metric_type over classifies
    compartments_observed: str[]             # post-import (Plan 3): sorted distinct compartment over all 3 DM edge types
```

- [ ] **Step 5: Lint-run the adapter smoke test to confirm schema still parses**

Run: `uv run python -c "from biocypher import BioCypher; bc = BioCypher(schema_config_path='config/schema_config.yaml', biocypher_config_path='config/biocypher_config.yaml'); print('schema OK')"`
Expected: `schema OK` (no YAML parse error, no BioCypher validation error)

- [ ] **Step 6: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema_config: declare Plan 3 post-import-computed DM rollup fields"
```

---

### Task 2: Extend Group 1 (indexes) with DM + Experiment.compartment indexes

**Purpose:** Add the 8 scalar indexes + 1 full-text index enumerated in slice spec §Post-import. Experiment.compartment index comes online now because Plan 2 Task 1 populated the property on every Experiment.

**Files:**
- Modify: `scripts/post-import.sh` (Group 1 block; current lines 29-105)
- Modify: `scripts/post-import.cypher` (byte-identical mirror, Task 9)

- [ ] **Step 1: Append new scalar + full-text indexes to Group 1**

In `scripts/post-import.sh`, find the line `// GeneCluster` in Group 1 (~line 102). Insert BEFORE that comment:

```cypher
// DerivedMetric scalar + full-text indexes
CREATE INDEX derived_metric_metric_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.metric_type);
CREATE INDEX derived_metric_value_kind_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.value_kind);
CREATE INDEX derived_metric_compartment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.compartment);
CREATE INDEX derived_metric_omics_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.omics_type);
CREATE INDEX derived_metric_treatment_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.treatment_type);
CREATE INDEX derived_metric_organism_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.organism_name);
CREATE INDEX derived_metric_experiment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.experiment_id);
CREATE FULLTEXT INDEX derivedMetricFullText IF NOT EXISTS
  FOR (dm:DerivedMetric) ON EACH [dm.name, dm.field_description];

// Experiment.compartment scalar index (adapter-emitted by Plan 2 Task 1)
CREATE INDEX experiment_compartment_idx IF NOT EXISTS FOR (e:Experiment) ON (e.compartment);
```

- [ ] **Step 2: Live-test the new index statements against the deployed graph**

Paste just the 9 new `CREATE INDEX`/`CREATE FULLTEXT INDEX` lines into cypher-shell:

```bash
docker exec -i deploy cypher-shell -u neo4j -p neo4j <<'CYPHER'
CREATE INDEX derived_metric_metric_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.metric_type);
CREATE INDEX derived_metric_value_kind_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.value_kind);
CREATE INDEX derived_metric_compartment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.compartment);
CREATE INDEX derived_metric_omics_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.omics_type);
CREATE INDEX derived_metric_treatment_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.treatment_type);
CREATE INDEX derived_metric_organism_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.organism_name);
CREATE INDEX derived_metric_experiment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.experiment_id);
CREATE FULLTEXT INDEX derivedMetricFullText IF NOT EXISTS FOR (dm:DerivedMetric) ON EACH [dm.name, dm.field_description];
CREATE INDEX experiment_compartment_idx IF NOT EXISTS FOR (e:Experiment) ON (e.compartment);
CYPHER
```

Expected: no errors; each statement returns `0 rows available`.

- [ ] **Step 3: Confirm all 9 indexes are ONLINE**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "SHOW INDEXES YIELD name, state WHERE name IN ['derived_metric_metric_type_idx','derived_metric_value_kind_idx','derived_metric_compartment_idx','derived_metric_omics_type_idx','derived_metric_treatment_type_idx','derived_metric_organism_idx','derived_metric_experiment_idx','derivedMetricFullText','experiment_compartment_idx'] RETURN name, state ORDER BY name;"
```

Expected: 9 rows, all `state = "ONLINE"`.

- [ ] **Step 4: Full-text index smoke test**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "CALL db.index.fulltext.queryNodes('derivedMetricFullText', 'periodic') YIELD node, score RETURN count(node) AS matches;"
```

Expected: `matches` ≥ 6 (all 6 boolean Biller 2018 DMs have "Periodic" in `name`).

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: add DerivedMetric + experiment_compartment indexes"
```

---

### Task 3: Post-import — DerivedMetric analysis-node rollups

**Purpose:** Parallels existing `ClusteringAnalysis.growth_phases` rollup. Each DerivedMetric gets `total_gene_count` (count of outgoing measurement edges, matching its `value_kind`) and `growth_phases` (union from parent Experiment). Lives in Group 2 (small-table aggregation) — there are only 7 DMs today, so no CALL IN TRANSACTIONS needed.

**Files:**
- Modify: `scripts/post-import.sh` (Group 2 block; between existing `ClusteringAnalysis growth_phases` and `OrganismTaxon clustering rollup`)

- [ ] **Step 1: Append DM analysis-node rollup block to Group 2**

In `scripts/post-import.sh`, find the line `// ClusteringAnalysis growth_phases (from linked experiments)` (~line 224). Insert AFTER the closing `SET ca.growth_phases = apoc.coll.sort(gps);` of that block (around line 228), BEFORE `// OrganismTaxon clustering rollup`:

```cypher
// DerivedMetric total_gene_count: count of outgoing measurement edges.
// Each DM emits exactly ONE of the 3 edge types based on its value_kind,
// so the union across types is unambiguous.
MATCH (dm:DerivedMetric)
OPTIONAL MATCH (dm)-[r:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(:Gene)
WITH dm, count(r) AS total
SET dm.total_gene_count = total;

// DerivedMetric growth_phases: union from parent Experiment
// (mirrors ClusteringAnalysis growth_phases; reads Experiment.growth_phases
// set earlier in Group 2).
MATCH (dm:DerivedMetric)
OPTIONAL MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm)
WITH dm, apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.growth_phases, [])) | s + t)) AS gps
SET dm.growth_phases = apoc.coll.sort(gps);
```

- [ ] **Step 2: Live-test against deployed graph**

```bash
docker exec -i deploy cypher-shell -u neo4j -p neo4j <<'CYPHER'
MATCH (dm:DerivedMetric)
OPTIONAL MATCH (dm)-[r:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(:Gene)
WITH dm, count(r) AS total
SET dm.total_gene_count = total;

MATCH (dm:DerivedMetric)
OPTIONAL MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm)
WITH dm, apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.growth_phases, [])) | s + t)) AS gps
SET dm.growth_phases = apoc.coll.sort(gps);
CYPHER
```

Expected: no errors; 7 + 7 rows updated.

- [ ] **Step 3: Verify computed values match expected**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric) RETURN dm.id AS id, dm.value_kind AS vk, dm.total_gene_count AS tgc, dm.growth_phases AS gps ORDER BY dm.id;"
```

Expected rows (7 DMs, Biller 2018):
- 6 boolean DMs: total_gene_count sum across all six = 4,160 (matches `Derived_metric_flags_gene` edge count); each DM's total is the per-metric row count
- 1 categorical DM: total_gene_count = 258 (matches `Derived_metric_classifies_gene` count)
- All `growth_phases` = `[]` (Biller 2018 Experiments don't set `growth_phase` on their DE edges, so parent `Experiment.growth_phases` is empty — inherited empty)

Cross-check sum:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric) RETURN dm.value_kind AS vk, sum(dm.total_gene_count) AS sum_tgc GROUP BY dm.value_kind;"
```

Expected: `boolean, 4160` and `categorical, 258`.

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: DerivedMetric total_gene_count + growth_phases rollups"
```

---

### Task 4: Post-import — Experiment DM rollups

**Purpose:** 5 new Experiment-level rollups per slice spec §Post-import. Follows existing `Experiment clustering rollup` pattern. Runs in Group 2 (166 Experiment nodes, no CALL IN TRANSACTIONS needed).

**Files:**
- Modify: `scripts/post-import.sh` (Group 2 block)

- [ ] **Step 1: Append Experiment DM rollup block to Group 2**

In `scripts/post-import.sh`, find the existing `// Experiment clustering rollup` block (~line 252). Insert AFTER its closing `SET e.clustering_analysis_count = ca_count,` (~line 261) and BEFORE the closing `CYPHER` of Group 2:

```cypher
// Experiment DM rollup defaults (empty-state; compute below overrides where children exist)
MATCH (e:Experiment)
SET e.reports_fold_change = 'false',
    e.reports_derived_metric_types = [],
    e.derived_metric_count = 0,
    e.derived_metric_value_kinds = [],
    e.derived_metric_gene_count = 0;

// Experiment reports_fold_change: 'true' iff outgoing Changes_expression_of exists
MATCH (e:Experiment)
WHERE EXISTS { (e)-[:Changes_expression_of]->() }
SET e.reports_fold_change = 'true';

// Experiment DM compute (overrides defaults)
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
WITH e,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds
SET e.derived_metric_count = dm_count,
    e.reports_derived_metric_types = apoc.coll.sort(metric_types),
    e.derived_metric_value_kinds = apoc.coll.sort(value_kinds);

// Experiment derived_metric_gene_count: distinct genes reachable via ANY child DM edge type
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(:DerivedMetric)
  -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH e, count(DISTINCT g) AS dmg_count
SET e.derived_metric_gene_count = dmg_count;
```

- [ ] **Step 2: Live-test against deployed graph**

Paste the 4 statements above into `docker exec -i deploy cypher-shell -u neo4j -p neo4j <<'CYPHER' … CYPHER`. Expected: no errors; row counts match (166 / N / 166 / 166).

- [ ] **Step 3: Verify rollup semantics**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (e:Experiment) WHERE e.derived_metric_count > 0 RETURN e.id AS id, e.reports_fold_change AS rfc, e.reports_derived_metric_types AS rdmt, e.derived_metric_count AS dmc, e.derived_metric_value_kinds AS dmvk, e.derived_metric_gene_count AS dmgc ORDER BY e.id;"
```

Expected: 4 Biller 2018 Experiments (the 4 distinct `experiment_id` values across the 7 DMs), each with:
- `reports_fold_change = "true"` (they also have DE edges)
- `reports_derived_metric_types` subset of `{"periodic_in_axenic_LD","periodic_in_axenic_extended_darkness","periodic_in_coculture_LD","periodic_in_coculture_extended_darkness","darkness_survival_class"}`
- `derived_metric_value_kinds = ["boolean"]` or `["boolean","categorical"]` for the experiment that carries S5

Empty-state check (should be 162 = 166 - 4):
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (e:Experiment) WHERE e.derived_metric_count = 0 RETURN count(e) AS without_dm, collect(DISTINCT e.derived_metric_value_kinds)[..3] AS vk_samples, collect(DISTINCT e.reports_derived_metric_types)[..3] AS mt_samples;"
```

Expected: `without_dm = 162`, `vk_samples = [[]]`, `mt_samples = [[]]` (defaults applied).

reports_fold_change check:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (e:Experiment) RETURN e.reports_fold_change AS rfc, count(e) AS cnt ORDER BY rfc;"
```

Expected: both `"true"` and `"false"` appear (most Experiments have DE; a small number — e.g., cluster-only analyses if any — may have `"false"`).

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: Experiment DM rollups (reports_fold_change, derived_metric_* 4 fields)"
```

---

### Task 5: Post-import — Publication + OrganismTaxon DM rollups

**Purpose:** Publication and OrganismTaxon get the same 5-field DM rollup shape. Combined into one task because the patterns are nearly identical — differ only in the source node and binding edge.

**Files:**
- Modify: `scripts/post-import.sh` (Group 2 block)

- [ ] **Step 1: Append Publication DM rollup block to Group 2**

In `scripts/post-import.sh`, find existing `// Publication clustering rollup` (~line 241). Insert AFTER its closing `SET p.cluster_count = total_clusters;` (~line 250):

```cypher
// Publication DM rollup defaults
MATCH (p:Publication)
SET p.derived_metric_count = 0,
    p.derived_metric_gene_count = 0,
    p.compartments = [],
    p.derived_metric_types = [],
    p.derived_metric_value_kinds = [];

// Publication DM compute
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
WITH p,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds
SET p.derived_metric_count = dm_count,
    p.derived_metric_types = apoc.coll.sort(metric_types),
    p.derived_metric_value_kinds = apoc.coll.sort(value_kinds);

// Publication compartments: from child Experiments
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
WITH p, [x IN collect(DISTINCT e.compartment) WHERE x IS NOT NULL] AS comps
SET p.compartments = apoc.coll.sort(comps);

// Publication derived_metric_gene_count
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(:DerivedMetric)
  -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH p, count(DISTINCT g) AS dmg_count
SET p.derived_metric_gene_count = dmg_count;
```

- [ ] **Step 2: Append OrganismTaxon DM rollup block to Group 2**

Insert immediately after the Publication block:

```cypher
// OrganismTaxon DM rollup defaults
MATCH (o:OrganismTaxon)
SET o.derived_metric_count = 0,
    o.derived_metric_gene_count = 0,
    o.compartments = [],
    o.derived_metric_types = [],
    o.derived_metric_value_kinds = [];

// OrganismTaxon DM compute
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
WITH o,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds,
     [x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL] AS comps
SET o.derived_metric_count = dm_count,
    o.derived_metric_types = apoc.coll.sort(metric_types),
    o.derived_metric_value_kinds = apoc.coll.sort(value_kinds),
    o.compartments = apoc.coll.sort(comps);

// OrganismTaxon derived_metric_gene_count
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
OPTIONAL MATCH (dm)-[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH o, count(DISTINCT g) AS dmg_count
SET o.derived_metric_gene_count = dmg_count;
```

- [ ] **Step 3: Live-test both blocks**

Paste both blocks into `docker exec -i deploy cypher-shell`. Expected: no errors.

- [ ] **Step 4: Verify Publication rollups**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (p:Publication {doi: '10.1128/mSystems.00040-18'}) RETURN p.derived_metric_count AS dmc, p.derived_metric_gene_count AS dmgc, p.derived_metric_types AS dmt, p.derived_metric_value_kinds AS dmvk, p.compartments AS comps;"
```

Expected for Biller 2018: `dmc = 7`, `dmgc` > 0 (distinct genes reached across 4,160 flags + 258 classifies — expected ~1,500-3,500), `dmt = ["darkness_survival_class","periodic_in_axenic_LD","periodic_in_axenic_extended_darkness","periodic_in_coculture_LD","periodic_in_coculture_extended_darkness"]` (5 types), `dmvk = ["boolean","categorical"]`, `comps = ["whole_cell"]`.

Empty-state: all other publications must have `dmc = 0`, `compartments = ["whole_cell"]` (every paper has whole_cell experiments), `dmt = []`, `dmvk = []`.

- [ ] **Step 5: Verify OrganismTaxon rollups**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (o:OrganismTaxon) WHERE o.derived_metric_count > 0 RETURN o.preferred_name AS org, o.derived_metric_count AS dmc, o.derived_metric_types AS dmt, o.derived_metric_value_kinds AS dmvk ORDER BY o.preferred_name;"
```

Expected: 2 organisms (counts verified from live Plan 2 rebuild on 2026-04-20):
- `"Prochlorococcus NATL2A"` — `dmc = 5` (4 boolean + 1 categorical: `periodic_in_axenic_LD`, `periodic_in_axenic_extended_darkness`, `periodic_in_coculture_LD`, `periodic_in_coculture_extended_darkness`, `darkness_survival_class`)
- `"Alteromonas macleodii MIT1002"` — `dmc = 2` (2 boolean from S4B: `periodic_in_coculture_LD`, `periodic_in_coculture_extended_darkness`)

Cross-check against direct edge join:

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o:OrganismTaxon) RETURN o.preferred_name AS org, count(dm) AS dm_count ORDER BY org;"
```

Rollup must match.

Empty-state: all other organisms have `dmc = 0`, `dmt = []`, `dmvk = []`, `compartments = []` (treatment/reference organisms without DM-carrying genes).

- [ ] **Step 6: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: Publication + OrganismTaxon DM rollups (5 fields each)"
```

---

### Task 6: Post-import — Numeric DM rank + percentile + bucket

**Purpose:** Per-DerivedMetric ranking of `Derived_metric_quantifies_gene` edges. Only runs when parent `rankable='true'` (boolean and categorical DMs emit the other two edge types, not quantifies — so `rankable='true'` on a categorical/boolean DM is an ingest bug, separately guarded by Plan 2 adapter validation). Lives in Group 3 (heavy writes).

**Pinned thresholds** (slice spec §Post-import): `top_decile ≥ 90`, `top_quartile ∈ [75, 90)`, `mid ∈ [25, 75)`, `low < 25`. **Percentile formula**: for N edges in a DM, edge at rank k (1 = highest value) has `percentile = 100 * (N-k) / (N-1)` when N > 1; `100.0` when N = 1.

**Files:**
- Modify: `scripts/post-import.sh` (Group 3 block)

- [ ] **Step 1: Append rank/percentile/bucket block to Group 3**

In `scripts/post-import.sh`, find the line `// closest_ortholog_group_size + closest_ortholog_genera` (~line 345). Insert BEFORE that block:

```cypher
// Numeric DM rank/percentile/bucket: ranks derived_metric_quantifies_gene edges
// grouped by DerivedMetric, only when parent DM has rankable='true'. Ties on
// value broken by Gene.locus_tag ascending (reproducibility).
// Percentile: rank 1 (highest value) -> 100.0; rank N (lowest) -> 0.0.
// Buckets pinned per slice spec §Post-import (thresholds must not drift).
MATCH (dm:DerivedMetric {rankable: 'true'})
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->(g:Gene)
  WITH r, r.value AS val, g.locus_tag AS lt
  ORDER BY val DESC, lt ASC
  WITH collect(r) AS edges, count(r) AS n
  UNWIND range(0, size(edges) - 1) AS i
  WITH edges[i] AS r, i, n,
       CASE WHEN n = 1 THEN 100.0
            ELSE 100.0 * toFloat(n - i - 1) / toFloat(n - 1)
       END AS pct
  SET r.rank_by_metric = i + 1,
      r.metric_percentile = pct,
      r.metric_bucket = CASE
        WHEN pct >= 90.0 THEN 'top_decile'
        WHEN pct >= 75.0 THEN 'top_quartile'
        WHEN pct >= 25.0 THEN 'mid'
        ELSE 'low'
      END
} IN TRANSACTIONS OF 10 ROWS;
```

- [ ] **Step 2: Live-test against current graph**

```bash
docker exec -i deploy cypher-shell -u neo4j -p neo4j <<'CYPHER'
MATCH (dm:DerivedMetric {rankable: 'true'})
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->(g:Gene)
  WITH r, r.value AS val, g.locus_tag AS lt
  ORDER BY val DESC, lt ASC
  WITH collect(r) AS edges, count(r) AS n
  UNWIND range(0, size(edges) - 1) AS i
  WITH edges[i] AS r, i, n,
       CASE WHEN n = 1 THEN 100.0
            ELSE 100.0 * toFloat(n - i - 1) / toFloat(n - 1)
       END AS pct
  SET r.rank_by_metric = i + 1,
      r.metric_percentile = pct,
      r.metric_bucket = CASE
        WHEN pct >= 90.0 THEN 'top_decile'
        WHEN pct >= 75.0 THEN 'top_quartile'
        WHEN pct >= 25.0 THEN 'mid'
        ELSE 'low'
      END
} IN TRANSACTIONS OF 10 ROWS;
CYPHER
```

Expected against current graph: no errors. Biller 2018 has zero `rankable='true'` DMs, so the outer MATCH returns 0 rows, inner CALL never executes — **0 edges updated** (correct no-op empty-state behavior). Full data validation happens after fixture-enabled rebuild (Task 11).

- [ ] **Step 3: Confirm no quantifies edges got unexpectedly touched**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH ()-[r:Derived_metric_quantifies_gene]->() RETURN count(r) AS total, count(r.rank_by_metric) AS with_rank;"
```

Expected: `total = 0`, `with_rank = 0` (no quantifies edges exist in current graph — Biller 2018 has no numeric DMs).

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: numeric DM rank/percentile/bucket (rankable=true only)"
```

---

### Task 7: Post-import — Numeric DM significance derivation

**Purpose:** Derive `significant` ('true'/'false' string enum) on `Derived_metric_quantifies_gene` edges. Only when parent DM has `has_p_value='true'` AND `p_value_threshold IS NOT NULL` AND edge has non-null `adjusted_p_value`. Mirrors existing DE `expression_status` derivation style.

**Files:**
- Modify: `scripts/post-import.sh` (Group 3 block)

- [ ] **Step 1: Append significance block to Group 3**

Insert AFTER the rank/percentile/bucket block from Task 6 (still BEFORE `// closest_ortholog_group_size` block):

```cypher
// Numeric DM significance: on derived_metric_quantifies_gene edges,
// only when parent DM has has_p_value='true' AND p_value_threshold IS NOT NULL
// AND the edge's adjusted_p_value is non-null. Left null otherwise.
MATCH (dm:DerivedMetric {has_p_value: 'true'})
WHERE dm.p_value_threshold IS NOT NULL
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->()
  WHERE r.adjusted_p_value IS NOT NULL
  SET r.significant = CASE
    WHEN r.adjusted_p_value < dm.p_value_threshold THEN 'true'
    ELSE 'false'
  END
} IN TRANSACTIONS OF 1000 ROWS;
```

- [ ] **Step 2: Live-test against current graph**

Paste block via `docker exec -i deploy cypher-shell`. Expected: no errors; 0 DMs match (current graph has no `has_p_value='true'` DMs — all 7 Biller 2018 DMs are `"false"`).

- [ ] **Step 3: Confirm no edges updated**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH ()-[r:Derived_metric_quantifies_gene]->() WHERE r.significant IS NOT NULL RETURN count(r) AS with_sig;"
```

Expected: `with_sig = 0` (correct — no quantifies edges exist yet).

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: numeric DM significance derivation (has_p_value=true only)"
```

---

### Task 8: Post-import — Gene DM routing counts

**Purpose:** Per-Gene routing counts and types-observed arrays, split by DM edge type so the MCP `gene_overview` tool can dispatch to the right detail tool. `compartments_observed` unions compartment across all 3 DM edge types. Lives in Group 3 (per-gene writes across 81K nodes) — use `CALL IN TRANSACTIONS`.

**Files:**
- Modify: `scripts/post-import.sh` (Group 3 block)

- [ ] **Step 1: Append Gene routing block to Group 3**

Insert AFTER the significance block from Task 7 (still BEFORE `// closest_ortholog_group_size` block):

```cypher
// Gene DM routing defaults + compute.
// Defaults run first (every gene gets 0 / []), then 4 OPTIONAL-MATCH passes
// override defaults only on genes reached by DM edges.

// Defaults
MATCH (g:Gene)
CALL {
  WITH g
  SET g.numeric_metric_count = 0,
      g.classifier_flag_count = 0,
      g.classifier_label_count = 0,
      g.numeric_metric_types_observed = [],
      g.classifier_flag_types_observed = [],
      g.classifier_label_types_observed = [],
      g.compartments_observed = []
} IN TRANSACTIONS OF 1000 ROWS;

// numeric_metric_count + numeric_metric_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_quantifies_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.numeric_metric_count = cnt,
      g.numeric_metric_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// classifier_flag_count + classifier_flag_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_flags_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.classifier_flag_count = cnt,
      g.classifier_flag_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// classifier_label_count + classifier_label_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_classifies_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.classifier_label_count = cnt,
      g.classifier_label_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// compartments_observed: union across all 3 DM edge types via parent DerivedMetric
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)
    -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g)
  WITH g, [x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL] AS comps
  SET g.compartments_observed = apoc.coll.sort(comps)
} IN TRANSACTIONS OF 1000 ROWS;
```

- [ ] **Step 2: Live-test all 5 passes against current graph**

Paste the whole Gene block via `docker exec -i deploy cypher-shell`. Expected: ~81K genes touched per pass; each pass completes in seconds.

- [ ] **Step 3: Verify routing counts against direct edge queries**

For Biller 2018 boolean evidence:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (g:Gene) WHERE g.classifier_flag_count > 0 RETURN count(g) AS genes_with_flags, sum(g.classifier_flag_count) AS total_flag_memberships;"
```

Compare to direct edge query:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric)-[r:Derived_metric_flags_gene]->(g:Gene) RETURN count(DISTINCT g) AS genes_with_flag_edges, count(DISTINCT [dm, g]) AS distinct_dm_gene_pairs;"
```

Expected: `genes_with_flags = genes_with_flag_edges`, and `total_flag_memberships = distinct_dm_gene_pairs` (each gene/DM pair contributes 1 to the rollup since we collect DISTINCT dm per gene).

Same check for categorical:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (g:Gene) WHERE g.classifier_label_count > 0 RETURN count(g) AS genes_with_labels; MATCH (dm:DerivedMetric)-[:Derived_metric_classifies_gene]->(g:Gene) RETURN count(DISTINCT g) AS genes_with_label_edges;"
```

- [ ] **Step 4: Verify compartments_observed**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (g:Gene) WHERE size(g.compartments_observed) > 0 RETURN g.compartments_observed AS comps, count(g) AS gene_count ORDER BY gene_count DESC;"
```

Expected: one group — `comps = ["whole_cell"]`, covering all genes reached by Biller 2018 DM edges. All other 80K+ genes have `compartments_observed = []`.

- [ ] **Step 5: Verify empty-state defaults**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (g:Gene) WHERE g.numeric_metric_count IS NULL OR g.classifier_flag_count IS NULL OR g.classifier_label_count IS NULL OR g.compartments_observed IS NULL OR g.numeric_metric_types_observed IS NULL OR g.classifier_flag_types_observed IS NULL OR g.classifier_label_types_observed IS NULL RETURN count(g) AS null_fields;"
```

Expected: `null_fields = 0`.

- [ ] **Step 6: Commit**

```bash
git add scripts/post-import.sh
git commit -m "post-import: Gene DM routing counts (numeric/flag/label counts + types_observed + compartments_observed)"
```

---

### Task 9: Mirror into post-import.cypher byte-identical; extend post-import-validate.sh

**Purpose:** CLAUDE.md invariant — `scripts/post-import.sh` and `scripts/post-import.cypher` must have identical Cypher logic. The `.cypher` file is a reference copy for non-Docker use (e.g., running post-import manually on a fresh local Neo4j). Also extend `post-import-validate.sh` with DM dump sections for future byte-identical regression diffs.

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import-validate.sh`

- [ ] **Step 1: Mirror Tasks 2-8 Cypher into post-import.cypher**

Open both files side-by-side. The `.cypher` file has NO bash wrapping (no `echo`, `time`, `neo4j start/stop`) — just raw Cypher statements organized in the same logical groups. Copy each block inserted in Tasks 2-8 into the corresponding position in `.cypher`, omitting the `cypher-shell <<'CYPHER'` / `CYPHER` wrappers.

Reference: existing `.sh` vs `.cypher` correspondence — any block wrapped in `time cypher-shell <<'CYPHER' … CYPHER` in `.sh` appears as raw Cypher (no wrapper) in `.cypher`.

- [ ] **Step 2: Diff-check the two files are logically identical**

Extract only the Cypher statements from `.sh` (strip the bash wrappers) and diff against `.cypher`:

```bash
# Extract Cypher-only from .sh
grep -vE "^(#|echo|time |neo4j |cypher-shell|CYPHER|for |if |fi|done|  Neo4j|  sleep|set |TIMEFORMAT)" scripts/post-import.sh | sed '/^$/d' > /tmp/sh_cypher.txt
grep -v "^//" scripts/post-import.cypher | sed '/^$/d' > /tmp/cypher_cypher.txt
diff /tmp/sh_cypher.txt /tmp/cypher_cypher.txt
```

Expected: empty diff (may need to tune the grep filter — goal is same Cypher logic in both).

If this filter heuristic proves unreliable, fall back to manual section-by-section review: both files organize logic into Group 1/2/3; walk each group and confirm every `CREATE`, `MATCH … SET`, etc. appears in both.

- [ ] **Step 3: Add DERIVEDMETRIC full-dump section to post-import-validate.sh**

In `scripts/post-import-validate.sh`, find the last `section "..."` call before the `======== END ========` printf (~line 228). Insert BEFORE the END printf:

```bash
section "DERIVEDMETRIC (full dump)"
CYPHER <<'CYPHER'
MATCH (dm:DerivedMetric)
RETURN dm.id AS id,
       dm.metric_type AS metric_type,
       dm.value_kind AS value_kind,
       dm.rankable AS rankable,
       dm.has_p_value AS has_p_value,
       dm.compartment AS compartment,
       dm.total_gene_count AS total_gene_count,
       apoc.coll.sort(coalesce(dm.growth_phases, [])) AS growth_phases
ORDER BY dm.id;
CYPHER

section "EXPERIMENT DM rollup dump"
CYPHER <<'CYPHER'
MATCH (e:Experiment)
WHERE e.derived_metric_count > 0
RETURN e.id AS id,
       e.reports_fold_change AS rfc,
       apoc.coll.sort(coalesce(e.reports_derived_metric_types, [])) AS metric_types,
       e.derived_metric_count AS dm_count,
       apoc.coll.sort(coalesce(e.derived_metric_value_kinds, [])) AS value_kinds,
       e.derived_metric_gene_count AS dm_gene_count
ORDER BY e.id;
CYPHER

section "PUBLICATION DM rollup dump"
CYPHER <<'CYPHER'
MATCH (p:Publication)
WHERE p.derived_metric_count > 0
RETURN p.id AS id,
       p.derived_metric_count AS dm_count,
       p.derived_metric_gene_count AS dm_gene_count,
       apoc.coll.sort(coalesce(p.compartments, [])) AS compartments,
       apoc.coll.sort(coalesce(p.derived_metric_types, [])) AS dm_types,
       apoc.coll.sort(coalesce(p.derived_metric_value_kinds, [])) AS dm_value_kinds
ORDER BY p.id;
CYPHER

section "ORGANISMTAXON DM rollup dump"
CYPHER <<'CYPHER'
MATCH (o:OrganismTaxon)
WHERE o.derived_metric_count > 0
RETURN o.id AS id,
       o.derived_metric_count AS dm_count,
       o.derived_metric_gene_count AS dm_gene_count,
       apoc.coll.sort(coalesce(o.compartments, [])) AS compartments,
       apoc.coll.sort(coalesce(o.derived_metric_types, [])) AS dm_types,
       apoc.coll.sort(coalesce(o.derived_metric_value_kinds, [])) AS dm_value_kinds
ORDER BY o.id;
CYPHER

section "DM EDGE AGGREGATES"
CYPHER <<'CYPHER'
MATCH ()-[r:Derived_metric_flags_gene]->()
RETURN 'flags' AS edge, count(r) AS total, count(r.value_flag) AS with_value,
       count(DISTINCT r.value_flag) AS distinct_value_flags
UNION ALL
MATCH ()-[r:Derived_metric_classifies_gene]->()
RETURN 'classifies' AS edge, count(r) AS total, count(r.value_text) AS with_value,
       count(DISTINCT r.value_text) AS distinct_value_flags
UNION ALL
MATCH ()-[r:Derived_metric_quantifies_gene]->()
RETURN 'quantifies' AS edge, count(r) AS total, count(r.value) AS with_value,
       count(r.rank_by_metric) AS distinct_value_flags
ORDER BY edge;
CYPHER

section "NUMERIC DM RANK SUBSAMPLE: top-5 per DerivedMetric by rank_by_metric"
CYPHER <<'CYPHER'
MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.rank_by_metric <= 5
RETURN dm.id AS dm_id,
       r.rank_by_metric AS rank,
       g.locus_tag AS locus_tag,
       r.value AS value,
       r.metric_percentile AS pct,
       r.metric_bucket AS bucket
ORDER BY dm_id, rank, locus_tag;
CYPHER

section "GENE DM ROUTING AGGREGATES"
CYPHER <<'CYPHER'
MATCH (g:Gene)
RETURN count(g) AS total_genes,
       sum(g.numeric_metric_count) AS sum_numeric_count,
       sum(g.classifier_flag_count) AS sum_flag_count,
       sum(g.classifier_label_count) AS sum_label_count,
       count(CASE WHEN size(g.compartments_observed) > 0 THEN 1 END) AS genes_with_compartment;
CYPHER
```

- [ ] **Step 4: Run extended validate.sh and diff DE-only sections against baseline**

```bash
scripts/post-import-validate.sh > /tmp/after_plan3.txt
# Filter out the new DM sections to compare only DE-scoped output
awk '/======== DERIVEDMETRIC/{skip=1} /======== GENE DM ROUTING AGGREGATES/{skip=0; next} !skip' /tmp/after_plan3.txt > /tmp/after_plan3_de_only.txt
awk '/======== DERIVEDMETRIC/{skip=1} /======== GENE DM ROUTING AGGREGATES/{skip=0; next} !skip' /tmp/before_plan3.txt > /tmp/before_plan3_de_only.txt  # noop on baseline (no DM sections yet), but symmetric
diff /tmp/before_plan3_de_only.txt /tmp/after_plan3_de_only.txt
```

Expected: empty diff on DE-scoped sections (indexes, Experiment/Publication/OrganismTaxon/ClusteringAnalysis/GeneCluster/BriteCategory dumps, Gene aggregates, rank subsamples). Any difference = regression; investigate.

Note: the INDEXES section WILL show the 9 new index names. That's expected. If no other DE-scoped field changes, the regression gate is green. Either accept the INDEXES diff as expected, or filter it from both files symmetrically.

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import-validate.sh
git commit -m "post-import: mirror DM logic into .cypher; extend validate.sh with DM sections"
```

---

### Task 10: Wire synthetic numeric-DM fixture into create_knowledge_graph.py + Docker rebuild #1

**Purpose:** Put the 100-row numeric-DM paperconfig fixture in Neo4j so Task 6-7 Cypher (rank/percentile/bucket/significance) has real data to compute on. Fixture lives at `tests/fixtures/non_de/paperconfig_files.txt` (Plan 2 Task 13). This wiring is TEMPORARY — reverted in Task 15.

**Files:**
- Modify: `create_knowledge_graph.py` (temporary; reverted in Task 15)

- [ ] **Step 1: Wire fixture into MultiObservationsAdapter config list**

In `create_knowledge_graph.py`, find the `observations_adapter = MultiObservationsAdapter(` block (~line 112). Extend its `config_list_file` list:

```python
    observations_adapter = MultiObservationsAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
            'tests/fixtures/non_de/paperconfig_files.txt',  # TEMP Plan 3 Task 10 — unwire before Task 15
        ],
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
```

Note: **ONLY** the observations_adapter gets the fixture path — NOT the MultiOMICSAdapter or MultiClusterAdapter, since the fixture paperconfig is a pure derived-metrics entry (no DE table, no clusters).

- [ ] **Step 2: Dry-run the adapter to confirm fixture is picked up**

```bash
uv run python -c "
from multiomics_kg.adapters.observations_adapter import MultiObservationsAdapter
a = MultiObservationsAdapter(
    config_list_file=[
        'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
        'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        'tests/fixtures/non_de/paperconfig_files.txt',
    ],
    genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
    test_mode=False,
)
print(f'{len(a.adapters)} adapters loaded')
nodes = list(a.get_nodes())
edges = list(a.get_edges())
print(f'{len(nodes)} DM nodes')
edge_by_type = {}
for e in edges:
    edge_by_type[e[3]] = edge_by_type.get(e[3], 0) + 1
for et, cnt in sorted(edge_by_type.items()):
    print(f'  {et}: {cnt}')
"
```

Expected:
- `len(a.adapters)` = number of paperconfigs with DM entries (Biller 2018 + synthetic = 2)
- `len(nodes)` = 7 (Biller 2018) + 3 (synthetic fixture metrics) = 10
- `derived_metric_flags_gene: 4160`, `derived_metric_classifies_gene: 258`, `derived_metric_quantifies_gene: 300` (100 rows × 3 numeric metrics), plus 10 each of the 3 binding edges = 30.

- [ ] **Step 3: Docker rebuild #1**

```bash
docker compose down
docker compose up -d --build
```

Wait for completion. Monitor with `docker compose logs -f post-process` — look for `Post-process complete` then check `deploy` is `Up`.

Expected duration: ~45-60 minutes (matches previous full rebuilds).

- [ ] **Step 4: Confirm new counts in live graph**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric) RETURN count(dm) AS dm_nodes; MATCH ()-[r:Derived_metric_quantifies_gene]->() RETURN count(r) AS quant; MATCH ()-[r:Derived_metric_flags_gene]->() RETURN count(r) AS flags; MATCH ()-[r:Derived_metric_classifies_gene]->() RETURN count(r) AS classifies; MATCH ()-[r:Changes_expression_of]->() RETURN count(r) AS ce;"
```

Expected:
- `dm_nodes = 10`
- `quant = 300` (100 × 3)
- `flags = 4160`
- `classifies = 258`
- `ce = 227361` (DE unchanged — fixture has no DE)

- [ ] **Step 5: Verify numeric Cypher blocks ran against fixture data**

Rank property coverage:
```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH ()-[r:Derived_metric_quantifies_gene]->() RETURN count(r) AS total, count(r.rank_by_metric) AS with_rank, count(r.metric_percentile) AS with_pct, count(r.metric_bucket) AS with_bucket, count(r.significant) AS with_sig;"
```

Expected (synthetic fixture has 2 rankable metrics + 1 non-rankable; 1 has_p_value metric + 2 without):
- `total = 300`
- `with_rank = 200` (only `fourier_score` and `peak_fit_r_squared` are `rankable='true'` → 2 × 100 = 200)
- `with_pct = 200`, `with_bucket = 200`
- `with_sig = 100` (only `fourier_score` has `has_p_value='true'` with threshold 0.05, applied to 100 rows)

- [ ] **Step 6: Commit (with explicit "TEMP" marker in message)**

```bash
git add create_knowledge_graph.py
git commit -m "create_knowledge_graph: TEMP wire synthetic numeric-DM fixture (Plan 3 Task 10; reverted in Task 15)"
```

---

### Task 11: Write tests/kg_validity/test_derived_metric.py (boolean/categorical + rollups)

**Purpose:** Parent-spec §Success-criteria KG validity assertions filtered to DerivedMetric scope. Covers binding cardinality, denormalization consistency, value constraints, rollup consistency — all exercised by real Biller 2018 data + synthetic fixture.

**Files:**
- Create: `tests/kg_validity/test_derived_metric.py`

- [ ] **Step 1: Create the test file**

Create `tests/kg_validity/test_derived_metric.py` with the following content:

```python
"""
KG validity tests for DerivedMetric nodes and the 6 new edge types
(3 binding + 3 measurement). Parent-spec §Success-criteria assertions
filtered to DerivedMetric scope.

Covers:
- Node presence and minimum counts
- Binding edge cardinality (1:many from parent to DM)
- Denormalized fields match parent Experiment
- Edge-target constraints (all measurement edges target :Gene)
- Per-edge-type value constraints (value_flag enum, value_text in allowed_categories)
- Rollup consistency (Experiment/Publication/OrganismTaxon/Gene) with direct query-time counts
- Empty-state defaults on nodes without DM children
- Index existence
"""

import pytest

pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Node presence
# ---------------------------------------------------------------------------

def test_derived_metric_nodes_exist(run_query):
    """DerivedMetric nodes must exist.

    Biller 2018 contributes 7 (6 boolean + 1 categorical); the synthetic
    fixture (when wired) adds 3 more. We assert the production floor (>=7)
    so this test passes both in fixture-enabled and fixture-free graphs.
    """
    result = run_query("MATCH (dm:DerivedMetric) RETURN count(dm) AS cnt")
    assert result[0]["cnt"] >= 7, (
        f"Expected at least 7 DerivedMetric nodes (Biller 2018 floor), got {result[0]['cnt']}"
    )


@pytest.mark.parametrize("prop", [
    "name",
    "metric_type",
    "value_kind",
    "experiment_id",
    "organism_name",
    "publication_doi",
    "compartment",
    "omics_type",
    "rankable",
    "has_p_value",
    "total_gene_count",
])
def test_derived_metric_required_properties(run_query, prop):
    """Every DerivedMetric must have these non-null properties."""
    result = run_query(
        f"MATCH (dm:DerivedMetric) WHERE dm.{prop} IS NULL RETURN count(dm) AS cnt"
    )
    assert result[0]["cnt"] == 0, f"Found DerivedMetric nodes with null {prop}"


def test_derived_metric_value_kind_enum(run_query):
    """value_kind must be one of numeric/boolean/categorical."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.value_kind IN ['numeric', 'boolean', 'categorical']
        RETURN collect(DISTINCT dm.value_kind) AS bad
    """)
    assert result[0]["bad"] == [], f"Unexpected value_kind values: {result[0]['bad']}"


def test_derived_metric_rankable_enum(run_query):
    """rankable must be string 'true' or 'false' (never bool, never null)."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.rankable IN ['true', 'false']
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


def test_derived_metric_has_p_value_enum(run_query):
    """has_p_value must be string 'true' or 'false'."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.has_p_value IN ['true', 'false']
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


def test_derived_metric_rankable_only_numeric(run_query):
    """rankable='true' must imply value_kind='numeric'."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE dm.rankable = 'true' AND dm.value_kind <> 'numeric'
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have rankable='true' but value_kind != 'numeric': {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Binding edge cardinality (1:many parent -> DM; exactly 1 parent per DM)
# ---------------------------------------------------------------------------

def test_every_dm_has_exactly_one_publication_parent(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (p:Publication)-[:PublicationHasDerivedMetric]->(dm)
        WITH dm, count(p) AS pc
        WHERE pc <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Publication parent: {result[0]['examples']}"
    )


def test_every_dm_has_exactly_one_experiment_parent(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm)
        WITH dm, count(e) AS ec
        WHERE ec <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Experiment parent: {result[0]['examples']}"
    )


def test_every_dm_has_exactly_one_organism(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (dm)-[:DerivedMetricBelongsToOrganism]->(o:OrganismTaxon)
        WITH dm, count(o) AS oc
        WHERE oc <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Organism parent: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Denormalized field consistency
# ---------------------------------------------------------------------------

def test_dm_denormalized_experiment_id_matches_parent(run_query):
    result = run_query("""
        MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WHERE dm.experiment_id <> e.id
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have experiment_id mismatch: {result[0]['examples']}"
    )


@pytest.mark.parametrize("field", [
    "organism_name", "compartment", "omics_type", "treatment",
    "light_condition", "experimental_context",
])
def test_dm_denormalized_scalar_matches_parent(run_query, field):
    result = run_query(f"""
        MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WHERE coalesce(dm.{field}, '') <> coalesce(e.{field}, '')
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have {field} != parent Experiment: {result[0]['examples']}"
    )


def test_dm_publication_doi_matches_parent(run_query):
    result = run_query("""
        MATCH (p:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        WHERE dm.publication_doi <> p.doi
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Edge target constraints: all measurement edges must target :Gene
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel_type", [
    "Derived_metric_quantifies_gene",
    "Derived_metric_flags_gene",
    "Derived_metric_classifies_gene",
])
def test_measurement_edges_target_gene(run_query, rel_type):
    result = run_query(f"""
        MATCH (dm:DerivedMetric)-[r:`{rel_type}`]->(t)
        WHERE NOT t:Gene
        RETURN count(*) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Per-edge-type value constraints
# ---------------------------------------------------------------------------

def test_flag_edges_value_flag_enum(run_query):
    """Every derived_metric_flags_gene edge has value_flag in {'true','false'}."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_flags_gene]->()
        WHERE r.value_flag IS NULL OR NOT r.value_flag IN ['true', 'false']
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_classify_edges_value_text_non_null(run_query):
    result = run_query("""
        MATCH ()-[r:Derived_metric_classifies_gene]->()
        WHERE r.value_text IS NULL OR r.value_text = ''
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_classify_edges_value_text_in_allowed_categories(run_query):
    """value_text on every classifies edge must be in parent DM's allowed_categories."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_classifies_gene]->()
        WHERE NOT r.value_text IN dm.allowed_categories
        RETURN count(r) AS bad, collect(r.value_text)[..5] AS bad_values
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} classifies edges have value_text outside allowed_categories: {result[0]['bad_values']}"
    )


def test_quantify_edges_value_non_null(run_query):
    """Every quantifies edge must have a non-null numeric value."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.value IS NULL
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_quantify_edges_metric_type_matches_parent(run_query):
    """Edge.metric_type must equal parent DM.metric_type on quantifies edges."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_type <> dm.metric_type
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# DerivedMetric total_gene_count consistency
# ---------------------------------------------------------------------------

def test_dm_total_gene_count_matches_edge_count(run_query):
    """total_gene_count equals actual count of outgoing measurement edges."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (dm)-[r:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->()
        WITH dm, dm.total_gene_count AS declared, count(r) AS actual
        WHERE declared <> actual
        RETURN count(dm) AS mismatched, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"{result[0]['mismatched']} DMs have total_gene_count != edge count: {result[0]['examples']}"
    )


def test_dm_emits_only_one_edge_type(run_query):
    """Each DerivedMetric emits exactly one of the 3 measurement edge types, matching its value_kind."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WITH dm,
             count { (dm)-[:Derived_metric_quantifies_gene]->() } AS n_quant,
             count { (dm)-[:Derived_metric_flags_gene]->() } AS n_flag,
             count { (dm)-[:Derived_metric_classifies_gene]->() } AS n_class
        WITH dm, dm.value_kind AS vk, n_quant, n_flag, n_class,
             CASE
               WHEN vk = 'numeric' THEN (n_flag > 0 OR n_class > 0)
               WHEN vk = 'boolean' THEN (n_quant > 0 OR n_class > 0)
               WHEN vk = 'categorical' THEN (n_quant > 0 OR n_flag > 0)
               ELSE true
             END AS wrong_type
        WHERE wrong_type
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs emit an edge type that doesn't match their value_kind: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Experiment DM rollup consistency
# ---------------------------------------------------------------------------

def test_experiment_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WITH e, e.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']} Experiments: {result[0]['examples']}"


def test_experiment_derived_metric_gene_count_matches_query(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
        WITH e, e.derived_metric_gene_count AS declared, count(DISTINCT g) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_experiment_reports_fold_change_enum(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT e.reports_fold_change IN ['true', 'false']
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0


def test_experiment_reports_fold_change_consistent(run_query):
    """reports_fold_change='true' iff Experiment has any Changes_expression_of edge."""
    result = run_query("""
        MATCH (e:Experiment)
        WITH e, e.reports_fold_change AS declared,
             EXISTS { (e)-[:Changes_expression_of]->() } AS has_de
        WHERE (declared = 'true') <> has_de
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# Publication DM rollup consistency
# ---------------------------------------------------------------------------

def test_publication_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        WITH p, p.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']}: {result[0]['examples']}"


def test_publication_derived_metric_gene_count_matches_query(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
        WITH p, p.derived_metric_gene_count AS declared, count(DISTINCT g) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_publication_compartments_matches_child_experiments(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
        WITH p, p.compartments AS declared,
             apoc.coll.sort([x IN collect(DISTINCT e.compartment) WHERE x IS NOT NULL]) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# OrganismTaxon DM rollup consistency
# ---------------------------------------------------------------------------

def test_organism_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
        WITH o, o.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(o) AS mismatched, collect(o.preferred_name)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']}: {result[0]['examples']}"


# ---------------------------------------------------------------------------
# Gene routing count consistency
# ---------------------------------------------------------------------------

def test_gene_classifier_flag_count_matches_edges(run_query):
    """On a sample of genes with flag edges, rollup equals distinct-DM count."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.classifier_flag_count > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_flags_gene]->(g)
        WITH g, g.classifier_flag_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_gene_classifier_label_count_matches_edges(run_query):
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.classifier_label_count > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_classifies_gene]->(g)
        WITH g, g.classifier_label_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_gene_compartments_observed_consistent(run_query):
    """compartments_observed equals union over all 3 DM edge types."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE size(g.compartments_observed) > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g)
        WITH g, g.compartments_observed AS declared,
             apoc.coll.sort([x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL]) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# Empty-state defaults on nodes without DM children
# ---------------------------------------------------------------------------

def test_experiment_dm_empty_state_defaults(run_query):
    """Experiments without DM children: all DM rollup fields are 0 / [], never null."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT EXISTS { (e)-[:ExperimentHasDerivedMetric]->() }
          AND (e.derived_metric_count IS NULL
               OR e.derived_metric_count <> 0
               OR e.derived_metric_value_kinds IS NULL
               OR size(e.derived_metric_value_kinds) <> 0
               OR e.reports_derived_metric_types IS NULL
               OR size(e.reports_derived_metric_types) <> 0
               OR e.derived_metric_gene_count IS NULL
               OR e.derived_metric_gene_count <> 0)
        RETURN count(e) AS bad, collect(e.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Experiments without DM children have wrong defaults: {result[0]['examples']}"
    )


def test_gene_dm_empty_state_defaults(run_query):
    """All genes have non-null DM routing defaults (never null)."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.numeric_metric_count IS NULL
           OR g.classifier_flag_count IS NULL
           OR g.classifier_label_count IS NULL
           OR g.numeric_metric_types_observed IS NULL
           OR g.classifier_flag_types_observed IS NULL
           OR g.classifier_label_types_observed IS NULL
           OR g.compartments_observed IS NULL
        RETURN count(g) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Index existence (extends test_post_import.EXPECTED_INDEXES)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("idx_name", [
    "derived_metric_metric_type_idx",
    "derived_metric_value_kind_idx",
    "derived_metric_compartment_idx",
    "derived_metric_omics_type_idx",
    "derived_metric_treatment_type_idx",
    "derived_metric_organism_idx",
    "derived_metric_experiment_idx",
    "experiment_compartment_idx",
    "derivedMetricFullText",
])
def test_plan3_index_exists(run_query, idx_name):
    result = run_query(
        f"SHOW INDEXES YIELD name, state WHERE name = '{idx_name}' "
        f"RETURN state"
    )
    assert len(result) == 1, f"Index {idx_name} not found"
    assert result[0]["state"] == "ONLINE", f"Index {idx_name} state = {result[0]['state']}"
```

- [ ] **Step 2: Run the new test file against live graph**

```bash
uv run pytest tests/kg_validity/test_derived_metric.py -v
```

Expected: all ~30 tests pass green (both fixture and Biller 2018 data present in graph).

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_derived_metric.py
git commit -m "test_derived_metric: KG validity assertions for DM nodes + 6 new edge types + rollups"
```

---

### Task 12: Write tests/kg_validity/test_numeric_derived_metric.py (fixture-driven)

**Purpose:** Dedicated file for numeric-DM assertions that need fixture data — rank contiguity, percentile range, bucket thresholds, significance gating. Separated from `test_derived_metric.py` because these assertions apply only to rankable/has_p_value DMs (the synthetic fixture's 2 rankable + 1 has_p_value metrics).

**Files:**
- Create: `tests/kg_validity/test_numeric_derived_metric.py`

- [ ] **Step 1: Create the test file**

```python
"""
KG validity tests for numeric DerivedMetric post-import computations:
rank_by_metric, metric_percentile, metric_bucket, significant.

These assertions require numeric DMs in the graph. The synthetic paperconfig
fixture at tests/fixtures/non_de/ provides deterministic expected values:
100 rows × 3 numeric metrics (fourier_score, peak_time_h, peak_fit_r_squared).

Fixture metric config (tests/fixtures/non_de/synthetic_paperconfig.yaml):
  fourier_score:      rankable=true,  has_p_value=true,  p_value_threshold=0.05
  peak_time_h:        rankable=false, has_p_value=false
  peak_fit_r_squared: rankable=true,  has_p_value=false

Expected distributions (100 rows each):
  Rankable metrics: rank 1..100 contiguous; percentile 100.0 → 0.0;
    buckets top_decile=10, top_quartile=15, mid=50, low=25 (given formula)
  Significance (fourier_score only): adj_p_value <= 0.05 for some subset
"""

import pytest

pytestmark = pytest.mark.kg


FIXTURE_DOI = "10.9999/synthetic-numeric-dm"


def _fixture_loaded(run_query):
    """Skip test file if fixture isn't in the graph."""
    result = run_query(
        "MATCH (p:Publication {doi: $doi}) RETURN count(p) AS cnt",
        doi=FIXTURE_DOI,
    )
    return result[0]["cnt"] > 0


@pytest.fixture(scope="module", autouse=True)
def skip_if_no_fixture(run_query):
    if not _fixture_loaded(run_query):
        pytest.skip(
            f"Synthetic numeric-DM fixture ({FIXTURE_DOI}) not in graph — "
            f"wire tests/fixtures/non_de/paperconfig_files.txt into create_knowledge_graph.py "
            f"and rebuild Docker"
        )


# ---------------------------------------------------------------------------
# Fixture presence sanity
# ---------------------------------------------------------------------------

def test_fixture_dm_nodes_present(run_query):
    result = run_query("""
        MATCH (p:Publication {doi: $doi})-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        RETURN dm.metric_type AS metric_type, dm.rankable AS rankable,
               dm.has_p_value AS has_p_value
        ORDER BY metric_type
    """, doi=FIXTURE_DOI)
    metric_types = {r["metric_type"] for r in result}
    assert metric_types == {"fourier_score", "peak_time_h", "peak_fit_r_squared"}


def test_fixture_quantifies_edge_count(run_query):
    """300 edges total: 100 rows × 3 numeric metrics."""
    result = run_query("""
        MATCH (p:Publication {doi: $doi})-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
          -[r:Derived_metric_quantifies_gene]->()
        RETURN count(r) AS cnt
    """, doi=FIXTURE_DOI)
    assert result[0]["cnt"] == 300


# ---------------------------------------------------------------------------
# Rank contiguity on rankable DMs
# ---------------------------------------------------------------------------

def test_rank_by_metric_only_on_rankable(run_query):
    """rank_by_metric non-null iff parent DM.rankable='true' (quantifies edges only)."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WITH
          count(CASE WHEN r.rank_by_metric IS NOT NULL AND dm.rankable <> 'true'
                     THEN 1 END) AS rank_on_non_rankable,
          count(CASE WHEN r.rank_by_metric IS NULL AND dm.rankable = 'true'
                     THEN 1 END) AS rankable_missing_rank
        RETURN rank_on_non_rankable, rankable_missing_rank
    """)
    row = result[0]
    assert row["rank_on_non_rankable"] == 0, (
        f"{row['rank_on_non_rankable']} edges have rank_by_metric set on non-rankable DM"
    )
    assert row["rankable_missing_rank"] == 0, (
        f"{row['rankable_missing_rank']} rankable-DM edges are missing rank_by_metric"
    )


def test_rank_contiguous_per_dm(run_query):
    """rank_by_metric should be 1..N contiguous per DerivedMetric."""
    result = run_query("""
        MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WITH dm.id AS dm_id, collect(r.rank_by_metric) AS ranks,
             size(collect(r)) AS n
        WITH dm_id, ranks, n, apoc.coll.sort(ranks) AS sorted
        WHERE sorted <> range(1, n)
        RETURN dm_id, n, sorted
    """)
    assert len(result) == 0, f"Non-contiguous ranks: {result}"


# ---------------------------------------------------------------------------
# Percentile and bucket consistency
# ---------------------------------------------------------------------------

def test_percentile_in_range(run_query):
    """metric_percentile ∈ [0, 100] on every quantifies edge where it's set."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_percentile IS NOT NULL
          AND (r.metric_percentile < 0.0 OR r.metric_percentile > 100.0)
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_bucket_matches_pinned_thresholds(run_query):
    """metric_bucket follows pinned thresholds: top_decile>=90, top_quartile>=75<90, mid>=25<75, low<25."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_bucket IS NOT NULL
        WITH r, r.metric_percentile AS pct, r.metric_bucket AS bucket,
             CASE
               WHEN pct >= 90.0 THEN 'top_decile'
               WHEN pct >= 75.0 THEN 'top_quartile'
               WHEN pct >= 25.0 THEN 'mid'
               ELSE 'low'
             END AS expected
        WHERE bucket <> expected
        RETURN count(r) AS bad, collect([pct, bucket, expected])[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have bucket != pinned-threshold match: {result[0]['examples']}"
    )


def test_rank_1_is_top_decile(run_query):
    """The top-ranked gene per rankable DM is always top_decile."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene {rank_by_metric: 1}]->()
        WHERE r.metric_bucket <> 'top_decile'
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_highest_value_has_rank_1(run_query):
    """For each rankable DM, the edge with highest value has rank_by_metric=1.

    Synthetic fixture: fourier_score has row 1 at value 0.98 (rank 1),
    peak_fit_r_squared has row 1 at value 1.0139 (rank 1).
    """
    result = run_query("""
        MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WITH dm, max(r.value) AS max_val
        MATCH (dm)-[r2:Derived_metric_quantifies_gene]->()
        WHERE r2.value = max_val
        WITH dm, r2
        WHERE r2.rank_by_metric <> 1
        RETURN dm.id AS dm_id, r2.rank_by_metric AS rank, r2.value AS value
    """)
    assert len(result) == 0, f"Max-value edge has rank != 1: {result}"


# ---------------------------------------------------------------------------
# Significance gating
# ---------------------------------------------------------------------------

def test_significant_only_on_has_p_value(run_query):
    """significant non-null only when parent DM.has_p_value='true' AND p_value_threshold IS NOT NULL
    AND edge.adjusted_p_value IS NOT NULL."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND (dm.has_p_value <> 'true'
               OR dm.p_value_threshold IS NULL
               OR r.adjusted_p_value IS NULL)
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_significant_computed_correctly(run_query):
    """significant='true' iff adjusted_p_value < threshold."""
    result = run_query("""
        MATCH (dm:DerivedMetric {has_p_value: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND dm.p_value_threshold IS NOT NULL
          AND r.adjusted_p_value IS NOT NULL
        WITH r, dm,
             CASE WHEN r.adjusted_p_value < dm.p_value_threshold THEN 'true' ELSE 'false' END AS expected
        WHERE r.significant <> expected
        RETURN count(r) AS bad, collect([r.adjusted_p_value, dm.p_value_threshold, r.significant, expected])[..3] AS examples
    """)
    assert result[0]["bad"] == 0, f"{result[0]['bad']}: {result[0]['examples']}"


def test_significant_enum(run_query):
    """When non-null, significant ∈ {'true','false'}."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND NOT r.significant IN ['true', 'false']
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_fourier_score_has_significance(run_query):
    """Synthetic fixture fourier_score (has_p_value='true', threshold=0.05) has 100 edges with significant set."""
    result = run_query("""
        MATCH (dm:DerivedMetric {metric_type: 'fourier_score'})
          -[r:Derived_metric_quantifies_gene]->()
        RETURN count(r) AS total, count(r.significant) AS with_sig
    """)
    assert result[0]["total"] == 100
    assert result[0]["with_sig"] == 100, (
        f"Expected all 100 fourier_score edges to have significant set; got {result[0]['with_sig']}"
    )


def test_peak_time_has_no_significance(run_query):
    """Synthetic fixture peak_time_h (has_p_value='false') should have no significant set."""
    result = run_query("""
        MATCH (dm:DerivedMetric {metric_type: 'peak_time_h'})
          -[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0
```

- [ ] **Step 2: Run the new test file against fixture-enabled graph**

```bash
uv run pytest tests/kg_validity/test_numeric_derived_metric.py -v
```

Expected: all tests green.

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_numeric_derived_metric.py
git commit -m "test_numeric_derived_metric: KG validity for rank/percentile/bucket/significance on fixture"
```

---

### Task 13: Update EXPECTED_INDEXES in test_post_import.py; run full pytest -m kg

**Purpose:** The existing `test_expected_indexes_present_and_online` test in `test_post_import.py` hard-codes the set of indexes post-import.sh creates. Adding new indexes without updating this set would silently pass but lose coverage.

**Files:**
- Modify: `tests/kg_validity/test_post_import.py`

- [ ] **Step 1: Add 9 new index names to EXPECTED_INDEXES**

In `tests/kg_validity/test_post_import.py`, find `EXPECTED_INDEXES = {` (~line 860). Add these 9 strings inside the set, following alphabetical convention:

```python
    # DerivedMetric scalar + full-text
    "derived_metric_metric_type_idx",
    "derived_metric_value_kind_idx",
    "derived_metric_compartment_idx",
    "derived_metric_omics_type_idx",
    "derived_metric_treatment_type_idx",
    "derived_metric_organism_idx",
    "derived_metric_experiment_idx",
    "derivedMetricFullText",
    # Experiment.compartment scalar (adapter-emitted)
    "experiment_compartment_idx",
```

- [ ] **Step 2: Run full KG validity suite**

```bash
uv run pytest -m kg -v
```

Expected: 100% green. New tests from Tasks 11-12 all pass; existing tests remain green.

- [ ] **Step 3: Run full non-kg suite to confirm no regressions**

```bash
uv run pytest -m "not slow and not kg" -v
```

Expected: 1,588+ tests pass (Plan 2 baseline was 1,588).

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/test_post_import.py
git commit -m "test_post_import: register 9 Plan-3 indexes in EXPECTED_INDEXES"
```

---

### Task 14: Extend /omics-edge-snapshot skill with DerivedMetric edge counts

**Purpose:** `/omics-edge-snapshot` currently tracks only `Changes_expression_of` edges. Extend it to also capture per-publication counts for all 3 DM measurement edges. Keeps the skill the single source of truth for "did the rebuild lose edges?"

**Files:**
- Modify: `.claude/skills/omics-edge-snapshot/omics_edge_snapshot.py`
- Modify: `.claude/skills/omics-edge-snapshot/SKILL.md`

- [ ] **Step 1: Add DM edge-count capture to the snapshot script**

In `omics_edge_snapshot.py`, find the `capture_snapshot()` function (~line 84). After the `per_publication_by_direction` block (~line 116), insert BEFORE the `By target gene organism` block:

```python
    # --- DM edge counts per publication + type ---
    dm_edge_queries = {
        "derived_metric_flags_gene": """
            MATCH (pub:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
              -[r:Derived_metric_flags_gene]->(g:Gene)
            RETURN pub.doi AS publication, count(r) AS edges
            ORDER BY publication
        """,
        "derived_metric_classifies_gene": """
            MATCH (pub:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
              -[r:Derived_metric_classifies_gene]->(g:Gene)
            RETURN pub.doi AS publication, count(r) AS edges
            ORDER BY publication
        """,
        "derived_metric_quantifies_gene": """
            MATCH (pub:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
              -[r:Derived_metric_quantifies_gene]->(g:Gene)
            RETURN pub.doi AS publication, count(r) AS edges
            ORDER BY publication
        """,
    }
    snapshot["dm_edges_per_publication"] = {}
    for edge_type, query in dm_edge_queries.items():
        rows = run_cypher(query)
        snapshot["dm_edges_per_publication"][edge_type] = {
            _strip_quotes(r[0]): int(r[1]) for r in rows
        }
    # Total DM edge counts (derived, for quick-look reports)
    for edge_type in dm_edge_queries:
        rows = run_cypher(f"MATCH ()-[r:{edge_type[0].upper() + edge_type[1:]}]->() RETURN count(r) AS total")
        # BioCypher capitalizes first letter of label_as_edge, matching actual Neo4j labels
        total_key = f"total_{edge_type}"
        snapshot[total_key] = int(rows[0][0]) if rows else 0
```

- [ ] **Step 2: Extend the comparison logic to report DM count changes**

`compare_snapshots()` in `omics_edge_snapshot.py:191` has signature `(old: dict, new: dict, old_name: str, new_name: str) -> int`. Add a DM totals section after the existing "By target organism" block (search for `print("\nBy target organism:")` ending around line 250, before the "Regressions" section starts):

```python
    # ---- DerivedMetric edge totals ----
    print("\nDerivedMetric edge totals:")
    for edge_type in ("derived_metric_flags_gene",
                      "derived_metric_classifies_gene",
                      "derived_metric_quantifies_gene"):
        total_key = f"total_{edge_type}"
        old_total = old.get(total_key, 0)
        new_total = new.get(total_key, 0)
        delta = new_total - old_total
        d_str = f"+{delta:,}" if delta >= 0 else f"{delta:,}"
        print(f"  {edge_type:<40} {old_total:>7,} -> {new_total:>7,}  ({d_str})")
        if delta < 0:
            regressions.append((f"<{edge_type} total>", old_total, new_total, -delta))
```

The `regressions` list already drives the regression banner + exit-code-1 logic below — appending DM-total losses keeps them in the same gate.

- [ ] **Step 3: Update SKILL.md**

In `.claude/skills/omics-edge-snapshot/SKILL.md`, update the `## What It Captures` table — add rows:

```markdown
| `total_derived_metric_flags_gene` | Total boolean-DM flag edges |
| `total_derived_metric_classifies_gene` | Total categorical-DM classify edges |
| `total_derived_metric_quantifies_gene` | Total numeric-DM quantify edges |
| `dm_edges_per_publication` | Nested dict: `{edge_type: {doi: count}}` per publication |
```

Also add a brief paragraph under `## Interpreting the Report`:

```markdown
### DerivedMetric evidence

New as of Plan 3 (non-DE evidence slice). Snapshot captures DM edge counts alongside `Changes_expression_of`. A rebuild that loses DM edges (regression) is reported analogously: IMPROVED for gained, LOST for lost.
```

- [ ] **Step 4: Run the skill against live graph, save snapshot "plan3_with_fixture"**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save plan3_with_fixture
```

Expected: snapshot file written to `.claude/skills/omics-edge-snapshot/snapshots/plan3_with_fixture.json`, includes new DM fields with non-zero values.

- [ ] **Step 5: Compare against Plan 2 baseline (after_biller2018_retrofit_removal)**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare after_biller2018_retrofit_removal
```

Expected: `Changes_expression_of` total unchanged (227,361); all DM edge counts IMPROVED (0 → > 0).

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py .claude/skills/omics-edge-snapshot/SKILL.md .claude/skills/omics-edge-snapshot/snapshots/plan3_with_fixture.json
git commit -m "omics-edge-snapshot: capture DerivedMetric edge counts + plan3_with_fixture snapshot"
```

---

### Task 15: Unwire synthetic fixture; Docker rebuild #2 (production, fixture-free)

**Purpose:** Remove the TEMP wiring from Task 10 so the production KG doesn't ship with synthetic paperconfig data. The regenerated `snapshot_data.json` (Task 16) must reflect production shape.

**Files:**
- Modify: `create_knowledge_graph.py` (revert Task 10 change)

- [ ] **Step 1: Revert the fixture path from create_knowledge_graph.py**

Remove the `'tests/fixtures/non_de/paperconfig_files.txt'` line added in Task 10 Step 1. Leave only the two production paperconfig paths:

```python
    observations_adapter = MultiObservationsAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        ],
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
```

- [ ] **Step 2: Docker rebuild #2**

```bash
docker compose down
docker compose up -d --build
```

Wait ~45-60 minutes.

- [ ] **Step 3: Confirm fixture is gone; DE + Biller 2018 DM counts unchanged**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (dm:DerivedMetric) RETURN count(dm) AS cnt; MATCH ()-[r:Derived_metric_quantifies_gene]->() RETURN count(r) AS quant; MATCH ()-[r:Derived_metric_flags_gene]->() RETURN count(r) AS flags; MATCH ()-[r:Derived_metric_classifies_gene]->() RETURN count(r) AS cls; MATCH ()-[r:Changes_expression_of]->() RETURN count(r) AS ce;"
```

Expected:
- `cnt = 7` (Biller 2018 only — fixture removed)
- `quant = 0` (fixture was the only numeric source)
- `flags = 4160`
- `cls = 258`
- `ce = 227361` (DE parity intact)

- [ ] **Step 4: Confirm fixture publication node is gone**

```bash
docker exec deploy cypher-shell -u neo4j -p neo4j --format plain "MATCH (p:Publication {doi: '10.9999/synthetic-numeric-dm'}) RETURN count(p) AS fixture_pub;"
```

Expected: `fixture_pub = 0`.

- [ ] **Step 5: Confirm `pytest -m kg` — test_numeric_derived_metric auto-skips**

```bash
uv run pytest -m kg -v 2>&1 | grep -E "test_numeric_derived_metric|passed|skipped|failed"
```

Expected: tests in `test_numeric_derived_metric.py` are SKIPPED (because the fixture publication DOI no longer matches). All other KG tests pass green.

- [ ] **Step 6: Commit (revert commit)**

```bash
git add create_knowledge_graph.py
git commit -m "create_knowledge_graph: unwire Plan 3 synthetic fixture (reverts Task 10)"
```

---

### Task 16: Extend generate_snapshot.py + regenerate snapshot_data.json

**Purpose:** The regression snapshot (`tests/kg_validity/snapshot_data.json`) samples anchor nodes + deterministic samples per label + edges. Extend to cover DerivedMetric nodes + all 6 new edge types. Regenerate against the production (fixture-free) graph.

**Files:**
- Modify: `tests/kg_validity/generate_snapshot.py`
- Regenerate: `tests/kg_validity/snapshot_data.json`

- [ ] **Step 1: Add DerivedMetric anchor nodes to generate_snapshot.py**

In `ANCHOR_NODES`, add a new entry (keep dict alphabetical/stable):

```python
    "DerivedMetric": [
        # Biller 2018 retrofitted boolean + categorical DMs — always present post-Plan-3
        "derived_metric:mSystems.00040-18:s4a_natl2a_axenic:periodic_in_axenic_LD",
        "derived_metric:mSystems.00040-18:s5_natl2a_survival:darkness_survival_class",
    ],
```

- [ ] **Step 2: Add DerivedMetric NODE_PROPERTIES entry**

```python
    "DerivedMetric": [
        "name", "metric_type", "value_kind", "experiment_id",
        "organism_name", "compartment", "omics_type", "total_gene_count",
        "rankable", "has_p_value",
    ],
```

- [ ] **Step 3: Add 6 new edge types to EDGE_PROPERTIES**

```python
    # DerivedMetric binding edges
    "PublicationHasDerivedMetric": [],
    "ExperimentHasDerivedMetric": [],
    "DerivedMetricBelongsToOrganism": [],
    # DerivedMetric measurement edges
    "Derived_metric_flags_gene": ["metric_type", "value_flag"],
    "Derived_metric_classifies_gene": ["metric_type", "value_text"],
    "Derived_metric_quantifies_gene": [
        "metric_type", "value", "adjusted_p_value",
        "rank_by_metric", "metric_percentile", "metric_bucket", "significant",
    ],
```

Note: `Derived_metric_quantifies_gene` will have 0 sample edges in a fixture-free production snapshot — that's fine; sample loop returns empty and it gets skipped. Keeping the entry means future builds that include numeric DMs will auto-populate samples.

- [ ] **Step 4: Regenerate snapshot against production graph**

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

Expected output:
```
Snapshot written to .../snapshot_data.json
  Nodes: <existing_count + ~7>  (DerivedMetric anchors + samples)
  Edges: <existing_count + ~15>  (5 samples × 3 binding edges + Biller 2018 flags/classifies samples)
```

- [ ] **Step 5: Run test_snapshot.py against regenerated snapshot**

```bash
uv run pytest tests/kg_validity/test_snapshot.py -v
```

Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add tests/kg_validity/generate_snapshot.py tests/kg_validity/snapshot_data.json
git commit -m "kg_validity: snapshot regeneration with DerivedMetric anchors + 6 new edge types"
```

---

### Task 17: docs/kg-changes/non-de-evidence-extension.md (downstream comms doc)

**Purpose:** Single-page handoff for the MCP / explorer team on the new DerivedMetric shape. Scoped to DerivedMetric — AbundanceAnalysis explicitly flagged as out-of-scope / follow-up.

**Files:**
- Create: `docs/kg-changes/non-de-evidence-extension.md`

- [ ] **Step 1: Write the doc**

Create `docs/kg-changes/non-de-evidence-extension.md`:

```markdown
# Non-DE evidence — DerivedMetric extension (Plan 3)

**Merged:** 2026-04-NN (Plan 3 completion)
**Scope:** DerivedMetric nodes + 6 new edge types. AbundanceAnalysis deferred to a follow-up slice.

## What changed

The knowledge graph now represents non-DE evidence from differential expression publications as **DerivedMetric** nodes, parallel to the existing `ClusteringAnalysis` pattern.

A DerivedMetric captures a column-level scalar summary per gene (periodicity flag, lag coefficient, fold-change category, etc.) that the paper reports **alongside or instead of** raw fold-change.

### New node type

| Node | Count (2026-04-NN rebuild) | Node ID format | Key properties |
|---|---|---|---|
| `DerivedMetric` | 7 (Biller 2018) | `derived_metric:{doi_short}:{entry_key}:{metric_type}` | `name`, `metric_type`, `value_kind ∈ {numeric, boolean, categorical}`, `rankable`, `has_p_value`, `p_value_threshold`, `unit`, `allowed_categories`, `field_description`, `total_gene_count` (post-import), `growth_phases` (post-import), plus denormalized `experiment_id` / `organism_name` / `publication_doi` / `compartment` / `omics_type` / `treatment` / `light_condition` / `experimental_context` / `treatment_type` / `background_factors` |

### New edge types

Three **binding** edges (1:many parent → DM):

| Edge | Source | Target | Cardinality |
|---|---|---|---|
| `PublicationHasDerivedMetric` | Publication | DerivedMetric | many DMs per Publication |
| `ExperimentHasDerivedMetric` | Experiment | DerivedMetric | many DMs per Experiment |
| `DerivedMetricBelongsToOrganism` | DerivedMetric | OrganismTaxon | one organism per DM |

Three **measurement** edges (DM → Gene). Each DerivedMetric emits exactly ONE edge type, chosen by its `value_kind`:

| Edge | `value_kind` | Edge properties |
|---|---|---|
| `Derived_metric_quantifies_gene` | `numeric` | `metric_type`, `value` (float), `p_value`, `adjusted_p_value`, `rank_by_metric` (post-import, only if parent `rankable="true"`), `metric_percentile` (post-import), `metric_bucket` (post-import: `top_decile`/`top_quartile`/`mid`/`low`), `significant` (post-import, only if parent `has_p_value="true"` and edge has non-null `adjusted_p_value`) |
| `Derived_metric_flags_gene` | `boolean` | `metric_type`, `value_flag ∈ {"true","false"}` |
| `Derived_metric_classifies_gene` | `categorical` | `metric_type`, `value_text` (must be in parent `allowed_categories`) |

### New `Experiment.compartment` property

Every Experiment now has `compartment` (string). Default `"whole_cell"`. Vocab: `whole_cell`, `vesicle`, `exoproteome`, `secretome` (parent spec §294).

## New post-import-computed rollups

Plan 3 post-import Cypher (`scripts/post-import.sh` + `scripts/post-import.cypher`) computes:

### Per DerivedMetric
- `total_gene_count` (int) — count of outgoing measurement edges
- `growth_phases` (str[]) — union from parent Experiment

### Per Experiment
- `reports_fold_change` (str `"true"`/`"false"`) — `"true"` iff has outgoing `Changes_expression_of`
- `reports_derived_metric_types` (str[])
- `derived_metric_count` (int)
- `derived_metric_value_kinds` (str[])
- `derived_metric_gene_count` (int) — distinct genes reachable via any DM edge type

### Per Publication
- `derived_metric_count` (int)
- `derived_metric_gene_count` (int)
- `compartments` (str[]) — from child Experiments
- `derived_metric_types` (str[])
- `derived_metric_value_kinds` (str[])

### Per OrganismTaxon
Same 5 fields as Publication.

### Per Gene (routing signals for MCP dispatch)
- `numeric_metric_count` / `classifier_flag_count` / `classifier_label_count` (int)
- `numeric_metric_types_observed` / `classifier_flag_types_observed` / `classifier_label_types_observed` (str[])
- `compartments_observed` (str[])

### Empty-state defaults

All new int properties default to `0`; all new str[] properties default to `[]`. **Never null** — downstream queries can skip null-guards.

## New indexes

- Scalar (8): `derived_metric_metric_type_idx`, `derived_metric_value_kind_idx`, `derived_metric_compartment_idx`, `derived_metric_omics_type_idx`, `derived_metric_treatment_type_idx`, `derived_metric_organism_idx`, `derived_metric_experiment_idx`, `experiment_compartment_idx`
- Full-text (1): `derivedMetricFullText` on `DerivedMetric(name, field_description)`

## Paperconfig surface

New supplementary-materials entry type: `derived_metrics_table`. See `.claude/skills/paperconfig/SKILL.md` for schema.

Biller 2018 (`10.1128/mSystems.00040-18`) is the first real-paper integration: 7 DerivedMetric nodes (6 boolean periodicity flags from Tables S4A/S4B across NATL2A axenic/coculture and MIT1002 coculture; 1 categorical darkness-survival class from Table S5). This retrofits evidence that was previously encoded as `GeneCluster` nodes pre-2026-04-20.

## Example Cypher

```cypher
// Find all genes flagged as periodic in NATL2A axenic L:D (boolean)
MATCH (dm:DerivedMetric {metric_type: 'periodic_in_axenic_LD'})
  -[r:Derived_metric_flags_gene]->(g:Gene)
WHERE r.value_flag = 'true'
  AND g.organism_name = 'Prochlorococcus NATL2A'
RETURN g.locus_tag, g.product;

// Top-decile fourier-periodicity genes across all papers (numeric, rankable)
MATCH (dm:DerivedMetric {metric_type: 'fourier_score'})
  -[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.metric_bucket = 'top_decile'
  AND r.significant = 'true'
RETURN g.organism_name, g.locus_tag, r.value, r.adjusted_p_value
ORDER BY r.value DESC;

// Classify genes by darkness-survival category (categorical)
MATCH (dm:DerivedMetric {metric_type: 'darkness_survival_class'})
  -[r:Derived_metric_classifies_gene]->(g:Gene)
RETURN r.value_text AS category, count(g) AS gene_count
ORDER BY gene_count DESC;
```

## Out of scope / follow-ups

- **AbundanceAnalysis** node (parent spec §1) — deferred to a later slice. Will cover per-sample abundance (spectral counts, copy-number-normalized transcript counts, community-proteomics fractions). Expected drivers: Biller 2022 vesicles, Oleza 2015/2017 exoproteome, Waldbauer 2012 paired RNA-seq + proteomics.
- **Metabolite layer** (parent spec Appendix) — own spec, own timeline.

## MCP surface

Unchanged in this slice. MCP-team TODO: implement `gene_derived_metrics(locus_tag)` / `derived_metric_ranked_genes(metric_type, bucket)` / `publication_derived_metrics(doi)` tools reading the post-import-computed properties above.
```

- [ ] **Step 2: Commit**

```bash
git add docs/kg-changes/non-de-evidence-extension.md
git commit -m "docs: non-DE evidence extension — DerivedMetric shape + post-import rollups"
```

---

### Task 18: Update paperconfig SKILL.md + cypher-queries SKILL.md

**Purpose:** The paperconfig wizard and query-template skills need to know about `derived_metrics_table` entries and the new query patterns. Combined into one task since both are doc-only modifications.

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md`
- Modify: `.claude/skills/cypher-queries/SKILL.md`

- [ ] **Step 1: Add `derived_metrics_table` section to paperconfig SKILL.md**

In `.claude/skills/paperconfig/SKILL.md`, find the existing table of supplementary_materials entry types (the one listing `csv`, `id_translation`, `annotation_gff`, `gene_clusters`). Add a new row:

```markdown
| `derived_metrics_table` | Column-level scalar summaries per gene (periodicity, lag, fold-change category, etc.); creates DerivedMetric nodes + one of 3 measurement edge types per metric |
```

Then add a new entry-type section AFTER the `gene_clusters` example (search for `**`gene_clusters` example**`):

````markdown
**`derived_metrics_table` example** — column-level derived metrics (periodicity flags, categories, numeric scores):
```yaml
# Boolean: periodicity flags per gene under axenic vs coculture conditions
s4a_natl2a_axenic:
  type: derived_metrics_table
  filename: "data/.../Biller 2018/table_s4a.csv"
  organism: "Prochlorococcus NATL2A"
  experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
  name_col: locus_tag
  id_columns:
    - column: locus_tag
      id_type: locus_tag
  metrics:
    - metric_type: periodic_in_axenic_LD
      name: "Periodic in NATL2A axenic L:D (Table S4A)"
      value_kind: boolean
      value_col: periodic_in_LD
      blank_policy: skip       # 'skip' | 'false' | 'true' — how to treat blank cells
      field_description: "Boolean flag: is this gene significantly periodic under 12:12 L:D light cycle in axenic NATL2A?"
    - metric_type: periodic_in_axenic_extended_darkness
      ...

# Categorical: survival classes
s5_natl2a_survival:
  type: derived_metrics_table
  filename: "data/.../Biller 2018/table_s5.csv"
  organism: "Prochlorococcus NATL2A"
  experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
  name_col: locus_tag
  id_columns:
    - column: locus_tag
      id_type: locus_tag
  metrics:
    - metric_type: darkness_survival_class
      name: "NATL2A darkness survival class (Table S5)"
      value_kind: categorical
      value_col: survival_class
      allowed_categories:
        - "strongly_upregulated"
        - "moderately_upregulated"
        - "unchanged"
        - "moderately_downregulated"
        - "strongly_downregulated"
      field_description: "Five-level survival-class assignment from Table S5"

# Numeric: rankable periodicity score with p-value
synthetic_fourier:
  type: derived_metrics_table
  filename: "..."
  organism: "..."
  experiment: some_experiment_key
  metrics:
    - metric_type: fourier_score
      value_kind: numeric
      value_col: fourier
      unit: ""
      rankable: "true"              # governs post-import rank_by_metric/percentile/bucket
      has_p_value: "true"
      p_value_col: p_value
      adjusted_p_value_col: adj_p_value
      p_value_threshold: 0.05       # used by post-import to derive `significant`
      field_description: "Fourier-transform periodicity score"
```

**Key rules:**
- `experiment:` field (top-level on the entry) links to the `experiments:` block — all metrics in this entry are attributed to that parent Experiment.
- Each metric emits ONE edge type based on `value_kind`: numeric → `Derived_metric_quantifies_gene`; boolean → `Derived_metric_flags_gene`; categorical → `Derived_metric_classifies_gene`.
- `rankable: "true"` and `has_p_value: "true"` are string enums (NEVER booleans — BioCypher booleans are unreliable). Only meaningful when `value_kind: numeric`.
- `allowed_categories` (str[]) REQUIRED when `value_kind: categorical`. Every `value_col` cell must match a category (validator enforces).
- `blank_policy` (`skip` / `false` / `true`) REQUIRED when `value_kind: boolean`. Governs how blank/NaN cells are interpreted.
- The entry key becomes the `entry_key` in the DerivedMetric node ID (`derived_metric:{doi_short}:{entry_key}:{metric_type}`). Keep short.
````

- [ ] **Step 2: Add `compartment` field doc to paperconfig SKILL.md**

Find the `experiments:` block documentation (search for `experimental_context:` or `treatment_type:`). After the `background_factors:` line, add:

```markdown
      compartment: "whole_cell"   # string, controlled vocab (optional; default "whole_cell")
                                  # Values: whole_cell | vesicle | exoproteome | secretome
                                  # Different compartments from one paper MUST be split into separate Experiments
```

- [ ] **Step 3: Add DM query templates to cypher-queries SKILL.md**

In `.claude/skills/cypher-queries/SKILL.md`, find the end of the existing query-templates section. Add a new section:

````markdown
## DerivedMetric queries (non-DE evidence)

### Boolean flags: genes flagged for a property

```cypher
MATCH (dm:DerivedMetric {metric_type: $metric_type})
  -[r:Derived_metric_flags_gene]->(g:Gene)
WHERE r.value_flag = 'true'
  AND g.organism_name = $organism
RETURN g.locus_tag, g.product
ORDER BY g.locus_tag;
```
Example: `$metric_type = "periodic_in_axenic_LD"`, `$organism = "Prochlorococcus NATL2A"` — lists genes flagged periodic under axenic L:D in NATL2A.

### Categorical labels: genes bucketed by category

```cypher
MATCH (dm:DerivedMetric {metric_type: $metric_type})
  -[r:Derived_metric_classifies_gene]->(g:Gene)
RETURN r.value_text AS category, count(g) AS gene_count
ORDER BY gene_count DESC;
```
Example: `$metric_type = "darkness_survival_class"` — gene count per survival class.

### Numeric rankings: top-N genes by rank_by_metric

```cypher
MATCH (dm:DerivedMetric {metric_type: $metric_type, rankable: 'true'})
  -[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.rank_by_metric <= $top_n
RETURN g.organism_name, g.locus_tag, r.rank_by_metric, r.value,
       r.metric_bucket, r.significant
ORDER BY r.rank_by_metric;
```

### Numeric filter by bucket + significance

```cypher
MATCH (dm:DerivedMetric {metric_type: $metric_type})
  -[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.metric_bucket IN ['top_decile', 'top_quartile']
  AND r.significant = 'true'
RETURN g.locus_tag, r.value, r.adjusted_p_value, r.metric_percentile
ORDER BY r.value DESC;
```

### All DerivedMetric evidence for one gene (routing query)

```cypher
MATCH (g:Gene {locus_tag: $locus_tag})
OPTIONAL MATCH (dm1:DerivedMetric)-[r1:Derived_metric_quantifies_gene]->(g)
OPTIONAL MATCH (dm2:DerivedMetric)-[r2:Derived_metric_flags_gene]->(g)
OPTIONAL MATCH (dm3:DerivedMetric)-[r3:Derived_metric_classifies_gene]->(g)
RETURN
  collect(DISTINCT {type: 'numeric', metric_type: dm1.metric_type, value: r1.value, bucket: r1.metric_bucket, significant: r1.significant}) AS numeric_metrics,
  collect(DISTINCT {type: 'boolean', metric_type: dm2.metric_type, value_flag: r2.value_flag}) AS flags,
  collect(DISTINCT {type: 'categorical', metric_type: dm3.metric_type, value_text: r3.value_text}) AS labels;
```
````

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md .claude/skills/cypher-queries/SKILL.md
git commit -m "skills: paperconfig derived_metrics_table + compartment; cypher-queries DM templates"
```

---

### Task 19: Update CLAUDE.md key-facts block

**Purpose:** `CLAUDE.md` is the first file every future conversation reads. It needs the DerivedMetric entry in its "Key graph facts" list and a mention of the new post-import computations.

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add DerivedMetric bullet to Key graph facts**

In `CLAUDE.md`, find the bullet list under `### Key graph facts` (search for `Gene↔Protein linkage` — that's the first bullet). After the `Clusteringanalysis_belongs_to_organism` / `Experiment_has_clustering_analysis` / `Gene_in_gene_cluster` bullets, add:

```markdown
- DerivedMetric nodes: 7 (Biller 2018 retrofitted evidence — 6 boolean + 1 categorical). Node IDs: `derived_metric:{doi_short}:{entry_key}:{metric_type}`. Properties: `name`, `metric_type`, `value_kind ∈ {numeric, boolean, categorical}`, `rankable`, `has_p_value`, `p_value_threshold`, `unit`, `allowed_categories`, `field_description`, `compartment`, plus denormalized parent-Experiment fields. Each DM emits exactly one of 3 measurement edge types: `Derived_metric_quantifies_gene` (numeric; ~0 edges currently, needs Waldbauer 2012 / zinser 2009 / Biller 2022 paper integrations), `Derived_metric_flags_gene` (boolean; ~4,160 edges in Biller 2018), `Derived_metric_classifies_gene` (categorical; ~258 edges in Biller 2018). Binding edges: `PublicationHasDerivedMetric`, `ExperimentHasDerivedMetric`, `DerivedMetricBelongsToOrganism`. See `docs/kg-changes/non-de-evidence-extension.md`.
- DerivedMetric computed properties (post-import): `total_gene_count` (int), `growth_phases` (str[] from parent Experiment).
- Numeric DM edge computed properties (post-import, only when parent `rankable="true"`): `rank_by_metric` (int, 1 = highest value), `metric_percentile` (float 0-100), `metric_bucket` (str: `top_decile`/`top_quartile`/`mid`/`low`, pinned thresholds ≥90/≥75/≥25). When parent `has_p_value="true"` AND threshold is set AND `adjusted_p_value` non-null: `significant` (str `"true"`/`"false"`).
- Experiment `compartment` property (str, adapter-emitted, default `"whole_cell"`). Vocab: `whole_cell`, `vesicle`, `exoproteome`, `secretome`. Computed Experiment rollups: `reports_fold_change` (str `"true"`/`"false"`), `reports_derived_metric_types` (str[]), `derived_metric_count` (int), `derived_metric_value_kinds` (str[]), `derived_metric_gene_count` (int).
- Publication/OrganismTaxon DM rollups (post-import): `derived_metric_count`, `derived_metric_gene_count`, `compartments` (str[]), `derived_metric_types` (str[]), `derived_metric_value_kinds` (str[]).
- Gene DM routing signals (post-import): `numeric_metric_count`, `classifier_flag_count`, `classifier_label_count` (ints); `numeric_metric_types_observed`, `classifier_flag_types_observed`, `classifier_label_types_observed` (str[]); `compartments_observed` (str[]).
```

- [ ] **Step 2: Extend Post-import indexes bullet**

Find the existing bullet `- Post-import indexes: scalar (...) + full-text (...)`. Inside the scalar list, insert before the closing `)`:

```
, `derived_metric_metric_type_idx`, `derived_metric_value_kind_idx`, `derived_metric_compartment_idx`, `derived_metric_omics_type_idx`, `derived_metric_treatment_type_idx`, `derived_metric_organism_idx`, `derived_metric_experiment_idx`, `experiment_compartment_idx`
```

And inside the full-text list, add `; `derivedMetricFullText` on DerivedMetric name, field_description`.

- [ ] **Step 3: Add paperconfig `derived_metrics_table` reference**

Find the `## Adding Omics Data from Publications` section. In the paperconfig supplementary_materials entry types table, add the same row as in paperconfig SKILL.md:

```markdown
| `derived_metrics_table` | Column-level scalar summaries per gene (DerivedMetric nodes + one of 3 measurement edge types per metric) |
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "CLAUDE.md: DerivedMetric + compartment + post-import rollups + derived_metrics_table entry type"
```

---

### Task 20: Final sweep — full tests + byte-identical DE regression

**Purpose:** Gate Plan 3 as complete. Every test passes; DE-only subset of post-import-validate is byte-identical against the step-0 baseline (no DE regression from new DM-scoped Cypher).

**Files:** None modified — verification only

- [ ] **Step 1: Full unit test suite**

```bash
uv run pytest -m "not slow and not kg" -v
```

Expected: 1,588+ tests pass (no regressions from Plan 2 baseline).

- [ ] **Step 2: Full KG validity suite**

```bash
uv run pytest -m kg -v
```

Expected: all pass — including new `test_derived_metric.py` (~30 tests); `test_numeric_derived_metric.py` tests SKIP cleanly (production graph has no fixture).

- [ ] **Step 3: DE-only regression diff against step-0 baseline**

The step-0 baseline (captured at the start of this slice effort, pre any DM work) is `scripts/post-import-validate.sh > step_0_baseline.txt` — should exist from Plan 1/2 or be captured fresh from the pre-Plan-3 rebuild state if not on disk.

If available:
```bash
scripts/post-import-validate.sh > /tmp/after_plan3_final.txt
# Filter to DE-only sections (excludes new DM sections + indexes list)
awk '
  /======== INDEXES ========/        {skip=1}
  /======== EXPERIMENT \(full/       {skip=0}
  /======== DERIVEDMETRIC/           {skip=1}
  /======== EXPRESSION_STATUS/       {skip=0}
  /======== END ========/            {skip=0}
  !skip
' /tmp/after_plan3_final.txt > /tmp/after_plan3_final_de_only.txt
# Apply same filter to baseline
# ... (adjust awk for baseline format)
diff /tmp/step_0_baseline_de_only.txt /tmp/after_plan3_final_de_only.txt
```

Expected: empty diff on DE-only sections (Experiment DE fields, Publication DE fields, OrganismTaxon, ClusteringAnalysis, GeneCluster, BriteCategory, Gene aggregates, rank subsamples).

- [ ] **Step 4: /omics-edge-snapshot final parity check**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save plan3_final_production
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare after_biller2018_retrofit_removal
```

Expected:
- `Changes_expression_of` total: 227,361 (unchanged from step-0 baseline)
- DerivedMetric flag edges: 0 → 4,160 (IMPROVED)
- DerivedMetric classifies edges: 0 → 258 (IMPROVED)
- DerivedMetric quantifies edges: 0 → 0 (unchanged — production has no numeric)
- No LOST / GONE publications

- [ ] **Step 5: Commit final snapshot**

```bash
git add .claude/skills/omics-edge-snapshot/snapshots/plan3_final_production.json
git commit -m "omics-edge-snapshot: plan3_final_production snapshot (production, fixture-free)"
```

- [ ] **Step 6: Mark plan complete**

Edit `docs/superpowers/plans/2026-04-20-plan3-non-de-evidence-biller-2018.md` — append a `## Completion status` section with a task→commit map (same format as Plan 2).

```bash
git add docs/superpowers/plans/2026-04-20-plan3-non-de-evidence-biller-2018.md
git commit -m "plan3: record completion status + task→commit map"
```

Plan 3 complete — slice DoD met:
- `pytest -m kg` green ✅
- `pytest -m "not slow and not kg"` green ✅
- `scripts/post-import-validate.sh` DE-only diff empty vs step-0 baseline ✅
- `/omics-edge-snapshot` reports expected DerivedMetric counts ✅
- Docs/skill files committed ✅

---

## Self-review

**Spec coverage check** against slice spec §Plan 3 (`docs/superpowers/specs/2026-04-19-non-de-evidence-biller-2018-slice.md` lines 208-226):

| Spec requirement | Plan 3 task(s) |
|---|---|
| Post-import Cypher additions (rank/percentile/bucket, significance) | Tasks 6, 7 |
| Post-import Cypher additions (analysis-node rollups) | Task 3 |
| Post-import Cypher additions (Experiment/Publication/OrganismTaxon rollups) | Tasks 4, 5 |
| Post-import Cypher additions (Gene routing counts) | Task 8 |
| Post-import Cypher additions (new scalar + full-text indexes) | Task 2 |
| Byte-identical `post-import.sh` ↔ `post-import.cypher` | Task 9 Step 1-2 |
| `post-import-validate.sh` DE-only diff empty vs baseline | Tasks 0, 9 Step 4, 20 Step 3 |
| KG validity tests per §Tests (filtered to DerivedMetric scope) | Tasks 11, 12, 13 |
| Regenerated `snapshot_data.json` | Task 16 |
| `/omics-edge-snapshot` extension (DM counts) | Task 14 |
| `docs/kg-changes/non-de-evidence-extension.md` | Task 17 |
| `.claude/skills/paperconfig/SKILL.md` update | Task 18 Step 1-2 |
| `.claude/skills/cypher-queries/SKILL.md` update | Task 18 Step 3 |
| `.claude/skills/omics-edge-snapshot/SKILL.md` update | Task 14 Step 3 |
| `CLAUDE.md` key-facts update | Task 19 |

Parent-spec §Success-criteria KG validity assertions (DerivedMetric scope) mapped to `test_derived_metric.py` / `test_numeric_derived_metric.py`:

| Assertion | Test |
|---|---|
| All 4 DM edges target :Gene | `test_measurement_edges_target_gene` (3 parametrized) |
| Each DM emits exactly one edge type matching value_kind | `test_dm_emits_only_one_edge_type` |
| Denormalized fields match parent Experiment | `test_dm_denormalized_scalar_matches_parent` (6 parametrized) + `test_dm_denormalized_experiment_id_matches_parent` + `test_dm_publication_doi_matches_parent` |
| `value_text` ∈ `allowed_categories` on classifies | `test_classify_edges_value_text_in_allowed_categories` |
| `value_flag` ∈ {'true','false'} on flags | `test_flag_edges_value_flag_enum` |
| Bucket ↔ percentile pinned thresholds | `test_bucket_matches_pinned_thresholds` + `test_rank_1_is_top_decile` |
| Significance gating + computation | `test_significant_only_on_has_p_value`, `test_significant_computed_correctly`, `test_significant_enum` |
| Rollup consistency (Experiment/Pub/Org/Gene) | 7 tests across both files |
| Empty-state defaults | `test_experiment_dm_empty_state_defaults`, `test_gene_dm_empty_state_defaults` |
| Binding edge 1:1 cardinality parent→DM | 3 tests |
| Index existence (9 new) | `test_plan3_index_exists` (parametrized) + `EXPECTED_INDEXES` update |
| `rank_by_metric` non-null iff rankable='true' | `test_rank_by_metric_only_on_rankable` |
| Rank contiguity | `test_rank_contiguous_per_dm` |
| String enums (`rankable`, `has_p_value`, `significant`, `value_flag`, `reports_fold_change`) | 5 separate tests |

No placeholders. All Cypher in tasks is complete and runnable.

## Risks + open questions at Plan 3 authoring

1. **Percentile formula edge case (N=1).** If a future DerivedMetric has only one gene, percentile = 100.0 → `top_decile`. Unavoidable but documented; a 1-gene DM is biologically uninformative anyway. Adapter-side minimum-rows validator could be added as follow-up if this becomes a real issue.

2. **BioCypher label capitalization** (binding edges CamelCase, measurement edges Caps_underscore) — already surfaced in live-graph check. Plan 3 Cypher uses exact Neo4j labels throughout. If a future BioCypher upgrade changes the convention, the Cypher blocks break loudly (test_derived_metric detects absent edges → fail rather than silent zero).

3. **Synthetic fixture wire-in/unwire dance** — mid-plan Docker rebuild with fixture, then rebuild without. The fixture DOI `10.9999/synthetic-numeric-dm` is clearly synthetic, so even if an "unwire" commit slipped, detection is trivial (`MATCH (p:Publication {doi: '10.9999/...'})` → >0 in production is a bug). Task 15 Step 4 explicitly checks this.

4. **Step-0 baseline may not exist on disk.** If Plan 1/2 didn't save `step_0_baseline.txt`, Task 20 Step 3 (DE-only byte-identical diff) can't run. Mitigation: capture baseline before Task 2 if not present. The live-graph incremental tests in Tasks 3-8 already detect most regressions; byte-identical diff is a stricter-still backstop.

5. **Per-block Cypher idempotency** — every SET statement is idempotent (re-running produces same result). CALL IN TRANSACTIONS OF N ROWS is safe to re-execute. So the per-task live-graph testing pattern (run the block against current graph, verify no error) doesn't leave the graph in a broken state even if subsequent tasks get delayed.

## Confirmed consistency with existing post-import code

- Group structure mirrored: indexes → small aggregations → heavy writes (via `CALL { … } IN TRANSACTIONS`)
- Empty-state defaults BEFORE compute (overrides) — same as existing `Experiment stats defaults`/`Experiment stats computation` pattern
- `apoc.coll.sort(apoc.coll.toSet(...))` + `coalesce(x, [])` patterns reused
- `label_in_input` → BioCypher-capitalized labels used throughout (verified live)
- String enums for boolean semantics (`'true'`/`'false'`) — matches existing `expression_status`, `has_cross_genus_members`, `level_is_best_effort`
- Timing: the 3 existing `time cypher-shell <<'CYPHER' … CYPHER` wrappers each get new statements appended; no new group added

## Follow-up refactors (out of scope for Plan 3)

- Once a paper with numeric DMs is integrated (Waldbauer 2012 / zinser 2009 / Biller 2022), regenerate `snapshot_data.json` so `Derived_metric_quantifies_gene` gets real sample edges in the snapshot (currently empty in production).
- AbundanceAnalysis slice — parent spec §1, triggered by Biller 2022 vesicles or Oleza 2015/2017 exoproteome. Separate plan.
- MCP tool surface (`gene_derived_metrics`, `derived_metric_ranked_genes`, `publication_derived_metrics`) — MCP team's follow-up, unblocked by this slice's comms doc.
