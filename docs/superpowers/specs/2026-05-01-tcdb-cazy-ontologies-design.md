# TCDB and CAZy Ontologies — Design Spec

**Date:** 2026-05-01
**Status:** Approved for plan-writing
**Scope:** Add Transporter Classification Database (TCDB) and CAZy carbohydrate-active-enzyme classification as first-class ontologies in the KG, surfaced via the existing MCP ontology tools. Bridge TCDB substrate annotations into the chemistry layer via Metabolite nodes.

## Motivation

Genes already carry flat array properties `transporter_classification: str[]` (TCDB IDs from eggNOG `KEGG_TC`) and `cazy_ids: str[]` (CAZy families from eggNOG `CAZy`). These are search-only — no graph traversal, no hierarchy navigation, no chemistry linkage.

Goals:

1. Promote both classification systems to graph nodes with hierarchical edges so MCP can surface them as ontologies (parallel to GO, EC, Pfam, BRITE, KEGG, COG).
2. Bridge TCDB to the chemistry layer: substrate annotations on TCDB families become `Metabolite` nodes (reusing existing where possible) connected by a new `Tcdb_family_transports_metabolite` edge.
3. Mirror the Pfam-normalization precedent: drop the flat Gene properties once the graph nodes/edges replace them.

## Non-goals

- Reaction-style modelling of transport (no `TransportReaction` nodes; no per-event stoichiometry; no direction). TCDB is treated as a hierarchical *classifier* of transport activity, structurally analogous to EC numbers — see "Why TCDB is an ontology, not a reaction" below.
- Full CAZy reference download (v1 stays observed-only — every CAZy family in the KG comes from a real eggNOG annotation across our configured strains).
- CAZy → EC linkage. Out of scope for v1.
- New MCP tool code. Existing ontology tools auto-pick-up the new node labels.

## Why TCDB is an ontology, not a reaction

TCDB has both flavours:

- *Ontology-like*: hierarchical classifier (5 levels), stable IDs, gene-annotation target, controlled vocabulary published as a curated scheme.
- *Reaction-like*: describes a biochemical process (something moves from somewhere to somewhere), has substrates.

The closest existing parallel is **EC numbers**: hierarchical classifier of catalytic activity, gene-annotation target, no chemistry directly attached — chemistry lives on separate `Reaction` nodes that *carry* EC numbers. EC alone has no per-reaction stoichiometry.

TCDB is structurally closer to EC than to Reaction: substrates are annotated **at the family level**, not per-reaction. A family says "transports {Na+, K+}" — no per-event stoichiometry, no direction (consistent with our spec 1.2 "no reaction direction"). There's no canonical "transport reaction" object in TCDB — only the family classification + the substrate set.

Treatment: TCDB is an ontology. Substrate links are *annotation edges from the ontology to chemistry* (`Tcdb_family_transports_metabolite`), not reaction edges. Genes attach via `Gene_has_tcdb_family` (membership verb, mirrors `Gene_has_pfam` / `Gene_in_cog_category` rather than `Gene_catalyzes_*`). No `TransportReaction` nodes. No fake stoichiometry.

If a future need arises for transport-reaction-level modelling (e.g. pairing with FBA, modelling co-transport stoichiometry of H+/Na+ symporters), `TransportReaction` nodes can be added later as a separate layer that *uses* `TcdbFamily` for classification — exactly how `Reaction` uses `EcNumber`.

## Schema additions

### Two new node labels

| YAML label | Neo4j label | Preferred ID prefix | Levels |
|---|---|---|---|
| `tcdb family` | `TcdbFamily` | `tcdb` (e.g. `tcdb:3.A.1.5.2`) | 5: `tc_class`, `tc_subclass`, `tc_family`, `tc_subfamily`, `tc_specificity` |
| `cazy family` | `CazyFamily` | `cazy` (e.g. `cazy:GH13`) | up to 3: `cazy_class`, `cazy_family`, `cazy_subfamily` |

Both follow the BRITE pattern: single label per ontology, distinguished by `level_kind` (str) + `level` (int, 0 = broadest = class). Level is set at adapter emit time — no post-import inference needed (mirrors the unified ontology level convention from `docs/kg-changes/ontology-level.md`).

```yaml
tcdb family:
  is_a: named thing
  represented_as: node
  preferred_id: tcdb
  properties:
    name: str                  # TCDB entry name; falls back to tcdb_id when source name is empty (typical for subclass/subfamily/specificity)
    tcdb_id: str               # bare "1.A.1.5.2" form for display
    level: int                 # 0=tc_class … 4=tc_specificity
    level_kind: str
    superfamily: str           # sparse, leaf-only
    tc_class_id: str           # post-import: pointer to root tc_class node ID (e.g. "tcdb:3"); self on tc_class nodes
    gene_count: int            # post-import
    organism_count: int        # post-import
    member_count: int          # post-import: direct child count
    metabolite_count: int      # post-import: distinct metabolites reachable in subtree

cazy family:
  is_a: named thing
  represented_as: node
  preferred_id: cazy
  properties:
    name: str                  # human-readable for class (e.g. "Glycoside Hydrolases"); falls back to cazy_id for family/subfamily
    cazy_id: str               # e.g. "GH13" or "GH13_5"
    level: int                 # 0=cazy_class, 1=cazy_family, 2=cazy_subfamily
    level_kind: str
    gene_count: int            # post-import
    organism_count: int        # post-import
```

### Five new edge labels

| YAML label | Neo4j label | Source → Target | Purpose |
|---|---|---|---|
| `gene has tcdb family` | `Gene_has_tcdb_family` | Gene → TcdbFamily | Gene annotation, attaches at exact level annotated |
| `tcdb family is a tcdb family` | `Tcdb_family_is_a_tcdb_family` | TcdbFamily → TcdbFamily | Parent edge within hierarchy |
| `tcdb family transports metabolite` | `Tcdb_family_transports_metabolite` | TcdbFamily → Metabolite | Substrate linkage; only on leaves with substrates |
| `gene has cazy family` | `Gene_has_cazy_family` | Gene → CazyFamily | Gene annotation, attaches at exact level annotated |
| `cazy family is a cazy family` | `Cazy_family_is_a_cazy_family` | CazyFamily → CazyFamily | Parent edge within hierarchy |

### Removed Gene properties

- `transporter_classification: str[]`
- `cazy_ids: str[]`

Mirrors Pfam normalization precedent. Query goes through the new graph edges. `gene_summary` / `geneFullText` already cover ID search needs.

## Substrate handling — TCDB → Metabolite

### Resolution rule

Every TCDB leaf substrate (format `CHEBI:NNNN;name`) is resolved through the MNX SQLite at step-6 build time:

1. Look up `chebi:NNNN` in `compound_aliases` → MNXM ID.
2. Derive primary Metabolite ID from the MNXM, following the existing `Metabolite` adapter rule:
   - If MNXM has a `kegg.compound:C*` alias → primary ID `kegg.compound:C*`.
   - Else if MNXM resolves to a ChEBI alias → primary ID `chebi:NNNN`.
   - Else → primary ID `mnx:MNXM*`.
3. Look up the primary ID in the existing pruned-KEGG compound set:
   - If present → reuse the existing Metabolite. Add `transport` to its `evidence_sources`.
   - If absent → create a new entry in `kegg_data.json` (under either `compounds` for `kegg.compound:C*` IDs, or a new `additional_compounds` section for non-KEGG IDs). Properties pulled from MNX: `name`, `formula`, `mass`, `inchikey`, `mnxm_id`, `chebi_id`. `evidence_sources: ['transport']`.

All 1,410 unique TCDB substrate strings resolve to MNXM (verified empirically). ~1,184 (84%) have structure (formula/mass/inchikey); 335 (16%) are MNX ontology classes (e.g. "molecule", "ion", "drug") that lack structure. Per user decision: **all** of them become Metabolite nodes — no slug-IDs needed because every substrate has a real CHEBI ID.

### Substrate collection rule

Genes annotate at any TCDB level. To collect substrates for pruning and edge emission:

- For each gene's annotated TCDB ID, walk **down** to all leaf descendants (`tc_specificity` level) in the TCDB hierarchy and union their substrate sets. A family-level annotation (e.g. `3.A.1`) thereby represents "could move any of these substrates."
- For pruning, also walk **up** from each annotated ID to `tc_class` so the hierarchy reaches the root and is queryable as ontology levels.

### Compound cache schema extension

`cache/data/kegg/kegg_data.json` (built by step 6) gains:

- New section `additional_compounds: {chebi_id_or_mnx_id: {...}}` for compounds without KEGG aliases.
- New per-entry property `evidence_sources: list[str]` on every entry in both `compounds` and `additional_compounds`. Values from `{"metabolism", "transport"}`. Compounds reachable via the existing KEGG/Reaction path get `metabolism`. Compounds reached via the TCDB substrate path get `transport`. A compound found by both gets both. Pre-existing entries' lists are unioned.

Adapter consequence: `metabolism_adapter.py` reads both sections, emits one Metabolite node per entry. New property `evidence_sources: str[]` lands on every Metabolite node.

### `Tcdb_family_transports_metabolite` edges

Emitted by the TCDB adapter only on leaf nodes (`tc_specificity`), one edge per (leaf, resolved Metabolite primary ID). The pre-resolved (leaf → metabolite primary IDs) mapping lives in `tcdb_pruned.json` (built by step 6) so the adapter does not re-resolve at build time.

## Pruning rule

### TCDB

Step 6 prunes the full ~13,642-entry TCDB hierarchy to the subtree reachable from gene annotations. Specifically:

- Seed set: every TCDB ID in any strain's `transporter_classification`.
- Walk **down** from each seed to all leaf descendants (so substrate linkage is complete).
- Walk **up** from each seed to `tc_class` (so hierarchy reaches root).
- Kept = union of seeds ∪ down-descendants ∪ up-ancestors.
- Persisted as `cache/data/tcdb/tcdb_pruned.json`. Adapter reads this file, never re-prunes.

Mirrors the BRITE pruning precedent (`docs/superpowers/specs/2026-04-14-brite-pruning-design.md`) but bidirectional (BRITE was bottom-up only).

### CAZy

No pruning step — `CazyFamily` nodes are inferred at adapter init from observed eggNOG annotations across configured strains. By construction every emitted node has at least one gene-annotation source.

## Build pipeline

### Pipeline ownership reorganization

| Script | Before | After |
|---|---|---|
| `scripts/refresh_mnx.sh` | Builds MNX SQLite + `tcdb_hierarchy.json` + `cazy_hierarchy.json` | MNX SQLite only |
| `prepare_data.sh` step 0 sub-step 6 | Downloads TCDB reference TSVs | Removed |
| `prepare_data.sh` step 6 (`build_kegg_metabolism_xrefs.py`) | Builds pruned KEGG cache | Builds pruned KEGG cache + TCDB download + TCDB hierarchy parse + TCDB pruning + substrate→metabolite resolution + `additional_compounds` fold-in |

`build_metabolite_resolver.py` is renamed to `build_mnx_resolver.py`. After this spec it does only one thing: build the heavy MNX SQLite (`metabolite_resolver.db`). All other artifacts move to step 6 or are deleted.

### Step 6 detailed flow (extended)

1. (existing) Walk every strain's `gene_annotations_merged.json` → collect gene-reachable {KOs, reactions, KEGG compounds, KEGG pathways}.
2. (new) Same walk → collect every TCDB ID into `known_tcdb_ids`, every CAZy ID into `known_cazy_ids` (CAZy collection is informational only — adapter handles its own pruning).
3. (new) Download TCDB reference TSVs via `download_metabolism_reference.download_all(sources=["tcdb"])` if missing or `--force`.
4. (new) Parse TCDB TSVs (logic moved from `build_metabolite_resolver.py`) → `cache/data/tcdb/tcdb_hierarchy.json` (raw, all ~13K entries, kept as a debugging artifact + `tcdb_utils.py` data source).
5. (new) Prune TCDB to gene-reachable (above + below). Write `cache/data/tcdb/tcdb_pruned.json` containing kept node set + per-leaf resolved metabolite primary IDs.
6. (new) For each kept leaf, resolve substrates via MNX → fold into `kegg_data.json` (`compounds` or `additional_compounds`), union `evidence_sources` lists.
7. (existing, extended) Write `kegg_data.json` with the `additional_compounds` section and `evidence_sources` tags.

### CAZy hierarchy: pure-Python data table

`cazy_hierarchy.json` is deleted. The hardcoded class table moves into `multiomics_kg/utils/cazy_utils.py`:

```python
CAZY_CLASSES = {
    "GH":  "Glycoside Hydrolases",
    "GT":  "GlycosylTransferases",
    "PL":  "Polysaccharide Lyases",
    "CE":  "Carbohydrate Esterases",
    "AA":  "Auxiliary Activities",
    "CBM": "Carbohydrate-Binding Modules",
}

_CAZY_FAMILY_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)(?:_(\d+))?$")

def parse_cazy_id(token: str) -> tuple[str, str | None] | None: ...
def is_valid_cazy(value: str) -> bool: ...
def cazy_ancestors(cazy_id: str) -> list[str]: ...
```

Pure-Python, no file I/O. The one external consumer (`annotation_transforms.py:265,275`) is also being removed in this spec — but `is_valid_cazy` stays in the utility module for tests / future use.

### Validation hooks dropped from gene annotation pipeline

Both TCDB and CAZy become passthrough fields, mirroring `kegg_reactions`:

- `config/gene_annotations_config.yaml`:
  - Remove `transform: validate_tcdb` from `transporter_classification`.
  - Remove `transform: validate_cazy` from `cazy_ids`.
- `multiomics_kg/download/utils/annotation_transforms.py`:
  - Drop `_tx_validate_tcdb`, `_tx_validate_cazy` functions.
  - Drop `"validate_tcdb"`, `"validate_cazy"` from the transform-name registry.
  - Drop `is_valid_tcdb`, `is_valid_cazy` imports.

`gene_annotations_merged.json` will then carry whatever eggNOG emitted. Unknown IDs are silently dropped at adapter / pruning time (with debug logs).

## Adapters

### File placement

Two new files, mirroring `brite_adapter.py`:

- `multiomics_kg/adapters/tcdb_adapter.py`
- `multiomics_kg/adapters/cazy_adapter.py`

### Two-class shape (mirrors Pfam adapter pattern)

Each module exposes:

- A **per-strain class** that loads `gene_annotations_merged.json`, exposes `get_all_*_ids()` for orchestrator pruning, and `get_edges()` yielding gene→ontology edges.
- A **multi-strain orchestrator** that owns ontology nodes + hierarchy edges + (TCDB only) substrate edges, and delegates per-strain gene-edge emission.

Naming convention: `XAnnotationAdapter` / `MultiXAnnotationAdapter` — matches Pfam, GO, EC, KEGG.

```python
class TcdbAnnotationAdapter:
    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None: ...
    def get_all_tcdb_ids(self) -> set[str]: ...
    def get_edges(self): ...   # yields Gene_has_tcdb_family

class MultiTcdbAnnotationAdapter:
    def __init__(self, genome_config_file, cache_root, test_mode, cache): ...
    def download_data(self, cache: bool = True) -> None: ...   # reads tcdb_pruned.json
    def get_nodes(self): ...   # TcdbFamily nodes
    def get_edges(self): ...   # Tcdb_family_is_a_tcdb_family + Gene_has_tcdb_family + Tcdb_family_transports_metabolite

class CazyAnnotationAdapter:
    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None: ...
    def get_all_cazy_ids(self) -> set[str]: ...
    def get_edges(self): ...   # yields Gene_has_cazy_family

class MultiCazyAnnotationAdapter:
    def __init__(self, genome_config_file, test_mode): ...
    def download_data(self, cache: bool = True) -> None: ...   # no-op (in-memory derivation)
    def get_nodes(self): ...   # CazyFamily nodes (class + family + subfamily)
    def get_edges(self): ...   # Cazy_family_is_a_cazy_family + Gene_has_cazy_family
```

### Wiring

Two new adapter blocks in `create_knowledge_graph.py`, slotted next to `brite_adapter`. Both depend on the metabolism adapter's `kegg_data.json` being complete (so transport-only Metabolites are present before TCDB tries to emit edges to them) — but since step 6 runs in `prepare_data.sh` (not in `create_knowledge_graph.py`), no in-pipeline ordering change is needed.

### Coverage estimate (rough, scaled from MED4)

- ~91/1976 MED4 genes have TCDB → ~1,500–2,500 `Gene_has_tcdb_family` edges across 25 strains.
- ~785 new transport-only Metabolite nodes (current chemistry layer ~2,188 → ~2,973).
- ~10K `Tcdb_family_transports_metabolite` edges from leaf substrate annotations.
- ~23 MED4 genes have CAZy → ~400–600 `Gene_has_cazy_family` edges across all strains.
- 6 `CazyFamily` class nodes + ~30–50 family/subfamily nodes (small ontology — Pro/Alteromonas have very little carbohydrate metabolism).

## Post-import additions

In `scripts/post-import.sh` and reference `scripts/post-import.cypher`:

### Indexes

- Scalar: `tcdb_family_level_idx`, `tcdb_family_level_kind_idx`, `tcdb_family_tcdb_id_idx`, `cazy_family_level_idx`, `cazy_family_level_kind_idx`, `cazy_family_cazy_id_idx`.
- Full-text: `tcdbFamilyFullText` on `name`, `tcdb_id`, `superfamily`; `cazyFamilyFullText` on `name`, `cazy_id`.

### Computed properties on TcdbFamily / CazyFamily

- `gene_count` (int): `COUNT(DISTINCT gene)` reachable via `Gene_has_*` edges. For non-leaf nodes, count includes genes attached anywhere in the subtree (descendant traversal).
- `organism_count` (int): same but distinct organisms.
- `member_count` (int, TCDB only): direct child count in the hierarchy.
- `metabolite_count` (int, TCDB only — per explorer ask TCDB-S4): distinct `Metabolite` nodes reachable via `Tcdb_family_transports_metabolite` edges in the subtree (mirrors the `gene_count` subtree-traversal pattern). On leaves equals direct count; on `tc_class` shows total substrate breadth across all descendants.
- `tc_class_id` (str | null, TCDB only — per explorer ask TCDB-S5): pre-computed pointer to the root `tc_class` ID for this node (e.g. `"tcdb:3"` for any node under class 3). Self-reference on `tc_class`-level nodes. Null only if a node somehow has no class ancestor (shouldn't happen). Avoids variable-length traversal in class-level filter queries. Cheap to populate since the pruning walk already touches every ancestor.

### Gene routing signals

Extend the existing `Gene.annotation_types` array post-import update to add `'tcdb'` / `'cazy'` when the gene has any `Gene_has_tcdb_family` / `Gene_has_cazy_family` edge.

Per explorer asks TCDB-S1 / TCDB-S2, also populate quantity-bearing rollups (mirror existing `expression_edge_count` / `numeric_metric_count` etc.):

- `Gene.tcdb_family_count` (int, default 0): `COUNT(*)` of `Gene_has_tcdb_family` edges per gene. Single-hop.
- `Gene.cazy_family_count` (int, default 0): `COUNT(*)` of `Gene_has_cazy_family` edges per gene. Single-hop.

### Coordination with chemistry-slice-1: `Gene.metabolite_count`

The chemistry-slice-1 KG-side asks (companion doc `multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md`, ask KG-A2 + cross-spec coordination point TCDB-S3) introduce `Gene.metabolite_count: int` defined as **distinct metabolites reachable via *any* gene-reaching path** — UNION of catalysis and transport.

Sequencing rule:

- If chemistry-slice-1 lands first, this spec extends its post-import block to add the transport arm — same property, additional UNION clause for `Gene → TcdbFamily → ... → Tcdb_family_transports_metabolite → Metabolite`.
- If this spec lands first, the post-import block ships with both arms baked in (catalysis + transport). The catalysis arm is identical to KG-A2's slice-1 form; the transport arm is the TCDB extension.

Either way the property's semantics is **UNION-across-paths**. No second property like `Gene.transport_metabolite_count` is introduced.

### Metabolite extended properties

- New: `transporter_count` (int) — distinct `TcdbFamily` nodes with `Tcdb_family_transports_metabolite` edges to this metabolite.
- Existing `gene_count` / `organism_count`: extend the materialization to UNION the transport path (`Gene → TcdbFamily → Metabolite`) with the existing catalysis path (`Gene → Reaction → Metabolite`).

### `evidence_sources` enum semantics

This spec defines two values: `"metabolism"` (compound participates in a Reaction) and `"transport"` (compound annotated as a TCDB family substrate). The enum is **open-ended** — explorer ask MET-M4 reserves `"metabolomics"` for a future metabolomics-DM spec (compound measured in a metabolomics experiment). Adapter / post-import logic must use set-union semantics so future contributors can add their own value without disturbing existing tags.

### Organism_has_metabolite extension

Existing post-import block creates these via `Gene → Reaction → Metabolite`. Add a UNION arm via `Gene → TcdbFamily → Metabolite`. Single MERGE deduplicates by `(organism, metabolite)`. Optional `via: str[]` property on the edge enumerating `{'metabolism', 'transport'}`.

## MCP surface

No MCP code changes. The existing ontology tools auto-pick-up the new node labels through Neo4j's schema introspection:

- `mcp__multiomics-kg__list_filter_values` and `kg_schema` will surface `TcdbFamily` and `CazyFamily` once they exist in the graph.
- `genes_by_ontology` works against any node label that has `Gene_has_*` edges + `name` + `level` — both new ontologies satisfy that contract.
- `gene_ontology_terms` returns all ontology nodes reachable from a gene — both will appear naturally.
- `ontology_landscape` aggregates by `level_kind` — works out of the box.
- `search_ontology` uses the `name` full-text indexes — `tcdbFamilyFullText` and `cazyFamilyFullText` plug in.

Side effect: the new substrate-attached Metabolite nodes (~785 transport-only) will show up in `gene_overview` / `genes_by_function` results. Their `evidence_sources: ['transport']` lets MCP-side rendering distinguish "this compound is annotated as a transport substrate" from "this compound is catalyzed by a reaction in this organism."

## Documentation updates

- `CLAUDE.md`:
  - Extend "Key Adapters" with `tcdb_adapter.py` and `cazy_adapter.py`.
  - Add `TcdbFamily` and `CazyFamily` to the node label list.
  - Add the five new edge labels.
  - Update prepare_data step 6 description.
  - Update mention of `build_metabolite_resolver.py` → `build_mnx_resolver.py`.
- `memory/MEMORY.md`: one-liner under "Recently Shipped" pointing to a new memory file `project_tcdb_cazy_ontologies.md`.
- `.claude/skills/cypher-queries/SKILL.md`: add `TcdbFamily`, `CazyFamily` to the Node Labels table; add the five new edge types; add 2–3 query templates.
- New doc `docs/kg-changes/tcdb-cazy-ontologies.md`: full integration narrative (data sources, schema additions, pruning rules, substrate-resolution rule, `evidence_sources` semantics, TCDB-as-ontology-not-reaction rationale).

## Tests

### Unit tests (`tests/`)

- `test_tcdb_adapter.py` — node emission per level; parent edge correctness; gene-edge attachment at exact level; transports-metabolite edges only on leaves; pruning correctness (above + below).
- `test_cazy_adapter.py` — class/family/subfamily emission; hardcoded class table coverage; parent edges; gene edges; malformed-token logging.
- `test_cazy_utils.py` — rewrite to test the pure-function API (no JSON fixture).
- `test_build_kegg_metabolism_xrefs.py` — extend with TCDB walk, substrate resolution, `additional_compounds` section emission, `evidence_sources` tagging logic.
- `test_build_metabolite_resolver.py` (renamed `test_build_mnx_resolver.py`) — drop `test_build_cazy_hierarchy`, drop TCDB-related tests, drop `cazy_hierarchy_entry_count` assertion.
- `test_prepare_data_step2_metabolism_smoke.py` — drop `cazy_hierarchy.json` path entry.
- `test_build_gene_annotations_metabolism.py`, `test_annotation_transforms_metabolism.py` — drop `cazy_hierarchy.json` bootstrap fixtures.
- `test_annotation_transforms.py` (or equivalent) — drop tests for `validate_tcdb` / `validate_cazy` transforms.

### KG validity tests (`tests/kg_validity/`)

- `test_tcdb_cazy.py` (new) — node count sanity ranges; gene-edge presence (no dangling); full-text indexes exist; level/level_kind invariants; orphan detection (every non-class-level node has a parent edge).
- `test_post_import.py` — extend to verify TCDB/CAZy `gene_count` / `organism_count` populated and consistent.
- `test_structure.py` — extend orphan checks to cover `TcdbFamily`, `CazyFamily`.
- `tests/kg_validity/snapshot_data.json` — regenerate after rebuild; sample a TCDB hierarchy chain, a CAZy family, a transport-only Metabolite with `evidence_sources: ['transport']`.

### Edge-count regression

`/omics-edge-snapshot` before and after — these changes don't touch `Changes_expression_of` so the snapshot must match exactly. Catches accidental regressions in the omics layer from schema-config edits.

## Migration sequencing — commit boundaries

Each commit is independently reviewable and leaves the build runnable.

1. **CAZy: pure-Python utils + drop validation + drop `cazy_hierarchy.json`.**
   - Rewrite `multiomics_kg/utils/cazy_utils.py` to pure-Python (regex + `CAZY_CLASSES` table).
   - Remove `_tx_validate_cazy`, `"validate_cazy"`, `is_valid_cazy` import from `annotation_transforms.py`.
   - Remove `transform: validate_cazy` from `gene_annotations_config.yaml`.
   - Drop CAZy logic from `build_metabolite_resolver.py` (functions, path entry, force-check, report counter).
   - Delete `cache/data/cazy/cazy_hierarchy.json`.
   - Update tests (drop fixtures, rewrite `test_cazy_utils.py`).
   - Verification: unit tests pass.

2. **TCDB: move download + parse to step 6, rename build script, drop validation.**
   - Move TCDB download from `download_genome_data.py` step 0 sub-step 6 → call into `download_metabolism_reference.download_all(sources=["tcdb"])` from `build_kegg_metabolism_xrefs.py`.
   - Move TCDB hierarchy parse from `build_metabolite_resolver.py` → step 6.
   - Rename `build_metabolite_resolver.py` → `build_mnx_resolver.py`. Update `scripts/refresh_mnx.sh` reference.
   - Drop `step6_tcdb_reference()` from `download_genome_data.py`; update arg help.
   - Update `scripts/prepare_data.sh` step 0 sub-step list.
   - Remove `_tx_validate_tcdb`, `"validate_tcdb"`, `is_valid_tcdb` import from `annotation_transforms.py`.
   - Remove `transform: validate_tcdb` from `gene_annotations_config.yaml`.
   - Verification: prepare_data step 6 still produces a valid `kegg_data.json` (no schema change yet); unit tests pass.

3. **Step 6: substrate resolution + `additional_compounds` + `evidence_sources` tagging.**
   - Add `evidence_sources: str[]` property on `metabolite` in `schema_config.yaml` (so the new property has a schema home before adapter starts emitting it).
   - Extend `build_kegg_metabolism_xrefs.py` to walk TCDB substrates, resolve via MNX, fold into `kegg_data.json` (new `additional_compounds` section).
   - Tag every compound entry with `evidence_sources: list[str]`.
   - Write `cache/data/tcdb/tcdb_pruned.json` (pruned kept-node set + per-leaf resolved metabolite primary IDs).
   - `metabolism_adapter.py` reads `additional_compounds` and emits Metabolite nodes (sets `evidence_sources` property).
   - No new node types yet — TcdbFamily/CazyFamily schema changes deferred to commit 4.
   - Update KG-validity Metabolite-count thresholds to accept the +785 transport-only growth (test_structure.py and any snapshot assertions on Metabolite cardinality).
   - Verification: `/omics-edge-snapshot` unchanged; Metabolite count grows by ~785; rebuild + KG validity pass.

4. **Schema additions in `schema_config.yaml` (TcdbFamily / CazyFamily).**
   - Add `tcdb family`, `cazy family` node labels and the five new edge labels.
   - Verification: schema validates against empty input (no adapter changes yet).

5. **TCDB adapter + wire-in.**
   - Create `multiomics_kg/adapters/tcdb_adapter.py` with `TcdbAnnotationAdapter` + `MultiTcdbAnnotationAdapter`.
   - Wire into `create_knowledge_graph.py`.
   - Drop `transporter_classification` from `Gene.properties` in `schema_config.yaml`.
   - Drop the `transporter_classification` field from `gene_annotations_config.yaml` (or keep for caching but don't surface to Gene — depends on whether other code reads it; default: drop).
   - Verification: unit tests pass; full build + KG validity pass; `Gene_has_tcdb_family` + `Tcdb_family_is_a_tcdb_family` + `Tcdb_family_transports_metabolite` edge counts within sanity ranges.

6. **CAZy adapter + wire-in.**
   - Create `multiomics_kg/adapters/cazy_adapter.py` with `CazyAnnotationAdapter` + `MultiCazyAnnotationAdapter`.
   - Wire into `create_knowledge_graph.py`.
   - Drop `cazy_ids` from `Gene.properties` in `schema_config.yaml`.
   - Drop the `cazy_ids` field from `gene_annotations_config.yaml` (same caveat).
   - Verification: unit tests pass; full build + KG validity pass.

7. **Post-import additions.**
   - Indexes (scalar + full-text).
   - `gene_count` / `organism_count` / `member_count` / `metabolite_count` computation for `TcdbFamily`; `gene_count` / `organism_count` for `CazyFamily`.
   - `TcdbFamily.tc_class_id` sparse pointer (per TCDB-S5).
   - Gene `annotation_types` extension to include `'tcdb'` / `'cazy'`.
   - `Gene.tcdb_family_count` / `Gene.cazy_family_count` rollups (per TCDB-S1 / TCDB-S2).
   - `Gene.metabolite_count` — extend (or initialize, if chemistry-slice-1 hasn't landed) with the UNION transport arm (per TCDB-S3 coordination).
   - Metabolite `transporter_count` + extended `gene_count` / `organism_count` (UNION transport path).
   - Extended `Organism_has_metabolite` materialization (UNION transport path).
   - Run `scripts/post-import-validate.sh` before/after to verify deterministic dump diff.
   - Verification: `pytest tests/kg_validity/ -v` passes.

8. **Snapshot regeneration + new validity tests + docs.**
   - Regenerate `tests/kg_validity/snapshot_data.json`.
   - Add `tests/kg_validity/test_tcdb_cazy.py`.
   - Update `CLAUDE.md`, `memory/MEMORY.md`, `.claude/skills/cypher-queries/SKILL.md`.
   - Add `docs/kg-changes/tcdb-cazy-ontologies.md`.
   - Verification: full test suite passes; `/omics-edge-snapshot` unchanged from commit 3.

## Open questions / deferred items

- **Arabinose / fructose / hexose substrate rescue**: ~30–50 named compounds resolve to MNX ontology classes (no structure) when looked up by ChEBI alone; their structured stereoisomer-specific MNXM lives one hop away. A name-lookup fallback in step 6 could rescue them, but adds complexity. **Deferred** — log them at build time, decide in a follow-up if it matters.
- **CAZy enrichment to full reference + EC linkage**: out of scope for v1. Reconsider if a use case emerges.
- **`TransportReaction` reified entities**: explicitly out of scope — see "Why TCDB is an ontology, not a reaction." Revisit if FBA / co-transport stoichiometry use cases emerge.
- **Renaming `build_metabolite_resolver.py` → `build_mnx_resolver.py`**: included in this spec (commit 2). Trivial, clarifies the script's purpose now that it does only one thing.

## Cross-spec coordination notes

- **`gene_details` MCP tool may surface removed Gene properties** — the explorer team flagged that `gene_details` uses a `g{.*}` pattern that would expose `transporter_classification` / `cazy_ids` raw. After commits 5 / 6 land, those properties are gone — the tool's output will simply omit them. No tool code change is required for this spec, but the explorer team will verify on the explorer side. Tracked in `multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md` under "Coordination note: removed `Gene.transporter_classification` array property".
- **Chemistry-slice-1 KG-side asks** (`KG-A1`, `KG-A2`, `KG-A3`, `KG-A4`) — adjacent chemistry-layer rollups owned by the explorer team's chemistry-slice-1 PR. Not in this spec's scope, except the `Gene.metabolite_count` UNION coordination above (TCDB-S3 / KG-A2). Both PRs commute; either can land first.
- **`evidence_sources` value `"metabolomics"`** — reserved for future metabolomics-DM spec per MET-M4. Documented in the enum-semantics section above.
