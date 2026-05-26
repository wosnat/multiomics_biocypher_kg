# Ontology patterns — the existing-ontology inventory, invariants & copy targets

The KG already carries **11 ontology node
types** (the explorer MCP exposes them as a fixed **12-value `ontology` enum** —
GO splits into `go_bp`/`go_mf`/`go_cc`, BRITE is scoped per-tree). Every one was
built with the same adapter → schema → post-import triad. A new ontology must
match this shape — **copy a template, don't invent.**

Line anchors verified against the repo on 2026-05-26; they drift, so they are
hedged with `~` and you should re-`grep` the symbol, not trust the number blindly.

## Shared invariants — a new ontology MUST satisfy all of these

- **Node id** is a CURIE with a per-ontology prefix (`cazy:`, `pfam:`, `pfam.clan:`,
  `tcdb:`, `kegg.brite:`, `go:`, …). Use `normalize_curie()` from
  `multiomics_kg/utils/curie_utils.py`.
- **`level: int`** on every node, **0 = root**, increasing with specificity. Set by
  the **adapter** — a `_classify`-style helper ([cazy_adapter.py:51](../../../../multiomics_kg/adapters/cazy_adapter.py#L51)),
  a pre-built cache entry (tcdb), or emit-time depth (brite/cyanorak) — **except GO**,
  whose `level` is computed post-import. For a flat tool, hardcode `level = 0`.
  The MCP routes every ontology by `level`, so it is mandatory.
- **Adapter-emitted `name`** (+ optional `short_name`); post-import adds `gene_count`,
  `organism_count`, and `member_count` (hierarchical only). See post-import-patterns.md.
- **`level_kind: str`** is OPTIONAL — set it on hierarchical (multi-level) ontologies,
  omit on flat. No MCP consumer filters on `level_kind` (only on `level`).
- **Hierarchy edge** `<x>_is_a_<x>`, **child → parent**, property-less (`{}`).
- **Gene→ontology edge** `gene_has_<x>` (a few use verbs: `gene_catalyzes_ec_number`,
  `gene_in_cog_category`, `gene_involved_in_biological_process`). **No gene→ontology
  edge in the graph today carries any edge properties — all yield `{}`** (confirmed
  across all 11). The *only* Gene-incident edge with rich properties is
  `Changes_expression_of` ([schema_config.yaml:~415](../../../../config/schema_config.yaml#L415):
  `log2_fold_change`, `adjusted_p_value`; post-import `rank_by_effect`/`rank_up`/
  `rank_down`). **A scored tool (psortb, signalp) will be the first scored ontology
  edge — borrow the edge-property + rank pattern from `Changes_expression_of`.**
- **Adapter split:** per-strain `<Tool>Adapter` reads the merged field and emits
  gene→ontology edges; `Multi<Tool>Adapter` **owns** the nodes + hierarchy edges and
  is the only class wired into `create_knowledge_graph.py` (instantiate → optional
  `download_data()` → `write_nodes()` → `write_edges()`). See
  [cazy_adapter.py:120](../../../../multiomics_kg/adapters/cazy_adapter.py#L120) +
  [create_knowledge_graph.py:233](../../../../create_knowledge_graph.py#L233).
- **Multi-call fan-out:** loop over the list-valued merged field, emit one edge per
  element, all converging on shared, de-duplicated nodes. 1:1 tools emit ≤1 edge.
- **`_clean_str`** (`'`→`^`, drop `|`) on every string node/edge property.

## Node-property contract (3A)

| Property | Type | Requirement |
|---|---|---|
| `name` | str | display; `_clean_str`-sanitized |
| `<tool>_id` (e.g. `cazy_id`) | str | bare identifier / CURIE local part |
| `level` | int | **0 = root**; adapter-emitted; mandatory |
| `level_kind` | str | optional — hierarchical only |
| `gene_count` | int | post-import |
| `organism_count` | int | post-import |
| `member_count` | int | post-import, hierarchical only |

## The merged-annotation field shape

- **Bare categorical (no score):** simple `str[]` list, like `cazy_ids`
  ([gene_annotations_config.yaml:540](../../../../config/gene_annotations_config.yaml#L540), `type: union`).
- **Categorical with per-call score(s):** a list of structured records (the
  "structured passthrough" shape). The merge writes well-formed JSON; the adapter
  reads it with `json.load`, so the BioCypher `|`/`'` CSV constraint does **not**
  apply to the merged file (only when the adapter *yields* strings to BioCypher).
  - **Genuinely multi-call** → a structured list-of-dicts via `passthrough` of the
    calls.json object.
  - **1:1** (psortb, signalp) → **parallel scalar fields** are simpler
    (`<tool>_call: str` + `<tool>_score: float`); no list.

## Copy targets

| New-tool shape | Copy | Why |
|---|---|---|
| Flat, multi-call, no score | `cazy_adapter` (observed-only, self-contained, `_classify`) | the canonical skeleton |
| Flat, single-call, **scored** (psortb, signalp) | `cazy_adapter` node/Multi skeleton **+** `Changes_expression_of` edge-property/rank pattern | no existing flat+scored ontology — compose the two |
| Hierarchical, pre-built cache | `tcdb_adapter` (5-level, `tcdb_pruned.json`, `seed_aliases`) or KEGG-KO 4-level | pruned-cache hierarchy |
| Hierarchical, computed closure | `MultiGoAnnotationAdapter` + `go_utils.compute_ancestry_closure` | DAG ancestry closure |
| Pruned-to-reachable | `brite_adapter` (`known_ko_ids` bottom-up mark) | reachable-subset pruning |

## Footguns found across the 11 — warn the operator

- **Edge-label casing is not uniform** (`gene_has_tcdb_family` vs
  `gene_catalyzes_ec_number` vs `gene_in_cog_category`) — read `schema_config.yaml`
  per-ontology; don't pattern-match the label name.
- **Node ownership splits differently:** most `Multi*` own nodes while the per-strain
  class emits gene edges, but `CogRoleAnnotationAdapter` emits COG/Cyanorak/TIGR edges
  with **no nodes** (Multi owns all three vocabularies), and **BRITE has no per-strain
  class at all**.
- **`level`/`level_kind` provenance varies** (adapter `_classify` / pre-built cache /
  emit-time depth / hardcoded-0 / post-import for GO) — pick the one matching your
  template, don't mix.

## What the MCP reads (the cross-node rollups Step 4 owns)

The ontology **query** tools (`ontology_landscape`, `search_ontology`,
`genes_by_ontology`, …) count genes/organisms **live** by traversing `Gene_has_<x>`
edges; they need only `level` + `name` + `is_uninformative` on the new node and
**nothing** on Gene/Organism/Publication. Registering the new ontology in the
explorer's `ONTOLOGY_CONFIG` + the 12-value enum is an **explorer-repo change, out
of scope**.

The **gene-routing / overview tools, however, read folded-in Gene rollups** — and
*that* is the cross-node work this skill's post-import owns (functional ontologies
only). See post-import-patterns.md → "Cross-node rollups" and decision-tree.md →
"Decision 2 — functional-vs-structural".

`OrganismTaxon` and `Publication` carry **no per-ontology rollup** — a new ontology
needs nothing on those.
