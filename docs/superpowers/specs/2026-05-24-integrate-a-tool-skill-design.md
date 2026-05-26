# `integrate-a-tool` skill — design spec

**Topic:** A new user-invocable Claude Code skill that performs the *second half* of per-strain
bioinformatics-tool integration — taking a Phase-1 tool's committed `<strain>.<tool>.calls.json`
artifacts and wiring them into the knowledge graph.

**Date:** 2026-05-24
**Status:** design approved; ready for implementation plan
**Sibling skill:** [`add-a-tool`](../../../.claude/skills/add-a-tool/SKILL.md) (Phase 1 — produces the calls.json)

---

## 1. Context & purpose

`add-a-tool` (Phase 1) scaffolds a `/<tool>-run` skill that runs an external predictor per strain
and writes inspectable, git-committed `<strain>.<tool>.calls.json` artifacts under
`cache/data/<org>/genomes/<strain>/<tool>/`. It explicitly **defers** all KG coupling to a
"Phase 2 hand-off" that, until now, had no skill behind it.

`integrate-a-tool` (Phase 2) is that missing skill. Invoked for a tool that already has Phase-1
artifacts, it walks the operator from `calls.json` to a validated, documented presence in the live
Neo4j graph: merged into `gene_annotations_merged.json`, surfaced as either a new ontology
(nodes + edges) or new Gene properties, rolled up + indexed in post-import, validated against a live
rebuild, and documented in a release-notes / what-changed doc. **Making the result queryable through the
explorer MCP is a separate, out-of-scope follow-up (§9).**

> **Phase-1 artifacts may predate this convention.** The per-tool run skills that exist (`psortb-run`,
> `tcdb-diamond`, `signalp-run`, `eggnog-run`) were written across the `add-a-tool` transition, so they are
> not all consistent with its calls.json contract — e.g. `psortb-run` and `tcdb-diamond` emit the standard
> `<strain>.<tool>.calls.json`, but `signalp-run` still writes raw `signalp/output.json`. (There is no `cazy`
> run-skill — CAZy annotations are eggNOG-derived.) Step 1 must inspect the actual artifact and, when it
> doesn't follow the calls.json convention, either normalize it (a fix that lives in the `/<tool>-run` skill)
> or parse the raw output in `load_<tool>()`. Treat this as the common case, not the exception.

### Design decisions (locked during brainstorming)

| Decision | Choice | Rationale |
|---|---|---|
| Skill name | `integrate-a-tool` | Reads as the Phase-2 partner to `add-a-tool`; triggers on "integrate X into the KG", "wire X into the graph", "Phase 2 for X". |
| Up-front design | **Full spec + plan first** | Graph-shape changes are high-stakes; matches the project's 5-phase workflow and the fact that every shipped tool already has a spec. |
| Attach target | **Gene-only (+ new ontology)** | All four existing tools (psortb/tcdb/cazy/signalp) are protein→Gene. Protein/Organism/region targets are an explicit future extension. |
| Definition of "done" | **Through live KG validation** | Final gate is a full Docker rebuild + post-import + `/omics-edge-snapshot` diff + `pytest -m kg` + snapshot regeneration. Not done until the new graph elements are verified in live Neo4j. |

### What is NOT changing (scope boundaries)

- **Phase 1 is untouched** — `integrate-a-tool` consumes calls.json; it never re-runs the external tool.
- **No refactor of `build_gene_annotations.py` into a fully config-driven merge.** The loader +
  `_get_raw`/`build_wide`/`build_merged` signatures are hand-wired for the existing sources; this skill
  *adds* a source following that hand-wired pattern. Generalizing the builder is an explicit non-goal.
- **Gene-only attach target.** Protein/Organism/per-region/per-genome-scalar attach targets are out of scope.
- **`integrate-a-tool` does not author a new `/<tool>-run` skill** — that is `add-a-tool`'s job. The durable
  artifacts here are the per-tool KG-integration spec, the code/schema/post-import edits, and the
  `docs/kg-changes/<tool>-extension.md` what-changed doc.

---

## 2. The canonical data path

Both integration tracks flow through the **same front door**: the tool's calls.json becomes a new
*source* in the gene-annotation merge, producing a field (or fields) in `gene_annotations_merged.json` —
exactly as eggNOG's output yields the `cazy_ids` / `transporter_classification` fields that
`cazy_adapter` / `tcdb_adapter` already read today.

```
<strain>.<tool>.calls.json                 (Phase 1 artifact, committed)
        │   new `sources:` + `fields:` entry in config/gene_annotations_config.yaml
        │   new load_<tool>() + builder edits in build_gene_annotations.py
        ▼
gene_annotations_merged.json[locus_tag]    field(s): list of calls, optionally with per-call score
        │
        ├── ONTOLOGY track (3A): <tool>_adapter.py reads the field → nodes + Gene_has_<X> edges
        │                        (edge carries score/tier as properties when present)
        │
        └── GENE-PROPERTY track (3B): schema_config gene property → cyanorak_ncbi_adapter auto-reads it
```

Merging first (rather than letting an adapter read calls.json directly) is **mandatory** because it is
what makes `contributing_sources` and the auto-generated `DataSource` node pick the tool up for free
([data_source_adapter.py](../../../multiomics_kg/adapters/data_source_adapter.py) derives everything from
the `logical_sources` block).

---

## 3. The categorical-vs-property decision tree (the crux)

The operator inspects the distinct values the tool emits across all strains (`jq` over the committed
calls.json) and routes:

```
├─ Small controlled vocabulary (≲ few hundred distinct values), shared across genes,
│  with hierarchy or cross-gene grouping semantics?
│        e.g. PSORTb localization {Cytoplasmic, OuterMembrane, …},
│             SignalP type {SP, LIPO, TAT, PILIN}, a TC/CAZy-style class tree
│   → ONTOLOGY TRACK (3A): new node type + Gene_has_<X> edges
│     (+ <X>_is_a_<X> hierarchy edges if there are levels). Copy the cazy_adapter skeleton.
│
├─ The categorical call ALSO carries a numeric score / probability / confidence / tier?
│        e.g. PSORTb localization + score, SignalP type + probability,
│             tcdb-diamond TC family + identity/qcov + tier
│   → STILL ONTOLOGY TRACK (3A); the score(s) ride as EDGE PROPERTIES on Gene_has_<X>
│     (the node de-duplicates the controlled vocabulary; per-gene evidence strength lives on the edge,
│      exactly like Changes_expression_of carries log2_fold_change + adjusted_p_value).
│
└─ Per-gene numeric / scalar / continuous value, OR free-text, OR a single boolean flag,
   with NO shared cross-gene vocabulary?
        e.g. a standalone probability, a count, GC%, a yes/no with no taxonomy
   → GENE-PROPERTY TRACK (3B): new property on the Gene node in schema_config;
     cyanorak_ncbi_adapter auto-reads it from gene_annotations_merged.json.
```

**Multiple calls per gene is the default for the ontology track.** A gene may carry many calls (one
protein in several TC families, several CAZy families, several localizations). The merged field is a
**list**; the adapter emits **one `Gene_has_<X>` edge per call**, fanning a single gene out to **N edges
pointing at N shared, de-duplicated nodes**. Nodes are owned once by the `Multi<Tool>Adapter`; many
genes converge on one node; one gene fans out to many.

---

## 3.5 The existing-ontology inventory — copy targets & invariants

The KG already carries **11 ontology node types** (the MCP exposes them as a fixed **12-value `ontology`
enum** — GO splits into `go_bp`/`go_mf`/`go_cc`, BRITE is scoped per-tree). Every one was built with the
same adapter → schema → post-import triad. A new ontology must match this shape — copy a template, don't
invent. (Operational detail — invariants, copy targets, footguns — lives in the skill's
[`references/ontology-patterns.md`](../../../.claude/skills/integrate-a-tool/references/ontology-patterns.md);
this section is the distilled rule set it was built from.)

**Shared invariants — a new ontology MUST satisfy all of these:**

- **Node id** is a CURIE with a per-ontology prefix (`cazy:`, `pfam:`, `pfam.clan:`, `tcdb:`, `kegg.brite:`,
  `go:`, …).
- **`level: int`** on every node, **0 = root**, increasing with specificity. Set by the adapter — a
  `_classify`-style helper (cazy), a pre-built cache entry (tcdb), or emit-time depth (brite/cyanorak) —
  **except GO**, whose `level` is computed post-import. For a flat tool, hardcode `level = 0`. The MCP
  routes every ontology by `level`, so it is mandatory (§4.2). `level_kind` is optional (§4.2).
- **Adapter-emitted `name`** (+ optional `short_name`); post-import adds `gene_count`, `organism_count`,
  and `member_count` (hierarchical only).
- **Hierarchy edge** `<x>_is_a_<x>`, **child → parent**, property-less.
- **Gene→ontology edge** `gene_has_<x>` (a few use verbs: `gene_catalyzes_ec_number`, `gene_in_cog_category`,
  `gene_involved_in_biological_process`). **No gene→ontology edge in the graph today carries any edge
  properties — all yield `{}`** (confirmed across all 11). The *only* Gene-incident edge with rich
  properties is `Changes_expression_of` (`schema_config.yaml` ~L406: `log2_fold_change`, `adjusted_p_value`,
  post-import `rank_by_effect`/`rank_up`/`rank_down`). **A scored tool (psortb, signalp) will be the first
  scored ontology edge — borrow the edge-property + rank pattern from `Changes_expression_of`, not from any
  ontology.**
- **Adapter split:** per-strain `<Tool>Adapter` reads the merged field and emits gene→ontology edges;
  `Multi<Tool>Adapter` **owns** the nodes + hierarchy edges and is the only class wired into
  `create_knowledge_graph.py` (instantiate → optional `download_data()` → `write_nodes()` → `write_edges()`).
- **Multi-call fan-out:** loop over the list-valued merged field, emit one edge per element, all converging
  on shared, de-duplicated nodes (one gene → N edges → N shared nodes). 1:1 tools simply emit ≤1 edge.
- **`_clean_str`** (`'`→`^`, drop `|`) on every string node/edge property.

**Copy targets:**

| New-tool shape | Copy | Why |
|---|---|---|
| Flat, multi-call, no score | `cazy_adapter` (observed-only, self-contained, `_classify`) | the canonical skeleton |
| Flat, single-call, **scored** (psortb, signalp) | `cazy_adapter` node/Multi skeleton **+** `Changes_expression_of` edge-property/rank pattern | no existing flat+scored ontology — compose the two (see §8) |
| Hierarchical, pre-built cache | `tcdb_adapter` (5-level, `tcdb_pruned.json`, `seed_aliases`) or KEGG-KO 4-level | pruned-cache hierarchy |
| Hierarchical, computed closure | `MultiGoAnnotationAdapter` + `go_utils.compute_ancestry_closure` | DAG ancestry closure |
| Pruned-to-reachable | `brite_adapter` (`known_ko_ids` bottom-up mark) | reachable-subset pruning |

**Inconsistencies to warn the operator about** (real footguns found across the 11):

- **Edge-label casing is not uniform** (`gene_has_tcdb_family` vs `gene_catalyzes_ec_number` vs
  `gene_in_cog_category`) — read `schema_config.yaml` per-ontology; don't pattern-match the label name.
- **Node ownership splits differently:** most `Multi*` own nodes while the per-strain class emits gene
  edges, but `CogRoleAnnotationAdapter` emits COG/Cyanorak/TIGR edges with **no nodes** (Multi owns all
  three vocabularies), and **BRITE has no per-strain class at all**.
- **`level`/`level_kind` provenance varies** (adapter `_classify` / pre-built cache / emit-time depth /
  hardcoded-0 / post-import for GO) — pick the one matching your template, don't mix.

---

## 3.6 What the MCP reads — the cross-node rollups the post-import must feed

The ontology **query** tools (`ontology_landscape`, `search_ontology`, `genes_by_ontology`,
`gene_ontology_terms`, `pathway_enrichment`, `cluster_enrichment`) count genes/organisms **live by
traversing `Gene_has_<x>` edges** (`genome_coverage` = live `count(DISTINCT g)` ÷ organism gene_count). They
need only `level` + `name` + `is_uninformative` on the new node and **nothing** on Gene/Organism/Publication.
A new ontology registers in the explorer's `ONTOLOGY_CONFIG` (`label` + `gene_rel` + `hierarchy_rels` (may be
`[]`) + `fulltext_index`) — but that is an explorer-repo change, **out of scope** (§9). The new node's own
`gene_count`/`organism_count` (§4.2) are for non-MCP/UI use + post-import sanity.

The **gene-routing / overview tools, however, read folded-in Gene rollups** — and *that* is the cross-node
work this skill's post-import owns:

| MCP surface | Stored Gene props it reads | Consequence for a new ontology |
|---|---|---|
| `gene_overview` | `annotation_types`, `informative_annotation_types` (whole arrays), `annotation_state`, `annotation_quality`, curated `*_count` signals (expression / ortholog / cluster / metric / `reaction_count` / `metabolite_count` / `tcdb_family_count`→`transporter_count`) — `queries_lib.py:505-524` | folding the source into the arrays surfaces it **automatically** (whole-array return — no MCP change) |
| `genes_by_function` | filters `g.annotation_quality >= min_quality` (`queries_lib.py:272`); envelope `by_annotation_state` + flattened `annotation_types` (`:408-414`) | if the source enters the `annotation_quality` bucket count it **shifts which genes pass `min_quality`** |
| `gene_details` | full `g{.*}` dump (`tools.py:1859`) | any new Gene property (a `<tool>_count`, a denormalized scalar) appears **for free** |

`OrganismTaxon` and `Publication` carry **no per-ontology rollup** (`list_organisms` surfaces only
gene/publication/experiment counts + treatment/omics + chemistry `reaction_count`/`metabolite_count`, which
are metabolism-layer specials, not a generic-ontology pattern). A new ontology needs **nothing** on those.

**The functional-vs-structural decision (Step 0) drives whether you touch any of this:**

- **Functional-annotation ontology** — a new functional *evidence* source (KO-like / domain-like / role-like):
  **fold the `Gene_has_<x>` edge into `Gene.annotation_types` + `Gene.informative_annotation_types`**
  (`post-import.cypher:609` / `:634`) so `gene_overview` reflects it. Then decide whether it is an
  *informative* source that also enters the **`annotation_state`/`annotation_quality` 8-bucket count**
  (`:548-576`) — the bucket list is *enumerated*, not auto-discovered, so adding one follows the maintenance
  procedure in `docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md` (both `.cypher` +
  `.sh`, CLAUDE.md, the bucket-count test) and **changes `genes_by_function` `min_quality` results**.
- **Structural / non-functional ontology** — psortb localization, signalp signal-peptide (*where the gene is*,
  not *what it does*): **do NOT** fold into `annotation_types`/`annotation_quality` (it would inflate
  functional-annotation coverage and skew `min_quality`). Optionally add a dedicated `<tool>_count` (or a
  denormalized scalar) Gene signal — it surfaces via `gene_details` for free; adding it to `gene_overview`'s
  curated signal list is an MCP-side change (out of scope, §9).

> **Correction to a tempting shortcut:** the ontology query tools traversing live does **not** mean
> "no cross-node rollup needed." The gene-routing tools (`gene_overview`, `genes_by_function`,
> `gene_details`) are the consumers of the folded-in Gene rollups — they are the reason Step 4 touches Gene.

---

## 4. Shapes & invariants

### 4.1 The merged-annotation field (clean JSON, not pipe-joined)

- **Bare categorical (no score):** a simple `str[]` list, like `cazy_ids`.
  ```json
  "cazy_ids": ["GH13", "GH13_5"]
  ```
- **Categorical with per-call score(s):** a **list of structured records** (the "structured passthrough"
  shape). The merge writes well-formed JSON; the adapter reads it with `json.load`, so the BioCypher
  `|`/`'` CSV constraint does **not** apply to the merged file itself (it applies only when the adapter
  *yields* string properties to BioCypher).
  ```json
  "<tool>_calls": [
    {"call": "FamilyA", "score": 9.5},
    {"call": "FamilyB", "score": 7.2}
  ]
  ```
  How the score travels through the merge is a per-tool Step-0 decision driven by **cardinality**:
  - **Genuinely multi-call** (one protein → several scored calls) → a structured list-of-dicts via
    `passthrough` of the calls.json object (preferred). The existing `union`/`single`/`passthrough_list`
    rules are string/list-oriented; carrying a list-of-dicts uses `passthrough` of an already-structured value.
  - **1:1** (one call per protein) — e.g. psortb, signalp (§8) → **parallel scalar fields** are simpler
    (`<tool>_call: str` + `<tool>_score: float`); no list needed.

### 4.2 The ontology node (3A)

Every new ontology node **must** carry, adapter-emitted:

| Property | Type | Requirement |
|---|---|---|
| `name` | str | human-readable display; `_clean_str`-sanitized |
| `<tool>_id` (e.g. `cazy_id`) | str | the bare identifier / CURIE local part |
| `level` | int | **0 = broadest/root**; increases with specificity. Mandatory — the MCP `ontology_landscape` / `search_ontology` tools rely on it. |
| `level_kind` | str | e.g. `psortb_class`; **optional** — set it on hierarchical (multi-level) ontologies; flat single-level ontologies may omit it. No MCP consumer filters on `level_kind` (only on `level`), so omitting it is safe. |
| `gene_count` | int | post-import |
| `organism_count` | int | post-import |
| `member_count` | int | post-import, only if hierarchical |

`level` / `level_kind` are set by the adapter from the hierarchy depth (see `_classify` in
`cazy_adapter`), **not** by post-import. Flat (single-level) ontologies set `level = 0` for every node
and may omit `level_kind`.

### 4.3 The `Gene_has_<X>` edge (3A)

- **Bare categorical:** property-less (`{}`), like `Gene_has_cazy_family`.
- **With score:** a real `properties:` block in `schema_config`, e.g. `score: float`, `probability: float`,
  `confidence_tier: str`. The adapter yields the per-call value in the edge-property dict. Floats stay
  numeric (no `_clean_str`); only string sub-fields (a tier label) are sanitized.

### 4.4 The Gene property (3B)

Two hand-edits are required — the read is **not** automatic over arbitrary merged fields:

1. Declare the property under `gene:` `properties:` in `schema_config.yaml`.
2. Add a matching member to the hand-maintained `GeneNodeField` enum in `cyanorak_ncbi_adapter.py`.
   `set_node_fields` only emits fields enumerated there; the `GeneEnumMeta` metaclass does **not** load
   from `schema_config`. **Numeric** properties additionally need their field name added to the hardcoded
   `int_fields` / `float_fields` sets in `_get_gene_nodes` — otherwise the value falls through the string
   branch (`unquote` → `clean_text` → `_split_field`) and lands as a string.

Type mapping: numeric → `float`/`int` (+ the `int_fields`/`float_fields` edit); flag → `str`
("true"/"false") per existing convention; free-text → `str`.

---

## 5. Skill identity & file layout

**Frontmatter**
- `name: integrate-a-tool`
- `user-invocable: true`
- `argument-hint: <tool-name> [--spec <path>]`
- `allowed-tools`: Read, Edit, Write, Grep, Glob, Bash(uv/python/pytest/docker/git/mkdir/bash/jq/diff),
  plus the `run_cypher` / `omics-edge-snapshot` invocation surface used in Step 5.
- `description`: triggers on "integrate signalp into the KG", "wire psortb into the graph", "add the X
  tool's predictions to the knowledge graph", "Phase 2 for <tool>", "turn <tool> calls.json into
  nodes/properties".

**`add-a-tool` change:** its "Phase 2 hand-off" section gets a one-line redirect: *"→ use
`/integrate-a-tool`."*

**Files the skill produces / edits per tool**

```
docs/superpowers/specs/<date>-<tool>-kg-integration-design.md   # Step 0 (spec-first)
docs/superpowers/plans/<date>-<tool>-kg-integration.md          # checkbox plan
config/gene_annotations_config.yaml                             # +source +fields +logical_sources        (Step 2)
multiomics_kg/download/utils/annotation_transforms.py           # +_tx_<name> + register in _TRANSFORMS   (Step 2, if needed)
multiomics_kg/download/build_gene_annotations.py                # +load_<tool>(), extend _get_raw/build_wide/build_merged sigs + row-fetch + contributing_sources (Step 2)
multiomics_kg/utils/<tool>.py                                   # pure parse/vocab helpers (reuse if exists)  (Steps 1,3)
multiomics_kg/adapters/<tool>_adapter.py                        # ONTOLOGY track only — copy cazy_adapter   (Step 3A)
config/schema_config.yaml                                       # +node +edges (3A)  OR  +gene property (3B)
multiomics_kg/adapters/cyanorak_ncbi_adapter.py                 # +GeneNodeField enum member (+int_fields/float_fields for numerics)  (Step 3B only)
create_knowledge_graph.py                                       # wire Multi<Tool>Adapter                   (3A only)
scripts/post-import.sh + scripts/post-import.cypher             # rollups + indexes, kept byte-identical     (Step 4)
tests/test_<tool>.py , tests/test_annotation_transforms.py      # unit tests                                 (Steps 2,3)
tests/kg_validity/test_<tool>.py (or extend existing)           # live-graph assertions                      (Step 5)
docs/kg-changes/<tool>-extension.md                             # release-notes / what-changed doc           (Step 6)
CLAUDE.md                                                       # Key graph facts + Actual Neo4j labels      (Step 6)
```

```
.claude/skills/integrate-a-tool/
├── SKILL.md
├── references/
│   ├── decision-tree.md              # expanded tree + worked PSORTb/SignalP examples
│   ├── ontology-patterns.md          # the 11-ontology inventory + invariants + copy-target table (§3.5), line-anchored
│   ├── gene-annotation-merge-recipe.md   # exact build_gene_annotations.py + transform edit points, line-anchored
│   └── post-import-patterns.md       # rollup + index snippets for both tracks, real cazy/tcdb examples
└── assets/
    ├── ontology_adapter_template.py          # per-strain + Multi* skeleton: level/level_kind, multi-call fan-out, edge-property payload
    ├── schema_node_edge_snippet.yaml         # node (+ level) + edge (+ optional properties block)
    ├── gene_annotations_source_snippet.yaml  # sources: + fields: + logical_sources: + transform stub
    ├── post_import_rollup_snippet.cypher     # gene_count/organism_count (*0.. hierarchical + flat), index block, annotation_types fold-in, optional edge-score rank
    └── what_changed_template.md              # docs/kg-changes/ template
```

---

## 6. The 7-step workflow

**Step 0 — Capture intent + write spec/plan.** Four bullets — (1) what the tool predicts and its value
space; (2) calls.json shape + key (WP_ vs locus_tag) + nulls/sentinels; (3) the categorical-vs-property
verdict from the decision tree (and, if categorical, whether a score rides on the edge); (4) the target
node label / edge label(s) / id prefix / hierarchy levels OR the target Gene property name + type. Expand
into `docs/superpowers/specs/<date>-<tool>-kg-integration-design.md` + a checkbox plan. The spec states the
merged field name(s), the merged-field shape (§4.1), and any new transform needed.

**Step 1 — Inspect the Phase-1 runner & calls.json.** Read the tool's `/<tool>-run` SKILL.md; `jq` the
committed calls.json to enumerate distinct values, the key, nulls/sentinels, and whether each call has an
associated score. **Reuse existing parse utils** — check `multiomics_kg/utils/<tool>.py` first. The
generic `multiomics_kg/utils/tool_calls_io.py` may not exist yet (as of 2026-05 each runner inlines its
own calls.json I/O; `add-a-tool` ships `assets/tool_calls_io_template.py` to copy on first use) — if
calls.json I/O isn't already centralized there, copy that template. Write new parsing only if absent, and
keep it pure (unit-testable, no filesystem).

**Step 2 — Wire the tool into the gene-annotation merge** (both tracks):
- `gene_annotations_config.yaml`: add a `sources:` entry (`type`, `path_pattern` to the calls.json,
  `join_to`, `logical_sources:` with `id`/`scope`/`provenance`) + one or more `fields:` rules.
- **If the raw value needs reshaping** (prefix, split, normalize, list-of-dicts assembly) → add a
  `_tx_<name>` to [annotation_transforms.py](../../../multiomics_kg/download/utils/annotation_transforms.py),
  register it in `_TRANSFORMS`, reference it via `transform:` in the field rule, and unit-test it in
  `tests/test_annotation_transforms.py`.
- `build_gene_annotations.py` (the hand-wired part — flagged explicitly): add `load_<tool>()`, extend
  `_get_raw` + `build_wide`/`build_merged` signatures + the `process_strain` row-fetch + the
  `_compute_contributing_sources` presence check.
- **Verify:** `bash scripts/prepare_data.sh --steps 2 --strains MED4 --force`, then `jq` that the field
  landed in `gene_annotations_merged.json`, that `contributing_sources` lists the tool, and that a new
  `DataSource` node will be emitted.

**Step 3 — Fork to the matching track.**
- **3A Ontology track:** vocab/hierarchy helper in `multiomics_kg/utils/<tool>.py` (shape of `cazy_utils`)
  → `schema_config.yaml` node (with mandatory `level`/`level_kind`, §4.2) + `gene_has_<x>` edge (with a
  `properties:` block iff a score rides, §4.3) + `<x>_is_a_<x>` hierarchy edge if levelled →
  `<tool>_adapter.py` with `<Tool>Adapter` (per-strain, reads the merged field, multi-call fan-out) +
  `Multi<Tool>Adapter` (owns nodes + parent edges), copied from `cazy_adapter`, using `_clean_str` on
  strings → wire in `create_knowledge_graph.py` → `tests/test_<tool>.py`.
- **3B Gene-property track:** add the property under `gene:` in `schema_config.yaml` (§4.4) **and** add the
  matching `GeneNodeField` enum member in `cyanorak_ncbi_adapter.py` (plus the `int_fields`/`float_fields`
  set for numeric properties, §4.4); then confirm the field lands on a built Gene node by spot-checking;
  unit-test the parse/transform.

**Step 4 — Post-import rollups + indexes** (edit **both** `post-import.sh` and `post-import.cypher`; keep
byte-identical):
- 3A: `gene_count` / `organism_count` (+ `member_count` if hierarchical) on the new node via the
  `MATCH … CALL { … } IN TRANSACTIONS` pattern (the cazy `*0..` subtree walk for hierarchical, no `*0..`
  for flat). **Cross-node rollups are conditional on the functional-vs-structural verdict (§3.6):**
  *functional* ontologies fold the new edge into Gene `annotation_types` + `informative_annotation_types`
  (`post-import.cypher:609`/`:634`) and decide on the `annotation_state`/`annotation_quality` 8-bucket count
  (`:548-576`, maintenance procedure applies — shifts `genes_by_function` `min_quality`); *structural*
  ontologies (psortb, signalp) **stay out** of those and at most add a dedicated `<tool>_count` Gene signal
  (surfaces via `gene_details`; `gene_overview`'s curated list is MCP-side/out-of-scope). `OrganismTaxon`/
  `Publication` need no per-ontology rollup. Scalar indexes on `level` (+ `level_kind` if hierarchical) + the
  id, and a full-text index on `name`. If a score rides on the edge, optionally compute a per-gene
  `rank_by_<score>` (mirrors expression `rank_by_effect`).
- 3B: scalar index on the new Gene property if it will be filtered; optional Gene rollup.
- Capture `scripts/post-import-validate.sh > baseline.txt` **before** rebuild; **after**, confirm
  additions-only (old properties unchanged, new properties present). (Not byte-identical — this is an
  addition, not a pure refactor — so the diff should show only the new lines.)

**Step 5 — Live KG validation gate** (definition of "done"):
`/omics-edge-snapshot` (before) → Docker rebuild (`docker compose up -d --build`, background, ~1h) →
post-import → `/omics-edge-snapshot` (after, diff: expression edges unchanged, new edges appear) → add
kg-validity assertions (node count > 0, edge presence, rollup sanity, no-orphan, numeric edge-property
range checks à la `test_expression.py`) → `pytest -m kg` → regenerate `snapshot_data.json` via
`uv run python tests/kg_validity/generate_snapshot.py`.

**Step 6 — Release-notes / what-changed doc + CLAUDE.md.** Write `docs/kg-changes/<tool>-extension.md` from the template
(new node/edge types + counts, property-additions table, indexes, vocabulary, conventions, validation
results, out-of-scope) → update CLAUDE.md "Key graph facts" + "Actual Neo4j labels (BioCypher PascalCase
output)" + link the doc with `See docs/kg-changes/<tool>-extension.md`.

**Step 7 — Register / hand-off.** Add the redirect line in `add-a-tool`'s "Phase 2 hand-off" section →
`/integrate-a-tool`; confirm `add-a-strain` step-10 ("loop back through prepare_data") covers the new
merged field so future strains pick the tool up automatically.

---

## 7. Reference checklist (final pass before merging)

The SKILL.md ends with this checklist (operator creates a TodoWrite item per line):

- [ ] Step 0 spec under `docs/superpowers/specs/<date>-<tool>-kg-integration-design.md` + checkbox plan
- [ ] calls.json inspected; distinct values + key + score-presence enumerated; existing parse utils reused
- [ ] Categorical-vs-property verdict recorded in the spec, with the §4.1 merged-field shape chosen
- [ ] `gene_annotations_config.yaml`: `sources:` + `fields:` + `logical_sources:` added
- [ ] New `_tx_<name>` added + registered in `_TRANSFORMS` + unit-tested (if reshaping was needed)
- [ ] `build_gene_annotations.py`: `load_<tool>()` + `_get_raw`/`build_wide`/`build_merged` + row-fetch + `contributing_sources`
- [ ] MED4 rebuild verifies the field in `gene_annotations_merged.json` + `contributing_sources` + DataSource node
- [ ] **3A:** schema node carries `level` (0=root); `level_kind` set when hierarchical (optional for flat); edge has a `properties:` block iff a score rides
- [ ] **3A:** `<tool>_adapter.py` (per-strain + Multi*) emits multi-call fan-out (one edge per call), `_clean_str` on strings, wired into `create_knowledge_graph.py`
- [ ] **3B:** Gene property declared in schema **and** added to the `GeneNodeField` enum (+ `int_fields`/`float_fields` for numerics); value confirmed on a built Gene node
- [ ] `tests/test_<tool>.py` (pure parsing/vocab) passes
- [ ] post-import.sh + post-import.cypher edited identically (`diff` clean); `post-import-validate.sh` baseline/after shows additions-only
- [ ] Scalar indexes (`level`, `level_kind` if hierarchical, id) + full-text (`name`) for the new node
- [ ] Cross-node rollups per §3.6: **functional** → Gene `annotation_types`/`informative_annotation_types` folded in (+ `annotation_quality` bucket decision, maintenance procedure); **structural** → kept out, optional `<tool>_count` only. `OrganismTaxon`/`Publication`: nothing
- [ ] `/omics-edge-snapshot` before/after: expression edges unchanged, new edges appear
- [ ] kg-validity assertions added; `pytest -m kg` passes; `snapshot_data.json` regenerated
- [ ] `docs/kg-changes/<tool>-extension.md` written; CLAUDE.md "Key graph facts" + "Actual Neo4j labels" updated + linked
- [ ] `add-a-tool` Phase-2 hand-off redirects here; `add-a-strain` step-10 covers the new field
- [ ] `pytest -m "not slow and not kg"` passes

---

## 8. Worked examples (implementation-ready)

PSORTb and SignalP are the **immediate targets** — both Phase-1 runners exist; integrate each as soon as its
calls land. Both resolve to **3A flat ontology + edge score**: the §3 decision tree routes a small shared
vocabulary that carries a per-call score to the ontology track with the score on the edge. Both are **1:1**
(one call per protein → ≤1 edge per gene, no fan-out) and both **skip the no-signal sentinel** (no node/edge
for `Unknown`/`Other`; the *absence* of an edge already encodes "no prediction").

> **Why ontology, not a Gene property?** Localization and SP-type are shared, closed vocabularies; a node
> makes "all OuterMembrane genes in organism X" a one-hop `MATCH` and de-duplicates the vocabulary. The 1:1
> cardinality doesn't change that. (3B Gene property, §8.3, is right only for a per-gene scalar with *no*
> shared vocabulary.)

| Tool | Track | Reference skeleton |
|---|---|---|
| `cazy` | 3A hierarchical, no score | the canonical property-less skeleton to copy |
| `tcdb` | 3A hierarchical, no score | 5-level + substrate bridge; pre-built-cache template |
| `psortb` (§8.1) | 3A flat + edge score | cazy node skeleton + `Changes_expression_of` edge pattern |
| `signalp` (§8.2) | 3A flat + edge score | same, + a calls.json-normalization prerequisite |

### 8.1 PSORTb → `SubcellularLocalization`

- **Phase-1 input:** `cache/data/<org>/genomes/<strain>/psortb/<strain>.psortb.calls.json` — a dict **keyed by
  WP_ accession**; record `{localization, score (2.0–10.0), is_unknown, secondary_localization, secondary_score}`.
  `is_multi_localized` is always false → strictly 1:1.
- **Vocabulary (5 nodes, skip `Unknown`):** `Cytoplasmic`, `CytoplasmicMembrane`, `Periplasmic`,
  `OuterMembrane`, `Extracellular`.
- **Node:** label `SubcellularLocalization`, id `psortb:<Class>` (e.g. `psortb:OuterMembrane`); props `name`
  (display), `psortb_id`, `level = 0` (flat — omit `level_kind`).
- **Edge:** `Gene_has_subcellular_localization` (Gene → node), **`properties: { score: float }`** — the first
  scored ontology edge; model the schema `properties:` block on `Changes_expression_of`.
- **Merge (Step 2):** `join_to: protein_id` (WP_, exactly like eggNOG). Two scalar fields via `passthrough`:
  `psortb_localization: str`, `psortb_score: float`; skip rows where `is_unknown`. `contributing_sources`
  picks up `psortb`; the `DataSource` node auto-emits.
- **Adapter (Step 3A):** `SubcellularLocalizationAdapter` (per-strain — reads the two fields, emits ≤1 scored
  edge/gene) + `MultiSubcellularLocalizationAdapter` (owns the 5 nodes at `level=0`). Copy `cazy_adapter`,
  then **delete the fan-out loop and the `*_is_a_*` hierarchy**, and **add `score` to the edge tuple**.
  `_clean_str` on `name`/`psortb_id`.
- **Post-import (Step 4):** flat `gene_count`/`organism_count` rollup (no `*0..`); scalar indexes on `level` +
  `psortb_id`, full-text on `name`; **optional** per-class `rank_by_score` (rank genes within each
  localization by `score` desc — the `rank_by_effect` UNWIND pattern). Do **not** fold into
  `annotation_types`/`annotation_quality` (structural, not functional evidence); optionally denormalize a
  `subcellular_localization` routing string onto Gene.
- **kg-validity (Step 5):** 5 nodes exist; no `Unknown` node; edge `score ∈ [2.0, 10.0]` (numeric-range
  assertion à la `test_expression.py`); no orphan gene→node edges; rollup sanity.
- **Release notes (Step 6):** `docs/kg-changes/psortb-extension.md`.

### 8.2 SignalP → `SignalPeptideType`  ⚠️ calls.json prerequisite

- **Phase-1 input — fix the artifact first:** `signalp-run` predates the calls.json convention and writes raw
  `signalp/output.json` (envelope + `.SEQUENCES` keyed by the **full FASTA header**), *not* a normalized
  `<strain>.signalp.calls.json`. **Prerequisite (Step 0 decision):** either extend `signalp-run` to emit the
  standard calls.json, or have `load_signalp()` parse `output.json` directly.
- **Per-record:** `Prediction` (winning type), `Likelihood[6]` (∈[0,1], parallel to `Protein_types[6]`),
  `CS_pos` (e.g. `"Cleavage site between pos. 26 and 27. Probability 0.862107"`, `""` when none).
- **Key:** WP_ accession = first whitespace token of the header → reuse the existing `_tx_first_token_space`
  transform.
- **Vocabulary (5 nodes, skip `Other`):** `SP` (Sec/SPI), `LIPO` (Sec/SPII), `TAT` (Tat/SPI), `TATLIPO`
  (Tat/SPII), `PILIN` (Sec/SPIII).
- **Node:** label `SignalPeptideType`, id `signalp:<TYPE>` (e.g. `signalp:SP`); props `name`, `signalp_id`,
  `level = 0`.
- **Edge:** `Gene_has_signal_peptide_type`, `properties: { probability: float, cleavage_site: int|null,
  cleavage_probability: float|null }` (`probability` = the winning class's `Likelihood`). The full 6-vector is
  intentionally **not** stored (winning prob only) — record this as a Step-0 decision; carry it as parallel
  fields only if a consumer needs it.
- **Merge (Step 2):** `join_to: protein_id` via the header-split transform. Fields: `signalp_type: str`,
  `signalp_probability: float`, `signalp_cleavage_site: int`, `signalp_cleavage_probability: float`; skip
  `Prediction == Other`. New transform `_tx_parse_cs_pos` (parse `CS_pos` → position + prob), registered in
  `_TRANSFORMS`, unit-tested in `tests/test_annotation_transforms.py`.
- **Adapter / post-import / kg-validity:** identical shape to §8.1 (flat, 1:1, scored). Validity:
  `probability ∈ [0,1]`; `cleavage_site` a positive int when present.
- **Release notes (Step 6):** `docs/kg-changes/signalp-extension.md`.

### 8.3 3B reference (no immediate target)

A per-gene numeric/scalar/free-text value with **no** shared cross-gene vocabulary (e.g. a GC% or a
codon-usage score) → Gene property per §4.4 (schema + `GeneNodeField` enum + `int_fields`/`float_fields`).
No current tool needs this; included so the decision tree's third branch has a worked reference.

---

## 9. Non-goals (explicit)

- Re-running Phase-1 tools, or authoring `/<tool>-run` skills (that is `add-a-tool`).
- **MCP / explorer integration.** `integrate-a-tool` stops at the live KG (nodes, edges, post-import rollups,
  kg-validity, snapshot) **+ the release-notes / what-changed doc** (`docs/kg-changes/<tool>-extension.md`).
  Making the new ontology queryable through the explorer MCP (`ontology_landscape` / `search_ontology` /
  `genes_by_ontology`) requires adding it to the MCP's 12-value `ontology` enum and its `level` handling — a
  separate change in the `multiomics_explorer` repo, tracked there. The release-notes doc is the hand-off
  artifact to that effort; the skill does not touch MCP code.
- **Normalizing a pre-`add-a-tool` runner's output** to the calls.json convention (e.g. signalp's
  `output.json`, §8.2). Flag it as a Step-1 prerequisite and choose how to consume it, but the durable fix
  belongs in that tool's `/<tool>-run` skill, not here.
- Refactoring `build_gene_annotations.py` into a fully config-driven merge.
- Attach targets other than Gene (Protein / Organism / per-region / per-genome-scalar).
- A machine-verified spot-check harness (the existing tools verify spot checks manually).
