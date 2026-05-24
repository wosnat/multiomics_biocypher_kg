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
rebuild, and documented for downstream MCP/explorer consumers.

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
  "psortb_localizations": [
    {"call": "OuterMembrane", "score": 9.5},
    {"call": "Extracellular",  "score": 7.2}
  ]
  ```
  How the score travels through the merge is a per-tool Step-0 decision: a **structured `passthrough`** of
  the calls.json object (preferred) vs. **parallel scalar fields**. The existing `union`/`single`/
  `passthrough_list` rules are string/list-oriented; carrying a list-of-dicts uses `passthrough` of an
  already-structured value.

### 4.2 The ontology node (3A)

Every new ontology node **must** carry, adapter-emitted:

| Property | Type | Requirement |
|---|---|---|
| `name` | str | human-readable display; `_clean_str`-sanitized |
| `<tool>_id` (e.g. `cazy_id`) | str | the bare identifier / CURIE local part |
| `level` | int | **0 = broadest/root**; increases with specificity. Mandatory — the MCP `ontology_landscape` / `search_ontology` tools rely on it. |
| `level_kind` | str | e.g. `psortb_class`; mandatory companion to `level` |
| `gene_count` | int | post-import |
| `organism_count` | int | post-import |
| `member_count` | int | post-import, only if hierarchical |

`level` / `level_kind` are set by the adapter from the hierarchy depth (see `_classify` in
`cazy_adapter`), **not** by post-import. Flat (single-level) ontologies set `level = 0` for every node.

### 4.3 The `Gene_has_<X>` edge (3A)

- **Bare categorical:** property-less (`{}`), like `Gene_has_cazy_family`.
- **With score:** a real `properties:` block in `schema_config`, e.g. `score: float`, `probability: float`,
  `confidence_tier: str`. The adapter yields the per-call value in the edge-property dict. Floats stay
  numeric (no `_clean_str`); only string sub-fields (a tier label) are sanitized.

### 4.4 The Gene property (3B)

Declared under `gene:` `properties:` in `schema_config`; `cyanorak_ncbi_adapter` auto-reads any field
present in `gene_annotations_merged.json`. Numeric → `float`/`int`; flag → `str` ("true"/"false") or
`bool` per existing convention; free-text → `str`.

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
create_knowledge_graph.py                                       # wire Multi<Tool>Adapter                   (3A only)
scripts/post-import.sh + scripts/post-import.cypher             # rollups + indexes, kept byte-identical     (Step 4)
tests/test_<tool>.py , tests/test_annotation_transforms.py      # unit tests                                 (Steps 2,3)
tests/kg_validity/test_<tool>.py (or extend existing)           # live-graph assertions                      (Step 5)
docs/kg-changes/<tool>-extension.md                             # what-changed doc                           (Step 6)
CLAUDE.md                                                       # Key graph facts + Actual Neo4j labels      (Step 6)
```

```
.claude/skills/integrate-a-tool/
├── SKILL.md
├── references/
│   ├── decision-tree.md              # expanded tree + worked PSORTb/SignalP examples
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
associated score. **Reuse existing parse utils** — check `multiomics_kg/utils/<tool>.py` and
`tool_calls_io.py` first; write new parsing only if absent, and keep it pure (unit-testable, no filesystem).

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
- **3B Gene-property track:** add the property under `gene:` in `schema_config.yaml` (§4.4); confirm
  `cyanorak_ncbi_adapter` auto-reads it by spot-checking a built Gene node; unit-test the parse/transform.

**Step 4 — Post-import rollups + indexes** (edit **both** `post-import.sh` and `post-import.cypher`; keep
byte-identical):
- 3A: `gene_count` / `organism_count` (+ `member_count` if hierarchical) on the new node via the
  `MATCH … CALL { … } IN TRANSACTIONS` pattern (the cazy `*0..` subtree walk for hierarchical, no `*0..`
  for flat); extend Gene `annotation_types` + `informative_annotation_types` to include the new edge, and
  add a `<tool>_count` Gene routing signal; scalar indexes on `level` + `level_kind` + the id, and a
  full-text index on `name`. If a score rides on the edge, optionally compute a per-gene `rank_by_<score>`
  (mirrors expression `rank_by_effect`).
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

**Step 6 — What-changed doc + CLAUDE.md.** Write `docs/kg-changes/<tool>-extension.md` from the template
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
- [ ] **3A:** schema node carries `level` (0=root) + `level_kind`; edge has a `properties:` block iff a score rides
- [ ] **3A:** `<tool>_adapter.py` (per-strain + Multi*) emits multi-call fan-out (one edge per call), `_clean_str` on strings, wired into `create_knowledge_graph.py`
- [ ] **3B:** Gene property declared in schema; `cyanorak_ncbi_adapter` auto-read confirmed on a built node
- [ ] `tests/test_<tool>.py` (pure parsing/vocab) passes
- [ ] post-import.sh + post-import.cypher edited identically (`diff` clean); `post-import-validate.sh` baseline/after shows additions-only
- [ ] Scalar indexes (`level`, `level_kind`, id) + full-text (`name`) for the new node; Gene `annotation_types` + `<tool>_count` extended
- [ ] `/omics-edge-snapshot` before/after: expression edges unchanged, new edges appear
- [ ] kg-validity assertions added; `pytest -m kg` passes; `snapshot_data.json` regenerated
- [ ] `docs/kg-changes/<tool>-extension.md` written; CLAUDE.md "Key graph facts" + "Actual Neo4j labels" updated + linked
- [ ] `add-a-tool` Phase-2 hand-off redirects here; `add-a-strain` step-10 covers the new field
- [ ] `pytest -m "not slow and not kg"` passes

---

## 8. Worked examples the references cite

| Tool | Track | Why | Edge payload |
|---|---|---|---|
| `cazy` | 3A hierarchical ontology | observed-only class/family/subfamily tree | none (`{}`) — the canonical skeleton to copy |
| `tcdb` | 3A hierarchical ontology | 5-level TC hierarchy + substrate bridge | none on `Gene_has_tcdb_family` |
| `psortb` | 3A flat ontology + edge score | localization vocab + per-call score | `score: float` on the edge — the canonical edge-property example |
| `signalp` | 3A flat ontology + edge score | type {SP,LIPO,TAT,PILIN} + probability | `probability: float` |
| (hypothetical numeric) | 3B Gene property | e.g. GC%, codon-usage score | n/a |

---

## 9. Non-goals (explicit)

- Re-running Phase-1 tools, or authoring `/<tool>-run` skills (that is `add-a-tool`).
- Refactoring `build_gene_annotations.py` into a fully config-driven merge.
- Attach targets other than Gene (Protein / Organism / per-region / per-genome-scalar).
- A machine-verified spot-check harness (the existing tools verify spot checks manually).
