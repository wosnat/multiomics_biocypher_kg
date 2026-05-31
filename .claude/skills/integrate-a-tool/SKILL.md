---
name: integrate-a-tool
description: Phase-2 partner to add-a-tool. Use when a per-strain tool already has committed `<strain>.<tool>.calls.json` artifacts and you need to wire them into the knowledge graph — triggers "integrate signalp into the KG", "wire psortb into the graph", "Phase 2 for <tool>", "turn <tool> calls.json into nodes/properties", "add the X tool's predictions to the knowledge graph". Routes the tool through one front door (a new source in the gene-annotation merge → a field in `gene_annotations_merged.json`), then forks to either an ONTOLOGY track (new node type + `Gene_has_<x>` edges, optionally carrying a per-call score) or a GENE-PROPERTY track (new Gene property), adds post-import rollups + indexes, validates against a live Docker rebuild (`/omics-edge-snapshot` + `pytest -m kg` + snapshot), and writes a release-notes / what-changed doc. Gene-only attach target. Stops at the live KG + release-notes doc — MCP/explorer integration is a separate, out-of-scope follow-up.
argument-hint: <tool-name> [--spec <path>]
user-invocable: true
allowed-tools: Read, Edit, Write, Grep, Glob, Skill, mcp__multiomics-kg__run_cypher, Bash(uv *), Bash(python *), Bash(pytest *), Bash(docker *), Bash(git *), Bash(mkdir *), Bash(jq *), Bash(diff *), Bash(bash *)
---

# Integrate a Tool

The **Phase-2 partner to [`add-a-tool`](../add-a-tool/SKILL.md)**. `add-a-tool` produces
committed, inspectable `<strain>.<tool>.calls.json` artifacts and explicitly defers all KG
coupling. This skill is that deferred half: it walks a tool from `calls.json` to a validated,
documented presence in the live Neo4j graph.

This SKILL.md is the operational front door; the deep, line-anchored detail lives in
[`references/`](references/) — read the one matching the step you're on. When in doubt, copy
[`cazy_adapter.py`](../../../multiomics_kg/adapters/cazy_adapter.py), the canonical skeleton.

| Reference | Read it at | What it carries |
|---|---|---|
| [`decision-tree.md`](references/decision-tree.md) | Step 0 | categorical-vs-property + functional-vs-structural routing; worked PSORTb/SignalP |
| [`ontology-patterns.md`](references/ontology-patterns.md) | Step 3A | the 11-ontology inventory, shared invariants, copy-target table, footguns |
| [`gene-annotation-merge-recipe.md`](references/gene-annotation-merge-recipe.md) | Step 2 | exact `build_gene_annotations.py` + config + transform edit points, line-anchored |
| [`post-import-patterns.md`](references/post-import-patterns.md) | Step 4 | rollup + index + fold-in snippets, real cazy/tcdb examples |

## When to use

- A tool already has committed `cache/data/<org>/genomes/<strain>/<tool>/<strain>.<tool>.calls.json`
  (or an equivalent Phase-1 artifact) and the user wants it in the graph.
- Triggers: "integrate <tool> into the KG", "wire <tool> into the graph", "Phase 2 for <tool>",
  "turn <tool> calls.json into nodes/properties".
- **Immediate targets:** PSORTb and SignalP (both runners exist — worked end-to-end in
  [`references/decision-tree.md`](references/decision-tree.md)).

**Do NOT use** to re-run a tool or author a `/<tool>-run` skill (that is `add-a-tool`), for
attach targets other than Gene, or to touch the explorer MCP (out of scope).

## The canonical data path

Both tracks flow through the **same front door**: the tool's calls.json becomes a new *source* in
the gene-annotation merge, producing a field in `gene_annotations_merged.json` — exactly as eggNOG
yields `cazy_ids` / `transporter_classification` that `cazy_adapter` / `tcdb_adapter` read today.

```
<strain>.<tool>.calls.json   →   new source+fields in gene_annotations_config.yaml
                                  + load_<tool>() in build_gene_annotations.py
        ▼
gene_annotations_merged.json[locus_tag]  →  field(s)
        ├── ONTOLOGY track (3A): <tool>_adapter.py → nodes + Gene_has_<x> edges (+ edge score)
        └── GENE-PROPERTY track (3B): schema_config gene property → cyanorak_ncbi_adapter reads it
```

Merging first is **mandatory** — it is what makes `contributing_sources` and the auto-generated
`DataSource` node pick the tool up for free.

> **Phase-1 artifacts may not follow the calls.json convention.** Some `/<tool>-run` skills predate
> `add-a-tool` (e.g. `signalp-run` still writes raw `signalp/output.json`). Step 1 inspects the real
> artifact and either normalizes it (a fix in the `/<tool>-run` skill) or parses the raw output in
> `load_<tool>()`. Treat this as common, not exceptional.

## Pick the track (decision tree — full version in [`references/decision-tree.md`](references/decision-tree.md))

```
Small controlled vocabulary, shared across genes, with grouping/hierarchy semantics?
  (psortb localizations, signalp types, a TC/CAZy-style tree)
   → ONTOLOGY TRACK (3A): new node type + Gene_has_<x> edges (+ <x>_is_a_<x> if levelled)

  …and the call ALSO carries a numeric score / probability / tier?
   → STILL 3A; the score rides as an EDGE PROPERTY on Gene_has_<x>
     (no existing ontology edge has properties — model the block on Changes_expression_of)

Per-gene numeric / scalar / free-text / single flag, NO shared cross-gene vocabulary?
   → GENE-PROPERTY TRACK (3B): new Gene property in schema_config
```

**Functional vs structural (drives cross-node rollups — see [`references/decision-tree.md`](references/decision-tree.md), Decision 2):** a *functional* evidence
source (KO-/domain-/role-like) folds into `Gene.annotation_types` + `informative_annotation_types`
(and maybe the `annotation_quality` bucket); a *structural* prediction (psortb localization,
signalp signal-peptide — *where* the gene is, not *what it does*) **stays out** of those.

## Templates (`assets/`)

Copy-paste starting points, grounded in `cazy_adapter` + the live schema / config / post-import:

| File | Use in | What it is |
|---|---|---|
| [`ontology_adapter_template.py`](assets/ontology_adapter_template.py) | Step 3A | `<Tool>Adapter` + `Multi<Tool>Adapter` skeleton — `level`, multi-call fan-out, edge-score payload, 1:1 variant |
| [`schema_node_edge_snippet.yaml`](assets/schema_node_edge_snippet.yaml) | Step 3 | node (+ `level`) + `gene_has_<x>` edge (+ optional `properties:` block) + hierarchy edge |
| [`gene_annotations_source_snippet.yaml`](assets/gene_annotations_source_snippet.yaml) | Step 2 | `sources:` + `fields:` + `logical_sources:` + transform stub |
| [`post_import_rollup_snippet.cypher`](assets/post_import_rollup_snippet.cypher) | Step 4 | flat + hierarchical rollups, indexes, `annotation_types` fold-in, optional `rank_by_score` |
| [`what_changed_template.md`](assets/what_changed_template.md) | Step 6 | the release-notes / what-changed doc skeleton |

## The 7-step workflow

Create a TodoWrite item per step. Each step links the matching `references/` file for full detail.

- **Step 0 — Capture intent + write spec/plan.** First **check for an existing Phase-1 design doc**
  (the `/<tool>-run` SKILL.md links it); its deferred "Phase 2" sketch may pick a *different* track than
  this skill's decision tree — follow the tree and note the supersession. Then four bullets: (1) what the
  tool predicts + value space; (2) calls.json shape + join key (WP_ vs locus_tag) + nulls/sentinels;
  (3) the categorical-vs-property verdict + functional-vs-structural verdict (routing:
  [`references/decision-tree.md`](references/decision-tree.md)); (4) target node/edge labels + id
  prefix + hierarchy levels OR target Gene property name+type. Expand into
  `docs/superpowers/specs/<date>-<tool>-kg-integration-design.md` + a checkbox plan.

- **Step 1 — Inspect the Phase-1 runner & calls.json.** Read the `/<tool>-run` SKILL.md; `jq` the
  committed calls.json to enumerate distinct values, the join key, nulls/sentinels, and per-call
  scores. (`allowed-tools` matches *bare* commands, so run **one `jq` per Bash call** — compound
  `a; b` lines and `for` loops over strains won't match the allow-list.) **Reuse existing parse utils**
  — check `multiomics_kg/utils/<tool>.py` first; the generic `multiomics_kg/utils/tool_calls_io.py` may
  not exist yet (copy `add-a-tool`'s template if needed). Keep new parsing pure (unit-testable, no filesystem).

- **Step 2 — Wire into the gene-annotation merge** (both tracks). Full line-anchored edit points in
  [`references/gene-annotation-merge-recipe.md`](references/gene-annotation-merge-recipe.md). In short:
  - `config/gene_annotations_config.yaml`: add a `sources:` + `fields:` + `logical_sources:` entry
    (only `logical_sources` is read by code; `type`/`path_pattern`/`join_to` are documentary).
  - `build_gene_annotations.py` (the hand-wired part): add `load_<tool>()`, then thread a new source
    dict through `_get_raw` → **all six** `_resolve_*` → `build_wide`/`build_merged` → the
    `process_strain` row-level `.get(<join_key>)` join, and add a `_compute_contributing_sources` check.
    Also patch `extract_first_match_in_sources` in `annotation_helpers.py` *iff* a field uses
    `extract_first_match`.
  - ⚠️ **Two test-breakage landmines** (recipe carries the fix): (1) make the new source arg
    *optional-defaulted* (`tool: dict | None = None`) — `test_build_gene_annotations.py` calls
    `build_*`/`_resolve_union` positionally ~120×, so a required positional detonates the suite;
    (2) update the exact DataSource node-count in **both** `test_data_source_adapter.py` (unit) and
    `tests/kg_validity/test_data_source.py` (live, Step 5) **and** add the tool to
    `data_source_adapter.py`'s `_name_for`/`_description_for` (else the DataSource node is named
    "Psortb" with an empty description).
  - If a raw value needs reshaping, add a `_tx_<name>` to `annotation_transforms.py` (register in
    `_TRANSFORMS`, unit-test it).
  - **Verify:** `bash scripts/prepare_data.sh --steps 2 --strains MED4 --force`, then `jq` the field
    landed + `contributing_sources` lists the tool; and `pytest tests/test_build_gene_annotations.py
    tests/test_data_source_adapter.py -q` is green.

- **Step 3 — Fork to the matching track.** (3A invariants + copy targets:
  [`references/ontology-patterns.md`](references/ontology-patterns.md).)
  - **3A Ontology:** vocab/hierarchy helper in `multiomics_kg/utils/<tool>.py` → `schema_config.yaml`
    node (mandatory `level: int`, `0`=root; `level_kind` optional/flat) + `gene_has_<x>` edge (add a
    `properties:` block **iff a score rides**) + `<x>_is_a_<x>` if levelled → `<tool>_adapter.py`
    (`<Tool>Adapter` per-strain + `Multi<Tool>Adapter` owns nodes), copied from `cazy_adapter`,
    `_clean_str` on strings → wire `Multi<Tool>Adapter` into `create_knowledge_graph.py` →
    `tests/test_<tool>.py`.
  - **3B Gene-property:** add the property under `gene:` in `schema_config.yaml` **AND** add the
    matching `GeneNodeField` enum member in `cyanorak_ncbi_adapter.py` (+ the `int_fields`/`float_fields`
    set for numerics — the read is **not** automatic over arbitrary merged fields); spot-check it on a
    built Gene node; unit-test the parse.

- **Step 4 — Post-import rollups + indexes** (edit **both** `post-import.sh` and `post-import.cypher`;
  keep the Cypher byte-identical). Snippets:
  [`references/post-import-patterns.md`](references/post-import-patterns.md). On the new node:
  `gene_count`/`organism_count` (+ `member_count` if
  hierarchical) via `CALL { … } IN TRANSACTIONS` (the cazy `*0..` subtree walk for hierarchical, no
  `*0..` for flat); scalar indexes on `level` (+ `level_kind` if hierarchical) + the id; full-text on
  `name`. **Cross-node rollups are conditional:** *functional* → fold the edge into Gene
  `annotation_types` + `informative_annotation_types` (+ `annotation_quality` bucket decision, which
  follows the bucket-maintenance procedure and shifts `genes_by_function` `min_quality`); *structural*
  → stay out, optionally add a `<tool>_count` Gene signal. `OrganismTaxon`/`Publication` need nothing.
  If a score rides, optionally compute a per-gene `rank_by_<score>` (mirrors expression
  `rank_by_effect`). Capture `scripts/post-import-validate.sh > baseline.txt` before rebuild; after,
  confirm **additions-only**.

- **Step 5 — Live KG validation gate** (definition of "done"): `/omics-edge-snapshot` (before) →
  `docker compose up -d --build` (background, ~1h) → post-import → `/omics-edge-snapshot` (after:
  expression edges unchanged, new edges appear) → add kg-validity assertions (node count > 0, edge
  presence, rollup sanity, no-orphan, numeric edge-property range checks à la `test_expression.py`) →
  `pytest -m kg` → regenerate `snapshot_data.json`.

- **Step 6 — Release-notes / what-changed doc + CLAUDE.md.** Write
  `docs/kg-changes/<tool>-extension.md` (new node/edge types + counts, property-additions table,
  indexes, vocabulary, conventions, validation results, out-of-scope) → update CLAUDE.md "Key graph
  facts" + "Actual Neo4j labels" + link the doc.

- **Step 7 — Register / hand-off.** Confirm `add-a-tool`'s "Phase 2 hand-off" redirects here; confirm
  `add-a-strain` step-10 ("loop back through prepare_data") covers the new merged field so future
  strains pick the tool up automatically. Record the MCP-registration follow-up (out of scope) in the
  release-notes doc.

## Worked examples

PSORTb (`SubcellularLocalization`) and SignalP (`SignalPeptideType`) are fully worked in
[`references/decision-tree.md`](references/decision-tree.md) — both are **3A flat ontology + edge
score, 1:1, skip the no-signal sentinel** (`Unknown`/`Other`). SignalP carries a
calls.json-normalization prerequisite. Copy those when implementing either.

## Reference checklist (final pass before merging)

- [ ] Step 0 spec under `docs/superpowers/specs/<date>-<tool>-kg-integration-design.md` + checkbox plan
- [ ] calls.json inspected; distinct values + join key + score-presence enumerated; parse utils reused
- [ ] Categorical-vs-property AND functional-vs-structural verdicts recorded; merged-field shape chosen
- [ ] `gene_annotations_config.yaml`: `sources:` + `fields:` + `logical_sources:` added
- [ ] New `_tx_<name>` added + registered in `_TRANSFORMS` + unit-tested (if reshaping was needed)
- [ ] `build_gene_annotations.py`: `load_<tool>()` + `_get_raw` + all six `_resolve_*` + `build_wide`/`build_merged` + the `process_strain` row-level `.get(<join_key>)` join + `contributing_sources` (+ `extract_first_match_in_sources` iff used)
- [ ] New source arg is **optional-defaulted** on `build_wide`/`build_merged`/`_resolve_union` (tests call them positionally); `pytest tests/test_build_gene_annotations.py` green
- [ ] DataSource node-count updated in **both** `test_data_source_adapter.py` (unit) **and** `tests/kg_validity/test_data_source.py` (live) **and** tool added to `data_source_adapter.py` `_name_for`/`_description_for`
- [ ] MED4 rebuild verifies the field in `gene_annotations_merged.json` + `contributing_sources` + DataSource node (with a real name, not `"Psortb"`)
- [ ] **3A:** schema node carries `level` (0=root); `level_kind` set when hierarchical; edge has a `properties:` block iff a score rides; node id is colon only for a **registered** bioregistry prefix, else underscore (`psortb_OuterMembrane`) — worked examples + test assertions match the real form
- [ ] **3A:** `<tool>_adapter.py` (per-strain + Multi*) emits multi-call fan-out (one edge per call, 1:1 for psortb/signalp), `_clean_str`, wired into `create_knowledge_graph.py`
- [ ] **3B:** Gene property declared in schema **and** added to `GeneNodeField` enum (+ `int_fields`/`float_fields` for numerics); confirmed on a built node
- [ ] `tests/test_<tool>.py` (pure parsing/vocab) passes
- [ ] post-import.sh + post-import.cypher edited identically (`diff` clean); `post-import-validate.sh` baseline/after shows additions-only
- [ ] Scalar indexes (`level`, `level_kind` if hierarchical, id) + full-text (`name`) for the new node
- [ ] Cross-node rollups: **functional** → Gene `annotation_types`/`informative_annotation_types` folded in (+ `annotation_quality` bucket decision); **structural** → kept out, optional `<tool>_count` only. `OrganismTaxon`/`Publication`: nothing
- [ ] `/omics-edge-snapshot` before/after: expression edges unchanged, new edges appear
- [ ] kg-validity assertions added; `pytest -m kg` passes; `snapshot_data.json` regenerated
- [ ] `docs/kg-changes/<tool>-extension.md` written; CLAUDE.md "Key graph facts" + "Actual Neo4j labels" updated + linked
- [ ] `add-a-tool` Phase-2 hand-off redirects here; `add-a-strain` step-10 covers the new field; MCP follow-up recorded as out-of-scope
- [ ] `pytest -m "not slow and not kg"` passes
