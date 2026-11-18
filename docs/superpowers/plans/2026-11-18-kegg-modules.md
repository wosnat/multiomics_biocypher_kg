# KEGG Modules Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the `KeggModule` ontology (modules + `ko00002` class nodes) and the per-organism `Organism_has_kegg_module` completeness edge to the KG.

**Architecture:** Step 6 downloads KEGG's module data and writes `modules` (class path inline) / `module_ec_to_kos` / `kegg_release` blocks into `kegg_data.json`. A pure-Python evaluator (`utils/kegg_modules.py`) parses module definitions and scores one organism's KOs (plus EC-only genes) against them. A new `kegg_module_adapter.py` emits nodes and the four edge types at build time; post-import adds rollups, indexes and `Schema_info.kegg_release`.

**Tech Stack:** Python 3, BioCypher, Neo4j / Cypher, pytest, KEGG REST.

**Spec:** `docs/superpowers/specs/2026-11-18-kegg-modules-design.md` (approved; §4–§5 and §11–§13 are binding).

**Deliberate deviations from spec §6 (implementation detail, no graph change):** no separate `module_classes` block (each module carries its `class_path`); no `link/ko/module` download (KO membership and step numbers come from parsing the definition, which is the authoritative source, and `link/ko/module` cannot give steps or optional subunits).

## Global Constraints

- Node ids: modules `kegg.module:M00531`; class nodes `kegg.module:<slug>[.<slug>[.<slug>]]`, slug = lowercase, runs of non-alphanumerics → `_`, stripped (e.g. `kegg.module:pathway_modules.energy_metabolism.nitrogen_metabolism`). Pass through `normalize_curie(raw) or raw` (verified: keeps both forms unchanged).
- `level` / `level_kind`: 0 `kegg_module_type`, 1 `kegg_module_class`, 2 `kegg_module_subclass`, 3 `kegg_module`.
- `module_type` closed: `pathway` | `signature` (spec §4: KEGG's 2 `complex` modules have no KO in the KG); dense on every `KeggModule` node (class nodes inherit their level-0 ancestor's). `kegg_data.json` keeps all three KEGG types; if the Task 3 dry run ever keeps a `complex` module, add the value to the vocabulary then.
- `status` closed: `complete` (0 missing) | `one_step_missing` (1) | `two_steps_missing` (2) | `partial` (≥ 3). `value_descriptions.complete` = "every step has an annotated gene; not shown to run".
- `component` closed: `required` | `optional`.
- Step numbering: 1-based over the steps that remain after joining multi-line definitions and **dropping whole-step `--` placeholders**; `steps` on `Kegg_term_in_kegg_module`, `missing_step_indices`, `ec_filled_step_indices` all use it. `definition` stays verbatim.
- EC fallback: only genes with **no live KO**; only **full** ECs (no `-`); never for a KO inside a `+` complex; at most one step per gene per module (maximum bipartite matching); only over steps KOs leave uncovered.
- `ec_genes` is parallel to `ec_filled_step_indices`; a module-reference step's slot holds the referenced module id (`M0xxxx`).
- `missing_kos`: KOs written directly in the missing steps, sorted, deduplicated; a missing `M0xxxx` step contributes none; **optional subunits (`-K…`) are left out** (they never complete a step).
- Module-reference steps: in the KO-only pass an `M0xxxx` term is true iff the referenced module is `complete` with `steps_present_via_ec == 0`; in the EC pass iff it is `complete` at all. A step made true by a reference is never sent to the gene matching (it consumes no gene).
- KEGG release skew: the six module raws must come from the same KEGG release as the existing `cache/data/kegg/raw/` files (all fetched 2026-11-18). If they are fetched on a later day, run step 6 with `--refetch-raw` instead, so `kegg_release` describes every KEGG input.
- Organism edge only when `steps_present >= 1`. No `completeness` property. No gene → module edge. No change to `annotation_types` / `annotation_quality` / BRITE.
- Every string property through the adapter's local `_clean_str` (`'` → `^`, `|` removed); never build values with `|`.
- `KeggModule.gene_count` / `organism_count` / `direct_gene_count` count KO-bridge genes only (`Gene_has_kegg_ko` → `Kegg_term_in_kegg_module`); `direct_gene_count` only on `level = 3`.
- Both `scripts/post-import.sh` and `scripts/post-import.cypher` get identical Cypher; every post-import property declared in `config/schema_config.yaml` with a `# post-import:` comment.

## Review Focus

1. **A module whose referenced module was itself completed through EC** — the referencing step must count in `steps_present_via_ec`, with `M0xxxx` in its `ec_genes` slot. Test in Task 1.
2. **A strain with no `gene_annotations_merged.json`, or with no KOs/ECs at all** — no edges, no crash, and registry comment lines (`#…`, e.g. the commented MIT1002 draft) never read. Test in Task 3.
3. **Definition shapes that break naive parsers**: leading optional subunit (`-K00003+K00001`), trailing optional (`K1+K2-K3`), nested complex `K1+(K2,K3)`, an optional subunit before a group `(K1,K2)-K3`. Tests in Task 1.
4. **A deleted KO named in a module definition** — never counted present, and no `Kegg_term_in_kegg_module` edge emitted for it (no `KeggTerm` node to land on). Test in Task 3.
5. **Step 6 run with `--force` on a checkout whose committed raw manifest predates the module raws** — the new raw files are fetched once and recorded, not reported as "missing raw input". Test in Task 2.

---

### Task 1: Module definition evaluator

**Files:**
- Create: `multiomics_kg/utils/kegg_modules.py`
- Test: `tests/test_kegg_modules.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) class Step: index: int; definition: str` (definition = the step's formula text, `--` steps excluded, index 1-based)
  - `parse_steps(definition: str) -> list[Step]` — multi-line (`\n`) definitions are joined as further steps.
  - `ko_roles(steps: list[Step]) -> dict[str, tuple[list[int], str]]` — KO → (sorted step indices, `required`|`optional`); raises `ValueError` if a KO is required in one step and optional in another.
  - `@dataclass(frozen=True) class ModuleResult` with fields `steps_total: int, steps_present: int, steps_present_via_ec: int, missing_step_indices: list[int], missing_step_definitions: list[str], missing_kos: list[str], ec_filled_step_indices: list[int], ec_genes: list[str]` and property `status -> str`.
  - `evaluate_organism(modules: dict[str, list[Step]], kos: set[str], ec_candidates: dict[str, set[str]]) -> dict[str, ModuleResult]` — `ec_candidates` maps KO → locus tags of EC-only genes standing in for it. Returns results for every module (callers drop `steps_present == 0`). Evaluates referenced modules first; raises `ValueError` on a reference cycle or a reference to a module not in `modules`.

- [ ] **Step 1: Write the failing tests**

```python
from multiomics_kg.utils.kegg_modules import parse_steps, ko_roles, evaluate_organism

NITRATE = "(K00367,K10534,K00372+K00360) (K00366,K17877)"           # M00531
GLYC = "K01803 (((K00134,K00150) K00927),K11389) (K01834,K15633,K15634,K15635) (K01689,K27394) (K00873,K12406)"

def ev(defn, kos, ec=None, extra=None):
    mods = {"M1": parse_steps(defn), **(extra or {})}
    return evaluate_organism(mods, set(kos), ec or {})["M1"]

def test_steps_and_dash_placeholder():
    assert [s.definition for s in parse_steps("K1 -- K2")] == ["K1", "K2"]   # '--' step dropped
    assert [s.index for s in parse_steps("K00001 -- K00002")] == [1, 2]
    assert len(parse_steps("K00001 K00002\nK00003")) == 3                     # multi-line joined

def test_or_and_complex_optional():
    assert ev(NITRATE, {"K00366", "K00372"}).steps_present == 1               # complex needs K00360 too
    assert ev(NITRATE, {"K00366", "K00372", "K00360"}).status == "complete"
    assert ev("K00001+K00002-K00003", {"K00001", "K00002"}).status == "complete"   # -K optional
    assert ev("K00001+(K00002,K00003)", {"K00001", "K00003"}).status == "complete"
    assert ev("(K00001,K00002)-K00003", {"K00002"}).status == "complete"
    assert ev("-K00003+K00001", {"K00001"}).status == "complete"               # leading optional

def test_missing_fields_and_status():
    r = ev(NITRATE, {"K00366"})
    assert (r.status, r.missing_step_indices) == ("one_step_missing", [1])
    assert r.missing_step_definitions == ["(K00367,K10534,K00372+K00360)"]
    assert r.missing_kos == ["K00360", "K00367", "K00372", "K10534"]

def test_ec_fallback_fills_plain_term():
    r = ev(GLYC, {"K01803", "K00134", "K00927", "K01834", "K01689"}, {"K00873": {"PMN2A_0295"}})
    assert (r.status, r.steps_present_via_ec, r.ec_filled_step_indices, r.ec_genes) == ("complete", 1, [5], ["PMN2A_0295"])

def test_ec_never_fills_complex_subunit():
    r = ev(NITRATE, {"K00366"}, {"K00372": {"g1"}, "K00360": {"g1"}})
    assert r.steps_present == 1 and r.ec_genes == []

def test_ec_one_step_per_gene_and_ko_wins():
    r = ev("K00001 K00002", set(), {"K00001": {"g1"}, "K00002": {"g1"}})
    assert r.steps_present == 1
    r = ev("K00001 K00002", {"K00001"}, {"K00001": {"g1"}, "K00002": {"g1"}})
    assert (r.steps_present, r.steps_present_via_ec, r.ec_filled_step_indices) == (2, 1, [2])

def test_module_reference_via_ec_propagates():
    sub = {"M2": parse_steps("K00001")}
    r = ev("M2 K00002", {"K00002"}, {"K00001": {"g1"}}, extra=sub)
    assert (r.status, r.steps_present_via_ec, r.ec_genes) == ("complete", 1, ["M2"])
    assert ev("M2 K00002", {"K00002"}, {}, extra=sub).missing_kos == []

def test_missing_kos_skip_optional_and_ko_pass_ignores_ec_completed_reference():
    assert ev("K00001+K00002-K00003", set()).missing_kos == ["K00001", "K00002"]
    sub = {"M2": parse_steps("K00001")}
    r = ev("(M2,K00009) K00002", {"K00002"}, {"K00001": {"g1"}}, extra=sub)
    assert (r.ec_filled_step_indices, r.ec_genes) == ([1], ["M2"])   # KO-only pass sees step 1 missing

def test_ko_roles_and_conflict():
    assert ko_roles(parse_steps(NITRATE))["K00360"] == ([1], "required")
    assert ko_roles(parse_steps("K00001+K00002-K00003"))["K00003"] == ([1], "optional")
    import pytest
    with pytest.raises(ValueError):
        ko_roles(parse_steps("K00001-K00003 K00003"))
```

- [ ] **Step 2: Run** `uv run pytest tests/test_kegg_modules.py -v` — expected: FAIL (module not found).

- [ ] **Step 3: Implement `multiomics_kg/utils/kegg_modules.py`.** Start from the dry-run evaluator (`tokenize` → `K\d{5}|M\d{5}|--|[(),+\- ]`; split top-level on space for steps, on `,` for OR, on space for AND inside groups, on `+`/`-` for complexes with `-` marking an optional subunit). Parse each step into a tree once; evaluate it twice: KO-only, then with `ec_candidates` admitted on non-complex KO terms. EC assignment: for the steps true with EC but false on KOs, candidate genes = union of `ec_candidates[k]` over non-complex KOs in the step; run augmenting-path bipartite matching (step ↔ gene); a step with no matched gene stays missing. A step that is only true through a referenced module whose own result has `steps_present_via_ec > 0` counts as EC-filled with the module id as its gene slot (it does not enter the matching).

- [ ] **Step 4: Run** `uv run pytest tests/test_kegg_modules.py -v` — expected: PASS.

- [ ] **Step 5: Commit** — `git add multiomics_kg/utils/kegg_modules.py tests/test_kegg_modules.py && git commit -m "feat(kegg): KEGG module definition evaluator with EC fallback"`

---

### Task 2: Step 6 — module raws and `kegg_data.json` blocks

**Files:**
- Modify: `multiomics_kg/utils/kegg_utils.py` (URLs, `_RAW_FILES`, new download + parsers)
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`_raw_inputs`, `_check_raw_inputs`, `main` → `kegg_rest`)
- Modify: `multiomics_kg/utils/metabolite_identity/pipeline.py` (pass the new keys through unchanged, the `deleted_kos` pattern)
- Test: `tests/test_kegg_utils_modules.py`, `tests/test_build_kegg_metabolism_xrefs.py`

**Interfaces:**
- Consumes: `parse_steps` (Task 1) — only to validate that every definition parses.
- Produces, in `kegg_data.json`:
  - `"modules": {"M00531": {"name": str, "module_type": "pathway"|"signature"|"complex", "definition": str, "class_path": [type, class, subclass], "pathways": ["ko00910", ...]}}` — **all** KEGG modules (583; the adapter prunes). `pathways` converted `map…` → `ko…`. Modules absent from `br:ko00002` get `class_path: []` and a warning.
  - `"module_ec_to_kos": {"2.7.1.40": ["K00873", ...]}` — `/link/ko/enzyme`, restricted to KOs named in some module definition.
  - `"kegg_release": "116.0+/11-18, Nov 25"` style string — the release field of `/info/kegg`'s `kegg` line.
- New raw files under `cache/data/kegg/raw/`: `list_module.txt`, `link_pathway_module.txt`, `link_ko_enzyme.txt`, `br_ko00002.json`, `get_modules.txt` (all `get/M…` responses concatenated, batches of 10), `info_kegg.txt`. All listed in `_raw_inputs`.
- `kegg_utils.download_module_raw(raw_dir: Path, force: bool) -> None`; `kegg_utils.parse_module_entries(text: str) -> dict[str, dict]` (NAME, DEFINITION multi-line kept with `\n`, type from the `ENTRY` line); `kegg_utils.parse_module_classes(brite_json: dict) -> dict[str, list[str]]`; `kegg_utils.parse_kegg_release(text: str) -> str`.

- [ ] **Step 1: Write failing tests** in `tests/test_kegg_utils_modules.py` on small fixture strings: `parse_module_entries` returns `{"M00031": {..., "definition": "K05826 ...\nK05827"}}` for a two-line DEFINITION and `module_type == "pathway"` from `ENTRY M00031 Pathway Module`; `parse_module_classes` maps a module under `Pathway modules / Amino acid metabolism / Lysine metabolism` to that 3-list; `parse_kegg_release` extracts the release from a fixture `info/kegg` text; `_require_parsed` raises on non-empty text that yields no entries (fail-loud parsers). In `tests/test_build_kegg_metabolism_xrefs.py`: a `--force` run whose committed manifest lacks the six module raw files downloads them (stubbed `download_module_raw`) instead of raising; the written `kegg_rest` carries `modules`, `module_ec_to_kos`, `kegg_release`; `module_ec_to_kos` drops KOs outside every definition.
- [ ] **Step 2: Run** `uv run pytest tests/test_kegg_utils_modules.py tests/test_build_kegg_metabolism_xrefs.py -v` — expected: new tests FAIL.
- [ ] **Step 3: Implement.** `download_module_raw` reads `list_module.txt` to get ids, then fetches `get/M…+M…` in batches of 10 (KEGG's `get` limit) only when `get_modules.txt` is missing or `force`. In `_check_raw_inputs`, only the six module raw files get the new-input treatment: when one is missing **and** the committed manifest has no entry for it, call `download_module_raw` before the checks, and leave such unpinned paths out of the `changed` comparison (today `pinned.get(k) != v` would flag them, since `None != sha`). Every other raw keeps today's fail-loud check, and a module raw the manifest *does* pin is checked like any other. The run's manifest write then pins them. Build the three blocks after `download_kegg_raw` and add them to the `kegg_rest` dict at `build_kegg_metabolism_xrefs.py:2115`; make `pipeline.py:512` keep them (`kegg_data = dict(kegg_rest, …)` already copies unknown keys — assert it in a test rather than change code if so). Every definition must `parse_steps` without error, else fail the step.
- [ ] **Step 4: Run** the two test files — expected: PASS. Then `bash scripts/prepare_data.sh --steps 6 --force` and check: `python -c "import json;d=json.load(open('cache/data/kegg/kegg_data.json'));print(len(d['modules']),len(d['module_ec_to_kos']),d['kegg_release'])"` → `583 <n> <release>`, and `git diff --stat cache/data/kegg/kegg_data.json` shows only the new blocks (+ the manifest file).
- [ ] **Step 5: Commit** — code, tests, `cache/data/kegg/kegg_data.json`, `cache/data/step6_raw_inputs.json`: `feat(step6): KEGG modules, module classes, EC→KO links and KEGG release`.

---

### Task 3: `kegg_module_adapter.py`

**Files:**
- Create: `multiomics_kg/adapters/kegg_module_adapter.py`
- Test: `tests/test_kegg_module_adapter.py`

**Interfaces:**
- Consumes: `parse_steps`, `ko_roles`, `evaluate_organism`, `ModuleResult` (Task 1); the Task 2 blocks via `kegg_utils.load_kegg_data(cache_root)`.
- Produces: `class MultiKeggModuleAdapter(genome_config_file: str, cache_root: Path, ko_node_ids: set[str], pathway_node_ids: set[str], test_mode: bool = False)` with `get_nodes()` / `get_edges()` (BioCypher tuples) and `module_node_ids() -> set[str]`.
  - `ko_node_ids` = bare KO ids the KEGG adapter emits (`MultiKeggAnnotationAdapter.all_ko_ids()`); `pathway_node_ids` = bare `ko…` pathway ids it emits (add a `MultiKeggAnnotationAdapter.all_pathway_ids() -> set[str]` if none exists, same style as `all_ko_ids`).
  - Node label input `kegg module`; edge labels `kegg_term_in_kegg_module`, `kegg_module_is_a_kegg_module`, `kegg_module_in_kegg_pathway`, `organism_has_kegg_module`.
  - Edge ids: `{ko}-in_module-{M}`, `{child_id}-is_a-{parent_id}`, `{M}-in_pathway-{ko_pathway}`, `{organism_id}-has_module-{M}`.
  - Organism id: `insdc.gcf:<ncbi_accession>` from the registry row (comment lines skipped, one per `data_dir`).

**Behaviour:**
- Per strain: `kos` = the strain's genes' `kegg_ko` (bare `K…` ids in the merged JSON) ∩ `ko_node_ids`; `ec_candidates[k]` = locus tags of genes whose live-KO list is empty and that carry a full EC with `k ∈ module_ec_to_kos[ec]`.
- Module node set = modules with ≥ 1 KO in `ko_node_ids` ∪ modules with ≥ 1 organism edge (dry run: 384 + 30 = 414); class nodes = ancestors of kept modules; modules with `class_path == []` get no `is_a` edge.
- `Kegg_term_in_kegg_module` only for KOs in `ko_node_ids` (dangling-proof); props `steps`, `component` from `ko_roles`.
- `Kegg_module_in_kegg_pathway` only for `pathway_node_ids`.
- Organism edge props exactly the spec §4 table (minus `completeness`): `status, steps_total, steps_present, steps_present_via_ec, missing_step_indices, missing_step_definitions, missing_kos`, plus `ec_filled_step_indices` / `ec_genes` only when non-empty.
- Node props per spec §4 (`name, level, level_kind, module_type`, + `module_id, definition, step_count` on modules). Logs: node counts by level, edge counts by type, EC-filled step total.

- [ ] **Step 1: Write failing tests** with a temp `cache_root` holding a 3-module `kegg_data.json` (M1 references M2; one complex step; one module under a deleted KO only) and a two-strain registry CSV (one row commented out, one strain dir without annotations): assert node ids/levels/`module_type` on class nodes; `steps`/`component` on KO edges; no KO edge for a KO outside `ko_node_ids`; the organism edge for an EC-only gene carries `ec_genes == ["LT_1"]`; the strain without annotations yields no edges; the commented row is not read; no edge with `steps_present == 0`; every property string clean of `'` and `|`.
- [ ] **Step 2: Run** `uv run pytest tests/test_kegg_module_adapter.py -v` — FAIL.
- [ ] **Step 3: Implement** per the Interfaces/Behaviour above.
- [ ] **Step 4: Run** — PASS. Dry-run against real data: `uv run python -c` instantiating the adapter with the real KEGG adapter's id sets; expected ≈ 414 module nodes, 57 class nodes, 1,723 KO edges, ≈ 26,005 organism edges, ≈ 6,987 `complete` (spec §2; small drift from a KEGG refresh is fine, record the actual numbers for Task 5). In the same dry run, list M00161 / M00163 status for **every** cyanobacterial strain (the spec measured only six): Task 5's "complete in every cyanobacterium" assertion is written from this list, and any strain that is not complete is reported to the user with its missing step before the test is written.
- [ ] **Step 5: Commit** — `feat(kegg): KeggModule adapter with per-organism module completeness`.

---

### Task 4: Schema, vocabularies, build wiring, post-import

**Files:**
- Modify: `config/schema_config.yaml` — node `kegg module` (`preferred_id: kegg.module`, properties from spec §4 incl. post-import `gene_count`, `direct_gene_count`, `organism_count`, `complete_organism_count` with `# post-import:` comments); the four edges (model `organism has metabolite` for the organism source label, `kegg term in brite category` for the KO edge).
- Modify: `config/controlled_vocabularies.yaml` — `KeggModule.module_type`, `KeggModule.level_kind`, `Kegg_term_in_kegg_module.component`, `Organism_has_kegg_module.status` (closed, with `value_descriptions`), plus numeric entries for `steps_total` / `steps_present` / `steps_present_via_ec` (`min_value: 0`), following `Gene_has_tcdb_family.tcdb_gate` / `.tier` as templates.
- Modify: `create_knowledge_graph.py` — after the BRITE block (~line 276): instantiate `MultiKeggModuleAdapter(... ko_node_ids=kegg_anno_adapter.all_ko_ids(), pathway_node_ids=kegg_anno_adapter.all_pathway_ids())`, write nodes/edges; next to `chebi_release.txt` (~line 469) write `kegg_release.txt` from `kegg_data.json["kegg_release"]` (warn and skip when absent).
- Modify: `scripts/post-import.sh` and `scripts/post-import.cypher` (identical logic): indexes `kegg_module_id_idx`, `kegg_module_level_idx`, `kegg_module_level_kind_idx`, `kegg_module_type_idx`, `kegg_module_module_id_idx` + full-text `keggModuleFullText` on `name` (index group); rollup in the `IN TRANSACTIONS` group modelled on the `MeropsFamily` block (`post-import.cypher` ~1155): subtree walk `(m)<-[:Kegg_module_is_a_kegg_module*0..]-(:KeggModule)<-[:Kegg_term_in_kegg_module]-(:KeggTerm)<-[:Gene_has_kegg_ko]-(g:Gene)` → `gene_count`, `organism_count` (distinct `g.organism_name`), `direct_gene_count` only `WHERE m.level = 3`; `complete_organism_count` = distinct organisms with `status = 'complete'` edges to any module in the subtree; `Schema_info.kegg_release` from `kegg_release.txt` exactly like `chebi_release` (`post-import.sh` ~2113, `.cypher` ~2040/2063).
- Test: `tests/test_controlled_vocab.py` (existing suite must pass). There is no general `.sh`/`.cypher` parity test, so `tests/test_kegg_module_schema.py` also asserts that every new index name and the `kegg_release` stamp appear in both post-import files (pattern: `tests/test_post_import_metabolite_identity.py`).

- [ ] **Step 1: Write failing test** — extend the vocabulary tests' expectations if they pin entry counts; add `tests/test_kegg_module_schema.py` asserting the schema declares the node, its 4 post-import properties and the 4 edges with the exact label strings from Task 3.
- [ ] **Step 2: Run** `uv run pytest tests/test_kegg_module_schema.py tests/test_controlled_vocab.py -v` — FAIL.
- [ ] **Step 3: Implement** the four files.
- [ ] **Step 4: Run** `uv run pytest -m "not slow and not kg" -q` (one heavy job at a time; do not run during a Docker build) — PASS. Then `uv run python create_knowledge_graph.py --test` completes and its output dir contains the `KeggModule` and `Organism_has_kegg_module` CSVs and `kegg_release.txt`.
- [ ] **Step 5: Commit** — `feat(kegg): KeggModule schema, vocabularies, build wiring and post-import rollups`.

---

### Task 5: Rebuild, KG tests, cross-check, docs

**Files:**
- Create: `tests/kg_validity/test_kegg_modules.py`
- Modify: `tests/kg_validity/test_schema_info.py` (`kegg_release` stamped, = `kegg_data.json`'s)
- Create: `docs/kg-changes/kegg-modules.md`
- Modify: `docs/kg-changes/vocabulary-contract.md` (R3 exception #2), `CLAUDE.md` (R3 sentence "One documented exception" → two; Key Adapters bullet; node + relationship lists; Key graph facts bullet; step 6 description), `CHANGELOG.md` (`### Added` + `### Highlights`), `plans/backlog.md` (mark done).

- [ ] **Step 1: Write the KG tests** (`@pytest.mark.kg`, `run_query` fixture):
  - counts within ±5% of the Task 3 dry run for module nodes, class nodes (57), KO edges, organism edges;
  - invariants on every organism edge: `size(missing_step_indices) = steps_total - steps_present`, `size(missing_step_definitions) = size(missing_step_indices)`, `coalesce(size(ec_genes),0) = coalesce(size(ec_filled_step_indices),0) = steps_present_via_ec`, `status` matches the missing count, `steps_total = m.step_count`;
  - biology: M00161 and M00163 `complete` in every organism whose `preferred_name` starts with Prochlorococcus / Synechococcus / Parasynechococcus / Cyanobium / Thermosynechococcus; not `complete` (no edge or another status) in every other genome strain; the spec §2 TCA / nitrate / glycolysis table for the six named strains; NATL2A M00002 has `ec_genes = ['PMN2A_0295']`;
  - one source of KOs: for edges to modules whose `definition` has no `+` and no nested group, `steps_present - steps_present_via_ec` = count of distinct `r.steps` values reached by `(o)<-[:Gene_belongs_to_organism]-(:Gene)-[:Gene_has_kegg_ko]->(:KeggTerm)-[r:Kegg_term_in_kegg_module]->(m)`; ≤ for all other edges;
  - no `KeggModule` node without `module_type`; every `level = 3` node except `class_path == []` ones has exactly one `is_a` parent.
- [ ] **Step 2: Rebuild.** Snapshot counts (`logs/pre_rebuild_counts_<date>.txt`, query as in the 2026-11-18 rebuild) and `tests/kg_validity/capture_annotation_state.py --save`; `docker compose stop app deploy`; `docker compose up --build build import post-process`; then `docker compose up -d --no-deps deploy app`.
- [ ] **Step 3: Verify.** `uv run pytest -m kg -q` — all pass; counts diff vs snapshot shows only the new labels/edge types and the new `ControlledVocabulary` entries; `capture_annotation_state.py --compare` shows no change; `/omics-edge-snapshot` unchanged.
- [ ] **Step 4: KEGG Mapper cross-check for MED4.** Write MED4's KO list (`MATCH (o:OrganismTaxon {preferred_name:'Prochlorococcus MED4'})<-[:Gene_belongs_to_organism]-(:Gene)-[:Gene_has_kegg_ko]->(k) RETURN DISTINCT k.id`) to a file, paste it at `kegg.jp/kegg/mapper/reconstruct.html` (manual; ask the user if no browser is available), compare module statuses; record agreement and every disagreement (expected only where EC fills a step) in `docs/kg-changes/kegg-modules.md`. If Mapper counts `--` steps, switch the rule in Task 1 to match and rerun Tasks 3–5.
- [ ] **Step 5: Docs** — `kegg-modules.md` covers: graph shape, counts, `kegg.module:` prefix use with the note that class ids are house-minted and fail bioregistry's `^M\d+$`, the EC fallback rules, `organism_count` vs edge count, `missing_kos` / `ec_genes` semantics (spec §13), Mapper cross-check result, `Schema_info.kegg_release`. Plus the CLAUDE.md / vocabulary-contract / CHANGELOG / backlog edits listed above.
- [ ] **Step 6: Commit** — `test+docs: KEGG modules KG checks and docs`.
