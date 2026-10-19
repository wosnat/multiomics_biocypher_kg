# TCDB diamond best-hit details — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish each diamond-called gene's best-hit TCDB system(s), their names, top-hit identity, and those systems' own substrates (ids + names) as six properties on the existing `Gene_has_tcdb_family` edge.

**Architecture:** The rule runs in the Phase-1 runner (`multiomics_kg/utils/tcdb_diamond.py`) and writes a per-protein `best_hit` block into the committed `<strain>.tcdb.calls.json`; the step-2 merge carries the light id list `tcdb_best_hit_ids` per gene. The TCDB name scraper keeps every system name on each family page and widens its scope to the best-hit systems. Step 6 resolves the best-hit systems' substrate strings to canonical ids and writes `best_hit_substrates` into `tcdb_pruned.json`. `tcdb_adapter` attaches the six properties.

**Tech Stack:** Python 3, pytest, BioCypher, diamond, Neo4j (Docker, built in the other clone).

**Spec:** `docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md`

## Global Constraints

- Gate: qcov ≥ 75 **and** scov ≥ 75; TC class ∈ {1,2,3,4}; reciprocal best hit over **gated** rows; window bitscore ≥ 0.9 × top, up to 3 distinct 5-part TCIDs, **any family**.
- Ordering / ties: hits by bitscore desc, then TCID asc; reciprocal reverse lookup by bitscore desc, then query id asc.
- Subject TCIDs longer than 5 parts are truncated to 5 (`lcl|Q39253-1.A.1.9.2.3` → `1.A.1.9.2`) — **only** inside the best-hit path; tier/consensus parsing is unchanged.
- Attach: the edge whose call is in the top hit's 3-part family; else the protein's highest-confidence call (`calls[0]`). Nothing dropped.
- Substrates are unfiltered (H⁺, Na⁺, classes, polymers stay) and kept whether or not a node exists; unresolved ChEBI → its primary `chebi:` CURIE.
- `best_hit_names` falls back to the TCID itself when unnamed.
- Diamond parameters unchanged (`--max-target-seqs 25 --evalue 0.001 --more-sensitive`). Existing tier calls, edges and all other properties must be unchanged.
- Rule version string: `best_hit_v1`.
- String props sanitized with `_clean_str` (no `'` or `|`).
- Stage only explicit paths in commits (concurrent sessions edit this clone).

## Review Focus

- A protein whose gated hits exist but whose top hit is claimed by a paralog → no `best_hit` block at all (not an empty one). Pinned in Task 1.
- A WP_ accession shared by two genes → both genes' edges carry the record. Pinned in Task 5.
- Two edges merged onto one (gene, TC) by `seed_aliases` remapping where only the non-base side carries `best_hit_*` → props survive the merge. Pinned in Task 5.
- A best-hit system with no substrates in TCDB → `best_hit_tcids`/`best_hit_names` present, all three substrate lists absent. Pinned in Task 5.
- Step 6 extra resolution must not change `kegg_data.json` or any existing `tcdb_pruned.json` key. Pinned in Task 4 (byte-compare).

---

### Task 1: Best-hit rule

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Test: `tests/test_tcdb_diamond.py`

**Interfaces:**
- Produces:
  - `parse_tcdb_subject_id(subject_id: str, *, truncate_to: int | None = None) -> tuple[str, str] | None` — with `truncate_to=5`, accepts 3–6+ part tails and truncates; default behaviour unchanged.
  - `BEST_HIT_RULE_VERSION = "best_hit_v1"`
  - `select_best_hits(rows: list[dict]) -> dict[str, dict]` — input = parsed diamond rows (`parse_diamond_row` dicts) for one strain; output `{query_id: {"tcids": [str, ...], "identity": float, "top_tcid": str}}`; proteins failing gate or reciprocal check are absent.
  - `attach_tcid(best: dict, candidates: list[dict]) -> str | None` — the candidate `tcid` whose 3-part family equals `best["top_tcid"]`'s (if several, the longest that is a prefix of `top_tcid`); else `candidates[0]["tcid"]`; `None` if no candidates.

- [ ] **Step 1: Write failing tests** (rows built with a small `_row(q, tcid, ident, qcov, scov, bits)` helper producing `lcl|ACC-<tcid>` subjects):
  - `test_parse_subject_truncates_six_part_only_when_asked`: default → `None` for `lcl|Q39253-1.A.1.9.2.3`; `truncate_to=5` → `("Q39253", "1.A.1.9.2")`.
  - `test_best_hit_gate_requires_both_coverages`: hit with qcov 80/scov 70 → protein absent.
  - `test_best_hit_drops_classes_5_to_9`: only hit `5.A.1.1.1` → absent; `8.A...` → absent.
  - `test_best_hit_reciprocal_fails_for_paralog`: P1 (bits 300) and P2 (bits 400) both hit `S1-3.A.1.16.2` → P1 absent, P2 present.
  - `test_best_hit_reciprocal_ignores_ungated_competitor`: P2 hits S1 with bits 400 but scov 50 → P1 present.
  - `test_best_hit_window_cross_family_top3`: hits bits 100 `3.D.3.1.1`, 95 `3.E.2.1.1`, 92 `3.D.3.2.1`, 91 `3.D.3.3.1`, 80 `x` → `tcids == ["3.D.3.1.1","3.E.2.1.1","3.D.3.2.1"]`, `top_tcid == "3.D.3.1.1"`, `identity` = top row's.
  - `test_best_hit_duplicate_tcid_keeps_best_row`: two rows same TCID → appears once.
  - `test_best_hit_ties_deterministic`: equal bitscores → TCID ascending.
  - `test_attach_prefers_top_family_else_first_call`: candidates `[{"tcid":"2.A.1"},{"tcid":"3.A.1.16"}]`, top `3.A.1.16.2` → `"3.A.1.16"`; top `1.A.11.2.5` → `"2.A.1"`.
- [ ] **Step 2:** `pytest tests/test_tcdb_diamond.py -v` → new tests FAIL (import errors).
- [ ] **Step 3:** Implement the three functions per the Interfaces block and Global Constraints.
- [ ] **Step 4:** `pytest tests/test_tcdb_diamond.py -v` → all PASS (old tests included).
- [ ] **Step 5:** Commit `feat(tcdb-diamond): best-hit selection rule (best_hit_v1)`.

### Task 2: Runner writes `best_hit`; rerun all strains

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py` (`build_strain_calls`)
- Modify: `.claude/skills/tcdb-diamond/SKILL.md` (Output Schema + Phase-2 notes)
- Test: `tests/test_tcdb_diamond.py`
- Data: `cache/data/*/genomes/*/tcdb/<strain>.tcdb.{tsv,calls.json,skill_summary.json}` (all 43)

**Interfaces:**
- Consumes: Task 1.
- Produces: `calls.json` record shape `{"calls": [...], "best_hit": {"rule": "best_hit_v1", "tcids": [...], "identity": float, "attach_tcid": str}}`. A protein with a best hit but no call still gets a record (`"calls": []`). `summary` gains `proteins_with_best_hit`, `best_hit_multi_tcid`, `best_hit_fallback_attach`.

- [ ] **Step 1: Failing test** `test_build_strain_calls_adds_best_hit_block` on a tmp TSV (cynD-like: top `3.A.1.16.2` bits 316 qcov 95 scov 94, plus `3.A.1.16.1` bits 246) → `calls[q]["best_hit"] == {"rule": "best_hit_v1", "tcids": ["3.A.1.16.2"], "identity": ..., "attach_tcid": calls[q]["calls"][0]["tcid"]}`; and the `calls` list equals what the function produced before the change (assert on tier/tcid).
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** In `build_strain_calls`, keep all parsed rows (not only tier-classified) for `select_best_hits`, then add the block and summary counts. Confirm `load_calls_json` / `iter_candidates` / `best_call` / adapter `_load_diamond_evidence` only read `calls` (no change needed; add a test if any does not).
- [ ] **Step 4:** Run → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** Snapshot old calls: copy every `cache/data/*/genomes/*/tcdb/*.tcdb.calls.json` into a scratch directory (outside the repo), then rerun: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --force`.
- [ ] **Step 6: Verify tier calls unchanged** — a short script comparing each strain's new `{q: rec["calls"]}` to the old file → 0 differences in all 43 strains. Report totals: proteins with best hit (expect ≈ 12.4K), multi-TCID, fallback attach (expect ≈ 188), MED4 PMM0372 → `["3.A.1.16.2"]`.
- [ ] **Step 7:** Update SKILL.md Output Schema with the `best_hit` block (rule text = Global Constraints gate/window). Commit code + SKILL.md + the 43×3 artifacts: `feat(tcdb-diamond): emit best_hit block; rerun all strains`.

### Task 3: Step-2 merge field `tcdb_best_hit_ids`

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py` (`load_tcdb_diamond`)
- Modify: `config/gene_annotations_config.yaml` (new passthrough field next to `tcdb_diamond_ids`)
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`_gene_reachable_sets`)
- Test: the existing `load_tcdb_diamond` / merge tests (`tests/test_tcdb_adapter.py` or the gene-annotation build tests), `tests/test_build_kegg_metabolism_xrefs.py`
- Data: every strain's `gene_annotations_merged.json`

**Interfaces:**
- Consumes: Task 2 calls.json `best_hit`.
- Produces:
  - `load_tcdb_diamond` returns `{wp: {"tcdb_ids": [...], "tcdb_best_hit_ids": [...]}}` (each key only when non-empty); a protein with a `best_hit` but no `calls` is included.
  - Merged gene field `tcdb_best_hit_ids: str[]` (`type: passthrough`, `source: tcdb_diamond`, `field: tcdb_best_hit_ids`), in window order. Sparse. Not added to `transporter_classification` (it must not seed TCDB pruning or create edges).
  - `_gene_reachable_sets(...)["tcdb_best_hit_ids"]: set[str]` — union over all strains' merged genes.

- [ ] **Step 1: Failing tests** — `test_load_tcdb_diamond_reads_best_hit_ids` (tmp calls.json with one protein carrying calls + best_hit, one carrying only best_hit) and `test_gene_reachable_sets_collects_best_hit_ids` (merged fixture gene with `tcdb_best_hit_ids`; assert they appear in `tcdb_best_hit_ids` and **not** in `tcdb_ids`).
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement per Interfaces; config comment states the field is a light id list (identity / attach stay in calls.json, adapter reads them there) and is deliberately outside `transporter_classification`.
- [ ] **Step 4:** Run → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** `bash scripts/prepare_data.sh --steps 2 --force`. Verify on MED4: PMM0372 has `tcdb_best_hit_ids == ["3.A.1.16.2"]`; a script comparing every strain's merged genes old vs new shows **only** the added key (save copies to a scratch directory first).
- [ ] **Step 6:** Commit code, config, tests, merged files: `feat(merge): tcdb_best_hit_ids from the tcdb_diamond source`.

### Task 4: TCDB system names + step-6 best-hit substrates

**Files:**
- Modify: `multiomics_kg/download/scrape_tcdb_names.py` (`assemble_names`, `main` scope)
- Modify: `multiomics_kg/utils/tcdb_diamond.py` (add collector over merged files)
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`main`)
- Modify: `multiomics_kg/utils/metabolite_identity/pipeline.py` (`Step6Inputs`, `run_identity_stages`)
- Test: `tests/test_scrape_tcdb_names.py`, `tests/test_tcdb_diamond.py`, `tests/test_step6_identity_pipeline.py` (+ `tests/step6_fixture.py` if needed)
- Data: `cache/data/tcdb/tcdb_names.json`, `cache/data/tcdb/tcdb_pruned.json`

**Interfaces:**
- Consumes: Task 3 merged `tcdb_best_hit_ids` (via `_gene_reachable_sets`).
- Produces:
  - `collect_best_hit_tcids(genome_dirs: Iterable[Path]) -> set[str]` in `tcdb_diamond.py`: union of merged `tcdb_best_hit_ids` over each `gene_annotations_merged.json` (used by the scraper; step 6 uses `_gene_reachable_sets`).
  - `assemble_names(...)` keeps **every** subfamily/system parsed on scraped family pages (drop the `kept_ids` filter for 4/5-part entries). `main` scope = `scoped_families(kept_ids | collect_best_hit_tcids(...))`. Top-up runs use the full scope; cached pages are not refetched, so only missing families hit the network. (Do **not** use `--families` for this: it rewrites the artifact from the subset.)
  - `Step6Inputs.best_hit_substrates: dict[str, list[tuple[str, str]]] = {}` — `{tcid5: [(primary_or_role_id, tcdb_name), ...]}`, pre-canonical.
  - `tcdb_pruned.json["best_hit_substrates"]: {tcid5: [[canonical_id, name], ...]}` sorted by id; `name` = final table entry's `name` when the id is a node, else `chemical_roles[id]["name"]` for roles, else the TCDB string's name part.

- [ ] **Step 1: Failing tests**
  - `test_assemble_names_keeps_unkept_systems`: a system id absent from `kept_ids` appears in `names` with `name` and `description`.
  - `test_collect_best_hit_tcids` on two tmp strain dirs with merged files.
  - `test_identity_stages_writes_best_hit_substrates` (fixture): an input primary that stage E merges onto another canonical id is rewritten to the canonical id and named from the table; a role id passes through named from `chemical_roles`; an id with no node keeps its id and TCDB name; a TC with no substrates is absent.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement. In step-6 `main`: take best-hit TCIDs from `catalysis_sets["tcdb_best_hit_ids"]`; take each id's own strings from `full_subtree_strings[tc]` (5-part = own substrates); split with `partition_role_strings`; resolve molecule strings with a **separate** `_resolve_substrate_strings` call using its own `rekeys` dict (never merged into the main `rekeys`, `compound_props` or `kegg_data`); map each to `(rekeys_bh.get(p, p), name)`; strings that resolve to nothing → `chebi:<primary n>` via `chebi.primary`. In `run_identity_stages`, rewrite via `identity.canonical_of` and name as above. Also in step-6 `main`, log a warning with the count and a sample when best-hit TCIDs have no entry in `tcdb_names.json` (`_warn_unscraped_families` pattern; report-only).
- [ ] **Step 4:** Run tests → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** Run the scraper (network, ~10 min): `uv run python -m multiomics_kg.download.scrape_tcdb_names`. Check: share of best-hit TCIDs named ≥ 95%; old names unchanged for previously named ids (diff the `names` dict on old keys → 0 changes).
- [ ] **Step 6:** Copy `kegg_data.json`, `tcdb_pruned.json`, `metabolite_identity_report.json` to the scratchpad; `bash scripts/prepare_data.sh --steps 6 --force`. Verify: `kegg_data.json` byte-identical; `tcdb_pruned.json` identical except the new `best_hit_substrates` key; report the number of TCIDs with substrates, ids with no node, ids in slice 1's target filter. MED4 `3.A.1.16.2` → cyanate + nitrite.
- [ ] **Step 7:** Commit code, tests, `tcdb_names.json`, `tcdb_pruned.json` (and any other step-6 output only if it changed and the change is explained): `feat(step6): best-hit system names and substrates`.

### Task 5: Adapter, schema, vocabulary

**Files:**
- Modify: `multiomics_kg/adapters/tcdb_adapter.py`
- Modify: `config/schema_config.yaml` (`gene_has_tcdb_family` properties)
- Modify: `config/controlled_vocabularies.yaml`
- Test: `tests/test_tcdb_adapter.py`

**Interfaces:**
- Consumes: calls.json `best_hit` (Task 2); `tcdb_pruned.json["best_hit_substrates"]`, `tcdb_names.json` (Task 4).
- Produces edge props (sparse, diamond-sourced edges only): `best_hit_tcids: str[]`, `best_hit_names: str[]` (aligned), `best_hit_identity: float`, `best_hit_substrate_tcids: str[]`, `best_hit_substrate_ids: str[]`, `best_hit_substrate_names: str[]` — long form, one entry per (window system, substrate), the three lists aligned; ordered by window rank then substrate id; a shared substrate repeats once per system.

- [ ] **Step 1: Failing tests**
  - `test_strain_edge_carries_best_hit_on_attach_edge_only`: protein with calls `2.A.1` and `3.A.1.16`, `attach_tcid = "3.A.1.16"` → only that edge has `best_hit_tcids` / `best_hit_identity`.
  - `test_shared_wp_both_genes_get_best_hit`.
  - `test_multi_adds_names_and_substrates`: names from `tcdb_names.json`, fallback to the TCID when missing; substrates from `best_hit_substrates` in long form — a window of `[A, B]` where A has {x, y} and B has {y} gives `substrate_tcids == [A, A, B]`, `substrate_ids == [x, y, y]`, names aligned; strings sanitized.
  - `test_best_hit_without_substrates_omits_substrate_props`.
  - `test_alias_merge_keeps_best_hit_from_non_base_side`: two edges collapsing to one (gene, TC) where the base side lacks `best_hit_*` → merged props carry them.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Strain adapter: load `best_hit` in `_load_diamond_evidence`'s pass (keyed by WP_); on the diamond-sourced edge whose source `tc == best_hit["attach_tcid"]` add `best_hit_tcids`, `best_hit_identity`. Multi adapter: load `tcdb_names.json` and `best_hit_substrates` in `download_data`; in the gene-edge loop add the four derived props (`best_hit_names` + the three substrate lists); in the alias-merge branch copy any `best_hit_*` key missing from `combined` from the other side. Schema: add the six props. Vocabulary: six entries modeled on `Gene_has_tcdb_family.identity` (`sparse: true`, `closed: false`; identity `min_value 0.0 / max_value 100.0`), descriptions researcher-facing (no build provenance in the text — `# provenance:` comment instead), each stating the caveat (raw best hit, long-form substrate lists, unfiltered substrates, ids may lack nodes, window may cross families).
- [ ] **Step 4:** `pytest tests/test_tcdb_adapter.py tests/test_controlled_vocab.py -v` → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** `uv run python create_knowledge_graph.py --test` → completes; grep the gene_has_tcdb_family CSV header for the six props.
- [ ] **Step 6:** Commit `feat(tcdb): best-hit properties on Gene_has_tcdb_family`.

### Task 6: Live rebuild, KG tests, validation, docs

**Files:**
- Create: `tests/kg_validity/test_tcdb_best_hit.py`
- Create: `scripts/tcdb_best_hit_validation.py`
- Create: `docs/kg-changes/tcdb-diamond-best-hit.md`
- Modify: `CLAUDE.md` (TCDB `Gene_has_tcdb_family` bullet), `CHANGELOG.md` (`### Added`), `plans/backlog.md` (transport section: project 2 done + follow-ups), `tests/kg_validity/snapshot_data.json` only if regenerated intentionally.

- [ ] **Step 1: Write KG tests** (`@pytest.mark.kg`): six props co-occur (all-or-none for the tcids/names/identity triple; the three substrate lists equal in length; every `best_hit_substrate_tcids` value is in `best_hit_tcids`); `best_hit_*` only on edges with `'tcdb_diamond' IN r.sources`; each gene has the record on ≤ 1 TCDB edge per protein; MED4 PMM0372 edge → `best_hit_tcids[0] = '3.A.1.16.2'` with a non-id name and a cyanate substrate; total edges with `best_hit_tcids` within 5% of the Task-2 protein count summed over genes.
- [ ] **Step 2:** Before rebuild: `/omics-edge-snapshot` save + `MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r)` and `TcdbFamily` / `Metabolite` counts noted.
- [ ] **Step 3:** Push the branch; the user pulls and rebuilds in the other clone (`multiomics_biocypher_kg`, `docker compose up -d`). Wait for the user.
- [ ] **Step 4:** After rebuild: `pytest -m kg -v` → PASS (new + existing); `/omics-edge-snapshot` compare → no change; TCDB edge / node / Metabolite counts identical to Step 2.
- [ ] **Step 5: Validation script** (`scripts/tcdb_best_hit_validation.py`, Cypher against the live graph; rule frozen, no thresholds tuned). Per strain and corpus-wide, MED4 reported separately as training data; hold-out = all other strains. For edges with best-hit details: eggNOG TC agreement (`same | best-hit descendant | best-hit ancestor | different family | no eggNOG`), KO in BRITE `ko02000`, best-hit family's Pfam bridge ∩ gene Pfams, Cyanorak transport role (Pro/Syn), GO-transport substrate overlap (`same | disjoint | no GO`), identity bins (<40 / 40–70 / ≥70), split calls, fallback attachments, substrate ids without a node, ids in slice 1's target filter; nitrogen sanity table (MED4, MIT9313, WH8102, one Alteromonas: N-source genes, best-hit substrates vs current tier-rollup substrate count). Also a seeded random sample of 30 eggNOG-different-family disagreements printed for manual review.
- [ ] **Step 6:** Write `docs/kg-changes/tcdb-diamond-best-hit.md`: what's new, the rule, the six properties (long-form substrate layout), caveat paragraph for the explorer (raw clue; identity guide; several TCIDs = near-tie, check substrate agreement; substrates unfiltered, may lack nodes; window may cross families; absent substrates = TCDB lists none), validation table, query examples. Update CLAUDE.md, CHANGELOG, backlog, memory `project_transport_programme_2026_10_18`.
- [ ] **Step 7:** Commit tests, script, docs: `docs(tcdb): best-hit details kg-changes, validation, CLAUDE.md`.
