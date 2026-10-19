# TCDB diamond best-hit details — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish each diamond-called gene's best-hit TCDB system(s), their names, top-hit identity, and those systems' own substrates (ids + names) as six properties on the existing `Gene_has_tcdb_family` edge.

**Architecture:** The rule runs in the Phase-1 runner (`multiomics_kg/utils/tcdb_diamond.py`) and writes a per-protein `best_hit` block into the committed `<strain>.tcdb.calls.json`, rebuilt from the committed `.tcdb.tsv` without rerunning diamond (`--from-tsv`); the step-2 merge carries the light id list `tcdb_best_hit_ids` per gene. The TCDB name scraper widens its scope to the best-hit systems and adds their names without touching existing entries. Step 6 resolves the best-hit systems' substrate strings to canonical ids and writes `best_hit_substrates` into `tcdb_pruned.json`. `tcdb_adapter` attaches the six properties.

**Tech Stack:** Python 3, pytest, BioCypher, diamond, Neo4j (Docker, built in the other clone).

**Spec:** `docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md`

## Global Constraints

- Gate: qcov ≥ 75 **and** scov ≥ 75; TC class ∈ {1,2,3,4}; reciprocal best hit over **gated** rows; window bitscore ≥ 0.9 × top, up to 3 distinct 5-part TCIDs, **any family**.
- Ordering / ties: hits by bitscore desc, then TCID asc; reciprocal reverse lookup by bitscore desc, then query id asc.
- Subject TCIDs longer than 5 parts are truncated to 5 (`lcl|Q39253-1.A.1.9.2.3` → `1.A.1.9.2`) — **only** inside the best-hit path; tier/consensus parsing is unchanged.
- Attach: the edge whose call is in the top hit's 3-part family (measured: always exists when the protein has any call); defensive fallback `calls[0]`, logged. The 188 best-hit proteins with **no** diamond call (HSP < 50 aa) get no edge record — their best hit lives in calls.json and the merge field `tcdb_best_hit_ids` only.
- Substrates are unfiltered (H⁺, Na⁺, classes, polymers stay) and kept whether or not a node exists; unresolved ChEBI → its primary `chebi:` CURIE.
- `best_hit_names` falls back to the TCID itself when unnamed.
- Diamond is not rerun; `.tcdb.tsv` files stay byte-identical; every strain's `calls[]` stays byte-identical. Existing edges, nodes and properties unchanged.
- Committed files expected to change: 43 × `.tcdb.calls.json` + `.skill_summary.json`, 43 × `gene_annotations_merged.json` and 43 × `gene_annotations_wide.json` (new key only), `tcdb_names.json` (new keys only), `tcdb_hierarchy.json` (non-kept `tc_specificity` entries gain `name`/`description`; kept entries byte-identical — step 6 copies names from `tcdb_names.json` into the hierarchy), `tcdb_pruned.json` (new key only), `cache/data/chebi/chebi_metabolite_xrefs.json` (more terms). Nothing else. KG-side: the vocabulary hash, six new `ControlledVocabulary` nodes, and `data_source:tcdb_diamond.info_types` gaining `tcdb_best_hit_ids` (derived from the config) move; nothing else.
- Rule version string: `best_hit_v1`.
- Every string element of the five list props passes `_clean_str` (no `'` or `|`); no element is ever empty (system name → TCID fallback, substrate name → id fallback).
- Stage only explicit paths in commits (concurrent sessions edit this clone).

## Review Focus

- A protein whose gated hits exist but whose top hit is claimed by a paralog → no `best_hit` block at all (not an empty one). Pinned in Task 1.
- A WP_ accession shared by two genes → both genes' edges carry the record. Pinned in Task 5.
- Two edges merged onto one (gene, TC) by `seed_aliases` remapping where only the non-base side carries `best_hit_*` → props survive the merge. Pinned in Task 5.
- A best-hit system with no substrates in TCDB → `best_hit_tcids`/`best_hit_names` present, all three substrate lists absent. Pinned in Task 5.
- Step 6 extra resolution must not change `kegg_data.json`, `metabolite_id_mapping.json`, the identity report, or any existing `tcdb_pruned.json` key; it must not share the main `rekeys` dict. Pinned in Task 4 (Step 6 byte-compare; the helper returns its own rekeys and `main` passes the original object to `Step6Inputs`).
- Reciprocal reverse lookup is keyed by the raw `subject_id` string (one TCDB accession can sit under several TCIDs).
- A best-hit-only protein (no calls) must not gain `tcdb_ids`, a `transporter_classification` entry, an edge, or a `contributing_sources` change. Pinned in Task 3.

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
  - `attach_tcid(best: dict, candidates: list[dict]) -> str | None` — the candidate `tcid` whose 3-part family equals `best["top_tcid"]`'s (calls are one per family, so at most one); else (defensive, never expected) `candidates[0]["tcid"]`; `None` if no candidates.

- [ ] **Step 1: Write failing tests** (rows built with a small `_row(q, tcid, ident, qcov, scov, bits, acc="ACC", length=200, evalue=1e-20)` helper producing `lcl|<acc>-<tcid>` subjects; the reciprocal lookup groups by raw `subject_id`):
  - `test_parse_subject_truncates_six_part_only_when_asked`: default → `None` for `lcl|Q39253-1.A.1.9.2.3`; `truncate_to=5` → `("Q39253", "1.A.1.9.2")`.
  - `test_best_hit_gate_requires_both_coverages`: hit with qcov 80/scov 70 → protein absent.
  - `test_best_hit_drops_classes_5_to_9`: only hit `5.A.1.1.1` → absent; `8.A...` → absent.
  - `test_best_hit_reciprocal_fails_for_paralog`: P1 (bits 300) and P2 (bits 400) both hit `S1-3.A.1.16.2` → P1 absent, P2 present.
  - `test_best_hit_reciprocal_ignores_ungated_competitor`: P2 hits S1 with bits 400 but scov 50 → P1 present.
  - `test_best_hit_window_cross_family_top3`: hits bits 100 `3.D.3.1.1`, 95 `3.E.2.1.1`, 92 `3.D.3.2.1`, 91 `3.D.3.3.1`, 80 `x` → `tcids == ["3.D.3.1.1","3.E.2.1.1","3.D.3.2.1"]`, `top_tcid == "3.D.3.1.1"`, `identity` = top row's.
  - `test_best_hit_duplicate_tcid_keeps_best_row`: two rows same TCID → appears once.
  - `test_best_hit_ties_deterministic`: equal bitscores → TCID ascending.
  - `test_attach_prefers_top_family_else_first_call`: candidates `[{"tcid":"2.A.1"},{"tcid":"3.A.1.16"}]`, top `3.A.1.16.2` → `"3.A.1.16"`; top `1.A.11.2.5` → `"2.A.1"`; `[]` → `None`.
  - `test_best_hit_includes_short_hsp_proteins`: a gated, reciprocal hit with length 29 is selected (the tier policy's 50-aa floor does not apply to the best-hit rule).
- [ ] **Step 2:** `pytest tests/test_tcdb_diamond.py -v` → new tests FAIL (import errors).
- [ ] **Step 3:** Implement the three functions per the Interfaces block and Global Constraints.
- [ ] **Step 4:** `pytest tests/test_tcdb_diamond.py -v` → all PASS (old tests included).
- [ ] **Step 5:** Commit `feat(tcdb-diamond): best-hit selection rule (best_hit_v1)`.

### Task 2: Runner writes `best_hit` (`--from-tsv`); regenerate all strains

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py` (`build_strain_calls`)
- Modify: `.claude/skills/tcdb-diamond/run_tcdb_diamond.py` (`--from-tsv` flag in `main` / `process_strain`)
- Modify: `.claude/skills/tcdb-diamond/SKILL.md` (Output Schema + Phase-2 notes)
- Test: `tests/test_tcdb_diamond.py`
- Data: `cache/data/*/genomes/*/tcdb/<strain>.tcdb.{calls.json,skill_summary.json}` (all 43; `.tsv` untouched)

**Interfaces:**
- Consumes: Task 1.
- Produces: `calls.json` record shape `{"calls": [...], "best_hit": {"rule": "best_hit_v1", "tcids": [...], "identity": float, "attach_tcid": str}}`. A protein with a best hit but no call gets `{"calls": [], "best_hit": {... no attach_tcid}}`. Existing `summary` keys keep their current meaning (computed over records with calls only); new keys `proteins_with_best_hit`, `best_hit_multi_tcid`, `best_hit_without_call`, `best_hit_fallback_attach`.
  - `run_tcdb_diamond.py --from-tsv [--strains ...]`: skips `resolve_tcdb_data_dir` / `build_tcdb_db` (no FASTA download, no `diamond makedb`), does not need `protein.faa`, ignores the calls.json-exists `SKIP_EXISTS` guard (always rewrites); for each strain reads the existing `<strain>.tcdb.tsv`, runs `build_strain_calls`, writes calls.json + summary. Statuses `OK` / `SKIP_NO_TSV`. `--from-tsv` with `--limit` is rejected.

- [ ] **Step 1: Failing test** `test_build_strain_calls_adds_best_hit_block` on a tmp TSV (cynD-like: top `3.A.1.16.2` bits 316 qcov 95 scov 94, plus `3.A.1.16.1` bits 246) → `calls[q]["best_hit"] == {"rule": "best_hit_v1", "tcids": ["3.A.1.16.2"], "identity": ..., "attach_tcid": calls[q]["calls"][0]["tcid"]}`; and the `calls` list equals what the function produced before the change (assert on tier/tcid). Plus `test_build_strain_calls_best_hit_without_call`: a protein whose only hit is 40 aa long (no tier call) gets `{"calls": [], "best_hit": {...}}` without `attach_tcid`, and `summary["proteins_with_call"]` does not count it. Plus `test_build_strain_calls_reciprocal_fail_has_no_best_hit_key`: a protein with a tier call whose top gated hit is claimed by a paralog → its record has no `best_hit` key.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** In `build_strain_calls`, keep all parsed rows (not only tier-classified) for `select_best_hits`, then add the block and summary counts. Add `--from-tsv` to the runner. Confirm `load_calls_json` / `iter_candidates` / `best_call` / `call_tc_families` / `call_tcids` / adapter `_load_diamond_evidence` / `load_tcdb_diamond` tolerate `"calls": []` records and only read `calls` (add a test for any that does not).
- [ ] **Step 4:** Run → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** Snapshot old calls + summaries into a scratch directory (outside the repo), then `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --from-tsv`.
- [ ] **Step 6: Verify** — `git diff --stat -- 'cache/data/*/genomes/*/tcdb/*.tsv'` empty; a script comparing each strain's new `{q: rec["calls"] for rec with calls}` to the old file → identical in all 43 strains, and no old protein lost. Report totals: proteins with best hit (expect 12,379), without call (expect 188), multi-TCID, fallback attach (expect 0), MED4 PMM0372 (look up its `protein_id` in MED4's merged annotations first; calls.json is keyed by WP_) → `["3.A.1.16.2"]`.
- [ ] **Step 7:** Update SKILL.md Output Schema with the `best_hit` block (rule text = Global Constraints gate/window). Document `--from-tsv`. Commit code + SKILL.md + the 43×2 artifacts: `feat(tcdb-diamond): emit best_hit block (--from-tsv)`.

### Task 3: Step-2 merge field `tcdb_best_hit_ids`

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py` (`load_tcdb_diamond`)
- Modify: `config/gene_annotations_config.yaml` (new passthrough field next to `tcdb_diamond_ids`)
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`_gene_reachable_sets`)
- Test: `tests/test_build_gene_annotations.py` (uses `AnnotationBuilder`, `build_merged`, `_compute_contributing_sources`), `tests/test_build_kegg_metabolism_xrefs.py`
- Data: every strain's `gene_annotations_merged.json` and `gene_annotations_wide.json`

**Interfaces:**
- Consumes: Task 2 calls.json `best_hit`.
- Produces:
  - `load_tcdb_diamond` returns `{wp: {"tcdb_ids": [...], "tcdb_best_hit_ids": [...]}}` (each key only when non-empty); a protein with a `best_hit` but no `calls` is included.
  - Merged gene field `tcdb_best_hit_ids: str[]` (`type: passthrough`, `source: tcdb_diamond`, `field: tcdb_best_hit_ids`), in window order. Sparse. Not added to `transporter_classification` (it must not seed TCDB pruning or create edges).
  - `_gene_reachable_sets(...)["tcdb_best_hit_ids"]: set[str]` — union over all strains' merged genes.

- [ ] **Step 1: Failing tests** — `test_load_tcdb_diamond_reads_best_hit_ids` (tmp calls.json with one protein carrying calls + best_hit, one carrying only best_hit: the latter gets `tcdb_best_hit_ids` and **no** `tcdb_ids`); `test_best_hit_only_gene_unchanged_otherwise` (merge a fixture gene whose only tcdb_diamond data is a best hit → no `tcdb_diamond_ids`, no `transporter_classification`, `contributing_sources` same as without the calls.json); and `test_gene_reachable_sets_collects_best_hit_ids` (merged fixture gene with `tcdb_best_hit_ids`; assert they appear in `tcdb_best_hit_ids` and **not** in `tcdb_ids`).
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement per Interfaces; config comment states the field is a light id list (identity / attach stay in calls.json, adapter reads them there) and is deliberately outside `transporter_classification`.
- [ ] **Step 4:** Run → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 0 (before any code change in this task):** baseline — `bash scripts/prepare_data.sh --steps 2 --force` on the Task-2 state; if any merged/wide file changes, stop and report the pre-existing drift to the user (do not commit it as part of this work); restore with `git checkout` on those files.
- [ ] **Step 5:** `bash scripts/prepare_data.sh --steps 2 --force`. Verify on MED4 (PMM0372 by locus tag in the merged file): PMM0372 has `tcdb_best_hit_ids == ["3.A.1.16.2"]`; a script comparing every strain's merged genes old vs new shows **only** the added key on any gene (save copies to a scratch directory first) — in particular `contributing_sources` unchanged for the 188 call-less genes; the wide files differ only by the new `tcdb_diamond_tcdb_best_hit_ids` key (and the call-less proteins' new rows carrying only it).
- [ ] **Step 6:** Commit code, config, tests, merged + wide files: `feat(merge): tcdb_best_hit_ids from the tcdb_diamond source`.

### Task 4: TCDB system names + step-6 best-hit substrates

**Files:**
- Modify: `multiomics_kg/download/scrape_tcdb_names.py` (merge-preserving `assemble_names`, `main` scope)
- Modify: `.claude/skills/tcdb-diamond/SKILL.md` ("Refreshing TCDB names (T6)": best-hit scope, top-up now merges)
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`main`)
- Modify: `multiomics_kg/utils/metabolite_identity/pipeline.py` (`Step6Inputs`, `run_identity_stages`)
- Test: `tests/test_scrape_tcdb_names.py`, `tests/test_build_kegg_metabolism_xrefs.py` (has an MNX sqlite fixture), `tests/test_step6_identity_pipeline.py` (+ `tests/step6_fixture.py` if needed)
- Data: `cache/data/tcdb/tcdb_names.json`, `cache/data/tcdb/tcdb_hierarchy.json`, `cache/data/tcdb/tcdb_pruned.json`, `cache/data/chebi/chebi_metabolite_xrefs.json`

**Interfaces:**
- Consumes: Task 3 merged `tcdb_best_hit_ids` (via `_gene_reachable_sets`).
- Produces:
  - `collect_best_hit_tcids(genome_dirs: Iterable[Path]) -> set[str]` in `scrape_tcdb_names.py` (not `tcdb_diamond.py`, whose contract is pure and never reads merged files): union of merged `tcdb_best_hit_ids`, reading each dir with `load_gene_annotations`; `main` gets the dirs from `load_genome_rows()` (`multiomics_kg/download/utils/cli.py`), as step 6 does. Step 6 itself uses `_gene_reachable_sets`.
  - `assemble_names(browse_names, per_family_results, wanted_ids, scraped_families, failed_families, existing: dict | None = None) -> dict`: starts from `existing["names"]` verbatim; adds an entry for each parsed 4/5-part id in `wanted_ids` **not already present** (never overwrites); browse-level names likewise added only when absent; `meta.scraped_families` = union with existing; `meta.scraped_at` = today (it records the last scrape). With `existing=None` the current behaviour is kept (fresh build).
  - `main` reads the current `tcdb_names.json` as `existing`; `wanted_ids = kept_ids | collect_best_hit_tcids(genome dirs)`; scope = `scoped_families(wanted_ids)`; `--families` filters that scope (fixes the pre-existing bug where a top-up rewrote the artifact from the listed families only, and where families with no kept id could never be selected).
  - `Step6Inputs.best_hit_substrates: dict[str, list[tuple[str, str]]] = field(default_factory=dict)` (placed after the other defaulted fields) — `{tcid5: [(primary_or_role_id, tcdb_name), ...]}`, pre-canonical.
  - `_best_hit_substrate_inputs(tcids: set[str], full_subtree_strings: dict[str, set[str]], conn, chebi, known_kegg: set[str]) -> tuple[dict[str, list[tuple[str, str]]], dict[str, str]]` in `build_kegg_metabolism_xrefs.py` — returns the inputs above and its own rekeys (for logging only).
  - `tcdb_pruned.json["best_hit_substrates"]: {tcid5: [[canonical_id, name], ...]}` sorted by id; `name` = final table entry's `name` when the id is a node, else `chemical_roles[id]["name"]` for roles, else the TCDB string's name part.

- [ ] **Step 1: Failing tests**
  - `test_assemble_names_adds_wanted_systems_without_touching_existing`: `existing` has `3.A.1.16.2` with an old name; a page re-parses it with new text and also lists `3.A.1.16.3` (wanted) and `3.A.1.16.9` (not wanted) → `3.A.1.16.2` keeps the old name, `3.A.1.16.3` added with `name` + `description`, `3.A.1.16.9` absent.
  - `test_collect_best_hit_tcids` on two tmp strain dirs with merged files.
  - `test_families_topup_merges_into_existing_artifact`: existing artifact with names for families F1 and F2; run `main(families=["F3"])` where F3 has only best-hit ids, against cached tmp pages (monkeypatch `PAGES_DIR`, `NAMES_FILE`, `PRUNED_FILE`, and `scrape_tcdb_names.collect_best_hit_tcids` to return the F3 ids; no network) → F1/F2 entries byte-identical, F3 systems added, `scraped_families` ⊇ {F1, F2, F3}.
  - `test_identity_stages_writes_best_hit_substrates` (fixture): an input primary that stage E merges onto another canonical id is rewritten to the canonical id and named from the table; a role id passes through named from `chemical_roles`; an id with no node keeps its id and TCDB name; a TC with no substrates is absent; all other `tcdb_pruned` keys equal a run without `best_hit_substrates` input.
  - `test_best_hit_substrate_inputs` (in `tests/test_build_kegg_metabolism_xrefs.py`, MNX fixture): a 5-part TCID's own molecule string resolves to its primary; a role string maps to its `chebi.role:` id; an unresolvable ChEBI string → `chebi:<n>`; a TCID absent from `full_subtree_strings` yields nothing (no KeyError). The no-shared-rekeys guarantee is pinned by Step 6's byte-compare of `metabolite_id_mapping.json` + identity report.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement. In step-6 `main`: take best-hit TCIDs from `catalysis_sets["tcdb_best_hit_ids"]` and call `_best_hit_substrate_inputs`: each id's own strings from `full_subtree_strings.get(tc, set())` (5-part = own substrates); split with `partition_role_strings`; resolve molecule strings with a **separate** `_resolve_substrate_strings` call using its own fresh `rekeys` dict (never merged into the main `rekeys`, `compound_props` or `kegg_data`), **tracked** (not inside `chebi.untracked()`, so the touched ChEBI terms join the committed cache and step 6 stays reproducible without raw `chebi.obo`); map each to `(primary, name)` (`_primary_for` already returns the KEGG id when it re-keys); strings that resolve to nothing → `chebi:<primary n>` via `chebi.primary`, or the raw `chebi:<n>` when `primary` returns `None`. Precondition: raw `cache/data/chebi/raw/chebi.obo.gz` present (only then is `chebi_metabolite_xrefs.json` rewritten — `chebi_from_raw`). In `run_identity_stages`, rewrite via `identity.canonical_of` and name as above. Also in step-6 `main`, log a warning with the count and a sample when best-hit TCIDs have no entry in `tcdb_names.json` (`_warn_unscraped_families` pattern; report-only).
- [ ] **Step 4:** Run tests → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** Copy the page cache from the other clone: `mkdir -p cache/data/tcdb/raw && cp -rn /home/osnat/github/multiomics_biocypher_kg/cache/data/tcdb/raw/pages cache/data/tcdb/raw/` (352 pages, gitignored). Run the scraper (≈ 53 new pages hit tcdb.org, ~2–3 min): `uv run python -m multiomics_kg.download.scrape_tcdb_names`. Check: every pre-existing `names` key byte-identical (script diff → 0 changes); share of best-hit TCIDs named ≥ 95% — if lower, report the unnamed families to the user before continuing (do not loosen the target).
- [ ] **Step 6:** Copy every step-6 output (`kegg_data.json`, `tcdb_pruned.json`, `metabolite_identity_report.json`, `metabolite_id_mapping.json`, the `brite_br08xxx.json` caches, `dominant_stereo_forms.json`, `chebi_metabolite_xrefs.json`, paper `_resolved.csv` files) to a scratch directory; `bash scripts/prepare_data.sh --steps 6 --force`. Verify with `git status`/diff: only `tcdb_pruned.json` (new `best_hit_substrates` key only, every other key identical), `tcdb_hierarchy.json` (only non-kept `tc_specificity` entries gain `name`/`description`; kept entries identical) and `chebi_metabolite_xrefs.json` (terms added, none changed) differ; `kegg_data.json`, `metabolite_id_mapping.json`, the identity report byte-identical; report the number of TCIDs with substrates, ids with no node, ids in slice 1's target filter. MED4 `3.A.1.16.2` → cyanate + nitrite.
- [ ] **Step 7:** Commit code, tests, SKILL.md, `tcdb_names.json`, `tcdb_hierarchy.json`, `tcdb_pruned.json`, `chebi_metabolite_xrefs.json`: `feat(step6): best-hit system names and substrates`.

### Task 5: Adapter, schema, vocabulary

**Files:**
- Modify: `multiomics_kg/adapters/tcdb_adapter.py`
- Modify: `config/schema_config.yaml` (`gene_has_tcdb_family` properties)
- Modify: `config/controlled_vocabularies.yaml`
- Test: `tests/test_tcdb_adapter.py`

**Interfaces:**
- Consumes: calls.json `best_hit` (Task 2); `tcdb_pruned.json["best_hit_substrates"]` and the system names that step 6 copied into `tcdb_hierarchy.json` (Task 4) — the Multi adapter already loads `self._hierarchy`, so names come from `self._hierarchy.get(tc, {}).get("name")`, no separate read of `tcdb_names.json`.
- Produces edge props (sparse, diamond-sourced edges only): `best_hit_tcids: str[]`, `best_hit_names: str[]` (aligned), `best_hit_identity: float`, `best_hit_substrate_tcids: str[]`, `best_hit_substrate_ids: str[]`, `best_hit_substrate_names: str[]` — long form, one entry per (window system, substrate), the three lists aligned; ordered by window rank then substrate id; a shared substrate repeats once per system.

- [ ] **Step 1: Failing tests**
  - `test_strain_edge_carries_best_hit_on_attach_edge_only`: protein with calls `2.A.1` and `3.A.1.16`, `attach_tcid = "3.A.1.16"` → only that edge has `best_hit_tcids` / `best_hit_identity`.
  - `test_shared_wp_both_genes_get_best_hit`.
  - `test_multi_adds_names_and_substrates`: names from the hierarchy fixture, fallback to the TCID when missing (existing fixtures without names must keep passing); substrates from `best_hit_substrates` in long form — a window of `[A, B]` where A has {x, y} and B has {y} gives `substrate_tcids == [A, A, B]`, `substrate_ids == [x, y, y]`, names aligned; strings sanitized.
  - `test_best_hit_without_substrates_omits_substrate_props`.
  - `test_best_hit_without_call_emits_nothing`: a calls.json record `{"calls": [], "best_hit": {...}}` → no edge, no error.
  - `test_best_hit_strings_sanitized_and_nonempty`: a substrate named `N,N'-diacetylchitobiose` → `N,N^-diacetylchitobiose`; an unnamed system → name == TCID.
  - `test_alias_merge_keeps_best_hit_from_non_base_side`: two edges collapsing to one (gene, TC) where the base side lacks `best_hit_*` → merged props carry them.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Strain adapter: load `best_hit` in `_load_diamond_evidence`'s pass (keyed by WP_); on the diamond-sourced edge whose source `tc == best_hit["attach_tcid"]` add `best_hit_tcids`, `best_hit_identity`. Multi adapter: read `best_hit_substrates` from the pruned file in `download_data` (optional key → `{}`); in the gene-edge loop add the four derived props (`best_hit_names` + the three substrate lists); in the alias-merge branch copy any `best_hit_*` key missing from `combined` from the other side. Schema: add the six props. Vocabulary: six entries modeled on `Gene_has_tcdb_family.identity` (`sparse: true`, `closed: false`, `values: []`; `value_type: string_array` for the five lists, `float` with `min_value 0.0 / max_value 100.0` for `best_hit_identity`), descriptions researcher-facing (no build provenance in the text — `# provenance:` comment instead), each stating the caveat (raw best hit, long-form substrate lists, unfiltered substrates, ids may lack nodes, window may cross families).
- [ ] **Step 4:** `pytest tests/test_tcdb_adapter.py tests/test_controlled_vocab.py -v` → PASS; `pytest -m "not slow and not kg" -q` → PASS.
- [ ] **Step 5:** `uv run python create_knowledge_graph.py --test` → completes; grep the gene_has_tcdb_family CSV header for the six props.
- [ ] **Step 6:** Commit `feat(tcdb): best-hit properties on Gene_has_tcdb_family`.

### Task 6: Live rebuild, KG tests, validation, docs

**Files:**
- Create: `tests/kg_validity/test_tcdb_best_hit.py`
- Create: `scripts/tcdb_best_hit_validation.py`
- Create: `docs/kg-changes/tcdb-diamond-best-hit.md`
- Modify: `CLAUDE.md` (TCDB `Gene_has_tcdb_family` bullet), `CHANGELOG.md` (`### Added`), `plans/backlog.md` (transport section: project 2 done + follow-ups), `tests/kg_validity/snapshot_data.json` only if regenerated intentionally.

- [ ] **Step 1: Write KG tests** (`@pytest.mark.kg`): six props co-occur (all-or-none for the tcids/names/identity triple; the three substrate lists equal in length; every `best_hit_substrate_tcids` value is in `best_hit_tcids`); `best_hit_*` only on edges with `'tcdb_diamond' IN r.sources`; each gene has the record on ≤ 1 TCDB edge per protein; MED4 spot checks — PMM0372 cynD and cynA → `best_hit_tcids[0] = '3.A.1.16.2'` with a non-id name and a cyanate substrate; amt1 → an ammonium substrate; the urt subunits → a urea substrate; count of edges with `best_hit_tcids` equals the number of genes whose protein has a Task-2 `best_hit` with `attach_tcid`, minus attach edges the adapter's prune drops (target not kept, no seed alias; expected 0 — compute and print it rather than assume), computed from calls.json + merged files.
- [ ] **Step 2:** Before rebuild: `/omics-edge-snapshot` save + `MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r)` and `TcdbFamily` / `Metabolite` counts noted.
- [ ] **Step 3:** Push the branch; the user pulls and rebuilds in the other clone (`multiomics_biocypher_kg`, `docker compose up -d`). Wait for the user.
- [ ] **Step 4:** After rebuild: `pytest -m kg -v` → PASS (new + existing); `/omics-edge-snapshot` compare → no change; TCDB edge / node / Metabolite counts identical to Step 2; `ControlledVocabulary` count +6, `Schema_info.controlled_vocabularies_hash` changed and `data_source:tcdb_diamond.info_types` gained `tcdb_best_hit_ids` (expected).
- [ ] **Step 5: Validation script** (`scripts/tcdb_best_hit_validation.py`, Cypher against the live graph; rule frozen, no thresholds tuned). Per strain and corpus-wide, MED4 reported separately as training data; hold-out = all other strains. For edges with best-hit details: eggNOG TC agreement (`same | best-hit descendant | best-hit ancestor | different family | no eggNOG`), KO in BRITE `ko02000`, best-hit family's Pfam bridge ∩ gene Pfams, Cyanorak transport role (Pro/Syn), GO-transport substrate overlap (`same | disjoint | no GO`), identity bins (<40 / 40–70 / ≥70), split calls, fallback attachments (expect 0), the 188 call-less best hits (merge field only), substrate ids without a node, ids in slice 1's target filter; nitrogen sanity table (MED4, MIT9313, WH8102, one Alteromonas: N-source genes, best-hit substrates vs current tier-rollup substrate count). Also from calls.json/TSV (not Cypher): genes gaining details, reciprocal-check failures by category (product/COG keyword list fixed in the script before running). KO-name substrate overlap is deferred (no structured KO→substrate source; say so in the doc). A seeded random sample of 30 eggNOG-different-family disagreements is printed and classified by hand (rule-wrong / comparator-wrong / unclear), result recorded in the kg-changes doc.
- [ ] **Step 6:** Write `docs/kg-changes/tcdb-diamond-best-hit.md`: what's new, the rule, the six properties (long-form substrate layout), caveat paragraph for the explorer (raw clue; identity guide; several TCIDs = near-tie, check substrate agreement; substrates unfiltered, may lack nodes; window may cross families; absent substrates = TCDB lists none; 188 call-less proteins carry the best hit only in the merged annotations; not the same as MEROPS `best_hit_id`), validation table, query examples. Update CLAUDE.md, CHANGELOG, backlog, memory `project_transport_programme_2026_10_18`.
- [ ] **Step 7:** Commit tests, script, docs: `docs(tcdb): best-hit details kg-changes, validation, CLAUDE.md`.
