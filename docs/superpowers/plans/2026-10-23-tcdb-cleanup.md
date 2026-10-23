# TCDB Cleanup (project 3a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the TCDB slice of the KG faithful to TCDB: substrates only on `tc_system` nodes, one gene edge per source (`tcdb_call`) stamped with `tcdb_gate`, best-hit systems as real nodes and edges, and every roll-up/workaround property removed.

**Architecture:** The change runs through the existing pipeline order. (1) The diamond artifact gains per-system best-hit statistics. (2) The step-6 identity pipeline swaps the subtree roll-up for per-system substrates. (3) Step 6 widens pruning (best-hit seeds plus substrate-bearing systems). (4) The names scraper widens scope. (5) The TCDB adapter is rewritten to emit the new nodes and edges. (6) Schema, vocabulary and post-import are rewired. (7) Tests follow, then the rebuild and docs. No new modules except one pure helper in `tcdb_utils.py`.

**Tech Stack:** Python 3 (uv), BioCypher adapters, Neo4j post-import Cypher (`scripts/post-import.sh` + identical `.cypher`), pytest (`-m "not slow and not kg"` unit, `-m kg` live).

**Spec:** `docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md`. The explorer contract lives in `multiomics_explorer/docs/kg-specs/2026-10-23-tcdb-cleanup-explorer-asks.md`: the KG replies there are binding. Read both.

## Global Constraints

- **`level_kind` values:** `tc_class · tc_subclass · tc_family · tc_subfamily · tc_system`. `tc_specificity` must not appear anywhere after this plan.
- **`tcdb_call` values:** `eggnog · diamond_tier · diamond_best_hit`. **`tcdb_gate` values:** `direct · annotated · below_gate`. `direction`: `uptake · export · unknown`. `direction_basis`: `abc_numbering · tcdb_name`.
- **Gate rules, verbatim from spec §5.3:**
  - `direct` = `diamond_best_hit` with `rank = 1` and `tier <= 2`, OR `eggnog` onto a `tc_system` (after alias remap), OR `diamond_tier` with `tier = 1`.
  - `annotated` = any other `eggnog` edge, OR `diamond_tier` with `tier = 2`.
  - `below_gate` = everything else.
- **Best-hit tier band:** `1` if identity ≥ 70, `2` if 40 ≤ identity < 70, `3` if identity < 40.
- **Gene edge id:** `{locus_tag}-has_tcdb-{tc_id}-{tcdb_call}`. `sources` is a one-element list: `['eggnog']` or `['tcdb_diamond']`. `evidence` is `family_inferred` for eggNOG and `homology` for both diamond kinds.
- **No native bool (R5)** on any new property. `_clean_str` applies to every string property.
- **Removals are by (label, property).** MEROPS keeps its `pfam_support`, `confidence_score`, `consensus_n`, `evidence_score`, `identity`, `qcov`, `evalue`, `tier` and `best_hit_kind`.
- **Post-import parity:** `scripts/post-import.sh` and `scripts/post-import.cypher` must keep identical Cypher logic.
- **`best_hit_v1` selection is unchanged:** same systems, same order. Only the stored statistics grow.
- **Two-clone workflow:** never `docker compose` in this clone. The rebuild runs in `/home/osnat/github/multiomics_biocypher_kg` after a push (`git push git@github.com:wosnat/multiomics_biocypher_kg.git HEAD:<branch>`). Other sessions edit this clone concurrently, so always `git add` explicit paths, never `-A`.

## Review Focus

1. **The first step-6 run after the change** reads a *previous* `tcdb_pruned.json` in the old format (`subtree_substrates`, `lumping_pin`, no `system_substrates`). It must not crash or pass legacy keys through. Test in Task 4.
2. **A protein with a best hit but no tier call** (the 188 call-less proteins) must still get a `diamond_best_hit` edge. Its gate follows its tier band. Test in Task 6.
3. **eggNOG aliasing:**
   - an eggNOG 5-part id re-anchored by `seed_aliases` to an ancestor is `annotated`, not `direct`;
   - an id aliased to `""` (dropped) produces no edge;
   - two eggNOG ids aliasing to one node produce one edge.
   Test in Task 6.
4. **A strain with no `calls.json`** still yields eggNOG edges and nothing from diamond. Test in Task 6.
5. **A kept `tc_system` with no scraped name** falls back to `name = tcdb_id` and keeps `direction` from its ancestors. Test in Task 6.

---

### Task 0: Branch and baselines

**Files:** none modified. Writes local baseline files only.

- [ ] **Step 1:** Create the branch: `git switch -c feat/tcdb-cleanup` from `main`.
- [ ] **Step 2:** Against the currently deployed graph, run `uv run python tests/kg_validity/capture_annotation_state.py --save` and the `/omics-edge-snapshot` save.
- [ ] **Step 3:** Copy `cache/data/kegg/kegg_data.json` and `cache/data/tcdb/tcdb_pruned.json` to `/tmp/claude-1001/-home-osnat-github-second-multiomics/d62c75b7-5c6e-451a-bc65-a5e09552f311/scratchpad/baseline/`. These are used for the before/after diff in Task 10.

### Task 1: Per-system best-hit statistics in the diamond artifact

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
  - `select_best_hits` (~l.194)
  - best-hit block in `build_strain_calls` (~l.426–440)
  - `_TIER_TO_LEVEL_KIND` (l.260)
  - add `best_hit_tier(identity: float) -> int`
- Test: `tests/test_tcdb_diamond.py`
- Regenerate: every `cache/data/*/genomes/*/tcdb/*.tcdb.calls.json` and `*.tcdb.skill_summary.json`

**Interfaces:**
- **Produces:** the calls.json `best_hit` block is
  `{"rule", "tcids", "identity", "top_tcid"?, "attach_tcid"?, "hits": [{"tcid", "identity", "qcov", "scov", "evalue", "bitscore"}]}`.
  - `hits` is aligned index-for-index with `tcids`, so `hits[0]` is rank 1.
  - Each hit's statistics come from the first row `select_best_hits` keeps for that TCID (its highest-bitscore subject row).
- **Produces:** `best_hit_tier(identity) -> int` (1/2/3, using the bands in Global Constraints).

- [ ] **Step 1: Write the failing tests.**
  - `test_best_hit_hits_aligned_with_tcids`: `best["hits"][i]["tcid"] == best["tcids"][i]` for every `i`; `len(hits) == len(tcids)`.
  - `test_best_hit_hits_use_each_systems_own_row`: for a fixture with two TCIDs whose top rows differ, assert each hit's `identity`/`qcov`/`scov`/`evalue`/`bitscore` equals that TCID's highest-bitscore row.
  - `test_best_hit_selection_unchanged`: the `tcids` order and membership are identical to the current output on the existing window/tie fixtures (`test_best_hit_window_cross_family_top3`, `test_best_hit_ties_deterministic`).
  - `test_best_hit_tier_bands`: `best_hit_tier(70.0) == 1`, `best_hit_tier(69.99) == 2`, `best_hit_tier(40.0) == 2`, `best_hit_tier(39.99) == 3`.
  - `test_tier_level_kind_uses_tc_system`: `_TIER_TO_LEVEL_KIND[1] == "tc_system"`. Update `test_build_strain_calls_full_pipeline`'s `level_kind` assertion to match.
- [ ] **Step 2:** Run `uv run pytest tests/test_tcdb_diamond.py -v`. Expected: the new tests FAIL.
- [ ] **Step 3: Implement.**
  - Keep the per-TCID first row in `select_best_hits`, and return it as `hits`.
  - Copy `hits` into the block in `build_strain_calls`.
  - Rename `tc_specificity` → `tc_system`.
  - Add `best_hit_tier`.
- [ ] **Step 4:** Run `uv run pytest tests/test_tcdb_diamond.py -v`. Expected: PASS.
- [ ] **Step 5:** Regenerate the artifacts: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --from-tsv`. Expected: `OK` for all 43 strains.
- [ ] **Step 6:** Check that the selection did not move. Compare the old and new calls.json for all strains: per protein, `best_hit.tcids` must be identical. Report any difference and stop if there is one.
- [ ] **Step 7: Commit.**
  ```bash
  git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
  git add 'cache/data/*/genomes/*/tcdb/*.tcdb.calls.json'
  git add 'cache/data/*/genomes/*/tcdb/*.tcdb.skill_summary.json'
  git commit
  ```
  Commit message: `feat(tcdb-diamond): per-system best-hit stats; tc_specificity -> tc_system`.

### Task 2: Per-system substrates in the metabolite-identity pipeline

**Files:**
- Modify: `multiomics_kg/utils/metabolite_identity/pipeline.py`
  - `Step6Inputs` (l.51)
  - `_best_hit_substrates` (l.342)
  - `run_identity_stages`: stage F (l.392), stage O (l.466), stage P (l.472), writer (l.536–539)
  - `lumping_report` (l.590/615)
- Modify: `canonical.py` (`rewrite_to_canonical`, l.367): rename the parameter to `system_substrates`; no logic change.
- Modify: `outputs.py` (`validate_build`, l.20): keyword `system_substrates` / `system_roles`; same dangling checks.
- Modify: `roles.py`: rename `build_subtree_roles` → `build_system_roles(system_strings, role_map) -> {tc: sorted role ids}`; same logic, non-empty only.
- Modify: `__init__.py` exports.
- Delete: `lumping.py`, `tests/test_mi_lumping.py`.
- Test: `tests/step6_fixture.py`, `tests/test_step6_identity_pipeline.py`, `tests/test_mi_roles.py`, `tests/test_mi_outputs.py`

**Interfaces:**
- **Produces:** `Step6Inputs` has `system_substrates: dict[str, list[str]]` (`{tc_system: [raw-resolved primary ids]}`) and `system_roles: dict[str, list[str]]`. Fields `subtree_substrates`, `subtree_roles`, `lumping_pin` and `best_hit_substrates` are removed.
- **Produces:** the `tcdb_pruned` written by stage P has exactly the keys `kept_tcdb_ids`, `system_substrates`, `system_roles`, `chemical_roles`, `seed_aliases`, `pfam_bridge`, `go_bridge`. The values come from `tcdb_static` plus the canonicalised maps.

- [ ] **Step 1: Write the failing tests.**
  - `tests/step6_fixture.py`: build inputs with `system_substrates` / `system_roles` keyed by 5-part ids only.
  - `test_pruned_keys_after_cleanup`: the written `tcdb_pruned` keys equal the set above.
  - `test_system_substrates_canonicalised`: a system listing a non-canonical id ends up with its canonical id after stage F.
  - `test_validate_flags_dangling_system_substrate` (renamed from `test_validate_flags_dangling_and_missing_assignment`).
  - `test_system_roles_sorted_nonempty_only`.
  - Update `test_rerun_with_committed_lock_and_pin_is_byte_identical` → `test_rerun_with_committed_lock_is_byte_identical`: drop the pin.
  - Delete the best-hit-key tests at l.366 and l.386.
- [ ] **Step 2:** Run `uv run pytest tests/test_step6_identity_pipeline.py tests/test_mi_roles.py tests/test_mi_outputs.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement** the renames. Remove stage O and the lumping report: the identity report loses its `lumping` key. Delete `lumping.py`.
- [ ] **Step 4:** Run the same command plus `uv run pytest tests/ -k "mi_ or step6" -q`. Expected: PASS. `test_build_kegg_metabolism_xrefs.py` still fails until Task 3.
- [ ] **Step 5: Commit.** Message: `refactor(step6): system-level substrates in identity pipeline; retire lumping pin`.

### Task 3: Step 6 — seeds, pruning, system substrates, bridges

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
  - `build_tcdb_hierarchy` (l.289; `tc_specificity` literal at l.376 → `tc_system`)
  - `_gene_reachable_sets` (l.447)
  - `_compute_subtree_substrate_strings` (l.798) → replace
  - `_best_hit_substrate_inputs` (l.1014) → delete
  - `_warn_unnamed_best_hits` (l.192) → generalise
  - `_tcdb_bridges` (l.1416)
  - `main` (l.1624–1872)
- Modify: `multiomics_kg/utils/tcdb_utils.py` (docstring l.142)
- Test: `tests/test_build_kegg_metabolism_xrefs.py`, `tests/test_tcdb_utils.py`

**Interfaces:**
- **Consumes:** `Step6Inputs(system_substrates=…, system_roles=…)` (Task 2).
- **Produces:**
  - `_gene_reachable_sets(raw)["tcdb_seed_ids"]` = ∪ of `transporter_classification` and `tcdb_best_hit_ids`.
  - `_add_range_systems(hierarchy: dict, kept: set[str]) -> tuple[set[str], set[str]]` returns `(kept_with_range, range_only)`. It adds every `tc_system` with non-empty `substrate_classes` whose 3-part family is in `kept`, plus its ancestors.
  - `_system_substrate_strings(hierarchy, kept) -> dict[str, set[str]]` = `{tc_system in kept: set(substrate_classes)}`, non-empty only.
  - `_warn_unnamed_systems(kept_ids, hierarchy) -> int`: logs and returns the number of kept `tc_system` nodes whose `name` is missing or equals the id.
- **The KEGG extension** uses `substrate_kegg_cpds` from `system_substrates` primaries. `transport_pws` is built from those, with the same rule as today.

- [ ] **Step 1: Write the failing tests.**
  - `test_gene_reachable_sets_seeds_include_best_hits`.
  - `test_prune_with_best_hit_seeds_keeps_their_ancestors`.
  - `test_add_range_systems_adds_substrate_systems_under_kept_families_only`:
    - a system with substrates in a kept family is added;
    - one in an unkept family is not;
    - one with no substrates is not;
    - its ancestors are added.
  - `test_system_substrate_strings_never_on_non_systems`: all keys have 5 parts.
  - `test_hierarchy_uses_tc_system` (update `test_step6_builds_tcdb_hierarchy` and `test_build_tcdb_hierarchy_seeds_5part_from_acc2tcid`).
  - `test_bridge_fallback_warns_when_previous_lacks_range_nodes`: with no map and a previous bridge missing some kept nodes, a warning is logged that names the count.
  - `test_main_tolerates_previous_pruned_in_old_format` (Review Focus 1): a previous `tcdb_pruned.json` with only `subtree_substrates` / `lumping_pin` / `best_hit_substrates` runs through `_stub_main` without error. Output has no legacy keys.
  - Update `test_main_passes_old_and_new_substrate_maps` so it still asserts exactly one MNX resolve pass.
  - Delete: `test_compute_subtree_substrate_strings_rolls_up_through_full_hierarchy`, `test_best_hit_substrate_inputs`, `test_warn_unnamed_best_hits_separates_unscraped_from_missing`.
- [ ] **Step 2:** Run `uv run pytest tests/test_build_kegg_metabolism_xrefs.py tests/test_tcdb_utils.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - Pruning: `_prune_tcdb(hierarchy, seeds)`, then `_add_range_systems`. Log the gene-anchored and range-only counts.
  - Substrates and roles: `_system_substrate_strings`, then role partition (`partition_role_strings`, `build_system_roles`), then **one** `_resolve_substrate_strings` call over every system string. Pass the result to `Step6Inputs`.
  - Bridges: `_tcdb_bridges` uses the full kept set. On the no-map fallback, warn when the previous bridge covers fewer kept nodes than now.
  - Remove the `lumping_pin` / `best_hit_substrates` plumbing.
- [ ] **Step 4:** Run the same pytest command. Expected: PASS. Then run `uv run pytest -m "not slow and not kg" -q -k "kegg or tcdb or mi_ or step6"`. Expected: PASS, apart from adapter tests (Task 6).
- [ ] **Step 5: Commit.** Message: `feat(step6): seed pruning from best hits, keep substrate-bearing systems, system-level substrates`.

### Task 4: Names scraper scope

**Files:**
- Modify: `multiomics_kg/download/scrape_tcdb_names.py`
  - `scoped_families` (l.209)
  - `main` (l.319)
- Test: `tests/test_scrape_tcdb_names.py`

**Interfaces:**
- **Produces:** `scoped_families(kept_ids)` also returns families kept only at 3-part level. Their query id is the family id itself.
- **Produces:** `main` scope = `kept_tcdb_ids` from `tcdb_pruned.json`, which now holds every kept system. The best-hit ids stay unioned.

- [ ] **Step 1: Write the failing tests.**
  - `test_scoped_families_includes_family_only_kept_ids`: `scoped_families({"2.A.1", "1.A.11.2.5"})` contains `"2.A.1"`.
  - `test_assemble_names_names_range_only_systems`: a system present on a page and in `kept_ids` gets a name.
- [ ] **Step 2:** Run `uv run pytest tests/test_scrape_tcdb_names.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4:** Run the same command. Expected: PASS.
- [ ] **Step 5: Commit.** Message: `feat(tcdb-names): scope covers family-only kept ids and range-only systems`.

### Task 5: Run step 6, the scraper, then step 6 again; commit the caches

**Files (regenerated and committed):**
- `cache/data/tcdb/{tcdb_hierarchy,tcdb_pruned,tcdb_names}.json`
- `cache/data/kegg/kegg_data.json`
- `cache/data/metabolomics/{metabolite_id_mapping,metabolite_identity_report}.json`
- regenerated paper `_resolved.csv` files, if any changed

Steps:
- [ ] **Step 1:** Run `bash scripts/prepare_data.sh --steps 6 --force`. Expected: completes. The log shows gene-anchored ≈ 3.5K and range-only ≈ 5.2K kept nodes. Requires `MNX_DATA_DIR`.
- [ ] **Step 2:** Run `uv run python -m multiomics_kg.download.scrape_tcdb_names`. Expected: about 112 new page fetches (~5 min) and names for the new systems. Network is required. If TCDB is down, stop and report.
- [ ] **Step 3:** Run `bash scripts/prepare_data.sh --steps 6 --force` again. Expected: the unnamed-systems count is logged and small.
- [ ] **Step 4:** Check against the Task 0 baseline and `docs/kg-changes/tcdb-cleanup/removed_ids.tsv` (the predicted id list the explorer is sweeping on). The `kegg_data.json` difference must equal spec §5.4 and that file:
  - −65 transport-only compounds;
  - 32 compounds that lose `'transport'` only;
  - −3 pathways: Phagosome, Focal adhesion, chondroitin sulfate GAG.

  `tcdb_pruned.json` keys must equal the Task 2 set. **Stop and report** if the difference is anything else.
- [ ] **Step 5:** Commit the regenerated caches, listing explicit paths. Message: `data(step6): system-level TCDB substrates, best-hit seeding, range-only systems`.

### Task 6: TCDB adapter rewrite

**Files:**
- Modify: `multiomics_kg/adapters/tcdb_adapter.py` (whole file)
- Modify: `multiomics_kg/utils/tcdb_utils.py`: add `tcdb_direction(tc_id: str, hierarchy: dict) -> tuple[str, str | None]`
- Modify: `create_knowledge_graph.py` (l.243 and l.332–342 keep their injected arguments; only the keys read change)
- Test: `tests/test_tcdb_adapter.py`, `tests/test_tcdb_utils.py`

**Interfaces:**
- **Consumes:**
  - calls.json `best_hit.hits` (Task 1);
  - `tcdb_pruned.json` keys (Task 2);
  - hierarchy `level_kind == "tc_system"` (Task 3);
  - `best_hit_tier` (Task 1).
- **Produces:** `tcdb_direction`, which implements spec §5.1:
  - ABC rule: `3.A.1.<n>` and below → `uptake` if `n < 100`, else `export`. Basis `abc_numbering`.
  - Otherwise: walk from the node up to family level (never class or subclass). The first name with a whole-word, case-insensitive match of exactly one of {uptake, import, influx} vs {export, efflux, extrusion} gives the direction. Basis `tcdb_name`.
  - Otherwise `("unknown", None)`.
- **Produces node properties:** `name`, `tcdb_id`, `level`, `level_kind`, `superfamily`?, `description`?, `direction`, `direction_basis`?. There is no `substrate_breadth`.
- **Produces gene edge properties per `tcdb_call`:** as in spec §5.3 / Global Constraints, plus `tcdb_gate`. `rank` = index in `best_hit.tcids` + 1.
- **Produces substrate and role edges:** only from `system_substrates` / `system_roles`; ids unchanged (`{tc}-transports-{primary}`, `t2r:{tc}:{role}`); no properties.

- [ ] **Step 1: Write the failing tests.** Replace the merge, depth, lumping and best-hit-on-diamond-edge tests listed in the code map with these.
  - `test_one_edge_per_tcdb_call`: a gene with eggNOG + tier-1 + best hit on one system → 3 edges with distinct ids ending `-eggnog` / `-diamond_tier` / `-diamond_best_hit`.
  - `test_edge_properties_by_call_kind`:
    - `eggnog`: no `tier` / `identity`;
    - `diamond_tier`: `tier`, `identity`, `qcov`, `scov`, `evalue`;
    - `diamond_best_hit`: `rank`, `tier`, `identity`, `qcov`, `scov`, `evalue`, `bitscore`, `bitscore_ratio`.
    - No `confidence_score`, `consensus_n` or `best_hit_*` anywhere.
  - `test_tcdb_gate_rules`, parametrised over:
    - best hit rank 1 at identity 85 → `direct`;
    - rank 1 at 39 → `below_gate`;
    - rank 2 at 85 → `below_gate`;
    - eggNOG 5-part → `direct`;
    - eggNOG 4-part → `annotated`;
    - tier 1 → `direct`;
    - tier 2 → `annotated`;
    - tier 3 → `below_gate`.
  - `test_eggnog_alias_to_ancestor_is_annotated` (Review Focus 3).
  - `test_eggnog_alias_dropped_emits_nothing` (Review Focus 3).
  - `test_two_eggnog_ids_aliasing_to_one_node_give_one_edge` (Review Focus 3).
  - `test_within_kind_duplicates_keep_higher_identity`.
  - `test_best_hit_without_call_emits_best_hit_edge` (Review Focus 2).
  - `test_missing_calls_json_yields_eggnog_only` (Review Focus 4).
  - `test_bitscore_ratio_is_one_at_rank_one`.
  - `test_substrate_edges_only_from_systems`: every substrate and role edge source has 5 parts.
  - `test_node_props`: `level_kind` ∈ the five values; no `substrate_breadth`; `direction` on every node.
  - `test_unnamed_system_falls_back_and_inherits_direction` (Review Focus 5).
  - `test_tcdb_direction`, in `test_tcdb_utils.py`:
    - `3.A.1.4.4` → `("uptake", "abc_numbering")`;
    - `3.A.1.122.2` → `("export", "abc_numbering")`;
    - a family named "… Efflux …" → `("export", "tcdb_name")`, and its child system inherits it;
    - a name with both words → `unknown`;
    - class `3` → `("unknown", None)`.
  - Keep the bridge tests, extended with a range-only node receiving a bridge edge.
  - Delete:
    - `test_missing_lumping_pin_fails`
    - `test_substrate_breadth`
    - `test_transport_substrate_resolution_reads_the_pinned_set`
    - `test_is_promiscuous_is_gone_from_post_import`, unless it still holds
- [ ] **Step 2:** Run `uv run pytest tests/test_tcdb_adapter.py tests/test_tcdb_utils.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - Per-strain `TcdbAnnotationAdapter.get_edges` yields raw `(tcdb_call, tc, props)` claims:
    - from `tcdb_eggnog_ids`;
    - from calls.json `calls` (tier);
    - from `best_hit.hits` — **all** proteins with a best hit, not only those with `attach_tcid`.
  - The Multi wrapper:
    - applies `seed_aliases`;
    - drops unkept targets;
    - collapses within kind (keeping higher `identity`);
    - computes `tcdb_gate` **after** the alias remap, using the target's level from the hierarchy;
    - logs collapsed and remapped counts.
  - Delete `_compute_substrate_depth`, `_with_best_hit_details` and the lumping-pin read.
- [ ] **Step 4:** Run the same pytest command. Expected: PASS. Then run `uv run pytest -m "not slow and not kg" -q`. Expected: PASS, apart from vocab and schema tests (Task 7).
- [ ] **Step 5: Commit.** Message: `feat(tcdb-adapter): one edge per tcdb_call with tcdb_gate; system-only substrates; direction`.

### Task 7: Schema, controlled vocabulary, vocabulary contract

**Files:**
- Modify: `config/schema_config.yaml`
  - Gene l.581–610
  - OrganismTaxon l.694
  - Metabolite l.906–914
  - ChemicalRole l.1318–1353
  - TcdbFamily l.1387–1407
  - gene edge l.1425–1474
  - substrate edge l.1484–1499
- Modify: `config/controlled_vocabularies.yaml`. Entries by line:
  - remove: 218, 236, 253, 268, 282, 296, 435, 449, 564, 2047, 2147, 2204, 2216–2295;
  - edit: 153, 168, 876, 894, 1001, 1286, 1985.
- Modify: `docs/kg-changes/vocabulary-contract.md` (R3 exceptions l.348, l.355–357)
- Test: `tests/test_controlled_vocab.py` (l.100–112, l.179–189), `tests/test_controlled_vocabulary_adapter.py`

**Interfaces:**
- **Produces** schema node `tcdb family`: `name, tcdb_id, level, level_kind, superfamily, description, direction, direction_basis, tc_class_id, gene_count, direct_gene_count, organism_count, member_count, tcdb_substrate_count`.
- **Produces** gene edge: `sources, evidence, tcdb_call, tcdb_gate, tier, rank, identity, qcov, scov, evalue, bitscore, bitscore_ratio`.
- **Produces** vocab entries:
  - Trust vocabularies with `value_descriptions` on every value (descriptions state the exact rule):
    - `Gene_has_tcdb_family.tcdb_call`
    - `Gene_has_tcdb_family.tcdb_gate`
    - `TcdbFamily.direction`
    - `TcdbFamily.direction_basis` (sparse)
  - Numeric entries:
    - `Gene_has_tcdb_family.rank` (min 1, max 3)
    - `.bitscore`
    - `.bitscore_ratio` (0.9–1.0)
    - `.scov`
    - `TcdbFamily.tcdb_substrate_count` (min 0)
  - Edited entries:
    - `Gene_has_tcdb_family.tier`: two derivations, not sparse on diamond kinds, absent on eggNOG, never compare across kinds.
    - `TcdbFamily.level_kind`: `tc_system`.
    - `Metabolite.evidence_sources`, `transport` = "TCDB lists it as a substrate of a system in the KG", not gene reach.
    - `Organism_has_metabolite.evidence_sources`, transport arm = `direct` edges.
    - `Gene.annotation_state`: describe the transporter bucket gate.
  - Descriptions carry no build provenance (`test_descriptions_carry_no_build_provenance`).

- [ ] **Step 1: Write the failing tests.**
  - In `test_evidence_ladder_placements`, replace the `attachment_depth` assertion with `tcdb_call == ["eggnog", "diamond_tier", "diamond_best_hit"]` and `tcdb_gate == ["direct", "annotated", "below_gate"]`.
  - In `test_trust_slice_carries_value_descriptions`, `trust_props` drops `attachment_depth`, `substrate_depth`, `pfam_support`, `go_support`, `source_agreement` and `substrate_breadth`. It adds `tcdb_call`, `tcdb_gate`, `direction` and `direction_basis`. MEROPS `pfam_support` stays covered.
  - New `test_tcdb_removed_entries_absent`: none of the removed `(label.property)` keys exist; every MEROPS key listed in Global Constraints still exists.
- [ ] **Step 2:** Run `uv run pytest tests/test_controlled_vocab.py tests/test_controlled_vocabulary_adapter.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - Edit both YAMLs.
  - In `vocabulary-contract.md`, mark the `substrate_breadth` R3 exception and the `attachment_depth` R3-exempt pair as retired. Add the `tcdb_gate` R3 exception with its reason.
- [ ] **Step 4:** Run the same command, then `uv run pytest -m "not slow and not kg" -q`. Expected: PASS, apart from the post-import static tests (Task 8).
- [ ] **Step 5: Commit.** Message: `feat(schema,vocab): TCDB cleanup contract`.

### Task 8: Post-import

**Files:** `scripts/post-import.sh` and `scripts/post-import.cypher`, kept identical.
- Test: `tests/test_annotation_quality_buckets.py`, `tests/test_post_import_metabolite_identity.py`, the post-import text assertions at `tests/test_tcdb_adapter.py:620-631`.

**Interfaces:**
- **Consumes:** the edge properties `tcdb_gate`, `tcdb_call` (Tasks 6–7).
- **Produces:** post-import properties defined in spec §5.1 and §7. Each block is rewritten as follows (`.sh` lines):
  - **728, 778, 838 (gates):** replace the gate predicate with `rt.tcdb_gate IN ['direct','annotated']`. Leave `has_any_edge` (733–735) as is.
  - **1042–1058 (TcdbFamily counts):**
    - `gene_count` / `direct_gene_count` / `organism_count` count only edges with `tcdb_gate IN ['direct','annotated']`.
    - Delete `metabolite_count`.
    - Add `tcdb_substrate_count` = distinct targets of `Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role` from `tc_system` nodes in the subtree (`*0..`), on all levels, with 0 when none.
  - **1481–1503:** delete (`attachment_depth`).
  - **1505–1552 (Gene transport arm):**
    - `tcdb_family_count` uses the spec §5.3 rule: gate-passing targets minus any node that is an ancestor of another.
    - `transported_metabolite_count` = distinct metabolites of the targets of `tcdb_gate = 'direct'` edges.
    - Delete `transport_substrate_resolution`.
  - **1554–1626:** delete the TCDB `source_agreement` / `pfam_support` / `go_support` / `evidence_score` / `tcdb_evidence_score_max` blocks. Leave the MEROPS blocks at 1163–1203 untouched.
  - **1677–1719 (Metabolite):**
    - `transporter_gene_count` = genes with a `direct` edge to a system carrying it.
    - `transporter_count` = distinct such systems.
  - **1741–1777:** the `Organism_has_metabolite` transport arm and `OrganismTaxon.transported_metabolite_count` use `direct` edges.
  - **1826–1839 and 1973–1982:** fix comments only.
  - **2051–2067 (ChemicalRole):**
    - `transporter_count` / `direct_transporter_count` = gene-direct systems.
    - **`tcdb_system_count` / `direct_tcdb_system_count`** (renamed from `tcdb_family_count` / `direct_tcdb_family_count`, explorer C4; also in schema, Task 7, and CHANGELOG `### Breaking`) = all kept systems carrying the role.
    - Keep the literal `r.direct_transporter_count`.
  - **1060–1078:** remove the stale `is_promiscuous` history comment.
  - **688–713 header comments:** update the gate description.

- [ ] **Step 1: Write the failing tests.**
  - Update the static text test at `tests/test_tcdb_adapter.py:620-631`: assert `tcdb_gate IN ['direct','annotated']` is present in both files, and that `substrate_breadth`, `attachment_depth`, `substrate_depth` and `tcdb_evidence_score_max` are absent.
  - Add `test_post_import_sh_and_cypher_tcdb_blocks_match`: extract the TCDB block text from both files and compare.
- [ ] **Step 2:** Run `uv run pytest tests/test_tcdb_adapter.py -k post_import tests/test_annotation_quality_buckets.py tests/test_post_import_metabolite_identity.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement** in both files.
- [ ] **Step 4:** Run the same command, then `uv run pytest -m "not slow and not kg" -q`. Expected: all PASS.
- [ ] **Step 5: Commit.** Message: `feat(post-import): TCDB gate, system-level transport arm, tcdb_substrate_count`.

### Task 9: KG validity tests and snapshot generator

**Files:**
- `tests/kg_validity/test_tcdb_cazy.py`
- `test_tcdb_best_hit.py`
- `test_annotation_trust.py`
- `test_chemistry_slice1_rollups.py`
- `test_metabolism.py`
- `test_metabolite_identity.py` + `fixtures/metabolite_identity_pre_counts.json`
- `test_controlled_vocabularies.py`
- `test_annotation_state.py`
- `test_post_import.py`
- `test_structure.py`
- `test_snapshot.py`
- `generate_snapshot.py`

**Tests to write.** Each spec §9 acceptance bullet becomes one kg test, in `test_tcdb_cazy.py` unless noted:
- `test_substrate_edges_leave_systems_only`: 0 substrate or role edges from `level_kind <> 'tc_system'`.
- `test_gene_edge_shape_by_tcdb_call`:
  - one `tcdb_call` and a one-element `sources` per edge;
  - `tier` on both diamond kinds only;
  - `rank` / `bitscore` / `bitscore_ratio` on best-hit edges only;
  - exactly one `rank = 1` per gene that has best hits;
  - `tcdb_gate` matches the rule recomputed in Cypher.
- `test_transport_arm_nests_in_annotation`: 0 genes with `transported_metabolite_count > 0` and without `'tcdb'` in `annotation_types`.
- `test_tcdb_substrate_count_equals_walk`: on a sample of nodes at every level, the property equals a `*0..` walk count.
- `test_tcdb_gene_count_uses_gate`.
- `test_tcdb_family_count_rule` (spec §5.3 examples, recomputed in Cypher).
- `test_med4_spot_checks`:
  - `PMM0970` has a `diamond_best_hit` edge with `rank = 1` to `tcdb:3.A.1.4.4`, which transports urea;
  - `PMM0710` best hit is `3.A.1.7.6`;
  - `tcdb:3.A.1` has 0 substrate edges and `direction = 'unknown'`;
  - `tcdb:3.A.1.4` is `uptake` / `abc_numbering`.
- `test_systems_are_named`: kept systems with `name = tcdb_id` ≤ the count step 6 logged.

**Rewrite or delete** each test the code map lists for removed properties. That covers the `test_tcdb_cazy.py` tests at 69, 103–130, 159, 176, 260–419, 451–587, plus the listed tests in the other files.

**Fixture updates:**
- `fixtures/metabolite_identity_pre_counts.json`: refresh the expectations to the measured deltas (−65 Metabolite nodes, −3 pathways).
- `test_metabolite_pathways_only_kg_evidenced`: tighten back. Transport evidence = system substrate edges.
- `test_snapshot_edge_exists`: when an edge spec carries `tcdb_call`, match on it, the way the existing time-point filter works.
- `generate_snapshot.py`:
  - l.252–259: drop `confidence_score` / `consensus_n`; add `tcdb_call` / `tcdb_gate`.
  - l.173: add `direction` to node properties.

Steps:
- [ ] **Step 1:** Write and rewrite the tests above.
- [ ] **Step 2:** Run `uv run pytest tests/kg_validity --collect-only -q`. Expected: collects with no errors. They run against the live graph in Task 10.
- [ ] **Step 3: Commit.** Message: `test(kg): TCDB cleanup acceptance tests`.

### Task 10: Rebuild and verify

- [ ] **Step 1:** Push: `git push git@github.com:wosnat/multiomics_biocypher_kg.git HEAD:feat/tcdb-cleanup`. Ask the user to pull and run `docker compose up -d --build` in `/home/osnat/github/multiomics_biocypher_kg`. **Wait** for confirmation.
- [ ] **Step 2:** Check `./output/import.status` and `import.report` in the other clone. Expected: 0 skipped relationships, and no duplicate-id drops for `Gene_has_tcdb_family`.
- [ ] **Step 3:** Run `uv run pytest -m kg -v`. Expected: all PASS. Fix and repeat on failure.
- [ ] **Step 4:** Run `uv run python tests/kg_validity/capture_annotation_state.py --compare`. Expected: movement only from the ~1,055 genes that gain `'tcdb'`. Anything else is a bug: stop and investigate.
- [ ] **Step 5:** Run the `/omics-edge-snapshot` compare. Expected: `Changes_expression_of` counts unchanged.
- [ ] **Step 6:** Run `uv run python tests/kg_validity/generate_snapshot.py`, then `uv run pytest tests/kg_validity/test_snapshot.py -v`. Expected: PASS. Commit `snapshot_data.json`.

### Task 11: Documentation

**Files:**
- `CLAUDE.md`: lines 82, 96, 400–419, 454, 519–523, 537–564, 568, 590–592, 599, per the code map.
- `CHANGELOG.md` `## [Unreleased]`:
  - `### Breaking`: removed properties; `tc_specificity` → `tc_system`; one edge per `tcdb_call`; transport arm ~30K → ~8.3K genes; −65 Metabolite nodes and −3 pathways; gated `TcdbFamily.gene_count`.
  - `### Added`: `tcdb_call`, `tcdb_gate`, `direction`, `tcdb_substrate_count`, best-hit edges with statistics, range-only systems.
  - `### Highlights`: one sentence.
- New `docs/kg-changes/tcdb-cleanup.md`: the what-changed doc for consumers, with measured live numbers from Task 10. Confirm or correct `docs/kg-changes/tcdb-cleanup/removed_ids.tsv` against the rebuilt graph, and drop its "Predicted" header line once confirmed.
- `.claude/skills/tcdb-diamond/SKILL.md`: the per-system `hits` block, `tc_system`, and the `--from-tsv` note.
- `plans/backlog.md`: tick 3a; close Project 2 follow-up "188 proteins" and Deferred-TCDB item (f).
- `multiomics_explorer/docs/kg-specs/2026-10-23-tcdb-cleanup-explorer-asks.md`: append "KG built" with live counts and the branch/commit.

Steps:
- [ ] **Step 1:** Write the docs, using measured numbers from Task 10.
- [ ] **Step 2:** Run `uv run pytest tests/test_controlled_vocab.py -q`. Expected: PASS (description hygiene).
- [ ] **Step 3: Commit** with explicit paths. Message: `docs: TCDB cleanup (project 3a)`.
