# GO Transport Substrates → Metabolite Edges Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish GO's transport-term substrate axioms (go-plus *has primary input*) as edges from GO term nodes to canonical `Metabolite` / `ChemicalRole` nodes, filtered of ion/class noise, with a separate GO-route organism rollup and a published exclusion list.

**Architecture:** prepare_data step 6 (which owns ChEBI/MNX canonicalisation) downloads go-plus, extracts transport axioms, resolves their ChEBI ids exactly like TCDB substrate strings, threads them through the identity pipeline (canonical rewrite, validation), applies the target filter before and after canonicalisation, and writes a `go_transport_substrates` block into `kegg_data.json`. `MultiGoAnnotationAdapter` emits the four edge types from that block, dangling-proofed against emitted GO / Metabolite / ChemicalRole ids. Post-import adds the organism rollup and stamps the exclusion list on `Schema_info`.

**Tech Stack:** Python 3 (uv), BioCypher, Neo4j 5 + APOC (post-import Cypher), pytest.

**Spec:** `docs/superpowers/specs/2026-10-18-go-transport-substrates-design.md` (read it with this plan). Explorer dialog: `multiomics_explorer/docs/kg-specs/2026-10-18-go-transport-substrates-explorer-asks.md`.

## Global Constraints

- Relation: only RO:0004009 *has primary input*; RO:0012001 *has small molecule activator* and PR / GO fillers are skipped.
- Transport term = logical-definition genus is, or `is_a`-descends from, GO:0005215 or GO:0006810.
- Role fillers: `GOCHE_<n>` → `CHEBI:<n>` (the GOCHE number is the ChEBI role id); role edges only to ChemicalRole ids already in `tcdb_pruned.json` `chemical_roles` — no new role nodes.
- Each GO term carries only its own axiom; no propagation (anti-goal A1).
- New `Metabolite.evidence_sources` value is exactly `go_transport`; `transport` keeps meaning TCDB (A3).
- No pathway pull-in: GO-only KEGG compounds join the compound set, never the pathway set.
- Target filter (canonical-id match): CHEBI:24870 ion, 24867 monoatomic ion, 23906 monoatomic cation, 23905 monoatomic anion, 25213 metal cation, 33515 transition element cation, 51143 nitrogen molecular entity, 33694 biomacromolecule, 24431 chemical entity, 36080 protein, 24834 inorganic anion, 15378 hydron (H⁺), 29101 sodium(1+). K⁺, Cl⁻, metals, iron, nutrient ions are kept.
- Edge labels (schema `label_as_edge`): `Molecular_function_transports_metabolite`, `Biological_process_transports_metabolite`, `Molecular_function_transports_chemical_role`, `Biological_process_transports_chemical_role`; one property `go_chebi_ids: str[]` (sorted); one edge per (GO term, canonical target).
- New post-import properties: `Metabolite.go_transport_organism_names` (str[], sorted `preferred_name`), `Metabolite.go_transport_organism_count` (int, `= size(names)`), sparse — present only on metabolites with ≥ 1 GO transport edge; `Schema_info.go_transport_excluded_targets` (str[]).
- All string properties through `_clean_str` (no `'` or `|`, never used as separators).
- `scripts/post-import.sh` and `scripts/post-import.cypher` carry identical Cypher.
- Staging: other sessions edit this clone — `git add` explicit paths only, never `-A`.
- Docker build runs in the other clone (`multiomics_biocypher_kg`): commit + push here, the user pulls and rebuilds there.
- Release date = go-plus `owl#versionInfo` (2026-07-26 for the measured file), never the download date.
- **One GO release, refreshed as a unit (Task 0).** `go-basic.obo` (the GO node set, `go_namespace_cache.json`) and `go-plus.json` (transport axioms) are downloaded **together** from the current URLs and accepted only if go-basic's `data-version: releases/R` equals go-plus's `owl#versionInfo` R (dated release URLs 404 on both purl.obolibrary.org and release.geneontology.org, checked 2026-10-18, so a release cannot be fetched by date). Both raw files stay in gitignored `cache/data/go_terms/raw/`; `go_release.txt` records R. Everything downstream reads those files and never downloads GO itself; a refresh happens only through `refresh_go_release(force=True)` (step 6 `--refetch-raw`). The KG build never refetches GO: `MultiGoAnnotationAdapter` no longer passes `force=not cache` to `load_go_data`, and it warns when `go_release.txt` ≠ `kegg_data["go_transport_release"]`. Baseline release: **2026-07-26** (current on 2026-10-18; the committed namespace cache was built from a March 2026 go-basic, 38,739 terms).
- **Order inside step 6** (each needs the previous): (1) gene annotations from step 2 (GO term seeds) → (2) GO files from Task 0's refresh (raw go-plus → `go_transport_substrates.json`, release checked against `go_release.txt`) → (3) ChEBI load (stage A; GO strings need ChEBI names and roles) → (4) TCDB role partition (stage C; GO role edges target its `chemical_roles`) → (5) resolve GO + filter strings with the TCDB molecule strings' resolver (stage D) → (6) extend the compound set (not pathways) **before** `build_pruned_kegg_data` → (7) `_fold_substrates_into_kegg_data` → (8) identity stages E–P rewrite / filter / validate.
- Step 6 skips when `kegg_data.json` exists without `--force` (existing rule): a new strain's GO terms reach the map only with `--steps 6 --force`, the same as its TCDB substrates.

## Review Focus

1. **A KG GO term absent from go-plus** (obsolete id, version skew between `go_namespace_cache.json` and go-plus) is skipped and counted, never raises — test in Task 4.
2. **A kept GO target that canonicalises onto a filtered node** (merge in stage E/F) is dropped by the post-canonical filter — test in Task 4.
3. **go-plus download failure** with the committed derived cache present reuses it; with it absent, step 6 fails loudly — test in Task 2.
4. **`--test` build** (metabolism adapter caps emitted metabolites at 100): GO edges to unemitted metabolites are pruned, no dangling — test in Task 5.
5. **A role filler whose ChemicalRole TCDB did not create** is dropped and counted, not emitted — test in Task 4.

---

### Task 0: Refresh GO to the current release and re-baseline the KG

The GO node-set refresh (March 2026 → 2026-07-26) moves GO node / gene→GO edge counts on its own. It ships and is rebuilt **first**, so every later gate compares against a graph that already carries the new GO release and sees additions only.

**Files:**
- Create: `multiomics_kg/download/refresh_go_release.py`
- Modify: `multiomics_kg/utils/go_utils.py` (`parse_obo` gains the header's `data-version`; add `GO_RELEASE_FILE = "go_release.txt"`, `read_go_release(cache_root) -> str | None`)
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:208` (`load_go_data(Path(cache_root))`, no `force`; docstring: GO is refreshed only by `refresh_go_release`)
- Modify: `.gitignore` (`cache/data/go_terms/raw/`)
- Test: `tests/test_refresh_go_release.py` (new), `tests/test_functional_annotation_adapter.py`

**Interfaces:**
- Produces: `refresh_go_release(cache_root: Path, force: bool = False) -> str` (returns R). Without `force` and with `raw/go-basic.obo`, `raw/go-plus.json`, `go_namespace_cache.json`, `go_release.txt` all present and agreeing → no download, returns R. With `force` (or any missing): download both current files to temp names; raise `GoReleaseMismatch` if their releases differ (nothing written); else move both into `raw/`, rebuild `go_namespace_cache.json` from the raw go-basic (`parse_obo`, same JSON format), write `go_release.txt`. CLI: `uv run python -m multiomics_kg.download.refresh_go_release [--force]`.
- `GO_BASIC_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"`, `GO_PLUS_URL = "http://purl.obolibrary.org/obo/go/extensions/go-plus.json"` (move `GO_PLUS_URL` here; Task 2 imports it).
- `build_namespace_cache` (lazy path, still used when the cache file is absent) also writes `go_release.txt` and keeps its OBO in `raw/` instead of deleting it.

- [ ] **Step 1: Write the failing tests** (fake downloader via monkeypatch): `test_refresh_writes_all_four_files_and_release` (both fakes stamped `2026-07-26` → returns it, files present, `read_go_release == "2026-07-26"`); `test_release_mismatch_writes_nothing` (go-basic `2026-07-26`, go-plus `2026-08-10` → raises `GoReleaseMismatch`, no file created); `test_no_force_no_download` (all files present → downloader never called); `test_no_cache_does_not_refetch_go` (adapter with `cache=False` → `load_go_data` called without `force=True`).
- [ ] **Step 2: Run** `uv run pytest tests/test_refresh_go_release.py tests/test_functional_annotation_adapter.py -v` — FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the same — PASS; then `uv run pytest -m "not slow and not kg" -q` — PASS.
- [ ] **Step 5: Refresh for real**: `uv run python -m multiomics_kg.download.refresh_go_release --force` → prints `2026-07-26` (or the then-current release — use whatever it returns throughout). Record `jq 'length'` of `go_namespace_cache.json` before (38,739) and after, and the added / removed term counts (`jq -r 'keys[]'` diff) in the scratchpad.
- [ ] **Step 6: Commit and push**

```bash
git add multiomics_kg/download/refresh_go_release.py multiomics_kg/utils/go_utils.py multiomics_kg/adapters/functional_annotation_adapter.py .gitignore tests/test_refresh_go_release.py tests/test_functional_annotation_adapter.py cache/data/go_terms/go_namespace_cache.json cache/data/go_terms/go_release.txt
git commit -m "feat(go): pin one GO release (go-basic + go-plus together); refresh to 2026-07-26"
git push
```

- [ ] **Step 7: Before-snapshot on the currently deployed (March-GO) graph**, taken before the user rebuilds: `/omics-edge-snapshot` (save) and GO counts per label (`MolecularFunction`, `BiologicalProcess`, `CellularComponent`) + the three gene→GO edge types, into the scratchpad.
- [ ] **Step 8: Rebuild** — the user pulls in `multiomics_biocypher_kg` and runs `docker compose up -d --build`; wait for confirmation that post-process finished.
- [ ] **Step 9: Verify the refresh build**: `/omics-edge-snapshot` compare (expression edges unchanged); GO counts after vs Step 7 (report the delta); `uv run pytest -m kg -v` — PASS (if `test_snapshot.py` fails only on GO terms the refresh obsoleted, regenerate `snapshot_data.json` and say so).
- [ ] **Step 10: Capture the baselines every later gate compares to**, on this rebuilt graph: `scripts/post-import-validate.sh > <scratchpad>/go_baseline.txt`; `uv run python tests/kg_validity/capture_annotation_state.py --save`; `/omics-edge-snapshot` (save); `MATCH (k:KeggTerm {level_kind:'pathway'}) RETURN count(k), sum(k.metabolite_count), sum(k.reaction_count)` and `MATCH (m:Metabolite) RETURN count(m)` into the scratchpad.
- [ ] **Step 11: Log it**: `CHANGELOG.md` `## [Unreleased]` → `### Changed`: "GO refreshed from the March 2026 release to 2026-07-26 (N terms added, M obsoleted; GO is now refreshed only with go-basic and go-plus together, never during a KG build)"; commit with `snapshot_data.json` if regenerated.

### Task 1: Shared "KG GO node set" utility (refactor, no behaviour change)

**Files:**
- Modify: `multiomics_kg/utils/go_utils.py` (add function after `compute_ancestry_closure`, ~line 233)
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:235-280` (`_all_seed_go_ids`, `all_go_terms`, `get_nodes`)
- Test: `tests/test_go_utils_kg_terms.py` (new), existing `tests/test_functional_annotation_adapter.py`

**Interfaces:**
- Produces: `kg_go_terms(go_data: dict[str, dict], genome_dirs: Iterable[Path]) -> dict[str, str]` — `{raw GO id ("GO:0015110"): schema label ("molecular function" | "biological process" | "cellular component")}`; seeds = every `go_terms` entry of each dir's merged annotations (via `multiomics_kg.utils.annotations_cache.load_merged_annotations`), closed with `compute_ancestry_closure`, filtered through `NAMESPACE_TO_LABEL`.
- `MultiGoAnnotationAdapter` gains `self.genome_dirs: list[Path]`; `all_go_terms()` returns `kg_go_terms(self.go_data, self.genome_dirs)`; `get_nodes()` iterates its keys.

- [ ] **Step 1: Write the failing tests**

`test_kg_go_terms_seeds_and_ancestors`: fixture `go_data` = `{"GO:1": {name, namespace "molecular_function", parents [["GO:0", "is_a"]]}, "GO:0": {..., parents []}, "GO:9": {namespace "biological_process", parents [["GO:1","regulates"]]}}` and two tmp genome dirs whose merged json carry `go_terms ["GO:1"]` and `["GO:9", "GO:404"]` → result `== {"GO:1": "molecular function", "GO:0": "molecular function", "GO:9": "biological process"}` (GO:404 absent from go_data is dropped; `regulates` parents are followed).
`test_adapter_all_go_terms_unchanged`: on the existing adapter test fixture, record `all_go_terms()` and the `get_nodes()` id set **before** the refactor (paste the literal expected dict/set into the test), assert both are equal after.

- [ ] **Step 2: Run** `uv run pytest tests/test_go_utils_kg_terms.py -v` — FAIL (`kg_go_terms` not defined).
- [ ] **Step 3: Implement `kg_go_terms`** and switch the adapter's three methods to it.
- [ ] **Step 4: Run** `uv run pytest tests/test_go_utils_kg_terms.py tests/test_functional_annotation_adapter.py tests/test_go_levels.py -v` — PASS.
- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/go_utils.py multiomics_kg/adapters/functional_annotation_adapter.py tests/test_go_utils_kg_terms.py tests/test_functional_annotation_adapter.py
git commit -m "refactor(go): one shared KG GO node-set function for adapter and step 6"
```

### Task 2: go-plus parser, target filter and committed derived cache

**Files:**
- Create: `multiomics_kg/utils/go_plus.py`
- Create: `config/go_transport_target_filter.yaml`
- Create: `tests/test_go_plus.py`, `tests/fixtures/go_plus_mini.json`

**Interfaces:**
- Produces (all in `go_plus.py`):
  - `parse_transport_axioms(graph: dict) -> dict[str, list[str]]` — `{"GO:0015110": ["CHEBI:29195"], ...}`, sorted lists, `GOCHE_n` mapped to `CHEBI:n`.
  - `go_plus_release(graph: dict) -> str` — the `owl#versionInfo` value.
  - `load_target_filter(path: Path) -> dict[str, str]` — `{"CHEBI:15378": "hydron", ...}`.
  - `build_go_transport_cache(cache_root: Path) -> dict` — reads Task 0's `cache_root/go_terms/raw/go-plus.json` (no download), checks `go_plus_release(graph) == read_go_release(cache_root)` (raise on mismatch: the raw file and the namespace cache must be one release), writes `cache_root/go_terms/go_transport_substrates.json` (`{"release": R, "md5": <md5 of the raw file>, "axioms": {...}}`, committed, indent=1, sorted keys) and returns it. If the raw file is absent (fresh checkout): call `refresh_go_release(cache_root)` (Task 0; no `force`) only when `go_release.txt` is also absent; if `go_release.txt` exists and the committed `go_transport_substrates.json` carries the same release, return the committed dict (offline rebuild); otherwise raise with "run refresh_go_release".
- YAML shape: `targets: [{chebi: "CHEBI:24870", name: "ion", reason: "<one line>"}, ...]` — the 13 ids in Global Constraints.

- [ ] **Step 1: Write the fixture** `tests/fixtures/go_plus_mini.json`: an obographs graph with `meta.basicPropertyValues` versionInfo `2026-07-26`; is_a edges GO_0022857→GO_0005215, GO_0015110→GO_0022857, GO_0015704→GO_0006810; logical definitions: GO_0015110 genus GO_0022857 + RO_0004009 CHEBI_29195; GO_0015704 genus GO_0006810 + RO_0004009 CHEBI_29195; GO_0005274 genus GO_0022857 + RO_0004009 CHEBI_15378 and CHEBI_15676; GO_0015891 genus GO_0006810 + RO_0004009 GOCHE_26672; GO_9999901 genus GO_0022857 + RO_0012001 CHEBI_29103 (activator); GO_9999902 genus GO_0022857 + RO_0004009 PR_000000001; GO_9999903 genus GO_0003824 (catalytic, not transport) + RO_0004009 CHEBI_15422.
- [ ] **Step 2: Write the failing tests**

```python
def test_parse_transport_axioms(mini):
    ax = parse_transport_axioms(mini)
    assert ax == {"GO:0015110": ["CHEBI:29195"], "GO:0015704": ["CHEBI:29195"],
                  "GO:0005274": ["CHEBI:15378", "CHEBI:15676"], "GO:0015891": ["CHEBI:26672"]}

def test_release(mini): assert go_plus_release(mini) == "2026-07-26"

def test_target_filter_contents():
    f = load_target_filter(Path("config/go_transport_target_filter.yaml"))
    assert set(f) == {f"CHEBI:{n}" for n in (24870, 24867, 23906, 23905, 25213, 33515, 51143,
                                              33694, 24431, 36080, 24834, 15378, 29101)}
    assert "CHEBI:29103" not in f and "CHEBI:17996" not in f        # K+, chloride kept

def test_build_cache_from_raw(tmp_path, mini): ...
    # raw go-plus = mini, go_release.txt = 2026-07-26 → file written, axioms as above, release matches
def test_build_cache_release_mismatch_raises(tmp_path, mini): ...
    # go_release.txt = 2026-03-01 → raises
def test_build_cache_offline_reuses_committed(tmp_path): ...
    # no raw file, go_release.txt + committed cache agree → returns committed dict, no download
```

- [ ] **Step 3: Run** `uv run pytest tests/test_go_plus.py -v` — FAIL (module missing).
- [ ] **Step 4: Implement** `go_plus.py` and the yaml (each entry's `reason` one line: "upper-level ChEBI class — transporter-ness, not a substrate" or, for H⁺/Na⁺, "coupling ion / bioenergetics — 5.5% / 4.0% of transport terms"). md5 over the raw file.
- [ ] **Step 5: Run** `uv run pytest tests/test_go_plus.py -v` — PASS.
- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/go_plus.py config/go_transport_target_filter.yaml tests/test_go_plus.py tests/fixtures/go_plus_mini.json
git commit -m "feat(go-plus): transport-axiom parser, target filter, committed derived cache"
```

### Task 3: Pure GO-target assembly helpers

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/go_transport.py`
- Test: `tests/test_go_transport_targets.py`

**Interfaces:**
- Consumes: `ChebiIndex` (`primary`, `term`), `role_node_id` from `metabolite_identity.roles`.
- Produces:
  - `go_substrate_strings(axioms: dict[str, list[str]], kg_terms: set[str], chebi: ChebiIndex) -> tuple[dict[str, list[tuple[str, str]]], list[str]]` — `({go_id: [(go_chebi_id, "CHEBI:<n>;<ChEBI name>")]}, skipped_go_ids)`; only `go_id ∈ kg_terms`; `skipped_go_ids` = KG terms with no axiom are NOT reported (normal) — it lists axioms whose ChEBI id `chebi.primary` cannot resolve.
  - `assemble_go_targets(term_strings, string_to_primary: dict[str, str], role_map: dict[str, str], role_ids: set[str], excluded_primaries: set[str]) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, list[str]]], Counter]` — `(go_substrates {go_id: {primary_id: [go_chebi_ids]}}, go_roles {go_id: {role_id: [go_chebi_ids]}}, drops)`; `drops` keys: `"filtered"`, `"role_without_node"`, `"unresolved"`.
  - `rewrite_go_targets(go_substrates, canonical_of: dict[str, str], excluded_canonical: set[str]) -> tuple[dict, int]` — rewrites target ids to canonical, merges `go_chebi_ids` lists (sorted, deduped) when two targets collapse, drops targets in `excluded_canonical`; returns `(rewritten, n_dropped_post_canonical)`.
  - `target_distribution(gene_terms_by_strain: dict[str, list[set[str]]], go_targets: dict[str, dict[str, list[str]]]) -> list[dict]` — rows `{"target", "genes", "gene_share", "terms", "term_share", "strains"}` sorted by `genes` desc.

- [ ] **Step 1: Write the failing tests** (a stub `ChebiIndex` built from a small terms dict, as in `tests/test_chebi_utils.py`):
  - `test_assemble_drops_filtered_and_keeps_k`: strings for H⁺, K⁺, cyanate; `excluded_primaries={"kegg.compound:C00080"}` → H⁺ absent, K⁺ and cyanate present, `drops["filtered"] == 1`.
  - `test_role_without_node_dropped`: role string mapping to `chebi.role:99999` not in `role_ids` → absent, `drops["role_without_node"] == 1`; `chebi.role:26672` in `role_ids` → present in `go_roles`.
  - `test_two_axioms_one_node_one_target`: term with CHEBI:29195 and its conjugate both resolving to `kegg.compound:C01417` → `go_substrates[go] == {"kegg.compound:C01417": ["CHEBI:28024", "CHEBI:29195"]}`.
  - `test_rewrite_merges_and_post_filter` (Review Focus 2): targets `A`, `B` with `canonical_of={"B": "A"}` → one target `A` with both id lists merged; a target whose canonical is in `excluded_canonical` → dropped, count 1.
  - `test_kg_term_without_go_plus_entry_is_skipped` (Review Focus 1): `kg_terms` includes `GO:0000000` absent from `axioms` → no entry, no exception.
  - `test_target_distribution`: two strains, three genes → exact rows.
- [ ] **Step 2: Run** `uv run pytest tests/test_go_transport_targets.py -v` — FAIL.
- [ ] **Step 3: Implement** the four functions.
- [ ] **Step 4: Run** — PASS.
- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/metabolite_identity/go_transport.py tests/test_go_transport_targets.py
git commit -m "feat(go-transport): pure GO-target assembly, canonical rewrite, distribution"
```

### Task 4: Step 6 wiring (stage D → identity pipeline → kegg_data block + report)

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`main` ~1441-1640; `_fold_substrates_into_kegg_data` ~1129)
- Modify: `multiomics_kg/utils/metabolite_identity/pipeline.py` (`Step6Inputs` ~50, `run_identity_stages` ~325-560)
- Test: `tests/test_build_kegg_metabolism_xrefs.py`, `tests/test_metabolite_identity_pipeline_go.py` (new; reuse `tests/metabolite_identity_helpers.py` builders)

**Interfaces:**
- Consumes: Task 1 `kg_go_terms`, Task 0 `refresh_go_release` / `read_go_release`, Task 2 `build_go_transport_cache` / `load_target_filter`, Task 3 helpers.
- `_fold_substrates_into_kegg_data(..., go_kegg_cpds: set[str] = frozenset(), go_primaries: set[str] = frozenset(), go_props: dict[str, dict] | None = None)` — appends `"go_transport"` to `evidence_sources` of compounds in `go_kegg_cpds`, adds non-KEGG `go_primaries` to `additional_compounds` (or appends the value to an existing entry).
- `Step6Inputs` new fields (defaulted so existing tests construct it unchanged): `go_substrates: dict = field(default_factory=dict)`, `go_roles: dict = field(default_factory=dict)`, `go_excluded_primaries: set = field(default_factory=set)`, `go_release: str | None = None`, `go_distribution_input: dict = field(default_factory=dict)`, `go_drops: dict = field(default_factory=dict)`.
- Produces in `kegg_data.json`: `go_transport_substrates: {go_id: {target_id: [go_chebi_ids]}}` (metabolite AND `chebi.role:` targets in one map), `go_transport_excluded_targets: [canonical ids]` (canonical node id when the filtered ChEBI id resolves to a written node, else `chebi:<n>`), `go_transport_release: str`. In the identity report: `report["go_transport"] = {"release", "axioms_in_corpus", "edges", "targets", "new_metabolites", "role_targets", "drops", "distribution": [...]}`.

- [ ] **Step 1: Write the failing tests**
  - `test_fold_tags_go_transport` (test_build_kegg_metabolism_xrefs.py): a compound in both `substrate_kegg_cpds` and `go_kegg_cpds` → `["transport", "go_transport"]`; a GO-only non-KEGG primary → `additional_compounds[...]["evidence_sources"] == ["go_transport"]`; an existing TCDB additional entry also in `go_primaries` → `["transport", "go_transport"]`.
  - `test_pipeline_rewrites_and_validates_go_targets` (new file): `Step6Inputs` from the helper builders with a `go_substrates` target that stage E merges into another node → output `kegg_data["go_transport_substrates"]` points at the canonical id with merged `go_chebi_ids`; a target not in the final table → `Step6ValidationError` naming `go_transport_substrates`.
  - `test_pipeline_no_pathway_pull_in`: a GO-only KEGG compound with a pathway outside `base_pathways` → its written entry's `pathways` excludes it, and `kegg_data["pathways"]` gains nothing.
- [ ] **Step 2: Run** the two test files — FAIL.
- [ ] **Step 3: Implement in `main`** following the Global-Constraints order: with `refetch_raw`, call `refresh_go_release(cache_root, force=True)` first (next to the TCDB download block; a `GoReleaseMismatch` or download failure with committed files present → warn and keep the committed release, TCDB outage style); then `go = build_go_transport_cache(cache_root)` **before** `_parse_raw_into_dict` / the resolver context; then, inside the resolver context after stage C and the TCDB `_resolve_substrate_strings` call and **before** `extended_sets` / `build_pruned_kegg_data`: `kg_terms = set(kg_go_terms(load_go_data(cache_root), [Path(r["data_dir"]) for r in load_genome_rows()]))`; `term_strings, _ = go_substrate_strings(go["axioms"], kg_terms, chebi)`; partition the distinct GO strings with `partition_role_strings`; resolve GO molecule strings and the filter's strings with `_resolve_substrate_strings` (same `chebi`, `rekeys`, `known_kegg`); `assemble_go_targets(..., role_ids=set(chemical_roles), excluded_primaries=<filter primaries>)`; add GO KEGG compounds to `extended_sets["cpds"]` only; pass `go_kegg_cpds` / `go_primaries` / props to `_fold_substrates_into_kegg_data`; collect `gene_terms_by_strain` (per gene, its GO term set ∩ GO-target terms) for the distribution; set the new `Step6Inputs` fields.
- [ ] **Step 4: Implement in `run_identity_stages`**: after stage F `rewrite_to_canonical`, `rewrite_go_targets(inp.go_substrates, identity.canonical_of, excluded_canonical)` where `excluded_canonical` = canonical of `inp.go_excluded_primaries`; roles pass through unchanged; add to `validate_build` inputs (every metabolite target ∈ table, every role ∈ `inp.chemical_roles`); write the three `kegg_data` keys and `report["go_transport"]` (distribution via `target_distribution`).
- [ ] **Step 5: Run** the two test files + `uv run pytest tests/test_build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs_metabolomics.py tests/test_post_import_metabolite_identity.py -q` — PASS.
- [ ] **Step 6: Run step 6 for real**: `bash scripts/prepare_data.sh --steps 6 --force` (no `--refetch-raw`: KEGG / TCDB raw stay as cached). `git diff --stat cache/data/go_terms/go_namespace_cache.json cache/data/go_terms/go_release.txt` → empty (Task 0 already refreshed GO; step 6 must not change it). Check `jq '.go_transport | {release, axioms_in_corpus, edges, targets, new_metabolites, drops}' cache/data/metabolomics/metabolite_identity_report.json` (release `2026-07-26`; edges > 0; drops.filtered > 0); `jq '.go_transport_substrates["GO:0015110"]' cache/data/kegg/kegg_data.json` → `{"kegg.compound:C01417": ["CHEBI:29195"]}`; `jq '.pathways | length' cache/data/kegg/kegg_data.json` equals the value before the run (record it first with `git show HEAD:cache/data/kegg/kegg_data.json | jq '.pathways|length'`).
- [ ] **Step 7: Commit** code, tests and the regenerated committed caches (review `git status` first; stage the step-6 outputs by path: `kegg_data.json`, `go_terms/go_transport_substrates.json`, the identity report, mapping, `tcdb_pruned.json` only if changed).

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py multiomics_kg/utils/metabolite_identity/pipeline.py tests/test_build_kegg_metabolism_xrefs.py tests/test_metabolite_identity_pipeline_go.py cache/data/kegg/kegg_data.json cache/data/go_terms/go_transport_substrates.json cache/data/metabolomics/metabolite_identity_report.json cache/data/metabolomics/metabolite_id_mapping.json
git commit -m "feat(step6): GO transport substrates resolved, filtered, canonicalised into kegg_data"
```

### Task 5: Schema + adapter emission + build wiring

**Files:**
- Modify: `config/schema_config.yaml` (4 edge types next to `tcdb family transports metabolite` ~1474; `metabolite:` gains `go_transport_organism_names: str[]`, `go_transport_organism_count: int` with `# post-import` comments)
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py` (`MultiGoAnnotationAdapter`)
- Modify: `create_knowledge_graph.py` (~line 231-238 and ~463)
- Test: `tests/test_functional_annotation_adapter.py`

**Interfaces:**
- Consumes: Task 4 `kegg_data.json` keys; Task 1 `all_go_terms()`.
- Produces: `MultiGoAnnotationAdapter.get_transport_edges(kegg_data_path: Path, metabolite_node_ids: set[str] | None, chemical_role_node_ids: set[str] | None) -> Iterator[tuple[None, str, str, str, dict]]` — yields `(None, go_node_id, target_id, label, {"go_chebi_ids": [...]})` with `label_in_input` ∈ `molecular_function_transports_metabolite | biological_process_transports_metabolite | molecular_function_transports_chemical_role | biological_process_transports_chemical_role` (by the term's `all_go_terms()` label and the target's prefix `chebi.role:`); skips terms not in `all_go_terms()` and targets not in the injected sets; `None` set → no edges of that target kind.
- Schema `label_as_edge` values exactly as in Global Constraints; `properties: {go_chebi_ids: str[]}`; comment on each: "GO defines this term's primary input as X — a term definition, not a gene fact; BP is broader than MF".

- Also: `get_transport_edges` logs a warning when `cache_root/go_terms/go_release.txt` ≠ `kegg_data["go_transport_release"]`.

- [ ] **Step 1: Write the failing tests**: `test_release_mismatch_warns` (caplog). Then: tmp `kegg_data.json` with `go_transport_substrates` for an MF term → metabolite, a BP term → `chebi.role:26672`, a term not in the adapter's GO set, and a metabolite not in `metabolite_node_ids` (Review Focus 4) → exactly two edges with the right labels and `go_chebi_ids`; `chemical_role_node_ids=None` → only the metabolite edge.
- [ ] **Step 2: Run** — FAIL.
- [ ] **Step 3: Implement** the method and schema entries; in `create_knowledge_graph.py` after `bc.write_edges(go_anno_adapter.get_edges())`: `bc.write_edges(go_anno_adapter.get_transport_edges(Path("cache/data/kegg/kegg_data.json"), metabolism_adapter.all_metabolite_ids(), set(json.loads(Path("cache/data/tcdb/tcdb_pruned.json").read_text())["chemical_roles"])))`; next to `chebi_release.txt`, write `go_transport_excluded_targets.json` (the kegg_data list) into `bc._output_directory`.
- [ ] **Step 4: Run** `uv run pytest tests/test_functional_annotation_adapter.py -v` and `uv run python create_knowledge_graph.py --test` (completes; the four edge CSVs exist in the output dir).
- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml multiomics_kg/adapters/functional_annotation_adapter.py create_knowledge_graph.py tests/test_functional_annotation_adapter.py
git commit -m "feat(go-transport): emit GO term → Metabolite / ChemicalRole edges"
```

### Task 6: Post-import rollup + Schema_info stamp + vocabulary

**Files:**
- Modify: `scripts/post-import.sh` (Metabolite rollup block after the KG-A8 `organism_names` recompute ~1973-1985; Group 4 ~2136-2215), `scripts/post-import.cypher` (same blocks, identical)
- Modify: `config/controlled_vocabularies.yaml` (`Metabolite.evidence_sources` description: add "go_transport (named as the primary input of a GO transport term a gene carries; a term definition, not a gene fact)")
- Test: `tests/test_post_import_go_transport.py` (new, static)

**Interfaces:**
- Cypher (both files):

```cypher
MATCH (m:Metabolite)
WHERE EXISTS { (:MolecularFunction|BiologicalProcess)-[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m) }
CALL {
  WITH m
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)
        -[:Gene_enables_molecular_function|Gene_involved_in_biological_process]->(t)
        -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m)
  WITH m, collect(DISTINCT o.preferred_name) AS names
  SET m.go_transport_organism_names = apoc.coll.sort(names),
      m.go_transport_organism_count = size(names)
} IN TRANSACTIONS OF 1000 ROWS;
```

- Group 4: the shell reads `/data/build2neo/go_transport_excluded_targets.json` (same flat/timestamped `ls -t` glob as `chebi_release.txt`), passes it as `-P "go_excluded => <list literal>"`, and the Schema_info `SET` adds `s.go_transport_excluded_targets = coalesce($go_excluded, [])`; warn when the file is missing.

- [ ] **Step 1: Write the failing test** `test_go_transport_rollup_identical_in_both_files`: extract the block between `// GO transport organism rollup` markers from both files → equal and present; `test_schema_info_sets_excluded_targets`: `.sh` contains `s.go_transport_excluded_targets`.
- [ ] **Step 2: Run** — FAIL.
- [ ] **Step 3: Implement** both files (copy the block verbatim) and the vocabulary edit.
- [ ] **Step 4: Run** `uv run pytest tests/test_post_import_go_transport.py tests/test_controlled_vocab.py -v` — PASS; `diff <(sed -n '/GO transport organism rollup/,/IN TRANSACTIONS/p' scripts/post-import.sh) <(sed -n '/GO transport organism rollup/,/IN TRANSACTIONS/p' scripts/post-import.cypher)` — empty.
- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher config/controlled_vocabularies.yaml tests/test_post_import_go_transport.py
git commit -m "feat(post-import): GO-route organism rollup + Schema_info exclusion list"
```

### Task 7: KG validity tests, rebuild and before/after gates

**Files:**
- Create: `tests/kg_validity/test_go_transport.py`
- Modify: `tests/kg_validity/snapshot_data.json` (regenerated), `tests/kg_validity/generate_snapshot.py` if it needs a GO-transport edge sampler

**Interfaces:**
- Consumes: the live graph at bolt://localhost:7687 (auto-skip when absent; `@pytest.mark.kg`).

- [ ] **Step 1: Baselines** are Task 0 Step 10's captures (taken on the GO-refreshed graph); re-capture only if another change was deployed in between.
- [ ] **Step 2: Write the kg tests**: four edge types each `> 0`; every edge's source is `MolecularFunction|BiologicalProcess` and target exists; `go:0015110` → `kegg.compound:C01417` with `'CHEBI:29195' IN r.go_chebi_ids`; `go:0015891` (siderophore transport) → `chebi.role:26672`; no GO transport edge targets an id in `Schema_info.go_transport_excluded_targets`, which contains the H⁺ and Na⁺ node ids; every edge has `size(go_chebi_ids) > 0`; PMM0402 (fadD) reaches no metabolite via the route; MED4 genes reaching ≥ 1 metabolite via the route ∈ [45, 75] (measured 59); `go_transport_organism_count = size(go_transport_organism_names)` on every carrier and absent elsewhere; pathway `KeggTerm` count and `sum(reaction_count)` equal the Step 1 values; `sum(metabolite_count)` equals the Step 1 value plus the GO-only compounds' links into kept pathways (6 measured in Task 4; read the exact number from the Task 4 report) — paste them as constants.
- [ ] **Step 3: Push and rebuild.** `git push`; the user pulls in `multiomics_biocypher_kg` and runs `docker compose up -d --build` (~1 h). Wait for the user's confirmation that post-process finished.
- [ ] **Step 4: Run gates**: `/omics-edge-snapshot` compare (expression edges unchanged); `scripts/post-import-validate.sh > after.txt` and `diff` against the baseline — additions only (the new properties / counts); `capture_annotation_state.py --compare` — no change; `uv run pytest -m kg -v` — PASS including `test_go_transport.py`; `cat output/import.report` shows no skipped GO-transport relationships.
- [ ] **Step 5: Regenerate the snapshot** `uv run python tests/kg_validity/generate_snapshot.py`, check it contains a GO-transport edge sample, rerun `uv run pytest tests/kg_validity/test_snapshot.py -v`.
- [ ] **Step 6: Commit**

```bash
git add tests/kg_validity/test_go_transport.py tests/kg_validity/snapshot_data.json tests/kg_validity/generate_snapshot.py
git commit -m "test(kg): GO transport edge validity + snapshot refresh"
```

### Task 8: Docs, changelog and explorer hand-off

**Files:**
- Create: `docs/kg-changes/go-transport-substrates.md`
- Modify: `CLAUDE.md`, `CHANGELOG.md`, `plans/backlog.md`, `multiomics_explorer/docs/kg-specs/2026-10-18-go-transport-substrates-explorer-asks.md`

- [ ] **Step 1: Write the kg-changes doc** from the live numbers (Task 7) and the step-6 report: what/why; the four edge types with counts; `go_chebi_ids`; the target filter table (13 ids, names, reasons) and drop counts; the target distribution top 25; the two new Metabolite properties and `Schema_info.go_transport_excluded_targets`; **one ready-to-copy caveat paragraph per edge type** (MF = "is a transporter for X"; BP = "involved in moving X", includes storage/regulation/accessory, ferritin example; "a term definition, not a gene fact"; "filtered targets return empty by curation"); the release date `2026-07-26`; query examples (gene → GO → metabolite; class reach via GO `is_a` and via `Metabolite_in_metabolite_class`); out-of-scope (direction, taxon constraints, Rhea bridge).
- [ ] **Step 2: Update `CLAUDE.md`**: relationships list (+4 labels), a Key-graph-facts bullet linking the doc, Data Locations entries (`cache/data/go_terms/go_transport_substrates.json`, `go_release.txt`, `go_namespace_cache.json` now refreshed by step 6, raw `go-plus.json` / `go-basic.obo`, `config/go_transport_target_filter.yaml`), a GO paragraph (`refresh_go_release`: go-basic + go-plus downloaded together, release-checked; step 6 refreshes only with `--refetch-raw`) + the step-6 GO-transport sub-step + "the KG build no longer refetches GO (`--no-cache` does not touch GO)", the CLI-flags table note for `--no-cache`, `Metabolite.evidence_sources` values, Metabolite post-import properties.
- [ ] **Step 3: `CHANGELOG.md` `## [Unreleased]`**: `### Added` (edges + counts, properties, Schema_info list, doc link); `### Highlights` ("Which genes carry a GO term defined as moving compound X — e.g. siderophores — is now a graph query, on the same metabolite nodes TCDB and KEGG use."); no `### Breaking` (additive; note the vocabulary hash moves).
- [ ] **Step 4: `plans/backlog.md`**: mark slice 1 done under the transport programme; file follow-ups (GO taxon-constraint flag, GO MF ↔ Rhea/EC/KEGG-reaction bridge, direction gap).
- [ ] **Step 5: Explorer asks doc**: append "KG built (date): dev graph carries slice 1; counts; doc path; vocabulary hash moved; GO release moved to 2026-07-26" so the explorer can absorb. Record in `plans/backlog.md` that the KG release carrying slice 1 sets `[tool.release-kg] mcp_min_version` in `pyproject.toml` to the explorer release that surfaces the GO arm (joint release, design §8) — done at `/release-kg` time, not in this plan.
- [ ] **Step 6: Run** `uv run pytest -m "not slow and not kg" -q` — PASS.
- [ ] **Step 7: Commit** (explorer repo separately, its doc path only)

```bash
git add docs/kg-changes/go-transport-substrates.md CLAUDE.md CHANGELOG.md plans/backlog.md
git commit -m "docs(go-transport): kg-changes doc, CLAUDE.md, CHANGELOG, backlog"
```
