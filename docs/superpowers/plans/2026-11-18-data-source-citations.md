# DataSource Citations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every `DataSource` node carries the citations its provenance unit needs, and `Schema_info.build_citations` carries the build framework's, so an explorer MCP tool can return `{source_id: [citations]}` from one query.

**Architecture:** A committed `config/data_source_citations.yaml` becomes the single source of per-source metadata (name, description, citations, dois, url, license, version). A small loader in `multiomics_kg/utils/data_source_citations.py` reads it and resolves versions; `DataSourceAdapter` consumes it (its hard-coded name/description dicts are deleted). Six reference databases join `reference_sources`. BioCypher's citation goes to `Schema_info` through the existing build-sidecar → post-import stamp pattern (`chebi_release.txt`).

**Tech Stack:** Python 3, PyYAML, BioCypher, pytest, cypher-shell (post-import).

**Spec:** `docs/superpowers/specs/2026-11-18-data-source-citations-design.md`

## Global Constraints

- Approach A only: no `applies_to` mapping, no DataSource→label edges.
- DataSource stays a provenance target: build software (BioCypher, Neo4j) is **not** a DataSource. Neo4j gets no citation.
- No DataSource nodes for DIAMOND, InterPro, NCBIfam, TIGRFAMs — they are entries in the `citations` lists of `tcdb_diamond` / `merops_diamond` / `interproscan`.
- New `reference_sources` ids exactly: `go`, `pfam`, `cazy`, `expasy_enzyme`, `metanetx`, `ncbi_taxonomy` — all `scope: global`, `provenance: download`. PubChem / Crossref only if their data lands in the graph (Task 2 checks).
- Curated text (every yaml string) contains no `'`, `"`, `|` or `\` (`_clean_str` maps `'`→`^` and deletes `|`; post-import swaps `"`→`'` when injecting JSON lists).
- `dois` is parallel to `citations` (same length, `""` where a reference has no DOI); DOIs bare (`10.xxxx/...`, no `https://doi.org/`).
- Every citation is taken from the tool's/database's own "how to cite" page — never from memory alone. Record the page URL in the yaml entry's `url`.
- Versions: `version_from: chebi_cache` (existing behaviour) and `version_from: go_release` (`cache/data/go_terms/go_release.txt`); literals for the rest; `interproscan` literal `"111.0"` must equal `INTERPRO_RELEASE` in `.claude/skills/interproscan-run/run_interproscan.py`; `merops_diamond` literal `"12.5"`.
- `Schema_info` is in `UNSCHEMED_LABELS` (`tests/kg_validity/test_schema_coverage.py:28`), so `build_citations` needs no `schema_config.yaml` entry (corrects spec §4's "declared in schema_config.yaml").
- Concurrent sessions edit this clone: `git add` explicit paths only.

## Review Focus

1. A citation string with an apostrophe (e.g. "Nature's") or a double quote → must fail the unit test, not reach the CSV as `^` or break the post-import `-P` list. (Task 1 test `test_curated_text_has_no_forbidden_chars`.)
2. A DataSource id with no yaml entry, or a yaml entry with no DataSource → test failure naming the id. (Task 1 `test_yaml_and_data_sources_match`.)
3. Fresh checkout with no `go_release.txt` / no ChEBI cache → adapter emits `version = ""` with a warning, never raises. (Task 1 `test_version_from_missing_file_is_empty`.)
4. Post-import with no `build_citations.json` (old build output) → `Schema_info.build_citations = []` with a warning, import continues. (Task 3 step: bash guard mirrors `go_excluded`.)
5. InterProScan rescan bumps `INTERPRO_RELEASE` but the yaml isn't updated → unit test fails. (Task 1 `test_interpro_version_matches_runner_pin`.)

## Citation research (used by Tasks 1, 2, 3)

For **every** yaml entry (each source and `build_citations`):
1. `WebSearch` (mode `standard`; `extended` if the result is thin) for `"<tool/database> how to cite"` / `"<tool> citation"`; open the project's own citation page with `WebFetch` (e.g. the InterPro "Citing InterPro" page, eggNOG-mapper's README, TCDB's "Citing TCDB").
2. Take the reference(s) the project itself asks for — when it names the most recent paper for the version we use, prefer that one. Format: `Authors (first author et al. when > 3). Year. Title. Journal Volume:Pages.` — no quotes.
3. Confirm each DOI by searching it (`"10.xxxx/..."`) and checking that title + year match.
4. Set `url` to the citation page you used.
5. Keep a working list (scratchpad, not committed) of `id → page URL + DOIs` for Task 2b.

---

### Task 1: Citations yaml + loader + adapter refactor (existing 11 sources)

**Files:**
- Create: `config/data_source_citations.yaml`
- Create: `multiomics_kg/utils/data_source_citations.py`
- Modify: `multiomics_kg/adapters/data_source_adapter.py` (delete `_name_for`, `_description_for`, `_version_for`, `_CHEBI_CACHE`; emit new props)
- Modify: `config/schema_config.yaml:494-507` (`data source` block)
- Test: `tests/test_data_source_adapter.py`, create `tests/test_data_source_citations.py`

**Interfaces:**
- Produces: `load_citations(path: str | Path = "config/data_source_citations.yaml") -> dict` returning `{"sources": {id: entry}, "build_citations": [str]}`; `resolve_version(entry: dict) -> str`; `FORBIDDEN_CHARS = ("'", '"', "|", "\\")`.
- Produces: DataSource node props `name`, `description`, `version`, `citations: list[str]`, `dois: list[str]`, `url: str`, `license: str` (omitted when absent).

yaml shape:
```yaml
build_citations:            # filled in Task 3
  - "..."
sources:
  eggnog:
    name: EggNOG-mapper
    description: "..."
    citations: ["<eggNOG-mapper v2 ref>", "<eggNOG 5 ref>", "<COG ref (Galperin)>"]
    dois: ["10.1093/molbev/msab293", "...", "..."]
    url: "<how-to-cite page>"
    version: "5.0.2"          # or version_from: go_release | chebi_cache
    license: "..."            # optional
```
Citation lists for the existing 11 per spec §2: `eggnog` (eggNOG-mapper v2, eggNOG 5, COG); `tcdb_diamond` (TCDB, DIAMOND); `merops_diamond` (MEROPS, DIAMOND); `interproscan` (InterProScan, InterPro, NCBIfam, TIGRFAMs — name "TIGRFAMs 15.0" in that citation's text); `kegg` (KEGG's requested papers; description widened to KO / pathways / BRITE / COMPOUND); `ncbi` (RefSeq/PGAP); `cyanorak`; `uniprot`; `psortb`; `signalp` (SignalP 6.0); `chebi` (keep "CC BY 4.0" in description, `license: CC BY 4.0`). Fix `interproscan` description to "InterProScan 6". Names unchanged from today's `_name_for` (tests assert "MEROPS (diamond)", "InterProScan").

- [ ] **Step 1: Write failing tests** in `tests/test_data_source_citations.py`:
  - `test_yaml_and_data_sources_match`: `set(load_citations()["sources"]) == data_source_ids(yaml.safe_load(open("config/gene_annotations_config.yaml")))`.
  - `test_every_entry_complete`: for each entry, `name`, `url` non-empty; `citations` non-empty list of non-empty strings; `len(dois) == len(citations)`; each non-empty doi matches `^10\.\d{4,9}/\S+$`.
  - `test_curated_text_has_no_forbidden_chars`: walk every string in the yaml (incl. `build_citations`); none contains any of `FORBIDDEN_CHARS`.
  - `test_interpro_version_matches_runner_pin`: regex `INTERPRO_RELEASE = "([^"]+)"` over `.claude/skills/interproscan-run/run_interproscan.py`; equals `sources["interproscan"]["version"]`.
  - `test_version_from_missing_file_is_empty(tmp_path, monkeypatch)`: point the go_release reader at a missing path → `resolve_version({"version_from": "go_release"}) == ""`; unknown `version_from` → `ValueError`.
  - `test_go_release_reader`: when `cache/data/go_terms/go_release.txt` exists, `resolve_version({"version_from": "go_release"})` equals its stripped content.
- [ ] **Step 2: Update `tests/test_data_source_adapter.py`**: every node has `citations` (non-empty list) and `url`; `nodes["kegg"]["version"] == ""` → assert it equals the yaml's kegg `version` (literal or `""`); keep the chebi-cache assertion. `test_emits_eleven_nodes` stays as-is in this task (Task 2 changes it).
- [ ] **Step 3: Run** `pytest tests/test_data_source_citations.py tests/test_data_source_adapter.py -v` → FAIL (module missing).
- [ ] **Step 4: Implement** the loader (module-level path constants for the two caches so tests can monkeypatch; missing file → `logger.warning` + `""`) and the yaml for the 11 sources. **Fill each entry by web search** (procedure in "Citation research" below); never write a citation from memory.
- [ ] **Step 5: Refactor the adapter**: `download_data` also loads citations; `get_nodes` takes name/description/version/citations/dois/url/license from the entry, `_clean_str` on every string and list item; raise `KeyError` naming the id if a DataSource has no entry.
- [ ] **Step 6: Schema**: add to `data source` properties `citations: str[]`, `dois: str[]`, `url: str`, `license: str` (comment: sparse); change the `scope` comment to `'gene_level' | 'organism_restricted' | 'global'`.
- [ ] **Step 7: Run** `pytest tests/test_data_source_citations.py tests/test_data_source_adapter.py tests/test_controlled_vocab.py -v` → PASS.
- [ ] **Step 8: Commit**
```bash
git add config/data_source_citations.yaml multiomics_kg/utils/data_source_citations.py multiomics_kg/adapters/data_source_adapter.py config/schema_config.yaml tests/test_data_source_citations.py tests/test_data_source_adapter.py
git commit -m "feat(data-source): citations yaml + loader; DataSource carries citations/dois/url/license"
```

### Task 2: Six new reference sources

**Files:**
- Modify: `config/gene_annotations_config.yaml:163-171` (`reference_sources`)
- Modify: `config/data_source_citations.yaml`
- Test: `tests/test_data_source_adapter.py`

**Interfaces:**
- Consumes: Task 1 loader + adapter (no code change expected).

New entries (all `scope: global`, `provenance: download`) with `info_types_extra`:
| id | name | info_types_extra | version |
|---|---|---|---|
| `go` | Gene Ontology | `[go_terms, go_hierarchy, go_transport_substrates]` | `version_from: go_release` |
| `pfam` | Pfam | `[pfam_domains, pfam_clans]` | literal (release of `Pfam-A.clans.tsv.gz` in use) |
| `cazy` | CAZy | `[cazy_families]` | `""` |
| `expasy_enzyme` | Expasy ENZYME | `[ec_numbers, ec_hierarchy]` | `""` unless `ec_data.json` records one |
| `metanetx` | MetaNetX | `[reaction_xrefs, metabolite_xrefs]` | literal MNX release used by `scripts/refresh_mnx.sh` |
| `ncbi_taxonomy` | NCBI Taxonomy | `[organism_taxonomy]` | `""` |

- [ ] **Step 1: Update tests**: rename `test_emits_eleven_nodes` → `test_emits_seventeen_nodes`, exact id set = the 11 + the 6 above; extend `test_reference_sources_have_info_types_and_chebi_version` loop to all 8 reference ids.
- [ ] **Step 2: Run** `pytest tests/test_data_source_adapter.py tests/test_data_source_citations.py -v` → FAIL.
- [ ] **Step 3: Check PubChem / Crossref**: read `multiomics_kg/utils/metabolite_identity/pipeline.py` + `__init__.py`. If any PubChem- or Crossref-derived value is written to a node/edge property, add that source too (same pattern) and to the test's id set; otherwise note "not cited: lookup only" in the commit message.
- [ ] **Step 4: Add** the `reference_sources` entries and yaml citation entries, each filled by web search ("Citation research" below).
- [ ] **Step 5: Run** `pytest tests/test_data_source_adapter.py tests/test_data_source_citations.py tests/test_controlled_vocab.py -v` → PASS.
- [ ] **Step 6: Commit**
```bash
git add config/gene_annotations_config.yaml config/data_source_citations.yaml tests/test_data_source_adapter.py
git commit -m "feat(data-source): GO, Pfam, CAZy, Expasy ENZYME, MetaNetX, NCBI Taxonomy as reference sources"
```

### Task 2b: Independent citation review

**Files:**
- Modify: `config/data_source_citations.yaml` (corrections only)

Run after Task 2 (and re-run on `build_citations` after Task 3 Step 0).

- [ ] **Step 1: Dispatch a fresh review agent** (Agent tool, `general-purpose`, no context from the implementer). Give it only: the path `config/data_source_citations.yaml`, the KG's versions in use (from the yaml `version` fields), and this brief: *"For every entry, independently web-search the tool/database's official citation page and verify: (a) the citations listed are the ones the project asks users to cite, for the version we use; nothing required is missing (e.g. a database paper alongside a tool paper); (b) authors, year, title, journal, volume, pages are correct; (c) each DOI resolves to that paper (search the DOI). Do not trust the yaml's `url`. Report per entry: OK, or each discrepancy with the source URL that shows the correct value. Do not edit files."*
- [ ] **Step 2: Triage** the report: fix confirmed discrepancies in the yaml; for disagreements, re-check the cited source yourself and record the call in the commit message.
- [ ] **Step 3: Run** `pytest tests/test_data_source_citations.py tests/test_data_source_adapter.py -v` → PASS.
- [ ] **Step 4: Commit**
```bash
git add config/data_source_citations.yaml
git commit -m "fix(data-source): citation corrections from independent review"
```

### Task 3: `Schema_info.build_citations`

**Files:**
- Modify: `config/data_source_citations.yaml` (`build_citations`: BioCypher, Lobentanzer et al. 2023, Nat Biotechnol)
- Modify: `create_knowledge_graph.py:~466` (next to `chebi_release.txt`)
- Modify: `scripts/post-import.sh:~2122-2181`, `scripts/post-import.cypher:~2040-2061`
- Test: `tests/test_data_source_citations.py`, `tests/kg_validity/test_data_source_citations.py` (create)

**Interfaces:**
- Consumes: `load_citations()["build_citations"]`.
- Produces: sidecar `<bc._output_directory>/build_citations.json` (JSON list of strings); `Schema_info.build_citations: str[]`.

- [ ] **Step 0: Fill `build_citations`** by web search ("Citation research"), then have the Task 2b reviewer check it (send it the one entry).
- [ ] **Step 1: Unit test** `test_build_citations_nonempty`: `load_citations()["build_citations"]` is a non-empty list, one entry mentions `BioCypher`.
- [ ] **Step 2: Write** the sidecar in `create_knowledge_graph.py` with `json.dumps(load_citations()["build_citations"])`.
- [ ] **Step 3: post-import.sh**: `KG_BUILD_CITATIONS="[]"`, locate `build_citations.json` with the same flat/one-level `ls -t` lookup as `go_transport_excluded_targets.json`, `tr -d '\n' | tr '"' "'"` (safe because Task 1 forbids `'`/`"` in curated text), warning when missing; `-P "build_citations => ${KG_BUILD_CITATIONS}"`; `s.build_citations = coalesce($build_citations, [])`.
- [ ] **Step 4: post-import.cypher**: `:param build_citations => []` and the same `SET` line (both files must carry identical logic).
- [ ] **Step 5: KG test** `tests/kg_validity/test_data_source_citations.py` (`@pytest.mark.kg`, existing `run_query` fixture pattern): `Schema_info.build_citations` non-empty; every `DataSource` has non-empty `citations` and `size(dois) = size(citations)`; `count(DataSource) == 17` (or 17 + PubChem/Crossref from Task 2).
- [ ] **Step 6: Run** `pytest tests/test_data_source_citations.py -v` → PASS; `bash -n scripts/post-import.sh` → no output.
- [ ] **Step 7: Commit**
```bash
git add config/data_source_citations.yaml create_knowledge_graph.py scripts/post-import.sh scripts/post-import.cypher tests/test_data_source_citations.py tests/kg_validity/test_data_source_citations.py
git commit -m "feat(schema-info): build_citations (BioCypher) stamped post-import"
```

### Task 4: Docs, full unit suite, rebuild verification

**Files:**
- Modify: `CLAUDE.md` (DataSource bullet: 17 nodes, the new ids, `citations`/`dois`/`url`/`license`, `config/data_source_citations.yaml` as source of truth, `Schema_info.build_citations`; Docker post-process paragraph lists the new sidecar)
- Modify: `CHANGELOG.md` (`## [Unreleased]` → `### Added`: one bullet)

- [ ] **Step 1: Edit** both files.
- [ ] **Step 2: Run** `pytest -m "not slow and not kg" -q` → all pass (report any pre-existing failures separately, unchanged).
- [ ] **Step 3: Commit**
```bash
git add CLAUDE.md CHANGELOG.md
git commit -m "docs: DataSource citations + Schema_info.build_citations"
```
- [ ] **Step 4: Rebuild + KG tests** — Docker build runs in the other clone (`multiomics_biocypher_kg`): push, user pulls + `docker compose up -d`, then `pytest tests/kg_validity/test_data_source_citations.py tests/kg_validity/test_schema_coverage.py -v` → PASS. Record the result here.
