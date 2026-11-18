# OrganismTaxon Assembly Quality Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every genome-backed `OrganismTaxon` node carries NCBI's assembly level, RefSeq status, contig stats and CheckM result.

**Architecture:** Step 0 sub-step 1 fetches each strain's NCBI Datasets v2 `dataset_report` into a committed `<strain dir>/assembly_report.json` (`{"fetched_at", "report"}`). A pure parser in a new `multiomics_kg/utils/assembly_report.py` turns one file into a property dict; `cyanorak_ncbi_adapter._get_organism_node` merges it. Nine `ControlledVocabulary` entries and 16 schema declarations publish the properties.

**Tech Stack:** Python 3, `requests`, pytest, BioCypher, Neo4j (Docker, other clone).

**Spec:** `docs/superpowers/specs/2026-11-18-organism-assembly-quality-design.md`

## Global Constraints

- Branch `feat/organism-assembly-quality` (spec already committed there). Stage explicit paths only — other sessions edit this clone (`git add -A` forbidden).
- Endpoint: `https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc,acc,…}/dataset_report?filters.assembly_version=all_assemblies&page_size=100`, GET, chunks of 40 accessions.
- File: `<data_dir>/assembly_report.json`, `json.dump(..., indent=1)`, shape `{"fetched_at": "YYYY-MM-DD", "report": <one NCBI report, verbatim>}`. Committed.
- Only accessions starting `GC` are fetched (registry has comment rows).
- Strings go through `_clean_str` semantics (`'`→`^`, `|` removed) before yield.
- Properties are absent (never `""`, never `None`) when NCBI gives no value; the one exception is `checkm_contamination = 0.0` when a `checkm_info` block exists without `contamination`.
- `assembly_level` / `assembly_status` values not in the vocabulary raise at build time via `VOCAB.check("OrganismTaxon", prop, value)`.
- Docker build/KG tests run in the other clone (`multiomics_biocypher_kg`): commit + push here, user pulls + rebuilds there.

## Review Focus

1. A strain with no `assembly_report.json` (new strain, fetch failed) → organism node built without assembly properties, adapter logs a warning, build does not fail. (Task 3 test `test_missing_report_file`.)
2. NCBI returns only a different version than requested (e.g. `.2` for a `.1` request) → no file written for the requested accession, warning names it. (Task 1 test `test_no_exact_accession_match`.)
3. `suppression_reason` without a `Reason:` clause → `"removed by RefSeq staff"`, not the boilerplate or empty. (Task 2 test.)
4. `--strains MED4` → only MED4's report is fetched/written. (Task 1 test `test_respects_genome_subset`.)
5. Second run without `--force` → no HTTP call for strains whose file exists. (Task 1 test `test_skips_existing_without_force`.)

---

### Task 1: Fetch assembly reports in step 0

**Files:**
- Modify: `multiomics_kg/download/download_genome_data.py` (new function near `_ncbi_download_genome`; call at end of `step1_ncbi`)
- Test: `tests/test_download_genome_data.py` (new `class TestFetchAssemblyReports`)

**Interfaces:**
- Produces: `ASSEMBLY_REPORT_FILENAME = "assembly_report.json"` (module constant in `multiomics_kg/utils/assembly_report.py` — create that module here with just the constant; Task 2 fills it). `_ncbi_fetch_assembly_reports(genomes: list[dict], force: bool, today: str | None = None) -> dict[str, str]` returning `{accession: "written" | "skipped" | "missing"}`.

- [ ] **Step 1: Write failing tests** (patch `requests.get`; `genomes` fixture with `data_dir` pointed at `tmp_path` subdirs)
  - `test_writes_wrapped_report`: mocked response `{"reports": [{"accession": "GCF_000011465.1", "assembly_info": {...}}]}` → file exists, `json["fetched_at"] == "2026-11-18"` (pass `today="2026-11-18"`), `json["report"]["accession"] == "GCF_000011465.1"`, status `"written"`.
  - `test_url_has_all_assemblies_filter_and_chunks`: 45 GC accessions → `requests.get` called twice; first URL contains 40 comma-joined accessions and `filters.assembly_version=all_assemblies`.
  - `test_no_exact_accession_match`: request `GCF_000179235.1`, response holds only `GCF_000179235.2` → no file, status `"missing"`, warning logged containing `GCF_000179235.1`.
  - `test_skips_existing_without_force`: file pre-exists → `requests.get` not called, status `"skipped"`; with `force=True` it is called and overwritten.
  - `test_respects_genome_subset`: genomes list of one strain → only that strain's file written.
  - `test_ignores_non_gc_accessions`: a row whose `ncbi_accession` is `"# comment"` is not in any URL.
  - `test_http_error_warns_and_continues`: `raise_for_status` raises on chunk 1 → statuses `"missing"` for that chunk, no exception.

- [ ] **Step 2: Run** `uv run pytest tests/test_download_genome_data.py -k AssemblyReports -v` — Expected: FAIL (ImportError).

- [ ] **Step 3: Implement** `_ncbi_fetch_assembly_reports` with `requests.get(url, timeout=120)`, `response.raise_for_status()`, match reports on `report["accession"] == requested` exactly; `today` defaults to `datetime.date.today().isoformat()`. Call it at the end of `step1_ncbi(genomes, force)` with an `log.info` summary of the counts.

- [ ] **Step 4: Run** the same command — Expected: PASS. Then `uv run pytest tests/test_download_genome_data.py -v` — all pass.

- [ ] **Step 5: Real fetch.** `uv run python -m multiomics_kg.download.download_genome_data --steps 1` (all strains; genome files are cached so only reports are fetched). Expected log: 118 written, 0 missing. Verify: `ls cache/data/*/genomes/*/assembly_report.json | wc -l` → 118; CB0101's file has `report.assembly_info.assembly_status == "previous"`.

- [ ] **Step 6: Commit** the module, tests, `multiomics_kg/utils/assembly_report.py`, and the 118 JSON files (`git add cache/data/*/genomes/*/assembly_report.json` — explicit glob, check `git status` shows only these).

---

### Task 2: Parser + vocabulary + schema

**Files:**
- Modify: `multiomics_kg/utils/assembly_report.py`
- Modify: `config/controlled_vocabularies.yaml` (after `OrganismTaxon.organism_type`)
- Modify: `config/schema_config.yaml` (`organism taxon.properties`, after `taxonomy_note`)
- Modify: `config/gene_annotations_config.yaml` (`ncbi` logical source: `info_types_extra: [assembly_quality]`)
- Test: `tests/test_assembly_report.py` (create), `tests/test_controlled_vocab.py`, `tests/test_data_source_adapter.py`
- Fixtures: `tests/fixtures/assembly_reports/{MED4,MIT9312,WH7803,CB0101}.json` — copied from the Task 1 files, trimmed to `fetched_at` + `report.accession`, `report.current_accession`, `assembly_info.{assembly_level,assembly_status,assembly_name,release_date,suppression_reason}`, `assembly_stats.{total_sequence_length,number_of_contigs,contig_n50,gc_percent}`, `checkm_info`.

**Interfaces:**
- Produces: `assembly_properties(wrapped: dict, accession: str) -> dict[str, str | int | float]` and `load_assembly_report(data_dir: str) -> dict | None` (returns `None` when the file is missing).

Field map (spec table; this is the contract):

| property | source | transform |
|---|---|---|
| `assembly_level` | `assembly_info.assembly_level` | `.lower().replace(" ", "_")`, then `VOCAB.check` |
| `assembly_status` | `assembly_info.assembly_status` | `VOCAB.check` |
| `assembly_suppression_reason` | `assembly_info.suppression_reason` | text after `"Reason:"` stripped of trailing `.`; if no `Reason:` → `"removed by RefSeq staff"`; only when status `suppressed` |
| `assembly_current_accession` | `report.current_accession` | only when `!= accession` |
| `assembly_name` | `assembly_info.assembly_name` | |
| `assembly_release_date` | `assembly_info.release_date` | |
| `assembly_report_date` | `wrapped.fetched_at` | |
| `genome_length` / `contig_count` / `contig_n50` | `assembly_stats.total_sequence_length` / `number_of_contigs` / `contig_n50` | `int()` |
| `gc_percent` | `assembly_stats.gc_percent` | `float()` |
| `checkm_completeness` | `checkm_info.completeness` | `float()` |
| `checkm_contamination` | `checkm_info.contamination` | `float()`; `0.0` if `checkm_info` present and key absent |
| `checkm_marker_set` / `checkm_marker_set_rank` / `checkm_version` | `checkm_info.*` | verbatim |

- [ ] **Step 1: Write failing tests** in `tests/test_assembly_report.py`:
  - `test_med4`: `assembly_level == "complete_genome"`, `assembly_status == "current"`, `checkm_completeness == 99.45`, `checkm_contamination == 0.0`, `checkm_marker_set_rank == "species"`, `"assembly_current_accession" not in props`, `"assembly_suppression_reason" not in props`, `assembly_report_date == "2026-11-18"`.
  - `test_mit9312_int_completeness`: `checkm_completeness == 100.0 and isinstance(..., float)`.
  - `test_wh7803_suppressed_no_checkm`: `assembly_status == "suppressed"`, `assembly_suppression_reason == "missing strain identifier"`, no key starting `checkm_`.
  - `test_cb0101_previous`: `assembly_level == "contig"`, `assembly_status == "previous"`, `assembly_current_accession == "GCF_000179235.2"`.
  - `test_reason_without_clause`: synthetic suppressed report whose reason is `"This record was removed by RefSeq staff. Please contact info@ncbi.nlm.nih.gov for further details."` → `"removed by RefSeq staff"`.
  - `test_unknown_status_raises`: status `"replaced"` → `pytest.raises(ValueError)`.
  - `test_unknown_level_raises`: level `"Plasmid Only"` → `ValueError`.
  - `test_load_missing_returns_none(tmp_path)`.
  - `test_no_empty_strings`: for every fixture, no value is `""` or `None`.
  - In `tests/test_controlled_vocab.py`, `test_assembly_quality_entries`: the 9 ids exist; `assembly_level.values == ["complete_genome", "chromosome", "scaffold", "contig"]`, `assembly_status.values == ["current", "suppressed", "previous"]`, both closed + sparse + with `value_descriptions`; `checkm_completeness` `(value_type, min_value, max_value) == ("float", 0, 100)`; `contig_count.min_value == 1`; `checkm_marker_set_rank.closed is False`. Add `"assembly_level", "assembly_status"` to the `trust_props` set in `test_trust_slice_carries_value_descriptions`.
  - In `tests/test_data_source_adapter.py`: the `ncbi` node's `info_types` contains `"assembly_quality"`.

- [ ] **Step 2: Run** `uv run pytest tests/test_assembly_report.py tests/test_controlled_vocab.py tests/test_data_source_adapter.py -v` — Expected: FAIL.

- [ ] **Step 3: Vocabulary entries** (`config/controlled_vocabularies.yaml`), all `applies_to: OrganismTaxon`, `applies_to_kind: node`, `sparse: true`, researcher-facing descriptions (no build provenance, no `'` or `|` — enforced by existing tests):
  - `assembly_level` (string, closed). Description: NCBI assembly level of the genome the organism's genes come from; for reference_proteome_match organisms it describes the reference assembly, not the sample. Value descriptions: `complete_genome` — "closed assembly, so a missing gene is not an assembly gap; necessary for a confident absence call, still check the annotation"; `chromosome`, `scaffold`, `contig` — one line each, saying gaps can hide genes.
  - `assembly_status` (string, closed). Value descriptions: `current`; `suppressed` — "RefSeq removed the record (see assembly_suppression_reason); the genome and its genes remain valid KG content"; `previous` — "a newer version of this assembly exists (assembly_current_accession)".
  - `checkm_completeness`, `checkm_contamination`, `gc_percent`: `float`, `min_value: 0`, `max_value: 100`, closed false, values `[]`. Contamination description: "0.0 = NCBI omitted the field (zero)". Completeness description: CheckM under-scores streamlined genomes (MED4 99.45, NS01 87.31 are closed genomes) — use assembly_level for absence calls; the reference-proteome caveat.
  - `genome_length`, `contig_count`, `contig_n50`: `int`, `min_value: 1`.
  - `checkm_marker_set_rank`: string, `closed: false`, values `[species, genus, family, order, domain]`; description: coarser rank means a weaker completeness estimate.

- [ ] **Step 4: Schema** — in `config/schema_config.yaml` `organism taxon.properties`, declare the 16 properties with types from the field map (`str` / `int` / `float`), comment `# NCBI Datasets assembly report (step 0); sparse`.

- [ ] **Step 5: DataSource** — `info_types_extra: [assembly_quality]` on the `ncbi` logical source.

- [ ] **Step 6: Implement** `assembly_properties` / `load_assembly_report` per the field map.

- [ ] **Step 7: Run** Step 2's command — Expected: PASS. Then `uv run pytest -m "not slow and not kg" -q` — all pass (hash-pinning tests, if any, updated deliberately).

- [ ] **Step 8: Commit** the six files + fixtures.

---

### Task 3: Organism node carries the properties

**Files:**
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` (`download_data` loads, `_get_organism_node` merges)
- Test: `tests/test_cyanorak_ncbi_adapter.py` (new `class TestOrganismAssemblyProperties`)

**Interfaces:**
- Consumes: `load_assembly_report(data_dir)`, `assembly_properties(wrapped, accession)` from Task 2.
- Produces: `self.assembly_report: dict | None` on `CyanorakNcbi`.

- [ ] **Step 1: Write failing tests** using the existing `adapter` fixture pattern (`temp_data_dir`, MED4, `GCF_000011465.1`):
  - `test_report_properties_on_node`: copy `tests/fixtures/assembly_reports/MED4.json` to `temp_data_dir/assembly_report.json` before `download_data()`; `_get_organism_node()[0][2]` has `assembly_level == "complete_genome"` and `checkm_contamination == 0.0`.
  - `test_missing_report_file`: no file → node built, no key starting `assembly_` or `checkm_`, a warning logged (`caplog`) naming the strain.
  - `test_strings_are_cleaned`: a report whose `assembly_name` contains `'` → value contains `^`.
  - Treatment organisms: `MultiCyanorakNcbi._get_treatment_organism_nodes` output has no `assembly_*` keys (existing fixture CSV).

- [ ] **Step 2: Run** `uv run pytest tests/test_cyanorak_ncbi_adapter.py -k Assembly -v` — Expected: FAIL.

- [ ] **Step 3: Implement**: in `download_data`, `self.assembly_report = load_assembly_report(self.data_dir)` (warn when `None`); in `_get_organism_node`, `properties.update({k: self.clean_text(v) if isinstance(v, str) else v for k, v in assembly_properties(...).items()})` when the report exists.

- [ ] **Step 4: Run** the adapter test file, then `uv run pytest -m "not slow and not kg" -q` — Expected: PASS.

- [ ] **Step 5: Smoke**: `uv run python create_knowledge_graph.py --test --output-dir <scratchpad>/kg-test` and grep the organism CSV header for `assembly_level` and `checkm_completeness`; MED4's row carries `complete_genome`.

- [ ] **Step 6: Commit** adapter + test.

---

### Task 4: KG validity tests, docs, rebuild

**Files:**
- Modify: `tests/kg_validity/test_organism.py`
- Modify: `CLAUDE.md` (OrganismTaxon computed-properties bullet; Data Locations; trust-vocabulary list + "63 nodes" → 65)
- Modify: `CHANGELOG.md` (`[Unreleased]` → `### Added` bullet + `### Highlights` line)
- Modify: `plans/backlog.md` (remove the assembly-quality item; add "migrate CB0101 to GCF_000179235.2" and "review the 10 RefSeq-suppressed assemblies (BS11, WH7803, HNLC1, HNLC2, A18-46.1, WH8103, A18-40, RCC307, RS9915, BOUM118)")

- [ ] **Step 1: Write KG tests** (`@pytest.mark.kg`, `run_query` fixture):
  - `test_assembly_properties_on_genome_organisms`: every node with `organism_type IN ['genome_strain','reference_proteome_match']` has non-null `assembly_level`, `assembly_status`, `assembly_report_date` (count of missing == 0).
  - `test_no_assembly_properties_on_treatment`: `organism_type = 'treatment'` and `assembly_level IS NOT NULL` → 0 rows.
  - `test_assembly_spot_checks`: MED4 `(complete_genome, current, 99.45, 0.0)`; WH7803 `assembly_status = 'suppressed'`; CB0101 `(contig, previous, 'GCF_000179235.2')`.
  - `test_assembly_level_distribution`: 118 nodes carry `assembly_level`; `complete_genome` count ≥ 70.

- [ ] **Step 2: Docs** — CLAUDE.md, CHANGELOG, backlog as listed above.

- [ ] **Step 3: Unit suite** `uv run pytest -m "not slow and not kg" -q` — PASS. Commit tests + docs, push the branch (`git push -u origin feat/organism-assembly-quality`).

- [ ] **Step 4: Rebuild (other clone, user-run)**: pull the branch in `multiomics_biocypher_kg`, `docker compose up -d`. Then `pytest tests/kg_validity/ -v` there — Expected: all pass, including `test_schema_coverage.py` and the four new tests; `test_import_report.py` shows no skipped rows.
