# Extraction Pipeline Cleanup + Quality Verification

**Date:** 2026-04-05
**Status:** Draft
**Scope:** Delete old pipeline, clean format, quality iteration loop

## Summary

Clean up the cluster extraction pipeline after the ad-hoc redesign. Delete ~1,500 lines of dead old pipeline code (visual/semantic/merge/synthesis/validation stages, RAG, RunManager). Promote the ad-hoc extraction script into a proper module at `multiomics_kg/extraction/cluster/extract.py`. Replace the hacky JSON format (fake stage3_validation verdicts) with a clean `clusters` dict — simplify the adapter to a single code path with no verdict gating. Build a diff-friendly quality report (`data/cluster_extraction_report.md`) with programmatic checks so we can iterate: extract → report → diff → tweak prompt → re-extract.

## Context

We redesigned the cluster extraction pipeline (see `2026-04-05-extraction-pipeline-redesign.md`) and populated all 115 clusters via an ad-hoc script (`scripts/run_extraction.py`). The KG now has 82 described clusters and 33 correctly-empty ones.

What remains:
- ~1,500 lines of dead old pipeline code still in the repo
- The ad-hoc script is in `scripts/` with sys.path hacks
- The JSON format fakes `stage3_validation` verdicts to pass the adapter's gate
- No quality verification has been done beyond the pilot's 3 ground-truth clusters
- No iteration workflow for improving extraction quality

## Deliverables

### 1. Delete old pipeline code + junk

**Delete:**
- `multiomics_kg/extraction/cluster/visual.py`
- `multiomics_kg/extraction/cluster/semantic.py`
- `multiomics_kg/extraction/cluster/merge.py`
- `multiomics_kg/extraction/cluster/synthesis.py`
- `multiomics_kg/extraction/cluster/validation.py`
- `multiomics_kg/extraction/cluster/pipeline.py`
- `multiomics_kg/extraction/cluster/migrate.py`
- `multiomics_kg/extraction/cluster/prompts.py`
- `multiomics_kg/extraction/cluster/run_manager.py`
- `multiomics_kg/extraction/rag.py`
- `multiomics_kg/review/cluster_review_app.py` (imports RunManager; rebuild later against new format)
- `multiomics_kg/review/review_components.py` (imports RunManager)
- `multiomics_kg/review/review_data.py` (imports RunManager)
- `scripts/run_extraction.py` (moved to `extract.py`)
- `scripts/pilot_extraction.py`
- `scripts/dry_run_prompts.py`
- `scripts/rag_experiment.py` (imports `rag.py`)
- `tests/test_cluster_extraction_migration.py` (tests `migrate.py` + `RunManager`)
- `tests/test_review_data.py` (tests review UI with `RunManager`)
- `tests/test_extraction_merge.py` (tests `merge.py`)
- `tests/test_run_manager.py` (tests `RunManager`)
- All `data/**/.extraction_cache/` directories
- All `data/**/pilot_comparison/` directories

**Keep:**
- `multiomics_kg/extraction/cluster/table.py` — Stage 0 CSV parsing, still used
- `multiomics_kg/extraction/cluster/__init__.py`
- `multiomics_kg/extraction/pdf_utils.py` — emergency base64 fallback
- `data/**/pilot_extraction_results.json` — ground truth (Tolonen clusters 1, 6, 7)

**Move** extraction results into a subfolder per paper:
- `data/<paper_dir>/cluster_extraction_*.json` → `data/<paper_dir>/cluster_extractions/`
- `data/<paper_dir>/cluster_extraction_*.md` → `data/<paper_dir>/cluster_extractions/`

This keeps paper directories clean (they already have CSVs, PDFs, paperconfig.yaml).

### 2. Clean extraction module

Move `scripts/run_extraction.py` → `multiomics_kg/extraction/cluster/extract.py` as the canonical extraction module.

**Cleanup:**
- Remove sys.path hack (proper package import)
- Pydantic schemas at module top
- Prompt template as module constant
- Functions organized into layers (see section 5 for full breakdown)

**CLI modes** (via `python -m multiomics_kg.extraction.cluster.extract`):
- Default: run extraction for all entries (or filtered by `--paper`/`--entry`)
- `--report`: generate `data/cluster_extraction_report.md` from existing JSONs (no API calls)
- `--verify`: run programmatic quality checks, append warnings to report
- `--dry-run`: show what would run
- `--force`: overwrite existing extractions
- `--flex`: use Flex processing (50% cheaper)
- `--model MODEL`: override default model

### 3. Clean JSON format + adapter

**New JSON format** — stored at `{paper_dir}/cluster_extractions/{entry_key}.json`:

```json
{
  "metadata": {
    "paper": "Tolonen 2006",
    "doi": "10.1038/msb4100087",
    "organism": "Prochlorococcus MIT9313",
    "entry_key": "mit9313_kmeans_nstarvation",
    "model": "gpt-4.1-mini",
    "extracted_at": "2026-04-05T13:25:22",
    "input_tokens": 62080,
    "output_tokens": 3141
  },
  "clusters": {
    "1": {
      "id": "mit9313_up_transport_binding",
      "name": "MIT9313 cluster 1 (up, transport and binding)",
      "functional_description": "...",
      "behavioral_description": "...",
      "peak_time_hours": 6.0,
      "period_hours": null,
      "direction": "up",
      "enrichment_category": "transport and binding",
      "enrichment_pvalue": 0.04,
      "enrichment_significant": true,
      "self_assessment": "high",
      "assessment_notes": "...",
      "confidence_notes": "...",
      "supporting_quotes": [{"quote": "...", "location": "..."}]
    }
  }
}
```

Changes from current hacky format:
- `clusters` dict replaces `stage2_results`
- No `stage3_validation` — no verdict gating
- No `stage4_review` — human edits go directly into clusters dict
- Self-assessment is informational, doesn't gate inclusion

**New utility** `multiomics_kg/extraction/cluster/extraction_utils.py`:

File I/O (knows about file locations and JSON structure):
- `load_extraction(paper_dir, entry_key)` → `dict[str, dict]` — reads `{paper_dir}/cluster_extractions/{entry_key}.json`, returns `data["clusters"]` dict. Returns `{}` if file missing. Logs warning if format is wrong.
- `save_extraction(paper_dir, entry_key, metadata, clusters)` → writes `{paper_dir}/cluster_extractions/{entry_key}.json` (creates dir if needed) + `{entry_key}.md` summary.
- `get_cluster_data(clusters, key)` → `dict` — returns `clusters.get(str(key), {})`.
- `list_extraction_files(paper_dir)` → list of entry_keys that have extraction JSONs.

Data loading (reads paperconfigs and CSVs):
- `find_all_entries()` → list of `(paper_dir, entry_key, table_config, pub_config)` — discovers all `gene_clusters` entries across all paperconfigs.
- `load_cluster_summaries(table_config)` → `dict[str, dict]` — reads cluster CSV, returns `{cluster_key: {gene_count, sample_genes}}`.

Matching (pure logic):
- `match_cluster_keys(parsed_clusters, expected_keys)` → `(matched: dict[str, dict], unmatched: list)` — maps model output names/ids to actual cluster keys. Tries name regex, id suffix, case-insensitive match, positional fallback.

Used by the adapter, `extract.py`, and the report/verify tooling.

**Adapter changes** in `cluster_adapter.py`:
- Replace `_load_extraction_json()` and `_get_extraction_cluster_data()` with calls to `extraction_utils.load_extraction()` and `extraction_utils.get_cluster_data()`.
- Delete: RunManager import, stage key mapping, verdict gating logic, all legacy code paths.

### 4. Quality verification

**Review report:** `data/cluster_extraction_report.md`

Diff-friendly design:
- Fixed sort order: papers alphabetical, clusters by key within analysis
- One cluster per section — diffs are localized
- Header: `### Cluster {key} | {direction} | {self_assessment}`
- Body: name, enrichment, functional, behavioral, notes
- No timestamps or token counts — only content that matters for quality
- `## Warnings` section at bottom with programmatic check results

```markdown
# Cluster Extraction Report

## Alonso-Saez 2023 / mit9301_softclusters_thermal_acclimation

### Cluster A | up | high
**Name:** MIT9301 cluster A (up, core metabolism)
**Enrichment:** core metabolism (p=0.01, sig=true)
**Functional:** This cluster contains genes involved in...
**Behavioral:** Genes show increased expression at higher temperatures...
**Notes:** Strong enrichment supported by Figure 2.

### Cluster B | down | medium
...

## Warnings

- [Tolonen 2006 / med4 / cluster 3] duplicate id: med4_up_regulation (also used by cluster 2)
- [Zinser 2009 / cluster 10] locus tag in description: "PMM0042"
- [Bernstein 2017 / bp1_oxygen / cluster 0] enrichment p-value 0.05 not found in Stage 0 table data
```

**Programmatic checks:**
1. Duplicate `id` fields within an analysis
2. Locus tags in descriptions (regex: `PMM\d+|PMT\d+|P9301_\d+|tll\d+|SY28_\d+` etc.)
3. Enrichment p-value cross-check against Stage 0 table data (where supplementary parsers exist)
4. Empty direction for clusters with non-empty descriptions
5. KG node data matches extraction JSON (query KG GeneCluster nodes, exact string match on `functional_description` — only run when `--verify` and KG is reachable, skip otherwise)

**Iteration workflow:**
1. `python -m multiomics_kg.extraction.cluster.extract --report --verify` — generate report from existing data
2. Read report, identify issues
3. Tweak prompt or fix specific clusters
4. `python -m multiomics_kg.extraction.cluster.extract --entry X --force` — re-extract
5. `python -m multiomics_kg.extraction.cluster.extract --report --verify` — regenerate
6. `git diff data/cluster_extraction_report.md` — see what changed
7. Repeat until satisfied

### 5. Module design for `extract.py`

The current script is a monolith. Decompose into focused functions with clear responsibilities:

```
multiomics_kg/extraction/cluster/extract.py
```

**Pydantic schemas** (at module top):
- `ClusterExtraction`, `AnalysisExtraction`, `SupportingQuote` — same as pilot

**Prompt layer** (pure string construction, easy to iterate):
- `DEVELOPER_MSG_TEMPLATE` — single constant
- `build_context_block(table_config)` → str
- `format_cluster_summaries(clusters)` → str

**LLM layer** (API calls, isolated for testing/mocking):
- `upload_pdf(client, pdf_path)` → file_id
- `extract_analysis(client, file_id, table_config, cluster_summaries, model, flex)` → (AnalysisExtraction, usage)
- `extract_paper(client, file_id, tables_and_summaries, model, flex)` → (AnalysisExtraction, usage)
  - Used for multi-organism papers (Tolonen) where per-analysis confuses the model

**Report layer** (reads JSONs via `extraction_utils`, no LLM):
- `generate_report(entries)` → markdown string
- `verify_quality(entries)` → list of warnings

**CLI** (`__main__` block):
- Parses args, uses `extraction_utils` for discovery/IO, calls LLM layer, handles paper-grouping loop

This decomposition means:
- Prompt changes = edit one constant in `extract.py`, re-extract
- Format/path changes = edit `extraction_utils.py` only
- New papers = just run the CLI
- Testing = mock the LLM layer, test everything else via `extraction_utils`

### 6. Tests

**extraction_utils tests** (unit, no LLM):
- `test_load_extraction`: write clean-format JSON to tmp dir, load via `load_extraction()`, verify cluster dict
- `test_load_wrong_format`: write old format with `stage2_results` key, verify returns `{}`
- `test_load_missing_file`: verify `{}` returned, no crash
- `test_save_load_roundtrip`: save via `save_extraction()`, load back, verify identical
- `test_list_extraction_files`: create 2 JSONs in `cluster_extractions/`, verify both entry_keys returned
- `test_match_cluster_keys_numeric`: keys "1"-"7", model outputs "cluster 1" through "cluster 7"
- `test_match_cluster_keys_alpha`: keys "HEG", "LEG", etc., model outputs "cluster HEG"
- `test_match_cluster_keys_composite`: keys "coculture_LD+coculture_darkness", model outputs numbered names → positional fallback
- `test_match_cluster_keys_unmatched`: model outputs garbage → returned in unmatched list
- `test_load_cluster_summaries`: load from real Tolonen CSV, verify 7 clusters with correct gene counts

**Adapter tests** (unit, no LLM):
- `test_cluster_properties_on_node`: verify adapter emits correct node properties (`functional_description`, `name`, `peak_time_hours`, etc.) from extraction data

**Extract module tests** (unit, no LLM):
- `test_build_context_block`: verify prompt includes organism, method, treatment
- `test_format_cluster_summaries`: verify output has all cluster keys
- `test_generate_report_stable_order`: generate report twice → identical output (no ordering jitter)

**Verification tests** (unit, no LLM):
- `test_detect_duplicate_ids`: two clusters with same id → warning
- `test_detect_locus_tags`: description containing "PMM0042" → warning
- `test_detect_empty_direction`: non-empty description with empty direction → warning

**Integration test** (requires existing extraction data, no LLM):
- `test_adapter_reads_all_entries`: load all 15 JSONs via adapter, verify 82 clusters have descriptions

### 7. Code review checklist

**Also update** (not delete):
- `tests/test_cluster_adapter.py` — remove RunManager-based tests, add tests against new format via `extraction_utils`

After implementation, verify:
- [ ] No imports of deleted modules anywhere in codebase
- [ ] `grep -r "RunManager\|stage4_review\|stage2_results\|stage3_validation\|stage1_merged" multiomics_kg/` returns nothing
- [ ] All 15 extraction JSONs are in new clean format under `cluster_extractions/`
- [ ] Report generates without errors
- [ ] Verification checks pass (or warnings are documented)
- [ ] `pytest -m "not slow and not kg"` passes
- [ ] KG build produces same 82 populated + 33 empty clusters

## Non-goals

- Review UI rebuild — the Streamlit app is deleted (depends on RunManager). Can be rebuilt later against `extraction_utils` if needed.
- Automated prompt optimization — manual iteration is fine for 15 entries
- RunManager versioning — one JSON per entry is sufficient for a one-time-per-paper task

## File structure after cleanup

```
multiomics_kg/extraction/
  cluster/
    __init__.py
    extract.py            # LLM extraction + CLI (was scripts/run_extraction.py)
    extraction_utils.py   # file I/O: load/save/list extraction JSONs (used by adapter + extract + report)
    table.py              # Stage 0 CSV parsing
  pdf_utils.py            # base64 fallback

multiomics_kg/adapters/
  cluster_adapter.py    # simplified: reads clusters dict, no verdict gating

data/<paper_dir>/cluster_extractions/
  {entry_key}.json                      # clean format
  {entry_key}.md                        # per-entry human-readable summary

data/
  cluster_extraction_report.md          # full report for review + diff
```
