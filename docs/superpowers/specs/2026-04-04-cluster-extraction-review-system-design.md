# Cluster Extraction Deployment & Review System

**Date:** 2026-04-04
**Status:** Approved

## Problem

The cluster extraction pipeline exists (`multiomics_kg/extraction/cluster/`) with a 3-path architecture (visual, semantic, table) but has only been run on 2 of 15 configured `gene_clusters` entries. Of 16 extracted clusters (both Tolonen 2006), only 3 have "pass" verdicts. The remaining 13 entries across 7 papers have no extraction at all — all cluster nodes in the KG use empty placeholder descriptions.

There is no way to:
- See where extraction is failing and why
- Classify failure modes to guide pipeline improvements
- Manually review and approve/edit/reject results before they enter the KG
- Iterate on the pipeline with visibility into what changed

## Goal

Improve cluster description extraction quality, scale to all 15 entries across 8 papers, with a human-in-the-loop review UI as the central diagnostic and approval tool.

## Key Principles

- Empty is always better than wrong — LLM agents downstream trust descriptions
- The review UI is the diagnostic instrument, not just a final gate
- Issue classification drives pipeline improvement — no guessing at fixes
- Full audit trail: every run preserved, every review decision tracked

## Phased Approach

### Phase 1: Infrastructure

1. **Refactor pipeline data organization** — new directory structure with versioned runs, separated pre-processing cache, and per-stage output files (replacing monolithic extraction JSON)
2. **Migrate existing Tolonen 2006 JSONs** — split `cluster_extraction_*.json` into stage1/2/3 files in an initial run directory
3. **Re-run Tolonen 2006** — fresh extraction on current pipeline version, validates new file structure and gives a clean baseline

### Phase 2: Review UI

4. **Build Streamlit review app** — develop against Tolonen 2006 results
5. **Manual review of Tolonen** — classify all 16 clusters, test the full review workflow

### Phase 3: Diagnose & Improve Pipeline

6. **Export and analyze issue report** — derive pipeline fixes from classified failures
7. **Pipeline improvements** — targeted fixes (prompts, source selection, merge logic, table parsing), iterate using review UI

### Phase 4: Scale

8. **Batch-run remaining 13 entries** across 7 papers
9. **Review queue** — work through unreviewed results in UI, iterate as needed

## Data Organization

### Run directory structure

```
{paper_dir}/.extraction_cache/{entry_key}/
  runs/
    2026-04-04T10-30-00/          # timestamped, immutable (stages 1-3)
      stage1_merged.json
      stage2_results.json
      stage3_validation.json
      stage4_review.json          # human review, copied forward on re-run
    2026-04-05T14-22-00/
      ...
  current -> 2026-04-05T14-22-00/ # symlink to active run
{paper_dir}/.extraction_cache/shared/
  pdf_text.json                   # pre-processing cache (shared across entries & runs)
  pages/                          # PDF page images for vision path
  embeddings.npz                  # text embeddings for RAG
  rag_chunks.json                 # RAG retrieved chunks
  tables/                         # parsed supplementary tables
```

Each entry key (e.g., `med4_kmeans_nstarvation`) gets its own `runs/` directory since Stage 1-3 outputs are per-entry. Pre-processing cache is shared at the paper level since PDF text, embeddings, and page images are the same for all entries from the same paper.

### Design rationale

- Each run directory holds the full pipeline trace for one extraction execution
- Stages 1-3 are immutable per run — never modified after creation
- Stage 4 (review) is copied forward when a new run is created
- Pre-processing cache lives outside `runs/`, shared across runs to avoid redundant PDF parsing, embedding, etc.
- Old runs are preserved for diff comparison and audit

### Change detection on re-run

Uses input hashing, not output comparison (LLM outputs vary in wording even with identical inputs):
- Per cluster: hash the Stage 1 merged inputs (quotes, confidence values, source fields)
- If hash unchanged from the run that was reviewed → review status carries forward as-is
- If hash changed → review status marked `stale` (highlighted in UI for re-review)
- Prompt-only changes (same inputs, different synthesis logic) require manual batch re-review since inputs haven't changed

### Migration of existing data

The two existing `cluster_extraction_*.json` files (Tolonen 2006) are split into stage1/2/3 files in an initial run directory. This is a one-time migration. After that, new extractions use the new structure exclusively.

## Review UI

### Technology

Streamlit app at `multiomics_kg/review/cluster_review_app.py`.

### Layout

**Left sidebar:**
- Paper → ClusteringAnalysis → Cluster hierarchy in dropdowns
- Color-coded by review status rollup:
  - Green: all clusters approved (reviewed in current run)
  - Light green: all clusters approved but review carried forward from a previous run
  - Yellow: some stale or partially reviewed
  - Red: any rejected or unreviewed
- Filter by: verdict (pass/warn/fail), review status (approved/stale/unreviewed/rejected)
- Summary stats: total clusters, pass/warn/fail counts, unreviewed count

**Center pane — merge-centric cluster view:**

Each cluster displayed as a card showing the full pipeline trace:

Per field (`functional_description`, `enrichment_category`, `direction`, `temporal_pattern`, etc.):

| Source | Value | Confidence |
|---|---|---|
| table | "Nitrogen transport and assimilation" | very_high |
| visual | "N uptake genes, amt transporters" | high |
| semantic | "Transport-related cluster" | medium |
| **-> synthesis** | **"Nitrogen transport and assimilation genes..."** | — |

Below the merge view:
- Supporting quotes — pooled from all paths, with source label, relevance score, page reference
- Stage 3 verdict badge (pass/warn/fail) + explanation text

**Review controls per cluster:**
- **Review status** (enum): `approve` / `edit` / `reject` / `flag-issue`
- **Issue classification** (multi-select enum): `wrong_cluster`, `hallucinated`, `low_info`, `cross_contamination`, `partial`, `source_missing`
- **Failing stage** (multi-select enum): `stage1_visual`, `stage1_semantic`, `stage1_table`, `stage2_synthesis`, `stage3_validation`, `merge`
- **Free text notes**: specifics, context, anything not captured by enums
- **Editable description fields**: active when status is `edit`

### Source access

- **Tabular files** (XLS/XLSX/CSV/TSV): rendered inline as Streamlit dataframes, searchable/sortable. Multi-sheet XLSX gets a sheet selector tab.
- **PDFs**: "open at page N" links via `file:///path#page=N` (opens in browser)
- **Everything else** (DOCX, HTML, TXT): `xdg-open` fallback
- **Per-quote links**: each supporting quote has a "-> source" link pointing to the originating file/page
- **Keyboard shortcut**: `p` to open main PDF for current paper

### Actions

- **Save reviews** — writes `stage4_review.json` in current run directory
- **Export issue report** — CSV summary across all papers: paper, entry, cluster, verdict, review_status, issues, failing_stages, notes
- **Re-run per cluster**: "Re-synthesize" (stage 2+3 only, uses cached Stage 1) or "Full re-extract" (all stages)
- **Re-run per entry**: with force mode dropdown:
  - Default: uses cached pre-processing, re-runs LLM stages
  - Force LLM: re-runs all LLM calls (stages 1-3), keeps pre-processing cache
  - Force all: clears pre-processing cache, full fresh extraction
- **Diff view**: after re-run, show old vs new extraction side-by-side with field-level diffs (green=improved, red=regressed, gray=unchanged)

### Review persistence format

Stored in `stage4_review.json` per run:

```json
{
  "1": {
    "status": "approve",
    "issues": ["partial"],
    "failing_stages": ["stage1_semantic"],
    "notes": "Cluster 3 genes attributed here instead",
    "edited_fields": {
      "functional_description": "corrected text..."
    },
    "reviewed_at": "2026-04-04T12:30:00",
    "reviewed_in_run": "2026-04-04T10-30-00",
    "input_hash": "abc123..."
  }
}
```

`reviewed_in_run` records which run the review was originally made in. When reviews are copied forward, `reviewed_at` and `reviewed_in_run` are preserved. The UI distinguishes:
- **Current review** (reviewed_in_run == current run): solid green/red badge
- **Carried-forward review** (reviewed_in_run != current run): faded/outlined badge with "reviewed on {date}" tooltip

This makes it immediately visible which clusters you've actually looked at in this run vs inherited from a previous one.

## Pipeline CLI

```bash
# Single entry
uv run python -m multiomics_kg.extraction.cluster.pipeline --entry med4_kmeans_nstarvation

# Single paper (all entries)
uv run python -m multiomics_kg.extraction.cluster.pipeline --paper "Biller 2018"

# All configured entries
uv run python -m multiomics_kg.extraction.cluster.pipeline --all

# Force modes
uv run python -m multiomics_kg.extraction.cluster.pipeline --all --force-llm
uv run python -m multiomics_kg.extraction.cluster.pipeline --all --force-all
```

## KG Graduation

The cluster adapter reads `current/stage2_results.json` + `current/stage4_review.json`. A cluster's extracted descriptions are used only when:
- Review status is `approve` or `edit`, AND
- Review is not `stale` (input hash matches current Stage 1)

Otherwise: empty descriptions (current behavior, safe default).

## Pre-processing Cache

| Operation | Cache file | Invalidation |
|---|---|---|
| PDF text extraction | `pdf_text.json` | Manual / `--force-all` |
| PDF page images (vision) | `pages/` | Manual / `--force-all` |
| Text embeddings (RAG) | `embeddings.npz` | If PDF text changes |
| RAG retrieved chunks | `rag_chunks.json` | If embeddings change |
| Parsed supplementary tables | `tables/` | If source files change |

Stage 1 results are cached per run in `stage1_merged.json`. Pre-processing caches are shared across runs.

## Scope Boundaries

**In scope:**
- Pipeline data organization refactor
- Streamlit review UI
- Pipeline CLI with force modes
- Per-cluster and per-entry re-run from UI
- Issue classification and export
- Diff view for comparing runs
- Adapter update to read new file structure + respect review status

**Out of scope:**
- New extraction paths (beyond existing visual/semantic/table)
- Changes to the extraction prompts (that's Phase 3, driven by review data)
- Changes to cluster adapter node/edge schema (already stable)
- Automated CI/pipeline triggers
