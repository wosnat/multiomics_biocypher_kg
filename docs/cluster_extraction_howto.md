# Cluster Extraction — How-To Guide

## Prerequisites

- OpenAI API key in `.env` (`OPENAI_API_KEY=...`)
- `uv sync` to install dependencies (includes `streamlit`, `openpyxl`)
- Papers with `type: gene_clusters` entries in their `paperconfig.yaml`

## 1. Review UI

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py
```

Opens in browser. Sidebar shows papers with gene_clusters entries. Select a paper and entry to review extraction results.

**Color coding:**
- Green: all clusters approved (current run)
- Yellow (light green): approved but carried forward from previous run
- Orange: some stale reviews
- Red: unreviewed or rejected

**Workflow:** Select cluster > review merge view (3 columns: table/visual/semantic) > check synthesis result > set review status (approve/edit/reject/flag-issue) > save.

## 2. Run Extraction Pipeline

### Full run (all stages)
```bash
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

### Single entry
```bash
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/.../paperconfig.yaml" --table-key med4_kmeans_nstarvation
```

### Re-run from a specific stage (uses cached earlier stages)
```bash
# Re-run synthesis + validation only (stages 2-3)
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/.../paperconfig.yaml" --from-stage 2

# Re-run validation only (stage 3)
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/.../paperconfig.yaml" --from-stage 3
```

### Force modes
```bash
# Re-run all LLM calls (stages 1-3), keep pre-processing cache
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/.../paperconfig.yaml" --force-llm

# Clear all caches and re-extract from scratch
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/.../paperconfig.yaml" --force-all
```

You can also trigger re-runs from the review UI sidebar.

## 3. Data Layout

```
{paper_dir}/
  paperconfig.yaml
  .extraction_cache/
    shared/                          # pre-processing cache (shared across entries & runs)
      pdf_text.json                  # extracted PDF text
      chunks.json                    # RAG text chunks
      embeddings.npy                 # chunk embeddings (OpenAI)
      pdf_content_parts.json         # base64-encoded PDF pages for vision API
    {entry_key}/                     # e.g., med4_kmeans_nstarvation
      runs/
        2026-04-04T18-36-35/         # timestamped run directory
          stage1_merged.json         # per-path extraction results merged
          stage2_results.json        # synthesized descriptions
          stage3_validation.json     # per-cluster validation verdicts
          stage4_review.json         # human review decisions
          metadata.json              # run metadata
          report.md                  # markdown summary
      current -> runs/2026-04-...    # symlink to latest completed run
```

- `current` symlink only updated after all stages complete (`finalize_run`)
- Review decisions (`stage4_review.json`) copied forward on re-runs; marked `stale` if inputs changed
- Shared cache persists across runs — delete to force re-extraction

## 4. Dry-Run Prompts (review without API calls)

```bash
# All stages, all clusters, all entries
uv run python scripts/dry_run_prompts.py "data/.../paperconfig.yaml"

# Specific entry + cluster + stage
uv run python scripts/dry_run_prompts.py "data/.../paperconfig.yaml" \
    --entry med4_kmeans_nstarvation --cluster 7 --stage synthesis

# Save to files for review
uv run python scripts/dry_run_prompts.py "data/.../paperconfig.yaml" \
    --output-dir /tmp/dry_run
```

## 5. RAG Query Experiments

Test what paper text the RAG retrieves for different queries, without calling extraction LLMs.

```bash
# Interactive mode (queries cached after first embedding)
uv run python scripts/rag_experiment.py "data/.../paper_dir" --interactive

# Test specific queries
uv run python scripts/rag_experiment.py "data/.../paper_dir" \
    --batch "MED4 cluster 1" "MED4 cluster 7 Translation downregulated" \
    --top-k 5
```

Caches chunks + embeddings to `.extraction_cache/shared/` — subsequent queries are instant.

## 6. Pipeline Stages

| Stage | What it does | Model | Input |
|---|---|---|---|
| 1a: Table | Parse cluster CSV + supplementary XLS for enrichment data | None (pandas) | CSV/XLS files |
| 1b: Visual | Extract cluster info from PDF figures/legends | gpt-4o | PDF pages (per cluster) |
| 1c: Semantic | RAG retrieve paper text, extract cluster info | gpt-5-nano | Embedded text chunks (per cluster) |
| 1d: Merge | Combine table/visual/semantic by confidence | None | Stage 1a-1c outputs |
| 2: Synthesis | Write descriptions from merged data | gpt-5-nano | Merged fields (per cluster) |
| 3: Validation | Verify descriptions against paper | gpt-4o | PDF pages + descriptions (per cluster) |

## 7. Export Issue Report

From the review UI sidebar: click "Export Issue Report" to download a CSV of all non-approved clusters with their issues, failing stages, and notes.

## 8. Adding New Papers

1. Add `type: gene_clusters` entry to the paper's `paperconfig.yaml` (see CLAUDE.md for format)
2. Run extraction: `uv run python -m multiomics_kg.extraction.cluster.pipeline "data/.../paperconfig.yaml"`
3. Review in UI
4. Approved clusters automatically feed into the KG via `cluster_adapter.py`

## 9. Clearing Cache

```bash
# Clear shared cache for a paper (forces re-extraction of PDF text + embeddings)
rm -rf "data/.../paper_dir/.extraction_cache/shared/"

# Clear all runs for an entry (start fresh)
rm -rf "data/.../paper_dir/.extraction_cache/entry_key/"

# Clear everything for a paper
rm -rf "data/.../paper_dir/.extraction_cache/"
```

## 10. Known Issues

See [extraction_issues_log.md](extraction_issues_log.md) for diagnosed issues and their status.
