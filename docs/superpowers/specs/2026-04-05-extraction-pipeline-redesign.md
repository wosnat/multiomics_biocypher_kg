# Cluster Extraction Pipeline Redesign

**Date:** 2026-04-05
**Status:** Draft (updated with pilot results)
**Scope:** Extraction pipeline refactor + Review UI improvements

## Problem

The current 4-stage extraction pipeline (visual → semantic → merge → synthesis → validation) produces low-quality results and is brittle:

- **Direction errors**: Cluster 6 (MIT9313 photosynthesis) extracted as "up" when the paper clearly shows "down" — the merge step lost the visual signal
- **JSON parse failures**: 10/16 Tolonen clusters fail validation because gpt-4o returns prose instead of JSON (Issue 6)
- **Thin behavioral descriptions**: The semantic path can't read figures, and the visual path's output gets diluted through merge + synthesis
- **Slow iteration**: 4 LLM calls per cluster, 4 prompts to tune — hard to diagnose which stage caused an error

## Pilot Results

A comparison pilot (`scripts/pilot_extraction.py`) tested 5 factors on Tolonen 2006 MIT9313 (7 clusters, ground truth for clusters 1, 6, 7). Full results at `data/Prochlorococcus/papers_and_supp/tolonen 2006/pilot_comparison/`.

### Factor 1: Granularity (per-cluster vs per-analysis vs per-paper)

All three achieved 100% accuracy on ground truth clusters. Cost differences are dramatic:

| Granularity | API calls | Input tokens | Projected cost (115 clusters, Flex) |
|---|---|---|---|
| per-cluster | 3 (for 3 clusters) | 185K | $1.48 |
| **per-analysis** | **1 (for 7 clusters)** | **62K** | **$0.23** |
| per-paper | 1 (for 16 clusters) | 62K | $0.11 |

Per-cluster sends the PDF once per cluster (3x for 3 clusters, ~115x at scale). Per-analysis sends it once for all clusters in an analysis. **3x fewer input tokens, same quality.** Per-paper is cheapest but mixes organisms in one output — per-analysis is the sweet spot.

**Decision: per-analysis extraction.** One call per `gene_clusters` entry (15 calls total for 115 clusters).

### Factor 2: RAG (with vs without text passages)

| Condition | Input tokens | Quality |
|---|---|---|
| per-cluster, no RAG | 184,698 | 100% |
| per-cluster, RAG | 187,386 | 100% |
| per-analysis, no RAG | 61,773 | 100% |
| per-analysis, RAG | 63,035 | 100% |

RAG adds ~1-2K tokens and produces identical quality. The model reads the full PDF via `input_file` (text + page images) — RAG passages are redundant.

**Decision: drop RAG.** Delete chunking, embedding, and retrieval infrastructure.

### Factor 3: LLM-as-judge (separate validation call)

| Condition | Calls | Input tokens | Findings |
|---|---|---|---|
| per-analysis, no judge | 1 | 62K | — |
| per-analysis, + judge | 8 | 492K | all pass, zero issues found |

The judge used 7 extra calls (430K tokens) re-reading the same PDF and found nothing wrong. Self-assessment fields (`confidence_notes`, `supporting_quotes`, `self_assessment`) in the extraction schema already provide validation signal.

**Decision: drop separate validation stage.** Use self-assessment fields + human review UI.

### Factor 4: Model (gpt-4.1 vs gpt-4.1-mini)

gpt-4.1 could not be tested — single request (~39K tokens for this PDF) exceeds the account's 30K TPM limit. gpt-4.1-mini handled all conditions successfully with 100% accuracy.

**Decision: default to gpt-4.1-mini.** Configurable via env var. Upgrade to gpt-4.1 when TPM limit is raised.

### Factor 5: PDF method (input_file vs base64)

Base64 via Chat Completions API failed (json_object mode requires "json" in message text). The Responses API `input_file` approach worked cleanly.

**Decision: use Responses API `input_file` exclusively.** Keep base64 `pdf_content_parts.json` as emergency fallback only.

### Cluster 6 Litmus Test

The key quality signal — cluster 6 direction must be "down" (photosynthesis repressed under N-starvation). The old pipeline got this wrong. **Every pilot condition got it right**, confirming that single-call extraction with the full PDF resolves the direction error.

## Goals

1. Higher extraction quality — especially behavioral descriptions from figures
2. Zero parse failures — strict structured outputs via Pydantic
3. Fast iteration loop — extract → review → fix prompt → re-run → compare
4. Scale to 115 clusters across 15 analyses in 8 papers
5. Cost under $0.50 per full run (with Flex processing)

## Non-Goals

- Polished production UI — this is a developer tool
- Automated prompt optimization — manual iteration guided by review insights
- Fine-tuning — not enough training data, prompt-based approach works
- RAG infrastructure — full PDF via input_file is sufficient

---

## Part 1: Extraction Pipeline

### Architecture: 2 Stages

| Stage | What | LLM? | Model | API | Output |
|---|---|---|---|---|---|
| **0: Table** | CSV parsing + enrichment parsers | No | — | — | `stage0_table.json` |
| **1: Extract** | Full PDF → all clusters per analysis | Yes | gpt-4.1-mini | Responses API | `stage1_extraction.json` |

Review data: `stage2_review.json` (human review, carried forward across runs).

No separate validation stage — self-assessment fields in the extraction schema + human review UI provide the quality gate.

### Stage 0: Table (no LLM)

Pure code. Reads from paperconfig + CSV + supplementary files.

**Input:** paperconfig `gene_clusters` entry + CSV file
**Output per cluster:**
- `gene_count`: int
- `genes`: list of gene IDs from CSV
- `enrichment_category`: str (from supplementary parser, if available)
- `enrichment_pvalue`: float
- `enrichment_significant`: bool
- `all_enrichments`: list of {category, pvalue, significant}

**Enrichment parsers:** Currently only `parse_tolonen_enrichment_xls`. Keep as-is — most papers won't have machine-readable enrichment data. The parser is invoked when matching XLS files exist in the paper directory; no config needed.

### Stage 1: Extract (one call per analysis)

**One Responses API call per `gene_clusters` entry** — extracts all clusters in the analysis simultaneously. The model sees the full paper + all cluster summaries, enabling cross-cluster disambiguation (e.g., "cluster 6 is down for photosynthesis while cluster 1 is up for transport").

**Input (all in one call):**
1. **Developer message**: Extraction instructions + paperconfig context + all cluster summaries
2. **PDF file**: Uploaded via Files API with `purpose="user_data"` (uploaded once per paper, reused across analyses)

No RAG passages — the PDF `input_file` provides full text + page images natively.

**Paperconfig context block** (injected into developer message):
```
Analysis: {name}
Organism: {organism}
Clustering: {cluster_method}, {n_clusters} clusters
Type: {cluster_type}
Treatment: {treatment}
Time points: {time_points}  # if available
Key figures: {figure_hint}  # if available
Experimental context: {experimental_context}
Omics type: {omics_type}
```

`cluster_type` drives behavioral expectations:
- `response_pattern` → `peak_time_hours` and `period_hours` should be null
- `diel` → expect periodic behavior, `peak_time_hours` relevant
- `periodic` → both `peak_time_hours` and `period_hours` relevant

**Output schema** (strict structured outputs via Pydantic):

```python
from pydantic import BaseModel
from typing import Literal, Optional

class SupportingQuote(BaseModel):
    quote: str
    location: str

class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    behavioral_description: str
    peak_time_hours: Optional[float]
    period_hours: Optional[float]
    direction: Literal["up", "down", "mixed", "not_described"]
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str

class AnalysisExtraction(BaseModel):
    clusters: list[ClusterExtraction]
```

Use `client.responses.parse(..., text_format=AnalysisExtraction)` to get a parsed object. The SDK auto-generates the correct JSON Schema (using `anyOf` for nullable types) and guarantees schema adherence.

**Field split:**
- **Node fields** (written to KG): `id`, `name`, `functional_description`, `behavioral_description`, `peak_time_hours`, `period_hours`
- **Review fields** (stored for debugging/review, not on node): `direction`, `enrichment_category`, `enrichment_pvalue`, `enrichment_significant`, `confidence_notes`, `supporting_quotes`, `self_assessment`, `assessment_notes`

**API configuration:**
- API: Responses API (`client.responses.parse()`)
- Model: `gpt-4.1-mini` (configurable via `CLUSTER_EXTRACTION_MODEL` env var; upgrade to `gpt-4.1` when TPM limits allow)
- `text_format`: `AnalysisExtraction` (Pydantic model)
- `temperature`: 0
- Responses stored by default (available in OpenAI dashboard for debugging)
- `service_tier`: `"flex"` when running in bulk mode (`--flex` flag) for 50% token discount
- `prompt_cache_key`: `"paper_{doi_short}"` to improve cache hits across re-runs of the same paper

**Prompt design principles:**
- Emphasize reading FIGURES for behavioral descriptions — "Figure panels showing heatmaps, time-course plots are the PRIMARY source for direction and timing" (the model receives PDF pages as images automatically via `input_file`)
- Provide cluster_type context so the model knows whether to expect periodic behavior
- "not described in paper" sentinel for unknown fields — never guess
- No locus tags in descriptions — only gene names from paper text
- No treatment conditions in descriptions — those live on the analysis node
- Max 3-5 named genes per description
- Include all cluster summaries so the model can disambiguate across clusters

### PDF Handling: Files API + Responses API `input_file`

Upload each paper's PDF(s) once at the start of a pipeline run:
1. Upload via Files API: `client.files.create(file=open("paper.pdf", "rb"), purpose="user_data")` → returns file object with `id`
2. Store `file.id` in run metadata
3. Reference in all analysis calls as an `input_file` content part:
   ```python
   {"type": "input_file", "file_id": file.id}
   ```

Benefits:
- Upload once, reference in multiple analysis calls for the same paper
- PDF pages automatically extracted as text + images by the Responses API
- Enables prompt caching (identical PDF prefix across analysis calls)

**Fallback to base64:** If Files API upload fails, use inline base64:
```python
{"type": "input_file", "filename": "paper.pdf", "file_data": f"data:application/pdf;base64,{base64_string}"}
```

**File size limit:** Each file must be under 50 MB; combined limit per request is 50 MB.

### Caching & Run Management

**Keep current RunManager pattern** (simplified stage numbering):
```
paper_dir/.extraction_cache/{entry_key}/
  runs/{timestamp}/
    stage0_table.json
    stage1_extraction.json
    stage2_review.json
    metadata.json
    report.md
  current -> runs/{latest}
shared/
  pdf_content_parts.json    # emergency base64 fallback only
```

**Removed from shared/:** `pdf_text.json`, `chunks.json`, `embeddings.npy` (RAG infrastructure no longer needed).

**Review forwarding:** When re-running, copy `stage2_review.json` from previous run. If `input_hash` matches → review carries forward. If hash differs → mark as "stale".

**Backward-compatible output:** Still write `cluster_extraction_{entry_key}.json` and `.md` for the adapter.

### Paperconfig Extensions

New optional fields on `gene_clusters` entries:

```yaml
mit9313_kmeans_nstarvation:
  type: gene_clusters
  # ... existing fields ...
  figure_hint: "Figure 3, Figure 5"       # which figures show cluster results
  time_points: [0, 3, 6, 12, 24, 48]      # explicit time points in hours
```

Both are optional. `figure_hint` is injected into the extraction prompt to focus the model. `time_points` helps the model interpret figure axes and validate behavioral descriptions.

### Model Configuration

| Env var | Default | Used by | Notes |
|---|---|---|---|
| `CLUSTER_EXTRACTION_MODEL` | `gpt-4.1-mini` | Stage 1 | Vision + structured output; upgrade to `gpt-4.1` when TPM allows |

All LLM calls use the **Responses API** (`client.responses.parse()` with Pydantic models). The SDK manages strict JSON schema generation and response parsing.

### Flex Processing for Bulk Runs

Two execution modes:

```python
# Iteration mode (default): sync, immediate results
response = client.responses.parse(model="gpt-4.1-mini", ...)

# Bulk mode (--flex flag): 50% token discount, variable latency
response = client.responses.parse(model="gpt-4.1-mini", service_tier="flex", ...)
```

Flex processing runs through the same Responses API but at 50% lower cost. Combined with `prompt_cache_key="paper_{doi}"`, re-runs after prompt tweaks get automatic cache hits on the PDF prefix.

**Full-run cost projection:** ~$0.23 at normal rate, ~$0.11 with Flex (based on pilot data).

---

## Part 2: Review UI

### Layout

**Left sidebar:**
- Paper list (only papers with extraction data)
- Each paper shows analysis entries underneath
- Entries colored by worst-case cluster status: green (all approved) / yellow (has issues) / red (has rejections) / gray (unreviewed)
- Click analysis → loads clusters in center pane

**Center pane:** Cluster list for selected analysis
- Collapsible rows, one per cluster
- Collapsed: `[6] mit9313_down_photosynthesis — approved (high)` (colored by status + self-assessment)
- Expanded: full detail view (see below)
- Sort by: cluster number (default) | status | direction
- Only show analyses that have extraction data in `.extraction_cache/`

**Expanded cluster detail:**
- **Extraction result**: id, name, functional_description, behavioral_description, direction, enrichment
- **Self-assessment**: high/medium/low confidence + assessment_notes
- **Confidence notes**: how the model arrived at its answer
- **Supporting quotes**: with paper locations
- **Review controls**: approve / edit / reject + issue tags + notes
- **History** (collapsed): previous run's result for this cluster, diff highlighting what changed

### Issue Taxonomy

When flagging an issue, tag with category:

| Tag | Meaning |
|---|---|
| `wrong_direction` | Up/down/mixed is incorrect |
| `hallucinated_genes` | Gene names not in the paper |
| `hallucinated_function` | Functional description not supported |
| `missing_behavioral` | Paper describes temporal pattern but extraction missed it |
| `vague_description` | Correct but too vague to be useful |
| `wrong_enrichment` | Enrichment category or p-value incorrect |
| `cross_contamination` | Info from a different cluster leaked in |
| `wrong_cluster_type` | Periodic vs response_pattern mismatch |

These tags enable **pattern detection**: "8 clusters across 3 papers have `wrong_direction`" → that's a prompt issue worth fixing.

### Cross-Run Comparison

- Per-cluster: side-by-side diff of current vs previous extraction (text diff of descriptions)
- Per-analysis summary: "re-run: 4 improved, 1 regressed, 2 unchanged"
- Status history: show status progression across runs (e.g., rejected → approved)

### Batch Operations

- "Approve all high-confidence" button (approves clusters where self_assessment=high and no issues flagged)
- "Re-run analysis" button (re-runs extraction for the selected analysis)
- Filter + select multiple → bulk approve/reject

### Edit Workflow

- Inline editing of `functional_description` and `behavioral_description`
- Pre-filled from extraction result — user corrects what's wrong
- Edited fields stored in `stage2_review.json` as `edited_fields`
- Adapter uses edited values when `status=edit`

---

## Part 3: File Structure Changes

### Files to delete
- `multiomics_kg/extraction/cluster/visual.py` — replaced by single-call extraction
- `multiomics_kg/extraction/cluster/semantic.py` — RAG retrieval no longer needed
- `multiomics_kg/extraction/cluster/synthesis.py` — replaced by single-call extraction
- `multiomics_kg/extraction/cluster/merge.py` — no longer needed
- `multiomics_kg/extraction/cluster/validation.py` — no separate validation stage
- `multiomics_kg/extraction/rag.py` — RAG infrastructure removed (pilot showed no benefit)

### Files to create
- `multiomics_kg/extraction/cluster/extract.py` — Stage 1: per-analysis extraction via Responses API
- `multiomics_kg/extraction/cluster/pdf_upload.py` — Files API upload (`purpose="user_data"`) + fallback

### Files to modify
- `multiomics_kg/extraction/cluster/pipeline.py` — rewire to 2-stage flow (table → extract)
- `multiomics_kg/extraction/cluster/prompts.py` — replace 4 prompts with 1 (extract)
- `multiomics_kg/extraction/cluster/table.py` — keep as-is (Stage 0)
- `multiomics_kg/extraction/cluster/run_manager.py` — update stage file names (stage0/1/2)
- `multiomics_kg/review/cluster_review_app.py` — new layout, remove validation display
- `multiomics_kg/review/review_components.py` — new components (cluster list, history, self-assessment display)
- `multiomics_kg/review/review_data.py` — adapt to new stage numbering
- `multiomics_kg/adapters/cluster_adapter.py` — adapt to new stage numbering

### Files unchanged
- `multiomics_kg/extraction/pdf_utils.py` — kept for emergency base64 fallback
- `multiomics_kg/extraction/cluster/run_manager.py` — mostly unchanged (stage file names update)

---

## Part 4: Migration

### Existing extraction data

Current runs under `.extraction_cache/{entry}/runs/` have the old stage numbering (stage1_merged, stage2_results, stage3_validation, stage4_review). Options:

**Approach: Fresh runs, keep history.** Don't migrate old data. Old runs stay in their directories with old filenames. New runs use new filenames (`stage0_table.json`, `stage1_extraction.json`, `stage2_review.json`). The `current` symlink points to the newest run. The adapter detects format by checking which filenames exist.

Review data (`stage2_review.json` / old `stage4_review.json`) is carried forward by RunManager when cluster keys match and input hashes are unchanged.

### Shared cache cleanup

Old RAG artifacts (`chunks.json`, `embeddings.npy`, `pdf_text.json`) remain in `shared/` but are no longer generated or used. Can be cleaned up manually.

### Backward-compatible JSON

Continue writing `cluster_extraction_{entry_key}.json` (the monolithic file the adapter falls back to). Generate from `stage1_extraction.json` — same content, just the node fields.

---

## Part 5: Testing

### Unit tests
- Stage 0: test CSV parsing + enrichment for Tolonen (existing tests, keep)
- Stage 1: test prompt construction from paperconfig context
- Stage 1: test Pydantic schema validation of mock responses
- Stage 1: test cluster key matching from model output names
- Run manager: test stage file renaming + review forwarding

### Integration test
- Run full pipeline on Tolonen MIT9313 (7 clusters) — compare with pilot results
- Verify adapter reads new format correctly
- Verify backward-compatible JSON is written

### Quality comparison
- After first full run: compare all 115 clusters on self-assessment distribution, description richness, parse success rate
- Regression test: Tolonen clusters 1, 6, 7 must match ground truth (direction, enrichment)

---

## Appendix: Pilot Data

### Workload

| Paper | Analyses | Clusters | Status |
|---|---|---|---|
| Tolonen 2006 | 2 (MED4 K=9, MIT9313 K=7) | 16 | ~4 populated |
| Zinser 2009 | 1 (MED4 diel K=18) | 18 | 0 |
| Steglich 2013 | 1 (MED4 expression level) | 6 | 0 |
| Lindell 2007 | 1 (MED4 phage groups) | 2 | 0 |
| Arias-Castro (mbio.03425-22) | 1 (MIT9301 thermal K=5) | 5 | 0 |
| Biller (mSystems.00040-18) | 3 (NATL2A + MIT1002 periodicity/darkness) | 18 | 0 |
| Weissberg (ycae131) | 2 (parental + darktolerant diel K=15) | 30 | 0 |
| Bernstein 2017 | 4 (BP1/M.ruber × light/oxygen) | 20 | 0 |
| **Total** | **15** | **115** | **~4 populated** |

### Cost comparison: old vs new pipeline

| | Old pipeline (spec v1) | New pipeline (pilot-validated) |
|---|---|---|
| Architecture | 4-stage per-cluster | 2-stage per-analysis |
| API calls per 115 clusters | ~230 (extract) + ~115 (validate) = ~345 | ~15 |
| PDF token waste | PDF sent ~345 times | PDF sent ~15 times |
| RAG infrastructure | Yes (embed, chunk, retrieve) | No |
| Projected cost (Flex) | ~$5+ | ~$0.11 |
| Parse failure rate | ~60% (gpt-4o prose) | 0% (strict structured outputs) |
| Cluster 6 direction | WRONG (up) | CORRECT (down) |
