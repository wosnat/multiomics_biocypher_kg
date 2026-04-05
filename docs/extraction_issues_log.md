# Cluster Extraction Pipeline — Issues Log

Tracks diagnosed issues from review UI sessions and the current state of the extraction system.

## Current State (2026-04-05)

### What was built (2026-04-04)

**Infrastructure:**
- `multiomics_kg/extraction/cluster/run_manager.py` — versioned run directories with `create_run()` / `finalize_run()`, stage file I/O, input hashing, review copy-forward
- `multiomics_kg/extraction/cluster/migrate.py` — one-time migration of legacy monolithic `cluster_extraction_*.json` into per-stage files
- Shared disk cache at `{paper_dir}/.extraction_cache/shared/` — `pdf_text.json`, `chunks.json`, `embeddings.npy`, `pdf_content_parts.json`
- Pipeline refactored (`pipeline.py`) to use RunManager, write per-stage files, finalize symlink after completion

**Review UI (Streamlit):**
- `multiomics_kg/review/cluster_review_app.py` — sidebar navigation with color-coded paper/entry dropdowns, verdict/review filters, export issue report, re-run triggers
- `multiomics_kg/review/review_components.py` — 3-column merge view (table|visual|semantic), synthesis result below, review controls (status/issues/failing stages/notes/editable fields), diff view, source access (inline tables, PDF links)
- `multiomics_kg/review/review_data.py` — paper scanning, summary stats, status color computation, issue report export
- Launch: `uv run streamlit run multiomics_kg/review/cluster_review_app.py`

**Tools:**
- `scripts/rag_experiment.py` — interactive RAG query testing with disk-cached embeddings. Usage: `uv run python scripts/rag_experiment.py "data/.../paper_dir" --interactive`
- `scripts/dry_run_prompts.py` — generates all LLM prompts without API calls. Usage: `uv run python scripts/dry_run_prompts.py "data/.../paperconfig.yaml" --output-dir /tmp/dry_run`

**Adapter updated:**
- `multiomics_kg/adapters/cluster_adapter.py` — reads from `.extraction_cache/{entry}/current/` structure, respects review status (approve/edit use descriptions, reject/stale/flag-issue use empty), applies `edited_fields` overrides

**Tests:** 53 targeted + 1356 full suite, all passing.

### Architecture: per-cluster LLM calls

All four LLM stages now run one call per cluster (no batch calls):

| Stage | Model | Per-cluster? | Shared context |
|---|---|---|---|
| Visual | gpt-4o | Yes | PDF pages (cached to `pdf_content_parts.json`) |
| Semantic | gpt-5-nano | Yes | RAG chunks (cached to `chunks.json` + `embeddings.npy`) |
| Synthesis | gpt-5-nano | Yes | Merged stage1 data (no PDF) |
| Validation | gpt-4o | Yes | PDF pages + CSV summary (same cache as visual) |

PDF prefix caching: OpenAI automatically caches identical prompt prefixes. Since all per-cluster calls within a stage share the same PDF pages, calls 2-N use cached input tokens (50% cheaper).

### Latest Tolonen run results (2026-04-04 ~23:00)

- **MED4**: 9 synthesized, **3 pass, 6 fail** — fails are JSON parse errors (see Issue 6)
- **MIT9313**: 7 synthesized, **1 pass, 2 warn, 4 fail** — same issue

Data is at `.extraction_cache/{entry}/current/` — viewable in the review UI.

---

## Open Issues

### Issue 6: Validator returns prose instead of JSON (BLOCKING)

**Diagnosed:** 2026-04-05
**Severity:** High — 6/9 MED4 and 4/7 MIT9313 clusters fail validation due to parse errors
**Status:** NOT FIXED

**Symptoms:** Validation LLM (gpt-4o) returns prose like "To validate the extracted description for Cluster 1, let's refer to the paper..." instead of the requested JSON. The `_validate_single_cluster` function logs "JSON parse error: Expecting value: line 1 column 1".

**Root cause:** The prompt asks for JSON but doesn't enforce it. gpt-4o sometimes "thinks out loud" before outputting JSON, or skips JSON entirely.

**Fix:** Add `response_format={"type": "json_object"}` to the `client.chat.completions.create()` call in `validation.py` `_validate_single_cluster()`. This forces the model to output valid JSON. May also want to add it to `_extract_single_cluster` in `visual.py` and `_extract_from_passages` in `semantic.py`.

**Where:** `multiomics_kg/extraction/cluster/validation.py` line ~82 (the `client.chat.completions.create` call).

### Issue 5: PDF sent as base64 in every API call — use Files API instead

**Diagnosed:** 2026-04-04
**Severity:** Medium — cost/performance optimization, not correctness
**Status:** NOT FIXED (logged for future)

**Current state:** Visual and validation send base64-encoded PDF pages in every per-cluster API call. Disk cache avoids re-encoding but doesn't reduce API payload size.

**Fix:** Use OpenAI Files API (`POST /v1/files` with `purpose="user_data"`) to upload PDFs once, get a `file_id`, and reference it in all subsequent calls. See https://developers.openai.com/api/docs/guides/file-inputs.

---

## Resolved Issues

### Issue 1: Semantic path RAG returns supplementary table rows instead of paper text (FIXED)

**Diagnosed:** 2026-04-04
**Fix:** Removed `scan_supplementary_text()` from semantic path. Table path handles supplementary data at `very_high` confidence.
**File:** `multiomics_kg/extraction/cluster/semantic.py`

### Issue 2: Synthesis prompt produces wrong content style (FIXED)

**Diagnosed:** 2026-04-04
**Fix:** Updated `SYNTHESIS_PROMPT` in `prompts.py`:
- No treatment in descriptions (belongs on analysis node)
- P-values: max 3 decimals or scientific notation
- Only named genes from paper text (not locus tags), max 3-5
- 2-3 sentences functional, 1-2 sentences behavioral
- Explicit instruction to ignore locus tag IDs (PMM*, PMT*)
**File:** `multiomics_kg/extraction/cluster/prompts.py`

### Issue 3: RAG query uses locus tags and misses cluster-specific text (FIXED)

**Diagnosed:** 2026-04-04
**Fix:**
1. Query changed to `"{organism_short} cluster {N} {enrichment_category} {direction}"` — gene IDs dropped
2. Post-retrieval tiered filtering: retrieve 2x, filter tier 1 (specific cluster mention) then tier 2 (any "cluster" mention), combine specific-first
**File:** `multiomics_kg/extraction/cluster/semantic.py`

### Issue 4: Validation conflates clusters when validating all at once (FIXED)

**Diagnosed:** 2026-04-04
**Fix:** Validate one cluster at a time. Each cluster gets its own API call with shared PDF pages + CSV context. 2s delay between calls to avoid rate limiting.
**File:** `multiomics_kg/extraction/cluster/validation.py`

### Issue 7: Bernstein 2017 cluster decomposition mismatch

**Diagnosed:** 2026-04-05
**Severity:** Medium — 20 clusters across 4 entries affected, mostly "Not discussed in paper."
**Status:** OPEN

**Symptoms:** `bp1_light_clusters` (5 clusters) all say "Not discussed in paper." despite the paper having detailed cluster descriptions with gene names and expression patterns (Figure 6, pages 6-7). `bp1_oxygen_clusters` clusters 0-2 also empty. `mruber_light_clusters` 0 and 4 empty. `mruber_oxygen_clusters` 0 and 4 empty.

**Root cause:** The paper defines clusters A-D (irradiance-responsive) and E-H (pO2-responsive), each containing genes from **both** T. elongatus and M. ruber jointly. Our paperconfig splits these into 4 per-organism × per-condition entries with numbered clusters 0-4. The model correctly cannot map the paper's joint "cluster A" to our per-organism "bp1_light cluster 0" — they're different decompositions.

The old extraction hallucinated descriptions by guessing which paper cluster maps to which CSV cluster. The new prompt correctly says "Not discussed" since the paper doesn't discuss per-organism numbered clusters.

**Options:**
- A) Accept "Not discussed" — the per-organism decomposition doesn't exist in the paper
- B) Add `figure_hint` or `cluster_mapping` to paperconfig telling the model the correspondence
- C) Restructure to joint organism clusters matching the paper (A-D, E-H)

**Scope:** Bernstein 2017 is currently the only paper with multi-organism clusters. Broader question: how should multi-organism clusters be modeled in the KG? Current schema ties ClusteringAnalysis to a single organism. Joint clusters would need a new pattern (multiple organism links, or a "community" clustering type).

**Decision:** Deferred. Needs design exploration before implementation.

### Issue 8: Model interprets rather than cites — hard to verify accuracy

**Diagnosed:** 2026-04-05
**Severity:** Medium — affects trust in descriptions downstream
**Status:** OPEN

**Symptoms:** Functional descriptions contain plausible-sounding content that may be the model's interpretation rather than direct paper statements. Example: Biller 2018 `mit1002_periodicity` cluster `coculture_LD` describes "Calvin cycle, glycolysis, fatty acid biosynthesis" — the supporting quote mentions these pathways but in the context of a *different* cluster category (extended darkness, not L:D only). The model applied the right information to the wrong cluster.

**Root cause:** The model synthesizes information from across the paper rather than strictly citing what the paper says about each specific cluster. Without human verification, it's hard to tell what's directly from the paper vs. model interpretation. The original 4-stage pipeline had a validation stage ("judge") that was supposed to catch this, but it was removed in the cleanup.

**Additional context:** Biller 2018 periodicity entries are *classification* (periodic Y/N across conditions), not clustering in the traditional sense. The paper discusses periodicity patterns and pathway enrichments for sets of genes, but the mapping to specific composite condition-clusters (e.g., `coculture_LD+coculture_darkness`) is indirect.

**Possible mitigations:**
- Strengthen prompt: "Only state what the paper explicitly says about this cluster. Do not synthesize or infer across clusters."
- Re-introduce a lightweight verification step (no separate LLM call — just a prompt section asking the model to self-check each quote against the cluster it's attributed to)
- Accept as known limitation and label descriptions as "LLM-interpreted, not human-verified"

**Decision:** Prompt improvement to reduce interpretation is worth trying. Full judge/verification is out of scope for now.

---

## Next Steps (TODO)

1. **Fix Issue 6** (JSON response format) — add `response_format={"type": "json_object"}` to validation, and optionally visual + semantic API calls
2. **Re-run Tolonen** with the JSON format fix — should get actual verdicts for all 16 clusters
3. **Review in UI** — classify any remaining issues, verify descriptions are accurate
4. **Iterate on prompts** if needed based on review findings
5. **Scale to other papers** — run extraction on remaining 13 `gene_clusters` entries across 7 papers
6. **Implement Files API** (Issue 5) — upload PDFs once, reference by file_id

## Spec & Plan

- **Spec:** `docs/superpowers/specs/2026-04-04-cluster-extraction-review-system-design.md`
- **Plan:** `docs/superpowers/plans/2026-04-04-cluster-extraction-review-system.md`
- **UI spec gaps** identified (2026-04-04): mostly fixed (reviewed_at, carried-forward badges, PDF button). Remaining: per-quote source links partial, diff view doesn't detect improvement/regression direction.
