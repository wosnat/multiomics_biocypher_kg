# Cluster Extraction Pipeline — Issues Log

Tracks diagnosed issues from review UI sessions and the current state of the extraction system.

## Current State (2026-04-06)

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
**Status:** MOSTLY RESOLVED (2026-04-06)

**Original problem:** CSV used numeric cluster IDs (0-4) that didn't match the paper's letter-based clusters (A-D, E-H). Model couldn't map between them.

**Fix applied (2026-04-06):** Changed `cluster_col` from `ClustID_light`/`ClustID_ox` to `Clust name_light`/`Clust name_ox` — the CSV has both columns. Now clusters are A-D (light) and E-H (oxygen), matching the paper. Added `figure_hint` and `extraction_notes` explaining the joint organism context.

**Result:** Bernstein clusters now have real descriptions with irradiance/oxygen response patterns from Figure 6. One remaining locus tag warning (tll1454 in bp1_oxygen H).

**Remaining:** The paper discusses clusters as joint T. elongatus + M. ruber, but our schema splits them per-organism. The multi-organism cluster modeling question (Issue 7 original scope) is still open for future design.

### Issue 8: Model interprets rather than cites — hard to verify accuracy

**Diagnosed:** 2026-04-05
**Severity:** Medium — affects trust in descriptions downstream
**Status:** MITIGATED (2026-04-06)

**Fix applied:** Added self-verification block to prompt — model checks each quote against the specific cluster before outputting. Combined with per-type prompts and stricter "only cite what paper says about THIS cluster" rules.

**Result:** Cross-cluster attribution reduced. Biller 2018 periodicity clusters now use appropriate classification language. Still a known limitation — LLM-generated descriptions are not human-verified.

### Issue 9: Undiscussed clusters get meaningless names and directions

**Diagnosed:** 2026-04-05
**Severity:** Low — cosmetic but noisy
**Status:** RESOLVED (2026-04-06)

**Fix:** Replaced `direction` (constrained enum) with `expression_dynamics` (free-text). Undiscussed clusters now get `expression_dynamics: "N/A"` and simpler names. Coe 2024 clusters correctly show N/A across the board.

### Issue 11: `behavioral_description` field name invites interpretation

**Diagnosed:** 2026-04-05
**Severity:** Low — prompt workaround in place
**Status:** RESOLVED in extraction (2026-04-06). KG schema rename deferred.

**Fix:** Renamed `behavioral_description` → `temporal_pattern` in extraction Pydantic schema. Combined with per-type prompt rules that tell the model exactly what to put in this field per cluster type. KG schema (`schema_config.yaml`) and adapter still use `behavioral_description` — rename deferred to separate task.

### Issue 10: Single prompt doesn't fit all cluster types

**Diagnosed:** 2026-04-05
**Severity:** Medium — affects extraction quality for periodicity and expression level clusters
**Status:** RESOLVED (2026-04-06)

**Fix:** Consolidated 7 cluster types → 4 (`time_course`, `diel`, `condition_comparison`, `classification`). Replaced single `DEVELOPER_MSG_TEMPLATE` with `SHARED_RULES` + `TYPE_RULES` dict + `SELF_VERIFICATION` block. Each type has specific guidance for `expression_dynamics` and `temporal_pattern` fields plus type-appropriate few-shot examples.

**Result:** Warnings dropped from 27 → 6. Biller periodicity uses classification language, Zinser uses diel timing, Tolonen uses time-course dynamics.

---

## Next Steps (TODO)

1. **Fix remaining 6 warnings** — 2 locus tags (PMM0958 in Tolonen, tll1454 in Bernstein), 2 filler phrases, 2 near-identical descriptions. Minor prompt tweaks or accept as-is.
2. **KG schema + adapter rename** — `behavioral_description` → `temporal_pattern` in `schema_config.yaml` and `cluster_adapter.py`. Deferred from this iteration.
3. **Multi-organism cluster modeling** (Issue 7 remainder) — Bernstein has per-organism splits of joint clusters. Design needed for how to represent joint clusters in KG.
4. **Human review** — spot-check extracted descriptions against papers for accuracy, especially Bernstein (new) and Biller classification entries.

## Spec & Plan

- **Spec (2026-04-06):** `docs/superpowers/specs/2026-04-06-extraction-per-type-prompts-design.md`
- **Plan (2026-04-06):** `docs/superpowers/plans/2026-04-06-extraction-per-type-prompts.md`
- **Spec (2026-04-05):** `docs/superpowers/specs/2026-04-05-extraction-quality-iteration-design.md`
- **Spec (2026-04-04):** `docs/superpowers/specs/2026-04-04-cluster-extraction-review-system-design.md`
- **Plan (2026-04-04):** `docs/superpowers/plans/2026-04-04-cluster-extraction-review-system.md`
