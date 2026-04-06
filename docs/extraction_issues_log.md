# Cluster Extraction Pipeline — Issues Log

Tracks diagnosed issues from review UI sessions and the current state of the extraction system.

## Current State (2026-04-06)

### Production extraction pipeline (single-call approach)

The production extraction uses `multiomics_kg/extraction/cluster/extract.py` — one OpenAI API call per `gene_clusters` entry (or per paper if `extraction.scope: paper`). Full PDF uploaded via Files API, Pydantic structured output.

**Architecture:**
- `SHARED_RULES` + `TYPE_RULES` dict (4 types: `time_course`, `diel`, `condition_comparison`, `classification`) + `SELF_VERIFICATION` block
- `build_prompt(table_config, cluster_summaries)` assembles type-specific prompts
- Paperconfig hints: top-level `extraction:` section (`scope`, `additional_pdfs`), per-entry `extraction_notes` and `figure_hint`
- Model: `gpt-4.1-mini` (default), configurable via `--model` or `CLUSTER_EXTRACTION_MODEL` env var

**Extraction fields:**
- `functional_description` — gene identity/pathway membership (2-3 sentences)
- `temporal_pattern` — expression dynamics prose (1-2 sentences). KG schema still uses `behavioral_description` — rename deferred.
- `expression_dynamics` — short free-text label (scaffolding, not in KG)
- Sentinel: `"N/A"` for undiscussed clusters

**Data:** 15 entries across 8 papers, 110 clusters (~60 described, ~50 N/A). Manually reviewed against papers and corrected (2026-04-06). 12 warnings remaining (all intentional low-confidence temporal patterns).

**Tests:** 1352 unit tests + 51 paperconfig validation tests, all passing.

### Legacy 4-stage pipeline (2026-04-04, not in production)

The 4-stage pipeline (`pipeline.py`, `run_manager.py`, review UI) from 2026-04-04 is still in the codebase but not used for production extraction. It was superseded by the simpler single-call approach which proved more reliable. Infrastructure includes `run_manager.py`, `migrate.py`, Streamlit review UI, RAG experiment tool, and dry-run prompt generator. See git history for details.

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
