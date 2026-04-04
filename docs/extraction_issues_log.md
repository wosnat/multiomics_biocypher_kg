# Cluster Extraction Pipeline — Issues Log

Tracks diagnosed issues from review UI sessions. Each entry includes the root cause, quick fix applied, and longer-term fix needed.

## Issue 1: Semantic path RAG returns supplementary table rows instead of paper text

**Diagnosed:** 2026-04-04, Tolonen 2006 review
**Severity:** High — affects most clusters across all papers
**Symptoms:** Supporting quotes from semantic path contain gene lists, p-value matrices, and spreadsheet fragments instead of paper discussion text. Example: cluster 4 quotes are gene rows like `PMT0411 4 Hypothetical No similarity` instead of the paper's description of cluster composition and function.

**Root cause:** `semantic.py` lines 66-76 mix paper PDF text with `scan_supplementary_text()` output (XLS/TXT table dumps) into a single RAG corpus. When querying for cluster-specific terms, tabular data with high keyword density (gene names, cluster numbers) scores higher than narrative text. The supplementary table rows drown out paper discussion.

**Quick fix (2026-04-04):** Remove `scan_supplementary_text()` from the semantic path. The table path already handles supplementary data at `very_high` confidence — the semantic path's job is finding paper narrative only.

**Longer-term fix:** Separate corpora approach:
- Embed paper text and supplementary text in separate vector stores
- Retrieve from each with different weights (paper text > supplementary)
- Or: detect and filter tabular chunks (high numeric density, tab-separated) before embedding
- Consider chunk quality scoring: penalize chunks that look like raw data vs prose

## Issue 2: Synthesis prompt produces wrong content style for cluster descriptions

**Diagnosed:** 2026-04-04, Tolonen 2006 review
**Severity:** Medium — affects all clusters
**Symptoms:** Synthesis output includes treatment conditions (which belong on the analysis, not individual clusters), overly precise p-values, long gene lists, and verbose descriptions.

**Required changes to synthesis prompt / post-processing:**

1. **No treatment info in cluster descriptions** — treatment is a property of the ClusteringAnalysis node, not individual clusters. The synthesis prompt should not include experimental conditions.
2. **P-values: max 3 decimal places or scientific notation** — e.g., `p=0.013` or `p=7.9e-10`, not `p=1.34e-02` with unnecessary precision.
3. **Gene lists: only named genes, 3-5 max** — only genes mentioned by name in the paper text. Not full gene lists from the CSV. E.g., "including urtA, cynA, and the nitrite permease" not a dump of all 5 genes.
4. **Mention gene categories/pathways/enrichments** — if the paper discusses functional enrichment for a cluster, include that.
5. **Short and to the point** — 2-3 sentences max per description field (functional_description, behavioral_description).

**Quick fix:** Update `SYNTHESIS_PROMPT` in `prompts.py` to specify these constraints explicitly.

**Longer-term fix:** Same prompt changes, plus post-processing validation that checks description length and flags violations.

## Issue 3: RAG query uses locus tags and misses cluster-specific text

**Diagnosed:** 2026-04-04, Tolonen 2006 RAG experiments
**Severity:** High — directly affects which paper text the LLM sees

**Root cause:** `build_cluster_query()` included gene locus tags from the CSV (e.g., PMM0970) which never appear in paper text. Papers use gene names (urtA, cynA). Also, no post-retrieval filtering — generic methodology chunks scored high.

**Quick fix (2026-04-04):**
1. Query changed to `"{organism_short} cluster {N} {enrichment_category} {direction}"` — all terms from table path that also appear in paper text. Gene IDs dropped.
2. Post-retrieval tiered filtering: retrieve 2x, then filter:
   - Tier 1: chunks mentioning this specific cluster (e.g., "cluster 7")
   - Tier 2: chunks mentioning "cluster" in general
   - Combine specific-first, capped at top_k
   - Fallback to unfiltered if no matches

**Longer-term fix:** Consider BM25 + embedding hybrid retrieval. Keyword match for "cluster N" is deterministic and cheap.

## Issue 4: Stage 3 validation conflates clusters when validating all at once

**Diagnosed:** 2026-04-04, Tolonen 2006 review
**Severity:** High — validator returns wrong verdicts, missing clusters

**Symptoms:** Stage 3 validation returns only 1 entry instead of 9, with key "4" but explanation describing cluster 9. Most clusters have no validation verdict.

**Root cause:** All 9 cluster descriptions sent to gpt-4o in a single call with 15 PDF pages. Too much context — the LLM conflates clusters in its JSON response, drops most entries.

**Quick fix (2026-04-04):** Validate one cluster at a time. Each cluster gets its own API call with the same PDF pages + CSV context but only its own description. More API calls but no cross-contamination, and every cluster gets a verdict.

**Longer-term fix:** Same per-cluster approach. Consider caching PDF page embeddings to reduce re-upload cost. Could also use a stronger model for validation since it's the quality gate.
