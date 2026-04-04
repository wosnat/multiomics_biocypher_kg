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
