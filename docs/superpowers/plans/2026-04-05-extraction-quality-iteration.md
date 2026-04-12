# Extraction Quality Iteration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve cluster extraction quality by updating the Pydantic schema, prompt, verification checks, and re-extracting all 15 entries.

**Architecture:** All changes are in `multiomics_kg/extraction/cluster/extract.py` (schema, prompt, context builder, verification) and one test file. Then re-extract all 15 entries and generate a baseline report.

**Tech Stack:** Python, Pydantic, OpenAI Responses API (gpt-4.1-mini)

**Spec:** `docs/superpowers/specs/2026-04-05-extraction-quality-iteration-design.md`

---

### Task 1: Update Pydantic schema

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:39-62`

- [ ] **Step 1: Remove `peak_time_hours` and `period_hours`, add `source_figures`**

In `multiomics_kg/extraction/cluster/extract.py`, replace the `ClusterExtraction` class:

```python
class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    behavioral_description: str
    direction: Literal["up", "down", "mixed", "not_described"]
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    source_figures: list[str]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str
```

Changes from current:
- Removed: `peak_time_hours: Optional[float]` and `period_hours: Optional[float]`
- Added: `source_figures: list[str]`

- [ ] **Step 2: Verify syntax**

Run: `uv run python -c "from multiomics_kg.extraction.cluster.extract import ClusterExtraction; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: update Pydantic schema — remove timing fields, add source_figures"
```

---

### Task 2: Update prompt with rules, few-shot examples, and cluster type guidance

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:68-92` (DEVELOPER_MSG_TEMPLATE)
- Modify: `multiomics_kg/extraction/cluster/extract.py:95-112` (build_context_block)

- [ ] **Step 1: Replace DEVELOPER_MSG_TEMPLATE**

Replace the entire `DEVELOPER_MSG_TEMPLATE` constant (lines 68-92) with:

```python
DEVELOPER_MSG_TEMPLATE = """\
You are extracting structured descriptions of gene expression clusters from \
a scientific paper.

{context_block}

For each cluster, extract all fields in the output schema.

## Rules

- functional_description: What genes/pathways are in this cluster?
  Can be: enrichment category with p-value, specific highlighted genes, both,
  or "Not discussed in paper." when not mentioned.
  2-3 sentences max. Only cite gene names from the paper (e.g., psbA, rbcLS),
  never locus tags (PMM*, PMT*, P9301_*, NATL2_*, PMN2A_RS*, tll*, SY28_*,
  BSR22_*, ALT831_RS*, MIT1002_*, SYNW*, sync_*, A9601_*, WP_*, cds-*).

- behavioral_description: How do genes in this cluster behave?
  Describe the expression dynamics — not just "up" or "down" but the pattern:
  timing (peak hours, periodicity), kinetics (rapid/gradual, transient/sustained),
  condition-dependence (increases with irradiance, decreases with oxygen).
  Include timing numbers when available from the paper. 1-2 sentences.

- If the paper does NOT describe a cluster's function or behavior, set the \
description to "Not discussed in paper." — an explicit statement. Do NOT \
invent, speculate, or generate generic descriptions. Each field is independent: \
a cluster can have a functional description but behavioral = \
"Not discussed in paper.", or vice versa. Partial descriptions are fine — \
describe what the paper says.

- Do NOT include treatment conditions in descriptions — those live on the analysis node.
- self_assessment: your confidence. assessment_notes: what you're uncertain about.
- Each cluster must have a unique id in snake_case: {{organism_short}}_{{direction}}_{{theme}}.
- name format: "{{Organism}} cluster {{KEY}} ({{direction}}, {{theme}})" — under 60 chars.
  Use the EXACT cluster key from the list below.
- For enrichment: use p-values from the paper/figures, max 3 decimal places or scientific notation.
- Max 3-5 named genes per cluster description.
- supporting_quotes: direct quotes from the paper that support your description.
- source_figures: list of figure/table references you used as evidence (e.g., "Figure 4A", "Table S4").

{cluster_type_guidance}

## Examples

Example 1 — Diel cluster with enrichment and timing:
{{"id": "med4_up_photosynthesis", "name": "Prochlorococcus cluster 1 (up, photosynthesis)", \
"functional_description": "Enriched for photosystem I and II components (p=1.5e-9). \
Includes genes psbA, psbD, and psaA involved in photosynthetic light reactions.", \
"behavioral_description": "Genes peak in expression near dawn (8.3 h) with 24h periodicity, \
coinciding with the onset of light and photosynthetic activity.", \
"direction": "up", "enrichment_category": "Photosystem I and II", \
"enrichment_pvalue": 1.5e-09, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", "confidence_notes": "", \
"supporting_quotes": [{{"quote": "Expression of approximately half of photosystem (PS) II genes, \
including reaction center genes psbA and psbD, peak in abundance at mid-day", "location": "Page 6"}}], \
"source_figures": ["Figure 4A", "Table S4"]}}

Example 2 — Stress time-course with dynamics:
{{"id": "mit9313_up_transport_binding", "name": "MIT9313 cluster 1 (up, transport and binding)", \
"functional_description": "Enriched for transport and binding (p=0.04). Contains nitrogen \
transport genes urtA and the nitrite permease, and hli genes hliS and hli7.", \
"behavioral_description": "Most rapidly and strongly upregulated cluster during nitrogen \
starvation, with genes responding within the first hours of the time course.", \
"direction": "up", "enrichment_category": "transport and binding", \
"enrichment_pvalue": 0.04, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", \
"confidence_notes": "Strong functional enrichment and early rapid upregulation.", \
"supporting_quotes": [{{"quote": "Cluster 1, the most rapidly and highly upregulated genes \
in each strain, contains N transport genes such as MED4 and MIT9313 urtA", "location": "Page 3"}}], \
"source_figures": ["Figure 2"]}}

Example 3 — Partially described:
{{"id": "med4_up_nitrogen_metabolism", "name": "Prochlorococcus cluster 7 (up, nitrogen metabolism)", \
"functional_description": "Includes nitrogen metabolism genes such as amt1.", \
"behavioral_description": "Genes peak near sunset (20.1 h), consistent with nitrogen uptake \
during the night.", \
"direction": "up", "enrichment_category": "Nitrogen metabolism", \
"enrichment_pvalue": 0.087, "enrichment_significant": false, \
"self_assessment": "medium", "assessment_notes": "Marginal enrichment; limited description in paper.", \
"confidence_notes": "Enrichment is marginal; functional role inferred from gene annotations and timing.", \
"supporting_quotes": [{{"quote": "Ammonium transport gene (amt1) and assimilation genes \
peak near sunset", "location": "Page 10"}}], \
"source_figures": ["Figure 7B"]}}

Example 4 — Not discussed in paper:
{{"id": "med4_down_cluster_10", "name": "Prochlorococcus cluster 10 (down, not discussed)", \
"functional_description": "Not discussed in paper.", \
"behavioral_description": "Not discussed in paper.", \
"direction": "down", "enrichment_category": "", \
"enrichment_pvalue": null, "enrichment_significant": false, \
"self_assessment": "low", "assessment_notes": "Paper does not discuss this cluster.", \
"confidence_notes": "", "supporting_quotes": [], "source_figures": []}}

CRITICAL: You MUST extract EXACTLY {n_clusters} clusters, one for each cluster key \
listed below. Use the EXACT cluster keys as they appear — do NOT renumber, skip, \
or merge clusters. Every key must appear exactly once in your output.

{cluster_summaries}
"""
```

- [ ] **Step 2: Add CLUSTER_TYPE_GUIDANCE mapping**

Add this constant after `DEVELOPER_MSG_TEMPLATE`:

```python
CLUSTER_TYPE_GUIDANCE = {
    "response_pattern": """\
Cluster type: response_pattern
Behavioral description style: Describe dynamics — rapid/gradual onset, transient/sustained,
timing relative to treatment. E.g. "Rapidly upregulated within 6h, then sustained." """,
    "diel_cycling": """\
Cluster type: diel_cycling
Behavioral description style: Describe timing — peak hours, periodicity, phase relative
to light/dark. E.g. "Peaks near dawn (8.3 h), 24h periodicity." """,
    "diel_expression_pattern": """\
Cluster type: diel_expression_pattern
Behavioral description style: Describe timing — peak hours, periodicity, phase relative
to light/dark. E.g. "Peaks near dawn (8.3 h), 24h periodicity." """,
    "periodicity_classification": """\
Cluster type: periodicity_classification
Behavioral description style: Describe which conditions show periodic expression and which
don't. E.g. "Periodic in coculture L:D and darkness, not in axenic." """,
    "expression_level": """\
Cluster type: expression_level
Behavioral description style: Describe constitutive vs. variable expression across
conditions. E.g. "Consistently highly expressed across all growth conditions." """,
    "expression_classification": """\
Cluster type: expression_classification
Behavioral description style: Describe presence/absence pattern across conditions.
E.g. "Transcripts present in both axenic and coculture during extended darkness." """,
    "expression_pattern": """\
Cluster type: expression_pattern
Behavioral description style: Describe expression dynamics across conditions — which
conditions drive changes and how. E.g. "Upregulated at cold temperatures during daytime." """,
}
```

- [ ] **Step 3: Update `build_context_block()` to include new fields and cluster type guidance**

Replace `build_context_block` function:

```python
def build_context_block(table: dict) -> str:
    """Build context block from paperconfig entry."""
    parts = [
        f"Analysis: {table.get('name', '')}",
        f"Organism: {table.get('organism', '')}",
        f"Clustering: {table.get('cluster_method', '')}",
        f"Type: {table.get('cluster_type', '')}",
        f"Treatment: {table.get('treatment', '')}",
    ]
    if table.get("treatment_type"):
        tt = table["treatment_type"]
        if isinstance(tt, list):
            parts.append(f"Treatment categories: {', '.join(tt)}")
        else:
            parts.append(f"Treatment categories: {tt}")
    if table.get("background_factors"):
        bf = table["background_factors"]
        if isinstance(bf, list):
            parts.append(f"Background factors: {', '.join(bf)}")
        else:
            parts.append(f"Background factors: {bf}")
    if table.get("experimental_context"):
        parts.append(f"Context: {table['experimental_context']}")
    if table.get("omics_type"):
        parts.append(f"Omics: {table['omics_type']}")
    if table.get("figure_hint"):
        parts.append(f"Key figures: {table['figure_hint']}")
    if table.get("time_points"):
        parts.append(f"Time points (hours): {table['time_points']}")
    return "\n".join(parts)


def get_cluster_type_guidance(table: dict) -> str:
    """Get behavioral description style guidance for this cluster type."""
    ct = table.get("cluster_type", "")
    return CLUSTER_TYPE_GUIDANCE.get(ct, f"Cluster type: {ct}")
```

- [ ] **Step 4: Update `extract_analysis()` to pass cluster type guidance**

In `extract_analysis()`, change the `dev_msg` construction:

```python
    dev_msg = DEVELOPER_MSG_TEMPLATE.format(
        context_block=ctx,
        n_clusters=len(cluster_summaries),
        cluster_summaries=summaries,
        cluster_type_guidance=get_cluster_type_guidance(table_config),
    )
```

- [ ] **Step 5: Update `extract_paper()` to pass cluster type guidance**

In `extract_paper()`, replace the loop and `dev_msg` construction. The existing loop populates `context_parts` and `all_summaries` — add `guidance_parts` alongside:

```python
    context_parts = []
    all_summaries = []
    guidance_parts = []
    total_clusters = 0
    for table_config, cluster_summaries in tables_and_summaries:
        context_parts.append(build_context_block(table_config))
        all_summaries.append(format_cluster_summaries(cluster_summaries))
        total_clusters += len(cluster_summaries)
        guidance_parts.append(get_cluster_type_guidance(table_config))

    dev_msg = DEVELOPER_MSG_TEMPLATE.format(
        context_block="\n\n".join(context_parts),
        n_clusters=total_clusters,
        cluster_summaries="\n\n".join(all_summaries),
        cluster_type_guidance="\n\n".join(guidance_parts),
    )
```

- [ ] **Step 6: Verify prompt renders**

Run: `uv run python -c "
from multiomics_kg.extraction.cluster.extract import DEVELOPER_MSG_TEMPLATE, get_cluster_type_guidance
table = {'cluster_type': 'diel_cycling'}
result = DEVELOPER_MSG_TEMPLATE.format(
    context_block='Test context',
    n_clusters=3,
    cluster_summaries='Cluster 1: 50 genes',
    cluster_type_guidance=get_cluster_type_guidance(table),
)
print('Lines:', len(result.splitlines()))
assert 'Not discussed in paper' in result
assert 'Example 4' in result
assert 'diel_cycling' in result
print('OK')
"`

Expected: `Lines: <number>` then `OK`

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: rewrite prompt with few-shot examples, cluster type guidance, and updated rules"
```

---

### Task 3: Update verification checks

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:242-277` (verify_quality function)

- [ ] **Step 1: Replace `verify_quality()` with expanded checks**

Replace the entire `verify_quality` function:

```python
def verify_quality(entries: list[tuple[Path, str, dict, dict]]) -> list[str]:
    """Run programmatic quality checks. Returns list of warning strings."""
    warnings = []

    locus_pat = re.compile(
        r"\b("
        r"PMM\d{3,}|PMT\d{3,}|P9301_\d+|tll\d{3,}|SY28_\d+|BSR22_\d+"
        r"|NATL[12]_\d+|PMN2A_RS\d+|ALT831_RS\d+|MIT1002_\d+"
        r"|SYNW\d{4}|sync_\d{4}|A9601_\d+"
        r"|WP_\d{6,}"
        r"|cds-[A-Z]{2}_\d+"
        r")\b"
    )

    filler_phrases = [
        "likely", "possibly", "not described but",
        "specific functions are not detailed",
        "not explicitly described",
    ]

    for paper_dir, entry_key, table_config, pub in entries:
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue
        paper = pub_name(pub)
        prefix = f"[{paper} / {entry_key}"

        # Check 1: locus tags in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "behavioral_description"):
                text = c.get(field, "")
                m = locus_pat.search(text)
                if m:
                    warnings.append(
                        f"{prefix} / cluster {key}] locus tag in {field}: {m.group()}"
                    )

        # Check 2: filler on low-confidence clusters
        for key, c in clusters.items():
            if c.get("self_assessment") == "low":
                for field in ("functional_description", "behavioral_description"):
                    text = c.get(field, "")
                    if text and text != "Not discussed in paper.":
                        warnings.append(
                            f"{prefix} / cluster {key}] low confidence but {field} "
                            f"is not 'Not discussed in paper.': {text[:60]}..."
                        )

        # Check 3: filler phrases in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "behavioral_description"):
                text = c.get(field, "").lower()
                for phrase in filler_phrases:
                    if phrase in text:
                        warnings.append(
                            f"{prefix} / cluster {key}] filler phrase '{phrase}' in {field}"
                        )
                        break

        # Check 4: near-identical descriptions within analysis
        descs = {}
        for key, c in clusters.items():
            fd = c.get("functional_description", "")
            if fd and fd != "Not discussed in paper." and len(fd) > 50:
                prefix_50 = fd[:50]
                if prefix_50 in descs:
                    warnings.append(
                        f"{prefix} / cluster {key}] near-identical functional_description "
                        f"as cluster {descs[prefix_50]}"
                    )
                else:
                    descs[prefix_50] = key

        # Check 5: empty direction with non-empty description
        for key, c in clusters.items():
            fd = c.get("functional_description", "")
            has_desc = fd and fd != "Not discussed in paper." and len(fd) > 20
            direction = c.get("direction", "")
            if has_desc and not direction:
                warnings.append(
                    f"{prefix} / cluster {key}] has description but empty direction"
                )

    return warnings
```

- [ ] **Step 2: Verify it runs on existing data**

Run: `uv run python -c "
from multiomics_kg.extraction.cluster.extract import verify_quality
from multiomics_kg.extraction.cluster.extraction_utils import find_all_entries
entries = find_all_entries()
warnings = verify_quality(entries)
print(f'{len(warnings)} warnings')
for w in warnings[:10]:
    print(f'  {w}')
"`

Expected: Many warnings (the old extraction data has lots of filler). This confirms the checks work.

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: expand verification — filler detection, near-identical check, broader locus tag regex"
```

---

### Task 4: Update report to show source_figures

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:210-239` (generate_report function)

- [ ] **Step 1: Add source_figures to report output**

In `generate_report()`, add after the `notes` block (after `lines.append(f"**Notes:** {notes}")`):

```python
            figs = c.get("source_figures", [])
            if figs:
                lines.append(f"**Sources:** {', '.join(figs)}")
```

- [ ] **Step 2: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: show source_figures in report"
```

---

### Task 5: Update test fixtures

**Files:**
- Modify: `tests/test_extraction_utils.py:39-57`

- [ ] **Step 1: Remove `peak_time_hours` and `period_hours` from roundtrip test**

In `tests/test_extraction_utils.py`, replace the `test_save_load_roundtrip` function:

```python
def test_save_load_roundtrip(tmp_path):
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction, save_extraction
    metadata = {"paper": "Test 2006", "model": "gpt-4.1-mini"}
    clusters = {
        "1": {
            "id": "test_up",
            "name": "Test cluster 1",
            "functional_description": "Genes involved in transport",
            "behavioral_description": "Upregulated early",
            "source_figures": ["Figure 2"],
        },
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)
    loaded = load_extraction(tmp_path, "test_entry")
    assert loaded == clusters
    md_path = tmp_path / "cluster_extractions" / "test_entry.md"
    assert md_path.exists()
    assert "Test cluster 1" in md_path.read_text()
```

- [ ] **Step 2: Run tests**

Run: `uv run pytest tests/test_extraction_utils.py -v`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add tests/test_extraction_utils.py
git commit -m "test: update extraction utils fixtures — remove timing fields, add source_figures"
```

---

### Task 6: Run full test suite

**Files:** None (verification only)

- [ ] **Step 1: Run unit tests**

Run: `uv run pytest -m "not slow and not kg" -v --tb=short 2>&1 | tail -30`
Expected: All tests pass. Pay attention to any failures in `test_cluster_adapter.py` — the adapter still reads `peak_time_hours`/`period_hours` via `.get()` which returns `None` when missing, so it should be fine.

- [ ] **Step 2: Commit if any fixes were needed**

Only if step 1 required fixes.

---

### Task 7: Re-extract all 15 entries

**Files:** None (runs extraction CLI)

**Prerequisites:** `OPENAI_API_KEY` set in `.env`.

- [ ] **Step 1: Dry run to verify all 15 entries found**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --dry-run`
Expected: 15 entries listed.

- [ ] **Step 2: Re-extract all with force**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --force 2>&1 | tee extraction_log.txt`

This takes ~5-10 minutes and costs ~$2-3. Watch for:
- `FAILED` lines (API errors)
- `Unmatched` / `Missing` lines (key matching issues)
- Final summary: "Done: N clusters, X in + Y out tokens"

Expected: ~115 clusters extracted across 15 entries.

- [ ] **Step 3: Generate report with verification**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --report --verify`
Expected: `data/cluster_extraction_report.md` written. Check warning count — should be significantly fewer than the pre-extraction check from Task 3 Step 2.

- [ ] **Step 4: Review report**

Read `data/cluster_extraction_report.md` and check:
1. Coe 2024 clusters — do undiscussed ones say "Not discussed in paper." instead of filler?
2. Biller 2018 periodicity — no more NATL2_* locus tags?
3. Zinser/Tolonen/Bernstein — descriptions still high quality?
4. Warnings section — what remains?

- [ ] **Step 5: Commit baseline**

```bash
git add data/*/papers_and_supp/*/cluster_extractions/*.json
git add data/*/papers_and_supp/*/cluster_extractions/*.md
git add data/cluster_extraction_report.md
git commit -m "extract: re-extract all 115 clusters with improved prompt (baseline v1)"
```

- [ ] **Step 6: Clean up**

```bash
rm -f extraction_log.txt
```

---

### Task 8: Iterate (if needed)

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py` (prompt tweaks)

This task is open-ended — repeat as needed based on report review from Task 7 Step 4.

- [ ] **Step 1: Identify issues from report**

Look at warnings and spot-check descriptions. Common issues to look for:
- Clusters that should say "Not discussed in paper." but have filler
- Locus tags that slipped through
- Wrong direction or enrichment
- Behavioral descriptions that don't match cluster type

- [ ] **Step 2: Tweak prompt and re-extract specific entries**

```bash
# Re-extract just one paper
uv run python -m multiomics_kg.extraction.cluster.extract --paper "Coe 2024" --force

# Regenerate report
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify

# Diff against last commit
git diff data/cluster_extraction_report.md
```

- [ ] **Step 3: Commit improvements**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git add data/*/papers_and_supp/*/cluster_extractions/*.json
git add data/*/papers_and_supp/*/cluster_extractions/*.md
git add data/cluster_extraction_report.md
git commit -m "extract: prompt iteration N — <describe what changed>"
```
