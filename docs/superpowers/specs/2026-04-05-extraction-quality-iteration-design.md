# Cluster Extraction Quality Iteration

**Date:** 2026-04-05
**Status:** Approved

## Goal

Improve cluster extraction quality by adding few-shot examples to the prompt, fixing the "filler text" problem for undiscussed clusters, removing unreliable structured timing fields, and re-extracting all 15 entries for a consistent baseline.

## Current Problems

1. **Filler text for undiscussed clusters** — The model generates vague, repetitive descriptions like "mixed, dark tolerance metabolism" for clusters the paper doesn't discuss. Should say "Not discussed in paper." instead.
2. **Locus tags in descriptions** — NATL2_*, PMN2A_RS*, PMT* appear in text. Regex only catches some patterns.
3. **No few-shot examples** — The prompt has rules but no concrete examples of good vs. empty output.
4. **`peak_time_hours` and `period_hours` are unreliable** — Only 37/116 clusters have peak_time, 70/116 have period_hours, but many are hallucinated (e.g., all 30 Coe 2024 clusters report `period_hours: 24.0` despite being dark-tolerance clusters with garbage filler descriptions). Real timing data exists mainly for Zinser 2009 (diel) and Alonso-Saez 2023. Timing is better expressed in behavioral_description prose where it can be verified.
5. **No cluster type guidance** — The prompt doesn't tell the model how to adapt behavioral descriptions to different cluster types (diel vs. stress response vs. expression level).

## Description Philosophy

Cluster types vary widely and descriptions should reflect that:

### Functional description

Can take multiple valid forms depending on what the paper provides:
- **Enrichment/pathway**: "Enriched for photosystem I and II components (p=1.5e-9)"
- **Highlighted genes**: "Contains nitrogen transport genes urtA, cynA, and the nitrite permease"
- **Both**: enrichment category + specific genes
- **Nothing**: "Not discussed in paper." — the cluster itself is still useful without a description

Partial descriptions are fine. No stats available is fine. The cluster existing in the KG has value on its own.

### Behavioral description

Varies by cluster type — it's not just "up/down", it's about dynamics:
- **Diel clusters**: "Peaks near dawn (8.3 h), 24h periodicity" — timing and period in prose
- **Stress time-course**: "Rapidly upregulated within 6h, then sustained" or "Initially upregulated then sharp drop after 24h"
- **Coculture response**: "Increased expression with irradiance, decreased with oxygen tension"
- **Expression level**: "Consistently highly expressed across all growth conditions"
- **Not discussed**: "Not discussed in paper."

Behavioral patterns are often shown in figures but not always described in text. Include timing (peak hours, periodicity) in the behavioral text when available from the paper — do NOT extract as separate structured fields.

### When the paper doesn't discuss a cluster

Set `functional_description` and/or `behavioral_description` to `"Not discussed in paper."` — an explicit statement, not an empty string and not vague filler. Each field is independent: a cluster can have a functional description but no behavioral, or vice versa. Do NOT generate generic or speculative text. This is the most important rule.

## Schema Changes

### Extraction-only fields (not in KG)

The following fields are extracted by the LLM but NOT stored on GeneCluster nodes in the KG:
`direction`, `enrichment_category`, `enrichment_pvalue`, `enrichment_significant`, `self_assessment`, `assessment_notes`, `confidence_notes`, `supporting_quotes`, `source_figures`.

These stay in the Pydantic schema and extraction JSON because they serve as scaffolding — forcing the model to extract structured facts and cite sources before writing prose descriptions improves quality (chain-of-thought effect). They also power verification checks (e.g., detecting filler where `enrichment_category` is empty but `functional_description` mentions a pathway).

**`source_figures`** (new field): list of figure/table references the model used as evidence (e.g., `["Figure 4A", "Table S4", "Figure 2B heatmap"]`). Especially important for behavioral descriptions which often come from plots rather than text. Forces the model to point at its source before writing the description.

### Remove `peak_time_hours` and `period_hours` from extraction

These structured fields are unreliable — the model hallucinates values for clusters it knows nothing about. When real timing data exists, it belongs in `behavioral_description` prose where it can be verified against the paper.

Remove from extraction only:
- `multiomics_kg/extraction/cluster/extract.py` — remove from Pydantic schema
- `tests/test_extraction_utils.py` — remove from test fixtures

Schema/adapter/KG changes (removing from GeneCluster node properties) are out of scope — tracked separately.

## Proposed Few-Shot Examples

### Example 1: Diel cluster with enrichment and timing

Source: Zinser 2009, cluster 1 — photosynthesis genes with strong enrichment, clear diel timing.

```json
{
  "id": "med4_up_photosynthesis",
  "name": "Prochlorococcus cluster 1 (up, photosynthesis)",
  "functional_description": "Enriched for photosystem I and II components (p=1.5e-9). Includes genes psbA, psbD, and psaA involved in photosynthetic light reactions.",
  "behavioral_description": "Genes peak in expression near dawn (8.3 h) with 24h periodicity, coinciding with the onset of light and photosynthetic activity.",
  "direction": "up",
  "enrichment_category": "Photosystem I and II",
  "enrichment_pvalue": 1.5e-09,
  "enrichment_significant": true,
  "self_assessment": "high",
  "assessment_notes": "",
  "confidence_notes": "",
  "supporting_quotes": [
    {"quote": "Expression of approximately half of photosystem (PS) II genes, including reaction center genes psbA and psbD, peak in abundance at mid-day", "location": "Page 6"}
  ],
  "source_figures": ["Figure 4A", "Table S4"]
}
```

**Why:** Shows a well-described diel cluster — timing embedded in behavioral prose, enrichment in functional, specific gene names from the paper. Quotes and figure references ground the description.

### Example 2: Stress time-course cluster with dynamics

Source: Tolonen 2006, MIT9313 cluster 1 — nitrogen starvation, rapid early response.

```json
{
  "id": "mit9313_up_transport_binding",
  "name": "MIT9313 cluster 1 (up, transport and binding)",
  "functional_description": "Enriched for transport and binding (p=0.04). Contains nitrogen transport genes urtA and the nitrite permease, and hli genes hliS and hli7.",
  "behavioral_description": "Most rapidly and strongly upregulated cluster during nitrogen starvation, with genes responding within the first hours of the time course.",
  "direction": "up",
  "enrichment_category": "transport and binding",
  "enrichment_pvalue": 0.04,
  "enrichment_significant": true,
  "self_assessment": "high",
  "assessment_notes": "",
  "confidence_notes": "Strong functional enrichment and early rapid upregulation.",
  "supporting_quotes": [
    {"quote": "Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA", "location": "Page 3"}
  ],
  "source_figures": ["Figure 2"]
}
```

**Why:** Shows a stress-response cluster — behavioral describes dynamics (rapid, early), not just direction. Functional has both enrichment and specific genes.

### Example 3: Partially described cluster

Source: Zinser 2009, cluster 7 — marginal enrichment, some detail.

```json
{
  "id": "med4_up_nitrogen_metabolism",
  "name": "Prochlorococcus cluster 7 (up, nitrogen metabolism)",
  "functional_description": "Includes nitrogen metabolism genes such as amt1.",
  "behavioral_description": "Genes peak near sunset (20.1 h), consistent with nitrogen uptake during the night.",
  "direction": "up",
  "enrichment_category": "Nitrogen metabolism",
  "enrichment_pvalue": 0.087,
  "enrichment_significant": false,
  "self_assessment": "medium",
  "assessment_notes": "Marginal enrichment; limited description in paper.",
  "confidence_notes": "Enrichment is marginal; functional role inferred from gene annotations and timing.",
  "supporting_quotes": [
    {"quote": "Ammonium transport gene (amt1) and assimilation genes peak near sunset", "location": "Page 10"}
  ],
  "source_figures": ["Figure 7B"]
}
```

**Why:** Shows it's OK to have short descriptions when the paper has limited detail. Still grounded — no speculation.

### Example 4: Not discussed in paper

Source: Based on Zinser 2009, cluster 10 — paper provides no functional discussion.

```json
{
  "id": "med4_down_cluster_10",
  "name": "Prochlorococcus cluster 10 (down, not discussed)",
  "functional_description": "Not discussed in paper.",
  "behavioral_description": "Not discussed in paper.",
  "direction": "down",
  "enrichment_category": "",
  "enrichment_pvalue": null,
  "enrichment_significant": false,
  "self_assessment": "low",
  "assessment_notes": "Paper does not discuss this cluster.",
  "confidence_notes": "",
  "supporting_quotes": [],
  "source_figures": []
}
```

**Why:** The key fix. Explicit "Not discussed in paper." instead of vague filler. Empty quotes and figures confirm nothing was found. `direction` still extracted when inferrable.

## Enriched Prompt Context

The current `build_context_block()` already reads most fields from the paperconfig (name, organism, cluster_method, cluster_type, treatment, experimental_context, omics_type, time_points). No KG query needed — the paperconfig is the source.

The key improvement is using `cluster_type` to **guide the description style**. Add a mapping in the prompt that tells the model what kind of behavioral description to write:

```
Cluster type guidance (match your behavioral_description style to the cluster type):
- response_pattern → Describe dynamics: rapid/gradual onset, transient/sustained, 
  timing relative to treatment. E.g. "Rapidly upregulated within 6h, then sustained."
- diel_cycling, diel_expression_pattern → Describe timing: peak hours, periodicity,
  phase relative to light/dark. E.g. "Peaks near dawn (8.3 h), 24h periodicity."
- periodicity_classification → Describe which conditions show periodic expression
  and which don't. E.g. "Periodic in coculture L:D and darkness, not in axenic."
- expression_level → Describe constitutive vs. variable expression across conditions.
  E.g. "Consistently highly expressed across all growth conditions."
- expression_classification → Describe presence/absence pattern across conditions.
  E.g. "Transcripts present in both axenic and coculture during extended darkness."
```

Also add to the context block (already available in paperconfig, just not passed to prompt currently):
- `treatment_type` list — helps the model understand what categories to look for
- `background_factors` list — experimental context factors
- Linked experiment names (from paperconfig `experiments` block) — helps ground what the comparison is

## Prompt Changes

### Updated rules

```
- functional_description: What genes/pathways are in this cluster?
  Can be: enrichment category with p-value, specific highlighted genes, both,
  or "Not discussed in paper." when not mentioned.
  2-3 sentences max. Only cite gene names from the paper (e.g., psbA, rbcLS), never locus
  tags (PMM*, PMT*, P9301_*, NATL2_*, PMN2A_RS*, tll*, SY28_*, BSR22_*, ALT831_RS*,
  MIT1002_*, SYNW*, sync_*, A9601_*, WP_*, cds-*).

- behavioral_description: How do genes in this cluster behave?
  Describe the expression dynamics — not just "up" or "down" but the pattern:
  timing (peak hours, periodicity), kinetics (rapid/gradual, transient/sustained),
  condition-dependence (increases with irradiance, decreases with oxygen).
  Include timing numbers when available from the paper. 1-2 sentences.

- If the paper does NOT describe a cluster's function or behavior, set the description
  to "Not discussed in paper." — an explicit statement. Do NOT invent, speculate,
  or generate generic descriptions. Each field is independent: a cluster can have
  a functional description but behavioral = "Not discussed in paper.", or vice versa.
  Partial descriptions are fine — describe what the paper says.

```

### Few-shot block

The four examples above added as `## Examples` section in the prompt, between rules and cluster list.

## Verification Improvements

### Expanded locus tag regex

```python
locus_pat = re.compile(
    r"\b("
    r"PMM\d{3,}|PMT\d{3,}|P9301_\d+|tll\d{3,}|SY28_\d+|BSR22_\d+"
    r"|NATL[12]_\d+|PMN2A_RS\d+|ALT831_RS\d+|MIT1002_\d+"
    r"|SYNW\d{4}|sync_\d{4}|A9601_\d+"
    r"|WP_\d{6,}"
    r"|cds-[A-Z]{2}_\d+"
    r")\b"
)
```

### Filler detection check

Flag clusters where:
- `self_assessment` is `"low"` but description is not `"Not discussed in paper."` (filler instead of explicit)
- Description contains filler phrases: "likely", "possibly", "not described but", "specific functions are not detailed", "not explicitly described"
- Multiple clusters in the same analysis have near-identical descriptions (substring match on first 50 chars, excluding "Not discussed in paper.")

## Execution Plan

1. Update Pydantic schema in `extract.py`: remove `peak_time_hours` and `period_hours`, add `source_figures: list[str]`
2. Update prompt in `extract.py`: new rules + few-shot examples + cluster type guidance
3. Update `build_context_block()`: add `treatment_type`, `background_factors`, cluster type guidance mapping
4. Update verification in `extract.py`: expanded locus tag regex + filler detection
5. Update test fixtures in `tests/test_extraction_utils.py`: remove `peak_time_hours`/`period_hours`
6. Run `pytest -m "not slow and not kg"` to confirm nothing breaks
7. Re-extract all 15 entries with `--force`
8. Generate report with `--report --verify`
9. Review report, iterate if needed (tweak prompt, re-extract specific entries)

Estimated cost: ~$2-3 with gpt-4.1-mini for all 15 entries.
