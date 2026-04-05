# Cluster Extraction Quality — Findings and Next Steps

**Date:** 2026-04-05
**Status:** Findings from first quality iteration session

## What Was Done

### Prompt improvements (committed)
1. Removed `peak_time_hours` and `period_hours` from Pydantic schema — unreliable, timing now in behavioral prose
2. Added `source_figures` field — forces model to cite figure/table sources
3. Added 4 few-shot examples (diel, stress, partial, not discussed)
4. Added cluster type guidance mapping (`CLUSTER_TYPE_GUIDANCE`) — tells model what behavioral style to use per cluster type
5. Enriched context block with `treatment_type`, `background_factors`
6. Expanded locus tag regex (13 patterns vs 6)
7. Added filler detection (low-confidence check, filler phrases, near-identical descriptions)
8. Rule: "Not discussed in paper." instead of filler text
9. Rule: don't describe membership stats (gene counts, sample IDs)
10. Rule: don't do LLM-driven enrichment or gene cherry-picking
11. Rule: no causal interpretation in behavioral descriptions (GOOD/BAD examples)
12. Rule: behavioral = temporal expression pattern only, not biological function
13. Rule: don't say what the paper didn't do — just say "Not discussed"
14. Rule: only state what paper explicitly says about THIS specific cluster

### Results
- **129 → 16-27 warnings** (depending on iteration) across 115 clusters
- **Coe 2024** (30 clusters): all correctly "Not discussed in paper." — was all filler before
- **Biller 2018 periodicity**: improved — no more locus tags, cleaner behavioral descriptions
- **Response pattern papers** (Tolonen, Zinser, Lindell, Alonso-Saez): good quality with real enrichment data and gene names
- **1 missing cluster**: Tolonen MED4 cluster 9 — model dropped it

### Remaining Issues (from issues log)

| Issue | Summary | Severity | Fix approach |
|-------|---------|----------|-------------|
| **7** | Bernstein 2017 multi-organism cluster mismatch | Medium | Needs KG schema design for joint clusters |
| **8** | Model interprets rather than cites (cross-cluster attribution) | Medium | Per-type prompts + lightweight self-verification |
| **9** | Undiscussed clusters get meaningless names/directions | Low | Default `not_discussed` direction, simpler name format |
| **10** | Single prompt doesn't fit all cluster types | Medium | Per-type prompt templates |
| **11** | `behavioral_description` field name invites interpretation | Low | Rename to `temporal_pattern` or `expression_pattern` |

## Proposed Next Steps

### Phase 1: Per-cluster-type prompt templates (Issues 9, 10, 11)

The single biggest improvement would be different prompt sections per cluster type. Current prompt is designed for `response_pattern` clusters and works well there. Other types need fundamentally different guidance.

**Cluster types and what they need:**

| Type | Papers | What to extract | Direction? | Timing? |
|------|--------|----------------|-----------|---------|
| `response_pattern` | Tolonen, Lindell, Bernstein | Enrichment, highlighted genes, dynamics (rapid/gradual) | Yes (up/down) | Relative to treatment |
| `diel_cycling` | Zinser | Enrichment, highlighted genes, peak phase | Peak phase (dawn/dusk/night) | Absolute (hours) |
| `diel_expression_pattern` | Coe 2024 | Not discussed — paper uses different decomposition | Not applicable | Not applicable |
| `periodicity_classification` | Biller 2018 (periodicity) | Which conditions show periodicity, pathway associations | Not applicable | Period (24h) |
| `expression_level` | Wang 2014 | COG category enrichment, evolutionary rates | Level (high/low/variable) | Not applicable |
| `expression_classification` | Biller 2018 (darkness) | Enrichment, presence/absence pattern | Enriched/depleted | Time window |
| `expression_pattern` | Alonso-Saez 2023 | Enrichment, highlighted genes, temperature response | Mixed per condition | Across conditions |

**Implementation:** `PROMPT_RULES_BY_TYPE` dict with type-specific rules, examples, and field guidance. Shared scaffolding fields stay the same. The `direction` field gets type-specific allowed values or is set to `not_applicable`.

**Schema changes needed:**
- Rename `behavioral_description` → `expression_pattern` (or `temporal_pattern`) in Pydantic schema, adapter, KG schema
- Add `not_applicable` to `direction` enum
- Update few-shot examples per type

### Phase 2: Lightweight self-verification (Issue 8)

Add a prompt section asking the model to self-check before outputting:
```
Before outputting, verify for each cluster:
1. Does each supporting_quote explicitly mention THIS cluster (by number/name)?
2. If a quote discusses a broader set of genes, is it correctly attributed?
3. Would the functional_description still be accurate if you removed the quotes — 
   or is it entirely dependent on them?
If a quote doesn't match this specific cluster, remove it and reconsider the description.
```

This is not a separate LLM call — just additional prompt instructions that encourage chain-of-thought verification.

### Phase 3: Multi-organism clusters (Issue 7)

Bernstein 2017 is the only paper with joint clusters (A-D/E-H containing genes from both T. elongatus and M. ruber). Options:
- A) Add `cluster_mapping` to paperconfig: `{csv_cluster_0: "paper_cluster_A", ...}`
- B) Restructure as joint organism entries (one `gene_clusters` entry per paper cluster, not per organism)
- C) New `multi_organism_clustering` type in schema with links to multiple organisms

Needs design exploration — affects KG schema, adapter, and how downstream agents query clusters.

### Phase 4: Full re-extract baseline

After Phase 1-2 changes, re-extract all 15 entries with `--force` and compare against current baseline. Track improvement per paper.

## Current Prompt State

The prompt in `multiomics_kg/extraction/cluster/extract.py` has accumulated many rules through iteration. Before Phase 1, it should be refactored — the rules are getting long and some may conflict. The per-type template approach naturally organizes this: shared rules + type-specific rules.

## Files Changed Today

- `multiomics_kg/extraction/cluster/extract.py` — Pydantic schema, prompt, context builder, verification, report
- `tests/test_extraction_utils.py` — updated fixtures
- `tests/test_extract.py` — replaced duplicate ID test with filler detection test
- `docs/extraction_issues_log.md` — Issues 7-11
- `docs/superpowers/specs/2026-04-05-extraction-quality-iteration-design.md` — original spec
- `docs/superpowers/plans/2026-04-05-extraction-quality-iteration.md` — implementation plan
- `data/*/papers_and_supp/*/cluster_extractions/*.json` — re-extracted data
- `data/cluster_extraction_report.md` — baseline report
