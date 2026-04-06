# Cluster Extraction: Per-Type Prompts, Hints, and Cleanup

**Date:** 2026-04-06
**Status:** Approved
**Builds on:** `2026-04-05-extraction-quality-iteration-design.md`, `2026-04-05-extraction-quality-next-steps.md`

## Goal

Get cluster extraction quality to deployment-ready by combining three changes:
1. Per-type prompt templates that give type-appropriate guidance instead of one-size-fits-all
2. Paperconfig hints that give the LLM better context per paper
3. Prompt cleanup consolidating yesterday's accumulated rules

## Scope

**In scope:**
- Consolidate 7 cluster types → 4
- Per-type prompt templates with type-specific rules, examples, and field guidance
- New `extraction:` top-level section in paperconfig (scope, additional_pdfs)
- Per-entry hints (`figure_hint`, `extraction_notes`)
- Rename extraction fields (`direction` → `expression_dynamics`, `behavioral_description` → `temporal_pattern`)
- Fix Bernstein `cluster_col` to use paper cluster names
- Prompt self-verification block
- Re-extract all 15 entries

**Out of scope:**
- KG schema + adapter field rename (`behavioral_description` → `temporal_pattern`) — deferred, separate task
- Multi-organism cluster modeling (Issue 7 — Bernstein joint clusters need separate design)
- 4-stage pipeline (`pipeline.py`, review UI) — this work stays on the single-call `extract.py`
- Removing `peak_time_hours`/`period_hours` from KG schema (already removed from extraction)

## 1. Cluster Type Consolidation

7 types → 4:

| New type | Replaces | Papers |
|----------|----------|--------|
| `time_course` | `response_pattern` | Tolonen 2006, Lindell 2007 |
| `diel` | `diel_cycling`, `diel_expression_pattern` | Zinser 2009, Coe 2024 |
| `condition_comparison` | `expression_pattern`, `expression_level` | Alonso-Saez 2023, Wang 2014, Bernstein 2017 |
| `classification` | `periodicity_classification`, `expression_classification` | Biller 2018 (all 3 entries) |

Note: Wang 2014 is about constitutive expression levels (high/low/variable) rather than condition-driven changes. The free-text `expression_dynamics` field accommodates this — prompt guidance tells the model to use level labels like "high across all", "variable".

Update `cluster_type` in all 15 paperconfig entries. The adapter reads `cluster_type` for ClusteringAnalysis node properties, so the KG vocabulary changes too.

## 2. Field Renames

### Extraction schema (Pydantic)

| Old field | New field | Type | Purpose |
|-----------|-----------|------|---------|
| `direction` | `expression_dynamics` | `str` (free-text) | Short label for the expression pattern. Type-specific — e.g., "early transient", "peaks at dawn", "periodic in L:D only", "variable across conditions". Not stored in KG. |
| `behavioral_description` | `temporal_pattern` | `str` | 1-2 sentence prose description of expression dynamics. Will be stored on GeneCluster node after the deferred KG rename. |

`expression_dynamics` becomes free-text (was `Literal` enum). The per-type prompt section tells the LLM what kind of label to write.

All other extraction fields unchanged: `id`, `name`, `functional_description`, `enrichment_category`, `enrichment_pvalue`, `enrichment_significant`, `confidence_notes`, `supporting_quotes`, `source_figures`, `self_assessment`, `assessment_notes`.

### KG schema + adapter (DEFERRED)

The rename of `behavioral_description` → `temporal_pattern` in `schema_config.yaml` and `cluster_adapter.py` is deferred to a separate task. For now, the extraction writes `temporal_pattern` in its JSON, and the adapter continues reading `behavioral_description`. The adapter will need a mapping when the KG rename happens.

`expression_dynamics` is extraction-only scaffolding — not in KG schema.

### "Not discussed" sentinel

Change from `"Not discussed in paper."` → `"N/A"` in extraction:
- Prompt rules and few-shot examples
- Filler detection in `verify_quality()`

Adapter fallback value update deferred with the other adapter changes.

## 3. Per-Type Prompt Templates

### Architecture

Replace the single `DEVELOPER_MSG_TEMPLATE` with assembled prompts:

```
SHARED_RULES          — rules for all types
TYPE_RULES[type]      — type-specific field guidance + examples
SELF_VERIFICATION     — quote/cluster attribution check
extraction_notes      — per-entry free-text (from paperconfig)
cluster list          — cluster keys + gene counts
```

Assembled by `build_prompt(table_config, cluster_summaries)` → single developer message string.

### Shared rules (all types)

- `"N/A"` for undescribed clusters — no filler, no speculation
- No locus tags in descriptions (PMM*, PMT*, P9301_*, etc.)
- No treatment conditions in descriptions (lives on analysis node)
- No membership stats or cluster definition restatement
- Only cite what paper explicitly says about THIS specific cluster
- Max 3-5 named genes per cluster, paper names only
- Enrichment: only paper-reported with p-values
- `source_figures` and `supporting_quotes` required
- `id` format: `{organism_short}_{dynamics}_{theme}` in snake_case
- `name` format: `"{Organism} cluster {KEY} ({theme})"` — under 60 chars
- `self_assessment`: confidence level. `assessment_notes`: what you're uncertain about

### Type-specific rules

**`time_course`** — treatment/stress time-series (Tolonen, Lindell)
- `expression_dynamics`: response timing label — "early transient", "late sustained", "rapid then declining", "gradual increase"
- `temporal_pattern`: describe when change begins, how fast, whether it persists or reverses. Include timepoints from paper.
- `functional_description`: enrichment category + highlighted genes from paper
- Few-shot: Tolonen MED4 cluster 1 (up, transport — rapid early response) + cluster with "N/A"

**`diel`** — 24h cycling (Zinser, Coe)
- `expression_dynamics`: peak phase label — "peaks at dawn", "peaks at dusk", "peaks midday", "peaks at night"
- `temporal_pattern`: describe peak timing (hours), periodicity, phase relative to light/dark cycle
- `functional_description`: enrichment + highlighted genes
- Few-shot: Zinser cluster 1 (photosynthesis, peaks near dawn)

**`condition_comparison`** — across discrete conditions (Bernstein, Alonso-Saez, Wang)
- `expression_dynamics`: condition-response label — "up with light", "down at cold", "variable", "high across all"
- `temporal_pattern`: describe which conditions drive changes and direction/magnitude
- `functional_description`: enrichment + highlighted genes
- Few-shot: Alonso-Saez cluster C (cold stress response, upregulated at Tmin)

**`classification`** — categorizing genes by condition/periodicity (Biller)
- `expression_dynamics`: category label — "periodic in L:D only", "not periodic", "present in darkness"
- `temporal_pattern`: describe which conditions/categories the genes fall into
- `functional_description`: pathway associations if paper reports them
- Few-shot: Biller periodicity cluster (periodic in coculture L:D and darkness)

### Self-verification block (all types)

Added to end of shared rules:

```
Before outputting, verify for each cluster:
1. Does each supporting_quote explicitly mention THIS cluster (by number/name)?
2. If a quote discusses a broader set of genes, is it correctly attributed?
3. Would the functional_description still be accurate if you removed the quotes?
If a quote doesn't match this specific cluster, remove it and set description to "N/A".
```

## 4. Paperconfig Extraction Hints

### New top-level `extraction:` section

```yaml
publication:
  papername: ...
  papermainpdf: ...
  experiments: ...
  supplementary_materials: ...

extraction:
  scope: paper              # "paper" | "analysis" (default: "analysis")
  additional_pdfs:
    - data/.../supplementary.pdf
```

All fields optional. The section itself is optional.

- `scope: paper` — one LLM call with all `gene_clusters` entries combined (uses `extract_paper()`). Better for papers where the discussion covers multiple organisms/conditions together.
- `scope: analysis` (default) — one LLM call per `gene_clusters` entry (uses `extract_analysis()`).
- `additional_pdfs` — supplementary PDFs uploaded alongside the main paper PDF. All PDFs included as file inputs in the LLM call.

### Per-entry fields on `gene_clusters` entries

- `figure_hint` (str, existing) — already wired in `build_context_block()`. E.g., `"Figure 6, Table S3"`. Tells the model which figures/tables to focus on.
- `extraction_notes` (str, new) — free-text guidance passed into prompt context block. E.g., `"Paper discusses clusters A-D (light-responsive) and E-H (oxygen-responsive) as joint T. elongatus + M. ruber clusters. This entry covers only the T. elongatus genes."`

`extraction_notes` is appended to the context block after all other fields.

## 5. Paper-Specific Paperconfig Changes

### Bernstein 2017

All 4 entries:
- Change `cluster_col`: `ClustID_light` → `Clust name_light`, `ClustID_ox` → `Clust name_ox`
- Change `cluster_type`: `response_pattern` → `condition_comparison`
- Add `figure_hint: "Figure 6"` (or appropriate figure per entry)
- Add `extraction_notes` explaining the joint organism context

This changes cluster node IDs from numeric (0, 1, 2...) to letter (A, B, C, D). The adapter handles string keys already.

### Coe 2024

Both entries:
- Change `cluster_type`: `diel_expression_pattern` → `diel`
- Add to publication level:
  ```yaml
  extraction:
    additional_pdfs:
      - data/Prochlorococcus/papers_and_supp/coe 2024/coe_isme_supp_final_rev_final_ycae131.pdf
  ```

### Tolonen 2006

- Change `cluster_type`: `response_pattern` → `time_course`
- Add to publication level:
  ```yaml
  extraction:
    scope: paper
  ```

### Lindell 2007

- Change `cluster_type`: `response_pattern` → `time_course`

### Zinser 2009

- Change `cluster_type`: `diel_cycling` → `diel`

### Alonso-Saez 2023

- Change `cluster_type`: `expression_pattern` → `condition_comparison`

### Wang 2014

- Change `cluster_type`: `expression_level` → `condition_comparison`

### Biller 2018

All 3 entries:
- Change `cluster_type`: `periodicity_classification` / `expression_classification` → `classification`

## 6. Code Changes Summary

### `multiomics_kg/extraction/cluster/extract.py`

- Rename Pydantic fields: `direction` → `expression_dynamics` (change to `str`), `behavioral_description` → `temporal_pattern`
- Replace `DEVELOPER_MSG_TEMPLATE` + `CLUSTER_TYPE_GUIDANCE` with `SHARED_RULES`, `TYPE_RULES` dict, `SELF_VERIFICATION`
- New `build_prompt(table_config, cluster_summaries)` function assembling the full prompt
- Update `build_context_block()` to include `extraction_notes`
- Update `extract_analysis()` and `extract_paper()` to read `extraction:` config and handle `additional_pdfs` (upload multiple files, pass as list of file inputs)
- Update `main()` CLI to respect `extraction.scope` when deciding per-paper vs per-analysis calls
- Update `verify_quality()`: `"N/A"` sentinel, update filler checks
- Update `generate_report()`: use new field names
- Update few-shot examples per type

### Paperconfig files (8 papers, 15 entries)

- Update `cluster_type` values per Section 5
- Add `extraction:` sections where needed
- Add `figure_hint` and `extraction_notes` where needed
- Change Bernstein `cluster_col` values

### `multiomics_kg/extraction/cluster/extraction_utils.py`

- Update `match_cluster_keys()` if needed for new field names in parsed output

### `.claude/skills/paperconfig/validate_paperconfig.py`

- Update `VALID_CLUSTER_TYPES` to the 4 new types: `time_course`, `diel`, `condition_comparison`, `classification`
- Add `extraction_notes` to allowed `gene_clusters` entry fields
- Accept new top-level `extraction:` section (with `scope` and `additional_pdfs` keys)
- Run validation on all paperconfigs after changes

### Tests

- Update extraction test fixtures for new field names
- Update paperconfig validation test fixtures for new cluster types

## 7. Execution Order

1. Extraction field renames (Pydantic schema in extract.py)
2. Update paperconfig validator (new cluster types, extraction section, extraction_notes)
3. Paperconfig changes (cluster_type consolidation, Bernstein cluster_col, hints)
4. Run paperconfig validation on all files
5. Prompt refactor (shared rules + per-type templates + self-verification)
6. Update `extract.py` CLI (additional_pdfs, scope logic)
7. Update verification and report
8. Update tests
9. Run full test suite
10. Re-extract all 15 entries with `--force`
11. Generate report, review warnings
12. Iterate on prompts if needed

## 8. Success Criteria

- Warnings drop from 27 to <10
- Bernstein clusters get real descriptions (not "N/A" for everything) — at least clusters A-D should have content
- Coe 2024 clusters correctly say "N/A" (not filler)
- Biller 2018 periodicity clusters use appropriate classification language
- No locus tags in any description
- No filler phrases in any described cluster
- All 115 clusters extracted, no missing keys
