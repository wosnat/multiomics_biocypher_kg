# Cluster Extraction — How-To Guide

## Prerequisites

- OpenAI API key in `.env` (`OPENAI_API_KEY=...`)
- `uv sync` to install dependencies
- Papers with `type: gene_clusters` entries in their `paperconfig.yaml`

## Quick Start

```bash
# See what would run (no API calls)
uv run python -m multiomics_kg.extraction.cluster.extract --dry-run

# Extract all entries
uv run python -m multiomics_kg.extraction.cluster.extract

# Generate quality report
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify
```

## 1. Run Extraction

### All entries
```bash
uv run python -m multiomics_kg.extraction.cluster.extract
```

### Single paper
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --paper "Tolonen 2006"
```

### Single entry
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --entry mit9313_kmeans_nstarvation
```

### Re-extract (overwrite existing)
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --entry mit9313_kmeans_nstarvation --force
```

### Cheaper bulk runs (50% discount via Flex processing)
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --flex
```

### Custom model
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --model gpt-4.1
```

Default model: `gpt-4.1-mini` (override via `CLUSTER_EXTRACTION_MODEL` env var).

## 2. Quality Review

### Generate report
```bash
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify
```

Creates `data/cluster_extraction_report.md` — a diff-friendly markdown file with all cluster descriptions and a Warnings section with programmatic quality checks.

### Quality checks performed
- Locus tags in descriptions (PMM*, PMT*, SY28_*, etc.)
- Low-confidence clusters with non-N/A descriptions (filler detection)
- Filler phrases ("likely", "possibly", "not explicitly described")
- Near-identical descriptions within an analysis
- No longer checked: empty direction (replaced by free-text `expression_dynamics`)

### Iteration workflow
```bash
# 1. Generate report
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify

# 2. Read report, identify issues
cat data/cluster_extraction_report.md

# 3. Tweak prompt (in extract.py: SHARED_RULES, TYPE_RULES, or SELF_VERIFICATION)

# 4. Re-extract specific entries
uv run python -m multiomics_kg.extraction.cluster.extract --entry X --force

# 5. Regenerate report, diff to see changes
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify
git diff data/cluster_extraction_report.md

# 6. Repeat until satisfied
```

## 3. Cluster Types

The extraction uses per-type prompt templates. Each `gene_clusters` entry must have a `cluster_type`:

| Type | Description | Papers |
|------|-------------|--------|
| `time_course` | Treatment/stress time-series | Tolonen 2006, Lindell 2007 |
| `diel` | 24h diel cycling | Zinser 2009, Coe 2024 |
| `condition_comparison` | Across discrete conditions | Alonso-Saez 2023, Wang 2014, Bernstein 2017 |
| `classification` | Categorizing genes by condition/periodicity | Biller 2018 |

## 4. Paperconfig Hints

### Publication-level `extraction:` section (optional)

```yaml
extraction:
  scope: paper              # "paper" (one LLM call for all entries) | "analysis" (default, one per entry)
  additional_pdfs:
    - data/.../supplementary.pdf   # uploaded alongside main PDF
```

### Per-entry fields on `gene_clusters` entries

- `figure_hint` (str) — tells the model which figures/tables to focus on (e.g., `"Figure 6, Table S3"`)
- `extraction_notes` (str) — free-text guidance passed into the prompt context (e.g., `"Paper discusses clusters jointly across organisms"`)

## 5. Extraction Fields

### Stored on GeneCluster nodes (via adapter)
- `functional_description` — gene identity and pathway membership (2-3 sentences max)
- `temporal_pattern` — observed expression dynamics (1-2 sentences)
- `expression_dynamics` — short label for expression response timing (e.g., "early transient", "peaks at dawn")

### Extraction-only scaffolding (not in KG)
- `expression_dynamics` — short free-text label (e.g., "early transient", "peaks at dawn", "periodic in L:D only")
- `enrichment_category`, `enrichment_pvalue`, `enrichment_significant`
- `supporting_quotes`, `source_figures`
- `self_assessment`, `assessment_notes`, `confidence_notes`

### Sentinel value
`"N/A"` — used when the paper doesn't discuss a cluster's function or behavior. Never filler text.

## 6. Data Layout

```
data/<paper_dir>/
  paperconfig.yaml
  cluster_extractions/               # extraction results (one JSON + MD per entry)
    {entry_key}.json                 # structured extraction data
    {entry_key}.md                   # human-readable summary

data/
  cluster_extraction_report.md       # full report across all papers (for review + diff)
```

### JSON format

```json
{
  "metadata": {
    "paper": "Tolonen 2006",
    "doi": "10.1038/msb4100087",
    "organism": "Prochlorococcus MIT9313",
    "entry_key": "mit9313_kmeans_nstarvation",
    "model": "gpt-4.1-mini",
    "extracted_at": "2026-04-06T10:08:07",
    "input_tokens": 62080,
    "output_tokens": 3141
  },
  "clusters": {
    "1": {
      "id": "mit9313_early_transport",
      "name": "MIT9313 cluster 1 (transport and binding)",
      "functional_description": "Enriched for transport and binding (p=0.04). Contains urtA and the nitrite permease.",
      "temporal_pattern": "Most rapidly upregulated cluster, responding within the first hours.",
      "expression_dynamics": "early rapid upregulation",
      "enrichment_category": "transport and binding",
      "enrichment_pvalue": 0.04,
      "enrichment_significant": true,
      "self_assessment": "high",
      "assessment_notes": "",
      "confidence_notes": "",
      "supporting_quotes": [{"quote": "...", "location": "Page 3"}],
      "source_figures": ["Figure 2"]
    }
  }
}
```

## 7. Architecture

The extraction uses **per-type prompt templates**:

```
SHARED_RULES          — rules for all cluster types
TYPE_RULES[type]      — type-specific field guidance + few-shot examples
SELF_VERIFICATION     — quote/cluster attribution self-check
extraction_notes      — per-entry free-text from paperconfig
cluster list          — cluster keys + gene counts from CSV
```

Assembled by `build_prompt(table_config, cluster_summaries)` into a single developer message.

**One API call per `gene_clusters` entry** (or per paper if `extraction.scope: paper`). The model receives the full PDF via OpenAI Files API and extracts all clusters simultaneously.

### Key modules

| File | Responsibility |
|---|---|
| `multiomics_kg/extraction/cluster/extract.py` | Prompt templates, LLM calls, report generation, CLI |
| `multiomics_kg/extraction/cluster/extraction_utils.py` | File I/O, cluster key matching, CSV loading, paperconfig discovery |
| `multiomics_kg/adapters/cluster_adapter.py` | Reads extraction JSONs, emits KG nodes/edges |
| `.claude/skills/paperconfig/validate_paperconfig.py` | Validates cluster_type vocabulary and extraction section |

## 8. Adding New Papers

1. Add `type: gene_clusters` entry to the paper's `paperconfig.yaml` (see CLAUDE.md for format)
2. Set `cluster_type` to one of: `time_course`, `diel`, `condition_comparison`, `classification`
3. Optionally add `figure_hint`, `extraction_notes`, and `extraction:` section
4. Run extraction:
   ```bash
   uv run python -m multiomics_kg.extraction.cluster.extract --paper "Author Year"
   ```
5. Review: `uv run python -m multiomics_kg.extraction.cluster.extract --report --verify`
6. For complex papers: manually review JSON against paper, fix misattributions
7. Rebuild KG — extracted descriptions automatically appear on `GeneCluster` nodes

## 9. Troubleshooting

### Model TPM limit exceeded
The PDF token count may exceed your account's TPM limit. Use `gpt-4.1-mini` (default) which has higher limits, or request a TPM increase for `gpt-4.1`.

### Cluster key mismatches
If the model renumbers or skips clusters, the matching logic tries: name regex, id suffix, case-insensitive, positional fallback. Check warnings in the extraction log.

### Cross-organism gene conflation
Papers discussing joint multi-organism clusters (e.g., Bernstein 2017) are prone to misattributing genes between organisms. Use `extraction_notes` to clarify which organism's genes are in scope, and manually verify the results.

### Self-assessment: low
Clusters with `self_assessment: "low"` typically mean the paper doesn't discuss them. These are correct — the adapter still creates the GeneCluster node (from CSV data) but with N/A descriptions.

## 10. Current State (2026-04-06)

- **15 analyses** across 8 papers, **110 clusters** total
- **~60 clusters** with populated descriptions (manually reviewed against papers)
- **~50 clusters** correctly N/A (paper doesn't discuss them)
- **12 warnings** — all intentional (low-confidence clusters with temporal patterns but N/A functional)
- Manually corrected: Tolonen cross-strain genes, Bernstein per-organism attributions, Biller organism confusion, Zinser timing, Alonso-Saez gene assignments, Wang spurious cluster removed
