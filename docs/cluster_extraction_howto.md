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
- Duplicate `id` fields within an analysis
- Locus tags in descriptions (PMM*, PMT*, etc.)
- Empty direction for clusters with non-empty descriptions

### Iteration workflow
```bash
# 1. Generate report
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify

# 2. Read report, identify issues
cat data/cluster_extraction_report.md

# 3. Tweak prompt (in multiomics_kg/extraction/cluster/extract.py, DEVELOPER_MSG_TEMPLATE)

# 4. Re-extract specific entries
uv run python -m multiomics_kg.extraction.cluster.extract --entry X --force

# 5. Regenerate report, diff to see changes
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify
git diff data/cluster_extraction_report.md

# 6. Repeat until satisfied
```

## 3. Data Layout

```
data/<paper_dir>/
  paperconfig.yaml
  cluster_extractions/               # extraction results (one JSON + MD per entry)
    {entry_key}.json                 # structured extraction data (clean format)
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
    "extracted_at": "2026-04-05T13:25:22",
    "input_tokens": 62080,
    "output_tokens": 3141
  },
  "clusters": {
    "1": {
      "id": "mit9313_up_transport_binding",
      "name": "MIT9313 cluster 1 (up, transport and binding)",
      "functional_description": "...",
      "behavioral_description": "...",
      "peak_time_hours": 6.0,
      "period_hours": null,
      "direction": "up",
      "enrichment_category": "transport and binding",
      "enrichment_pvalue": 0.04,
      "enrichment_significant": true,
      "self_assessment": "high",
      "assessment_notes": "...",
      "confidence_notes": "...",
      "supporting_quotes": [{"quote": "...", "location": "..."}]
    }
  }
}
```

## 4. Architecture

The extraction pipeline has two stages:

| Stage | What | LLM? | Output |
|---|---|---|---|
| **0: Table** | CSV parsing + enrichment parsers | No | Per-cluster gene counts, enrichment data |
| **1: Extract** | Full PDF + all cluster summaries → structured descriptions | Yes (per analysis) | `cluster_extractions/{entry_key}.json` |

**One API call per `gene_clusters` entry** — the model receives the full PDF via `input_file` and extracts all clusters in the analysis simultaneously. This enables cross-cluster disambiguation (e.g., "cluster 6 is down while cluster 1 is up").

### Key modules

| File | Responsibility |
|---|---|
| `multiomics_kg/extraction/cluster/extract.py` | LLM calls, prompts, report generation, CLI |
| `multiomics_kg/extraction/cluster/extraction_utils.py` | File I/O, cluster key matching, CSV loading |
| `multiomics_kg/extraction/cluster/table.py` | Stage 0 CSV parsing + enrichment parsers |
| `multiomics_kg/adapters/cluster_adapter.py` | Reads extraction JSONs via `extraction_utils`, emits KG nodes |

## 5. Adding New Papers

1. Add `type: gene_clusters` entry to the paper's `paperconfig.yaml` (see CLAUDE.md for format)
2. Run extraction:
   ```bash
   uv run python -m multiomics_kg.extraction.cluster.extract --paper "Author Year"
   ```
3. Review: `uv run python -m multiomics_kg.extraction.cluster.extract --report --verify`
4. Rebuild KG — extracted descriptions automatically appear on `GeneCluster` nodes

## 6. Troubleshooting

### Model TPM limit exceeded
The PDF token count may exceed your account's TPM limit. Use `gpt-4.1-mini` (default) which has higher limits, or request a TPM increase for `gpt-4.1`.

### Cluster key mismatches
If the model renumbers or skips clusters, the matching logic tries: name regex, id suffix, case-insensitive, positional fallback. Check warnings in the extraction log. For multi-organism papers (e.g., Tolonen), use `extract_paper()` instead of `extract_analysis()` — see the Tolonen one-off in the git history.

### Self-assessment: low
Clusters with `self_assessment: "low"` typically mean the paper doesn't discuss them. These are correct — the adapter still creates the GeneCluster node (from CSV data) but without populated descriptions.

## 7. Current State

- **15 analyses** across 8 papers, **115 clusters** total
- **82 clusters** with populated descriptions
- **33 clusters** correctly empty (paper doesn't discuss them)
- **5 warnings**: 1 locus tag, 4 duplicate IDs on undescribed Zinser diel clusters
- Ground truth validated: Tolonen MIT9313 clusters 1, 6, 7 all correct
