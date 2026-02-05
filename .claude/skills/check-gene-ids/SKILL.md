---
name: check-gene-ids
description: Validate gene ID matching between paper CSV data and gene nodes in the knowledge graph. Reports per-paper match rates, categorizes mismatches, and suggests fix strategies. Use to diagnose dangling edges.
argument-hint: [paper-name or "all"]
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *)
---

# Gene ID Validation Skill

Check whether gene IDs used in paper supplementary tables match existing gene nodes in the knowledge graph.

## Quick Start

```bash
# Check all papers against latest biocypher output
uv run python .claude/skills/check-gene-ids/check_gene_ids.py

# Check a single paper
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"

# Specify biocypher output dir
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --biocypher-dir biocypher-out/20260205200505
```

## What It Checks

For each publication and each supplementary table/analysis:

1. Reads the CSV file using `name_col` and `skip_rows` from the paperconfig
2. Takes a sample of gene IDs and checks if `ncbigene:<id>` matches any Gene node `:ID`
3. If not matched, determines **why** and suggests a **fix strategy**

Skips tRNA/ncRNA/rRNA IDs automatically (these are intentionally not in gene nodes).

## Fix Strategies

| Strategy | Meaning | Action |
|----------|---------|--------|
| **MATCH** | IDs match gene nodes | No action needed |
| **LOAD_ORGANISM** | Organism genome not loaded | Add genome to `cyanobacteria_genomes.csv` |
| **CHANGE_NAME_COL** | Another CSV column has matching IDs | Change `name_col` in paperconfig.yaml |
| **CREATE_MAPPING** | IDs exist in GFF/GBK or alt gene fields | Build ID mapping from genome files |
| **UNRELATED_IDS** | IDs don't match anything | Manual investigation needed |

## Workflow

When invoked with a paper name (e.g., `/check-gene-ids "Biller 2018"`):

1. Find the paperconfig.yaml in `data/Prochlorococcus/papers_and_supp/$ARGUMENTS/`
2. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/$ARGUMENTS/paperconfig.yaml"`
3. Review the output
4. For CHANGE_NAME_COL: read the CSV to verify the suggested column, then update paperconfig.yaml
5. For CREATE_MAPPING: examine the GFF field identified, plan a mapping approach
6. For LOAD_ORGANISM: note the organism for future genome loading

When invoked with "all" or no arguments:

1. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py`
2. Review the summary table at the end
3. Prioritize fixes: CHANGE_NAME_COL (easy) > CREATE_MAPPING (moderate) > LOAD_ORGANISM (needs genome data)
