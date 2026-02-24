---
name: eggnog-run
description: Run eggNOG-mapper on genome protein FASTA files to generate functional annotations (COG, GO, KEGG, EC, PFAMs). Reads genomes from cyanobacteria_genomes.csv, skips already-annotated strains unless --force is given. Results are written to cache/<organism>/genomes/<strain>/eggnog/.
argument-hint: [--strain <name> | --force | --cpu <n>]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(python *)
---

# eggNOG-mapper Run Skill

Run eggNOG-mapper on all (or selected) genome protein FASTA files.

## Quick Start

```bash
# Run all genomes (skip already-annotated)
uv run python .claude/skills/eggnog-run/run_eggnog.py

# Run a single strain
uv run python .claude/skills/eggnog-run/run_eggnog.py --strain MED4

# Force re-run even if output exists
uv run python .claude/skills/eggnog-run/run_eggnog.py --force

# Use more CPUs
uv run python .claude/skills/eggnog-run/run_eggnog.py --cpu 8
```

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` (skips comment lines)
2. For each genome, checks for `<data_dir>/protein.faa`
3. Skips if `<data_dir>/eggnog/<strain>.emapper.annotations` already exists (unless `--force`)
4. Runs `emapper.py` with `--data_dir $EGGNOG_DATA_DIR` (from `.env`)
5. Prints a status table: skipped / ran / failed per strain

## Output Location

`cache/<organism>/genomes/<strain>/eggnog/<strain>.emapper.annotations`

## Dependencies

- `EGGNOG_DATA_DIR` set in `.env` (e.g. `~/tools/eggnog-mapper`)
- eggnog-mapper installed (`uv sync`)
- `protein.faa` downloaded (run `create_knowledge_graph.py` or genome download step first)

## Workflow When Invoked

When invoked (e.g. `/eggnog-run`):

1. Run: `uv run python .claude/skills/eggnog-run/run_eggnog.py`
2. Review the status table — note any FAILED strains
3. For failed strains, check that `protein.faa` exists and `EGGNOG_DATA_DIR` is correct
