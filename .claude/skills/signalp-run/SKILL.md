---
name: signalp-run
description: Run SignalP 6.0 (via Docker) on genome protein FASTA files to predict signal peptides (SP, LIPO, TAT, PILIN). Reads genomes from cyanobacteria_genomes.csv, skips already-predicted strains unless --force is given. Results are written to cache/<organism>/genomes/<strain>/signalp/.
argument-hint: [--strain <name> | --force]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(python *), Bash(docker *), Bash(bash *)
---

# SignalP 6.0 Run Skill

Run SignalP 6.0 signal peptide prediction on all (or selected) genome protein
FASTA files. Uses a Docker container to handle SignalP's Python ≤3.10 /
PyTorch <2.0 requirements; image lives in Docker's daemon storage (not in
the repo).

| File | Purpose |
|------|---------|
| `SKILL.md` | This file — invocation guide |
| `INSTALL.md` | One-time setup (download tarball, build image) |
| `Dockerfile` | Image recipe; mirrors the [official install steps](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md) |
| `build.sh` | Stages `~/tools/signalp-6.0i/` as the build context and runs `docker build` |
| `run_signalp.py` | Batch driver: reads `cyanobacteria_genomes.csv`, runs SignalP per strain |

## Quick Start

```bash
# One-time setup — see INSTALL.md
bash .claude/skills/signalp-run/build.sh

# Run all genomes (skip already-predicted)
uv run python .claude/skills/signalp-run/run_signalp.py

# Run a single strain
uv run python .claude/skills/signalp-run/run_signalp.py --strain HOT1A3

# Force re-run even if output exists
uv run python .claude/skills/signalp-run/run_signalp.py --force
```

## What It Does

1. Checks that the `signalp6` Docker image exists (build with `bash .claude/skills/signalp-run/build.sh`)
2. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`
3. For each genome, checks for `<data_dir>/protein.faa`
4. Skips if `<data_dir>/signalp/prediction_results.txt` already exists (unless `--force`)
5. Runs `docker run --rm --user $(id -u):$(id -g) -v <genome_dir>:/data signalp6 --organism other --mode fast ...`
6. Prints a status table with prediction counts per strain

## Output Location

`cache/<organism>/genomes/<strain>/signalp/prediction_results.txt`

Tab-delimited columns: protein ID, prediction type (`SP`/`LIPO`/`TAT`/`TATLIPO`/`PILIN`/`OTHER`), per-type probabilities, predicted cleavage site. SignalP also writes `processed_entries.fasta`, `output.gff3`, `region_output.gff3`, and `output.json` to the same directory. See the [official output reference](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#output-interpretation).

## GPU Support

CPU-only on this host (no NVIDIA hardware). The `signalp6_convert_models gpu`
upgrade path for an NVIDIA box is documented in [INSTALL.md](./INSTALL.md#gpu-upgrade-path-nvidia-hosts-only).

## Prerequisites

- Docker installed and `signalp6` image built — see [INSTALL.md](./INSTALL.md)
- `protein.faa` downloaded per strain (run `bash scripts/prepare_data.sh` first)

## Workflow When Invoked

When invoked (e.g. `/signalp-run`):

1. Run: `uv run python .claude/skills/signalp-run/run_signalp.py`
2. Review the status table — note any `FAILED` strains
3. For failed strains, check that `protein.faa` exists and Docker is running
