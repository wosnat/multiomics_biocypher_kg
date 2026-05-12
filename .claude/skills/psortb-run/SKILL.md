---
name: psortb-run
description: Run PSORTb v3.0.3 (Gram-negative bacterial subcellular-localization predictor) via the brinkmanlab Docker image against each strain's protein.faa, emitting per-protein localization calls. Phase 1 — produces inspectable `<strain>.psortb.calls.json` artifacts; KG integration deferred to Phase 2.
argument-hint: [--strain <name> | --force | --prepare-image | --refresh-image | --limit N]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(docker *)
---

# PSORTb-Run Skill

Per-strain PSORTb v3.0.3 subcellular-localization prediction via the official
`brinkmanlab/psortb_commandline:1.0.2` Docker image. All KG strains are
Gram-negative, so `--negative` is hardcoded.

## Quick Start

```bash
# One-time: pull the ~2 GB Docker image
uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image

# Run all genome strains (skip already-done)
uv run python .claude/skills/psortb-run/run_psortb.py

# Run a single strain
uv run python .claude/skills/psortb-run/run_psortb.py --strain MED4

# Force re-run even if calls.json exists
uv run python .claude/skills/psortb-run/run_psortb.py --force

# Re-pull the image (e.g. after upstream release)
uv run python .claude/skills/psortb-run/run_psortb.py --refresh-image

# Smoke test: first 100 proteins of one strain (~1-2 min)
uv run python .claude/skills/psortb-run/run_psortb.py --strain MED4 --limit 100
```

## One-Time Setup

1. Install Docker and confirm it is on PATH: `docker --version`.
2. Pull the PSORTb image once: `uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image`.
   - Downloads `brinkmanlab/psortb_commandline:1.0.2` (~2 GB) into Docker's
     daemon storage, which is shared across all projects on the host.
   - Prints the image digest, then exits.
3. Normal runs refuse to start until the image is present locally — this
   avoids accidentally pulling 2 GB inside a long batch.

No `.env` configuration is required.

## What It Does

1. Loads genome rows via `multiomics_kg.download.utils.cli.load_genome_rows`
   (every row in `cyanobacteria_genomes.csv` with both `strain_name` and
   `data_dir`). Currently 32 rows: 30 `genome_strain` + 2
   `reference_proteome_match` (HP15, Alt_MarRef). Rows lacking `protein.faa`
   are skipped at run time with a `MISSING_INPUT` status.
2. For each strain, runs PSORTb against `<data_dir>/protein.faa` via Docker:
   ```
   docker run --rm --user $(id -u):$(id -g) \
     -v <data_dir>:/tmp/results \
     brinkmanlab/psortb_commandline:1.0.2 \
     /usr/local/psortb/bin/psort --negative --output terse -i /tmp/results/protein.faa
   ```
   Notes on the invocation:
   - The brinkmanlab `psort` Perl wrapper hardcodes its output directory to
     `/tmp/results/` inside the container and writes a timestamped file
     `<YYYYMMDDHHMMSS>_psortb_gramneg.txt` there — not to stdout. The skill
     mounts the strain's data_dir at `/tmp/results` so the file lands on the
     host, then renames it to `<strain>.psortb.terse.tsv`.
   - `-i` is mandatory; the wrapper does not accept a positional file argument.
   - `--user` makes produced files host-owned (no `sudo` needed for cleanup).
3. Renames the timestamped file to `<data_dir>/psortb/<strain>.psortb.terse.tsv`
   and writes container stdout/stderr to `logs/psortb_<strain>.log`
   (auto-gitignored via `*.log`).
4. Parses the 3-column terse TSV (`SeqID`, `Localization`, `Score`).
5. Writes `<data_dir>/psortb/<strain>.psortb.calls.json` (per-protein records,
   keyed by WP_ accession) and `<strain>.psortb.skill_summary.json`
   (per-strain stats including image digest, wallclock time, localization
   distribution).
6. Prints a status table to stdout — one row per strain.

## Output Schema (`<strain>.psortb.calls.json`)

Keyed by NCBI protein_id (WP_ accession):

```json
{
  "WP_011131900.1": {
    "localization": "OuterMembrane",
    "score": 9.97,
    "secondary_localization": null,
    "secondary_score": null,
    "is_multi_localized": false,
    "is_unknown": false
  }
}
```

Field semantics:
- `localization` ∈ {`Cytoplasmic`, `CytoplasmicMembrane`, `Periplasmic`,
  `OuterMembrane`, `Extracellular`, `Unknown`}.
- `score` ∈ [0, 10]; PSORTb assigns `Unknown` when no class scores ≥ 7.5
  (`is_unknown=true`).
- Multi-localization handling: the parser splits slash-separated values into
  primary + secondary fields. **However, with `--output terse` PSORTb emits a
  single `Final_Localization` per protein even when two classifiers both
  score ≥ 7.5**, so `is_multi_localized` is effectively always `false` in
  Phase 1. To expose true multi-localization, re-run with `--output long`
  and parse the per-class score block (deferred to Phase 2 if needed).

## Observed batch results (32-strain run, 2026-05-11)

96,418 proteins classified across all 32 rows. Cross-strain distribution:

| Localization | Count | % |
|---|---:|---:|
| Cytoplasmic | 40,415 | 41.9% |
| CytoplasmicMembrane | 21,231 | 22.0% |
| Periplasmic | 1,647 | 1.7% |
| OuterMembrane | 1,742 | 1.8% |
| Extracellular | 928 | 1.0% |
| Unknown | 30,455 | 31.6% |

Per-strain wallclock: **8.7 min** (SS120, 1858 proteins) → **30 min**
(KT2440, 5452 proteins). Total batch: ~8 h sequential. PSORTb startup
overhead dominates small-strain time; throughput scales sublinearly with
proteome size. OuterMembrane density tracks lifestyle — Prochlorococcus
~0.7% vs heterotrophs ~3% (richer OMP repertoire for nutrient scavenging).

## Phase 2 (Future)

Phase 1 artifacts sit in the strain cache for inspection. Phase 2 (separate
spec) will integrate them into `gene_annotations_merged.json` and the KG.
See [docs/superpowers/specs/2026-05-10-psortb-localization-design.md §6](../../../docs/superpowers/specs/2026-05-10-psortb-localization-design.md).

## Workflow When Invoked

1. Verify `docker` is on PATH (`docker --version`).
2. If first run on this host, prepare the image: `uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image`.
3. Run the skill: `uv run python .claude/skills/psortb-run/run_psortb.py`.
   Each strain takes ~9-30 min depending on proteome size; the batch runs
   sequentially (~8 h for all 32 strains). Use `nohup ... &` for unattended
   batch runs since stdout is line-buffered until each strain finishes.
4. Review the status table — note any FAILED or MISSING_INPUT strains.
5. Inspect `<data_dir>/psortb/<strain>.psortb.calls.json` for spot checks.
6. For batch failures, check `logs/psortb_<strain>.log`.
