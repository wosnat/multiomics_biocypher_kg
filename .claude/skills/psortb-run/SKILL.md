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

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`, filters to
   `organism_type='genome_strain'` (the 25 strains in the KG). Reference
   proteome match rows (Marinobacter, Alteromonas MarRef v6) have no
   `data_dir` and are skipped.
2. For each strain, runs PSORTb against `<data_dir>/protein.faa` via Docker:
   `docker run --rm -v <data_dir>:/tmp/psortb brinkmanlab/psortb_commandline:1.0.2 /usr/local/psortb/bin/psort --negative --output terse /tmp/psortb/protein.faa`.
3. Captures raw terse-format stdout to `<data_dir>/psortb/<strain>.psortb.terse.tsv`
   and container stderr to `logs/psortb_<strain>.log` (auto-gitignored via `*.log`).
4. Parses the 3-column terse TSV (`SeqID`, `Localization`, `Score`),
   splits multi-localization values (`Cytoplasmic/CytoplasmicMembrane`,
   `9.83/9.71`) into primary + secondary fields.
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
- When PSORTb returns a multi-localization (e.g.
  `Cytoplasmic/CytoplasmicMembrane` with two ≥ 7.5 scores), the
  highest-scoring class becomes `localization` + `score`, the second goes
  into `secondary_localization` + `secondary_score`, and
  `is_multi_localized=true`.

## Phase 2 (Future)

Phase 1 artifacts sit in the strain cache for inspection. Phase 2 (separate
spec) will integrate them into `gene_annotations_merged.json` and the KG.
See [docs/superpowers/specs/2026-05-10-psortb-localization-design.md §6](../../../docs/superpowers/specs/2026-05-10-psortb-localization-design.md).

## Workflow When Invoked

1. Verify `docker` is on PATH (`docker --version`).
2. If first run on this host, prepare the image: `uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image`.
3. Run the skill: `uv run python .claude/skills/psortb-run/run_psortb.py`.
   Each strain takes ~5-30 min depending on proteome size; the batch runs
   sequentially.
4. Review the status table — note any FAILED or MISSING_INPUT strains.
5. Inspect `<data_dir>/psortb/<strain>.psortb.calls.json` for spot checks.
6. For batch failures, check `logs/psortb_<strain>.log`.
