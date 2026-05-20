---
name: <tool>-run
description: <one sentence: what the tool predicts + per-strain output location + Phase 1 marker. Inform from Step 0 Q1 + Q4>
argument-hint: "[--strains <name> ... | --force | --limit N | --threads N | --prepare-<external> | --refresh-<external>]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(docker *), Bash(<external> *)
---

# <Tool>-Run Skill

<one-paragraph overview: what the tool predicts, what input it consumes,
which Phase this skill covers, what the produced artifact looks like. Pull
straight from Step 0 Q1 + Q2.>

## Quick Start

```bash
# Run all genome strains (skip already-done)
uv run python .claude/skills/<tool>-run/run_<tool>.py

# Run one or more strains by name
uv run python .claude/skills/<tool>-run/run_<tool>.py --strains <NEW_STRAIN>
uv run python .claude/skills/<tool>-run/run_<tool>.py --strains MED4 MIT9313

# Force re-run even if calls.json exists
uv run python .claude/skills/<tool>-run/run_<tool>.py --force

# Smoke test: first 100 records of one strain (output is gitignored)
uv run python .claude/skills/<tool>-run/run_<tool>.py --strains MED4 --limit 100
```

## One-Time Setup

<placeholder — filled in Step 2 while installing.

Flavor A: `--prepare-image` command + image size + digest verification + `.env` vars.
Flavor B: `--refresh-<external>` command + reference DB size + URL + `<TOOL>_DATA_DIR` + which system CLI tool must be on PATH and how to install it.
Flavor C: license URL + tarball placement + `bash build.sh` + auxiliary DBs under `~/tools/<TOOL>/`; point at sibling INSTALL.md for verbose detail.

State explicitly:
- Disk footprint (image size, DB size)
- `.env` entry and whether it's required or optional
- Any portability gotchas
- Verification command (`docker images | grep <tool>`, `ls ~/tools/<TOOL>/DB/`)
>

## What It Does

<placeholder — Step 3.

Plain-prose walkthrough of the orchestrator: loads genome rows → invokes
external tool → parses → writes JSON → status table. Cite the exact docker
command / CLI invocation it builds. ~10–20 lines.>

## Output Schema

<placeholder — Step 3.

Field-by-field for `<strain>.<tool>.calls.json` + `<strain>.<tool>.skill_summary.json`.
Be explicit about nulls, sentinels, multi-valued fields. Document the
per-strain QC fields here so future users know what to look at.>

## QC

### Per-strain QC — what to inspect in `<strain>.<tool>.skill_summary.json`

<placeholder — Step 3 (initial); Step 5 (thresholds calibrated to real ranges).

Document each QC field with its expected range. Example template:

- `sentinel_rate < 0.40` — anything higher means the tool is failing to classify most proteins; check input encoding and external-tool version.
- `parse_failures == 0` — non-zero means raw output drifted from what the parser expects.
- `calls_made == input_proteins` — anything else means the runner skipped records; investigate.
>

### Cross-strain QC narrative

<placeholder — Step 5 (after batch).

Sanity benchmarks observed across strains. Example:

> OuterMembrane% ranged 0.5–1.5% in Pro strains, 2.5–4% in heterotrophs —
> matches the lifestyle prediction. Cross-tool sanity: of N proteins called
> `OuterMembrane` by psortb, M overlap with tcdb-diamond's outer-membrane
> TC families.
>

### Spot checks

<placeholder — Step 3 (table drafted); Step 4 (MED4 verified); Step 5 (all strains verified).

3–5 known-correct proteins per strain with expected calls + a `jq` one-liner. Example template:

| Strain | Protein ID | Expected | Why this is the ground truth |
|---|---|---|---|
| MED4 | WP_011131900.1 | OuterMembrane (score ≥ 9) | Pal lipoprotein, canonical OM marker |
| MIT9313 | WP_011131123.1 | CytoplasmicMembrane | UniProt-annotated membrane protein |
| MIT1002 | WP_NNNNN.1 | Extracellular | Secreted protease, signal peptide present |

```bash
jq '."WP_011131900.1".localization' \
  cache/data/Prochlorococcus/genomes/MED4/<tool>/MED4.<tool>.calls.json
# Expected: "OuterMembrane"
```
>

## Observed batch results

<placeholder — Step 5.

Cross-strain distribution table + per-strain wallclock range. Canonical shape:

## Observed batch results (32-strain run, YYYY-MM-DD)

N proteins classified. Cross-strain distribution:

| <Primary value> | Count | % |
|---|---:|---:|
| ... | ... | ... |

Per-strain wallclock: **X min** (SS120, N proteins) → **Y min** (KT2440, M proteins).
>

## Phase 2 (Future)

<placeholder — Step 3.

One paragraph + link to the design spec under `docs/superpowers/specs/`.
Be explicit that Phase 1 artifacts sit in the strain cache for inspection
and are NOT yet wired into `gene_annotations_merged.json` or any KG adapter.>

## Workflow When Invoked

<placeholder — Step 4.

Numbered checklist of what a user does when invoking the skill: verify
prerequisite (`docker --version`, `diamond --version`), run install if
needed, run the batch, inspect the status table, verify spot checks,
dig into FAILED rows in `logs/<tool>_<strain>.log`. ~5 lines.>
