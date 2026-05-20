---
name: add-a-tool
description: Scaffold a new per-strain bioinformatics tool integration (Phase 1) following the psortb-run / tcdb-diamond / eggnog-run / signalp-run template. Step 0 captures intent (what it predicts, input shape, install flavor, trigger phrasings), then six operational steps — install script (+ SKILL.md skeleton) → actually install (+ fill One-Time Setup) → runner with QC fields (+ fill What It Does / Output Schema / QC fields / spot-checks) → smoke test via the skill (+ validate Quick Start) → full batch via the skill (+ fill Observed batch results / cross-strain QC narrative) → register on new strains. The SKILL.md grows in lockstep with each step so the documentation is validated by the next action. Produces a `/tool-run` skill that loops `cyanobacteria_genomes.csv`, runs the external tool on each strain's `protein.faa` (or genome / GFF), writes inspectable `<strain>.<tool>.calls.json` artifacts committed under git, and ships a SKILL.md with single-strain invocation (run on new strain), fresh-install-on-new-machine instructions, and QC including spot checks. Use whenever the user wants to add a new prediction tool (SignalP, InterProScan, DeepLoc, AntiSMASH, dbCAN3, BUSCO, …), integrate a new external annotation source per-strain, or says "add a tool for X", "wire up SignalP", "we need per-protein <prediction> for all strains". Also use when generalising an ad-hoc one-off script into a reusable per-strain skill.
argument-hint: <tool-name> [--external "<tool description>"] [--phase 1|2]
user-invocable: true
allowed-tools: Read, Edit, Write, Grep, Glob, Bash(uv *), Bash(python *), Bash(mkdir *), Bash(git *), Bash(docker *), WebFetch
---

# Add a Tool

Per-strain external prediction tool → fully-scaffolded `/<tool>-run` skill,
following the as-built patterns from `psortb-run`, `tcdb-diamond`,
`eggnog-run`, and `signalp-run`. Phase-1 only — inspectable per-strain JSON
artifacts under git, no KG-side coupling until a separate Phase-2 spec lands.

When in doubt, read [`psortb-run/SKILL.md`](../psortb-run/SKILL.md),
[`tcdb-diamond/SKILL.md`](../tcdb-diamond/SKILL.md), and
[`signalp-run/SKILL.md`](../signalp-run/SKILL.md) and copy the patterns. This
skill encodes their lessons so the next tool doesn't drift.

## When to use

- User mentions adding any external bioinformatics tool that runs per-strain over `protein.faa` / `genome.fna` / `genomic.gff`.
- Concrete recent examples: PSORTb, TCDB-diamond, eggNOG-mapper, SignalP.
- Plausible near-future examples: InterProScan, AntiSMASH, DeepLoc, dbCAN3, BUSCO, GhostKOALA.

If the tool produces **graph-shape data** (new node type or new edge type
going directly into BioCypher), that's a Phase-2 adapter and outside this
skill's scope — write a spec under `docs/superpowers/specs/` and integrate via
the adapter pattern in `multiomics_kg/adapters/`.

## Why Phase 1 is its own step

Splitting per-strain prediction from KG integration is intentional:

1. **Inspection cycle.** Phase-1 artifacts are committed JSON the user can `jq`, `diff`, and code-review before committing to a schema. PSORTb's `is_multi_localized` field surfaced an interesting `--output terse` vs `--output long` nuance only after Phase-1 artifacts existed for all 32 strains.
2. **Wallclock.** Tool runs are slow (hours per batch). Decoupling lets the artifacts persist across many Phase-2 schema iterations without re-running.
3. **Cache locality.** Outputs live in `<data_dir>/<tool>/` next to `protein.faa` — same place every adapter and Phase-2 integration will look.

Phase 2 (deferred) integrates the calls.json into either
`gene_annotations_merged.json` (`prepare_data.sh` step 2 reads it) or a
dedicated KG adapter — separate design spec.

## Before you start

Every tool that's shipped has a design spec under
`docs/superpowers/specs/<YYYY-MM-DD>-<tool>-design.md` and an implementation
plan under `docs/superpowers/plans/<YYYY-MM-DD>-<tool>-skill.md`. They aren't
numbered steps below because they're scaffolding for the implementation, not
the operational workflow — but write them. The shortest sufficient version
of the spec lives in Step 0 below; the formalized markdown spec then expands
on those answers. The spec answers what the tool
predicts, the install footprint, the per-strain wallclock, the output schema
field-by-field, the tier/scoring policy (if applicable), and a Phase-2
sketch. The plan is a checkbox task list. Canonical examples:
[`2026-05-10-psortb-localization-design.md`](../../../docs/superpowers/specs/2026-05-10-psortb-localization-design.md),
[`2026-05-10-tcdb-diamond-augmentation-design.md`](../../../docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md),
[`2026-05-10-tcdb-diamond-skill.md`](../../../docs/superpowers/plans/2026-05-10-tcdb-diamond-skill.md).

---

## SKILL.md authoring map

The tool-run SKILL.md grows alongside each operational step rather than
being written in one block. Every step has a "Document this" sub-task that
fills (or refines) a section while the work is fresh. By step 5 the
SKILL.md is complete and battle-tested — every documented invocation has
been used, every QC threshold has been calibrated against real data.

| SKILL.md section | First drafted in | Refined / finalized in | Why this timing |
|---|---|---|---|
| Frontmatter | Step 1 | — | Tool name + flag surface decided up front |
| Quick Start | Step 1 | Step 4 (validated by smoke) | Documented invocations are the runner's contract |
| One-Time Setup | Step 2 | Step 2 | Capture every install command as you run it |
| What It Does | Step 3 | — | Written while implementing the runner |
| Output Schema | Step 3 | Step 4 (if smoke reveals tweaks) | Co-evolves with the parser |
| QC — per-strain field docs | Step 3 | Step 5 (thresholds calibrated to real ranges) | Document the fields you're emitting |
| QC — spot checks table | Step 3 | Step 4 (MED4 verified) → Step 5 (all strains verified) | Pick known proteins before any data lands |
| QC — cross-strain narrative | (placeholder Step 3) | Step 5 (after batch) | Needs cross-strain data to write |
| Observed batch results | (placeholder Step 3) | Step 5 | After batch completes |
| Phase 2 (Future) | Step 3 | — | Link to the design spec; deferred |
| Workflow When Invoked | Step 4 | — | Stable once Quick Start is validated |

Step 0 below produces the answers that drive every column. Don't skip it —
it's three minutes of writing that saves a half day of rework.

---

## Step 0 — Capture intent

Before any code, write down answers to these four questions. They become the
spec; the rest of the steps consume them. Three minutes of writing here
prevents an afternoon of rewriting later.

### Q1 — What does the tool predict?

What does one call from the tool say about one protein / gene / region?
Pick the shape that matches:

- **Per-protein call** — one record per protein, with a primary value + score. Examples: PSORTb (`localization` + score), SignalP (`signal_peptide_type` + probability), DeepLoc.
- **Per-protein list-of-candidates** — multiple ranked hits per protein. Example: tcdb-diamond emits one candidate per 3-part TC family.
- **Per-genome region / interval** — annotates positions on the contig. Examples: AntiSMASH BGC clusters, CRISPRCasFinder, prophage finders.
- **Per-gene metric** — numeric or categorical value attached to a gene. Examples: BUSCO completeness, GC content per gene, codon-usage scores.
- **Cross-protein relationship** — output is a graph or pairwise table. Examples: domain-domain interactions, paralog clustering.

This decision picks the `<strain>.<tool>.calls.json` schema shape (single
object keyed by WP_ vs nested candidates list vs intervals list).

### Q2 — What input does the tool need?

Pick the file(s) the tool consumes:

- `protein.faa` — most prediction tools (PSORTb, SignalP, eggNOG, tcdb-diamond)
- `genome.fna` — nucleotide-based predictors (AntiSMASH, prophage, CRISPR)
- `genomic.gff` — annotation-aware tools that need gene calls
- pre-computed annotation table (e.g., `gene_annotations_merged.json`) — meta-tools that consume other tools' output
- Other — call it out; you may need a custom loader

This decides what `run_<tool>.py` reads inside the per-strain loop.

### Q3 — Install flavor?

Pick one — the install script in Step 1 follows the matching pattern:

- **Flavor A** — Public Docker image (psortb-run). `docker pull` only; no build.
- **Flavor B** — Public CLI + downloadable reference DB (tcdb-diamond, eggnog-run). Tool on PATH; DB lands at `~/tools/<TOOL>/`.
- **Flavor C** — Build from source / academic license (signalp-run). `Dockerfile` + `build.sh` + `INSTALL.md` bundled with the runner.

Picking now decides the file layout (Dockerfile? `build.sh`? `INSTALL.md`?)
and the runner's setup flag (`--prepare-image` vs `--refresh-<external>` vs
`bash build.sh`).

### Q4 — What user phrasings should `/<tool>-run` trigger on?

List 3–5 realistic phrasings a user might type. These inform the SKILL.md
frontmatter `description:` field and the optional `/skill-creator` description
optimization in step 5.5. Examples for a hypothetical SignalP:

- "predict signal peptides for all strains"
- "run signalp on MED4"
- "where are the secreted proteins in MIT9313"
- "we need SP predictions for the new strain"

If the tool's name is obvious-matching (e.g., user always types the tool's
exact name), 1–2 phrasings is fine. If the tool serves a semantic need users
phrase many ways (subcellular location, transporter classification),
broader coverage matters.

### What you produce from Step 0

A four-bullet answer block, e.g.:

```
- Predicts: per-protein call (signal_peptide_type ∈ {SP, LIPO, TAT, PILIN, OTHER} + probability)
- Input: protein.faa
- Install flavor: C (build from source, academic license)
- Triggers: "signal peptides for X", "secreted proteins", "run signalp on Y", "SP predictions for new strain"
```

These four answers drive Steps 1–5: Q1 shapes the calls.json schema, Q2
shapes the runner's input handling, Q3 shapes the install script + One-Time
Setup section, Q4 shapes the frontmatter description.

Optionally formalize the answers into a `docs/superpowers/specs/<date>-<tool>-design.md`
spec — the four existing tools all have one — but Step 0's four bullets are
the minimum viable input to start Step 1.

---

## Step 1 — Install script + SKILL.md skeleton

### Work

The install script is what gets the external tool onto a fresh host. Pick the
flavor that matches the tool:

**Flavor A — Public Docker image (psortb-run pattern).** Tool ships as a
public Docker image (e.g., `brinkmanlab/psortb_commandline:1.0.2`). No build
needed; just `docker pull`. Bundle into the runner as a `--prepare-image`
flag that pulls + prints the digest + exits.

**Flavor B — Public CLI + downloadable reference DB (tcdb-diamond / eggnog-run pattern).** Tool itself is already on PATH (`diamond`, `emapper.py`). What needs installing is the reference DB. Bundle into the runner as `--refresh-<external>` / first-run-auto-download. DB lands under `~/tools/<TOOL>/` (outside the repo — see "Installation lives outside the repo" below).

**Flavor C — Build-from-source / academic-license tarball (signalp-run pattern).** Tool requires an academic license to download the source tarball, then a Dockerfile build to produce a local image. Bundle three files alongside the runner:

```
.claude/skills/<tool>-run/
├── SKILL.md
├── run_<tool>.py
├── Dockerfile       # builds local image
├── build.sh         # one-shot: docker build → tag as <tool>:local
└── INSTALL.md       # license URL + tarball placement instructions for the user
```

For all three flavors, the runner **refuses to start a normal run if the
install is missing** and prints the exact next command — PSORTb's pattern is
gold: the `--prepare-image` step prints the digest and exits, so the user
can't accidentally launch an 8-hour batch with no image. Don't auto-install
silently — installs are slow and license-sensitive; the user wants to see
them happen.

### Document this — write the SKILL.md skeleton

The skeleton lives at
[`assets/tool_run_skill_template.md`](assets/tool_run_skill_template.md) —
copy it into the new tool's directory and fill in the Step-0 answers plus
the flag surface you just decided:

```bash
mkdir -p .claude/skills/<tool>-run
cp .claude/skills/add-a-tool/assets/tool_run_skill_template.md \
   .claude/skills/<tool>-run/SKILL.md
```

Then edit the new SKILL.md:

- Replace every `<tool>` / `<TOOL>` / `<NEW_STRAIN>` literal.
- Fill the **frontmatter** `description:` from Step 0 Q1 + Q4 (what it predicts + trigger phrasings).
- Fill the **one-paragraph overview** from Step 0 Q1 + Q2.
- Fill the **Quick Start** block — the exact flags you'll implement (`--prepare-<external>` for Flavor A, `--refresh-<external>` for B, etc.). It's OK that the runner doesn't exist yet — Quick Start is your contract for what it'll accept.
- Leave the placeholder blocks (`<placeholder — Step N>`) for the later steps to fill in.

Commit the skeleton. You now have a stable target the rest of the workflow
fills out.

## Step 2 — Actually install it + fill One-Time Setup

### Work

Run the install on the current host. This is a real step, not a footnote —
the install is what proves the script works, and without it the rest of the
workflow blocks.

```bash
# Flavor A
uv run python .claude/skills/<tool>-run/run_<tool>.py --prepare-image

# Flavor B
uv run python .claude/skills/<tool>-run/run_<tool>.py --refresh-<external>

# Flavor C
bash .claude/skills/<tool>-run/build.sh
# Follow INSTALL.md for the license tarball + license-key placement.
```

Verify by checking what landed under `~/tools/<TOOL>/` (or `docker images |
grep <tool>`).

### Document this — fill in **One-Time Setup**

While you're running the install, write down every command and observation
in the SKILL.md's **One-Time Setup** section. A user who just cloned the
repo on a fresh machine should be able to follow this section and end up
with a working tool. Cover whichever flavor applies:

- **Flavor A**: `--prepare-image` command + image size + digest verification + which `.env` vars to set
- **Flavor B**: `--refresh-<external>` command + reference DB size + URL + `<TOOL>_DATA_DIR` env var + which system CLI tool must be on PATH (`diamond`, `emapper.py`, etc.) and how to install it
- **Flavor C**: license URL + tarball placement + `bash build.sh` + which auxiliary DBs (if any) go under `~/tools/<TOOL>/` — point to a sibling `INSTALL.md` for verbose detail

State explicitly:

- Disk footprint (image size, DB size)
- The `.env` entry (`<TOOL>_DATA_DIR=~/tools/<TOOL>`) and whether it's required or optional
- Any portability gotchas (extractTCDB.pl dash-vs-bash, PSORTb's `/tmp/results` hardcode, Ubuntu-specific package names, etc.)
- Verification command (`docker images | grep <tool>`, `ls ~/tools/<TOOL>/DB/`, etc.) so the user knows the install worked

Writing this **while installing** captures the gotchas you actually hit, not
the idealized version you'd write a week later.

## Step 3 — Runner script (with QC fields) + fill data-free SKILL.md sections

### Work — the runner

Write the runner at `.claude/skills/<tool>-run/run_<tool>.py`. It should:

1. Parse CLI flags (see "Canonical CLI surface" below) — match what your Step 1 Quick Start documented.
2. Iterate strains via `multiomics_kg.download.utils.cli.load_genome_rows`.
3. For each strain, invoke the external tool against the input from Step 0 Q2 (`<data_dir>/protein.faa`, `<data_dir>/genome.fna`, etc.); capture full stdout+stderr to `logs/<tool>_<strain>.log`.
4. Parse the tool's raw output, keyed by **the natural identifier for one record of what the tool predicts** (see "Choosing the calls.json key" below). For all four existing tools that's WP_ protein_id (NCBI accession) because they all consume `protein.faa`; tools consuming `genomic.gff` or `genome.fna` will pick a different key.
5. Write `<data_dir>/<tool>/<strain>.<tool>.calls.json` and `<strain>.<tool>.skill_summary.json` — with **per-strain QC fields baked into skill_summary** (see below).
6. Emit a status table to stdout — one row per strain.

### Choosing the calls.json key

Driven by Step 0 Q1 (prediction shape) and Q2 (input file):

| Q1 prediction shape | Q2 input | calls.json key | Example |
|---|---|---|---|
| Per-protein call (single primary value) | `protein.faa` | WP_ `protein_id` from gene_mapping.csv | psortb, signalp, eggnog |
| Per-protein list-of-candidates | `protein.faa` | WP_ `protein_id` → `{shared fields, calls: [...]}` | tcdb-diamond |
| Per-gene metric | `genomic.gff` or merged annotations | `locus_tag` | hypothetical codon-usage-per-gene tool |
| Per-genome region / interval | `genome.fna` | top-level list `[{contig, start, end, ...}, ...]` — no per-record key | AntiSMASH BGC clusters, prophage finders |
| Per-genome scalar | any | top-level scalar `{value, metadata}` | BUSCO completeness |
| Cross-protein relationship | `protein.faa` | composite key (e.g., `"WP_A|WP_B"`) or list of pairs | domain-domain interaction predictors |

WP_ is the right key for any tool consuming `protein.faa` because it's the
join key into `gene_mapping.csv` — works across strains regardless of locus
prefix differences. Tools consuming GFF/genome inputs pick the identifier
that naturally indexes the prediction shape; document the choice in the
Output Schema section so consumers know what they're joining on.

```python
from multiomics_kg.download.utils.cli import load_genome_rows

rows = load_genome_rows(strains=args.strains)
for row in rows:
    data_dir = Path(row["data_dir"])
    protein_faa = data_dir / "protein.faa"
    if not protein_faa.exists():
        # Emit MISSING_INPUT status row; continue.
        continue
    ...
```

**Split pure-Python parsing from the orchestrator.** Put parsing + scoring +
post-step logic in `multiomics_kg/utils/<tool>.py` so it's unit-testable
without subprocess or filesystem. Cover it with `tests/test_<tool>.py`.
tcdb-diamond does this cleanly; psortb-run inlined the parsing because the
post-step was trivial — judgement call.

### Shared I/O utils — write once, reused by runner + QC + Phase 2

Don't inline path resolution + `json.dump` / `json.load` in the runner. Put
them in **`multiomics_kg/utils/tool_calls_io.py`** — a single generic
module reused by every runner, by any QC verification scripts, and by the
future Phase-2 adapter that consumes calls.json into the KG.

As of 2026-05 this module doesn't exist yet (each of the four existing
runners inlines its own I/O). The starter file is at
[`assets/tool_calls_io_template.py`](assets/tool_calls_io_template.py) —
copy it once on the first new tool that needs it:

```bash
cp .claude/skills/add-a-tool/assets/tool_calls_io_template.py \
   multiomics_kg/utils/tool_calls_io.py
```

Then write `tests/test_tool_calls_io.py` covering path resolution +
round-tripping + iteration across a fixture genome CSV. The tests pin the
file-naming convention (`.limited_<N>.` infix, indented JSON, sort_keys)
so a future schema drift breaks loudly.

Subsequent tools just `from multiomics_kg.utils.tool_calls_io import save_calls, iter_strain_calls`.
Later, opportunistically migrate psortb-run / tcdb-diamond / eggnog-run /
signalp-run to use it (deferred — not blocking).

What this buys you:

- Runner becomes one-liner-per-write: `save_calls(data_dir, "psortb", strain, calls)` instead of inline path-joining + `json.dump`.
- Cross-strain QC scripts (Step 5 narrative) and any future spot-check verifier use the same iteration: `for strain, dd, calls in iter_strain_calls("psortb"): ...`.
- Phase-2 adapter walks every strain's calls.json with the same iterator — no per-tool I/O code duplicated.
- Schema invariants (indent, sort_keys, file naming) live in one place.

Schema-aware iteration (e.g., "yield every record from a per-protein calls
dict OR a per-region calls list") stays in the **per-tool**
`multiomics_kg/utils/<tool>.py` because it depends on the Step 0 Q1
prediction shape. The generic module only handles paths + JSON load/save
+ cross-strain walks.

### Work — per-strain QC fields in skill_summary.json

Every runner emits QC-flagged fields in
`<strain>.<tool>.skill_summary.json` alongside the existing stats. These
are the per-strain health signals that gate "did this strain run cleanly?"
without anyone needing to read the raw output:

| Field | What it catches |
|---|---|
| `input_proteins` | Total proteins in `protein.faa` (vs `calls_made`) |
| `calls_made` | Proteins with a call (anything other than `MISSING_INPUT` / parse fail) |
| `parse_failures` | Records the tool produced that the parser couldn't read |
| `distribution` | Histogram of the primary call value — PSORTb's `localization`, tcdb-diamond's tier — shape at a glance |
| `sentinel_rate` | Fraction of calls that hit the tool's "no call" sentinel (PSORTb's `Unknown`, score-below-threshold). High = tool struggling with this strain's proteins |
| `wallclock_s` | Anomalies (10x slowdown) usually signal a bad input |
| `tool_version` / `image_digest` | Reproducibility provenance |

tcdb-diamond's `filter_action_distribution` + `pfam_agreement_distribution`
are gold. PSORTb's `localization_distribution` captures the shape. eggNOG's
summary is minimal — when adding a new tool, prefer the tcdb-diamond level
of detail over the eggNOG minimum.

### Document this — fill the data-free SKILL.md sections

Fill in every SKILL.md section that doesn't need cross-strain data yet:

**What It Does** — Plain-prose walkthrough of the orchestrator: loads
genome rows → invokes external tool → parses → writes JSON → status table.
Cite the exact docker command / CLI invocation it builds. ~10–20 lines.

**Output Schema** — Field-by-field for `<strain>.<tool>.calls.json` +
`<strain>.<tool>.skill_summary.json`. Be explicit about nulls, sentinels,
multi-valued fields. **Document the per-strain QC fields** here so future
users know what to look at.

**QC — per-strain field docs** — Document each QC field with its expected
range. Initial thresholds will be calibrated against real data in step 5;
write your best guesses now:

> - `sentinel_rate < 0.40` — anything higher means the tool is failing to classify most proteins; check input encoding and external-tool version.
> - `parse_failures == 0` — non-zero means raw output drifted from what the parser expects.
> - `calls_made == input_proteins` — anything else means the runner skipped records; investigate.

**QC — spot checks** — A handful of known-correct proteins per strain with
expected calls + a `jq` one-liner. These are inline regression checks; when
a future re-run breaks a spot check, something upstream changed:

```markdown
| Strain | Protein ID | Expected | Why this is the ground truth |
|---|---|---|---|
| MED4 | WP_011131900.1 | OuterMembrane (score ≥ 9) | Pal lipoprotein, canonical OM marker |
| MIT9313 | WP_011131123.1 | CytoplasmicMembrane | UniProt-annotated membrane protein |
| MIT1002 | WP_NNNNN.1 | Extracellular | Secreted protease, signal peptide present |
```

```bash
jq '."WP_011131900.1".localization' \
  cache/data/Prochlorococcus/genomes/MED4/<tool>/MED4.<tool>.calls.json
# Expected: "OuterMembrane"
```

Pick spot checks from the literature, UniProt annotations, or existing tool
agreement. Aim for 3–5 per tool, covering at least one Pro strain, one Syn
strain, one heterotroph.

**QC — cross-strain narrative** — Leave as a placeholder; filled in step 5
after the batch.

**Phase 2 (Future)** — One paragraph + link to the design spec under
`docs/superpowers/specs/`. Be explicit that Phase 1 artifacts sit in the
strain cache for inspection and are NOT yet wired into
`gene_annotations_merged.json` or any KG adapter.

## Step 4 — Smoke test via the skill, validate Quick Start + MED4 spot check

### Work

**Invoke the smoke test via the documented Quick Start command — not via
raw python.** This dual-validates the runner AND the documentation. Copy
the command character-for-character from the SKILL.md's Quick Start:

```bash
uv run python .claude/skills/<tool>-run/run_<tool>.py --strains MED4 --limit 100
```

MED4 is the right default — smallest proteome (~1800 proteins) and almost
every paper uses it. The smoke test should finish in seconds-to-minutes
and produce a calls.json you can `jq` to verify the schema.

Then verify the MED4 spot check(s) from the SKILL.md's QC section pass.
If a spot check fails on a smoke test, either the expected value is wrong
(literature didn't actually claim that) or the tool / parser is broken —
fix before proceeding.

`--limit` produces files named `<strain>.<tool>.limited_<N>.calls.json`
that are auto-gitignored.

**Don't start step 5 until smoke + MED4 spot checks pass.** A bug in the
runner means you've wasted hours of batch wallclock.

### Document this — refine Quick Start + Output Schema; fill Workflow When Invoked

If the smoke test surfaced any tweak to the Quick Start (an extra flag, a
warning the user should see), update the SKILL.md now while it's fresh. If
the Output Schema needed a small correction based on real output, fix it.

**Workflow When Invoked** — Numbered checklist of what a user does when
invoking the skill: verify prerequisite (`docker --version`,
`diamond --version`), run install if needed, run the batch, inspect the
status table, **verify spot checks**, dig into FAILED rows in
`logs/<tool>_<strain>.log`. ~5 lines.

## Step 5 — Full batch via the skill, then finalize SKILL.md

### Work

Run the documented full-batch invocation from the SKILL.md's Quick Start:

```bash
nohup uv run python .claude/skills/<tool>-run/run_<tool>.py > logs/<tool>_batch.log 2>&1 &
```

Monitor with `tail -f logs/<tool>_batch.log`. Most tools take hours wallclock.

After completion:

1. Inspect the status table for `FAILED` rows; chase down via `logs/<tool>_<strain>.log`.
2. Verify the cross-strain spot checks from the SKILL.md's QC section all pass. Failures here mean either the spot check was wrong or the tool drifted on certain strains — investigate before committing.
3. Commit the per-strain `<strain>.<tool>.calls.json` + `<strain>.<tool>.skill_summary.json` artifacts (committed by convention — see "Outputs are committed").

### Document this — finalize the SKILL.md

With real cross-strain data in hand, fill in the last two placeholders and
calibrate the QC thresholds:

**Observed batch results** — Cross-strain distribution table + per-strain
wallclock range. PSORTb's section is the canonical shape:

```markdown
## Observed batch results (32-strain run, YYYY-MM-DD)

N proteins classified. Cross-strain distribution:

| <Primary value> | Count | % |
|---|---:|---:|
| ... | ... | ... |

Per-strain wallclock: **X min** (SS120, N proteins) → **Y min** (KT2440, M proteins).
```

**QC — cross-strain narrative** — sanity benchmarks observed across
strains. Example:

> OuterMembrane% ranged 0.5–1.5% in Pro strains, 2.5–4% in heterotrophs —
> matches the lifestyle prediction. SS120 sat at 0.7%, in family. Cross-tool
> sanity: of N proteins called `OuterMembrane` by psortb, M overlap with
> tcdb-diamond's outer-membrane TC families.

**Calibrate per-strain QC thresholds** in the QC field-docs section. The
initial guesses from step 3 may be too loose or too tight; replace with
ranges that match what you actually observed across all 32 strains.

## Step 5.5 — *Optional* — description optimization via `/skill-creator`

The produced `<tool>-run/SKILL.md`'s `description:` field is what makes
`/<tool>-run` actually fire when a user types something like "predict signal
peptides on MED4" or "wire up subcellular localization for the new strain".
A good description triggers on varied phrasings without firing on adjacent
requests; a sloppy one either undertriggers or hijacks unrelated prompts.

For tools whose name doesn't obviously match common user phrasings (or where
you want extra confidence the trigger surface is clean), run
`/skill-creator`'s description-optimization loop:

```bash
# Build a small eval set of should-trigger / should-not-trigger queries:
# 20 entries, ~10 of each, realistic prompts (see /skill-creator §"Generate trigger eval queries")

# Then optimize:
python -m scripts.run_loop \
  --eval-set <path-to-trigger-eval.json> \
  --skill-path .claude/skills/<tool>-run \
  --model claude-opus-4-7 \
  --max-iterations 5 \
  --verbose
```

It evaluates the current description, proposes improvements, iterates up to
5 times. The output JSON has a `best_description` — paste that into the
SKILL.md frontmatter, commit, done.

Skip this step if the tool's name and description already match what users
naturally say (most pure-prediction tools — psortb, signalp, tcdb-diamond —
trigger fine on their default descriptions).

The rest of `/skill-creator` (interview, test-case eval loop, packaging) is
**not** relevant here — `add-a-tool` already prescribes the structure
`/skill-creator`'s interview would otherwise discover, and tool-run skills
live in the repo rather than shipping as standalone `.skill` files.

## Step 6 — Register to be run on new strains

Edit [`add-a-strain/SKILL.md`](../add-a-strain/SKILL.md) step 3 ("Per-strain
tools — background, parallel") and append a `nohup` invocation for the new
runner:

```bash
# <Tool description>
nohup uv run python .claude/skills/<tool>-run/run_<tool>.py --strains <NEW_STRAIN> > logs/<tool>_<NEW_STRAIN>.log 2>&1 &
```

This is what makes future brand-new strains pick up the tool automatically.
Don't skip this — a tool that isn't in the batch list silently gets stale on
new strains.

If the tool's Phase 2 will later be consumed by `prepare_data.sh --steps 1 2`
(eggNOG pattern — output read by `build_gene_annotations.py`), also flag that
in the SKILL.md so `add-a-strain` step 12 ("loop back through prepare_data")
applies. Phase-1-only tools (psortb, tcdb, signalp as of 2026-05) don't need
the loop-back.

---

## Reference material

### File layout to produce

For a new tool called `<tool>` (e.g., `signalp`):

```
multiomics_kg/utils/<tool>.py                              # pure-Python parsing + scoring; no filesystem
tests/test_<tool>.py                                       # unit tests for the utility module
.claude/skills/<tool>-run/SKILL.md                         # workflow doc
.claude/skills/<tool>-run/run_<tool>.py                    # orchestrator: load_genome_rows → external tool → calls.json
.claude/skills/<tool>-run/{Dockerfile,build.sh,INSTALL.md} # IF building from source (Flavor C) — see signalp-run
docs/superpowers/specs/<YYYY-MM-DD>-<tool>-design.md       # design spec (see "Before you start")
docs/superpowers/plans/<YYYY-MM-DD>-<tool>-skill.md        # checkbox implementation plan
~/tools/<TOOL>/                                            # OUTSIDE the repo — installation + reference DBs; host-shared across checkouts; never committed
cache/data/<org>/genomes/<strain>/<tool>/                  # per-strain outputs (COMMITTED — see "Outputs are committed")
  <strain>.<tool>.raw.{tsv,gff,json}                       # raw tool output preserved
  <strain>.<tool>.calls.json                               # post-processed; key shape depends on Step 0 Q1 (WP_ for per-protein tools)
  <strain>.<tool>.skill_summary.json                       # per-strain stats
logs/<tool>_<strain>.log                                   # full stdout+stderr (auto-gitignored via *.log)
```

### Outputs are committed

Per-strain `<strain>.<tool>.calls.json` + `<strain>.<tool>.skill_summary.json`
go **under git** alongside the existing eggnog/tcdb/psortb/signalp artifacts.
The repo's `.gitignore` has an explicit comment locking this in:

```
# Full-run outputs (<strain>.tcdb.{tsv,calls.json,skill_summary.json}) are
# intentionally committed alongside eggnog/ outputs.
cache/data/*/genomes/*/tcdb/*.limited_*
```

Only smoke-test artifacts (`*.limited_*` from `--limit N`) are ignored.
Reasons to commit full-run outputs:

- Reproducibility — anyone can `git pull` and inspect a strain's calls without re-running an 8-hour batch
- PR review — schema changes show up as readable JSON diffs
- Cross-tool joins — Phase-2 integrations can read sibling tool outputs without coordinating runtimes

When adding a new tool, append a one-line `*.limited_*` rule to `.gitignore`
for its directory (e.g., `cache/data/*/genomes/*/<tool>/*.limited_*`) so
smoke tests don't pollute commits.

### Installation lives outside the repo

`~/tools/<TOOL>/` is the canonical location for everything that isn't a
per-strain artifact: the tool binary (if not on PATH), reference databases,
diamond indices, downloaded FASTAs, model weights, license tarballs, Docker
image build context that's been baked into a local image. **Never** commit
any of this. Reasons:

- Size — eggNOG's DB is 39 GB, MNX is ~4 GB; would blow up the repo
- Licensing — PSORTb / SignalP / eggNOG / MetaNetX have non-permissive licenses the repo can't redistribute
- Host-sharing — a single download serves every checkout on the same machine

Make the path overridable via env var in `.env` (`SIGNALP_DATA_DIR`,
`TCDB_DATA_DIR`, etc.) parsed with `python-dotenv`, defaulting to
`~/tools/<TOOL>/` when unset.

### Canonical CLI surface

Existing tools diverge slightly (psortb / eggnog / signalp use `--strain`;
tcdb-diamond uses `--strains` plural; signalp uses `--cpu` while tcdb uses
`--threads`). The standard going forward is below — match it for new tools
and migrate old ones opportunistically.

Required flags every `run_<tool>.py` must support:

| Flag | Purpose | Notes |
|---|---|---|
| `--strains <a> <b> ...` | Run one or more strains by name | Plural canonical (`nargs='+'`); accepts a single value too. Existing `--strain` (singular) skills work — keep both forms via `argparse`. |
| `--force` | Re-run even if `<strain>.<tool>.calls.json` already exists | Default is to skip strains that already have output (idempotency). |
| `--limit N` | Smoke test: only first N proteins (or records) per strain | Output filenames get `.limited_<N>.` infix and are auto-gitignored. |

Common optional flags:

| Flag | Purpose | Notes |
|---|---|---|
| `--threads N` | Parallelism within the external tool | Standardize on `--threads` (tcdb-diamond pattern). Avoid `--cpu` going forward; signalp-run uses it for legacy reasons. |
| `--prepare-<external>` | One-shot: download DB / pull image / build local image, then exit | E.g., `--prepare-image` (psortb), `--prepare-db`. |
| `--refresh-<external>` | Re-download the reference DB / re-pull the Docker image | E.g., `--refresh-image` (psortb), `--refresh-tcdb` (tcdb-diamond). |

**About `--test`:** the project-wide `create_knowledge_graph.py --test`
means "stop each adapter after 100 items" — a different semantic from a
tool runner. Don't introduce `--test` to tool runners; `--limit N` is the
explicit, unambiguous equivalent.

**About `--qc`:** there is no `--qc` flag. Per-strain QC fields in
`skill_summary.json` are emitted on **every** run, not behind a flag —
they're always on. Spot checks live in the SKILL.md as a table + `jq`
one-liner; the user verifies them manually after the batch (Steps 4 + 5
both call for it explicitly). The four existing tool runners follow this
pattern. If you ever want machine-verified spot checks, that's a separate
follow-up — currently out of scope.

Status table to stdout — columns like `strain | n_proteins | n_calls |
wallclock_s | status`. Status values: `OK` / `SKIPPED` (calls.json exists, no
`--force`) / `MISSING_INPUT` / `FAILED`. Full external stdout+stderr goes to
`logs/<tool>_<strain>.log`.

### Output schema conventions

**`<strain>.<tool>.calls.json`** — shape driven by Step 0 Q1 (see "Choosing the calls.json key" in Step 3). For per-protein tools (the common case, all four existing tools), it's a top-level dict keyed by WP_ accession:

```json
{
  "WP_011131900.1": {
    "<tool>_call": "<primary_value>",
    "score": 0.97,
    "<tool>_specific_fields": "..."
  }
}
```

For per-region / per-genome tools, it's a top-level list of records or a single scalar object — see the key-choice table in Step 3.

Conventions independent of key choice:

- Key on `protein_id` (WP_) for protein-input tools, **not** locus_tag — many strains share locus prefixes across assemblies, and WP_ is the join key in `gene_mapping.csv`.
- Sentinel values for "the tool ran but didn't make a call" (PSORTb's `Unknown`, score-below-threshold, etc.) — be explicit; don't conflate with missing-key.
- If the tool emits multiple candidates per record (tcdb-diamond's per-3-part-TC-family candidates), nest them in a `calls: [...]` list sorted by confidence descending — friendlier to downstream filtering than collapsing to a single winner.
- Include the raw scores even when you also emit a tier/bucket — Phase 2 may want to re-threshold.

**`<strain>.<tool>.skill_summary.json`** — per-strain stats:

```json
{
  "strain": "MED4",
  "tool_version": "...",
  "image_digest": "sha256:...",
  "input_proteins": 1858,
  "calls_made": 1842,
  "wallclock_s": 522,
  "distribution": {"<primary_value_1>": 800, ...}
}
```

The summary is what you cite in the SKILL.md "Observed batch results" table.

## Patterns to follow / avoid

### Follow

- **Pure logic in `multiomics_kg/utils/<tool>.py`, orchestration in `.claude/skills/<tool>-run/run_<tool>.py`, generic filesystem I/O in `multiomics_kg/utils/tool_calls_io.py`.** Three-tier split keeps the testable surface large and ensures runner + QC + Phase-2 adapter all read/write through the same code.
- **`Dockerfile` + `build.sh` + `INSTALL.md` next to the runner for Flavor C tools.** signalp-run is the template.
- **Write the SKILL.md skeleton BEFORE the smoke test.** Then invoke smoke + batch via the documented Quick Start commands — dual-validates the runner and the documentation.
- **Bake QC fields into skill_summary.json**, don't bolt them on later. `sentinel_rate`, `parse_failures`, `distribution`, etc. — present from the first run so anomalies are visible from the first smoke test.
- **Ship 3–5 spot checks in the SKILL.md.** Known-correct proteins with expected calls + a `jq` one-liner. They're inline regression checks anyone can run; they catch upstream drift across tool/DB updates.
- **Status table to stdout, full output to log file.** Humans read tables; debugging needs the log.
- **WP_ accession as the calls.json key for protein-input tools.** Join key into `gene_mapping.csv` regardless of locus-prefix differences across strains. Tools consuming GFF / genome.fna pick the natural identifier for their prediction shape — see Step 3's "Choosing the calls.json key" table.
- **`--prepare-<external>` flag that exits after setup.** Prevents accidental 8-hour runs with a stale reference DB.
- **Smoke test (step 4) before batch (step 5).** Treat them as separate gates.
- **Document gotchas in the SKILL.md.** Portability bugs (extractTCDB.pl, PSORTb's `/tmp/results` hardcode) are exactly what the next user will hit.

### Avoid

- **Don't put the new tool's output into `gene_annotations_merged.json` yet.** Phase 2 — separate spec, separate review.
- **Don't write graph-shape data (nodes/edges) from the runner.** That's an adapter, not a tool runner.
- **Don't hardcode a per-strain protein FASTA path.** Use `row["data_dir"] / "protein.faa"`.
- **Don't skip strains silently when input is missing.** Emit a `MISSING_INPUT` status row so the cause is visible.
- **Don't gitignore per-strain calls.json + skill_summary.json.** Committed by convention; only `*.limited_*` gets ignored.
- **Don't put the install under the repo.** `~/tools/<TOOL>/` is outside the repo. Always.
- **Don't reinvent strain iteration.** Use `load_genome_rows` from `multiomics_kg.download.utils.cli`.
- **Don't forget step 6.** A tool not in `add-a-strain` is a regression waiting to happen.

## Reference checklist

Final pass before merging:

- [ ] Step 0 intent captured — 4-bullet answer block (prediction shape, input file, install flavor, trigger phrasings)
- [ ] Design spec under `docs/superpowers/specs/<date>-<tool>-design.md` (optional but recommended; expands on Step 0)
- [ ] Implementation plan under `docs/superpowers/plans/<date>-<tool>-skill.md`
- [ ] Install script in place (Dockerfile + build.sh + INSTALL.md for Flavor C; `--prepare-image` flag for A; `--refresh-<external>` for B)
- [ ] Install actually run on this host; `~/tools/<TOOL>/` or `docker images` confirms
- [ ] `multiomics_kg/utils/<tool>.py` + `tests/test_<tool>.py` — pure Python parsing/scoring, unit-tested
- [ ] `multiomics_kg/utils/tool_calls_io.py` exists (created by the first tool to need it) with `calls_path` / `skill_summary_path` / `load_calls` / `save_calls` / `load_skill_summary` / `save_skill_summary` / `iter_strain_calls` — covered by `tests/test_tool_calls_io.py`
- [ ] `.claude/skills/<tool>-run/run_<tool>.py` uses `load_genome_rows` + `tool_calls_io` helpers, supports `--strains` / `--force` / `--limit`
- [ ] Runner emits per-strain QC fields in `<strain>.<tool>.skill_summary.json` — at minimum `input_proteins`, `calls_made`, `parse_failures`, `distribution`, `sentinel_rate`, `wallclock_s`, `tool_version` (or `image_digest`)
- [ ] Tool-run SKILL.md skeleton written BEFORE smoke — has **Quick Start** (single-strain invocation), **One-Time Setup** (fresh install on a new machine; disk footprint + `.env` entry + verification command), **What It Does**, **Output Schema**, **QC** (per-strain field docs + cross-strain narrative placeholder + spot-checks table), **Observed batch results** (placeholder), **Phase 2 (Future)**, **Workflow When Invoked**
- [ ] QC section has at least 3–5 **spot checks** — known-correct proteins per strain with expected calls + `jq` verification one-liner
- [ ] Smoke test on MED4 with `--limit 100` invoked via the documented Quick Start command (not raw python); MED4 spot check(s) pass
- [ ] Full-batch run completed via the documented Quick Start full-batch command; status table clean
- [ ] Cross-strain spot checks all pass after the full batch
- [ ] SKILL.md finalized after batch — Observed batch results table populated; cross-strain QC narrative filled in; per-strain QC thresholds adjusted to match observed ranges
- [ ] Per-strain calls.json + skill_summary.json committed
- [ ] `.gitignore` gets a `cache/data/*/genomes/*/<tool>/*.limited_*` rule
- [ ] `add-a-strain/SKILL.md` step 3 lists the new runner
- [ ] `CLAUDE.md` updated if the tool changes Architecture / Adapter section or adds a prepare_data step
- [ ] `pytest -m "not slow and not kg"` passes

## Worked examples to copy from

| Tool | Install flavor | External interface | Per-strain wallclock | Pattern highlights |
|---|---|---|---|---|
| `psortb-run` | A (public Docker pull) | `brinkmanlab/psortb_commandline:1.0.2` | ~9-30 min | `--prepare-image`, `--user $(id -u):$(id -g)`, container's `/tmp/results` mounted to host data_dir |
| `tcdb-diamond` | B (CLI + downloaded DB) | `diamond` + downloaded TCDB FASTA | ~30 s setup + seconds per strain | Tier policy + multi-call schema + filter_action annotations; tests cover the policy |
| `eggnog-run` | B (CLI + 39 GB DB) | `emapper.py` | ~5-15 min cached | Output consumed by `prepare_data` step 2 (Phase-2 integration done) |
| `signalp-run` | C (build from source) | local image `signalp:local` | ~? per strain | `build.sh` + `Dockerfile` + `INSTALL.md` next to runner; license tarball goes under `~/tools/SignalP/` |

## Phase 2 hand-off (when ready, separately)

When a tool moves from Phase 1 to Phase 2, it grows one of these surfaces:

- **Merged into `gene_annotations_merged.json`** — e.g., eggNOG. The merge happens in `multiomics_kg/download/build_gene_annotations.py`, configured via `config/gene_annotations_config.yaml`'s `logical_sources` block. Add a new logical source for the tool's output, with field-merge rules.
- **New KG adapter** — e.g., PSORTb's future `localization_adapter.py` that creates `SubcellularLocation` nodes + `Protein_located_in` edges. Write the adapter; add to `create_knowledge_graph.py`; update `config/schema_config.yaml`; update `CLAUDE.md` "Actual Neo4j labels" + "Key graph facts"; add KG-validity tests.

Phase 2 is **out of scope** for this skill — write a new design spec under
`docs/superpowers/specs/` when ready.
