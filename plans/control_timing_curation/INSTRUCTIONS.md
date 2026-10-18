# Control-timing curation — instructions for curation agents

> **Historical (2026-10-15 pass).** Since 2026-10-18 `time_matched` is `matched` (same sampling time
> or same growth state) and `statistical_model_id` no longer exists. Current rules: the `/paperconfig`
> skill, "Control timing and design notes".

Plan: `plans/per_arm_sampling_and_shared_model.md` (read § Design and § Curation pass first).
You curate a **batch of papers**. For each paper you (1) read the paper and all its material,
(2) add up to three new fields to its `paperconfig.yaml`, (3) write a report so a human can
review every decision and every uncertainty.

## Scope of a paper

- Work only on experiments that have ≥ 1 `statistical_analyses` entry (search
  `supplementary_materials.*.statistical_analyses[].experiment`). Those are "DE experiments".
  Experiments with no DE contrast (metabolomics, derived metrics, clusters only) get **no** new
  field.
- Change nothing else in the paperconfig. If you notice another problem (wrong `growth_phase`,
  wrong `timepoint`, wrong `control_condition`, mislabelled experiment), **report it** under
  "Other findings" — do not fix it.
- Do not commit, do not touch other papers, do not run prepare_data / Docker.

## What to read (all of it that bears on sampling and statistics)

In the paper's directory (`dirname` of the paperconfig):
1. The main PDF (`papermainpdf`) — **Methods** in full (culture set-up, sampling schedule, RNA /
   protein extraction times, statistics / DE section), plus Figure 1 / the experimental-design
   figure and its legend, plus the Results paragraphs that introduce each comparison.
   Use the Read tool with `pages` (≤ 20 pages per call) for PDFs.
2. Supplementary PDFs / docx / `legends of supp tables.txt` / README.md — supplementary methods
   often hold the sampling schedule and the DESeq2/edgeR/limma design.
3. The supplementary table headers (`head -3` of the CSVs referenced by the paperconfig) — column
   names like `d18_vs_d14` or `T6_vs_T0` reveal what was compared.
4. The paperconfig itself — `timepoint` / `timepoint_hours` / `growth_phase` on each analysis,
   `treatment_condition` / `control_condition`, `test_type`.

## The three fields

### 1. `control_timing` — required for every DE experiment (closed vocab)

| value | meaning |
|---|---|
| `time_matched` | Every contrast in the experiment compares treatment and control arms **sampled at the same time** (treated vs untreated at 1 h; coculture vs axenic both on day 11; steady-state chemostats/turbidostats of two conditions sampled in parallel; acclimated cultures of two conditions harvested at the same growth state). |
| `baseline` | Every contrast compares against a control **sampled at another time** — typically one early reference point (each time t vs t0 / pre-shock / before transfer; "late vs early growth phase"; "day 18 vs day 7" within one culture series). Growth-phase / time-course-vs-start comparisons are `baseline`. |
| `mixed` | Some timepoints of the experiment are time-matched, others are not. Requires the experiment to span ≥ 2 distinct `timepoint` labels. |

Placement:
- `publication.control_timing` = the value most of the paper's DE experiments take. Put it right
  after `papermainpdf` (or after `doi` if that comes later) in the `publication:` block.
- An experiment whose value differs gets its own `control_timing:` line (override). Put override
  lines at the end of the experiment block (after `table_scope` / `table_scope_detail` if present,
  otherwise after the last existing key). Do **not** add an override that equals the paper default.

`time_matched` is a claim checked in the methods, never a fallback. If the methods do not state
when the control was sampled, pick the value the design implies, say so in `design_notes`, and
mark it `low` confidence in the report.

Hints: "control" = the arm the fold change is relative to (denominator). For a time course where
each time point is compared with time 0 of the *same* culture → `baseline`. For a time course
where each time point compares treated vs mock/untreated sampled at that same time → `time_matched`.
Watch for phage infection papers (often each time vs uninfected at t0, or vs mock at each t —
check), stress shocks (vs pre-shock = baseline), diel series (vs a reference time = baseline).

### 2. `statistical_model_id` — sparse, only for shared fitted models

Set **only** when ≥ 2 DE experiments of this paper are contrasts extracted from **one fitted
model** over shared samples (one DESeq2/edgeR design fitted on all samples, one limma lmFit with a
contrast matrix, one GLM with interaction terms). All experiments sharing an id must have the same
`organism` and `omics_type`. Value: a short snake_case key unique within the paper, e.g.
`med4_rnaseq_deseq2` (the adapter later prefixes the DOI — do not add it).

Do **not** set it when: each comparison was fitted separately (pairwise DESeq2 runs, per-comparison
t-tests, Cyber-T per pair, per-comparison proteomics tests); all contrasts from the model sit
inside one experiment (nothing to link); or the paper does not say how DE was fitted (→ no id,
flag it in the report as an open).

### 3. `design_notes` — sparse free text

**Required** when `control_timing` is `baseline` or `mixed` (on that experiment, even if `baseline`
came from the paper default). Recommended on experiments with a `statistical_model_id`. Optional
elsewhere — only when it records something a reader needs (e.g. "methods do not state control
sampling time; time-matched assumed from parallel cultures").

Writing rules:
- Use the exact `timepoint` labels from the analyses (`"day 18"`, `"T6"`).
- Refer to other experiments by their paperconfig key verbatim.
- State the control sampling time and its growth phase using the `growth_phase` vocabulary
  (`exponential`, `nutrient_limited`, `stationary`, `death`, … — see validator line ~275 of
  `scripts/validate_paperconfig.py` for the list; use `unknown` if not determinable).
- For a shared model: what the model is (tool + design formula if given) and any exact linear
  relation between contrasts.
- No `|` and no `'` (apostrophes) anywhere — write "do not", "cells", not "don't", "cell's".
- Use a YAML folded block:
  ```yaml
      design_notes: >-
        day 11: coculture vs axenic, both sampled on day 11. day 18: coculture sampled on day 18 vs
        axenic sampled on day 14 (death) - not time-matched.
  ```
- Short: 1–4 sentences.

## Editing the YAML

- Use the Edit tool with small, exact insertions. Never rewrite or reformat the file; keep
  comments, quoting style and key order. Match the file's indentation.
- After editing, check it parses and the values are what you meant:
  `uv run python -c "import yaml,sys; y=yaml.safe_load(open(sys.argv[1]))['publication']; print(y.get('control_timing')); [print(k, v.get('control_timing'), v.get('statistical_model_id'), bool(v.get('design_notes'))) for k,v in y['experiments'].items()]" "<path>"`
- Also run `uv run python scripts/validate_paperconfig.py "<path>"` and confirm you introduced no
  new errors (it does not know the new fields yet; it must simply not get worse).

## Report — `plans/control_timing_curation/<batch>.md`

The report is the review surface. Be explicit about evidence and uncertainty. Format:

```markdown
# Batch <X> — control-timing curation

## Summary
- Papers: N · DE experiments: M · baseline: a · mixed: b · time_matched: c · shared-model ids: k
- Needs human review: <list of paper/experiment keys with a one-line reason each>

## <papername>
- Sources read: main PDF pp. x–y (Methods, Fig. 1); Supp methods pp. …; CSV headers …
- Paper default: `time_matched`
- Sampling: <one or two sentences describing the sampling schedule as the paper states it>
- DE fitting: <tool, one model or per-comparison, design formula if stated>

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `key` | day 11, day 18 | mixed (override) | med4_rnaseq | high | p. 4 Methods: "…short quote…" |

Opens / uncertainties:
- `key`: <what is uncertain, what you chose, what would resolve it>

Other findings (not fixed):
- <e.g. growth_phase on analysis X says nutrient_limited, paper says death phase (p. 5)>
```

Confidence: `high` = the paper states it explicitly; `medium` = strongly implied by the design /
column names; `low` = guessed — always listed under Opens.

When done, reply with only: the report path, the summary block, and the "Needs human review" list.
