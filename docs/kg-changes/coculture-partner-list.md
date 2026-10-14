# KG Change: `coculture_partner` becomes a dense list; partner role moves to the edge

**Date:** 2026-10-13
**Plan:** `plans/coculture_partner_background.md`
**Dialog with the explorer:** `multiomics_explorer/docs/kg-specs/2026-10-13-coculture-partner-list-explorer-asks.md`
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P1 item 2
(friction log 2026-10-06: `list_experiments(coculture_partner="Alteromonas")` missed the
Weissberg 2025 proteomics because the partner was the culture's background, not its treatment).

## What was wrong

`Experiment.coculture_partner` (a string) and the `Tests_coculture_with` edge were set only from the
paperconfig's `treatment_organism`, i.e. only when coculture *was* the treatment. Measured over
every wired paperconfig on 2026-10-13:

| Role of coculture | Partner recorded | Experiments |
|---|---|---|
| treatment (`'coculture' IN treatment_type`) | yes | 76 |
| background (`'coculture' IN background_factors` only) | **no** | **35** |
| background | yes (`treatment_organism` on a growth-phase / light contrast) | 3 |

The 35: Moreno 2023 (20), Barreto 2022 (5), Weissberg 2025 (4), Biller 2018 (2), Coe 2024 (2),
Hennon 2017 (2). Moreno's cultures contain **two** identified heterotrophs at once, so no scalar
could have represented them.

## The contract now

| | Before | After |
|---|---|---|
| `Experiment.coculture_partner` | `str`, sparse (treatment partner only) | **`str[]`, dense**: sorted `preferred_name`s of every organism sharing the culture (treatment, background, phage); `[]` when axenic |
| who writes it | adapter, from `treatment_organism` | adapter seeds it; **post-import overwrites it from the edges** (`apoc.coll.sort(collect(o.preferred_name))`), so list == edge set by construction |
| `Tests_coculture_with` | one edge, treatment partner only, no properties | **one edge per partner**, `partner_role ∈ {treatment, background}` on every edge (phage = treatment) |
| role derivation | — | **edge only.** `'coculture' IN treatment_type` is NOT a role rule: 8 phage experiments are `viral`, and Biller 2016's growth-phase contrasts have a background partner |
| paperconfig | `treatment_organism` / `treatment_taxid` / `treatment_assembly_accession` | + optional `background_organisms: [{organism, taxid \| assembly_accession}]` |
| validator | warned when a partner was unaccounted | **errors** when `coculture` is in `background_factors` and no partner is named; warns when a `treatment_organism` is read as background |
| `Publication.organisms` | rolled up the scalar | rolls up the list (6 publications gain organisms) |
| `OrganismTaxon` | 48 | **47** — `ncbitaxon:277` *Meiothermus ruber* retired (see below) |

Adapter rule (`omics_adapter._experiment_partners`): a `treatment_organism` is `treatment` iff
`treatment_type ∩ {coculture, viral} ≠ ∅`, else `background`; every `background_organisms` entry is
`background`; an organism named on both sides keeps `treatment`. Edge ids are
`{doi}_coculture_{experiment_key}_{organism node id}` (were per-experiment, now per partner).

Structural edges (`Has_experiment`, `Tests_coculture_with`) are now emitted before the
"no supplementary materials" early return in `get_edges`, so an experiment-bearing paperconfig
without tables still links its Experiment nodes.

## Data corrections that rode along

- **Moreno 2023, 40 heterotroph-side experiments (curator-confirmed).** All three sides of each
  condition compare "Light/Dark + glucose" vs "no glucose"; coculture is the culture state, not the
  treatment. `treatment_type: ['carbon','coculture'] → ['carbon']`, `coculture` added to
  `background_factors` on the 20 Alteromonas-side entries that lacked it, partner moved to
  `background_organisms`. Every one of the 60 Moreno experiments now lists the other two members
  of its culture (cyanobacterium / Marinobacter (MarRef v6) / Alteromonas (MarRef v6)). Visible
  effect: `treatment_type = ['coculture']` loses 40 experiments that were never coculture-vs-axenic
  contrasts. **No experiment mixes partner roles afterwards** (a KG validity test asserts it).
- **Barreto 2022, 5 entries** (4 pooled pCO₂ main effects + 1 interaction term over a design that
  spans axenic and coculture arms): partners recorded as `background` so the experiments are
  discoverable; `background_factors = ['axenic','coculture']` already marks the design as pooled.
  The `background` wording is therefore "present in the culture the measured samples came from",
  not "present in every arm" (corrected 2026-10-14 after explorer review).
- **Biller 2016 (2) and Bernstein 2017 (1)**: `treatment_organism` → `background_organisms`.
- **Bernstein 2017's partner edge retargeted** from `ncbitaxon:277` to the genome node
  `insdc.gcf:GCF_000836395.1` (MruberA *is* "M. ruber strain A"), per explorer ask; the treatment
  row is retired from `treatment_organisms.csv` (comment kept), "Meiothermus ruber" added to the
  validator's canonical genome-organism set. Two `OrganismTaxon` nodes no longer share a name, which
  the new list invariant needs.

## Consumers

- **Explorer** (absorbs on the dev build, ships `0.1.0-alpha.6`; KG release sets
  `mcp_min_version = 0.1.0a6`, a `warn` verdict in `kg_release_info`, not a hard stop): the two
  `CONTAINS` filter sites become `ANY(p IN e.coculture_partner WHERE toLower(p) CONTAINS …)`; four
  Pydantic models typed `str | None` must accept a list; new `coculture_role` filter reads
  `EXISTS { (e)-[r:Tests_coculture_with]->(o) WHERE r.partner_role = $role AND … }`.
- Anyone reading `coculture_partner` as a string: it is a list now; the role is on the edge.

## Verification

- Unit: `tests/test_omics_adapter_experiment.py::TestCoculturePartners` (7: treatment / phage /
  background-via-treatment_organism / `background_organisms` / two partners / combined / axenic),
  `tests/test_paperconfig_validation.py::TestBackgroundOrganisms` (6) + the parametrized
  validation of every wired paperconfig.
- KG validity (`tests/kg_validity/test_post_import.py`): dense list with no `''`; list == sorted
  edge partners; every coculture/viral experiment names ≥ 1 partner; non-coculture experiments carry
  `[]` and no edge; `partner_role` non-null on every edge, phage = treatment, no mixed roles.
  `test_organism.py`: 47 nodes, 4 treatment organisms.
