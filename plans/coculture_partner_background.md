# Coculture partner when coculture is the background, not the treatment

**Status:** implemented 2026-10-13 (adapter, validator, post-import, vocab, 78 experiment blocks in 8 paperconfigs, Moreno correction curator-confirmed, `ncbitaxon:277` retired); dev rebuild + `pytest -m kg` in progress; explorer absorb next.
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P1 item 2
(friction log 2026-10-06 "`coculture_partner` filter misses experiments where coculture is
the background, not the treatment").

## Problem

`Experiment.coculture_partner` and the `Tests_coculture_with` edge are set only from the
paperconfig's `treatment_organism`, i.e. only when coculture *is* the treatment. An
experiment run *in* coculture whose treatment is something else (N starvation, glucose,
pCO₂) carries `background_factors: [coculture, …]` but no partner anywhere in the graph. The
explorer's `list_experiments(coculture_partner=…)` filters on the string property, so those
experiments are invisible to it; the analysis first concluded "no *Prochlorococcus*
proteomics in coculture with *Alteromonas*", which was wrong.

Measured 2026-10-13 over every wired paperconfig (114 experiments touch coculture):

| Role of coculture | Partner recorded | Experiments |
|---|---|---|
| treatment | yes | 76 |
| background | **no** | **35** |
| background | yes | 3 |

The 35, by paper, with the partner the experiment text names:

| Paper | n | Measured organism(s) | Partner(s) in the culture |
|---|---|---|---|
| Moreno 2023 | 20 | MED4, SS120, WH7803, WH8102, BL107 (proteomics) | **two at once**: Marinobacter (MarRef v6) + Alteromonas (MarRef v6) — "non-axenic" cultures whose heterotrophs were identified by reference-proteome match |
| Barreto 2022 | 5 | EZ55, MIT9312, CC9311 | EZ55 ↔ MIT9312 / CC9311 — but some entries are a pCO₂ × coculture *interaction* term over a design that spans axenic and coculture arms; "coculture in background" may mean "design includes coculture", not "the culture was a coculture" — **check per entry** |
| Weissberg 2025 | 4 | MED4, HOT1A3 (RNA-seq + proteomics) | MED4 ↔ HOT1A3 |
| Biller 2018 | 2 | NATL2A, MIT1002 | NATL2A ↔ MIT1002 |
| Coe 2024 | 2 | NATL2A, MIT1002 | NATL2A ↔ MIT1002 |
| Hennon 2017 | 2 | MIT9312, EZ55 | MIT9312 ↔ EZ55 |

Two constraints fall out of the table: a background coculture can have **more than one
partner** (Moreno), and the partner's node can be a genome strain, a treatment organism, or a
reference-proteome-match organism — exactly the three id forms the existing edge code already
resolves (`treatment_taxid` → `ncbitaxon:` / strain name → `insdc.gcf:`).

## Scope

**Changing**
1. **paperconfig** — new optional experiment field, a list because of Moreno:
   ```yaml
   background_organisms:
     - organism: "Alteromonas macleodii HOT1A3"
       taxid: 28108          # optional, same resolution rules as treatment_taxid
   ```
   `treatment_organism` / `treatment_taxid` keep their meaning (the partner that *is* the
   treatment). `scripts/validate_paperconfig.py`: accept the field; error when
   `background_factors` contains `coculture` but neither `treatment_organism` nor
   `background_organisms` is set (so the 35 can never silently recur); warn on the reverse.
2. **omics_adapter**:
   - `Tests_coculture_with` emitted for every partner (treatment ∪ background ∪ phage), carrying
     **`partner_role ∈ {treatment, background}`** on every edge (R5 two-state string, new vocab
     entry with `value_descriptions`).
   - `Experiment.coculture_partner: str[]` — **dense**, every partner name sorted, `[]`
     when axenic (post-import re-materialises `[]` like `treatment_type`). Same name, new type.
3. **post-import**: writes `Experiment.coculture_partner` from the edges (sorted `preferred_name`s,
   `[]` when none); `Publication.organisms` / `OrganismTaxon.*` rollups `UNWIND` the list.
4. **Data**: fill `background_organisms` on the 35 experiments (6 paperconfigs), after the
   Barreto per-entry check; migrate the 3 background `treatment_organism` entries (Biller 2016,
   Bernstein 2017); apply the Moreno correction (40 heterotroph-side experiments) once confirmed.
   Every partner must resolve to an existing OrganismTaxon node.
5. **Vocabulary / schema / docs / CHANGELOG / `docs/kg-changes/`**: `Tests_coculture_with.partner_role`,
   `Experiment.coculture_partner` (type change → `### Breaking`, `mcp_min_version` bump); CLAUDE.md Experiment bullet; backlog tick.
6. **Explorer ask** (`multiomics_explorer/docs/kg-specs/`): `list_experiments(coculture_partner=…)`
   should filter `ANY(p IN e.coculture_partner WHERE toLower(p) CONTAINS …)`; the MCP
   docs for `coculture_partner` must say which role it covers.

**Not changing**
- Experiments where coculture is the treatment (76) — their partner, edge and property stay
  exactly as they are.
- `background_factors` semantics; the `coculture` token stays the flag that a culture was
  mixed.
- No new OrganismTaxon nodes — every partner above already has one.

## Decision (2026-10-13, user + explorer dialog)

**`coculture_partner` becomes a dense `str[]` in place, same name.** No second property, no
semicolon-joined string. **Role lives on the edge only** (`Tests_coculture_with.partner_role ∈
{treatment, background}`, set on every edge, phage included) — the node-level "derive from
`treatment_type`" rule was dropped after the explorer showed it misroles the 8 phage experiments
(`treatment_type = ['viral']`, partner `Phage`) and cannot tie a name to a role on a multi-partner
list. **Post-import is the single writer of the list**:
`e.coculture_partner = apoc.coll.sort(collect(o.preferred_name))` over the experiment's edges,
`[]` when none — so the explorer's invariant (list == sorted edge partner names; 87/87 today)
holds by construction. Dialog record + explorer-side change list:
`multiomics_explorer/docs/kg-specs/2026-10-13-coculture-partner-list-explorer-asks.md`
(explorer responses R0.1–R0.3, CCP-001…005 + the KG reply). The explorer absorbs it on the dev
build and ships `0.1.0-alpha.6`; the KG release sets `mcp_min_version = 0.1.0a6` (a `warn`
verdict in `kg_release_info`, not a hard stop — 3 explorer tools fail on the new build until then,
because their Pydantic models type the field `str | None`).

**Adapter role rule** (deterministic, no new semantics on existing fields): a `treatment_organism`
is `treatment` iff `treatment_type ∩ {coculture, viral} ≠ ∅`, else `background`; every
`background_organisms` entry is `background`. Validator warns when a `treatment_organism` is read
as background (today: Biller 2016 ×2, Bernstein 2017 — migrate those to `background_organisms`).

**Moreno 2023 correction (pending curator confirmation).** All 40 heterotroph-side experiments
compare glucose vs no glucose exactly like the cyanobacterium side, yet carry `coculture` in
`treatment_type`. Fix: `treatment_type → ['carbon']`, `coculture` into `background_factors` on the
20 Alteromonas-side entries that lack it, `treatment_organism → background_organisms`. Every one
of the 60 Moreno experiments then lists the other two culture members as `background`; the corpus
has no mixed-role experiment afterwards. Visible consequence (→ `### Data` bullet):
`treatment_type = ['coculture']` filters lose 40 Moreno heterotroph experiments that were never
coculture-vs-axenic comparisons.

## Acceptance criteria

- `list_experiments`-style query `MATCH (e:Experiment) WHERE 'coculture' IN e.background_factors
  AND size(e.coculture_partner) = 0 RETURN count(e)` → 0 (was 35).
- Moreno cyanobacterial experiments carry two `Tests_coculture_with` edges each, both
  `partner_role = background`, to the two MarRef organisms.
- The 56 treatment cocultures outside Moreno + the 8 phage experiments keep exactly one
  `Tests_coculture_with` edge, `partner_role = treatment`, list == `[old scalar]`.
- KG validity: on every Experiment, `coculture_partner` == sorted `preferred_name` over its
  `Tests_coculture_with` edges (`[]` when none); `partner_role` non-null on every edge.
- `validate_paperconfig.py` errors on a fresh paperconfig with `coculture` in
  `background_factors` and no partner.
- `pytest -m "not slow and not kg"` and, after rebuild, `pytest -m kg` green;
  `/omics-edge-snapshot` shows no expression-edge change.
