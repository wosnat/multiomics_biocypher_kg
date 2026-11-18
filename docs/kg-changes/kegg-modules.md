# KEGG modules and per-organism module completeness (2026-11-19)

Spec: `docs/superpowers/specs/2026-11-18-kegg-modules-design.md` (explorer-approved, §10–§13).
Plan: `docs/superpowers/plans/2026-11-18-kegg-modules.md`.

KEGG publishes 583 **modules**: small functional units, each a formula over KOs
(M00531 *Assimilatory nitrate reduction* = `(K00367,K10534,(K00372-K00360)) (K00366,K17877,(K26139+K26138),K00361)`).
The KG now carries them as an ontology and, for every genome, says which modules are complete,
which steps are missing, and which steps were filled only through an EC number.

## Graph shape

| Element | Count | Notes |
|---|---|---|
| `KeggModule` level 3 (`level_kind = kegg_module`) | 414 | 382 `pathway`, 32 `signature`; 384 with a KO in the KG + 30 reached only by an organism edge (EC-only steps or a module reference) |
| `KeggModule` levels 0–2 (`kegg_module_type` / `_class` / `_subclass`) | 2 / 12 / 43 | the `ko00002` class tree, pruned to the kept modules |
| `Kegg_module_is_a_kegg_module` | 469 | module → subclass → class → type |
| `Kegg_term_in_kegg_module` (KO → module) | 1,722 | `steps` (int[], 1-based), `component` (`required` \| `optional`) |
| `Kegg_module_in_kegg_pathway` (module → map) | 1,254 | pruned to the pathway maps the KG keeps (global maps such as `ko01100` included) |
| `Organism_has_kegg_module` | 26,004 | 118 organisms; one edge per (organism, module) with ≥ 1 step covered |

Counts are from the build of 2026-11-19 (KEGG release below).

### Ids and the `kegg.module:` prefix

- Modules: `kegg.module:M00531` — a real KEGG module CURIE (bioregistry `kegg.module`, `^M\d+$`).
- Class nodes: `kegg.module:<slug>[.<slug>[.<slug>]]`, e.g.
  `kegg.module:pathway_modules.energy_metabolism.nitrogen_metabolism`. **House-minted** under the
  same prefix (TigrRole mainrole precedent) so one prefix routes to one label. They do **not**
  match bioregistry's `^M\d+$` and will fail external CURIE validation; only
  `level_kind = 'kegg_module'` ids are real KEGG CURIEs. Slugs, not BRITE positions, because
  positions shift when KEGG inserts a class.

### `KeggModule` properties

`name`, `level`, `level_kind`, `module_type` (dense: `pathway` | `signature`; class nodes carry
their tree's type, so `module_type = 'pathway'` filters at any level), and on modules only
`module_id`, `definition` (KEGG's formula, verbatim; a multi-line definition keeps its `\n`) and
`step_count`.

Post-import:
- `gene_count`, `organism_count` — **subtree**, through the KO bridge only
  (`Gene_has_kegg_ko` → `Kegg_term_in_kegg_module`, then `is_a*0..`), the set
  `genes_by_ontology` reads. EC-only genes are not counted (no KO edge reaches them).
- `direct_gene_count` — modules only; equals `gene_count` there.
- `complete_organism_count` — organisms with a `complete` edge to a module in the subtree.

**`organism_count` is not the number of `Organism_has_kegg_module` edges.** It counts organisms
with a gene whose KO is in the module. The edge count also includes organisms that reach the
module only through EC-only genes or a module reference, so it can be larger, and the 30
KO-less modules have `organism_count = 0` with real edges. `complete_organism_count` counts edges.

### `Organism_has_kegg_module`

| Property | Meaning |
|---|---|
| `status` | `complete` \| `one_step_missing` \| `two_steps_missing` \| `partial` (≥ 3 missing). R3 exception #2 (see `vocabulary-contract.md`). `complete` = every step has an annotated gene; not shown to run. |
| `steps_total` | = the module's `step_count` |
| `steps_present` | steps covered, KO or EC |
| `steps_present_via_ec` | of those, steps covered **only** through the EC fallback; `steps_present - steps_present_via_ec` is the KO-only count |
| `missing_step_indices` | 1-based, same numbering as `Kegg_term_in_kegg_module.steps`; `[]` when complete |
| `missing_step_definitions` | the formula of each missing step, parallel to `missing_step_indices` |
| `missing_kos` | KOs written **directly** in the missing steps, alternatives and complex subunits flattened, optional subunits left out. **One KO from this list may not complete a step** (a `+` complex needs all its subunits). A missing module-reference step (`M0xxxx`) contributes nothing; `missing_step_definitions` shows the module id. |
| `ec_filled_step_indices` | sparse; steps counted only through EC |
| `ec_genes` | sparse; **parallel to `ec_filled_step_indices`** (same length, same order): the locus tag the matching assigned to that step, or, for a module-reference step, the referenced module id (`M0xxxx`). When several EC-only genes could fill a step, only the matched one is listed. |

KO-only missing steps = `missing_step_indices ∪ ec_filled_step_indices`. An organism with no step
covered has no edge; "absent" is a left join against all modules.

## How completeness is computed

Build time, `multiomics_kg/utils/kegg_modules.py` (evaluator) + `kegg_module_adapter.py`, from
each strain's `gene_annotations_merged.json` KOs filtered to the KOs the KEGG adapter emits (so it
sees exactly the `Gene_has_kegg_ko` set; deleted KOs never count).

- Space separates steps; `,` is OR; parentheses group, a space inside means AND; `+` joins complex
  subunits; `-K…` is an optional subunit (ignored).
- `--` (a reaction with no KO) is always a whole step in this release and is **dropped** from the
  step list, otherwise those 14 modules could never be complete. A step written as one optional
  subunit on its own (`-K08080` in M00955; 7 modules) is dropped the same way and its KO gets no
  `Kegg_term_in_kegg_module` edge.
- A multi-line definition's lines are further steps.
- A step naming another module is covered when that module is complete; if that module needed an
  EC-filled step, the referencing step counts in `steps_present_via_ec` with the module id in its
  `ec_genes` slot. A KO serving several steps fills all of them (Mapper parity).

**EC fallback.** A gene with **no live KO** and a **full** EC number (no `-`) can stand in for a
KO that KEGG links to that EC (`/link/ko/enzyme`), except:
1. never for a KO inside a complex (`+` or `-` group), because one EC often names every subunit;
2. a gene is a candidate for a step only if it completes the step on its own (with the
   organism's KOs); a step that needs two EC-only genes stays missing;
3. one step per gene per module, assigned by maximum bipartite matching;
4. only over steps the KO evidence leaves uncovered — EC never competes with or replaces a KO.

Measured: 3,470 EC-filled steps; 6,528 → 6,985 complete pairs (`one_step_missing` 5,864, `two_steps_missing` 4,788, `partial` 8,367). Example: NATL2A's `pyk`
(PMN2A_0295, EC 2.7.1.40, no KO) completes M00002 Glycolysis core.

## Cross-check against KEGG

KEGG Mapper's web form does not accept scripted calls, so the evaluator was checked against
KEGG's own per-genome module lists instead: run on KEGG's own KO assignments
(`/link/ko/<org>`), it reproduces KEGG's module calls (`/link/module/<org>`) for MED4 (48/48),
*E. coli* (112/112), *P. putida* KT2440 (92/92), WH8102 (53/53), human (131/131), *Arabidopsis*
(126/126) and yeast (85/85): every module KEGG lists is complete here. The only extras are
"module set" signature modules built only from module references (M00611, M00364, M00366,
M00367), which KEGG lists per gene and so never for a module with no genes of its own. **The `--` rule is confirmed:** every
module KEGG lists that contains a `--` step (8 across human, *Arabidopsis* and yeast) is complete
here with `--` dropped. On the KG's own MED4 KOs (1,014, vs KEGG's 982 for `pmm`) 49 modules
are complete; 3 of KEGG's 48 are not and 4 others are, all from KO-assignment differences.

## What `complete` does not mean

The route runs; the substrate is used; the genes are expressed together; or that a missing step is
absent from the genome (unannotated KO, or a draft assembly — see `OrganismTaxon.assembly_level`).

Known annotation gaps (one-step photosystem complexes, so a missing subunit KO means no edge at
all): M00161 Photosystem II is not complete in *S. elongatus* PCC 7942 (no psbD K02706),
Synechococcus CB0205 (psbD, psbC K02707) and Prochlorococcus HNLC1 (psbB K02705, psbD; a scaffold
draft); M00163 Photosystem I is not complete in HNLC1 (K02693). Every other cyanobacterium has
both complete; no heterotroph does.

## KEGG release

Module definitions change between KEGG releases (M00531 gained an optional subunit between the
spec and the build). Step 6 records `/info/kegg`'s per-database update dates as
`kegg_data.json` `kegg_release`, stamped post-import on **`Schema_info.kegg_release`**, e.g.
`pathway 2026-11-18, brite 2026-11-17, module 2026-11-09, ko 2026-11-18, compound 2026-11-18, reaction 2026-11-18, enzyme 2026-11-16`.
KEGG REST publishes no release number, so the dates are the release.

## Data flow

Step 6 (`build_kegg_metabolism_xrefs`) fetches six raws into gitignored `cache/data/kegg/raw/`
(`list_module.txt`, `link_pathway_module.txt`, `link_ko_enzyme.txt`, `br_ko00002.json`,
`get_modules.txt` — `get/M…` batched 10 per call — and `info_kegg.txt`), pinned in
`cache/data/step6_raw_inputs.json`, and writes three `kegg_data.json` blocks: `modules` (all 583,
with `class_path` and `pathways` as `ko…` ids), `module_ec_to_kos` and `kegg_release`. A step-6
`--force` run on a checkout whose committed manifest predates these raws fetches them once
instead of failing. KO membership and step numbers come from parsing each definition, not from
`link/ko/module` (which has neither).

## Not in this slice

No gene → module edge (reach a module through the gene's KO); eggNOG's `KEGG_Module` field stays
unused. No reaction modules (`RM*`). No change to `annotation_types` / `annotation_quality` /
BRITE (`ko00002` is not a BRITE tree).

## Explorer

`ONTOLOGY_CONFIG["kegg_module"] = {label: KeggModule, gene_rel: Gene_has_kegg_ko,
bridge: {node_label: KeggTerm, edge: Kegg_term_in_kegg_module},
hierarchy_rels: [Kegg_module_is_a_kegg_module]}`. Enrichment at `level = 3` by default; never mix
levels in one test. The completeness tool reads `Organism_has_kegg_module`.
