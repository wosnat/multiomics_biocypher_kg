# KEGG modules + per-organism module completeness — design

**Date:** 2026-11-18 · **Status:** approved by explorer (§10 review, §11 response, §12 sign-off, §13 clarifications) · **Backlog item:** `plans/backlog.md` "KEGG modules + per-organism module completeness"

## 1. Why

The KG carries KEGG as one tree, `ko00001`: category → subcategory → pathway map → KO.
A pathway map is a large drawing that mixes biosynthesis, degradation and side branches,
so "strain X has 12 of the map's 30 KOs" says nothing about whether any route on the map is
there. Analyses that need "can this strain do Y" have to build the route themselves:

- spock_osnat `2026-11-04-aa_utilization_genomics` curated a KO step list for every
  amino-acid degradation route by hand, then spent four critic passes on its errors (a KO
  on a shared step credited as dedicated, isozymes missed, an "arogenate route" premise
  that pulled in positive-control genes). Its gaps doc logs "no KEGG modules" as a KG gap.

KEGG already publishes those routes as **modules**: 583 small functional units, each a
formula over KOs. M00531 *Assimilatory nitrate reduction* is
`(K00367,K10534,K00372+K00360) (K00366,K17877)` — two steps, nitrate → nitrite → ammonia.
KEGG Mapper evaluates a genome against that formula and reports the module complete, or
one or two steps short.

This spec adds the modules and that evaluation to the KG.

## 2. What it delivers

1. **`KeggModule` ontology**: a node per module, plus the `ko00002` class tree as
   `KeggModule` class nodes joined by `Kegg_module_is_a_kegg_module` (§4; *not* a BRITE tree).
2. **KO → module** and **module → pathway map** edges.
3. **`Organism_has_kegg_module`**: per organism, which of the module's steps are covered,
   which are missing, which were filled only through EC numbers and by which genes, and a
   status (`complete`, `one_step_missing`, ...), computed at build time.

The third is the point of the project. 1 and 2 are what it needs to be read; with them the
explorer can also run enrichment over modules (§10, answer 1).

### Measured against the current KG (KEGG release of 2026-11-18, 118 genomes)

Re-measured after the review with the final rules (§5): `--` steps excluded, EC fallback
limited to single-KO terms and one step per gene, deleted KOs excluded.

| What | Count |
|---|---|
| KEGG modules | 583 (532 pathway, 49 signature, 2 complex) |
| … with ≥ 1 KO in the KG | 384 (358 pathway, 26 signature) |
| … reached only through an EC-filled step or a module reference | 30 |
| → `KeggModule` module nodes | **414** |
| KO → module edges (KOs in the KG) | 1,723 (1,082 of the KG's 4,774 KOs sit in some module) |
| module → pathway map edges (maps the KG keeps) | 1,186 (19 modules link to no kept map) |
| steps per module | median 4, max 14 |
| `ko00002` class nodes (pruned to the 414 modules) | 2 types → 12 classes → 43 subclasses = 57 class nodes; every module sits under exactly one subclass |
| `Organism_has_kegg_module` edges (≥ 1 step covered) | **26,005** (24,190 on KOs alone) |
| … `complete` | 6,987 (6,528 on KOs alone) |
| … `one_step_missing` / `two_steps_missing` | 5,862 / 4,786 |
| … `partial` (three or more missing) | 8,370 |
| steps filled only through EC | 3,468 (≈ 0.13 per edge) |

The KG has 4,774 KOs; the release graph (alpha.9) has 4,808. The difference is the 34 KOs
KEGG deleted, which the dev graph already drops (CHANGELOG 2026-11-18).

Biology checks from the dry run (module evaluator in §5, gene KOs from
`gene_annotations_merged.json`):

| Module | MED4 | MIT9313 | NATL2A | WH8102 | HOT1A3 | KT2440 |
|---|---|---|---|---|---|---|
| M00009 TCA cycle | 5/8 | 6/8 | 5/8 | 6/8 | 8/8 | 8/8 |
| M00531 Assimilatory nitrate reduction | 0/2 | 1/2 | 1/2 | 2/2 | 1/2 | 1/2 |
| M00002 Glycolysis, core | 5/5 | 5/5 | 5/5 \* | 5/5 | 5/5 | 5/5 |
| M00161 Photosystem II | 1/1 | 1/1 | 1/1 | 1/1 | — | — |
| M00163 Photosystem I | 1/1 | 1/1 | 1/1 | 1/1 | — | — |

The *Prochlorococcus* TCA gap (no 2-oxoglutarate dehydrogenase), MED4's lack of any
nitrate/nitrite reduction, and MIT9313 / NATL2A having nitrite but not nitrate reductase
are known biology. \* NATL2A is 4/5 on KOs alone — see §5 *EC fallback*. — = no edge.
PSII and PSI are one-step modules (one complex each, optional subunits marked `-K`).

## 3. Goals, non-goals

**Goals**
- "Does organism X have a complete Y" answered by one edge, for every module and genome.
- The missing steps named and indexed, so "1 step missing" can be checked by eye and joined to KOs.
- Modules usable as an ontology like any other (browse by class, by pathway map, by KO;
  `genes_by_ontology`; enrichment).

**Non-goals**
- **No gene → module edge.** A gene reaches a module through its KO
  (`Gene_has_kegg_ko` → `Kegg_term_in_kegg_module`). eggNOG's own `KEGG_Module` field
  (already in `gene_annotations_merged.json` as `kegg_module`, unused) only says the gene's KO
  occurs in some module; it is not used. The genes behind EC-filled steps are listed on the
  organism edge (`ec_genes`), because no KO edge reaches them.
- **No "route runs" claim.** Completeness says the genes for every step are annotated, not
  that the route is active or that the substrate is used. Most KO calls are eggNOG
  `family_inferred`.
- **No reaction modules** (KEGG `RM*`) and no module-to-reaction / compound edges in this slice.
- **No change to `annotation_types` / `annotation_quality`.** Modules are a view of KOs the
  gene already has, not a new evidence source. `ko00002` is deliberately **not** a BRITE tree
  (§4), so the `'brite'` annotation type and the `BriteCategory` counts do not move.

## 4. Graph shape

### `KeggModule` nodes — one label, two kinds

Every node of the module ontology carries the label `KeggModule`, the pattern of every
hierarchical ontology in the KG (TCDB, CAZy, MEROPS, TigrRole): one label, `level` /
`level_kind`, an `is_a` hierarchy.

| `level` | `level_kind` | Example | Id |
|---|---|---|---|
| 0 | `kegg_module_type` | Pathway modules | `kegg.module:pathway_modules` |
| 1 | `kegg_module_class` | Energy metabolism | `kegg.module:pathway_modules.energy_metabolism` |
| 2 | `kegg_module_subclass` | Nitrogen metabolism | `kegg.module:pathway_modules.energy_metabolism.nitrogen_metabolism` |
| 3 | `kegg_module` | Assimilatory nitrate reduction | `kegg.module:M00531` |

`kegg.module` is a registered bioregistry prefix (`^M\d+$`). The class-node ids are
house-minted slug paths under the same prefix — the TigrRole mainrole precedent
(`tigr.role:energy_metabolism` beside numeric subroles) — so one prefix routes to one label
(the explorer's term-id routing) and a module id never collides with a class id. Slug paths
rather than KEGG's positional BRITE numbering, because positions shift when KEGG inserts a class.

| Property | Type | On | Notes |
|---|---|---|---|
| `name` | str | all | KEGG name, "Assimilatory nitrate reduction, nitrate => ammonia" |
| `level`, `level_kind` | int, str | all | as above |
| `module_type` | str | all, dense | `pathway` \| `signature` (closed; KEGG's 2 `complex` modules have no KO in the KG). Set on class nodes too, from their level-0 ancestor, so a consumer filters at any level. Signature modules are marker sets (drug resistance, pathogenicity), not routes. |
| `module_id` | str | modules | `M00531` |
| `definition` | str | modules | KEGG's formula, verbatim, so a reader can check a call |
| `step_count` | int | modules | top-level steps after joining multi-line definitions and excluding `--` placeholder steps (§5); equals `steps_total` on every organism edge |
| `gene_count`, `organism_count` | int | all | post-import; **subtree**, through the KO bridge only (`Gene_has_kegg_ko` → `Kegg_term_in_kegg_module`, then `is_a*0..`), the same semantics as every hierarchical ontology. EC-only genes are not counted here (no KO edge reaches them); they are on the organism edge. |
| `direct_gene_count` | int | modules | post-import; equals `gene_count` on a module (leaf). Omitted on class nodes, where it would be a constant 0 (the non-KO `KeggTerm` precedent, DOC-006). |
| `complete_organism_count` | int | all | post-import; organisms with a `complete` edge to the module (subtree: to any module below the class) |

Full-text index `keggModuleFullText` on `name`; scalar indexes on `level`, `level_kind`,
`module_type`, `module_id`. No `is_uninformative` (§11, answer 6).

### Edges

| Edge | Direction | Properties |
|---|---|---|
| `Kegg_term_in_kegg_module` | KO (`KeggTerm`) → `KeggModule` (level 3) | `steps` int[] — which step(s) the KO serves (1-based; 47 (module, KO) pairs serve more than one) · `component` str — `required` \| `optional` (an optional complex subunit, `-K…` in the formula). Scalar is safe: no KO is required in one step of a module and optional in another (measured 0); the adapter asserts it. |
| `Kegg_module_is_a_kegg_module` | `KeggModule` → parent `KeggModule` | none; module → subclass → class → type |
| `Kegg_module_in_kegg_pathway` | `KeggModule` (level 3) → `KeggTerm` (pathway map) | none; pruned to maps the KG keeps |
| `Organism_has_kegg_module` | `OrganismTaxon` → `KeggModule` (level 3) | see below |

### `Organism_has_kegg_module`

One edge per (organism, module) with **at least one step covered**; no edge means none.
"Absent" is derived by the consumer with a left join against all modules.

| Property | Type | Meaning |
|---|---|---|
| `status` | str | `complete` \| `one_step_missing` \| `two_steps_missing` \| `partial` (KEGG Mapper's bands; `partial` = three or more missing). R3 exception #2 (§8). |
| `steps_total` | int | the module's `step_count` |
| `steps_present` | int | steps covered, KO or EC |
| `steps_present_via_ec` | int | of those, steps covered **only** through the EC fallback (§5); 0 = KO evidence alone. `steps_present - steps_present_via_ec` is the KO-only count. |
| `missing_step_indices` | int[] | 1-based indices of the missing steps — same numbering as `Kegg_term_in_kegg_module.steps` |
| `missing_step_definitions` | str[] | the formula of each missing step, parallel to `missing_step_indices`, e.g. `["(K00367,K10534,K00372+K00360)"]` |
| `missing_kos` | str[] | every KO named in the missing steps, alternatives and complex subunits flattened. **One KO from this list may not complete a step**: a `+` complex needs all its subunits. |
| `ec_filled_step_indices` | int[] | sparse; indices of the steps counted only through EC (`size = steps_present_via_ec`) |
| `ec_genes` | str[] | sparse; locus tags of the EC-only genes that fill those steps. Each fills exactly one step of the module (§5). |

`completeness` (a ratio) is dropped: it restated two published counts and nothing read it.

With `ec_filled_step_indices` the KO-only view is fully recoverable: KO-only missing steps =
`missing_step_indices ∪ ec_filled_step_indices`.

## 5. Completeness semantics

### Evaluating a definition (KEGG Mapper convention)

- **space** separates steps; a step is covered when its expression is true.
- **comma** separates alternatives (OR).
- **parentheses** group; inside them a space means all of (AND).
- **`+`** joins subunits of one complex (all required); **`-K…`** marks an optional
  subunit (ignored).
- **`--`** is a reaction KEGG has no KO for. In today's release it only ever appears as a
  whole step (15 steps in 14 modules, never as an alternative). Such steps are **excluded**
  from `step_count` / `steps_total`, otherwise those 14 modules could never be complete. The
  MED4 cross-check against KEGG Mapper (§7) confirms Mapper does the same; if it does not, the
  rule follows Mapper.
- A definition on several lines (37 modules) is the module's separate parts, all required —
  the lines are joined as further steps (M00031 lists the LysW carrier protein on its own
  line). `step_count` counts the joined steps.
- A step that names another module (2 definitions) is covered when that module is `complete`
  for the organism. If that module needed an EC-filled step, the referencing step counts
  toward `steps_present_via_ec` (and its index goes into `ec_filled_step_indices`), so the
  KO-only count stays exact.

### Which genes count

All of an organism's gene → KO assignments, whatever their `evidence`. Filtering to stronger
evidence would leave almost nothing: KO calls are eggNOG `family_inferred`. Deleted KOs
(dropped from the KG on 2026-11-18) do not count.

### EC fallback

eggNOG leaves some genes with an EC number but no KO: 8,613 genes (2.6%). NATL2A's `pyk`
(PMN2A_0295, EC 2.7.1.40) has no KO, so on KOs alone NATL2A is 4/5 on glycolysis core and
"lacks" pyruvate kinase.

**Rule.** A gene with **no KO** and a **full** EC number (no `-`) can cover a KO term when
KEGG links that EC to the KO (`/link/ko/enzyme`), with two limits from the review:

1. **Not a complex subunit.** A KO inside a `+` complex is never filled through EC, because one
   EC often names every subunit (both nitrite-reductase subunits share EC 1.7.1.4), so one gene
   would fill the whole complex. Effect: −23 complete pairs.
2. **One step per gene per module.** An EC that KEGG links to KOs in several steps fills at
   most one of them; steps are assigned to genes by maximum bipartite matching, so the count
   is the most the EC-only genes can explain, never more. Effect: −4 complete pairs, −122 steps.

Measured: 3,468 steps filled, **6,528 → 6,987 complete** pairs (+459). Those steps are
counted in `steps_present_via_ec` and listed in `ec_filled_step_indices` / `ec_genes`.

A KO serving several steps still fills all of them (Mapper parity): that is KEGG's own
statement that one enzyme does both, unlike an EC → several-KO link.

### What `complete` does not mean

- The route runs, or the organism uses the substrate.
- The genes are in one operon or expressed together.
- An absent step is absent from the genome: the KO may be unannotated, or the genome a
  draft (OrganismTaxon assembly level / CheckM completeness is being added separately).

## 6. Where it is built

| Piece | Where |
|---|---|
| Raw downloads (`list/module`, `link/ko/module`, `link/pathway/module`, `link/ko/enzyme`, `get/M…` batched 10 per call, `br:ko00002`, `info/kegg`) | step 6, into gitignored `cache/data/kegg/raw/`, covered by the step-6 raw-input manifest |
| Module table (name, type, definition, steps, KO→step map, pathways) + class tree | step 6 → new blocks `modules` / `module_classes` in `kegg_data.json` |
| KEGG release | step 6 records `info/kegg`'s release line in `kegg_data.json`; post-import stamps `Schema_info.kegg_release` (module definitions change between releases) |
| Per-organism evaluation → `Organism_has_kegg_module`; `KeggModule` nodes and the `is_a` / KO / pathway edges | build time, new `kegg_module_adapter.py`: reads each strain's KOs + ECs from `gene_annotations_merged.json` and the step-6 blocks. Python, not post-import Cypher — the formula needs a parser. The module node set (384 with a KO + 30 reached by an organism edge) is fixed after evaluation, so no organism edge can dangle. |
| `gene_count`, `direct_gene_count`, `organism_count`, `complete_organism_count` | post-import (`.sh` + `.cypher`), declared in `schema_config.yaml` |

The evaluator is a small pure-Python module (`multiomics_kg/utils/kegg_modules.py`) with
unit tests on every operator and on known modules. It reads KOs from
`gene_annotations_merged.json` filtered to the KOs the KEGG adapter emits, so it sees exactly
the `Gene_has_kegg_ko` edge set; a KG test checks this (§7).

## 7. Validation

- Unit: the evaluator on each operator, nested groups, `--` steps, multi-line and
  module-referencing definitions (including an EC-completed referenced module); the EC
  fallback's complex and one-step-per-gene limits; the step-6 blocks on a fixture.
- **One source of KOs** (KG test): for every edge to a module whose steps are plain KO
  alternatives (no `+`, no nested AND), `steps_present - steps_present_via_ec` must equal the
  number of distinct steps reached by `(organism)<-[:Gene_belongs_to_organism]-(g)-[:Gene_has_kegg_ko]->(ko)-[r:Kegg_term_in_kegg_module]->(m)`
  over `r.steps`; for every other edge it is an upper bound. Catches the merged-JSON vs
  `Gene_has_kegg_ko` drift the review raises.
- Biology spot checks in KG tests, per named module after reading its definition:
  `complete` in every cyanobacterium for M00161 (PSII) and M00163 (PSI); not `complete` (edge
  absent or partial) in every heterotroph; the §2 TCA / nitrate / glycolysis table.
- Cross-check against KEGG Mapper's own reconstruction for MED4 (paste its KOs at
  `kegg.jp/kegg/mapper/reconstruct.html`): statuses agree except where §5's EC fallback adds
  a step; settles the `--` rule.
- `/omics-edge-snapshot` and the annotation-state baseline: no change expected (nothing
  touches BRITE or the annotation buckets).

## 8. Vocabulary and docs

- Controlled vocabularies (with `value_descriptions`): `KeggModule.module_type`,
  `KeggModule.level_kind` (4 new values, `kegg_`-namespaced per R1b), `Kegg_term_in_kegg_module.component`,
  `Organism_has_kegg_module.status` — `complete`: "every step has an annotated gene; not shown to run".
- **R3 exception #2**: `Organism_has_kegg_module.status` is a band over the edge's own
  `steps_total - steps_present`, the same shape as `tcdb_gate`. Declared in
  `docs/kg-changes/vocabulary-contract.md` as the second sanctioned exception, and CLAUDE.md's
  "one documented exception" sentence updated.
- `docs/kg-changes/kegg-modules.md` lists the new prefix use (`kegg.module:` for modules and
  class slugs) for the explorer's term-id routing, the node / edge counts, and the
  `Schema_info.kegg_release` stamp.
- CLAUDE.md key-facts bullet, relationship list, CHANGELOG `### Added` + `### Highlights`
  ("which strains have a complete route for Y").

## 9. Questions for the explorer review

1. **Tool surface.** Is `Organism_has_kegg_module` enough on its own, or do you want
   module completeness in an existing tool (`list_organisms`, `ontology_landscape`,
   `genes_by_ontology` with `ontology="kegg_module"`) or a new one (e.g.
   `module_completeness(organisms, modules | class)`)?
2. **EC fallback.** Count EC-only steps as covered (proposed, with the separate count), or
   keep completeness KO-only and report the EC steps as a hint?
3. **Status bands.** KEGG's `complete / 1 / 2 missing` plus `partial` — or only the counts
   and let the explorer label them?
4. **Missing steps.** Are `missing_steps` (formula text) and `missing_kos` both useful, or
   would you rather resolve missing KOs to names yourselves?
5. **Organisms with no step covered.** No edge (proposed) — or an explicit `absent` edge for
   modules relevant to the organism's group?
6. **Signature modules.** Same node label with `module_type = 'signature'` (proposed) — or
   leave them out?
7. **Gene-level view.** "Which genes cover step 2 of M00531 in MIT9313" is a join
   (organism's genes → KO → `steps` on the KO→module edge). Fine as a query, or worth a
   stored list on the organism edge?

(Original questions, kept so §10's answers read in context; §4–§8 above already carry the resolutions.)

---

## 10. Explorer review (2026-11-18)

Checked against the live KG (118 genomes with genes out of 123 `OrganismTaxon`; 4,808 KOs, all with
genes, vs the 4,774 quoted in §2, so re-measure before freezing) and the explorer's `ONTOLOGY_CONFIG`.
Overall: the design is good and §3's non-goals are the right ones. The blocking items are B1, B2 and B3 below.

### Blocking / must clarify

**B1. `ko00002` as a 13th BRITE tree breaks the explorer's BRITE path and the "no annotation-state change" claim.**
The explorer reads every BRITE tree as `Gene -[Gene_has_kegg_ko]-> KeggTerm(ko) -[Kegg_term_in_brite_category]-> BriteCategory`
(`ONTOLOGY_CONFIG["brite"].bridge`). The post-import `BriteCategory.gene_count` / `member_ko_count` (post-import.cypher ~1030)
and `Gene.annotation_types += 'brite'` (~798) read the same path. ko00002 is type → class → subclass → **module** → KO.
The spec has modules hang off subclasses via a new `Kegg_module_in_brite_category`. If the BRITE adapter then emits no
KO→category edges for this tree, `genes_by_ontology(ontology='brite', tree='modules')`, `ontology_landscape` and
the `BriteCategory` counts all read 0 for it, and the tree is unusable without new code. If the adapter *does* emit KO → subclass edges
(its default leaf behaviour), the counts work, but any gene whose only BRITE coverage is ko00002 gains `'brite'` in
`annotation_types`. That would be a change to the annotation-state baseline, which §3 and §7 say won't happen. Please state which one you choose. Preference: emit
**both** `KeggModule → subclass` (the module hierarchy) and the derived `KO → subclass` `Kegg_term_in_brite_category`
edges, so the existing BRITE consumers work unchanged. Then measure the `annotation_types` delta and either accept it
in §3 or exclude `tree='modules'` from the `'brite'` check.

**B2. The EC fallback makes the organism edge impossible to reproduce from the graph.**
(a) The gene-level join in Q7 (organism genes → KO → `steps`) cannot see EC-filled steps, because those genes have no KO edge.
So "which genes cover M00531 in X" will disagree with the organism edge on exactly the 3,073 steps that the
fallback fills. (b) `missing_steps` lists only the steps still missing, so a KO-only consumer can compute the KO-only
*count* (`steps_present - steps_present_via_ec`) but cannot tell *which* steps were EC-filled. That makes it impossible to
rebuild a KO-only `missing_steps`. Ask: add `ec_filled_steps` (int[], step indices) and `ec_genes` (str[], locus tags)
to the edge. They are small (EC fill ≈ 0.12 steps/edge) and close both gaps.
(c) **Complex steps:** a single EC often maps to all subunits of a `+` complex (e.g. both nitrite-reductase
subunits share an EC). One EC-only gene should not fill a whole `K00372+K00360` step. Proposed rule: the EC fallback only
fills a single-KO term; each `+` subunit needs its own gene. Please also confirm that a broad EC that maps to many KOs across *different* steps
fills only one of them, or say why that is acceptable.

**B3. `status` is a second R3 exception and needs to be declared as one.** The CLAUDE.md rule says there is *one* documented
exception (`tcdb_gate`). `status` is a band over the edge's own `steps_total - steps_present`, which is the same shape as `tcdb_gate`, so
the R3 justification in §4 is sound. Even so, it has to go into `vocabulary-contract.md` as exception #2, not just be asserted.
`completeness` is a plain ratio of two published counts with no threshold. It is harmless, but drop it unless something reads it
(the explorer would compute it anyway).

### Answers to §9

1. **Tool surface.** The organism edge alone is not enough; the explorer will add surface. Plan, explorer-side:
   - A new **`kegg_module` ontology** (the 18th), using the same pattern as BRITE: `label: KeggModule`, `gene_rel: Gene_has_kegg_ko`,
     `bridge: KeggTerm -[Kegg_term_in_kegg_module]->`, and the hierarchy through ko00002. This gives `search_ontology`,
     `genes_by_ontology`, `gene_ontology_terms`, `ontology_term_details`, and **`pathway_enrichment` / `cluster_enrichment`
     over modules** with no extra KG work. Modules are much better ORA units than maps, which is arguably the biggest analytical payoff.
     It needs `KeggModule` to carry the standard ontology props: `level`, `level_kind` (e.g. `module`), `direct_gene_count`,
     `is_uninformative` where relevant (signature modules?), and `keggModuleFullText` (already proposed).
   - A new tool, `kegg_module_completeness(organisms?, module_ids?, classes?, status?)`, reading `Organism_has_kegg_module`
     with `by_status` / `by_class` envelopes.
   - Rollups: `list_organisms` gets `complete_module_count`. `ontology_term_details(organism=…)` on a module gets that organism's
     status. These are explorer-side and need nothing new from the KG.
2. **EC fallback.** Keep it on, with B2's additions. The NATL2A `pyk` case shows it fixes real false negatives.
3. **Status bands.** Keep them (with the B3 declaration). LLM consumers read a named band far more reliably than they compute one.
   Add `value_descriptions`. For `complete`, write "every step has an annotated gene; not shown to run".
4. **Missing steps.** Keep both. Also add **`missing_step_indices` (int[])**, 1-based to match `Kegg_term_in_kegg_module.steps`,
   so the explorer can join missing steps to KO edges without parsing formula text. `missing_kos` flattens away
   complex structure, so document that one KO from it may not be enough to complete the step.
5. **No-step organisms.** No edge, as proposed. The explorer derives "absent" by a left join against all modules. An `absent` edge
   for "relevant" modules would bring in a relevance judgement that the KG shouldn't make.
6. **Signature modules.** Same label, as proposed. Suggestion: mark them `is_uninformative`, or have the explorer default
   `module_type='pathway'` on completeness and enrichment, because "complete" for a marker set reads differently. Explorer will document either way.
7. **Gene-level view.** A join is fine for KO-covered steps. The only stored list needed is `ec_genes` (B2).

### Smaller points

- **`component` is a single string, but a KO can serve several `steps`.** A KO that is required in one step and optional (`-K`) in another
  can't be represented. Make `component` a str[] parallel to `steps`, or emit one edge per (KO, step).
- **`--` placeholders.** If they count in `steps_total`, those modules can never be `complete`. Say whether Mapper excludes them
  (the MED4 cross-check in §7 will show it) and make `step_count` / `steps_total` agree with that.
- **Multi-line definitions.** Confirm `step_count` on the node counts the joined lines, the same as `steps_total` on the edge.
- **Module-referencing steps.** If the referenced module was completed via EC, does the referencing step count toward
  `steps_present_via_ec`? It should, so the KO-only count stays correct.
- **One source of KOs.** The adapter reads `gene_annotations_merged.json`, but every explorer join reads `Gene_has_kegg_ko`. After
  gene-mapping drops and remaps (e.g. the −423 genes of 2026-08-29) the two can disagree. Add a post-import check that
  recomputes KO-only `steps_present` in Cypher for a sample of edges and compares.
- **§7 photosynthesis checks.** Write them as "`complete` in every cyanobacterium" per named module after checking each
  definition: PSII / PSI modules list subunits that some *Prochlorococcus* may lack or that eggNOG misses. For heterotrophs,
  assert "not `complete`" rather than "no edge", because a shared KO can produce a 1-step partial.
- **Vocabulary / bookkeeping.** New prefix `kegg.module:`: the explorer's term-ID routing (`ontology_term_details` mixed
  batches) needs it, so please list it in the KG-changes doc. `level_kind` gains a value if modules get one (R1b namespacing). Add the KEGG
  release date to `Schema_info` / release info, because module definitions change between KEGG releases.
- The comment in the explorer's `ONTOLOGY_CONFIG["kegg"]` already expects "module rows" on `KeggTerm`. Confirm modules are a
  **separate label** (`KeggModule`), not a new `KeggTerm` `level_kind`. Separate is the better choice: it keeps `ko00001` clean.

### Explorer follow-up (once the KG lands)

The ontology is an 18th `ONTOLOGY_CONFIG` entry, using the Mode-B template from the TCDB/CAZy rollout. The new completeness tool goes through
`add-or-update-tool`. The doc touch-points are the "17 ontologies" mentions across tool docs plus a `docs://ontologies/kegg_module` page.

---

## 11. KG response to the review (2026-11-18)

All numbers in §2–§5 were re-measured with the final rules. The 4,808 vs 4,774 KO gap is the
34 KOs KEGG deleted: the review read the release graph (alpha.9), the dev graph already drops them.

| Item | Resolution | Where |
|---|---|---|
| **B1** ko00002 as BRITE | **Different from the suggested fix, same goal.** `ko00002` is not a BRITE tree at all: its type / class / subclass levels become `KeggModule` class nodes (`level` 0–2) above the modules (`level` 3), joined by `Kegg_module_is_a_kegg_module`. The explorer's `kegg_module` ontology then has its hierarchy inside one label, like TCDB / CAZy / MEROPS; the BRITE path, `BriteCategory` counts and the `'brite'` annotation type are untouched, so there is no annotation-state delta to measure or exclude. Emitting KO → subclass `Kegg_term_in_brite_category` edges as well would have put each KO in two hierarchies and made BRITE counts double-report module membership. | §3, §4 |
| **B2a/b** EC steps not reproducible | Added `ec_filled_step_indices` (int[]) and `ec_genes` (str[]); renamed the formula list `missing_step_definitions`, parallel to the new `missing_step_indices`. KO-only missing steps = `missing_step_indices ∪ ec_filled_step_indices`. | §4 |
| **B2c** EC filling a complex / several steps | Both rules adopted: no EC fill inside a `+` complex (−23 complete pairs); one step per EC-only gene per module, by bipartite matching (−4 pairs, −122 steps). Net: 6,528 → 6,987 complete (+459), 3,468 EC steps. | §5 |
| **B3** `status` is an R3 exception | Declared as exception #2 in `vocabulary-contract.md` + CLAUDE.md. `completeness` dropped. | §4, §8 |
| Answer 1: ontology props | `level`, `level_kind`, `direct_gene_count` (modules only), `keggModuleFullText`, `module_type` dense on every node. `gene_count` counts KO-bridge genes only, so it agrees with `genes_by_ontology`. | §4 |
| Answer 4: `missing_step_indices` | Added, 1-based, same numbering as `Kegg_term_in_kegg_module.steps`. `missing_kos` documented as "one KO may not be enough". | §4 |
| Answer 6: signature modules | **Not** `is_uninformative`: that presence marker means "says nothing about function", and a drug-resistance marker set does. `module_type` is dense on every node, so the explorer's default `module_type = 'pathway'` filter works at any level. | §4 |
| `component` vs several `steps` | Measured: no KO is required in one step and optional in another within a module (0 cases); scalar kept, adapter asserts the invariant and fails the build if KEGG ever breaks it. | §4 |
| `--` placeholders | Always a whole step (15 steps, 14 modules); excluded from `step_count` / `steps_total`; the MED4 Mapper cross-check confirms. | §5, §7 |
| Multi-line definitions | `step_count` counts the joined lines; equals `steps_total` on every edge. | §4, §5 |
| Module-referencing steps via EC | Counted in `steps_present_via_ec` and `ec_filled_step_indices`. | §5 |
| One source of KOs | Evaluator reads KOs filtered to the emitted KO set; KG test recomputes KO-only `steps_present` from `Gene_has_kegg_ko` (exact for plain-alternative modules, upper bound otherwise). | §6, §7 |
| Photosynthesis checks | Per named module: M00161 / M00163 `complete` in every cyanobacterium (both 1/1 in the six checked strains); heterotrophs "not `complete`". | §2, §7 |
| Prefix / release | `kegg.module:` for modules and class slugs, listed in `docs/kg-changes/kegg-modules.md`; `Schema_info.kegg_release` from `info/kegg`. | §4, §6, §8 |
| Separate label | Confirmed: `KeggModule`, not a `KeggTerm` `level_kind`. | §4 |
| Node set | New, from re-measuring: 30 modules are reached by an organism edge without any KO in the KG (EC-only steps or a module reference). They get nodes (414 total) so no edge dangles; their `gene_count` is 0 and their genes are in `ec_genes`. | §2, §6 |

**For the explorer to confirm:** `ONTOLOGY_CONFIG["kegg_module"]` as `label: KeggModule`,
`gene_rel: Gene_has_kegg_ko`, `bridge: KeggTerm -[Kegg_term_in_kegg_module]->`,
`hierarchy: Kegg_module_is_a_kegg_module`, with class nodes and modules under one label and
`level_kind` telling them apart. Enrichment should run at `level = 3` (modules) by default;
class levels are ORA-able too but, as everywhere, never mixed with modules in one test.

## 12. Explorer sign-off (2026-11-18)

**Approved.** The §11 resolutions close B1–B3. The B1 alternative is better than what the review suggested: one label with
an `is_a` hierarchy avoids double membership in BRITE and leaves the annotation state alone.

**Config confirmed.** `ONTOLOGY_CONFIG["kegg_module"] = {label: KeggModule, gene_rel: Gene_has_kegg_ko,
bridge: {node_label: KeggTerm, edge: Kegg_term_in_kegg_module}, hierarchy_rels: [Kegg_module_is_a_kegg_module]}`.
This is exactly the BRITE code path: genes enter at the bridged level-3 modules, and the `*0..` walk lifts them to the class nodes.
Explorer-side additions:
- `links_out`: KO → module (`Kegg_term_in_kegg_module`, membership) and module → map (`Kegg_module_in_kegg_pathway`).
- Enrichment docs recommend `level=3`.
- The completeness tool defaults to `module_type='pathway'`; enrichment does not filter signature modules, it only exposes `module_type`.

Minor asks. None of them block the build; they are clarifications for `kegg-modules.md`:

1. **`missing_kos` for a module-reference step.** When the missing step is `M0xxxx`, what goes in `missing_kos`: nothing,
   the referenced module's missing KOs, or the module id? Please pick one and document it. Our preference: nothing in `missing_kos`;
   `missing_step_definitions` already shows `M0xxxx`.
2. **Make `ec_genes` parallel to `ec_filled_step_indices`.** The matching already gives one gene per step, so a parallel list
   answers "which gene fills step 3" without a second lookup. A module-reference step has no gene in this module, so put the
   referenced module id in that slot. If several EC-only genes could fill a step, document that only the matched one is listed.
3. **Bipartite matching scope.** Confirm it runs only over steps the KO evidence leaves uncovered, i.e. that EC never competes with
   or replaces a KO-covered step.
4. **`organism_count` ≠ edge count.** On a module, `organism_count` (KO bridge) and the number of `Organism_has_kegg_module`
   edges differ: EC-only strains, and the 30 KO-less modules with `gene_count = 0` but real edges. Please add one line to the
   node docs so the two are not read as the same quantity.
5. **Slug ids under `kegg.module`.** Class ids don't match bioregistry's `^M\d+$`. That is fine given the TigrRole precedent, but say so
   in `kegg-modules.md`, since external CURIE validation will reject them.
6. **Still open:** re-measure the `ko00002` class-node counts on the final module set (§2).

**Explorer follow-up** (after the build lands, via `absorb-kg-rebuild` and then `add-or-update-tool`):
- the `kegg_module` ontology across the ontology and enrichment tools ("17 → 18 ontologies" doc sweep, `docs://ontologies/kegg_module`);
- a new `kegg_module_completeness` tool;
- `list_organisms.complete_module_count`;
- organism-scoped status in `ontology_term_details`;
- surfacing `Schema_info.kegg_release` in `kg_release_info`.

## 13. KG answers to the sign-off asks (2026-11-18)

All six adopted as asked; each goes into `docs/kg-changes/kegg-modules.md` and the schema descriptions.

1. **`missing_kos` for a module-reference step: nothing.** `missing_kos` lists KOs named directly
   in the missing steps; a missing `M0xxxx` step contributes no KO, and `missing_step_definitions`
   shows the module id.
2. **`ec_genes` is parallel to `ec_filled_step_indices`** (same length, same order). Each slot holds
   the locus tag the matching assigned to that step; for a module-reference step it holds the
   referenced module id (`M0xxxx`). When several EC-only genes could fill a step, only the matched
   one is listed.
3. **Matching scope: confirmed.** The bipartite matching runs only over steps the KO evidence leaves
   uncovered. EC never competes with or replaces a KO-covered step (that is how the dry run computed it).
4. **`organism_count` ≠ edge count.** Node docs will say: `organism_count` counts organisms with a
   gene whose KO is in the module (KO bridge, like `gene_count`); the number of
   `Organism_has_kegg_module` edges also counts organisms reaching the module only through EC-only
   genes or a module reference, so it can be larger, and the 30 KO-less modules have
   `organism_count = 0` with real edges. `complete_organism_count` counts edges, not the bridge.
5. **Slug ids.** `kegg-modules.md` states that class ids (`kegg.module:pathway_modules…`) are
   house-minted, do not match bioregistry's `^M\d+$`, and will fail external CURIE validation;
   only `level_kind = 'kegg_module'` ids are real KEGG module CURIEs.
6. **Class-node counts re-measured** on the 414 modules: 2 types, 12 classes, 43 subclasses
   (57 class nodes); §2 updated.

