# TCDB cleanup (transport programme, project 3a)

**Date:** 2026-10-23 · **Branch:** `feat/tcdb-cleanup` (commits `a46f83b5..f46c084e`, docs in the commit after)
**Spec:** [`docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md`](../superpowers/specs/2026-10-22-tcdb-cleanup-design.md)
**Plan:** `docs/superpowers/plans/2026-10-23-tcdb-cleanup.md`
**Explorer dialog:** `multiomics_explorer/docs/kg-specs/2026-10-23-tcdb-cleanup-explorer-asks.md`
**Removed / changed ids:** [`tcdb-cleanup/removed_ids.tsv`](tcdb-cleanup/removed_ids.tsv) (measured on the rebuilt graph)
**Breaking:** yes. No backwards-compatibility aliases.

All numbers below were measured on the rebuilt graph (Docker build of `c129ffc8`, kg suite
1,346 passed / 0 failed after `f46c084e`).

## In one paragraph

The TCDB slice now mirrors TCDB. `TcdbFamily` is TCDB's hierarchy (class → subclass → family →
subfamily → system) with TCDB's names. Substrate edges leave **system nodes only**, as TCDB
curates them, so a family no longer "transports" hundreds of compounds through its own edges. To
get a family's range, walk down to its systems. Each source's evidence is its **own**
`Gene_has_tcdb_family` edge, which carries only that source's numbers. Every scored, rolled-up or
reconciled TCDB property is gone. Reconciling the sources into a per-gene transporter call is
project 3b.

## Nodes — `TcdbFamily`

| | Before | After |
|---|---|---|
| Nodes | 1,515 | **9,389**: 3,527 gene-anchored + 5,862 range-only |
| By level | class 7 · subclass 34 · family 592 · subfamily 596 · specificity 286 | class 7 · subclass 34 · family 592 · subfamily 1,432 · **system 7,324** |
| With `gene_count > 0` | all | 2,190 (class 7 · subclass 31 · family 375 · subfamily 687 · system 1,090) |

- **`level_kind` `tc_specificity` → `tc_system`.** "System" is TCDB's term for a 5-part TCID.
- **Kept set.** It holds every TC id a gene edge points at (eggNOG ids, diamond tier calls of all
  tiers, and diamond best-hit systems, which are new seeds), plus their ancestors. On top of that
  it holds every `tc_system` with at least one curated substrate under a kept family, plus its
  ancestors (**range-only**). Range-only nodes have `gene_count = 0`. They exist so the walk from
  a family returns TCDB's full substrate range. **Ontology listings should filter
  `gene_count > 0`.**
- **Names.** The scraper's scope widened to every kept system: 592 of 592 family pages, 9,622
  names. 134 kept systems have no upstream name and fall back to the TCID (`t.name = t.tcdb_id`).
- **New `direction`** (dense: `uptake` 956 · `export` 1,367 · `unknown` 7,066) and sparse
  **`direction_basis`**:
  - `abc_numbering` (414 uptake + 510 export): in the ABC superfamily, subfamilies `3.A.1.1–99`
    are uptake and `3.A.1.100+` are export, by TCDB's numbering convention.
  - `tcdb_name` (542 + 857): the nearest node, going up to family level, whose TCDB name says
    import/uptake/influx or export/efflux/extrusion. Only one of the two may appear. Matching is
    whole-word and includes agent nouns ("Exporter").
  - ABC numbering applies to the subfamily and everything below it and is never overridden.
    Otherwise the nearest named node wins: the node's own name first, then its ancestors up to
    family level. So a system named as an uptake porter inside an exporter family reads uptake;
    the direction is not simply inherited downward. Class/subclass names and free-text
    descriptions are never used.
  - Known false positive: `2.A.86.1.7` reads uptake from its descriptive name ("may function
    with a heme uptake transporter") although its family is an exporter family. No genes carry
    it.
  - `3.A.1` itself is `unknown`; `3.A.1.4` is `uptake` / `abc_numbering`.
- **New `tcdb_substrate_count`** (post-import, dense, all five levels): distinct Metabolite +
  ChemicalRole targets of the system substrate edges in the node's subtree, i.e. TCDB's range for
  the node. On a system it is the system's own count (0 when none). Medians: class 96, subclass 5,
  family 1, subfamily 1, system 1. `3.A.1` = 551.
- **`gene_count` / `direct_gene_count` / `organism_count`** keep subtree semantics but count
  **gate-passing edges only** (`tcdb_gate IN ['direct','annotated']`), the same set as the
  explorer's TERM2GENE. A system reached only through near-tie or tier-3 calls reads 0.
- **Removed:** `metabolite_count` (use `tcdb_substrate_count`) and `substrate_breadth`. The pinned
  lumping set and its R3 exception go with it.

## Substrate edges

| | Before | After |
|---|---|---|
| `Tcdb_family_transports_metabolite` | 10,790, from every level (subtree roll-up) | **10,778, from `tc_system` only** |
| `Tcdb_family_transports_chemical_role` | 392 | **367, from `tc_system` only** |
| `substrate_depth` | on both | **removed** |

There is one edge per (system, canonical target), straight from TCDB's `substrates.tsv` through
the ChEBI → canonical identity pipeline. Odd curation (e.g. Pst → triphosphate, roles,
macromolecules) is reported as-is. The ABC superfamily node has 0 own substrate edges. MED4's
urea system `3.A.1.4.4` ("The high-affinity (<1 μM) urea porter") carries urea.

## Gene edges — `Gene_has_tcdb_family`

There is now **one edge per (gene, TC node, `tcdb_call`)**. The adapter no longer merges sources.
Edge id: `{locus_tag}-has_tcdb-{tc_id}-{tcdb_call}`. In total there are **75,823 edges** over
30,725 genes.

| `tcdb_call` | Edges | `sources` / `evidence` | Target | Properties |
|---|---|---|---|---|
| `eggnog` | 17,069 | `['eggnog']` / `family_inferred` | each eggNOG `KEGG_TC` id after the `seed_aliases` remap (114 re-anchored to an ancestor) | — |
| `diamond_tier` | 41,242 | `['tcdb_diamond']` / `homology` | tier-truncated call: tier 1 → system, 2 → subfamily, 3 → family | `tier`, `identity`, `qcov`, `scov`, `evalue` |
| `diamond_best_hit` | 17,512 | `['tcdb_diamond']` / `homology` | each `best_hit_v1` system (top + near-ties ≥ 0.9 × top bitscore, ≤ 3) | `rank`, `tier`, `identity`, `qcov`, `scov`, `evalue`, `bitscore`, `bitscore_ratio`, each from that system's **own** alignment |

- **`rank`**: 1 = top hit, unique per gene.
- **`tier` has two derivations.** On `diamond_tier` it is the consensus tier. On
  `diamond_best_hit` it is the identity band (≥ 70 → 1, 40–70 → 2, < 40 → 3). Never compare tiers
  across kinds.
- **Units.** `identity`, `qcov` and `scov` are percent. `qcov` used to be mis-declared as a
  fraction.
- **Call-less proteins.** The ~188 proteins with a best hit but no tier call (top HSP < 50 aa)
  now get a `diamond_best_hit` edge.
- **Duplicates.** Within one call kind, duplicate (gene, node) pairs collapse to the
  higher-identity edge.
- **Agreement.** The same (gene, node) can carry up to three edges. That overlap is the agreement
  signal, read by query; nothing about it is stored.

**`tcdb_gate`** (dense and ordered) stamps which downstream rule each edge feeds. Every consumer
reads rows and counts from one set:

| `tcdb_gate` | Edges | Rule | Feeds |
|---|---|---|---|
| `direct` | 11,158 (best hit 7,059 · eggNOG 3,907 · tier 192) | `diamond_best_hit` with `rank = 1` and `tier <= 2`; `eggnog` onto a `tc_system`; `diamond_tier` with `tier = 1`. Always onto a system, after the alias remap | the transport-arm metabolite rollups, plus everything `annotated` feeds |
| `annotated` | 19,866 (eggNOG 13,162 · tier 6,704) | other `eggnog` edges; `diamond_tier` with `tier = 2` | the `'tcdb'` / `transporter` annotation gate, `Gene.tcdb_assignment_count`, `TcdbFamily.gene_count`, TERM2GENE |
| `below_gate` | 44,799 (tier 34,346 · best hit 10,453) | everything else: tier-3 calls, near-tie best hits, rank-1 tier-3 best hits | nothing computed; evidence only |

Filter with `r.tcdb_gate IN ['direct','annotated']`; do not re-derive the rule. `tcdb_gate` is
the sanctioned R3 exception (`vocabulary-contract.md`).

**Removed from the edge:** `evidence_score`, `source_agreement`, `pfam_support`, `go_support`,
`attachment_depth`, `confidence_score`, `consensus_n`, and the six `best_hit_*` properties. The
best-hit substrates are now the system node's own substrate edges. The MEROPS edge keeps its
same-named properties.

**Example — MED4 cynD (PMM0372):**

| `tcdb_call` | Target | `tcdb_gate` | Statistics |
|---|---|---|---|
| `diamond_best_hit` | `3.A.1.16.2` "Bispecific cyanate/nitrite transporter." | `direct` | rank 1, tier 2, 59.2% |
| `diamond_tier` | `3.A.1` | `below_gate` | tier 3 |
| `eggnog` | `3.A.1.16.1` "Four component nitrate/nitrite porter." | `direct` | — |

Its counts: `tcdb_assignment_count = 2` (two sibling systems), `transported_metabolite_count = 3`
(cyanate + nitrite from `3.A.1.16.2`, nitrate + nitrite from `3.A.1.16.1`). MED4 urtA (PMM0970) has a rank-1 best hit
on `3.A.1.4.4` (`direct`), a tier-2 call on `3.A.1.4` (`annotated`), and eggNOG edges on
`3.A.1.4.4` / `3.A.1.4.5` (`direct`). MED4 pstS (PMM0710) has its best hit on `3.A.1.7.6`.

## Gene / metabolite / organism rollups

- **Annotation gate.** `'tcdb'` in `annotation_types` and `transporter` in
  `informative_annotation_types` mean "has a gate-passing edge". The rank-1 best-hit arm is new,
  so these go from 14,745 to **15,800 genes (+1,055)**. `annotation_state` moved for 8 genes
  (`informative_single` → `informative_multi`); nothing else moved (`capture_annotation_state.py
  --compare`).
- **`Gene.tcdb_assignment_count`** takes the TC nodes targeted by the gate-passing edges, drops any
  node that is an ancestor of another node in that set (whatever the call kind), and counts the
  rest. A tier-2 call on `3.A.1.4` plus a rank-1 best hit on `3.A.1.4.4` = 1; near-tie best hits
  never add an assignment. 15,800 genes have it, max 9.
  - *Renamed after 3a (2026-10-24, `feat/post-3a-fixes`):* this shipped in 3a as
    `Gene.tcdb_family_count`. The name misled — the counted targets are mostly systems and
    subfamilies, not families — so it is now `tcdb_assignment_count`. Same meaning, same
    values; no alias. `ChemicalRole.tcdb_system_count` is unrelated.
- **Transport arm = direct systems only.** A gene reaches the substrates of the `tc_system` nodes
  it points at through `tcdb_gate = 'direct'` edges, with no ancestor or descendant walk.
  - **6,995 genes** have `transported_metabolite_count > 0` (was ~30K through the family
    roll-up; p90 3, max 20). The design estimate of ~8.3K counted genes with a direct edge, but
    some direct systems have no curated substrates.
  - Genes known only at family/subfamily level, or only through a near-tie or identity < 40 best
    hit, now read 0. They keep their edges, and their family's compounds stay reachable by
    walking gene → family → systems → substrates.
  - `Metabolite.transporter_gene_count`, `OrganismTaxon.transported_metabolite_count`, the
    transport arm of `Organism_has_metabolite` and `Metabolite.organism_names` / `organism_count`
    are projections of the same set.
- **`Metabolite.transporter_count`** = distinct systems that carry the metabolite and are the
  target of at least one `direct` edge. Range-only systems never count. 460 metabolites, median 1,
  max 88.
- **`ChemicalRole`.** `transporter_count` / `direct_transporter_count` follow the same rule.
  `tcdb_family_count` / `direct_tcdb_family_count` are **renamed** `tcdb_system_count` /
  `direct_tcdb_system_count`: kept systems curated with the role, range-only included.
- **`Organism_has_metabolite`: 88,579 → 52,648.**
- **Removed from Gene:** `tcdb_evidence_score_max` and `transport_substrate_resolution`.
- **`'transport' IN Metabolite.evidence_sources`** means "TCDB lists it as a substrate of a kept
  system", not "a gene here reaches it". It is on 1,285 metabolites, and 825 of them have
  `transporter_gene_count = 0`. Route on `transporter_gene_count > 0`.

## Metabolites, roles, classes and pathways that left

Substrates now come only from kept systems, so the old class/subclass roll-up is gone. That
roll-up reached systems in families no gene touches. The full list is in
[`removed_ids.tsv`](tcdb-cleanup/removed_ids.tsv):

| `change` | Rows | Notes |
|---|---|---|
| `metabolite_removed` | 65 | transport-only, off-target for marine bacteria: prostaglandins, statins, ouabain (`kegg.compound:C01443`), vitamin D3, phosphoinositides, retinal, propidium iodide, collagen, double-stranded RNA, D-threonine, … Metabolite 3,243 → 3,178 |
| `metabolite_rekeyed` | 1 | icosanoid `chebi:23899` → `chebi:62937` (icosanoid anion): the conjugate group's canonical member flips once TCDB stops reaching 23899, and the GO route keeps 62937. A rename, not a loss |
| `metabolite_loses_transport_tag` | 31 | stay as nodes via catalysis / metabolomics / GO, without `'transport'` |
| `chemical_role_removed` | 9 | 5 lost their only route + 4 parents left with no children. ChemicalRole 111 → 102 |
| `chemical_role_loses_transport_route` | 3 | hormone, ionophore, nutrient stay as nodes |
| `metabolite_class_removed` | 11 | classes whose only members were removed. MetaboliteClass 445 → 434 |
| `pathway_removed` | 3 | `ko04145` Phagosome, `ko04510` Focal adhesion, `ko00532` chondroitin / dermatan sulfate GAG biosynthesis (all gene-less). Pathways 452 → 449 |

Related edge counts:
- `Metabolite_in_pathway` 9,348 → 9,241.
- `Metabolite_in_metabolite_class` 1,153 → 1,117.

Side effect: in Biller 2022, the assay row "Threonine" changes `stereo_assignment` from
`dominant_form_assumed` to `as_reported`. It still resolves to `kegg.compound:C00188`. The flag
depends on whether stereo siblings exist as KG nodes, and D-threonine is one of the removed nodes.
See the backlog.

## Bridges

`Tcdb_family_has_pfam_domain` (5,772), `Tcdb_family_involved_in_biological_process` (5,780),
`Tcdb_family_enables_molecular_function` (5,545) and `Tcdb_family_located_in_cellular_component`
(3,451) are now emitted for **all kept nodes, range-only systems included**. Before, they covered
gene-anchored nodes only (1,856 / 2,588 / 2,266 / 1,708). More of them land on `tc_system` nodes
instead of rolling up to families. They are clues for the systems project, and nothing on gene
edges is computed from them any more. The rule is unchanged: nearest kept node, targets shallower
than `tc_family` dropped, `curated_tcids` on the edge.

## Unchanged

- `Changes_expression_of`: 323,849 edges.
- No other node or edge type.
- The diamond best-hit rule `best_hit_v1` and its selected systems. The calls.json `best_hit`
  block gained a per-system `hits` list, regenerated with `run_tcdb_diamond.py --from-tsv`, with
  no diamond re-run.

## Vocabulary

| Change | Entries |
|---|---|
| New | `Gene_has_tcdb_family.tcdb_call`, `.tcdb_gate`, `.rank`, `.bitscore`, `.bitscore_ratio`, `.scov`; `TcdbFamily.direction`, `.direction_basis`, `.tcdb_substrate_count` |
| Removed, by (label, property) | the properties listed above |
| Reworded | `tier` (two derivations); `TcdbFamily.level_kind` value `tc_system` |
| Renamed | `ChemicalRole` `tcdb_system_count` / `direct_tcdb_system_count` |

`controlled_vocabularies_hash` moves. See the TCDB cleanup addendum in
[`vocabulary-contract.md`](vocabulary-contract.md).

## `tcdb_pruned.json`

- **Keys now:** `kept_tcdb_ids`, `system_substrates`, `system_roles`, `chemical_roles`,
  `seed_aliases`, `pfam_bridge`, `go_bridge`.
- **Retired:** `subtree_substrates`, `subtree_roles`, `lumping_pin`, `best_hit_substrates`.
