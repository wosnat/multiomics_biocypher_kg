# TCDB cleanup — design

Transport programme, project 3a (of 3a TCDB cleanup → 3b gene transport layer → systems).
Status: design, approved in conversation 2026-10-19…22. Breaking changes allowed (no backwards compatibility).

> **As built (2026-10-23, `feat/tcdb-cleanup` `a46f83b5..f46c084e`).** Measured on the rebuilt
> graph: 9,389 `TcdbFamily` nodes (3,527 gene-anchored + 5,862 range-only; 134 systems unnamed
> upstream); 75,823 `Gene_has_tcdb_family` edges (eggnog 17,069 · diamond_tier 41,242 ·
> diamond_best_hit 17,512; gate direct 11,158 · annotated 19,866 · below_gate 44,799);
> substrate edges 10,778 + 367 role edges, all from `tc_system`; `'tcdb'` 14,745 → 15,800;
> **6,995** genes with a transport arm (the ~8.3K estimate below counted genes with a direct
> edge — some direct systems have no curated substrates); Metabolite 3,243 → 3,178 (65 removed,
> **31** lose `'transport'` only, **1** re-keyed: icosanoid `chebi:23899` → `chebi:62937`);
> ChemicalRole 111 → 102 (5 of the predicted 8 removed + 4 parents left childless; 3 keep their
> node); MetaboliteClass 445 → 434; pathways 452 → 449. The draft name `call_kind` shipped as
> `tcdb_call`. Consumer doc: `docs/kg-changes/tcdb-cleanup.md`; measured id list:
> `docs/kg-changes/tcdb-cleanup/removed_ids.tsv`.

## 1. Intent

The TCDB part of the KG has grown a layer of rolled-up substrates and workaround
properties (`substrate_depth`, `substrate_breadth`, `attachment_depth`,
`evidence_score`, `source_agreement`, `pfam_support`, `go_support`,
`transport_substrate_resolution`, `tcdb_evidence_score_max`, …) that the user
cannot follow and that make the obvious query wrong. A gene placed in the ABC
superfamily inherits ~550 substrates; one gene carries 2–9 overlapping edges at
different depths, some marked "superseded".

This project makes the TCDB slice **faithful to TCDB and readable**:

- the ontology is TCDB's hierarchy, with TCDB's names, descriptions and curated
  substrates, attached exactly where TCDB curates them (systems);
- each source's evidence is its own edge carrying only that source's numbers;
- nothing is rolled up, scored or reconciled across sources here.

Reconciling sources into a per-gene answer ("is this a transporter, of what
kind, how sure") is **project 3b, the gene transport layer**, which sits on top
of this cleanup and on the rest of the KG (KEGG BRITE, GO, Pfam/InterPro,
Cyanorak roles). Substrate assignment per transport system is a later
**systems** project. This document covers 3a only.

### Why now

Two analyses that tried to build transport systems from the KG
(`multiomics_analysis/analyses/2026-10-06-pro_n_import_coculture`,
`…/2026-07-06-alteromonas_coculture_carbon_sources`) spent most of their effort
working around this slice: the superfamily substrate lumping, "most_specific"
rows on superfamilies, merged multi-source edges. The recent diamond best hit
(project 2, `docs/kg-changes/tcdb-diamond-best-hit.md`) names the specific
characterised system a protein resembles, which makes a system-level model
possible.

## 2. Goals

1. `TcdbFamily` = TCDB's hierarchy (class → subclass → family → subfamily →
   system) using TCDB's standard terms and names.
2. Substrate edges leave **system nodes only**, exactly as TCDB curates them.
   Walking down from any family gives TCDB's real substrate range for it.
3. One `Gene_has_tcdb_family` edge **per source** (eggNOG, diamond tier call,
   diamond best hit), each with only its own properties. Enrichment filters by
   `tcdb_call` and `tier`.
4. Every derived/advisory TCDB property listed in §6 is removed.
5. Direction (uptake / export) on TCDB nodes where TCDB itself states it, with
   its basis.

## 3. Non-goals

- Any per-gene transporter call, confidence, subunit role, purpose or substrate
  answer (project 3b / systems).
- Reconciling the three sources, flagging their agreement or conflict (3b).
- Using KEGG, GO, Pfam, InterPro or Cyanorak to filter or score TCDB evidence (3b).
- Correcting TCDB's curation. TCDB substrates are reported as-is, including odd
  ChEBI choices (Pst → triphosphate), ChEBI roles, macromolecules and drugs.
- MCP / explorer changes (the explorer absorbs the breaking changes separately).
- Re-running diamond or changing the best-hit rule `best_hit_v1`.

## 4. Anti-goals (would make this a failure even if it ran)

- A family node that again "transports" hundreds of compounds via its own edges.
- A gene edge carrying another source's properties (sparse property soup).
- Any stored score or flag that combines sources.
- Hand-curated substrate corrections.

## 5. The model

### 5.1 Nodes — `TcdbFamily`

| Property | Notes |
|---|---|
| `name` | TCDB name (scraped, `tcdb_names.json`); falls back to `tcdb_id` |
| `tcdb_id` | bare id, e.g. `3.A.1.4.4` |
| `level` | 0–4 |
| `level_kind` | `tc_class` · `tc_subclass` · `tc_family` · `tc_subfamily` · **`tc_system`** (renamed from `tc_specificity`; "system" is TCDB's term for a 5-part TCID) |
| `superfamily` | sparse, as today |
| `description` | sparse, `tc_system` only — TCDB's curated system text from the cached family pages (`cache/data/tcdb/raw/pages/`), capped as today |
| `direction` | `uptake` · `export` · `unknown` (dense, every level; inherited downward from the node where the rule fires, so systems carry it) |
| `direction_basis` | sparse: `abc_numbering` · `tcdb_name` |
| `tcdb_substrate_count` | dense on **all five levels** (explorer R5; on a system = its own substrate count, 0 when none): distinct canonical targets (Metabolite + ChemicalRole nodes) of the `tc_system` substrate edges in this node's subtree. Computed **post-import from the edges**, so it equals the walk by construction; because §5.4 step 3 keeps every substrate-bearing system under kept families, that is TCDB's full range for the node (unresolvable substrate strings excluded, as everywhere) |
| `gene_count`, `direct_gene_count`, `organism_count` | post-import, subtree semantics as today, but counted over **gate-passing edges only** (`tcdb_gate IN ['direct','annotated']`), so term sizes match the explorer's TERM2GENE; a system reached only by near-tie or tier-3 calls reads 0 |
| `member_count` | post-import, as today |
| `tc_class_id` | as today |

Direction rule (TCDB-stated only):
1. ABC superfamily: subfamilies `3.A.1.1–3.A.1.99` are uptake systems,
   `3.A.1.100+` are exporters — a TCDB numbering convention. Applies to the
   subfamily and everything below it. Basis `abc_numbering`.
2. Otherwise, walk up from the node through its ancestors, stopping at the
   family level: the first node whose TCDB name says uptake/import/influx vs
   export/efflux/extrusion (exactly one of the two) gives the direction. Basis
   `tcdb_name`. Class and subclass nodes are never used (their names are
   mechanism classes, not direction).
3. Otherwise `unknown`. Because rules 1–2 look at the node and its ancestors,
   a descendant inherits its ancestor's direction unless its own (or a nearer
   ancestor's) name states one. Name matching is whole-word, case-insensitive.
   Computed in step 6 after names are applied (see §8 run order). **Free-text descriptions are not used** (the POC showed
   them misleading, e.g. a type-VI secretion system read as "uptake").

### 5.2 Substrate edges

`Tcdb_family_transports_metabolite` and `Tcdb_family_transports_chemical_role`
**leave `tc_system` nodes only**, one per (system, canonical target), from
TCDB's `substrates.tsv` through the existing ChEBI → canonical identity
pipeline. No `substrate_depth` property; no edges from levels 0–3.

### 5.3 Gene evidence edges — `Gene_has_tcdb_family`

One edge per (gene, TC node, call kind). Same relationship type; the adapter
**no longer merges** sources onto one edge.

| `tcdb_call` | `sources` | Target | Properties |
|---|---|---|---|
| `eggnog` | `['eggnog']` | each TC id eggNOG's `KEGG_TC` names (after `seed_aliases` remap) | — |
| `diamond_tier` | `['tcdb_diamond']` | the tier-truncated call (tier 1 → system, 2 → subfamily, 3 → family) | `tier`, `identity`, `qcov`, `scov`, `evalue` |
| `diamond_best_hit` | `['tcdb_diamond']` | each best-hit system (`best_hit_v1`, top + near-ties within 0.9×, ≤ 3) | `rank`, `tier`, `identity`, `qcov`, `scov`, `evalue`, `bitscore`, `bitscore_ratio` — each system's **own** alignment |

All edges also carry **`tcdb_gate`** — one dense, ordered string (explorer R2 +
F1; R5 string, not bool) that stamps which downstream rule the edge passes, so
consumers read rows and counts from the same set by construction:

| `tcdb_gate` | Edges | Used by |
|---|---|---|
| `direct` | `diamond_best_hit` with `rank = 1` and `tier <= 2`; `eggnog` onto a `tc_system`; `diamond_tier` with `tier = 1` | transport-arm rollups (§7) **and** everything `annotated` feeds |
| `annotated` | the remaining annotation-gate edges: other `eggnog` edges; `diamond_tier` with `tier = 2` | `'tcdb'` annotation gate, `tcdb_family_count`, TERM2GENE |
| `below_gate` | everything else: tier-3 calls, near-tie best hits, rank-1 best hits with `tier = 3` | nothing computed; evidence/clue only |

Ordered: `direct` ⊂ passes-annotation-gate = {`direct`, `annotated`}. The nesting
holds per edge (every direct arm is also an annotation arm). Computed **in the
adapter** (it has `tcdb_call`, `rank`, `tier` and the target's level), **after**
the `seed_aliases` remap — an eggNOG 5-part id re-anchored to an ancestor is
`annotated`, not `direct`. Registered in `controlled_vocabularies.yaml` as a
trust vocabulary with `value_descriptions` stating each value's exact rule, and
listed in `docs/kg-changes/vocabulary-contract.md` as a sanctioned **R3
exception** (it restates thresholds over the edge's own `tcdb_call`/`rank`/
`tier` on purpose, so rows and counts come from one set — explorer R2/F1).

**`Gene.tcdb_family_count` precisely (explorer F2):** take the TC nodes targeted by
the gene's edges with `tcdb_gate IN ['direct','annotated']`; drop any node that is
an ancestor of another node in that set, whatever the call kind; count the
distinct remaining nodes. So a `diamond_tier` edge on `3.A.1.4` plus a rank-1 best
hit on `3.A.1.4.4` = **1**; an eggNOG edge to the sibling system `3.A.1.4.5` on top
of that = **2**; eggNOG and a best hit on the same system = **1**. No per-edge
"counted" marker (that would re-create `attachment_depth`).

**Best-hit statistics (decision, user 2026-10-23).**
- `rank`: 1 = top; unique per gene (order = bitscore desc, then TCID, then
  subject — deterministic).
- `tier` on a best-hit edge is the **identity band** with the tier-call cutoffs:
  1 = identity ≥ 70, 2 = 40 ≤ identity < 70, 3 = < 40 (exactly 70 → 1, exactly
  40 → 2).
- Per-system statistics come from that system's highest-bitscore subject row
  (`select_best_hits` dedupes by TCID; take the first row it keeps). No consensus step (unlike tier calls) —
  the vocabulary entry says so. It gives one gate, `tier <= 2`, for both diamond
  kinds.
- `bitscore_ratio` = bitscore ÷ the top hit's bitscore (1.0 at rank 1, ≥ 0.9 for
  near-ties).
- Reciprocity is not stored: every best hit passed it by definition.
- Today `calls.json` `best_hit` keeps only the top hit's identity. Its block is
  extended to per-system statistics (`hits: [{tcid, identity, qcov, scov,
  evalue, bitscore}]`), regenerated from the saved diamond TSVs with
  `run_tcdb_diamond.py --from-tsv` (no diamond re-run). This is a format
  extension; the selection rule `best_hit_v1` and its output systems are
  unchanged.
- Measured (2026-10-23): among genes whose tier call is only tier 3, 619 have a
  best hit ≥ 70% identity and 2,545 at 40–70% — tier calls come from a
  consensus collapse, not identity alone, so the best-hit band is a separate,
  real signal.

- **Edge id** = `{locus_tag}-has_tcdb-{tc_id}-{tcdb_call}` (today:
  `{locus_tag}-has_tcdb-{tc_id}`, merged). Including the kind is what keeps the
  edges distinct at import; the same (gene, node) can legitimately carry up to
  three edges (e.g. eggNOG, tier-1 call and best hit all on one system) — that
  overlap is the agreement signal, read by query.
- `sources` stays a one-element list for the KG-wide trust surface (R2: joins
  `DataSource` via `data_source:` + value). It cannot tell the two diamond kinds
  apart (both come from `tcdb_diamond`), hence `tcdb_call` (R1 snake_case,
  closed vocabulary, registered in `controlled_vocabularies.yaml`).
- `evidence`: `family_inferred` for eggNOG, `homology` for both diamond kinds.
- Within one call kind, duplicates of the same (gene, TC id) are collapsed to
  one edge (e.g. two diamond candidates truncating to the same family, or two
  eggNOG ids aliasing to the same node): keep the higher-identity one's
  properties (today's adapter keeps the higher `confidence_score`; that property
  is gone).

- **Tier-3-only diamond calls are kept** as ordinary edges (decision: keep
  evidence; the transport layer decides whether to use it). `tier` lets
  enrichment exclude them.
- The six `best_hit_*` properties and `confidence_score`, `consensus_n` leave
  the diamond edge. Best-hit substrates are no longer copied onto edges: they
  are the system node's own substrate edges.
- Proteins with a best hit but no tier call (188, top HSP < 50 aa) now get a
  `diamond_best_hit` edge (they had none).
- Agreement between sources is a query (do two edges' targets share a family?),
  not a stored property.

### 5.4 Pruning

Kept TCDB node set (step 6):

1. **Seeds** = every TC id a gene edge points at: eggNOG ids, diamond tier calls
   (all tiers), **and diamond best-hit systems (new)**.
2. **+ ancestors** of the seeds (as today; no downward walk from seeds).
3. **+ every `tc_system` with ≥ 1 curated substrate under a kept family (new)**,
   with its ancestors. These nodes have `gene_count = 0`; they exist so that
   walking a family's hierarchy returns TCDB's full substrate range (§5.1
   `tcdb_substrate_count` therefore equals the walk).

Measured on the 2026-10-19 data: today 1,515 nodes (286 systems); best-hit
seeding adds ~2,012 (1,789 systems); substrate-bearing systems under kept
families not already kept add **5,249** more (re-measured 2026-10-23; the
earlier ~6.8K counted systems best-hit seeding already adds) → roughly **8.8K
TCDB nodes** plus any subfamily ancestors step 3 adds, most with no genes.
Ontology listings should filter `gene_count > 0`. Without step 3 the walk shows a
median 50% of a family's range and nothing at all for 389 of 1,081
families/subfamilies.

Seeding moves from `transporter_classification` to an explicit seed set read in
step 6 (eggNOG ∪ diamond tiers ∪ best hits); `transporter_classification` stays
a merge field if other code needs it.

Two kinds of kept node result: **gene-anchored** (steps 1–2) and **range-only**
(step 3). Nothing downstream needs the split stored — bridges cover all kept
nodes (§5.5) and `gene_count` already tells consumers which nodes have genes —
so step 6 only logs the two counts.

Metabolite / ChemicalRole nodes: the compound set's transport part becomes **the
substrates of all kept `tc_system` nodes** (gene-anchored + range-only). The
range-only systems add no new compound nodes. Compared with today it **loses the
class/subclass-only roll-up** (decision (a), user 2026-10-23): today the class
nodes union whole TCDB classes, reaching systems in families no gene touches.
The predicted ids are listed in `docs/kg-changes/tcdb-cleanup/removed_ids.tsv`
(explorer C1). Measured 2026-10-23 on `tcdb_pruned.json`: class roll-up 1,382 primaries vs
1,285 via kept families → **97 primaries + 8 ChemicalRoles** lose their transport
route. Of the 97, **65 Metabolite nodes disappear** (transport-only; eukaryote
lipids, prostaglandins, statins, ouabain, vitamin D3, phosphoinositides,
retinal, propidium iodide, collagen, dsRNA, …) and 32 stay via catalysis /
metabolomics but lose `'transport'` from `evidence_sources` (**as built: 31 lose the
tag and 1 — icosanoid `chebi:23899` — is re-keyed to `chebi:62937`**). Gene-less pathways
that lose their only support: **Phagosome, Focal adhesion, Glycosaminoglycan
biosynthesis – chondroitin sulfate / dermatan sulfate** (closes backlog
Deferred-TCDB item (f)). The pathway extension otherwise works as today (§7).
Semantics: `'transport' IN Metabolite.evidence_sources` means "TCDB lists it as
a substrate of a system in the KG" — **not** "a gene here reaches it" (gene
reach = the transport-arm rollups, §7), so a metabolite can carry the tag with
`transporter_gene_count = 0`. Say so in the vocabulary description.

`tcdb_pruned.json` keys after 3a: `kept_tcdb_ids`, **`system_substrates`** /
**`system_roles`** (`{tc_system: [canonical id, …]}`, replacing
`subtree_substrates` / `subtree_roles`; still canonicalised by the identity
pipeline's stage F and checked by `outputs.validate_build`), `seed_aliases`,
`pfam_bridge`, `go_bridge`, `chemical_roles`. Retired: `subtree_substrates`,
`subtree_roles`, `lumping_pin`, `best_hit_substrates` (and its separate resolve
`bh_rekeys`) — best-hit substrates are now the system's own `system_substrates`.

### 5.5 TCDB cross-ontology bridges

`Tcdb_family_has_pfam_domain` and the three GO bridges stay (they are TCDB's own
published maps) and are emitted for **all kept TCDB nodes, gene-anchored and
range-only alike** (decision, user 2026-10-23): they are clues for the later
system calling (e.g. "the binding protein's Pfam matches TCDB's curated members
of system X"), so they must exist on systems no gene points at yet. Rule as
today: each curated 5-part TCID in the map rolls to its nearest kept node
(itself when kept), targets shallower than `tc_family` dropped, Pfam/GO targets
pruned to existing nodes, `curated_tcids` kept on the edge.

Measured before Pfam/GO node pruning (2026-10-23 data): Pfam pairs 3,630 →
6,666, GO pairs 13,849 → 22,476 when range-only systems are included, i.e.
~1.6–1.8× the gene-anchored output; more edges now land on `tc_system` nodes
rather than rolling up to families. They are no longer used to compute anything
on gene edges. TCDB-outage fallback (bridges filtered from the previous
`tcdb_pruned.json` to the new kept set) must **warn** that range-only bridges are
missing when the previous file predates them.

## 6. Removed

| Removed | Where | Replacement |
|---|---|---|
| `Tcdb_family_transports_*` edges from levels 0–3 (subtree roll-up) | step 6 `subtree_substrates`, adapter | system edges + hierarchy walk |
| `substrate_depth` | substrate edges | not needed (systems only) |
| `TcdbFamily.substrate_breadth` (R3 exception, pinned lumping set) | node | `tcdb_substrate_count`; sunsets the R3 exception |
| `TcdbFamily.metabolite_count` | node | `tcdb_substrate_count` |
| `attachment_depth` | gene edge | query |
| `evidence_score`, `source_agreement`, `pfam_support`, `go_support` | gene edge | 3b |
| `confidence_score`, `consensus_n` | gene edge | dropped (tool-internal) |
| `best_hit_*` (6 props) | diamond edge | `diamond_best_hit` edges + system substrate edges |
| `Gene.tcdb_evidence_score_max` | gene | 3b |
| `Gene.transport_substrate_resolution` | gene | 3b / systems |
| `lumping_pin`, `subtree_substrates`, `subtree_roles` | `tcdb_pruned.json` | — |

Their `ControlledVocabulary` entries go too; the vocabulary hash moves.
**Remove by (label, property), not by name:** MEROPS shares `pfam_support`,
`confidence_score`, `consensus_n`, `evidence_score`, `best_hit_*`
(`Gene_has_merops_family`), which stay. The R3 exception for
`substrate_breadth` sunsets; `tcdb_gate` takes a new, documented one (§5.3).

New vocabulary entries: `tcdb_call`, `tcdb_gate`, `TcdbFamily.direction`,
`direction_basis`, `tcdb_substrate_count`, and numeric entries for `rank`,
`bitscore`, `bitscore_ratio`, `scov`; `tier` reworded for its two derivations
(consensus tier on `diamond_tier`, identity band on `diamond_best_hit` — never
compare across kinds); `level_kind` value `tc_specificity` → `tc_system`.

## 7. Downstream that reads this slice (must be rewired or retired)

- **Gene / annotation surface:** `annotation_types` /
  `informative_annotation_types` `'tcdb'`/`transporter` and the
  `annotation_quality` `transporter` bucket use the gate — gene has an `eggnog`
  edge, a `diamond_tier` edge with `tier <= 2`, **or a `diamond_best_hit` edge
  with `rank = 1` and `tier <= 2` (new, 2026-10-23)**. The best-hit arm makes the
  transport arm nest inside the annotated set by construction (explorer R3):
  measured 1,055 genes gain `'tcdb'` (14,745 → 15,800), and 0 genes would have a
  transport arm without the annotation. `has_any_edge` keeps counting
  `Gene_has_tcdb_family`. `Gene.tcdb_family_count` = distinct TC nodes over the
  gene's gate-passing edges (the three arms above), deepest per lineage — near-tie
  best hits never count as a second family.
- **Transport-arm metabolite rollups** — `Gene.transported_metabolite_count`,
  `Metabolite.transporter_gene_count` / `transporter_count`,
  `OrganismTaxon.transported_metabolite_count`, the transport arm of
  `Organism_has_metabolite` and `Metabolite.organism_names`/`organism_count`,
  `ChemicalRole.transporter_count` / `tcdb_family_count`. All are defined today
  through the roll-up and `attachment_depth`. **Decision (user, 2026-10-22):
  direct systems only.** A gene reaches the substrates of the `tc_system` nodes
  it points at directly — its `diamond_best_hit` edge with `rank = 1` **and
  `tier <= 2` (identity ≥ 40; decision, user 2026-10-23)**, an `eggnog` edge to a
  5-part id, or a `diamond_tier` edge with `tier = 1`. Near-tie best hits
  (`rank > 1`) and weak best hits (identity < 40) are excluded; no ancestor or
  descendant walk. Those edges carry `tcdb_gate = 'direct'` (§5.3). Every
  rollup above is redefined on that one (gene, metabolite) set, so the gene-side
  and metabolite-side counts stay two projections of it.
  Exact definitions of the system-side counts:
  - `Metabolite.transporter_count` = distinct `tc_system` nodes that carry the
    metabolite **and** are the target of ≥ 1 `tcdb_gate = 'direct'` edge
    (gene-anchored systems only; range-only systems never count).
  - `ChemicalRole.transporter_count` / `direct_transporter_count` = the same rule
    for role edges, over the role's `is_a` subtree / the role itself.
  - `ChemicalRole.tcdb_system_count` / `direct_tcdb_system_count` (**renamed**
    from `tcdb_family_count` / `direct_tcdb_family_count`, explorer C4) = the
    **TCDB-side** count: distinct kept `tc_system` nodes with a role edge
    (range-only included), over the subtree / the role itself — "how many TCDB
    systems are curated with this role". Roles are still never propagated to
    genes. Measured expectation (2026-10-23, with the
  identity gate): **~8.3K genes** keep a transport arm (**as built: 6,995** with
  `transported_metabolite_count > 0`) (was ~30K; ~13.2K without
  the gate); genes known only at family/subfamily level, or only through a weak
  best hit, read 0 until the systems project. CHANGELOG `### Breaking` must say so.
  Enrichment gene sets (explorer side) use edges with
  `tcdb_gate IN ['direct','annotated']` — the same set `TcdbFamily.gene_count`
  counts.
  Why the drop is acceptable (measured 2026-10-23, ~30.7K genes with TCDB
  evidence): 12,579 tier-3-only and 1,944 tier-2-only genes lose the arm, only
  8% / 3% transporter-named (regulators, chaperones, fold look-alikes that today
  carry a whole family's list); 2,974 eggNOG family/subfamily-only genes lose it,
  46% transporter-named. ~2,400 real transporters (e.g. MED4 phnD on `3.A.1.9`)
  therefore read 0 in the counts — **but keep their edges, and their family's
  compounds stay reachable by walking** gene → family → systems → substrates
  (complete because of §5.4 step 3). User confirmed strict (no narrow-family
  fallback in the counts). eggNOG still feeds counted compounds where it names a
  system: 2,270 genes (2,129 with curated substrates), 771 of them with no best
  hit.
- **Step 6 KEGG extension — same rule, new source (decision, user 2026-10-23).**
  Every KEGG compound among the substrates of kept `tc_system` nodes
  (`system_substrates`, gene-anchored + range-only) extends the compound set and
  pulls in its pathways, as today's rule does from `subtree_substrates`. The only
  outcome change is the measured class-roll-up loss of §5.4 (65 transport-only
  Metabolite nodes, 8 ChemicalRoles, 3 gene-less pathways: Phagosome, Focal
  adhesion, chondroitin sulfate GAG biosynthesis). Builder checks the
  before/after `kegg_data.json` diff equals exactly that list. Pathway scope is handled separately — see `plans/backlog.md`
  (KEGG pathway category filter; catalysis-only pathway creation, needs
  scoping). Measured 2026-10-23 for those items: 66 of 452 pathways have no
  gene; a direct-systems rule would remove only ~13 (Ca²⁺/Na⁺/K⁺/Cl⁻,
  amino acids and antibiotics keep the rest), so narrowing here was dropped.
- **Code that reads the roll-up or removed names (all must change):**
  `tcdb_adapter.py`; step 6 `build_kegg_metabolism_xrefs.py` (TCDB section,
  best-hit resolve); `metabolite_identity/pipeline.py` (`subtree_substrates` /
  `subtree_roles` inputs, stage O `update_lumping_pin`, the `tcdb_pruned`
  writer), `canonical.rewrite_to_canonical`, `outputs.validate_build`,
  `roles.build_subtree_roles`, `metabolite_identity/lumping.py` (+ `__init__`
  export; retired); `tcdb_utils.py`; `tcdb_diamond.py` (`_TIER_TO_LEVEL_KIND`
  writes `tc_specificity` into committed calls.json; best-hit per-system stats);
  `scrape_tcdb_names.py` (scope); `create_knowledge_graph.py` (bridge / adapter
  wiring); `scripts/post-import.sh` **and** `.cypher` (kept identical);
  `scripts/tcdb_best_hit_validation.py`; `tests/kg_validity/generate_snapshot.py`;
  `prepare_data.sh` comment; `.claude/skills/tcdb-diamond/SKILL.md`; schema;
  `controlled_vocabularies.yaml`; `vocabulary-contract.md`; CLAUDE.md.
- **Tests:** kg — `test_tcdb_cazy.py`, `test_tcdb_best_hit.py`,
  `test_annotation_trust.py`, `test_chemistry_slice1_rollups.py`,
  `test_metabolism.py`, `test_metabolite_identity.py` (+
  `fixtures/metabolite_identity_pre_counts.json`), `test_controlled_vocabularies.py`,
  snapshot; unit — `test_tcdb_adapter.py`, `test_tcdb_diamond.py`,
  `test_scrape_tcdb_names.py`, `test_controlled_vocab.py`,
  `test_controlled_vocabulary_adapter.py`, `test_mi_lumping.py` (retired),
  `test_mi_roles.py`, `test_mi_outputs.py`, `test_step6_identity_pipeline.py` /
  `step6_fixture.py`, `test_build_gene_annotations.py`,
  `test_build_kegg_metabolism_xrefs.py`.

## 8. Constraints

- House rules R1–R5 (`docs/kg-changes/vocabulary-contract.md`): `tcdb_call` and
  `direction` values are R1 snake_case; no native bool.
- String sanitization (`_clean_str`) on every string property, incl. scraped
  descriptions.
- Step 6 stays outage-tolerant for TCDB (committed `tcdb_hierarchy.json`).
- **Names and descriptions for the new systems must be generated.** Measured
  2026-10-23: none of the 5,249 range-only systems (nor best-hit systems outside
  today's scrape scope) has a name in `tcdb_names.json` or a name/description in
  `tcdb_hierarchy.json` — both cover only today's kept nodes. All of them sit
  under kept families, whose pages are already cached in
  `cache/data/tcdb/raw/pages/`, so the scraper's scope is widened to every kept
  system and re-run **from the page cache** (no network; merge-preserving). Any
  system still unnamed falls back to the TCID by design; the step-6 unnamed
  warning should report the count. **Correction (review 2026-10-23):** 112 of
  the kept families holding substrate systems (763 systems) have **no cached
  page** — families kept today only at 3-part level, which the scraper's
  `scoped_families` skips — and `fetch_page` has no offline mode. So the run
  fetches ~112 pages (~5 min at the 2.5 s cadence). The scraper reads its scope
  from the committed `tcdb_pruned.json`, and step 6 applies names when it builds
  the hierarchy, so the **run order** is: step 6 (new kept set) → scraper (new
  scope) → step 6 `--force` (names + `direction` rule 2, which needs names).

## 9. Acceptance

- No substrate edge leaves a node with `level_kind <> 'tc_system'`.
- Every `Gene_has_tcdb_family` edge has exactly one `tcdb_call` and a
  one-element `sources`; `tier` on both diamond kinds and never on `eggnog`;
  `rank`/`bitscore`/`bitscore_ratio` only on `diamond_best_hit`; exactly one
  `rank = 1` per gene with best hits; `tcdb_gate` on every edge, and its counts match the arms in §5.3. Edge ids
  are unique per (gene, TC id, call kind); `import.report` shows no
  duplicate-id drops for this type.
- No gene has a transport arm (`transported_metabolite_count > 0`) without
  `'tcdb'` in `annotation_types`.
- Edge counts by source match the merge inputs (eggNOG ids, diamond calls of
  all tiers, best-hit systems incl. near-ties and the 188 call-less proteins)
  **after** the within-kind collapse and alias collisions of §5.3; the adapter
  logs how many it collapsed.
- Every gene edge target exists (0 dangling in `import.report`).
- For every kept node (all five levels), `tcdb_substrate_count` = distinct
  substrate targets found by walking its kept subtree (a system: its own).
- Every kept `tc_system` has a non-TCID `name` except a reported remainder that
  is unnamed upstream.
- `kegg_data.json` compound and pathway sets identical before/after **except**
  the measured class-roll-up loss (§5.4): 65 transport-only Metabolite nodes,
  8 ChemicalRoles, 3 pathways; 31 compounds lose `'transport'` only and 1 is re-keyed
  (icosanoid, as built).
- MED4 spot checks: urtA → `diamond_best_hit` to `3.A.1.4.4` (Urea on the
  system); pstS best hit `3.A.1.7.6`; the ABC superfamily node has 0 own
  substrate edges and `direction` unset at its own level (`unknown`), while
  `3.A.1.4` is `uptake`/`abc_numbering`.
- Annotation-state distribution vs baseline
  (`capture_annotation_state.py --compare`): the only movement is from the
  ~1,055 genes the best-hit arm adds to `'tcdb'`; anything else is a bug.
- Removed properties absent from schema, graph, vocab and CLAUDE.md.

## 10. Decisions taken with the user

- Split into 3a (this cleanup) and 3b (gene transport layer); substrates per
  system come later.
- Report TCDB substrates as-is; no correction list.
- Keep tier-3-only diamond edges (with `tier`); filtering is 3b's / the
  analysis's job.
- Seed pruning from best hits; keep every substrate-bearing system under kept
  families (option A) so the hierarchy walk gives TCDB's full range.
- Transport-arm rollups: direct systems only, strict — no narrow-family
  fallback; family-level genes reach compounds by walking the graph (§7).
- KEGG compound and pathway sets keep today's rule; pathway scope goes to the
  backlog (§7). The class-roll-up-only compounds are accepted as lost (option a,
  2026-10-23): 65 off-target Metabolite nodes, 8 roles (as built: 9 role nodes removed, 3
  kept without a transport route), 3 gene-less pathways.
- `TcdbFamily.gene_count` and the explorer's TERM2GENE both use gate-passing
  edges.
- Best hits keep TC classes 1–4 only (no 5/8/9).
- Best-hit edges carry per-system diamond statistics and `tier` as an identity
  band; "direct" and the annotation gate use rank-1 best hits with `tier <= 2`.
- Pfam/GO bridges on all kept TCDB nodes, including range-only systems (§5.5).

Explorer dialog: `multiomics_explorer/docs/kg-specs/2026-10-23-tcdb-cleanup-explorer-asks.md`
(responses R1–R9, F1–F2 and the KG replies). From it: `call_kind` renamed
`tcdb_call`; `tcdb_gate` marker on edges (direct · annotated · below_gate);
exact `tcdb_family_count` rule; `tcdb_substrate_count` on all levels;
`direction` dense and inherited downward.

## 11. Calls made by Claude (for the user to check)

- One relationship type, one edge per call kind (`tcdb_call` + one-element
  `sources`), id includes the kind; rather than three relationship types or one
  merged edge.
- `evidence` rung per source: eggNOG `family_inferred`, both diamond kinds
  `homology`.
- `rank` on best-hit edges (near-ties kept as edges).
- `tcdb_family_count` counts gate-passing edges (incl. rank-1 best hits with
  `tier <= 2`), never near-ties.
- `tcdb_gate: direct | annotated | below_gate` as the R5-compliant form of the
  explorer's requested boolean.
- `call_kind` renamed `tcdb_call` (explorer R4 suggestion).
- Direction from names + ABC numbering only; descriptions not used.
- `tcdb_substrate_count` computed post-import from edges (canonical targets,
  metabolites + chemical roles), not in step 6.
