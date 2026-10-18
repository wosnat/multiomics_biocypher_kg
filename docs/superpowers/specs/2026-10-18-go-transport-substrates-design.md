# GO transport substrates → Metabolite edges — design

**Date:** 2026-10-18 · **Status:** draft for review · **Slice 1 of the transport programme**

## 1. Why, and where this sits

### The problem

Transport annotation is hard, and the KG makes it harder to use correctly than it
should be. The 2026-10-06 *Prochlorococcus* N-import analysis
(`multiomics_analysis/analyses/2026-10-06-pro_n_import_coculture/`, its
`gaps_and_friction.md` and `methods/notebook.md`) spent seven rule versions and 21
decisions building its own transporter and substrate rules on top of the KG. Every
hand-built rule had a flaw that only a critic pass or a QC check caught:

- non-transporters carrying TCDB substrates (`fadD`, helicases, response regulators,
  GSTs, ferritin);
- subfamily substrate lists read as gene facts (`sul1` "transports" amino acids,
  `ktrA` ammonium);
- the wrong sibling system marked confident (`cynD` reads nitrate/nitrite,
  `transport_substrate_resolution = 'resolved'`; it is the cyanate system).

### The framing

No route through the graph gives *the* transporters or *the* substrates. Every
annotation route is biased by design:

| Route | Built-in bias |
|---|---|
| TCDB family transfer (eggNOG) | inherits the subfamily's substrate list; can name the wrong sibling system |
| TCDB best hit (diamond) | one curated reference at whatever identity; no direction |
| GO substrate terms | only as good as whoever assigned the GO term; correlated with InterPro/KO naming |
| Metabolite class (KEGG BRITE) | lists *representative* compounds, so "unlisted" ≠ "not in the class" |
| Cyanorak / BRITE transporter categories | coarse class, no compound |

So the KG should publish **separate, cleanly queryable routes, each with its caveat
stated where a consumer meets it**, not fuse them into flags. The consumer — in
practice an LLM agent working through the explorer — answers a question by reading
several routes side by side.

That is less bleak than it sounds. On MED4, for the 46 genes where both GO and the
diamond best hit name a substrate, after ChEBI normalisation 25 name the same
compound, 11 a related one, and 10 disagree (mostly naming granularity). Agreement is
the confident core; the disagreements are few enough to inspect. Two caveats on that
number: the routes are **not independent** — GO terms often arrive via InterPro2GO or
eggNOG, which share homology with the TCDB calls — so agreement is weaker evidence than
it reads; and the comparison used the diamond best hit, which is not in the KG until
project 2.

### The three questions the routes must serve

1. What are the transporters, and what do they transport?
2. Are up-regulated genes enriched for a transport family?
3. What transports X (a compound or a class, e.g. amino acids), if anything — and in
   which direction?

### The programme (four projects, in order)

1. **GO transport substrates → Metabolite edges** — this document.
2. **Diamond best-hit details** — keep each gene's best-hit TCDB systems and their own
   substrate lists, which the current tier/consensus collapse discards. Stored in the
   smallest form that survives project 3 (hit records on the existing diamond edge /
   merge field, no new node types).
3. **TCDB nodes/edges redesign** — the hard one. Goal: the obvious query is the correct
   one. Likely retires `substrate_depth`, `substrate_breadth`,
   `transport_substrate_resolution` (breaking is acceptable).
4. **A transport-analysis recipe** for LLM consumers — several answers per question,
   each route's caveat, how to read agreement, what an empty answer means. Probably a
   `docs://analysis/transport` page on the explorer side.

A per-gene `Gene → Metabolite` transport edge was discussed and is **deferred until
after project 3**.

## 2. What this slice delivers

GO defines its transport terms logically. "cyanate transmembrane transporter activity"
(GO:0015110) ≡ *transmembrane transporter activity* (GO:0022857) **and** *has primary
input* (RO:0004009) **some** CHEBI:29195. These axioms ship in GO's `go-plus` release.
This slice turns them into edges from GO term nodes to the KG's `Metabolite` (and
`ChemicalRole`) nodes, so the route

```
Gene —[Gene_enables_molecular_function | Gene_involved_in_biological_process]→ GO term
     —[<GO>_transports_metabolite]→ Metabolite
```

becomes a two-hop query, joined to the same canonical Metabolite nodes that TCDB,
KEGG reactions and metabolomics use.

Measured 2026-10-18 (go-plus downloaded that day; scripts in the session scratchpad):

- 1,593 GO terms are transport terms (genus under GO:0005215 *transporter activity*
  or GO:0006810 *transport*) with a ChEBI primary input; 579 distinct ChEBI
  substrates. (A raw-id match onto existing nodes gave 386, but that undercounts:
  node `chebi_id` is often the conjugate partner of GO's id — cyanate `C01417` carries
  28024, GO says 29195 — so it is not a baseline. The step-6 identity report gives the
  real figure.) A further 42 transport axioms name a *role* through GO's `GOCHE_<n>`
  "substance with X role" classes (siderophore, xenobiotic, drug, toxin, vitamin,
  acid) — see §4 role substrates.
- Gene reach over 43 strains (TCDB, for scale: 30,547 genes), counted three ways:
  6,290 genes via ChEBI-class substrates only (first pass; 5,999 also carry TCDB, 291
  GO-only) → **6,786** once role substrates (`GOCHE_n`) are included — the figure §4
  uses before filtering → **5,436 after the §4 target filter** (2,777 MF+BP, 2,122
  BP-only, 537 MF-only; MED4 59). 5,436 is the expected-build figure.
- MED4 spot checks: cynA → cyanate; cynD → cyanate + nitrate; sul1/sul3 → sulfate;
  ktrA → K⁺; amt1 → ammonium; urtA → urea; dpp → peptide; pstS → phosphate;
  fadD, rpaA, salY → none (correct). ftn → iron ion, through the BP term "iron ion
  transport" (ferritin stores iron — what BP's broader meaning admits; see §4 MF vs BP).

## 3. Goals, non-goals, anti-goals

**Goals**
- G1. Every GO transport term in the KG that has a ChEBI primary-input axiom links to
  the canonical Metabolite (or ChemicalRole) node for that substrate.
- G2. Substrate identity is canonical: ChEBI secondary ids, conjugate acid/base forms
  and existing identity merges resolve exactly as they do for TCDB substrates
  (step 6's ChEBI/MNX machinery), so GO and TCDB meet on one node.
- G3. Dangling-proof: no edge to a GO term or Metabolite node that is not emitted.
- G4. The route's caveat is published with it (vocabulary / schema description) so an
  LLM reading the graph meets it.

**Non-goals (this slice)**
- Gene-level substrate calls, agreement flags, or any change to TCDB nodes/edges
  (projects 2–3).
- GO taxon-constraint flagging (measured: 1,912 of 1.4M gene→GO annotations violate a
  GO taxon constraint) — backlog item.
- GO molecular-function → Reaction / EC bridges via Rhea / EC / KEGG reaction xrefs —
  backlog item, unmeasured.
- GO metabolic-process input/output → Metabolite (2,646 terms; on MED4 mostly generic
  classes such as "nucleic acid", "lipid"; low value) — not planned.
- Transport direction (GO target start/end location; too sparse — 5 MED4 import
  annotations, all on one gene). **No project in the programme currently answers
  question 3's "in which direction"**; that is an open gap, not a deferral.

**Anti-goals** — it fails even if it "works" when:
- A1. GO substrates are rolled up or propagated (to ancestor GO terms, to gene
  properties, to organism reach). The TCDB subtree rollup is the complexity this
  programme is undoing; each edge states exactly one axiom.
- A2. The edge is presented as "this gene transports X". It says *this GO term is
  defined as transporting X*; whether a gene deserves the term is the gene→GO edge's
  `sources` / `evidence`.
- A3. GO and TCDB substrates are blended into one property (e.g. widening
  `Metabolite.evidence_sources = 'transport'` to mean "TCDB or GO").

## 4. Graph shape

### New edges

| Edge label | Source → target | Meaning |
|---|---|---|
| `Molecular_function_transports_metabolite` | MolecularFunction → Metabolite | GO defines this MF term's primary input as this compound |
| `Biological_process_transports_metabolite` | BiologicalProcess → Metabolite | same, for transport process terms ("cyanate transport") |
| `Molecular_function_transports_chemical_role` | MolecularFunction → ChemicalRole | the axiom's ChEBI entity is a role (e.g. xenobiotic) |
| `Biological_process_transports_chemical_role` | BiologicalProcess → ChemicalRole | same |

Naming follows `Tcdb_family_transports_metabolite` /
`Tcdb_family_transports_chemical_role`. Cellular-component terms carry no such axioms.

**Relation.** Always RO:0004009 *has primary input*, in both namespaces (measured
2026-10-18 over transport-genus axioms: 965 MF and 882 BP chemical fillers use it; no
transport term uses RO:0002313 *transports or maintains localization of*). The only
other chemical relation there, RO:0012001 *has small molecule activator* (23 MF terms,
ligand-gated channels), names an activator, not a substrate, and is excluded.

**Edge property: `go_chebi_ids`** (str[], sorted) — the ChEBI id(s) exactly as GO's
axiom states them (for role substrates, the ChEBI role id behind the `GOCHE_n` class),
before KG canonicalisation. One edge per (GO term, canonical target); a list because
canonicalisation can fold two of a term's axioms (e.g. a conjugate acid/base pair, or a
secondary id) onto one node — no id is lost and no parallel edges appear (explorer F2).
It is the only way to tell whether a GO↔TCDB disagreement comes from canonicalisation
or from the annotation (explorer TPS-004). No other properties; the edge type is the
route.

**MF vs BP.** Both are emitted, as separate edge types, because they mean different
things and the caveat must say so: MF *transporter activity* is close to "this protein
moves X"; BP *X transport* also covers storage, regulation and accessory subunits
(MED4 `ftn` reaches iron through BP "iron ion transport", GO:0006826 — ferritin stores
iron). MF-only was considered and rejected: of the 5,436 genes reaching a GO substrate
after the target filter, 2,122 do so only through BP (before the filter: 2,395 of
6,786).

### Which GO terms

Only terms that are KG GO nodes (seed annotations + ancestor closure, i.e.
`MultiGoAnnotationAdapter.all_go_terms()`), and only their **own** axioms — no
inheritance from ancestors, no propagation to descendants (A1). A gene annotated to
"L-glutamate transmembrane transporter activity" reaches L-glutamate through that
term's edge; reaching "amino acid" is the consumer's choice (walk GO `is_a`, or walk
`Metabolite_in_metabolite_class`), documented in the recipe.

Transport term = the logical definition's genus is, or descends (GO `is_a`) from,
GO:0005215 *transporter activity* or GO:0006810 *transport*; the substrate is the
RO:0004009 filler when it is a CHEBI class. Axioms whose filler is PR (protein) or GO
are skipped.

### Role substrates

GO names role-typed substrates through its own `GOCHE_<n>` classes ("substance with
siderophore role"), each defined as *chemical entity* **and** *has role* CHEBI:n — the
GOCHE number **is** the ChEBI role id, so the mapping is `GOCHE_n → chebi.role:n`
(passed through `ChebiIndex.primary` for secondary ids; no MNX, no name matching).
Measured 2026-10-18: 12 GOCHE classes on 42 transport terms; 15 of those terms are KG
GO nodes, using 5 roles — siderophore (`siderophore transport`, 290 genes direct),
xenobiotic (5 terms, up to 168 genes), vitamin, toxin, neurotransmitter
(`neurotransmitter:sodium symporter activity`, 19 genes — the bacterial NSS family;
the role name is a eukaryote-centred caveat). **All 5 already exist as `ChemicalRole`
nodes** (created by TCDB). The other 7 GOCHE roles (auxin, androgen, growth hormone,
…) sit only on GO terms outside the corpus.

Decision: role edges target **existing ChemicalRole nodes only**; this slice creates no
role nodes, so ChemicalRole ownership (TCDB side: `tcdb_pruned.json` →
`tcdb_adapter`) is unchanged. A future GO role with no node is dropped and counted in
the identity report. The emitting adapter reads the role id set from
`tcdb_pruned.json`'s `chemical_roles` (it runs before the TCDB adapter is built).

### Target filter (decided 2026-10-18)

Measured over the 6,786 corpus genes that reach a GO substrate (224 distinct targets
after merging ChEBI secondary ids and conjugate forms), the noise is two specific
kinds, and a frequency cut-off cannot separate them from signal: no target reaches
more than 18% of genes, and real substrates sit among the noise (peptide 10.1%, amino
acid 6.5% of genes, between "monoatomic cation" 17.8% and "metal cation" 8.7%).

1. **Upper-level ChEBI classes** say "moves ions / molecules", which is
   transporter-ness evidence, not a substrate call. No edge to: ion (CHEBI:24870),
   monoatomic ion (24867), monoatomic cation (23906), monoatomic anion, metal cation
   (25213), transition element cation, nitrogen molecular entity (51143; 14.3% of
   genes through 0.1% of terms — one broad BP term), biomacromolecule (33694),
   chemical entity (24431), protein (36080), inorganic anion.
2. **H⁺ and Na⁺.** They stand out on the share of *terms* (H⁺ 5.5%, Na⁺ 4.0%; next K⁺
   1.8%) as the coupling ions of every symporter/antiporter, and most genes that reach
   them as the sole input are bioenergetic machinery (ATP synthase, Na⁺-NQR,
   oxaloacetate decarboxylase, Na⁺/H⁺ and cation:H⁺ antiporters) — "moves protons" says
   nothing about nutrients. Dropping them loses "what transports Na⁺" from the GO route
   only; TCDB still answers it.

**Kept:** everything else, including K⁺ (94% of its genes reach it as the sole input —
Trk, Ktr, K⁺ channels), Cl⁻ (CLC channels, ClcA), iron (Fe²⁺/Fe³⁺, siderophore role,
iron coordination entity), Cu, Zn, Mn, Ni, Co, Mg, Ca, Hg, and the nutrient ions
(nitrate, ammonium, phosphate, sulfate, phosphonate, molybdate). Measured check on
metals: of 1,124 genes reaching an iron/metal target, ~85% are transporters (MF term
and/or a TCDB class 1–3 call with a transporter product; 15 NRAMP / FeoB genes are
caught by GO alone); the ~170 non-transporters (bacterioferritin, ferritin, Cu
resistance CopC/D, AhpC) arrive through BP terms and are what the BP caveat covers.

The filter is an explicit list in config (e.g. `config/go_transport_target_filter.yaml`,
one entry per ChEBI id with name and reason), applied in step 6 **after**
canonicalisation and matched on canonical ids (so conjugate / secondary forms of a
listed entity are caught). Dropped targets are counted per id in the identity report.

**The filter must be visible to the reader** (explorer F1): "what transports Na⁺ via
GO" returns empty by curation, not biology, and an LLM cannot otherwise tell
"filtered" from "no transporter". So post-import Group 4 stamps
**`Schema_info.go_transport_excluded_targets`** (str[], the canonical Metabolite node
ids the filter removed, sorted) — readable through `kg_release_info` — and the
kg-changes doc carries the same list with names and reasons as caveat text. The
explorer's GO arm warns when a requested target is on it. Post-import reads the list
from the committed derived cache (step 6 writes the resolved exclusion ids next to
`go_transport_substrates`), not from the raw config, so it matches what was built.

### Generic / class substrates

Classes below the filter ("amino acid", "carbohydrate", "peptide", "lipid") get exactly
TCDB's treatment: step 6 has no class/generic filter for molecule substrates — roles
are split out (`partition_role_strings`), everything else MNX can resolve becomes a
Metabolite node (generic ones carry `structure_kind = 'generic'`), and only ids MNX
cannot resolve are dropped.

### New Metabolite nodes

A GO substrate with no existing Metabolite node gets one, **only** if its GO term is a
KG GO node (gene-reachable), through the same two step-6 paths TCDB substrates use:
a KEGG primary joins the extended compound set that `build_pruned_kegg_data` prunes
to (`compounds`); a non-KEGG primary (`chebi:` / `mnx:`) goes to
`additional_compounds`. Its provenance is a new `Metabolite.evidence_sources` value
**`go_transport`** (open vocabulary; the existing `transport` keeps meaning "TCDB
substrate" — A3). A metabolite that both routes reach carries both values.

**No pathway pull-in** (decided 2026-10-18). Unlike TCDB substrates (whose KEGG
pathways extend the kept pathway set — the mechanism behind backlog item (f),
Phagosome / Focal adhesion), a compound reached only through GO does **not** add
KEGG pathways: its `Metabolite_in_pathway` edges are limited to pathways already kept
by the metabolism and TCDB arms. So this slice adds no `KeggTerm` node and changes no `reaction_count`;
`KeggTerm.metabolite_count` on already-kept pathways can rise where a GO-only compound
is in that pathway (6 such links measured 2026-10-18, ruling during implementation).

### Post-import

One new, separate rollup (decided 2026-10-18, explorer TPS-002):
**`Metabolite.go_transport_organism_names`** (str[], sorted distinct
`OrganismTaxon.preferred_name`) and **`Metabolite.go_transport_organism_count`** (int,
size invariant), over the path Gene → (MF | BP GO term) → Metabolite, organism via
`Gene_belongs_to_organism`. Dense on metabolites with a GO transport edge, absent
otherwise (sparse; the `organism_names` rule for absent vs `[]` decides which —
match it). Without it a GO-only metabolite reads `organism_count = 0` although its gene
path exists, and organism-scoped queries silently miss it.

Everything else is **unchanged**: `Organism_has_metabolite`,
`Metabolite.organism_count` / `organism_names` (still TCDB / KEGG / assay only — no
blending, A3), `Metabolite.transporter_count` / `transporter_gene_count`, and the Gene
metabolite counts. Other GO-route rollups belong to project 3 / the recipe. Both
`post-import.sh` and `post-import.cypher` get the block (identical logic, per
CLAUDE.md).

## 5. Where it is built

Follows the existing split: identity work in prepare_data step 6, emission in an
adapter, injected node sets for dangling-proofing.

1. **One GO release, pinned in prepare_data** (step 6, first sub-step, before anything
   reads GO). Today the GO node set comes from a committed `go_namespace_cache.json`
   built lazily at KG-build time from whatever `go-basic.obo` was current (March 2026),
   and `create_knowledge_graph.py --no-cache` re-downloads it mid-build — so the axiom
   map and the node set could come from different GO releases, and terms new in one
   would silently lose their edges. So GO is refreshed **as a unit, first, and on its
   own** (plan Task 0): `refresh_go_release` downloads the current `go-basic.obo` and
   `go-plus.json` (≈138 MB, `http://purl.obolibrary.org/obo/go/extensions/go-plus.json`;
   the `.obo` URL returns 403) together and accepts them only if go-basic's
   `data-version` equals go-plus's `owl#versionInfo` — dated release URLs 404 on both
   purl.obolibrary.org and release.geneontology.org (checked 2026-10-18), so a release
   cannot be fetched by date. It keeps both raw files (gitignored
   `cache/data/go_terms/raw/`), rebuilds `go_namespace_cache.json`, and writes
   `go_terms/go_release.txt`. Step 6 reads the raw go-plus and checks its release against
   `go_release.txt`; it refreshes GO only with `--refetch-raw`. The KG build stops
   refetching GO (`MultiGoAnnotationAdapter` no longer forces a rebuild on `--no-cache`)
   and warns if `go_release.txt` and `kegg_data.json`'s GO release disagree. Baseline:
   **2026-07-26** (current on 2026-10-18; the committed cache was March 2026). The
   refresh is committed, rebuilt and re-baselined **before** the GO transport work, so
   its GO count changes are reported on their own and every later gate sees additions
   only.
2. **Derived, committed cache** `cache/data/go_terms/go_transport_substrates.json`:
   `{go_id: [chebi_id, ...]}` for all transport axioms (not pruned to the corpus — a
   small file, and the corpus changes), plus a header with the go-plus release date
   and md5.
3. **Shared GO node set** (refactor, `multiomics_kg/utils/go_utils.py`). Today the
   "which GO terms become KG nodes" rule is split: `compute_ancestry_closure` and
   `load_go_data` live in `go_utils`, but seed collection (every strain's
   `go_terms`) and the namespace filter are private to `MultiGoAnnotationAdapter`
   (`_all_seed_go_ids`, `all_go_terms`, `get_nodes`). Lift the whole rule into one
   `go_utils` function — seed from all strains' `gene_annotations_merged.json`
   (given the genome config) → `compute_ancestry_closure` (every parent relation,
   including `regulates`) → namespace filter → `{go_id: schema label}` — and make
   the adapter's `all_go_terms()` / `get_nodes()` and step 6 both call it. The
   adapter's emitted node set must not change (assert in a test against the current
   behaviour before switching).
4. **Canonicalisation** (step 6). For the GO terms in that shared set, resolve each ChEBI id
   through the same `ChebiIndex` + MNX path TCDB substrates use, partition roles with
   the same `ChebiIndex.is_role`, add missing compounds to `additional_compounds` with
   `evidence_sources: ['go_transport']`, and write the resolved map into
   `kegg_data.json` as a new top-level block `go_transport_substrates:
   {go_id: [metabolite_or_role_node_id, ...]}`. The identity report records counts:
   axioms seen, in-corpus, resolved, new nodes, roles, unresolved (with ids).
5. **Emission** (adapter). `MultiGoAnnotationAdapter` (it owns GO nodes and knows each
   term's label) emits the four edge types from `kegg_data.json`'s block, pruned to
   its own `all_go_terms()` and to injected `metabolite_node_ids` /
   `chemical_role_node_ids` (`None` → emit nothing, the `pfam_node_ids` contract).
   `create_knowledge_graph.py` already builds the metabolism adapter first; inject
   `metabolism_adapter.all_metabolite_ids()`, and the role ids read from
   `tcdb_pruned.json` (§4 role substrates).

   *Builder's choice:* a small dedicated adapter instead, if wiring into
   `MultiGoAnnotationAdapter` is awkward — the contract (inputs, pruning, `None`
   behaviour) is what matters.

## 6. Validation

- **Unit tests:** axiom parsing (genus filter incl. `is_a` descent; RO:0004009 only;
  PR/GO fillers skipped); ChEBI secondary-id + conjugate resolution onto existing
  nodes (urea 9888 → urea's canonical node; ammonium 28938 ↔ the ammonia/ammonium
  node; cyanate); role partition; dangling-proof pruning; `None` injection → no edges;
  `_clean_str` on any string; the `go_transport` value lands only on GO-reached
  compounds.
- **Before/after:** `/omics-edge-snapshot` and
  `tests/kg_validity/capture_annotation_state.py --save/--compare` (expected: no
  change); `Metabolite` node count and `evidence_sources` distribution before/after
  (report the new-node count); pathway `KeggTerm` count and every existing pathway's
  `metabolite_count` / `reaction_count` unchanged (no pathway pull-in).
- **KG validity (`-m kg`):** new edge types present; every edge's endpoints exist;
  spot checks — GO:0015110 → cyanate (`kegg.compound:C01417`), GO:0015204 urea
  transporter → urea, GO:0015112 → nitrate; `fadD` (PMM0402) reaches no metabolite
  through this route; MED4 genes reaching ≥ 1 metabolite via the route within a
  stated band of the measured 65.
- **Snapshot:** add a GO→Metabolite edge sample to `snapshot_data.json`.
- **Target filter:** unit tests that each listed id (and a conjugate / secondary form of
  one) is dropped, and that K⁺, Cl⁻, Fe²⁺, nitrate, ammonium, phosphate survive; a
  `-m kg` check that no GO transport edge targets a filtered id or H⁺ / Na⁺, and that
  `Schema_info.go_transport_excluded_targets` equals the set of filtered canonical ids.
- **`go_chebi_ids`:** a unit case where two axioms of one term canonicalise to one node
  yields one edge with both ids; `-m kg` check that every GO transport edge has a
  non-empty `go_chebi_ids`.
- **Organism rollup:** `go_transport_organism_count = size(go_transport_organism_names)`;
  for a sample metabolite it equals the distinct organisms of genes carrying a GO term
  with an edge to it.
- **Report** (identity report every build, and in the kg-changes doc): edges per type,
  distinct metabolites, new Metabolite nodes, generic/role share, genes and organisms
  reached, unresolved ChEBI ids, dropped-by-filter counts, and the **target
  distribution table** (per canonical target: gene share, term share, strain count) so
  a new strain that shifts it gets noticed.

## 7. Vocabulary and docs

- Release: `### Added` + `mcp_min_version` raised to the explorer release that
  surfaces the GO route (§8); stamp the go-plus **release date recorded in the derived
  cache header** (not the download date) in the kg-changes doc.
- `config/controlled_vocabularies.yaml`: add `go_transport` to the
  `Metabolite.evidence_sources` description (open vocab; description-only edit plus
  the new value — tell the explorer, it moves `controlled_vocabularies_hash`); register
  `Metabolite.go_transport_organism_names` / `_count` if the coverage test requires it.
- `config/go_transport_target_filter.yaml` (new) — the §4 target filter, each entry with
  its reason.
- `schema_config.yaml` / post-import Group 4: `Schema_info.go_transport_excluded_targets`
  (str[]).
- `schema_config.yaml`: the four edge types, with a comment carrying the caveat (A2).
- New `docs/kg-changes/go-transport-substrates.md` (what / why / caveat / measured
  numbers / query examples).
- `CLAUDE.md`: Key-graph-facts bullet; Data Locations entry for
  `go_transport_substrates.json` and raw `go-plus.json`; relationships list; step-6
  description gains the sub-step.
- `CHANGELOG.md` `### Added` (+ a `### Highlights` line: "which GO-annotated genes move
  compound X" is now a graph query).
- `plans/backlog.md`: file the two follow-ups (GO taxon-constraint flag; GO MF ↔
  Rhea/EC/KEGG-reaction bridge) and mark the transport programme's four projects.

## 8. Joint release with the explorer

Dialog: `multiomics_explorer/docs/kg-specs/2026-10-18-go-transport-substrates-explorer-asks.md`.
Slice 1 ships as a joint release (explorer accepted 2026-10-18): the explorer adds
`go_transport` as a third `evidence_source` on `genes_by_metabolite` /
`metabolites_by_gene` (rows carry the GO term id + name, its namespace, and the gene→GO
edge's `evidence` / `sources`), a `go_transport` value in `gene_overview`'s per-gene
`evidence_sources`, and role inputs (`chebi.role:`) reaching the role edges; the route
caveats live in `kg_schema` edge-type descriptions plus tool docs until the explorer's
`docs://analysis/transport` recipe page lands. Order: KG dev build → explorer absorbs +
surfaces → explorer MCP release → KG release with `mcp_min_version` raised to it (§7).

The caveat text is a KG deliverable: the kg-changes doc carries, ready to copy, one
paragraph per edge type (MF vs BP meaning; "term definition ≠ gene fact"; the target
filter with its list). `kg_schema` has no edge-type descriptions today, so the explorer
adds a hand-authored edge-descriptions block (next to its `schema_baseline.yaml`) filled
from that text — explorer work, no KG graph change (explorer F4).

## 9. Left to the builder

- Exact sub-step placement inside step 6 and the identity-report field names.
- Name and signature of the shared `go_utils` function (§5.3).
- Adapter placement (§5.5).
- Test fixture shape (a tiny go-plus JSON fixture is expected).
