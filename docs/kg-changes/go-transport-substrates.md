# GO transport substrates: GO term to Metabolite / ChemicalRole edges

**Date:** 2026-10-18
**Spec:** [`docs/superpowers/specs/2026-10-18-go-transport-substrates-design.md`](../superpowers/specs/2026-10-18-go-transport-substrates-design.md)
**Plan:** `docs/superpowers/plans/2026-10-18-go-transport-substrates.md`
**Track:** ontology to ontology bridge, slice 1 of the transport programme
**GO release:** `2026-07-26` (go-plus `owl#versionInfo` = go-basic `data-version`, recorded in `cache/data/go_terms/go_release.txt`)

## What changed

GO defines its transport terms logically: "cyanate transmembrane transporter activity"
(GO:0015110) is *transmembrane transporter activity* **and** *has primary input* (RO:0004009)
**some** CHEBI:29195. The axioms ship in GO's `go-plus` release. They are now edges from GO
term nodes to the same canonical `Metabolite` / `ChemicalRole` nodes that TCDB, KEGG
reactions and metabolomics use, so

```
Gene -[Gene_enables_molecular_function | Gene_involved_in_biological_process]-> GO term
     -[<GO>_transports_metabolite]-> Metabolite
```

is a two-hop query. Nothing about TCDB nodes or edges, gene properties or existing
organism rollups changed (spec anti-goals A1 to A3: no propagation, no gene-level claim, no
blending with the TCDB `transport` evidence value).

Also changed, and logged on its own in CHANGELOG `### Changed`: GO itself was refreshed from the
March 2026 release to 2026-07-26 before this work (MF nodes 2,899 to 2,890, BP 3,433 to 3,322, CC
unchanged at 451; `Gene_enables_molecular_function` 363,298 to 362,439,
`Gene_involved_in_biological_process` 546,971 to 533,186, `Gene_located_in_cellular_component`
154,722 to 155,459). GO is now refreshed as a unit (go-basic + go-plus, release-checked) by
`refresh_go_release` and by step 6 only with `--refetch-raw`; the KG build no longer refetches GO and
`--no-cache` does not touch it.

## The four edge types

| Edge | Source to target | Edges | Distinct targets |
|---|---|---|---|
| `Molecular_function_transports_metabolite` | MolecularFunction to Metabolite | 257 | 171 |
| `Biological_process_transports_metabolite` | BiologicalProcess to Metabolite | 246 | 182 |
| `Molecular_function_transports_chemical_role` | MolecularFunction to ChemicalRole | 7 | 5 |
| `Biological_process_transports_chemical_role` | BiologicalProcess to ChemicalRole | 8 | 5 |
| **Total** | | **518** | 214 (208 metabolites + 6 roles) |

The relation is always RO:0004009 *has primary input*. One edge per (GO term, canonical target).
Single edge property **`go_chebi_ids`** (str[], sorted): the ChEBI id(s) exactly as GO's axiom states
them (for roles, the ChEBI role id behind the `GOCHE_<n>` class), before canonicalisation. It is a
list because canonicalisation can fold two axioms of one term (a conjugate acid/base pair, a
secondary id) onto one node; it is how to tell whether a GO vs TCDB disagreement comes from
canonicalisation or from the annotation.

Only terms that are KG GO nodes get edges, and only their **own** axioms: no inheritance from
ancestors, no propagation to descendants. Transport term = the logical definition's genus is, or
descends from, GO:0005215 *transporter activity* or GO:0006810 *transport*; PR and GO fillers and the
RO:0012001 *has small molecule activator* relation are skipped.

Built from the committed `cache/data/go_terms/go_transport_substrates.json` (1,624 transport terms
with substrates, 1,847 GO-ChEBI pairs in go-plus, not pruned to the corpus). Of those, 514 terms
(560 axioms) are KG GO nodes. Step 6 resolves them through the same ChEBI + MNX path as TCDB
substrates: 0 unresolved ChEBI ids, 0 not in ChEBI.

## Target filter (13 ids, 88 axioms dropped)

`config/go_transport_target_filter.yaml`, applied in step 6 **after** canonicalisation on canonical
ids (so conjugate and secondary forms of a listed entity are caught). A frequency cut-off cannot
separate the noise from real substrates (peptide 12.6% of genes, between "monoatomic cation" and
"metal cation"), so the list is explicit.

| Target | Canonical id | Axioms dropped | Reason |
|---|---|---|---|
| ion | `chebi:24870` | 1 | upper-level ChEBI class: transporter-ness, not a substrate |
| monoatomic ion | `chebi:24867` | 7 | same |
| monoatomic cation | `chebi:23906` | 14 | same |
| monoatomic anion | `chebi:23905` | 6 | same |
| metal cation | `chebi:25213` | 4 | same |
| transition element cation | `chebi:33515` | 2 | same |
| nitrogen molecular entity | `chebi:51143` | 1 | same (14.3% of genes through one broad BP term) |
| biomacromolecule | `chebi:33694` | 1 | same |
| chemical entity | `chebi:24431` | 5 | same |
| protein | `chebi:36080` | 1 | same |
| inorganic anion | `chebi:24834` | 1 | same |
| hydron (H+) | `kegg.compound:C00080` | 23 | coupling ion of every symporter/antiporter and of bioenergetic machinery (ATP synthase, Na+-NQR); "moves protons" says nothing about nutrients |
| sodium(1+) | `kegg.compound:C01330` | 22 | same |

Kept on purpose: K+ (94% of its genes reach it as the sole input: Trk, Ktr, K+ channels), Cl-, iron
(Fe2+/Fe3+, siderophore role), Cu, Zn, Mn, Ni, Co, Mg, Ca, Hg, and the nutrient ions (nitrate,
ammonium, phosphate, sulfate, phosphonate, molybdate). The list is published on
`Schema_info.go_transport_excluded_targets` (str[], the 13 canonical ids; also written to
`go_transport_excluded_targets.json` and `kegg_data.json`) so a reader can tell "filtered by
curation" from "no transporter exists". Querying what transports Na+ or H+ through GO therefore
returns empty by curation, not biology; TCDB still answers it.

## Target distribution (top 25 of 214, after the filter)

From `cache/data/metabolomics/metabolite_identity_report.json` `.go_transport.distribution`; genes
counted over 43 strains (5,436 genes reach at least one target in total). Gene share is the share of
those genes.

| Target | Name | Genes | Gene share | GO terms | Strains |
|---|---|---|---|---|---|
| `kegg.compound:C00012` | Peptide | 686 | 12.6% | 5 | 43 |
| `chebi.role:35703` | xenobiotic (role) | 443 | 8.2% | 5 | 27 |
| `chebi:35238` | amino-acid zwitterion | 440 | 8.1% | 9 | 39 |
| `chebi:16646` | carbohydrate | 374 | 6.9% | 7 | 43 |
| `kegg.compound:C00238` | Potassium cation | 367 | 6.8% | 12 | 43 |
| `chebi.role:26672` | siderophore (role) | 305 | 5.6% | 2 | 15 |
| `chebi:24875` | iron cation | 282 | 5.2% | 5 | 43 |
| `kegg.compound:C00009` | Orthophosphate | 244 | 4.5% | 5 | 43 |
| `kegg.compound:C01356` | Lipid | 216 | 4.0% | 3 | 29 |
| `kegg.compound:C00338` | Lipopolysaccharide | 204 | 3.8% | 4 | 43 |
| `kegg.compound:C00032` | Heme | 168 | 3.1% | 6 | 43 |
| `chebi:22918` | branched-chain amino acid | 163 | 3.0% | 2 | 40 |
| `chebi:29067` | carboxylic acid anion | 142 | 2.6% | 4 | 19 |
| `kegg.compound:C00086` | Urea | 140 | 2.6% | 4 | 27 |
| `chebi:64709` | organic acid | 137 | 2.5% | 2 | 23 |
| `chebi:63299` | carbohydrate derivative | 135 | 2.5% | 2 | 14 |
| `chebi.role:23888` | drug (role) | 120 | 2.2% | 1 | 18 |
| `kegg.compound:C00098` | Oligopeptide | 119 | 2.2% | 5 | 31 |
| `chebi:33892` | iron coordination entity | 115 | 2.1% | 1 | 43 |
| `chebi:23378` | copper cation | 111 | 2.0% | 4 | 41 |
| `kegg.compound:C00059` | Sulfate | 103 | 1.9% | 5 | 38 |
| `kegg.compound:C00151` | L-Amino acid | 96 | 1.8% | 4 | 27 |
| `kegg.compound:C00698` | Cl- | 96 | 1.8% | 6 | 43 |
| `kegg.compound:C00014` | Ammonia | 83 | 1.5% | 2 | 43 |
| `kegg.compound:C00420` | Polysaccharide | 79 | 1.4% | 2 | 31 |

The report is regenerated by every step-6 run, so a new strain that shifts the distribution shows up
as a diff in that file.

## New Metabolite properties and nodes

- **`Metabolite.evidence_sources`** gains the value **`go_transport`** (open vocabulary). It means "a
  GO transport term the KG carries names this compound as its primary input". It is separate from
  `transport` (TCDB substrate); a metabolite both routes reach carries both.
- **+29 Metabolite nodes** (3,214 to 3,243): GO substrates with no existing node, with
  `evidence_sources` containing `go_transport` (KEGG primaries via the compound set, non-KEGG via
  `additional_compounds`). Side effects of the 29: `Metabolite_in_pathway` +6, `Metabolite_is_stereo_form_of`
  +2, `Metabolite_is_enantiomer_of` +1, `Metabolite_in_metabolite_class` +1. **No new `KeggTerm`**: a
  GO-only compound never pulls pathways in. The 6 new pathway links land in already-kept pathways
  (C00855 D-methionine to ko00470, ko01100, ko02010; C21066 to ko02010; C21382 to ko02024, ko02026),
  so those pathways' `metabolite_count` can rise by a few; pathway `KeggTerm` count stays 452 and no
  `reaction_count` changes.
- **Antimonite** stays `chebi:30297`: GO brings in its conjugate acid `chebi:49870` (antimonous acid)
  and the canonical rule would have renamed the existing node and 14 TCDB substrate lists, so the
  canonical id is pinned in `config/metabolite_identity_overrides.yaml` (chloride / hydrochloric acid
  precedent). `tcdb_pruned.json` is byte-identical to before.
- **`Metabolite.go_transport_organism_names`** (str[], sorted distinct `OrganismTaxon.preferred_name`)
  and **`go_transport_organism_count`** (int, `= size(names)`): organisms with a gene annotated to a
  GO term that has an edge to this metabolite. Sparse: present on the 198 metabolites that some gene reaches through a
  DIRECT gene-to-GO annotation (count min 1, median 10, max 43); 10 of the 208 metabolites with a GO
  transport edge (e.g. fructose, disaccharide, oligosaccharide, alcohol, icosanoid) have edges only from
  GO terms present as ancestors with no directly annotated gene and carry no organism properties. Absent
  therefore means "no directly annotated gene path", not "no GO transport edge": to test for an edge, test the edge. `Metabolite.organism_count` /
  `organism_names`, `Organism_has_metabolite`, `transporter_count` and the Gene metabolite counts are
  unchanged (still TCDB / KEGG / assay only: no blending).

**The rollup uses direct gene-to-GO annotations only**, with no `is_a` propagation (anti-goal A1).
That is why 208 metabolite targets exist but only 198 carry the organism properties: 10 targets are
reached only by ancestor GO terms that no gene carries directly. To include inheritance, walk GO
`is_a` in the query.

- **Gene reach** (genes reaching at least one metabolite through the route, after the filter):
  KT2440 390, DSS-3 332, Marinobacter (MarRef v6) 252, AltDE 200, W3-18-1 198; MED4 59 (TCDB for
  scale: 30,547 genes over the corpus). Of 5,436 genes (43 strains): 2,777 reach a substrate through
  both MF and BP, 2,122 through BP only, 537 through MF only (the reason BP edges are kept: MF-only
  was rejected).

## Agreement with TCDB

Measured on MED4 (design section 1): of 46 genes where both GO and the TCDB diamond best hit name a
substrate, after ChEBI normalisation 25 name the same compound, 11 a related one, 10 disagree (mostly
naming granularity). Two caveats: the routes are not independent (GO terms often arrive through
InterPro2GO or eggNOG, which share homology with the TCDB calls), so agreement is weaker evidence than
it reads; and the comparison used the best hit, which is not in the KG yet (programme project 2).

## Caveat text (ready to copy into `kg_schema` edge descriptions)

**`Molecular_function_transports_metabolite`**
> GO defines this molecular-function term as a transporter activity whose primary input
> (RO:0004009) is this compound. Read it as close to "a gene product with this term moves X". It is
> a term definition, not a gene fact: whether a gene deserves the term is the gene-to-GO edge's
> `sources` / `evidence`. Only the term's own axiom is stated (no ancestor inheritance, no
> propagation). Upper-level ChEBI classes (ion, monoatomic cation, chemical entity, ...), H+ and Na+
> are deliberately not targets (list in `Schema_info.go_transport_excluded_targets`): an empty
> answer for those is curation, not biology. `go_chebi_ids` holds the ChEBI id(s) as GO states them.

**`Biological_process_transports_metabolite`**
> GO defines this biological-process term (for example "iron ion transport", GO:0006826) as
> involving the movement of this compound. BP terms are broader than MF terms: they also cover
> storage, regulation and accessory subunits, so a gene carrying one is not necessarily a
> transporter (MED4 `ftn`, ferritin, reaches iron through "iron ion transport"). It is a term
> definition, not a gene fact: whether a gene deserves the term is the gene-to-GO edge's `sources` /
> `evidence`. Only the term's own axiom is stated. Filtered targets (H+, Na+, upper-level ion and
> molecule classes) return empty by curation, not biology.

**`Molecular_function_transports_chemical_role`**
> GO defines this transporter-activity term as moving something GO calls `<role>` (for example
> siderophore, xenobiotic, drug, toxin, vitamin). The target is a `ChemicalRole`; roles are not
> propagated to metabolites, so this does not say which compound is moved. A term definition, not a
> gene fact: check the gene-to-GO edge's `sources` / `evidence`. The neurotransmitter role names a
> eukaryote-centred class: in bacteria the matching terms (neurotransmitter:sodium symporter activity)
> are mostly amino-acid transporters.

**`Biological_process_transports_chemical_role`**
> GO defines this transport-process term as involving something GO calls `<role>` (siderophore,
> xenobiotic, drug, toxin, vitamin, ...). The broader BP meaning applies (storage, regulation,
> accessory subunits are admitted). Roles are not propagated to metabolites. A term definition, not a
> gene fact: check the gene-to-GO edge's `sources` / `evidence`.

## Query examples

```cypher
// Which genes carry a GO term defined as moving urea (GO route), with the gene->GO evidence?
MATCH (m:Metabolite {id: 'kegg.compound:C00086'})<-[t:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]-(go)
MATCH (g:Gene)-[r:Gene_enables_molecular_function|Gene_involved_in_biological_process]->(go)
RETURN g.locus_tag, go.id, go.name, type(t) AS route, r.evidence, r.sources;

// Class reach, route 1: walk GO is_a (term inherits from descendants), e.g. amino-acid transporters
MATCH (go:MolecularFunction)-[:Molecular_function_is_a_molecular_function*0..]->(p:MolecularFunction)
WHERE p.name = 'amino acid transmembrane transporter activity'
MATCH (go)-[:Molecular_function_transports_metabolite]->(m:Metabolite)
RETURN DISTINCT m.id, m.name;

// Class reach, route 2: walk the chemical class instead
// Amino acids (br08001); members attach at child classes, so walk the class subtree; coverage is partial
MATCH (c:MetaboliteClass {id: 'kegg.brite:br08001.A5.B2'})<-[:Metabolite_class_is_a_metabolite_class*0..]-(d)<-[:Metabolite_in_metabolite_class]-(m:Metabolite)
MATCH (m)<-[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]-(go)
RETURN m.id, count(DISTINCT go) AS go_terms;

// Which organisms carry the GO route for a metabolite
MATCH (m:Metabolite {id: 'chebi:24875'}) RETURN m.go_transport_organism_names;
```

## Out of scope (see `plans/backlog.md`, transport programme)

- **Direction**: GO target start/end location is too sparse (5 MED4 import annotations, one gene).
  No project in the programme currently answers question 3's "in which direction".
- **Taxon constraints**: 1,912 of 1.4M gene-to-GO annotations violate a GO taxon constraint; not flagged.
- **GO MF to Rhea / EC / KEGG reaction xrefs**: unmeasured bridge.
- **GO metabolic-process input/output** to Metabolite: not planned (mostly generic classes).
- **Obsolete gene GO annotations**: ~330,594 gene `go_terms` in the merged annotations reference ids
  not live in 2026-07-26 (139K obsolete with no replacement, 145K with a single `replaced_by`, 45K
  `consider`-only); 23,462 alt-id annotations are remapped to the live primary id, the rest are dropped
  at build as before. Pre-existing; the source annotations (eggNOG, UniProt, Cyanorak) carry old ids.

## Vocabulary and release notes

Adding `go_transport` to the `Metabolite.evidence_sources` description moves
`Schema_info.controlled_vocabularies_hash`; consumers pinning the hash must refresh. The KG release
that carries this slice sets `[tool.release-kg] mcp_min_version` to the explorer release that surfaces
the GO arm (joint release; done at `/release-kg` time).
