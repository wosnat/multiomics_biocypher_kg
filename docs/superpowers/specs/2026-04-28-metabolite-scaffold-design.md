# Metabolite scaffold design

**Date:** 2026-04-28
**Status:** Phase 1 detailed; Phases 2 and 3 sketched
**Driver goal:** identify nutrient exchanges between *Prochlorococcus* and *Alteromonas*

## Goal

Add metabolites and the surrounding metabolism scaffold to the KG so that we can answer, in increasing depth:

1. Which compounds are predicted as Pro→Alt or Alt→Pro nutrient exchanges from genome content alone? (Phase 1)
2. Of those, which have actually been measured in published experiments — especially in cocultures and exometabolome samples? (Phase 2)
3. What does a query layer over both look like? (Phase 3)

Phase 1 builds the deterministic scaffold. Phase 2 layers paper measurement evidence onto it. Phase 3 adds tooling and optional augmentations.

## Architecture overview

Three independent phases, each independently shippable + testable:

| Phase | Scope | Inputs | Output |
|---|---|---|---|
| **1** | Metabolism scaffold from existing eggNOG-mapper output | MNX, TCDB, CAZy, eggNOG annotations already on disk | Metabolite + Reaction + TransporterClass + CazyFamily nodes; gene→reaction/transporter/CAZy edges; post-import organism→metabolite edges |
| **2** | Paper measurement evidence | Marquee metabolomics papers (Capovilla 2023, Kujawinski 2023, Szul 2019, Fang 2019) | MetaboliteAssay nodes + measurement edges; augments existing organism→metabolite edges |
| **3** | Query tooling + GEM-quality augmentation | (depends on Phases 1+2) | MCP tools; optional KofamScan / gapseq / Pfam-bridge / BiGG-GPR augmentations |

## Phase 1 — Metabolism scaffold from eggNOG

### Identifier backbone

**MetaNetX MNXref (release 4.5+).** The canonical Rosetta-stone for chemistry and reaction identifiers, actively maintained, sources updated 2024–2025. Cross-references KEGG, ChEBI, BiGG, MetaCyc, ModelSEED, HMDB, Reactome, Rhea, LipidMaps, SwissLipids, VMH, enviPath. We pull the TSV files (`chem_prop.tsv`, `chem_xref.tsv`, `reac_prop.tsv`, `reac_xref.tsv`).

- Metabolite IDs: `mnx:MNXM<NNNNNN>`
- Reaction IDs: `mnx:MNXR<NNNNNN>`
- Cross-refs stored as properties (`kegg_compound`, `chebi`, `bigg`, `metacyc`, `modelseed`, `hmdb`, `inchikey` for compounds; `kegg_reaction`, `rhea`, `metacyc`, `bigg`, `modelseed` for reactions).

### Pipeline changes (`prepare_data.sh`)

Existing steps are renumbered to absorb new prerequisite steps. The new step 4 (formerly step 2) is MNX-aware: it resolves raw eggNOG IDs into structured MNXR / MNXM / TC-leaf / CAZy-leaf forms during the gene-annotations merge, so adapters read pre-resolved data.

| Step | Purpose | Source / Output |
|---|---|---|
| 0 | Genome data download | NCBI, Cyanorak, UniProt, eggNOG-mapper, `gene_mapping.csv` |
| 1 | Build protein annotations per taxid | `protein_annotations.json` |
| **2 (new)** | Download MNX, TCDB, CAZy reference data | `cache/data/mnx/{chem,reac}_{prop,xref}.tsv`, `cache/data/tcdb/families.tsv`, `cache/data/cazy/families.tsv` |
| **3 (new)** | Build metabolite resolver | `cache/data/mnx/metabolite_resolver.db` (SQLite) + `metabolite_id_mapping_report.json` |
| **4 (was 2)** | Build gene annotations — MNX-aware | reads `protein_annotations.json` + eggNOG output + `metabolite_resolver.db` + TCDB/CAZy lookups; writes enriched `gene_annotations_merged.json` |
| **5 (was 3)** | Build extended gene ID mapping | unchanged |
| **6 (was 4)** | Resolve paper gene IDs | unchanged |
| **7 (was 5)** | Build OG descriptions | unchanged |
| **8 (new)** | Build pruned metabolism scaffold | reads ALL strain `gene_annotations_merged.json` to compute referenced MNXR/MNXM/TC/CAZy IDs (with parent closure), writes `cache/data/mnx/scaffold_pruned.json` |

CLI breakage: existing `--steps N M` invocations need updating to the new numbering. Documented in the commit message.

### Resolver — SQLite

```
cache/data/mnx/metabolite_resolver.db
├── compounds (mnxm_id PK, name, formula, inchikey, smiles, charge, mass)    ~340K rows
├── aliases   (source, value, mnxm_id)  INDEX(value), INDEX(source, value)   ~1.6M rows
└── names     (name_normalized, mnxm_id)  INDEX(name_normalized)              ~1M rows
```

API mirrors `gene_id_utils.resolve_row()`:

```python
mnxm_id, method = resolve_metabolite(value, conn)
# method ∈ {"xref:exact", "xref:ambiguous", "name:normalized", "ambiguous", "unresolved"}
```

Used in step 4 to resolve eggNOG `KEGG_Reaction` R-numbers into MNXR IDs, and in Phase 2 to resolve paper-reported metabolite names/IDs.

### eggNOG columns ingested in step 4

| Column | Used as | Result |
|---|---|---|
| `KEGG_Reaction` (col 15) | per-gene direct reaction list | resolved → MNXR IDs in `gene_annotations_merged.json` |
| `KEGG_TC` (col 18) | per-gene transporter classification | validated against TC hierarchy → list of TC leaf IDs |
| `CAZy` (col 19) | per-gene CAZy family | validated against CAZy → list of CAZy leaf IDs |
| `KEGG_ko`, `EC`, `BRITE`, `PFAMs`, `eggNOG_OGs`, `COG_category`, `GOs`, `Description`, `Preferred_name` | already ingested | unchanged |

`BiGG_Reaction` (col 20) is **not** ingested in Phase 1 — see "Out of scope".

`KEGG_Pathway`, `KEGG_Module`, `KEGG_rclass` are not ingested in Phase 1; can be revisited.

### New `gene_annotations_merged.json` fields per gene

```json
{
  "WP_011045886.1": {
    "...existing fields...": "...",
    "kegg_reactions":   ["MNXR101234", "MNXR101567"],
    "reactions":        ["MNXR101234", "MNXR101567"],
    "reaction_evidence": {
      "MNXR101234": ["kegg_reaction"],
      "MNXR101567": ["kegg_reaction"]
    },
    "transporter_classifications": ["1.A.1.1.1"],
    "cazy_families":               ["GH13_1"]
  }
}
```

Per-strain step-4 report: number of genes with each new annotation, resolution success rate.

### Pruned scaffold cache (step 8 output)

```
cache/data/mnx/scaffold_pruned.json
{
  "metabolites": {
    "MNXM41": { "name": "D-glucose", "formula": "C6H12O6", "inchikey": "...", "aliases": [...], ... },
    ...
  },
  "reactions": {
    "MNXR101234": {
      "name": "...",
      "equation_text": "...",
      "direction_source": "directional",
      "substrates": [["MNXM3", 1], ["MNXM41", 1]],
      "products":   [["MNXM7", 1], ["MNXM58", 1]],
      "kegg_terms": ["K00844", "K12407"],
      "ec_numbers": ["2.7.1.1", "2.7.1.2"],
      "aliases": {...}
    },
    ...
  },
  "transporter_classes": {
    "1": {"name": "Channels and Pores", "level": 0, "level_kind": "tc_class", "parent": null, ...},
    "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1", ...},
    "1.A.1.1.1": {"level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.1", "substrate_classes": ["sugar"], ...},
    ...
  },
  "cazy_families": { ... }
}
```

### Pruning rule

**Inclusion:** node appears in the KG iff at least one gene references it OR it is an ancestor of a referenced leaf. Same pattern as BRITE pruning.

- **Reaction**: any MNXR present in `kegg_reactions` of any gene, plus all reactions reachable via `Reaction_classified_by → EcNumber` for ECs already in the KG, plus all reactions reachable via `Reaction_catalyzed_by → KeggTerm` for KOs already in the KG.
- **Metabolite**: any MNXM appearing as substrate or product of any retained reaction.
- **TransporterClass**: any TC value in `transporter_classifications` of any gene, plus its ancestor classes up to the root.
- **CazyFamily**: any CAZy value in `cazy_families` of any gene, plus its ancestor classes up to the root.

Estimated pruned counts:

| Node type | Count |
|---|---|
| Metabolite | 3,000–8,000 |
| Reaction | 5,000–15,000 |
| TransporterClass | ~1,500 (full TCDB ~12K) |
| CazyFamily | ~400 (full CAZy ~700) |

### New nodes — full property lists

#### Metabolite

```
metabolite:
  is_a: chemical entity
  represented_as: node
  preferred_id: id   # mnx:MNXM<NNNNNN>
  label_in_input: metabolite
  properties:
    name: str
    formula: str
    inchikey: str
    smiles: str
    charge: int
    mass: float
    is_cofactor: str  # "true" | "false"
    kegg_compound: str
    chebi: str
    bigg: str
    metacyc: str
    modelseed: str
    hmdb: str
    aliases: str[]   # all sources collapsed for full-text
    synonyms: str[]
```

#### Reaction

```
reaction:
  is_a: biochemical reaction
  represented_as: node
  preferred_id: id   # mnx:MNXR<NNNNNN>
  label_in_input: reaction
  properties:
    name: str
    equation_text: str
    direction_source: str  # "reversible" | "directional" | "unknown"
    is_transport_reaction: str  # "true" | "false"
    kegg_reaction: str
    rhea: str
    metacyc: str
    bigg: str
    modelseed: str
    aliases: str[]
```

#### TransporterClass

```
transporter class:
  is_a: ontology term
  represented_as: node
  preferred_id: tcdb
  label_in_input: transporter_class
  properties:
    name: str
    description: str
    level: int               # 0 = class, 4 = specificity
    level_kind: str          # tc_class | tc_subclass | tc_family | tc_subfamily | tc_specificity
    substrate_classes: str[] # e.g. ["sugar", "amino acid"]
```

#### CazyFamily

```
cazy family:
  is_a: ontology term
  represented_as: node
  preferred_id: cazy
  label_in_input: cazy_family
  properties:
    name: str
    description: str
    level: int               # 0 = class, 1 = family, 2 = subfamily
    level_kind: str          # cazy_class | cazy_family | cazy_subfamily
    class: str               # GH | GT | PL | CE | AA | CBM
```

### New adapter-emitted edges

```
Gene -[catalyzes_reaction {evidence_source}]-> Reaction
Gene -[has_transporter_classification]-> TransporterClass
Gene -[in_cazy_family]-> CazyFamily

Reaction -[catalyzed_by]-> KeggTerm                          # MNX-derived KO link
Reaction -[classified_by]-> EcNumber                         # MNX-derived EC link
Reaction -[has_substrate {stoichiometry}]-> Metabolite       # MNX reac_prop.tsv
Reaction -[has_product {stoichiometry}]-> Metabolite

Transporter_class -[is_a]-> Transporter_class                # 5-level hierarchy
Cazy_family -[is_a]-> Cazy_family                            # 3-level hierarchy
```

`Gene_catalyzes_reaction.evidence_source ∈ {"kegg_reaction"}` in Phase 1. Reserved enum values for future: `bigg_reaction`, `gapseq_curated`, `gapseq_gapfilled`, `kofamscan`, `pfam_bridge`.

Edges always point at the most specific TC / CAZy leaf, never at interior class nodes (per existing convention with KEGG/BRITE).

### Post-import-derived edge

```
OrganismTaxon -[has_metabolite {properties}]-> Metabolite
```

| Property | Phase 1 | Phase 2 |
|---|---|---|
| `as_substrate_count` | int — # genes catalyzing reactions where metabolite is substrate | — |
| `as_product_count` | int — # genes catalyzing reactions where metabolite is product | — |
| `transporter_count` | int — # genes with TC for this metabolite (sparse until Phase 2/3 TC→compound) | — |
| `cazy_count` | int — # genes in CAZy families acting on this compound (sparse) | — |
| `evidence_sources` | str[] — `["kegg_reaction"]` in Phase 1 | future enum extensions |
| `via_ortholog_group_only` | — (no OG inheritance in Phase 1) | added in Phase 3 if/when OG inheritance lands |
| `measured_assay_count` | — | int (Phase 2) |
| `measured_compartments` | — | str[] (Phase 2) |
| `measured_value_kinds` | — | str[] (Phase 2) |
| `measured_paper_count` | — | int (Phase 2) |

Reaction direction is honoured permissively. MNX `direction_source = "reversible"` causes substrates and products to *both* contribute to `as_substrate_count` AND `as_product_count`. `direction_source = "directional"` is strict (substrate → as_substrate_count only, product → as_product_count only). Rationale: in vivo direction is condition-dependent, and biological reversibility is the conservative default.

### Post-import properties on existing nodes

- **Metabolite**: `producer_organism_count` (int), `consumer_organism_count` (int), `producer_genera` (str[]), `consumer_genera` (str[]), `predicted_pro_alt_exchange_flag` (str — `"true"` if both Prochlorococcus and Alteromonas are in the producer/consumer combinations), `total_gene_links` (int)
- **Reaction**: `organism_count` (int), `genera` (str[]), `gene_link_count` (int)
- **Gene** (routing signals, parallel to existing `expression_edge_count`): `reaction_count` (int), `transporter_count` (int), `cazy_count` (int)
- **EcNumber**: `reaction_count` (int)
- **KeggTerm**: `reaction_count` (int)

### Indexes

Scalar: `metabolite_inchikey_idx`, `metabolite_name_idx`, `reaction_id_idx`, `reaction_kegg_idx`, `transporter_class_level_idx`, `transporter_class_name_idx`, `cazy_family_level_idx`, `cazy_family_class_idx`.

Full-text:
- `metaboliteFullText` on `name + aliases + synonyms`
- `reactionFullText` on `equation_text + name`
- `transporterClassFullText` on `name + description`
- `cazyFamilyFullText` on `name + description`

### Adapter changes

- **`metabolism_adapter.py`** (new file). `MultiMetabolismAdapter` wrapper to match existing pattern. Reads only `cache/data/mnx/scaffold_pruned.json`. Emits Metabolite, Reaction, TransporterClass, CazyFamily nodes plus inter-scaffold edges (`Reaction_catalyzed_by`, `_classified_by`, `_has_substrate`, `_has_product`, hierarchy edges). No external file parsing.
- **`functional_annotation_adapter.py`** (extended). Three new edge generators reading new fields in `gene_annotations_merged.json`: `Gene_catalyzes_reaction`, `Gene_has_transporter_classification`, `Gene_in_cazy_family`. No MNX/TCDB/CAZy parsing in the adapter.

Both adapters are thin (file-read → property-fill → yield). Matches the existing pattern of doing resolution upstream in `prepare_data.sh` and adapters consuming pre-resolved data.

### `schema_config.yaml` changes

- **Drop** the dead `compound` node entry (lines ~526–541) and its `compound_targets_protein` edge entry.
- **Add** the new node and edge entries described above.
- Bump any version comment.

### `scripts/post-import.cypher` + `scripts/post-import.sh` changes

New blocks:

1. Compute `Organism_has_metabolite` edges with all Phase-1 properties (single `MERGE` per (org, metabolite) pair, rich property assignment via `CALL ... IN TRANSACTIONS`).
2. Compute Metabolite + Reaction + Gene + EC + KO rollup properties.
3. Add new scalar and full-text indexes listed above.
4. Reference Cypher (`scripts/post-import.cypher`) kept byte-identical to the bash version, per existing convention.

### KG validity tests

New file `tests/kg_validity/test_metabolism.py`:

- Metabolite node count ≥ minimum threshold (e.g., ≥ 2,000)
- Reaction node count ≥ minimum threshold (e.g., ≥ 3,000)
- Every Reaction has ≥1 substrate AND ≥1 product (no orphan reactions)
- Every `Gene_catalyzes_reaction` edge points to an existing Reaction
- Every `Reaction_has_substrate` / `_has_product` edge points to an existing Metabolite
- TransporterClass hierarchy is connected (no orphan classes; root `level=0` exists for each top-level class)
- CazyFamily hierarchy is connected; classes ∈ {GH, GT, PL, CE, AA, CBM}
- `Organism_has_metabolite` edges exist for every (organism, metabolite) where genes catalyse a connecting reaction
- Both Prochlorococcus and Alteromonas have `Organism_has_metabolite` edges (sanity check for the driver query)
- `predicted_pro_alt_exchange_flag = "true"` for at least N metabolites (sanity threshold for the driver goal)

### Estimated KG impact

| Quantity | Estimate |
|---|---|
| Metabolite nodes | 3,000–8,000 |
| Reaction nodes | 5,000–15,000 |
| TransporterClass nodes | ~1,500 |
| CazyFamily nodes | ~400 |
| `Gene_catalyzes_reaction` edges | ~30K–80K |
| `Gene_has_transporter_classification` edges | ~15K–40K |
| `Gene_in_cazy_family` edges | ~10K–30K |
| `Reaction_*_metabolite` edges | ~25K–80K |
| `Reaction_catalyzed_by → KeggTerm` edges | ~10K |
| `Reaction_classified_by → EcNumber` edges | ~10K |
| Hierarchy edges (TC + CAZy `is_a`) | ~2K |
| `Organism_has_metabolite` edges (post-import) | ~80K–150K |
| **Total new edges** | **~150K–300K** |

Comparable scale to existing OrthologGroup work (~223K edges, ~45K nodes).

### Driver-goal queries enabled by Phase 1 alone

**Pro↔Alt exchange candidates from genome content:**

```cypher
MATCH (pro:OrganismTaxon)-[r1:has_metabolite]->(m:Metabolite)
      <-[r2:has_metabolite]-(alt:OrganismTaxon)
WHERE pro.preferred_name STARTS WITH 'Prochlorococcus'
  AND alt.preferred_name STARTS WITH 'Alteromonas'
  AND r1.as_product_count > 0       // Pro has biosynthesis genes
  AND r2.as_substrate_count > 0     // Alt has consumption genes
RETURN m, r1.as_product_count, r2.as_substrate_count
```

**Pro genes that could produce a specific compound:**

```cypher
MATCH (m:Metabolite {name: "D-glucose"})<-[:has_product]-(r:Reaction)
      <-[:catalyzes_reaction]-(g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.preferred_name STARTS WITH 'Prochlorococcus'
RETURN o.preferred_name, collect(DISTINCT g.locus_tag), collect(DISTINCT r.name)
```

**Alt transporters for a substrate class:**

```cypher
MATCH (g:Gene)-[:has_transporter_classification]->(tc:TransporterClass)
      -[:is_a*0..]->(parent:TransporterClass)
WHERE 'sugar' IN parent.substrate_classes
  AND g.organism_name STARTS WITH 'Alteromonas'
RETURN g.locus_tag, tc.name
```

## Phase 2 — Paper measurement evidence (sketch)

Detailed design deferred until Phase 1 lands. Key shape decisions made:

- **`MetaboliteAssay`** as its own node type (separate from DerivedMetric — different terminology, different domain), but reusing every property name that maps cleanly: `name`, `value_kind`, `unit`, `field_description`, `compartment`, `metric_type`, `rankable`, `has_p_value`, `p_value_threshold`, denormalized block (`experiment_id`, `organism_name`, `publication_doi`, `omics_type`, `treatment_type`, `background_factors`, `treatment`, `light_condition`, `experimental_context`).
- **Three measurement edge types** mirroring the DerivedMetric pattern: `Assay_quantifies_metabolite` (numeric), `Assay_flags_metabolite` (boolean / presence), `Assay_classifies_metabolite` (categorical).
- **Binding edges**: `Experiment_has_metabolite_assay`, `Publication_has_metabolite_assay`, `MetaboliteAssayBelongsToOrganism`.
- **Edge properties on `Assay_quantifies_metabolite`** parallel to `Changes_expression_of`: `value`, `value_se`, `condition_label`, `time_point`, `time_point_order`, `time_point_hours`, `log2_fold_change`, `adjusted_p_value`. Post-import `rank_by_metric`, `metric_percentile`, `metric_bucket`, `significant` (when applicable, parallel to DM numeric edges).
- **Augments existing `Organism_has_metabolite` edges** with measurement properties (`measured_assay_count`, `measured_compartments`, `measured_value_kinds`, `measured_paper_count`).
- **Resolution** uses the Phase-1 `metabolite_resolver.db`. Unresolved compounds are logged in a per-paper `metabolite_resolution_report.json`. Policy for unresolved compounds (drop / keep as `unresolved_metabolite_feature` node / manual override) is decided per paper.
- **Compartment vocab** decided here once we see actual paper data; likely additions include `exometabolome` and possibly `particulate` / `lipid_fraction`.

Marquee paper candidates for the first integration (decided when Phase 2 starts):

- **Capovilla 2023** — chitin utilization, MIT9303/MIT9313 + GlcNAc, paired DE genes + intracellular metabolites. Probable first integration: paired DE+metabolites in CSV form already.
- **Kujawinski 2023** — Pro ecotype metabolite diversity. Profiling-style data, exercises the `Assay_flags_metabolite` (presence) edge type. Currently only PDF + docx supplementary; needs extraction.
- **Szul 2019** — Pro N-limitation metabolomics. Exercises the multi-condition `condition_label` property and isotopologue time-course (deferred or treated specially).
- **Fang 2019** — viral lysis DOM products. Likely exometabolome compartment. Needs extraction from PDF.

## Phase 3 — Tooling and augmentation (sketch)

- **MCP tools**: `metabolite_overview`, `predicted_exchange_candidates`, `measured_exchange_candidates`, `metabolite_response_profile`, `genes_predicted_to_produce_metabolite`, `genes_predicted_to_consume_metabolite`, `metabolites_in_compartment`.
- **OrthologGroup inheritance for reactions**: post-import `Reaction_catalyzed_by_ortholog_group` edges with `support` property (fraction of OG members with the catalyzing KO). Lifts effective gene→reaction reachability from ~52% to ~85–95% via OG bridging. Decided after Phase 1 lands — see Open questions.
- **KofamScan**: parallel KO source via KEGG-official HMM scoring. Adds ~5–15% KO coverage. New step in `prepare_data.sh`.
- **gapseq**: full genome-scale model reconstruction per strain. Output is per-strain SBML with explicit gene-reaction associations and predicted exchange reactions. Heaviest, but directly answers "predicted Pro–Alt nutrient exchanges" with a curated network. Output parsed into additional `Gene_catalyzes_reaction` edges with `evidence_source ∈ {gapseq_curated, gapseq_gapfilled}`.
- **CarveMe** (alternative to gapseq): faster ML-based GEM reconstruction, BiGG namespace, smaller compound coverage.
- **Pfam → reaction bridges**: derived from UniProt-curated Pfam ↔ Rhea links. Adds another precision tier (~78% gene coverage). Less precise than KO; included only if Phase 2 reveals coverage gaps that affect the driver query.
- **antiSMASH**: biosynthetic gene clusters + secondary-metabolite class predictions. Useful if marine NRPS/PKS compounds become a topic.
- **BiGG model GPR ingestion**: parse BiGG model JSON files to extract gene-reaction rules from cited model genes; resolves the eggNOG `BiGG_Reaction` column properly. Likely subsumed by gapseq for our purposes.
- **TC → Metabolite direct edges**: for TC entries with explicit substrate compounds in TCDB or KEGG TRANSPORTER. Supplements the `substrate_classes` text property on TransporterClass nodes with compound-level edges.

## Open questions (revisit based on data)

These are explicitly left open in the spec and revisited after Phase 1 (or later phase) lands and we observe actual data:

1. **Direct vs KO-mediated `Gene → Reaction` edges.** Phase 1 keeps both: `Gene_catalyzes_reaction` (direct, from eggNOG `KEGG_Reaction` column, ~25% gene coverage, per-gene precision) AND `Reaction_catalyzed_by → KeggTerm` (MNX-derived, ~52% gene coverage via existing `Gene_has_kegg_ko`, class-level). Trade-off: redundancy for KEGG-sourced reactions vs. precision/coverage tier flexibility. Revisit after observing query-pattern usage; may simplify to one path.
2. **OrthologGroup inheritance for reactions.** Deferred from Phase 1. Potential `Reaction_catalyzed_by_ortholog_group` post-import edge with `support` property would lift gene→reaction reachability from ~52% to ~85–95% by inheriting from KO-annotated OG members. Decision postponed until Phase 1 lands and we see how thin direct + KO-mediated coverage feels in practice.
3. **Compartment vocabulary for Phase 2.** Existing Experiment values: `whole_cell` (default), `vesicle`, `exoproteome`, `secretome`. Phase 2 likely adds `exometabolome` (dissolved extracellular). Final vocab decided once we examine actual paper data.
4. **Phase 2 marquee paper choice.** First integrated paper TBD when Phase 2 starts; Capovilla 2023 is the leading candidate (paired DE+metabolites, intracellular concentrations in fg/cell, treatment/control structure).
5. **Unresolved-name policy for Phase 2.** Per-paper decision: drop, keep as `unresolved_metabolite_feature` node, or manual override. Default proposal: drop with logging in per-paper resolution report.
6. **TransporterClass → Metabolite compound-level edges.** Phase 1 keeps TC substrate info as a text property (`substrate_classes`). Phase 2 or 3 may add direct `TransporterClass_transports → Metabolite` edges where TCDB or KEGG TRANSPORTER has explicit compound data. Sparse data; deferred.

## Out of scope (for Phase 1)

- **`BiGG_Reaction` ingestion.** The eggNOG output column contains gene references in BiGG models (`<model_id>.<gene_id>` format) rather than reaction IDs; converting to reactions requires BiGG model JSON GPR parsing. Coverage is also marginal beyond what `KEGG_Reaction` already captures (~14 BiGG-only genes for MED4, ~88 for EZ55). Listed as a Phase 3 candidate, but likely subsumed by gapseq.
- **`KEGG_Pathway` and `KEGG_Module` nodes.** Functional groupings rather than reactions; partial overlap with BRITE. May be added later if a query pattern requires them.
- **`KEGG_rclass` ingestion.** Abstract reaction patterns, not connectable to specific compounds.
- **Lipid-specific extensions.** MNX includes LipidMaps and SwissLipids, so lipid metabolites are reachable, but lipid-class hierarchies and LM-specific properties are not modelled in Phase 1.
- **Thermodynamics, charge balancing, stoichiometric coefficients beyond integer count.** Phase 1 stoichiometry is the integer coefficient from MNX. Sub-integer or charge-balanced coefficients deferred.
- **Paper measurement integration of any kind.** Strictly Phase 2.

## Acceptance criteria for Phase 1

- All Phase 1 nodes and edges present in the deployed graph at expected magnitudes (within 50% of the estimates table).
- New KG validity tests pass.
- The Pro↔Alt driver query (above) returns ≥1 result with measurement-free predicted exchange candidates.
- For at least 5 known biological exchange compounds (e.g. glycine betaine, DMSP, glucose, alanine, glycolate, ammonium), the corresponding `Metabolite` nodes exist with `predicted_pro_alt_exchange_flag = "true"`.
- `prepare_data.sh` documentation in `CLAUDE.md` updated with new step ordering.
- `omics-edge-snapshot` skill confirms no regression in existing expression-edge counts.
