# Metabolism Chemistry Layer (Reaction + Metabolite nodes)

**Added:** 2026-04-30 (Phase 1.2 + 1.2.1 + 1.2.2)

## What changed

Added `Reaction` and `Metabolite` nodes plus 5 new edge types, ingesting the gene-reachable subset of KEGG reactions and compounds enriched with MetaNetX (MNX), Rhea, ChEBI, and HMDB cross-references. This makes the chemistry layer queryable: which genes catalyze which reactions, which metabolites participate in those reactions, which pathways contain them, and which organisms can produce or consume each metabolite.

The KEGG primary IDs are the canonical identifiers (`kegg.reaction:R*`, `kegg.compound:C*`); MNX is used as a deduplication hub and cross-reference enrichment source — not as the primary key. Reaction direction is **not modelled** — `Reaction_has_metabolite` is direction-agnostic (substrates and products are not distinguished).

## New node types

### Reaction (~2,349 nodes)

Primary ID: `kegg.reaction:R*` (e.g. `kegg.reaction:R00200` = pyruvate kinase).

| Property | Type | Notes |
|---|---|---|
| `kegg_reaction_id` | str | Raw KEGG R-number, e.g. `"R00200"` |
| `name` | str | KEGG reaction name |
| `ec_numbers` | str[] | EC numbers from MNX `classifs`, e.g. `["2.7.1.40"]` |
| `kegg_pathway_ids` | str[] | KO-prefixed pathway IDs, e.g. `["ko00010"]` |
| `mnxr_id` | str (nullable) | MNX cross-reference ID; ~93/2349 reactions lack one |
| `rhea_ids` | str[] | Rhea cross-references (sorted) |
| `reaction_class` | str | `"chemical"` \| `"transport"` (currently all `chemical` in deployed graph) |
| `mass_balance` | str | `"balanced"` \| `"unbalanced"` (1,922 balanced / 427 unbalanced) |
| `gene_count` | int | post-import: distinct genes catalyzing |
| `organism_count` | int | post-import |
| `organisms` | str[] | post-import: distinct `organism_name` values |

### Metabolite (~2,188 nodes)

Primary ID: `kegg.compound:C*` (with `chebi:` / `mnx:` fallback for compounds without a KEGG entry; in the deployed graph **100% have a `kegg.compound:` primary**).

| Property | Type | Notes |
|---|---|---|
| `kegg_compound_id` | str (nullable) | Raw KEGG C-number; nullable for ChEBI/MNX-only metabolites |
| `name` | str | Primary KEGG compound name (first synonym before `; `) |
| `formula` | str | Chemical formula from MNX |
| `mass` | float | Monoisotopic mass |
| `inchikey` | str | InChI key |
| `smiles` | str | SMILES string |
| `mnxm_id` | str (nullable) | MNX cross-reference (100% populated in deployed graph) |
| `chebi_id` | str (nullable) | ChEBI ID (lowest-sorted alias when multiple); ~1,887/2,188 populated |
| `hmdb_id` | str (nullable) | HMDB ID; ~1,175/2,188 populated |
| `gene_count` | int | post-import |
| `organism_count` | int | post-import |

## New edge types

| Edge type | Source | Target | Count | Direction-agnostic? |
|---|---|---|---|---|
| `Gene_catalyzes_reaction` | Gene | Reaction | ~52,742 | n/a |
| `Reaction_has_metabolite` | Reaction | Metabolite | ~10,050 | **YES** — substrate vs product not distinguished |
| `Reaction_in_kegg_pathway` | Reaction | KeggTerm (`kegg.pathway:ko*`) | ~6,349 | n/a |
| `Metabolite_in_pathway` | Metabolite | KeggTerm (`kegg.pathway:ko*`) | ~8,095 | n/a |
| `Organism_has_metabolite` | OrganismTaxon | Metabolite | ~37,010 | n/a |

`Organism_has_metabolite` is **materialized post-import** as a 2-hop save: `(Organism)<-[:Gene_belongs_to_organism]-(Gene)-[:Gene_catalyzes_reaction]->(Reaction)-[:Reaction_has_metabolite]->(Metabolite)`. This makes the multi-evidence-source pattern (Phase 2 will add paper-derived metabolomics evidence) queryable by a single edge type.

`Metabolite_in_pathway` and `Reaction_in_kegg_pathway` reuse the existing `KeggTerm` pathway nodes (`kegg.pathway:ko*`) — no new pathway node type was introduced.

### Pathway pruning rule (Option B)

A pathway gets `Metabolite_in_pathway` / `Reaction_in_kegg_pathway` edges **only if the pathway is also reachable from KOs or reactions tied to genes in the KG**. Pathways that have only compound annotations (e.g. Leishmaniasis, MAPK signalling) are dropped — they're biologically inappropriate for marine bacteria. KG validity test [tests/kg_validity/test_metabolism.py](tests/kg_validity/test_metabolism.py#L122-L136) enforces this invariant.

## New properties on existing nodes

### OrganismTaxon (post-import)

| Property | Type | Source |
|---|---|---|
| `reaction_count` | int | distinct reactions catalyzed by genes belonging to this organism |
| `metabolite_count` | int | distinct metabolites linked via `Organism_has_metabolite` |

Top organisms by metabolic coverage in the deployed graph: *Pseudomonas putida* KT2440 (1,449 reactions / 1,490 metabolites), *Ruegeria pomeroyi* DSS-3 (1,377 / 1,468), *Alteromonas macleodii* EZ55 (1,348 / 1,428).

## New indexes

| Index | Type | On |
|---|---|---|
| `reaction_id_idx` | scalar | `Reaction.id` |
| `reaction_kegg_id_idx` | scalar | `Reaction.kegg_reaction_id` |
| `reaction_mnxr_idx` | scalar | `Reaction.mnxr_id` |
| `metabolite_id_idx` | scalar | `Metabolite.id` |
| `metabolite_kegg_id_idx` | scalar | `Metabolite.kegg_compound_id` |
| `metabolite_mnxm_idx` | scalar | `Metabolite.mnxm_id` |
| `metabolite_chebi_idx` | scalar | `Metabolite.chebi_id` |
| `reactionFullText` | full-text | `Reaction.name` |
| `metaboliteFullText` | full-text | `Metabolite.name`, `Metabolite.formula` |

## Pipeline changes

A new step 6 in `prepare_data.sh` runs [multiomics_kg/download/build_kegg_metabolism_xrefs.py](multiomics_kg/download/build_kegg_metabolism_xrefs.py): walks every strain's `gene_annotations_merged.json` to identify gene-reachable {KOs, reactions, compounds, pathways}, prunes raw KEGG to that subset, and enriches reactions/compounds with MNX/Rhea/ChEBI/HMDB cross-refs. Output: a single `cache/data/kegg/kegg_data.json` (~3-4 MB, indented JSON for git-friendly diffs).

Both the build-time `metabolism_adapter` and the `kegg_annotation_adapter` are pure file readers of `kegg_data.json` — no SQLite or KEGG REST calls at build time. The 2.6 GB MNX SQLite resolver is opened only at prepare-data time, mirroring the step-5 `build_og_descriptions.py` pattern. The MNX resolver itself is built by the standalone `scripts/refresh_mnx.sh` (rerun only when MNX releases).

## Implementation files

- [multiomics_kg/download/build_kegg_metabolism_xrefs.py](multiomics_kg/download/build_kegg_metabolism_xrefs.py) — step 6: prune-then-enrich
- [multiomics_kg/adapters/metabolism_adapter.py](multiomics_kg/adapters/metabolism_adapter.py) — `MultiMetabolismAdapter`
- [multiomics_kg/utils/kegg_utils.py](multiomics_kg/utils/kegg_utils.py) — KEGG endpoint parsers
- [multiomics_kg/utils/metabolite_utils.py](multiomics_kg/utils/metabolite_utils.py) — `mnxm_to_primary_id()`, `mnxr_to_primary_id()`
- [config/schema_config.yaml](config/schema_config.yaml#L637-L671) — `reaction` + `metabolite` schema entries (lines 637-671)
- [config/schema_config.yaml](config/schema_config.yaml#L1108-L1146) — 5 association entries (lines 1108-1146)
- [scripts/post-import.sh](scripts/post-import.sh#L137-L148) — indexes
- [scripts/post-import.sh](scripts/post-import.sh#L723-L762) — rollups + `Organism_has_metabolite` materialization
- [tests/test_build_kegg_metabolism_xrefs.py](tests/test_build_kegg_metabolism_xrefs.py) — step 6 unit tests
- [tests/test_metabolism_adapter.py](tests/test_metabolism_adapter.py) — adapter unit tests
- [tests/kg_validity/test_metabolism.py](tests/kg_validity/test_metabolism.py) — KG validity (counts, primary-ID coverage, Option B pruning)

## Querying

```cypher
// All reactions a gene catalyzes
MATCH (g:Gene {locus_tag: 'PMM0807'})-[:Gene_catalyzes_reaction]->(r:Reaction)
RETURN r.kegg_reaction_id, r.name, r.ec_numbers, r.mass_balance

// All metabolites in glycolysis (ko00010) for MED4
MATCH (g:Gene {organism_name: 'Prochlorococcus MED4'})
      -[:Gene_catalyzes_reaction]->(r:Reaction)
      -[:Reaction_in_kegg_pathway]->(p:KeggTerm {id: 'kegg.pathway:ko00010'})
MATCH (r)-[:Reaction_has_metabolite]->(m:Metabolite)
RETURN DISTINCT m.kegg_compound_id, m.name, m.formula

// Find which organisms can make/consume a given metabolite
MATCH (m:Metabolite {kegg_compound_id: 'C00031'})  // D-glucose
      <-[:Organism_has_metabolite]-(o:OrganismTaxon)
RETURN o.preferred_name, o.metabolite_count, o.reaction_count
ORDER BY o.metabolite_count DESC

// Cross-reference: KEGG → MNX → ChEBI for ATP
MATCH (m:Metabolite {kegg_compound_id: 'C00002'})
RETURN m.name, m.mnxm_id, m.chebi_id, m.hmdb_id, m.formula, m.mass

// Pyruvate kinase spot-check (Spec 1.2 acceptance test)
MATCH (r:Reaction {id: 'kegg.reaction:R00200'})-[:Reaction_has_metabolite]->(m)
RETURN r.name, r.ec_numbers, collect(m.name) AS metabolites
```

## What is *not* in this layer

- **Reaction direction** — substrates vs products are not distinguished; queries traverse `Reaction_has_metabolite` symmetrically.
- **Per-paper metabolite measurements** — paper-derived metabolomics integration (peak intensities, fold changes, etc.) is Phase 2 and will add paper-evidenced edges onto the same Metabolite nodes.
- **OrthologGroup-mediated reaction inheritance** — Phase 3 will let an OG inherit reactions from its members so genes without direct KEGG annotation can still surface in reaction queries.
- **Transporter (TCDB) and CAZy family nodes** — Spec 1.3 will add `TransporterClass` + `CazyFamily` nodes; the cached TCDB tables and CAZy hierarchy are already downloaded but not yet ingested.
