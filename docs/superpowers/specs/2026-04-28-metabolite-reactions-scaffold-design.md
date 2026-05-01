# Phase 1.2 — Reactions scaffold (Metabolite + Reaction nodes; driver query)

**Date:** 2026-04-28
**Status:** ⚠️ **SUPERSEDED** by [`2026-04-30-metabolite-reactions-scaffold-revised.md`](2026-04-30-metabolite-reactions-scaffold-revised.md)
**Parent overview:** [2026-04-28-metabolite-scaffold-design.md](2026-04-28-metabolite-scaffold-design.md)
**Depends on:** [Spec 1.1 — Foundation](2026-04-28-metabolite-foundation-design.md) (must ship first)

> **REVISION NOTICE (2026-04-30):** After the conversations recorded in
> `memory/spec_1_2_design_decisions.md`, the design pivoted on three load-bearing
> points that contradict the body of this spec:
>
> 1. **Primary IDs are KEGG-native**, not MNX-native. KEGG R-numbers and
>    C-numbers are the primary node IDs. MNX is reserved for (a) cross-reference
>    enrichment and (b) Phase 2 paper integration. Mixed-prefix fallback
>    (`kegg.{compound,reaction}` → `chebi:` → `mnx:`) only when a KEGG entry
>    doesn't exist.
> 2. **Reaction direction is not modelled.** No `Reaction_has_substrate` /
>    `Reaction_has_product`. Single direction-agnostic edge
>    `Reaction_has_metabolite`. Substrate/product split was determined to be
>    positional (KEGG editorial convention), not physiological — unreliable.
> 3. **No build-time resolution.** The 2.6 GB MNX SQLite resolver only opens
>    during prepare-data. Build-time adapter reads pre-computed JSON only —
>    same pattern as `build_og_descriptions.py` (step 5).
>
> See the revised spec for the corrected design. This document is preserved
> for history; do not implement against it.

## Goal

Land the chemistry layer of the metabolism scaffold. After this spec ships, the Pro↔Alt nutrient-exchange driver query returns useful candidates from genome content alone.

What this spec adds:

- `Metabolite` and `Reaction` nodes with full property sets.
- Intra-scaffold edges: `Reaction_has_substrate`, `Reaction_has_product`, `Reaction_catalyzed_by → KeggTerm`, `Reaction_classified_by → EcNumber`.
- Per-gene edge: `Gene_catalyzes_reaction` with `evidence_source` property.
- Post-import edge: `Organism_has_metabolite` (substrate/product columns populated; transporter/cazy columns reserved-but-zero, populated by Spec 1.3).
- Rollup properties on Metabolite, Reaction, Gene, EcNumber, KeggTerm.
- New `prepare_data.sh` step 8 (pruned scaffold cache) emitting only metabolite + reaction sections (extended to TC + CAZy in Spec 1.3).
- New adapter `metabolism_adapter.py`.
- Schema, post-import, and KG validity test additions.

What this spec does **not** include:

- TransporterClass and CazyFamily nodes / edges → Spec 1.3.
- Paper measurement integration → Phase 2.
- OrthologGroup inheritance for reactions → Phase 3.

## Step 8 — Pruned scaffold cache

`multiomics_kg/download/build_metabolism_scaffold.py` (new module).

### Inputs

- All strains' `gene_annotations_merged.json` from Spec 1.1.
- `cache/data/mnx/metabolite_resolver.db` from Spec 1.1.
- `cache/data/mnx/{chem,reac}_{prop,xref}.tsv` from Spec 1.1.
- Existing per-KG state (KEGG KOs and EC numbers already declared in current adapters; read by glob over `gene_annotations_merged.json` for any `kegg_ko` or `ec_numbers` field).

### Pruning rule

A node appears in the scaffold cache iff it is **referenced** (directly used by a gene or required as part of a referenced reaction). Same pattern as BRITE pruning.

| Type | Inclusion rule |
|---|---|
| Reaction | Any `MNXR` in any gene's `kegg_reactions` field. |
| Reaction (extended) | Plus reactions whose `kegg_reaction` xref points at a KO already in the KG (via gene `kegg_ko` field). |
| Reaction (extended) | Plus reactions whose `ec_number` xref points at an EC already in the KG (via gene `ec_numbers` field). |
| Metabolite | Any `MNXM` appearing as substrate or product of any retained Reaction. |

The KO-mediated and EC-mediated extensions exist so that `Reaction_catalyzed_by → KeggTerm` and `Reaction_classified_by → EcNumber` edges are useful — they bring in reactions that no gene's `KEGG_Reaction` column listed but which are reachable via existing KO/EC annotations.

Estimated retained counts:

| Type | Count |
|---|---|
| Metabolite | 3,000–8,000 |
| Reaction | 5,000–15,000 |

### Output (`cache/data/mnx/scaffold_pruned.json`)

```json
{
  "schema_version": "1.2",
  "metabolites": {
    "MNXM41": {
      "name": "D-glucose",
      "formula": "C6H12O6",
      "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
      "smiles": "...",
      "charge": 0,
      "mass": 180.063,
      "is_cofactor": "false",
      "kegg_compound": "C00031",
      "chebi": "CHEBI:17234",
      "bigg":  "glc__D",
      "metacyc":   "GLC",
      "modelseed": "cpd00027",
      "hmdb":     "HMDB0000122",
      "aliases":  ["C00031", "CHEBI:17234", "glc__D", "GLC", "cpd00027", "HMDB0000122"],
      "synonyms": ["D-Glucose", "α-D-Glucopyranose", "Dextrose"]
    }
  },
  "reactions": {
    "MNXR101234": {
      "name":             "...",
      "equation_text":    "1 MNXM3 + 1 MNXM41 -> 1 MNXM7 + 1 MNXM58",
      "direction_source": "directional",
      "is_transport_reaction": "false",
      "substrates":       [["MNXM3", 1], ["MNXM41", 1]],
      "products":         [["MNXM7", 1], ["MNXM58", 1]],
      "kegg_terms":       ["K00844", "K12407"],
      "ec_numbers":       ["2.7.1.1", "2.7.1.2"],
      "kegg_reaction":    "R00299",
      "rhea":             "16332",
      "metacyc":          "GLUCOKIN-RXN",
      "bigg":             "HEX1",
      "modelseed":        "rxn00148",
      "aliases":          ["R00299", "16332", "GLUCOKIN-RXN", "HEX1", "rxn00148"]
    }
  }
}
```

`is_cofactor` is set to `"true"` for a hardcoded short list of common cofactors (ATP/ADP/AMP, NAD(P)(H), CoA, water, protons, ortho-/diphosphate, CO2, O2, etc.), keyed by MNXM. Hardcoded list lives in `multiomics_kg/utils/metabolite_id_utils.py:COFACTOR_MNXM_IDS` and is documented inline. Used downstream by Phase 3 query tooling.

`is_transport_reaction` is set to `"true"` when the reaction's substrate set and product set are equal as MNXM ID multisets (a transport reaction shuttles the same compound between compartments). Coarse heuristic; refined later if needed.

`schema_version` is bumped when the cache shape changes; Spec 1.3 bumps it to `"1.3"` after adding `transporter_classes` and `cazy_families` keys.

## Schema changes (`config/schema_config.yaml`)

### Drop dead entries

- The current `compound` node entry (~lines 526–541) and its `compound_targets_protein` edge entry are removed. (They were placeholders for an earlier design and have never been emitted.)

### New nodes

```yaml
metabolite:
  is_a: chemical entity
  represented_as: node
  preferred_id: id
  label_in_input: metabolite
  properties:
    name: str
    formula: str
    inchikey: str
    smiles: str
    charge: int
    mass: float
    is_cofactor: str
    kegg_compound: str
    chebi: str
    bigg: str
    metacyc: str
    modelseed: str
    hmdb: str
    aliases: str[]
    synonyms: str[]

reaction:
  is_a: biochemical reaction
  represented_as: node
  preferred_id: id
  label_in_input: reaction
  properties:
    name: str
    equation_text: str
    direction_source: str
    is_transport_reaction: str
    kegg_reaction: str
    rhea: str
    metacyc: str
    bigg: str
    modelseed: str
    aliases: str[]
```

### New edges

```yaml
gene catalyzes reaction:
  is_a: gene to chemical entity association
  represented_as: edge
  source: gene
  target: reaction
  properties:
    evidence_source: str    # "kegg_reaction" in Phase 1

reaction has substrate:
  is_a: chemical to chemical association
  represented_as: edge
  source: reaction
  target: metabolite
  properties:
    stoichiometry: int

reaction has product:
  is_a: chemical to chemical association
  represented_as: edge
  source: reaction
  target: metabolite
  properties:
    stoichiometry: int

reaction catalyzed by:
  is_a: ontology relationship
  represented_as: edge
  source: reaction
  target: kegg term
  properties: {}

reaction classified by:
  is_a: ontology relationship
  represented_as: edge
  source: reaction
  target: ec number
  properties: {}

organism has metabolite:
  is_a: organism to chemical association
  represented_as: edge
  source: organism taxon
  target: metabolite
  properties:
    as_substrate_count:    int
    as_product_count:      int
    transporter_count:     int     # 0 in Spec 1.2; populated in Spec 1.3
    cazy_count:            int     # 0 in Spec 1.2; populated in Spec 1.3
    evidence_sources:      str[]
    measured_assay_count:  int     # default 0; populated in Phase 2
    measured_compartments: str[]   # default []; populated in Phase 2
    measured_value_kinds:  str[]   # default []; populated in Phase 2
    measured_paper_count:  int     # default 0; populated in Phase 2
```

The Phase-2 properties are declared up front to avoid a schema migration when Phase 2 ships; they default to empty/zero in Phase 1.

## Adapter changes

### New file: `multiomics_kg/adapters/metabolism_adapter.py`

```
MetabolismAdapter(scaffold_path, test_mode):
  download_data() — load scaffold_pruned.json into memory
  get_nodes()     — yield Metabolite then Reaction nodes
  get_edges()     — yield reaction-internal edges:
                      Reaction_has_substrate (per (rxn, sub, stoich))
                      Reaction_has_product   (per (rxn, prod, stoich))
                      Reaction_catalyzed_by  (per (rxn, KO))
                      Reaction_classified_by (per (rxn, EC))

MultiMetabolismAdapter(...):
  Wrapper matching existing patterns; instantiated by create_knowledge_graph.py.
```

The adapter does **no** parsing of MNX/TCDB/CAZy — that all happens in step 8. The adapter only reads the pre-pruned JSON cache.

`Reaction_catalyzed_by → KeggTerm` and `Reaction_classified_by → EcNumber` edges target only KOs and ECs that already exist as nodes in the KG. Step 8 enforces this (`kegg_terms` and `ec_numbers` arrays are filtered to in-KG IDs before writing).

### Extended file: `multiomics_kg/adapters/functional_annotation_adapter.py`

Add a new `GeneCatalyzesReactionAdapter` class + `MultiGeneCatalyzesReactionAdapter` wrapper. For each gene's `kegg_reactions` list, emit one `Gene_catalyzes_reaction` edge with `evidence_source = "kegg_reaction"`. (Property is kept on the edge for future Phase-3 sources; only `"kegg_reaction"` ships in Phase 1.)

Edge ID format: `gene_catalyzes_reaction:{locus_tag}:{mnxr_id}` (deterministic for `omics-edge-snapshot` regression detection).

### `create_knowledge_graph.py`

Instantiate `MultiMetabolismAdapter` and `MultiGeneCatalyzesReactionAdapter` alongside existing adapters. Both call `download_data()` then `write_nodes()` / `write_edges()`.

## Post-import (`scripts/post-import.sh` + `scripts/post-import.cypher`)

### New blocks (in this order, appended to existing post-import flow)

**Block A — `Organism_has_metabolite` edges (substrate/product columns).** Single MERGE per `(organism, metabolite)` pair:

```cypher
CALL { } IN TRANSACTIONS OF 5000 ROWS
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
MATCH (g)-[:Gene_catalyzes_reaction]->(r:Reaction)
MATCH (r)-[:Reaction_has_substrate]->(m:Metabolite)
WITH o, m, count(DISTINCT g) AS as_sub_genes
MERGE (o)-[e:Organism_has_metabolite]->(m)
ON CREATE SET e.as_substrate_count = as_sub_genes,
              e.as_product_count   = 0,
              e.transporter_count  = 0,
              e.cazy_count         = 0,
              e.evidence_sources   = ["kegg_reaction"],
              e.measured_assay_count = 0,
              e.measured_compartments = [],
              e.measured_value_kinds = [],
              e.measured_paper_count = 0
ON MATCH SET  e.as_substrate_count = as_sub_genes;
```

A second pass mirrors the above for products. Reaction direction is honoured **permissively**: when `r.direction_source = "reversible"`, the substrate-side gene count is also added to `as_product_count` and vice versa. Rationale: in vivo direction is condition-dependent; biological reversibility is the conservative default. `direction_source = "directional"` is strict (substrate counts to substrate column only). `direction_source = "unknown"` is treated as `"reversible"` for the rollup (most permissive).

The reversibility logic is encoded as a third pass that adds reversed-direction contributions only where `r.direction_source IN ['reversible', 'unknown']`.

**Block B — Metabolite rollups:**

```cypher
MATCH (m:Metabolite)<-[e:Organism_has_metabolite]-(o:OrganismTaxon)
WITH m,
     collect(DISTINCT CASE WHEN e.as_product_count   > 0 THEN o.preferred_name END) AS producers,
     collect(DISTINCT CASE WHEN e.as_substrate_count > 0 THEN o.preferred_name END) AS consumers,
     sum(e.as_substrate_count + e.as_product_count) AS total_links
WITH m, [x IN producers WHERE x IS NOT NULL] AS prod, [x IN consumers WHERE x IS NOT NULL] AS cons, total_links
SET m.producer_organism_count = size(prod),
    m.consumer_organism_count = size(cons),
    m.producer_genera = [x IN apoc.coll.toSet([n IN prod | head(split(n, ' '))]) ],
    m.consumer_genera = [x IN apoc.coll.toSet([n IN cons | head(split(n, ' '))]) ],
    m.predicted_pro_alt_exchange_flag =
        CASE WHEN ('Prochlorococcus' IN m.producer_genera AND 'Alteromonas' IN m.consumer_genera)
              OR ('Alteromonas'      IN m.producer_genera AND 'Prochlorococcus' IN m.consumer_genera)
             THEN 'true' ELSE 'false' END,
    m.total_gene_links = total_links;
```

`producer_genera` and `consumer_genera` derive from the first whitespace-separated token of `preferred_name` (matches existing convention with `genera` on OrthologGroup).

**Block C — Reaction rollups:**

```cypher
MATCH (g:Gene)-[:Gene_catalyzes_reaction]->(r:Reaction)
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WITH r, count(DISTINCT g) AS gc, collect(DISTINCT o.preferred_name) AS orgs
SET r.gene_link_count = gc,
    r.organism_count  = size(orgs),
    r.genera          = [x IN apoc.coll.toSet([n IN orgs | head(split(n, ' '))]) ];
```

**Block D — Gene routing-signal updates:**

```cypher
MATCH (g:Gene)
OPTIONAL MATCH (g)-[r:Gene_catalyzes_reaction]->()
WITH g, count(r) AS rc
SET g.reaction_count = rc;
```

(`g.transporter_count` and `g.cazy_count` are added in Spec 1.3.)

**Block E — Ontology rollups:**

```cypher
MATCH (k:KeggTerm)<-[:Reaction_catalyzed_by]-(r:Reaction)
WITH k, count(DISTINCT r) AS rc SET k.reaction_count = rc;

MATCH (e:EcNumber)<-[:Reaction_classified_by]-(r:Reaction)
WITH e, count(DISTINCT r) AS rc SET e.reaction_count = rc;
```

### New indexes

Scalar:
- `metabolite_inchikey_idx` on `Metabolite(inchikey)`
- `metabolite_name_idx`     on `Metabolite(name)`
- `reaction_id_idx`         on `Reaction(id)`
- `reaction_kegg_idx`       on `Reaction(kegg_reaction)`

Full-text:
- `metaboliteFullText` on `Metabolite(name, aliases, synonyms)`
- `reactionFullText`   on `Reaction(equation_text, name)`

### Reference Cypher parity

`scripts/post-import.cypher` is kept byte-identical to the bash version per existing convention. Validation: run `scripts/post-import-validate.sh > before.txt` against the currently deployed graph, then `diff before.txt after.txt` after rebuild — must be byte-identical for any non-functional refactor.

## KG validity tests

New file: `tests/kg_validity/test_metabolism_reactions.py` (marker `@pytest.mark.kg`).

- `test_metabolite_count_minimum` — ≥ 2,000 nodes.
- `test_reaction_count_minimum` — ≥ 3,000 nodes.
- `test_no_orphan_reactions` — every Reaction has ≥1 substrate AND ≥1 product.
- `test_no_orphan_metabolites` — every Metabolite has at least one Reaction edge (substrate or product).
- `test_gene_catalyzes_reaction_targets_exist` — every `Gene_catalyzes_reaction` targets an existing Reaction.
- `test_reaction_substrate_targets_exist` — every `Reaction_has_substrate` and `Reaction_has_product` targets an existing Metabolite.
- `test_reaction_kegg_term_targets_exist` — every `Reaction_catalyzed_by` targets an existing KeggTerm in the KG.
- `test_reaction_ec_targets_exist` — every `Reaction_classified_by` targets an existing EcNumber in the KG.
- `test_organism_has_metabolite_present` — both Prochlorococcus and Alteromonas have ≥100 `Organism_has_metabolite` edges each.
- `test_predicted_exchange_flag_threshold` — `count(Metabolite WHERE predicted_pro_alt_exchange_flag = "true") ≥ 50`.
- `test_known_exchange_compounds_present` — for each of glycine betaine, DMSP, D-glucose, L-alanine, glycolate, ammonium, the corresponding Metabolite node exists with `predicted_pro_alt_exchange_flag = "true"`.
- `test_metabolite_rollup_consistency` — `producer_organism_count = size(producer_genera_origins)` etc. (sanity).
- `test_reaction_rollup_consistency` — `gene_link_count > 0 ⇔ organism_count > 0`.

Update `test_snapshot.py` fixture: regenerate after first successful build via `uv run python tests/kg_validity/generate_snapshot.py` so subsequent rebuilds are checked against the new baseline.

## Driver-goal queries enabled by this spec

**Pro↔Alt exchange candidates from genome content alone:**

```cypher
MATCH (pro:OrganismTaxon)-[r1:Organism_has_metabolite]->(m:Metabolite)
      <-[r2:Organism_has_metabolite]-(alt:OrganismTaxon)
WHERE pro.preferred_name STARTS WITH 'Prochlorococcus'
  AND alt.preferred_name STARTS WITH 'Alteromonas'
  AND r1.as_product_count   > 0   // Pro biosynthesis
  AND r2.as_substrate_count > 0   // Alt consumption
RETURN m.name, m.kegg_compound,
       r1.as_product_count   AS pro_producers,
       r2.as_substrate_count AS alt_consumers
ORDER BY pro_producers + alt_consumers DESC
LIMIT 50;
```

**Pro genes that could produce a specific compound:**

```cypher
MATCH (m:Metabolite {name: "D-glucose"})<-[:Reaction_has_product]-(r:Reaction)
      <-[:Gene_catalyzes_reaction]-(g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.preferred_name STARTS WITH 'Prochlorococcus'
RETURN o.preferred_name, collect(DISTINCT g.locus_tag), collect(DISTINCT r.name);
```

## Estimated KG impact

| Quantity | Estimate |
|---|---|
| Metabolite nodes | 3,000–8,000 |
| Reaction nodes | 5,000–15,000 |
| `Gene_catalyzes_reaction` edges | 30K–80K |
| `Reaction_has_substrate` edges | 12K–40K |
| `Reaction_has_product` edges | 12K–40K |
| `Reaction_catalyzed_by → KeggTerm` edges | ~10K |
| `Reaction_classified_by → EcNumber` edges | ~10K |
| `Organism_has_metabolite` edges (post-import) | 50K–100K |
| **Total new edges** | **~120K–280K** |

## Acceptance criteria

- All Phase 1.2 nodes and edges present in the deployed graph at expected magnitudes (within 50% of the estimates table).
- New KG validity tests pass.
- The Pro↔Alt driver query (above) returns ≥1 result.
- For at least 5 of {glycine betaine, DMSP, D-glucose, L-alanine, glycolate, ammonium}, the corresponding Metabolite node exists with `predicted_pro_alt_exchange_flag = "true"`.
- `omics-edge-snapshot` confirms no regression in existing expression-edge counts.
- `scripts/post-import.cypher` and `scripts/post-import.sh` produce byte-identical computed properties when invoked against the same graph state.
- CLAUDE.md updated with the new step 8, new node/edge types, new indexes, and the post-import additions.
- The deployed graph still exits cleanly (no `import.report` rejections) — verified by `test_import_report.py`.

## Out of scope (handled in later sub-specs)

- TransporterClass and CazyFamily nodes / edges → Spec 1.3.
- Backfill of `Organism_has_metabolite.transporter_count` and `cazy_count` → Spec 1.3.
- MetaboliteAssay (paper measurements) → Phase 2.
- OrthologGroup inheritance for reactions → Phase 3.
