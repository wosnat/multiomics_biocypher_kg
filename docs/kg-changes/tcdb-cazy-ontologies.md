# TCDB and CAZy Ontologies (TcdbFamily, CazyFamily nodes)

**Status:** LANDED 2026-05-02 — design at [docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md](../superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md). Counts below are measured against the post-rebuild graph.

**Updated 2026-05-03:** Substrate edges (`Tcdb_family_transports_metabolite`) are rolled up from `tc_specificity` leaves to every ancestor, flattening the descendants walk into edges so all explorer queries become single-hop at any TcdbFamily level. Recover leaf-only semantics via `WHERE source.level_kind = 'tc_specificity'`. Same date: chemistry-slice-1 follow-up (KG-A5..A8) added pathway/organism rollup arrays onto Metabolite — see "Metabolite — chemistry-slice-1 follow-up rollups" below.

**Updated 2026-05-08:** Substrate rollup moved from adapter time to **step-6 time, computed pre-prune over the FULL TCDB hierarchy**. Previously the rollup walked only the pruned hierarchy, so kept ancestors reached only by `walk_up` from a deeper gene seed undercounted: substrates from non-gene-annotated sibling specificities under the same ancestor were missed (~950 substrate edges across ~50 ancestor families). Step 6 now (1) computes the subtree rollup over the full hierarchy as raw `CHEBI:NNNN;name` strings, (2) prunes structurally, (3) resolves only the surviving distinct substrate strings via MNX. Result lives in `tcdb_pruned.json` as `subtree_substrates: {tc_id: [primary_id, ...]}`; the adapter emits one edge per pair without traversing the hierarchy. Total edges: 21,530 → **22,476**. Affected ancestor `metabolite_count`/`is_promiscuous` accordingly.

**Proposed:** 2026-05-01

## What's changing

Two new ontologies join GO, EC, Pfam, BRITE, KEGG, COG, CyanorakRole, TigrRole as first-class hierarchical classifications in the graph:

- **TcdbFamily** — Transporter Classification Database. ~13K curated transport-protein families across 5 levels of specificity. Each leaf carries a curated list of substrates the family is known to move.
- **CazyFamily** — Carbohydrate-Active enZymes. 6 classes (GH, GT, PL, CE, AA, CBM) of glycoside hydrolases, glycosyltransferases, etc., with optional subfamilies.

Both surface through the existing MCP ontology tools (`genes_by_ontology`, `gene_ontology_terms`, `ontology_landscape`, `search_ontology`) without any MCP code changes.

The biggest change for explorers: **TCDB substrates become Metabolite nodes**. **837 new transport-evidence Metabolites** join the chemistry layer (2,188 → 3,025). Every Metabolite gains an `evidence_sources: str[]` property indicating whether it was reached via metabolism (existing reactions) or transport (new TCDB substrate annotations) or both.

A second consequence (post-merge refactor): **transport-reachable KEGG pathways** also enter the graph. When a TCDB substrate resolves to a `kegg.compound:*` Metabolite, that compound's KEGG pathway memberships pull the corresponding pathways into the gene-reachable pws set — even if no local gene catalyzes a reaction in those pathways. The pathway count grows from 377 (catalysis-only) to 442 (catalysis ∪ transport-reachable). 65 new pathways enter, mostly transporter-rich (e.g. ABC transporters ko02010, two-component systems ko02020). 15 nameless KEGG meta-classification maps (ko010** plant biosynthesis overviews and ko07*** drug-classification maps with no KOs/reactions) are filtered out — they have no biological meaning for marine bacteria.

The flat Gene properties `transporter_classification` and `cazy_ids` are removed — query through the new graph edges instead. Same precedent as the Pfam normalization.

## New node label: TcdbFamily

Node IDs follow the bioregistry pattern: `tcdb:1.A.1.5.2`

Properties:

- `name` (str) — TCDB entry name, e.g. "The Voltage-gated Ion Channel (VIC) Superfamily". Falls back to `tcdb_id` when the source name is empty (typical for `tc_subclass`/`tc_subfamily`/`tc_specificity`).
- `tcdb_id` (str) — bare hierarchical ID, e.g. `"1.A.1.5.2"`.
- `level` (int) — hierarchy depth, 0 = broadest.
- `level_kind` (str) — one of `tc_class` (0), `tc_subclass` (1), `tc_family` (2), `tc_subfamily` (3), `tc_specificity` (4).
- `superfamily` (str, sparse) — leaf-only superfamily annotation, e.g. `"VIC Superfamily"`.
- `tc_class_id` (str, post-import) — pointer to root `tc_class` node ID (e.g. `"tcdb:3"` for any node under class 3); self on `tc_class`-level nodes. Lets class-level filters skip variable-length traversals: `WHERE tf.tc_class_id = 'tcdb:3'`.
- `gene_count` (int, post-import) — distinct genes reachable via subtree.
- `organism_count` (int, post-import) — distinct organisms among those genes.
- `member_count` (int, post-import) — direct child count in the hierarchy.
- `metabolite_count` (int, post-import) — distinct Metabolite nodes reachable via `Tcdb_family_transports_metabolite` (single-hop after the 2026-05-03 rollup; substrate edges exist directly on every ancestor). On a `tc_class` node, shows total substrate breadth across all descendants.

## New node label: CazyFamily

Node IDs: `cazy:GH13`, `cazy:GH13_5`, etc.

Properties:

- `name` (str) — human-readable for class (e.g. "Glycoside Hydrolases"). Falls back to `cazy_id` for family/subfamily levels.
- `cazy_id` (str) — e.g. `"GH13"` or `"GH13_5"`.
- `level` (int) — 0 = class, 1 = family, 2 = subfamily.
- `level_kind` (str) — `cazy_class` | `cazy_family` | `cazy_subfamily`.
- `gene_count` (int, post-import) — distinct genes reachable via subtree.
- `organism_count` (int, post-import) — distinct organisms.

## New edge types

| Edge type | Source | Target | Count | Notes |
|---|---|---|---|---|
| `Gene_has_tcdb_family` | Gene | TcdbFamily | **10,576** | Attaches at exact level annotated by eggNOG (no walk-up at edge time) |
| `Tcdb_family_is_a_tcdb_family` | TcdbFamily (child) | TcdbFamily (parent) | **4,838** | Hierarchy parent edge (one per non-root TcdbFamily) |
| `Tcdb_family_transports_metabolite` | TcdbFamily | Metabolite | **13,641** | Substrate linkage rolled up from leaves to every ancestor (tc_class: 1,497, tc_subclass: 1,510, tc_family: 2,040, tc_subfamily: 2,832, tc_specificity: 5,762). Filter `WHERE source.level_kind = 'tc_specificity'` for leaf-only (1,097 distinct primary IDs across 3,095 leaves). |
| `Gene_has_cazy_family` | Gene | CazyFamily | **1,181** | |
| `Cazy_family_is_a_cazy_family` | CazyFamily (child) | CazyFamily (parent) | **58** | |

## TCDB pruning

The full TCDB hierarchy has ~13,642 entries; 99% of those are not relevant to our 25 genome strains. The graph keeps only the **subhierarchy reachable from gene annotations**:

- Seed set: every TCDB ID present in any strain's `transporter_classification` from eggNOG.
- For each seed, walk **down** to all leaf descendants (so substrate linkage is complete) and **up** to `tc_class` (so the hierarchy reaches root).
- Kept = union of seeds + descendants + ancestors.

This is the same approach as BRITE pruning (one-way), made bidirectional because TCDB substrates live at the leaves but eggNOG often annotates at family/subfamily level.

## CAZy is observed-only

CAZy hierarchy is inferred at adapter init from observed eggNOG annotations across configured strains. Every CazyFamily node corresponds to a real gene annotation. No external download — the 6 class names (GH, GT, PL, CE, AA, CBM) are hardcoded in the adapter.

## TCDB substrate → Metabolite bridge

This is the most consequential semantic change. Every TCDB family with substrate annotations gets `Tcdb_family_transports_metabolite` edges to actual Metabolite nodes — not properties, not strings.

How substrates are resolved:

1. TCDB curates substrates as `CHEBI:NNNN;name` strings (e.g. `CHEBI:9314;sucrose`, `CHEBI:3308;calcium(2+)`).
2. Each ChEBI ID is resolved through MetaNetX (MNX) to its canonical compound.
3. The resulting compound is matched against existing Metabolite nodes (by `kegg.compound:C*` ID, then `chebi:NNNN`, then `mnx:MNXM*`).
4. If the substrate's primary ID is `kegg.compound:*` (645 of 1,097 distinct primaries), the compound is added to the gene-reachable cpds set so it rides through the regular MNX-enrichment pipeline and picks up its KEGG pathway memberships natively. Reuses the existing Metabolite node if it was already catalysis-reachable; otherwise becomes a transport-only KEGG Metabolite.
5. If the substrate's primary ID is non-KEGG (`chebi:*` or `mnx:MNXM*` — the remaining 452), it lands in `kegg_data.json`'s `additional_compounds` bucket as a new transport-only Metabolite (covers compounds like tetracycline that TCDB knows but KEGG doesn't curate).

End result: **3,025 Metabolites** = 1,928 metabolism-only + 260 metabolism+transport overlap + 385 transport-only KEGG + 452 transport-only non-KEGG (additional_compounds).

## Metabolite — new and extended properties

### New: `evidence_sources: str[]`

Every Metabolite (including pre-existing ones) now carries an `evidence_sources` array. Values:

- `"metabolism"` — at least one Reaction in the KG involves this compound.
- `"transport"` — at least one TcdbFamily is annotated as transporting this compound.

A compound found by both paths gets `["metabolism", "transport"]`. A pre-existing pure-metabolism compound stays `["metabolism"]`. The 837 new transport-evidence Metabolites carry just `["transport"]`.

| Bucket | Count | `evidence_sources` |
|---|---|---|
| Metabolism-only (catalysis-reachable, not transport) | 1,928 | `['metabolism']` |
| Both | 260 | `['metabolism', 'transport']` |
| Transport-only KEGG (substrate `kegg.compound:*`, no local catalysis) | 385 | `['transport']` |
| Transport-only non-KEGG (substrate `chebi:*` or `mnx:MNXM*`) | 452 | `['transport']` |

This lets you distinguish "this organism has a known metabolic reaction with this compound" from "this organism's transporters are annotated as moving this compound" — both are evidence, but they answer different biological questions.

The enum is **open-ended**: a future metabolomics-DM spec is expected to add `"metabolomics"` for compounds measured in metabolomics experiments. Filter logic should use set membership (`'transport' IN m.evidence_sources`) rather than equality so it stays robust to additions.

### New: `transporter_count: int`

Distinct **`tc_specificity` leaf** `TcdbFamily` nodes pointing at this Metabolite via `Tcdb_family_transports_metabolite`. Filtered to leaves so the count reflects "actual transporter systems" rather than ancestors-via-rollup (post-import filters `level_kind = 'tc_specificity'`). Useful as a quick "how many transporter families curate this compound as a substrate" signal.

### Extended semantics: `gene_count` / `organism_count`

Both already exist on Metabolite, populated by post-import as a 2-hop count via `Gene → Reaction → Metabolite`. After this spec lands, the materialization adds a UNION arm via `Gene → TcdbFamily → Metabolite` (single-hop after the 2026-05-03 rollup). Same property, broader semantics:

> Distinct genes / organisms reachable via *any* chemistry path — catalysis OR transport.

A Metabolite with `evidence_sources: ['transport']` (transport-only) will show non-zero `gene_count` and `organism_count` from the transport path alone. Pre-existing pure-metabolism Metabolites grow their counts only if any of their consumer organisms also has a transporter for them.

**2026-05-03 fix.** Pre-rollup, the transport arm was implemented at post-import time as `(g)-[:Gene_has_tcdb_family]->(:TcdbFamily)<-[:Gene_has_tcdb_family]-(...)` matching only when the gene was annotated **at the `tc_specificity` leaf**. Genes annotated at higher levels (tc_family, tc_subfamily) were missed, so `gene_count` and `organism_count` undercounted relative to the materialized `Organism_has_metabolite` edge. After the rollup, both properties derive from single-hop traversal of the rolled-up substrate edges and are now consistent (verified: `size(m.organism_names) == m.organism_count` invariant holds for all 3,025 Metabolites; was failing for 952 pre-fix).

### Extended materialization: `Organism_has_metabolite` edges

Already exist as 2-hop saves via `Gene → Reaction → Metabolite`. Post-import adds a UNION arm via `Gene → TcdbFamily → Metabolite` (single-hop after the 2026-05-03 rollup; previously walked descendants via `*0..` to a `tc_specificity` leaf). Single dedup by `(organism, metabolite)`. Optional `via: str[]` property on the edge enumerating which paths produced it (`{'metabolism', 'transport'}`).

## Metabolite — chemistry-slice-1 follow-up rollups (2026-05-03, KG-A5..A8)

Four denormalized arrays / scalars added to every Metabolite to collapse `EXISTS`-subquery filter clauses and per-row edge traversals in `list_metabolites`:

| Property | Type | Source |
|---|---|---|
| `pathway_ids` | `list[str]` | Distinct sorted `KeggTerm.id` reachable via `Metabolite_in_pathway`. Empty list when no memberships. |
| `pathway_names` | `list[str]` | `KeggTerm.name`, **index-aligned** with `pathway_ids` (sorted by `KeggTerm.id`). |
| `pathway_count` | `int` | `size(pathway_ids)`. |
| `organism_names` | `list[str]` | Distinct sorted `OrganismTaxon.preferred_name` reachable via `Organism_has_metabolite` (UNION of catalysis + transport). Invariant: `size(organism_names) == organism_count`. |

**Filter pattern (before):**

```cypher
WHERE EXISTS {
  MATCH (m)-[:Metabolite_in_pathway]->(p:KeggTerm)
  WHERE p.id IN $pathway_ids
}
```

**Filter pattern (after):**

```cypher
WHERE ANY(p IN $pathway_ids WHERE p IN m.pathway_ids)
```

Same shape for `organism_names` against `Organism_has_metabolite`. Roughly 2-3× faster on filtered detail queries spanning thousands of metabolites.

## New Gene routing-signal properties

To match the existing rollup-based routing pattern (`expression_edge_count`, `numeric_metric_count`, etc.) the post-import script populates two new counts on every Gene:

- `tcdb_family_count` (int, default 0) — number of `Gene_has_tcdb_family` edges (i.e. distinct TCDB family memberships across all levels for this gene).
- `cazy_family_count` (int, default 0) — number of `Gene_has_cazy_family` edges.

The pre-existing `annotation_types: str[]` array also gains `'tcdb'` / `'cazy'` values when the corresponding count is > 0.

`Gene.metabolite_count` (the property introduced by the chemistry-slice-1 chemistry-layer asks) is defined as **distinct Metabolite nodes reachable via *any* gene-reaching path** — UNION of catalysis (`Gene → Reaction → Metabolite`) and transport (`Gene → TcdbFamily → Metabolite`). On TCDB landing the count grows to include transport substrates automatically; consumers read one property regardless of source path. **2026-05-03**: post-import computes the transport arm via single-hop edges (no descendants walk) thanks to the substrate-edge rollup.

## Removed Gene properties (breaking change)

Gone:

- `transporter_classification: str[]`
- `cazy_ids: str[]`

Migration: existing queries using these properties must rewrite to traverse the new graph edges. Examples:

```cypher
// Before:
MATCH (g:Gene)
WHERE 'GH13' IN g.cazy_ids
RETURN g

// After:
MATCH (g:Gene)-[:Gene_has_cazy_family]->(cf:CazyFamily {cazy_id: 'GH13'})
RETURN g

// Before:
MATCH (g:Gene)
WHERE any(tc IN g.transporter_classification WHERE tc STARTS WITH '3.A.1')
RETURN g

// After (single-level):
MATCH (g:Gene)-[:Gene_has_tcdb_family]->(tf:TcdbFamily {tcdb_id: '3.A.1'})
RETURN g

// After (any descendant of 3.A.1):
MATCH (g:Gene)-[:Gene_has_tcdb_family]->(tf:TcdbFamily)
MATCH (tf)-[:Tcdb_family_is_a_tcdb_family*0..]->(:TcdbFamily {tcdb_id: '3.A.1'})
RETURN g
```

`gene_summary` and `geneFullText` continue to surface raw IDs for text-search needs.

## What this enables

### Substrate-driven gene queries

```cypher
// All genes in MED4 whose transporters move calcium.
// Single-hop after the 2026-05-03 substrate-edge rollup: substrate edges exist
// on every TcdbFamily ancestor, so the gene's annotation level no longer
// matters and no descendants walk is needed.
MATCH (org:OrganismTaxon {strain_name: 'MED4'})
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(org)
MATCH (g)-[:Gene_has_tcdb_family]->(:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
WHERE m.name CONTAINS 'calcium'
RETURN DISTINCT g.locus_tag, g.product, m.name
```

### Hierarchy traversal

```cypher
// What does TCDB class "1: Channels and Pores" cover in MED4?
// (Using the post-import-computed tc_class_id pointer to skip the variable-length traversal)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(:OrganismTaxon {strain_name: 'MED4'})
MATCH (g)-[:Gene_has_tcdb_family]->(tf:TcdbFamily {tc_class_id: 'tcdb:1'})
RETURN tf.tcdb_id, tf.name, tf.level_kind, count(g) AS gene_count
ORDER BY tf.level, tf.tcdb_id
```

### Transport vs metabolism distinction

```cypher
// Compounds the organism is annotated to TRANSPORT but has no metabolism for
MATCH (m:Metabolite)
WHERE 'transport' IN m.evidence_sources
  AND NOT 'metabolism' IN m.evidence_sources
RETURN m.name, m.formula, m.evidence_sources
LIMIT 50
```

### Ontology landscape via MCP

After the work lands, `ontology_landscape` will show TCDB and CAZy alongside GO, EC, Pfam, BRITE, KEGG, COG, CyanorakRole, TigrRole. `search_ontology` will return TCDB families and CAZy families when searching by chemical name or family code.

### Genes by carbohydrate activity

```cypher
// All glycosyltransferase genes across organisms
MATCH (g:Gene)-[:Gene_has_cazy_family]->(:CazyFamily)
      -[:Cazy_family_is_a_cazy_family*0..]->(:CazyFamily {cazy_id: 'GT'})
RETURN g.organism_name, g.locus_tag, g.product
ORDER BY g.organism_name
```

### Routing-signal-driven Gene shortlists

```cypher
// Genes with rich TCDB annotation (≥2 family memberships across the hierarchy)
MATCH (g:Gene {organism_name: 'Prochlorococcus MED4'})
WHERE g.tcdb_family_count >= 2
RETURN g.locus_tag, g.product, g.tcdb_family_count
ORDER BY g.tcdb_family_count DESC
LIMIT 20
```

### Most chemically diverse TCDB families

```cypher
// Top transporter families by substrate breadth
MATCH (tf:TcdbFamily {level_kind: 'tc_family'})
WHERE tf.metabolite_count > 0
RETURN tf.tcdb_id, tf.name, tf.metabolite_count, tf.gene_count
ORDER BY tf.metabolite_count DESC
LIMIT 20
```

## Pipeline changes (for awareness only)

If you re-run prepare_data:

- TCDB reference TSV download moves from `prepare_data.sh` step 0 sub-step 6 → step 6 (alongside KEGG metabolism caches). One less script to think about; everything chemistry-related lives in step 6.
- `cache/data/cazy/cazy_hierarchy.json` is deleted — the CAZy hierarchy is now a hardcoded Python table; no file artifact.
- `build_metabolite_resolver.py` is renamed to `build_mnx_resolver.py` (does only one thing now).
- `cache/data/kegg/kegg_data.json` gains an `additional_compounds` section for non-KEGG compounds (transport-only chemistry).

## What does NOT change

- `Changes_expression_of` edge counts and properties — omics layer untouched.
- Existing Reaction nodes, Gene_catalyzes_reaction edges, Reaction_has_metabolite edges — unchanged.
- Direction of transport — explicitly NOT modeled (consistent with spec 1.2 "no reaction direction"). `Tcdb_family_transports_metabolite` is direction-agnostic, same as `Reaction_has_metabolite`.
- TransportReaction nodes — explicitly NOT introduced. TCDB is treated as an ontology (like EC), not a reaction layer (rationale in the design spec).
- Existing Metabolite IDs, names, formulas, etc. — unchanged. The new `evidence_sources` property is additive.

## Graph-size impact (measured)

| Quantity | Before | After | Delta |
|---|---|---|---|
| Gene nodes | 81,458 | 81,458 | 0 |
| Metabolite nodes | 2,188 | 3,025 | **+837** |
| TcdbFamily nodes | 0 | 4,844 | **new** (pruned from 13,643 in raw TCDB) |
| CazyFamily nodes | 0 | 64 | **new** |
| KeggTerm pathway-level nodes | 377 | 442 | **+65** (transport-reachable extension; 15 nameless meta-classification maps filtered out) |
| `Changes_expression_of` edges | 232,439 | 232,439 | 0 |
| `Gene_has_tcdb_family` edges | 0 | 10,576 | **new** |
| `Gene_has_cazy_family` edges | 0 | 1,181 | **new** |
| `Tcdb_family_transports_metabolite` edges | 0 | 13,641 | **new** (rolled up from 5,762 leaf-only edges to all 5 levels via 2026-05-03 adapter change) |
| `Tcdb_family_is_a_tcdb_family` edges | 0 | 4,838 | **new** |
| `Cazy_family_is_a_cazy_family` edges | 0 | 58 | **new** |
| `Organism_has_metabolite` edges | 37,010 | 56,898 | **+19,888** (transport arm contributing) |
| `Metabolite_in_pathway` edges | 8,095 | 9,444 | **+1,349** (transport-extended pathway memberships) |
| `Reaction_has_metabolite` edges | 10,050 | 10,050 | 0 |
| `Reaction_in_kegg_pathway` edges | 6,349 | 6,349 | 0 |
| `Gene_catalyzes_reaction` edges | 52,742 | 52,742 | 0 |

## Property changes summary

At-a-glance reference for what to update on the explorer side:

| Node type | Property | Status | Notes |
|---|---|---|---|
| Gene | `transporter_classification: str[]` | **REMOVED** | Use `(g)-[:Gene_has_tcdb_family]->(tf)` instead |
| Gene | `cazy_ids: str[]` | **REMOVED** | Use `(g)-[:Gene_has_cazy_family]->(cf)` instead |
| Gene | `tcdb_family_count: int` | NEW | Routing signal (per TCDB-S1) |
| Gene | `cazy_family_count: int` | NEW | Routing signal (per TCDB-S2) |
| Gene | `metabolite_count: int` | EXTENDED | UNION across catalysis + transport paths (per TCDB-S3); the underlying property is introduced by chemistry-slice-1 KG-A2 |
| Gene | `annotation_types: str[]` | EXTENDED | Now also contains `'tcdb'` / `'cazy'` |
| Metabolite | `evidence_sources: str[]` | NEW | `metabolism` / `transport`; `metabolomics` reserved |
| Metabolite | `transporter_count: int` | NEW | Distinct **`tc_specificity` leaf** TCDB families pointing at this metabolite (filtered to leaves post 2026-05-03 rollup so the count keeps "actual transporter systems" meaning) |
| Metabolite | `gene_count: int` | EXTENDED | Now UNION'd with transport path; **2026-05-03**: derives from rolled-up substrate edges, fixes pre-existing transport-arm undercount |
| Metabolite | `organism_count: int` | EXTENDED | Now UNION'd with transport path; **2026-05-03**: derives from materialized `Organism_has_metabolite` edge so `size(organism_names) == organism_count` is invariant |
| Metabolite | `pathway_ids: list[str]` | NEW (KG-A5, 2026-05-03) | Distinct sorted `KeggTerm.id` reachable via `Metabolite_in_pathway`; empty when none |
| Metabolite | `pathway_names: list[str]` | NEW (KG-A6, 2026-05-03) | Index-aligned with `pathway_ids` |
| Metabolite | `pathway_count: int` | NEW (KG-A7, 2026-05-03) | `size(pathway_ids)` |
| Metabolite | `organism_names: list[str]` | NEW (KG-A8, 2026-05-03) | Distinct sorted `OrganismTaxon.preferred_name` reachable via `Organism_has_metabolite` |
| TcdbFamily | (entire node type) | NEW | See properties section above |
| CazyFamily | (entire node type) | NEW | See properties section above |

## See also

- [Design spec](../superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md) — implementation details, commit-by-commit migration plan, cross-spec coordination notes.
- [Chemistry-slice-1 KG-side asks](../../../multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md) — adjacent chemistry-layer rollups (`Gene.reaction_count`, `Gene.metabolite_count`, `Metabolite.elements`, `KeggTerm.reaction_count` / `metabolite_count`) owned by the explorer team's chemistry-slice-1 PR. Coordinated via TCDB-S3 (`Gene.metabolite_count` UNION semantics).
- [BRITE Categories](brite-categories.md) — closest precedent for ontology-with-pruning.
- [Ontology Level](ontology-level.md) — unified `level` property convention shared by all ontologies.
- [Metabolism Chemistry Layer](metabolism-chemistry-layer.md) — the existing chemistry layer that TCDB substrates fold into.
