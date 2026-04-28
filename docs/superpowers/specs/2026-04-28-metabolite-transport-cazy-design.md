# Phase 1.3 — Transport + CAZy scaffold

**Date:** 2026-04-28
**Status:** Spec ready for implementation plan
**Parent overview:** [2026-04-28-metabolite-scaffold-design.md](2026-04-28-metabolite-scaffold-design.md)
**Depends on:** [Spec 1.1 — Foundation](2026-04-28-metabolite-foundation-design.md), [Spec 1.2 — Reactions scaffold](2026-04-28-metabolite-reactions-scaffold-design.md) (both must ship first)

## Goal

Add the protein-family scaffold orthogonal to reactions: TransporterClass (TCDB hierarchy) and CazyFamily (CAZy hierarchy) nodes + per-gene edges. Backfill `Organism_has_metabolite` with TC/CAZy counts so the existing driver query benefits without changes.

What this spec adds:

- `TransporterClass` and `CazyFamily` nodes with full property sets and 5-level / 3-level hierarchies.
- Per-gene edges: `Gene_has_transporter_classification`, `Gene_in_cazy_family`.
- Hierarchy edges: `Transporter_class_is_a → Transporter_class`, `Cazy_family_is_a → Cazy_family`.
- Step 8 extension: scaffold cache now also emits `transporter_classes` and `cazy_families` sections.
- Post-import backfill of `Organism_has_metabolite.transporter_count` and `cazy_count`.
- New routing-signal properties on Gene: `transporter_count`, `cazy_count`.
- Schema, post-import, and KG validity test additions.

What this spec does **not** include:

- Direct `TransporterClass → Metabolite` substrate-compound edges (kept as text property `substrate_classes`; deferred to Phase 3).
- Paper measurement integration → Phase 2.

## Step 8 — Scaffold cache extension

The step-8 module from Spec 1.2 (`build_metabolism_scaffold.py`) is extended to emit two additional top-level sections in `cache/data/mnx/scaffold_pruned.json`. Schema version bumped to `"1.3"`.

### Pruning rule

| Type | Inclusion rule |
|---|---|
| TransporterClass | Any TC value in any gene's `transporter_classification` field, plus all ancestor classes up to the root. |
| CazyFamily | Any CAZy value in any gene's `cazy_ids` field, plus all ancestor classes up to the root. |

Ancestors are derived from `tcdb_hierarchy.json` and `cazy_hierarchy.json` (built in Spec 1.1). A leaf at level 4 (e.g. `1.A.1.1.1`) drags in its ancestors `1`, `1.A`, `1.A.1`, `1.A.1.1`.

Estimated retained counts: ~1,500 TC nodes (full TCDB ~12K), ~400 CAZy nodes (full CAZy ~700).

### Output additions

```json
{
  "schema_version": "1.3",
  "metabolites": {...},                 // unchanged
  "reactions":   {...},                 // unchanged
  "transporter_classes": {
    "1":         {"name": "Channels and Pores", "level": 0, "level_kind": "tc_class",       "parent": null,        "substrate_classes": []},
    "1.A":       {"name": "α-Type Channels",    "level": 1, "level_kind": "tc_subclass",    "parent": "1",         "substrate_classes": []},
    "1.A.1":     {"name": "VIC Family",         "level": 2, "level_kind": "tc_family",      "parent": "1.A",       "substrate_classes": []},
    "1.A.1.1":   {"name": "Kv subfamily",       "level": 3, "level_kind": "tc_subfamily",   "parent": "1.A.1",     "substrate_classes": []},
    "1.A.1.1.1": {"name": "Shaker channel",     "level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.1",   "substrate_classes": ["potassium ion"]}
  },
  "cazy_families": {
    "GH":      {"name": "Glycoside Hydrolases", "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "GH"},
    "GH13":    {"name": "GH13 family",          "level": 1, "level_kind": "cazy_family",    "parent": "GH",    "class": "GH"},
    "GH13_1":  {"name": "GH13 subfamily 1",     "level": 2, "level_kind": "cazy_subfamily", "parent": "GH13",  "class": "GH"}
  }
}
```

`substrate_classes` for TC entries is propagated only at the leaf level (the TCDB level where it's declared); interior nodes get `[]`. Phase 3 may compute aggregated substrate-class arrays at interior levels, but Phase 1.3 keeps the leaf-only convention.

## Schema changes (`config/schema_config.yaml`)

### New nodes

```yaml
transporter class:
  is_a: ontology term
  represented_as: node
  preferred_id: tcdb
  label_in_input: transporter_class
  properties:
    name: str
    description: str
    level: int
    level_kind: str         # tc_class | tc_subclass | tc_family | tc_subfamily | tc_specificity
    substrate_classes: str[]

cazy family:
  is_a: ontology term
  represented_as: node
  preferred_id: cazy
  label_in_input: cazy_family
  properties:
    name: str
    description: str
    level: int
    level_kind: str         # cazy_class | cazy_family | cazy_subfamily
    class: str              # GH | GT | PL | CE | AA | CBM
```

ID conventions: `tcdb:1.A.1.1.1` and `cazy:GH13_1`. (Bioregistry-style, matching existing `pfam:` / `kegg.brite:` conventions.)

### New edges

```yaml
gene has transporter classification:
  is_a: gene to ontology association
  represented_as: edge
  source: gene
  target: transporter class
  properties: {}

gene in cazy family:
  is_a: gene to ontology association
  represented_as: edge
  source: gene
  target: cazy family
  properties: {}

transporter class is a transporter class:
  is_a: ontology relationship
  represented_as: edge
  source: transporter class
  target: transporter class
  properties: {}

cazy family is a cazy family:
  is_a: ontology relationship
  represented_as: edge
  source: cazy family
  target: cazy family
  properties: {}
```

Per existing convention with KEGG/BRITE, `Gene_has_transporter_classification` and `Gene_in_cazy_family` edges always point at the most specific leaf, never at interior class nodes. Hierarchy traversal happens through `is_a*0..` patterns at query time.

## Adapter changes

### Extended file: `multiomics_kg/adapters/metabolism_adapter.py`

`MetabolismAdapter.get_nodes()` is extended to also yield TransporterClass and CazyFamily nodes from `scaffold_pruned.json["transporter_classes"]` and `["cazy_families"]`.

`MetabolismAdapter.get_edges()` is extended to also yield:

- `Transporter_class_is_a` edges (one per non-root TC entry, child → parent).
- `Cazy_family_is_a` edges (one per non-root CAZy entry, child → parent).

### Extended file: `multiomics_kg/adapters/functional_annotation_adapter.py`

Two new edge generators (alongside the `GeneCatalyzesReactionAdapter` from Spec 1.2):

- `GeneHasTransporterClassificationAdapter` — for each gene's `transporter_classification` list, emit one `Gene_has_transporter_classification` edge per leaf TC.
- `GeneInCazyFamilyAdapter` — for each gene's `cazy_ids` list, emit one `Gene_in_cazy_family` edge per leaf CAZy family.

Edge ID formats:
- `gene_has_tc:{locus_tag}:{tc_id}`
- `gene_in_cazy:{locus_tag}:{cazy_id}`

Both adapters get Multi* wrappers and are instantiated by `create_knowledge_graph.py`.

## Post-import (`scripts/post-import.sh` + `scripts/post-import.cypher`)

### New blocks (appended after Spec 1.2 blocks)

**Block F — `Organism_has_metabolite.transporter_count` backfill.** Sparse data: only TC entries with non-empty `substrate_classes` contribute. Phase 1.3 uses the substrate-class **text** as a soft join — if the TC's `substrate_classes` text matches a `Metabolite.synonyms` or `Metabolite.name` (case-insensitive, after normalization), the gene is counted.

```cypher
CALL { } IN TRANSACTIONS OF 5000 ROWS
MATCH (g:Gene)-[:Gene_has_transporter_classification]->(tc:TransporterClass)
WHERE size(tc.substrate_classes) > 0
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
UNWIND tc.substrate_classes AS sc
WITH g, o, toLower(sc) AS scn
MATCH (m:Metabolite)
WHERE toLower(m.name) = scn OR scn IN [s IN m.synonyms | toLower(s)]
WITH o, m, count(DISTINCT g) AS tc_genes
MERGE (o)-[e:Organism_has_metabolite]->(m)
ON CREATE SET e.as_substrate_count = 0,
              e.as_product_count   = 0,
              e.transporter_count  = tc_genes,
              e.cazy_count         = 0,
              e.evidence_sources   = ["transporter_classification"],
              e.measured_assay_count = 0,
              e.measured_compartments = [],
              e.measured_value_kinds = [],
              e.measured_paper_count = 0
ON MATCH SET  e.transporter_count  = tc_genes,
              e.evidence_sources   = apoc.coll.toSet(e.evidence_sources + "transporter_classification");
```

This is intentionally permissive — sparse data plus name-based matching is expected to recall a small but useful subset (low-double-digit-percent of metabolites). Phase 3 will replace this with explicit `TransporterClass_transports → Metabolite` edges sourced from TCDB / KEGG TRANSPORTER, at which point the sparse text join is dropped.

**Block G — `Organism_has_metabolite.cazy_count` backfill.** Same shape, with the join on CazyFamily `description` and a hardcoded substrate hint table (e.g. `GH13` family → "starch", "glycogen", "amylose"). Hardcoded table lives in `multiomics_kg/utils/cazy_substrates.py:CAZY_SUBSTRATE_HINTS` and is documented inline as a stop-gap; Phase 3 replaces with curated annotations from CAZy substrate-binding records.

```cypher
CALL { } IN TRANSACTIONS OF 5000 ROWS
MATCH (g:Gene)-[:Gene_in_cazy_family]->(cf:CazyFamily)
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WITH g, o, cf
WITH g, o, $cazy_substrate_hints[cf.id] AS hints
WHERE hints IS NOT NULL
UNWIND hints AS hint
WITH g, o, toLower(hint) AS hn
MATCH (m:Metabolite)
WHERE toLower(m.name) = hn OR hn IN [s IN m.synonyms | toLower(s)]
WITH o, m, count(DISTINCT g) AS cazy_genes
MERGE (o)-[e:Organism_has_metabolite]->(m)
ON CREATE SET e.as_substrate_count = 0, e.as_product_count = 0,
              e.transporter_count  = 0, e.cazy_count       = cazy_genes,
              e.evidence_sources   = ["cazy_family"], ...
ON MATCH SET  e.cazy_count         = cazy_genes,
              e.evidence_sources   = apoc.coll.toSet(e.evidence_sources + "cazy_family");
```

The `$cazy_substrate_hints` parameter is passed in by `post-import.sh` after loading from a small JSON file (`config/cazy_substrate_hints.json`) shipped in-repo.

**Block H — Gene routing-signal updates:**

```cypher
MATCH (g:Gene)
OPTIONAL MATCH (g)-[r:Gene_has_transporter_classification]->()
WITH g, count(r) AS tc
SET g.transporter_count = tc;

MATCH (g:Gene)
OPTIONAL MATCH (g)-[r:Gene_in_cazy_family]->()
WITH g, count(r) AS cz
SET g.cazy_count = cz;
```

### New indexes

Scalar:
- `transporter_class_level_idx` on `TransporterClass(level)`
- `transporter_class_name_idx`  on `TransporterClass(name)`
- `cazy_family_level_idx`       on `CazyFamily(level)`
- `cazy_family_class_idx`       on `CazyFamily(class)`

Full-text:
- `transporterClassFullText` on `TransporterClass(name, description, substrate_classes)`
- `cazyFamilyFullText`       on `CazyFamily(name, description)`

### Reference Cypher parity

`scripts/post-import.cypher` kept byte-identical to the bash version. Same byte-diff validation as Spec 1.2.

## KG validity tests

New file: `tests/kg_validity/test_metabolism_transport.py` (marker `@pytest.mark.kg`).

- `test_transporter_class_count` — ≥ 500 nodes.
- `test_cazy_family_count` — ≥ 100 nodes.
- `test_transporter_hierarchy_connected` — every non-root TC has a parent that exists; every root has `level = 0` and `parent IS NULL`.
- `test_cazy_hierarchy_connected` — same shape; class set ⊆ {GH, GT, PL, CE, AA, CBM}.
- `test_gene_tc_targets_exist` — every `Gene_has_transporter_classification` targets an existing TransporterClass.
- `test_gene_tc_targets_are_leaves` — every target TC has no child via `is_a*1..` (or equivalently, no other TC has it as `parent`).
- `test_gene_cazy_targets_exist` — every `Gene_in_cazy_family` targets an existing CazyFamily.
- `test_gene_cazy_targets_are_leaves` — same shape.
- `test_organism_has_metabolite_transporter_backfill` — at least 50 `Organism_has_metabolite` edges have `transporter_count > 0` (sanity check; the join is sparse but non-empty).
- `test_organism_has_metabolite_cazy_backfill` — at least 20 with `cazy_count > 0`.
- `test_gene_routing_signals` — for every gene, `transporter_count = count((g)-[:Gene_has_transporter_classification]->())` and same for cazy.

Snapshot fixture regenerated: include 5 spot TransporterClass and 3 CazyFamily entries.

## Estimated KG impact

| Quantity | Estimate |
|---|---|
| TransporterClass nodes | ~1,500 |
| CazyFamily nodes | ~400 |
| `Gene_has_transporter_classification` edges | 15K–40K |
| `Gene_in_cazy_family` edges | 10K–30K |
| `Transporter_class_is_a` edges | ~1,500 |
| `Cazy_family_is_a` edges | ~400 |
| Backfilled `Organism_has_metabolite.transporter_count > 0` edges | ~500–2,000 |
| Backfilled `Organism_has_metabolite.cazy_count > 0` edges | ~200–1,000 |
| **Total new edges** | **~30K–75K** |

## Acceptance criteria

- All Phase 1.3 nodes and edges present in the deployed graph at expected magnitudes (within 50% of the estimates table).
- New KG validity tests pass.
- Backfill of `Organism_has_metabolite.transporter_count` and `cazy_count` is non-zero for both Prochlorococcus and Alteromonas (sanity check that the sparse join produces something).
- The Phase 1.2 driver query still runs and returns the same or larger candidate set (TC/CAZy fields are additive — they don't filter anything out).
- `omics-edge-snapshot` confirms no regression in existing expression-edge counts.
- `scripts/post-import.cypher` and `scripts/post-import.sh` produce byte-identical computed properties.
- CLAUDE.md updated with the new node/edge types and post-import additions.
- `test_import_report.py` still green (no rejected relationships).

## Open question — sparse text-join for substrate matching

The TC/CAZy → Metabolite text-join (Blocks F + G) is intentionally a stop-gap. The cleaner Phase 3 design is direct `TransporterClass_transports → Metabolite` edges sourced from TCDB and KEGG TRANSPORTER's explicit substrate-compound links. Phase 1.3 lands the text-join because:

1. The Phase-3 source data is sparse and would need its own download + parse pipeline.
2. The driver query's value is from `as_substrate_count` + `as_product_count`, not from `transporter_count`. The TC/CAZy backfill is augmentation, not foundation.

If the text-join produces more noise than signal in practice (e.g. spurious matches on common substrate words like "sugar"), drop Blocks F + G entirely and let `transporter_count` / `cazy_count` stay at 0 until Phase 3. Decision deferred until first deployment.

## Out of scope (handled in later phases)

- Direct `TransporterClass_transports → Metabolite` edges → Phase 3.
- Direct `CazyFamily_acts_on → Metabolite` edges → Phase 3.
- Aggregated `substrate_classes` at interior TC levels → Phase 3.
- MetaboliteAssay integration → Phase 2.
