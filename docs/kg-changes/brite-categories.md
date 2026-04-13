# KEGG BRITE Functional Hierarchies (BriteCategory nodes)

**Added:** 2026-04-13

## What changed

Added `BriteCategory` nodes and two new edge types to the knowledge graph, ingesting 12 KEGG BRITE functional hierarchy trees. This provides hierarchical functional classification for ~1,800 KO terms that are present in the graph but have no KEGG pathway membership ("pathway-orphan KOs").

## New node type: BriteCategory

Node IDs use a positional scheme: `kegg.brite:{tree_code}.{A#}.{B#}...`

Example: `kegg.brite:ko02000.A1.B1.C1` (Phosphate transport system within ABC Transporters)

Properties:
- `name` (str) — category label from KEGG
- `tree` (str) — snake_case canonical tree name (e.g. `"transporters"`)
- `tree_code` (str) — KEGG tree ID (e.g. `"ko02000"`)
- `level` (int) — hierarchy depth, 0 = broadest (no parent edge emitted)
- `level_kind` (str) — `brite_class` | `brite_subclass` | `brite_family` | `brite_subfamily`
- `member_ko_count` (int, post-import) — recursive count of KO leaves under this node
- `gene_count` (int, post-import) — distinct genes reachable via KO leaves
- `organism_count` (int, post-import) — distinct organisms among those genes

## New edge types

| Edge type | Source | Target | Count |
|-----------|--------|--------|-------|
| `Brite_category_is_a_brite_category` | BriteCategory (child) | BriteCategory (parent) | ~8,400 |
| `Kegg_term_in_brite_category` | KeggTerm (KO level) | BriteCategory | ~17,000 |

## 12 configured trees

| Tree code | Name |
|-----------|------|
| ko01000 | enzymes |
| ko02000 | transporters |
| ko01002 | peptidases |
| ko03000 | transcription_factors |
| ko02044 | secretion |
| ko02022 | two_component |
| ko02048 | defense |
| ko03110 | chaperones |
| ko03011 | ribosome |
| ko03012 | translation_factors |
| ko03016 | trna_biogenesis |
| ko03032 | dna_replication |

## Level semantics

No synthetic tree-root node is emitted. A-level entries are the broadest nodes:

| BRITE level | `level` int | `level_kind` | Has parent edge |
|-------------|-------------|--------------|-----------------|
| A-level | 0 | `brite_class` | No |
| B-level | 1 | `brite_subclass` | Yes |
| C-level | 2 | `brite_family` | Yes |
| D-level (non-KO) | 3 | `brite_subfamily` | Yes |
| KO leaf | — | (not a node) | (edge only) |

## Implementation files

- `multiomics_kg/utils/brite_utils.py` — tree download/cache, `compute_level_kind`
- `multiomics_kg/adapters/brite_adapter.py` — `MultiBriteAdapter`
- `config/schema_config.yaml` — schema entries for `brite category`, `brite category is a brite category`, `kegg term in brite category`
- `scripts/post-import.sh` / `scripts/post-import.cypher` — indexes + computed properties
- `tests/test_brite_utils.py` — 16 unit tests for utilities
- `tests/test_brite_adapter.py` — 29 unit tests for adapter
- `tests/kg_validity/test_brite.py` — 13 KG validity tests

## Querying

```cypher
// Find all transporters classification categories for pstB
MATCH (g:Gene {locus_tag: 'PMM1416'})-[:Gene_has_kegg_ko]->(ko:KeggTerm)
MATCH (ko)-[:Kegg_term_in_brite_category]->(b:BriteCategory)
RETURN ko.name, b.name, b.level, b.level_kind, b.tree

// Find all genes in the ribosome tree
MATCH (b:BriteCategory {tree_code: 'ko03011'})
MATCH (ko:KeggTerm)-[:Kegg_term_in_brite_category]->(b)
MATCH (g:Gene)-[:Gene_has_kegg_ko]->(ko)
RETURN DISTINCT g.locus_tag, g.organism_name, b.name

// Navigate up the hierarchy from a specific category
MATCH (leaf:BriteCategory {name: 'Phosphate transport system'})
MATCH path = (leaf)-[:Brite_category_is_a_brite_category*0..]->(root:BriteCategory {level: 0})
RETURN [n IN nodes(path) | n.name] AS hierarchy
```
