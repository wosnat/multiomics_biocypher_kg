# KG Changes for Category/Role Ontology Support

Spec for `multiomics_biocypher_kg` changes needed to support COG category,
Cyanorak role, and TIGR role in the explorer's ontology MCP tools.

This is the KG-side companion to the explorer's
`plans/redefine_mcp_tools/add_category_role_ontologies.md`.

---

## 1. Rename `description` → `name` on CyanorakRole and TigrRole nodes

**Problem:** All existing ontology nodes (BiologicalProcess, MolecularFunction,
CellularComponent, EcNumber, KeggTerm) use `name` as their human-readable text
property. CyanorakRole and TigrRole use `description` instead. CogFunctionalCategory
already has `name`.

**Fix:** Rename `description` → `name` on CyanorakRole and TigrRole nodes.

**Current state:**
```
CogFunctionalCategory: {id, code, name, preferred_id}              ✓ has name
CyanorakRole:          {id, code, description, preferred_id}        ✗ no name
TigrRole:              {id, code, description, preferred_id}        ✗ no name
```

**After fix:**
```
CogFunctionalCategory: {id, code, name, preferred_id}
CyanorakRole:          {id, code, name, preferred_id}
TigrRole:              {id, code, name, preferred_id}
```

### Files to change

- `config/schema_config.yaml`: `description: str` → `name: str` for both node types
- `multiomics_kg/adapters/functional_annotation_adapter.py`: property dict key `"description"` → `"name"` in CyanorakRole and TigrRole node emission
- `tests/test_cog_role_annotation_adapter.py`: update property assertions
- `tests/kg_validity/test_functional_annotation.py`: update property-presence checks

## 2. Create fulltext indexes

The explorer's `search_ontology` tool uses Neo4j fulltext indexes for Lucene-syntax
search. No indexes exist for these three node types.

Add to `scripts/post-import.sh` and `scripts/post-import.cypher`:

```cypher
CREATE FULLTEXT INDEX cogCategoryFullText IF NOT EXISTS
  FOR (n:CogFunctionalCategory) ON EACH [n.name];

CREATE FULLTEXT INDEX cyanorakRoleFullText IF NOT EXISTS
  FOR (n:CyanorakRole) ON EACH [n.name];

CREATE FULLTEXT INDEX tigrRoleFullText IF NOT EXISTS
  FOR (n:TigrRole) ON EACH [n.name];
```

## 3. Add TigrRole hierarchy edges

**Problem:** CyanorakRole has `Cyanorak_role_is_a_cyanorak_role` hierarchy edges
(154 edges linking sub-roles to main roles). TigrRole has the same main/sub
structure encoded in its descriptions but no hierarchy edges.

**Current TigrRole structure** (~114 nodes, all flat):
```
tigr.role:100  "Central intermediary metabolism / Amino sugars"
tigr.role:102  "Central intermediary metabolism / Other"
tigr.role:108  "Energy metabolism / Aerobic"
tigr.role:120  "Energy metabolism / TCA cycle"
```

The text before ` / ` is the main category (19 unique), after is the sub-category.

**Fix:** Create main-category TigrRole nodes and `Tigr_role_is_a_tigr_role` edges:

1. Parse unique main categories from the ` / ` split of existing descriptions
2. Create a TigrRole node for each main category (e.g. `tigr.main:energy_metabolism`)
   with `name` = "Energy metabolism"
3. Add `Tigr_role_is_a_tigr_role` edge from each sub-role to its main category

**Example result:**
```
(tigr.role:120 "Energy metabolism / TCA cycle")
  -[:Tigr_role_is_a_tigr_role]->
(tigr.main:energy_metabolism "Energy metabolism")
```

### Files to change

- `config/schema_config.yaml`: add `tigr role is a tigr role` edge type
- `multiomics_kg/adapters/functional_annotation_adapter.py`:
  - `_all_tigr_codes()`: also collect main categories from ` / ` split
  - `get_nodes()`: emit main-category TigrRole nodes with `tigr.main:` prefix
  - `get_edges()`: emit `tigr_role_is_a_tigr_role` hierarchy edges
- `tests/test_cog_role_annotation_adapter.py`: test TIGR hierarchy
- `tests/kg_validity/test_functional_annotation.py`: test TIGR hierarchy edges

## 4. Explorer integration

Suggested `ONTOLOGY_CONFIG` entries:

```python
"cog_category": {
    "label": "CogFunctionalCategory",
    "gene_rel": "Gene_in_cog_category",
    "hierarchy_rels": [],          # flat — 26 single-letter categories
    "fulltext_index": "cogCategoryFullText",
},
"cyanorak_role": {
    "label": "CyanorakRole",
    "gene_rel": "Gene_has_cyanorak_role",
    "hierarchy_rels": ["Cyanorak_role_is_a_cyanorak_role"],
    "fulltext_index": "cyanorakRoleFullText",
},
"tigr_role": {
    "label": "TigrRole",
    "gene_rel": "Gene_has_tigr_role",
    "hierarchy_rels": ["Tigr_role_is_a_tigr_role"],
    "fulltext_index": "tigrRoleFullText",
},
```

### Node ID formats

- COG: `cog.category:J` (single letter)
- CyanorakRole: `cyanorak.role:B.5.1` (hierarchical code)
- TigrRole sub-roles: `tigr.role:120` (numeric code from data)
- TigrRole main categories: `tigr.main:energy_metabolism` (slugified name)

### Expected counts

| Entity | Count |
|--------|-------|
| CogFunctionalCategory nodes | 25 (hardcoded) |
| CyanorakRole nodes | ~172 (full tree) |
| TigrRole nodes (sub-roles) | ~114 (from gene data) |
| TigrRole nodes (main categories) | ~19 (from ` / ` split) |
| Gene_in_cog_category edges | ~20,000 |
| Gene_has_cyanorak_role edges | ~10,000 |
| Gene_has_tigr_role edges | ~10,000 |
| Cyanorak_role_is_a_cyanorak_role edges | ~154 |
| Tigr_role_is_a_tigr_role edges | ~106 |

## Verification

After rebuild:

```cypher
-- Step 1: name property exists on all role/category nodes
MATCH (t:CyanorakRole) WHERE t.name IS NULL RETURN count(t)
-- Expected: 0

MATCH (t:TigrRole) WHERE t.name IS NULL RETURN count(t)
-- Expected: 0

-- Step 2: fulltext indexes work
CALL db.index.fulltext.queryNodes('cogCategoryFullText', 'energy')
YIELD node, score RETURN node.name, score LIMIT 5

CALL db.index.fulltext.queryNodes('cyanorakRoleFullText', 'DNA')
YIELD node, score RETURN node.name, score LIMIT 5

CALL db.index.fulltext.queryNodes('tigrRoleFullText', 'metabolism')
YIELD node, score RETURN node.name, score LIMIT 5

-- Step 3: TIGR hierarchy
MATCH (child:TigrRole)-[:Tigr_role_is_a_tigr_role]->(parent:TigrRole)
RETURN parent.name, count(child) AS sub_roles
ORDER BY sub_roles DESC
```
