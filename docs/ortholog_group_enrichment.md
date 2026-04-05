# OrthologGroup Node Enrichment

## Summary

OrthologGroup nodes carry pre-computed aggregate properties derived from their member genes, plus functional descriptions from eggNOG and derived annotations. These enable richer queries and fulltext search without traversing membership edges.

## Properties on OrthologGroup Nodes

| Property | Type | Source | Description | Example |
|---|---|---|---|---|
| `consensus_product` | string | Member genes | Most common product (excludes "hypothetical protein" when possible) | `"DNA polymerase III, beta subunit"` |
| `consensus_gene_name` | string | Member genes | Most frequent non-null gene name | `"dnaN"` |
| `description` | string (nullable) | eggNOG local DB | Functional narrative from eggNOG 5.0. Null for Cyanorak groups. | `"DNA-directed DNA polymerase activity"` |
| `functional_description` | string (nullable) | Derived from member gene roles/categories | Majority-vote CyanorakRole (hierarchical) + CogFunctionalCategory names | `"DNA metabolism > DNA replication, recombination, and repair; Replication, recombination and repair"` |
| `member_count` | integer | Member genes | Total number of genes in the group | `12` |
| `organism_count` | integer | Member genes | Number of distinct organism strains represented | `8` |
| `genera` | string[] | Member genes | Sorted list of distinct genera | `["Alteromonas", "Prochlorococcus"]` |
| `has_cross_genus_members` | string | Member genes | `"cross_genus"` or `"single_genus"` | `"cross_genus"` |

### `description` data source

Loaded from the local eggNOG SQLite database (`$EGGNOG_DATA_DIR/eggnog.db`). The adapter parses the OG node ID (e.g., `eggnog:COG0592@2` -> og=`COG0592`, level=`2`) and queries the `og` table. Cyanorak groups have null description (no external source).

### `functional_description` derivation

Built from member genes' `cyanorak_Role` and `cog_category` fields (both lists). Majority-vote: a code passes if >50% of member genes have it. Cyanorak roles use full hierarchical names (e.g., `"Photosynthesis and respiration > Photosystem II"`). Filtered out: COG "S" (Function unknown), Cyanorak "R"/"R.2"/"R.4" (hypothetical). "Other" as leaf is kept when parent is informative.

## Annotation Edges

| Type | Direction | Description |
|------|-----------|-------------|
| `Og_has_cyanorak_role` | OrthologGroup -> CyanorakRole | Created if >50% of member genes have that Cyanorak role |
| `Og_in_cog_category` | OrthologGroup -> CogFunctionalCategory | Created if >50% of member genes have that COG category |

## Indexes

| Index | Type | Fields |
|---|---|---|
| `ortholog_group_rank_idx` | scalar | `OrthologGroup.specificity_rank` |
| `orthologGroupFullText` | fulltext | `OrthologGroup(consensus_product, consensus_gene_name, description, functional_description)` |

## Example Queries

### Fulltext search across all text fields
```cypher
CALL db.index.fulltext.queryNodes('orthologGroupFullText', 'photosynthesis')
YIELD node AS og, score
RETURN og.id AS group_id, og.consensus_gene_name,
       og.consensus_product, og.description,
       og.functional_description, score
ORDER BY score DESC
LIMIT 10
```

### Traverse from group to functional category
```cypher
MATCH (og:OrthologGroup {id: 'cyanorak:CK_00000570'})-[:Og_has_cyanorak_role]->(cr:CyanorakRole)
RETURN cr.id AS role_id, cr.name AS role_name
```

### Find ortholog groups that bridge Prochlorococcus and Alteromonas
```cypher
MATCH (og:OrthologGroup)
WHERE og.has_cross_genus_members = 'cross_genus'
RETURN og.name, og.consensus_product, og.consensus_gene_name,
       og.member_count, og.organism_count, og.genera
ORDER BY og.member_count DESC
LIMIT 20
```

### Find large conserved groups with known function
```cypher
MATCH (og:OrthologGroup)
WHERE og.member_count > 10
  AND og.consensus_product IS NOT NULL
  AND og.consensus_product <> 'hypothetical protein'
RETURN og.name, og.source, og.consensus_product, og.member_count
ORDER BY og.member_count DESC
```

### Find all groups for a specific gene and their consensus annotations
```cypher
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup)
RETURN og.name, og.source, og.taxonomic_level, og.specificity_rank,
       og.consensus_product, og.consensus_gene_name,
       og.member_count, og.organism_count, og.has_cross_genus_members
```

### Filter by specificity rank (curated > family > order > domain)
```cypher
MATCH (og:OrthologGroup)
WHERE og.specificity_rank <= 1
  AND og.organism_count >= 5
RETURN og.name, og.consensus_gene_name, og.organism_count, og.genera
ORDER BY og.organism_count DESC
```

## Existing Properties (unchanged)

| Property | Type | Description |
|---|---|---|
| `name` | string | Raw OG identifier (e.g., `"CK_00000364"`, `"COG0592@2"`) |
| `source` | string | `"cyanorak"` or `"eggnog"` |
| `taxonomic_level` | string | `"curated"`, `"Prochloraceae"`, `"Bacteria"`, etc. |
| `taxon_id` | integer | NCBI taxon ID of the taxonomic level |
| `specificity_rank` | integer | 0=curated, 1=family, 2=order, 3=domain |
