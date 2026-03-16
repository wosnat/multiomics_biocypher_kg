# OrthologGroup Node Enrichment

## Summary

OrthologGroup nodes now carry pre-computed aggregate properties derived from their member genes. These properties enable richer queries without traversing membership edges for basic group characterization.

## New Properties on OrthologGroup Nodes

| Property | Type | Description | Example |
|---|---|---|---|
| `consensus_product` | string | Most common product annotation among member genes (excludes "hypothetical protein" when possible) | `"DNA polymerase III, beta subunit"` |
| `consensus_gene_name` | string | Most frequent non-null gene name among members | `"dnaN"` |
| `member_count` | integer | Total number of genes in the group | `12` |
| `organism_count` | integer | Number of distinct organism strains represented | `8` |
| `genera` | string[] | Sorted list of distinct genera | `["Alteromonas", "Prochlorococcus"]` |
| `has_cross_genus_members` | string | `"cross_genus"` or `"single_genus"` | `"cross_genus"` |

## New Index

| Index | Property |
|---|---|
| `ortholog_group_rank_idx` | `OrthologGroup.specificity_rank` |

## Example Queries

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
