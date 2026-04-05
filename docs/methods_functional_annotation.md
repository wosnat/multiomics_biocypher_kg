# Methods: Functional Annotation as Graph Edges

## Overview

Gene nodes in the knowledge graph store functional annotations (GO terms, EC numbers, KEGG orthology assignments, COG categories, protein domain identifiers, and role classifications) as array properties derived from multiple annotation sources. While property-based storage permits fast filtering (e.g., `WHERE 'GO:0005737' IN g.go_terms`), it does not support graph traversal. An LLM agent asking "which genes are involved in DNA replication" cannot follow edges from a Gene Ontology term to the genes it annotates; instead, it must know to search within array properties on gene nodes, which is fragile and limits cross-annotation queries.

To address this, functional annotations were converted from gene node properties into explicit graph edges connecting genes to ontology and classification nodes. The original array properties were retained for backward compatibility and fast filtering, while the new edges enable traversal-based queries across annotation systems.

## Annotation Sources

Functional annotations originate from three upstream sources, with coverage varying by organism:

**Cyanorak database** (Doré et al. 2020): For *Prochlorococcus* and *Synechococcus* strains, the Cyanorak v2.1 database provides expert-curated COG orthologous group assignments, a three-level functional role hierarchy specific to cyanobacteria (Cyanorak Roles), and automated TIGR role codes from the IGS Annotation Engine. These annotations are embedded in the per-strain `gene_annotations_merged.json` during the genome data preparation pipeline.

**eggNOG-mapper v2**: Automated functional annotation applied to all 13 strains, producing COG functional categories, KEGG Orthology (KO) assignments, KEGG pathway memberships, GO terms, EC numbers, and Pfam domain identifiers. For *Prochlorococcus* and *Synechococcus* strains, eggNOG-mapper supplements the Cyanorak annotations; for the three *Alteromonas* strains (MIT1002, EZ55, HOT1A3), which are not represented in Cyanorak, eggNOG-mapper is the primary source of functional classification.

**UniProt**: GO terms and EC numbers from UniProt proteome records, merged during the gene annotation build step.

When annotations from multiple sources overlap, Cyanorak data takes priority for *Prochlorococcus*/*Synechococcus* strains, as it reflects expert manual curation. eggNOG-mapper fills annotation gaps not covered by Cyanorak. For *Alteromonas* strains, eggNOG-mapper and UniProt are the sole sources.

## Graph Modeling

Four adapter modules convert the per-gene annotation arrays into graph nodes and edges. Each follows the project's Multi* wrapper pattern: a per-strain adapter reads `gene_annotations_merged.json` and yields edges, while a multi-strain wrapper aggregates nodes across all strains (deduplicating by identifier) and yields hierarchy edges where applicable.

### Gene Ontology

GO terms referenced in gene annotations were resolved against the GO OBO ontology to determine their namespace (biological process, cellular component, or molecular function), and edges were created with namespace-specific relationship types: `gene_involved_in_biological_process`, `gene_located_in_cellular_component`, and `gene_enables_molecular_function`.

Rather than importing the full Gene Ontology (~47,000 terms), only the subset of terms directly referenced in at least one gene annotation was collected across all strains, and an ancestry closure was computed: all transitive parent terms reachable via `is_a` and `part_of` relationships were added. This produces a compact GO subgraph (typically ~5,000 terms) sufficient for hierarchical queries such as "all genes annotated with any descendant of GO:0006259 (DNA metabolic process)" while avoiding the overhead of unused ontology branches. The GO-GO hierarchy edges (including `is_a`, `part_of`, `regulates`, `positively_regulates`, and `negatively_regulates` relationships) were emitted for terms within this closure.

### Enzyme Commission Numbers

EC number nodes were imported from the full Expasy enzyme classification hierarchy, with `ec_number_is_a_ec_number` edges encoding the four-level numerical hierarchy (e.g., EC:2.7.7.7 → EC:2.7.7.- → EC:2.7.- → EC:2.-). Gene-to-EC edges (`gene_catalyzes_ec_number`) were created from the `ec_numbers` array in gene annotations.

### KEGG Orthology and Pathway Hierarchy

KEGG annotations were modeled as a unified `KeggTerm` node type with a `level` property discriminating four hierarchy levels, retrieved from the KEGG REST API and cached locally:

1. **ko** (KO, e.g., K02338) — the functional ortholog unit
2. **pathway** (e.g., ko00230 "Purine metabolism") — metabolic or signaling pathways
3. **subcategory** (BRITE B-level, e.g., "Nucleotide metabolism") — functional groupings of pathways
4. **category** (BRITE A-level, e.g., "Metabolism") — top-level functional classes

Gene-to-KO edges (`Gene_has_kegg_ko`) were created from the `kegg_ko` array in gene annotations. The hierarchy edges (`Kegg_term_is_a_kegg_term`) were derived from the KEGG BRITE orthology classification. Pathway names are extracted from the BRITE hierarchy C-level node labels. Only KO entries, pathways, subcategories, and categories reachable from KOs present in the gene annotations were included.

### COG Functional Categories, Cyanorak Roles, and TIGR Roles

Three classification systems with differing scope and granularity were modeled:

**COG functional categories** (25 single-letter codes, e.g., "L" for "Replication, recombination and repair") provide a coarse functional classification present across all 13 strains. The 25 category nodes are hardcoded; `gene_in_cog_category` edges were created from the `cog_category` array. COG categories are the only functional classification that spans both cyanobacteria and heterotrophic bacteria in the knowledge graph, making them the natural basis for cross-phylum functional comparisons.

**Cyanorak Roles** are a cyanobacteria-specific three-level functional hierarchy maintained by the Cyanorak database. The hierarchy contains approximately 172 nodes organized as a tree (e.g., "B" → "B.5" Pigments → "B.5.1" Carotenoids). Parent codes are derivable from child codes by stripping the last dot-separated component. The full tree was parsed from the Cyanorak role definition file and imported as `cyanorak_role` nodes with `cyanorak_role_is_a_cyanorak_role` hierarchy edges, enabling queries at any level of specificity (e.g., "all genes involved in cofactor biosynthesis" traverses from root "B" to all descendants). Cyanorak Roles exist only for *Prochlorococcus* and *Synechococcus* strains; *Alteromonas* strains silently yield no edges for this classification.

**TIGR Roles** are automated functional role codes assigned by the IGS Annotation Engine during the original genome annotation. They are modeled as flat nodes (no hierarchy) because the role codes are opaque integers without a derivable parent structure. TIGR Roles are retained despite the availability of the more refined Cyanorak Roles because they cover approximately 13% of *Prochlorococcus*/*Synechococcus* genes that lack a Cyanorak Role assignment, and they represent the upstream automated annotation that Cyanorak experts subsequently refined.

## Dual Representation: Properties and Edges

Functional annotations are stored both as array properties on gene nodes and as explicit graph edges. This dual representation serves two complementary access patterns:

- **Property-based filtering** is efficient for queries that test membership within a known gene set (e.g., `WHERE 'K02338' IN g.kegg_ko`). LLM agents generating Cypher queries naturally produce this pattern.
- **Edge-based traversal** is necessary for queries that start from a functional category and seek genes (e.g., "which genes are in KEGG pathway ko03030?"), for hierarchical queries (e.g., "all genes in any sub-pathway of Metabolism"), and for cross-annotation joins (e.g., "genes with both COG category L and GO:0006260").

The two representations are generated from the same underlying data (`gene_annotations_merged.json`) and are therefore consistent by construction.

Not all annotation types warrant edge modeling. Pfam protein domain identifiers, COG orthologous group IDs (COG####), and eggNOG taxonomic-level OGs are retained as array properties on gene nodes but are not converted to edges. For these annotations, the most common query patterns (e.g., "does this gene have domain PF00712?") are well served by property filtering, and the cross-gene grouping queries that edge traversal would enable (e.g., "all genes sharing COG0592") are lower priority given that five classification hierarchies already provide multiple axes for functional grouping.

## Coverage and Scope

Coverage varies across annotation systems and organisms. COG functional categories provide the broadest coverage (present for nearly all protein-coding genes across all 13 strains), while Cyanorak Roles and TIGR Roles are limited to the 8 *Prochlorococcus* and 2 *Synechococcus*/*Parasynechococcus* strains present in the Cyanorak database. GO term annotations are generally richer for *Alteromonas* strains (where eggNOG-mapper and UniProt contribute larger GO term sets typical of well-studied Gammaproteobacteria) than for *Prochlorococcus* strains.

The number of nodes created per classification system reflects the scope of each hierarchy:

- **GO**: ~5,000 terms (ancestry closure of referenced terms)
- **EC**: Full Expasy hierarchy (~8,000 nodes including incomplete entries)
- **KEGG**: ~2,000 KO nodes, ~400 pathways, plus subcategories and categories reachable from annotated KOs
- **COG functional categories**: 25 nodes (fixed)
- **Cyanorak Roles**: ~172 nodes (full tree)
- **TIGR Roles**: Variable (only codes present in at least one gene annotation)
