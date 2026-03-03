# Methods: Heterogeneous Gene Identifier Harmonization

## Overview

Integrating transcriptomic, proteomic, and metabolomic data from diverse publications into a unified knowledge graph requires resolving the heterogeneous gene identifiers used across studies. Different research groups and time periods use distinct identifier vocabularies: NCBI locus tags, JGI Integrated Microbial Genomes (IMG) catalog IDs, microarray probeset IDs, old-style locus tags from earlier genome annotations, UniProt entry names, RefSeq protein accessions, and human-readable gene symbols. We developed a three-tier identifier harmonization system that builds a pre-computed mapping for each genome strain and resolves publication CSV files to canonical locus tags prior to knowledge graph construction.

## Data Sources

Gene identifier information was collected from four categories of sources per strain:

1. **NCBI genome annotations** (GFF3 format): Provided locus tags, NCBI RefSeq locus tags (RS-format), old locus tags from earlier annotation versions, RefSeq protein accessions (WP_ format), and gene symbols.

2. **Cyanorak database annotations** (GFF/GBK format): Provided Cyanorak-specific locus tag identifiers and cross-references for *Prochlorococcus* and *Synechococcus* strains.

3. **UniProt proteome records**: Provided UniProt accession numbers and entry names (e.g., `DNAA_PROM0`) linked to RefSeq protein accessions.

4. **Publication-specific identifier translation tables** (`id_translation` entries in paperconfig.yaml): Provided JGI IMG OID numbers, microarray probeset IDs, and annotation tables mapping non-standard identifiers used by individual research groups to cross-reference identifiers present in the genome annotation.

## Three-Tier Identifier Classification

Gene identifiers were classified into three tiers based on their expected uniqueness within a single organism:

**Tier 1 (gene-unique identifiers)**: Expected to uniquely identify exactly one gene within an organism. Includes: `locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `old_locus_tag`, `alternative_locus_tag`, `jgi_id` (JGI IMG catalog OIDs), `probeset_id` (microarray probeset identifiers), and `uniprot_entry_name` (with the `_ORGANISM` suffix stripped; e.g., `DNAA_PROM0` → `DNAA`). Tier 1 identifiers participate in the convergence graph as one-to-one mappings. When two different locus tags share a Tier 1 identifier, this is treated as a data error and recorded in a `conflicts` registry for manual inspection.

**Tier 2 (protein-level identifiers)**: Can legitimately match multiple paralogous genes within one organism. Includes: `protein_id` (RefSeq WP_ accessions) and `uniprot_accession`. Paralogs sharing a single protein accession are biologically expected for gene families with conserved protein sequences. These identifiers participate in the convergence graph for transitive linking but are stored in a `multi_lookup` registry (one-to-many) rather than the one-to-one `specific_lookup`. A Tier 2 identifier is used for resolution only when it maps to exactly one gene in the organism (singleton case), permitting unambiguous assignment.

**Tier 3 (generic identifiers)**: Many-to-many by design, such as `gene_name` (e.g., `dnaA`), `gene_synonym`, and `gene_oln`. These are added to `multi_lookup` directly during the seeding phase and are similarly used for resolution only as singletons.

## Iterative Convergence Graph

For each genome strain, a `GeneIdGraph` was constructed and populated in two phases:

**Seeding phase**: For each gene in `gene_annotations_merged.json` (a per-strain annotation table built from NCBI, Cyanorak, and UniProt data), the canonical locus tag was registered as an anchor node. All alternative identifiers were added to the graph according to their tier classification.

**Paper-source processing phase**: Each row from all `id_translation`, `annotation_gff`, and `csv` supplementary table entries in the paperconfig files was processed. A row was represented as a list of (identifier value, identifier type) pairs. The graph searched for an anchor (canonical locus tag) reachable from any identifier in the row, using a priority order: Tier 1 direct match → whitespace-token split (for compound cells such as `"dnaA PMM0001"`) → Tier 2 singleton. When an anchor was found, all remaining identifiers in that row were linked to that locus tag.

All paper sources were processed simultaneously in iterative passes until no new identifier-to-locus-tag mapping was added (convergence). Typically two to three passes were required. This approach is order-independent: if identifier A from paper X co-occurs with identifier B in one row, and identifier B co-occurs with identifier C in another row (possibly from a different paper), then A, B, and C all resolve to the same locus tag after convergence, regardless of the order in which the papers appear in the configuration.

The resulting mappings were stored in a per-strain `gene_id_mapping.json` file (schema version 2) containing: a `specific_lookup` dictionary (Tier 1 alt-IDs → locus tag, 1:1), a `multi_lookup` dictionary (Tier 2+3 IDs → list of locus tags), and a `conflicts` dictionary (Tier 1 IDs with conflicting assignments).

## Row-Level Resolution Strategy

For each row in a publication's differential expression CSV file, locus tag resolution was performed using the `resolve_row` function, which attempts the following passes in order across all configured identifier columns (primary `name_col` first, then any explicit `id_columns` from the paperconfig):

**Pass 1 — Tier 1 specific lookup with list expansion**: Cell values were split on commas and semicolons (for list-valued cells such as `"PMM0001, PMM0002"`), and each token was looked up in `specific_lookup`. Normalization heuristics were also applied: trailing asterisks (footnote artifacts) were stripped, and numeric suffixes were zero-padded when the original failed to match (e.g., `MIT1002_0001` → `MIT1002_00001`). The first successful match was returned.

**Pass 2 — Tier 2+3 multi-lookup, singletons only**: Each expanded token was looked up in `multi_lookup`. A match was accepted only when the identifier resolved to exactly one gene within the organism. When `multi_lookup` returned multiple candidates (ambiguous), the next configured column was tried before the row was flagged.

Each resolved row was annotated with a `resolution_method` string (e.g., `tier1:Gene`, `heuristic:Gene`, `multi:JGI_ID`) describing which pass and column produced the result. Unresolved rows were annotated with a reason code (`unresolved`, `ambiguous`, `tier1_conflict`) and retained in the output file with a null locus tag, ensuring no data was silently discarded. Per-table diagnostic reports listing all unresolved rows with their raw values and attempted columns were written alongside each resolved CSV.

## Diagnostic Output and Quality Monitoring

After building `gene_id_mapping.json`, a diagnostic report (`gene_id_mapping_report.json`) was written for each strain. This report enumerates, for each identifier type, the fraction of identifiers that resolved uniquely (in `specific_lookup`), matched multiple genes (in `multi_lookup`), or produced Tier 1 conflicts. Reclassification warnings were emitted when an identifier type declared as Tier 1 exhibited a multi-match rate exceeding 10%, suggesting the type should be moved to Tier 2.

Regression testing used a snapshot-based approach: `Affects_expression_of` edge counts were recorded per publication before and after each pipeline change. Any publication showing a reduction in edge count triggered investigation of the corresponding `_resolved_report.txt` and re-examination of the paperconfig `id_columns` configuration.

## Phased Validation

The redesigned pipeline was validated strain by strain, beginning with *Prochlorococcus* MIT9312. For each strain, the pipeline proceeded as: (1) rebuild `gene_id_mapping.json` from all annotation sources and paperconfigs; (2) re-resolve all publication CSVs; (3) compare expression edge counts in the knowledge graph against the pre-rebuild snapshot; (4) investigate any publication showing reduced edge counts. This phased approach ensured that improvements in identifier resolution were verified before proceeding to the next strain.
