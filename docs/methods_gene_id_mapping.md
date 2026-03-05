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

## Cross-Assembly Identifier Bridging

Some organisms in the knowledge graph lack a single canonical genome assembly. When different publications used gene annotations from different assemblies or annotation pipelines, identifiers from each study belong to independent namespaces that share no common IDs. In these cases, a cross-assembly bridge must be constructed before the standard three-tier resolution can proceed.

### Case study: *Alteromonas macleodii* EZ55

EZ55 exemplifies this challenge. Two independent annotations exist:

1. **IMG/JGI annotation** (genome 2785510739): A draft assembly consisting of three scaffolds (`scf7180000008833`, `scf7180000008834`, `scf7180000008835`), totaling approximately 4.9 Mb with 5,016 predicted genes. Gene identifiers use the prefix `AEZ55_` with 4-digit numbering (e.g., `AEZ55_0520`). This annotation was used by Hennon et al. (2017), who aligned reads to the draft genome (SRA accession SRX022631; BioProject PRJNA50031, sequenced at JCVI). The draft assembly was never deposited in GenBank as a genome record — only the raw reads exist in NCBI.

2. **NCBI RefSeq annotation** (GCF_901457815.2 / GCA_901457815.2): A complete genome assembly submitted to ENA by the Institute for Chemistry and Biology of the Marine Environment, comprising two contigs (OX359237.1, 4.7 Mb chromosome; OX359238.1, 203 kb plasmid) with 4,324 annotated genes. The GCA annotation uses `EZ55_NNNNN` (5-digit) locus tags and was originally produced by IMG (`source:img_core_v400` in CDS Note fields). NCBI PGAP reannotated the assembly with `ALTBGP6_RS*` locus tags (GCF). This annotation was used by Barreto et al. (2022).

The `AEZ55_NNNN` and `EZ55_NNNNN` identifiers are not related by zero-padding or string transformation — they originate from independent gene-calling runs on different assemblies. No NCBI cross-reference exists between the two.

### Protein sequence-based bridging

To bridge the AEZ55 and EZ55 identifier namespaces, we obtained the original protein FASTA and GenBank flat files for the three draft contigs directly from the corresponding author (G. Hennon, personal communication). This FASTA contained 4,930 predicted protein sequences with `AEZ55_NNNN` headers, covering all CDS features in the draft annotation (the remaining 86 of the 5,016 gene features are tRNAs, rRNAs, and ncRNAs). We then matched these against the 4,096 RefSeq proteins from the current assembly (GCF_901457815.2) using a three-phase approach:

**Phase 1 — Exact sequence matching.** Each draft protein sequence was compared against all RefSeq sequences. An exact string match was accepted when the draft protein mapped to exactly one RefSeq accession, which was then chained to a locus tag via `gene_mapping.csv` (RefSeq protein accession → locus tag). This phase resolved 2,971 proteins (60.3%).

**Phase 2 — Subsequence matching.** For unmatched proteins, we tested whether one sequence was a contiguous subsequence of the other (capturing annotation boundary differences where the same CDS was predicted with different start or stop codons). A match was accepted when the shorter sequence was contained within the longer and covered ≥95% of the longer sequence's length. This phase resolved an additional 136 proteins (cumulative 63.0%).

**Phase 3 — Diamond blastp similarity search.** Remaining unmatched proteins were searched against the full RefSeq protein database using Diamond v2.1.9 (Buchfink et al., 2021) in `--very-sensitive` mode, requiring ≥80% sequence identity and ≥60% query coverage. Subject coverage was intentionally not filtered because the draft assembly's fragmented scaffolds introduced frameshifts that caused the draft annotation pipeline to split single NCBI genes into multiple shorter open reading frames. In these cases, each draft ORF aligns well across its own length (high query coverage) but covers only a fraction of the larger RefSeq protein (low subject coverage). This phase resolved an additional 946 proteins (cumulative 82.2%, 4,053 of 4,930).

A reciprocal best-hit (RBH) filter was not applied. RBH is standard for cross-species ortholog detection, where paralogous hits confound one-to-one assignment. In this within-organism context, both annotations derive from the same genome and the proteins are near-identical (≥95% identity for the vast majority of Phase 3 hits). The risk of a false paralog assignment at ≥80% identity with ≥60% query coverage is negligible for same-genome comparisons. Furthermore, the fragment deduplication step (below) already resolves the many-to-one direction that RBH would otherwise address.

**Fragment deduplication.** When multiple draft ORFs mapped to the same RefSeq locus tag (indicating a single gene split into fragments by the draft annotation), only the longest fragment was retained in the final mapping. The rationale is that RNA-seq read counts scale with transcript length; the longest fragment captured the most reads and therefore best represents the gene's expression level. Retaining all fragments would create duplicate expression edges with partial, potentially conflicting fold-change values. This deduplication removed 496 shorter fragments.

The final mapping was written as a two-column CSV (`aez55_id`, `locus_tag`) and added to the Hennon 2017 paperconfig as an `id_translation` entry with `id_type: old_locus_tag` (Tier 1), enabling the standard three-tier resolution pipeline to resolve Hennon's differential expression tables. Resolution rates were 345/422 (81.8%) for the coculture DE table and 98/125 (78.4%) for the axenic DE table. The 104 unresolved genes comprise 47 draft-specific ORFs with no detectable homolog in the current assembly (the draft predicted 834 more CDS than the finished genome) and 49 genes discarded by the fragment deduplication step, with the remainder falling below the Diamond identity or coverage thresholds. All unresolved rows are annotated as such in the resolved CSV and produce no edges in the knowledge graph, avoiding dangling references.

The mapping script (`scripts/map_img_to_ncbi_proteins.py`) is provided for reproducibility and can be adapted for other strains with similar cross-assembly challenges. It requires Diamond (Buchfink et al., 2021) for the similarity search phase.

### Barreto 2022 annotation

Barreto et al. (2022) provided a separate GTF annotation file (`EZ55.exon.fixed2.gtf`) based on the CABDXN contig accessions (ENA WGS accessions for GCA_901457815.2). This annotation uses `EZ55_NNNNN` gene IDs identical to the GCA locus tags, so no cross-assembly bridging is needed — the standard `locus_tag` pass resolves these directly at ~95%.

### Case study: *Alteromonas macleodii* MIT1002

MIT1002 presents a similar cross-assembly challenge. Two independent genome assemblies exist:

1. **GCF_001077695.1** (MIT/Chisholm Lab, 2015): A draft Illumina assembly of two contigs (NZ_JXRW01000001.1, NZ_JXRW01000002.1) with approximately 4,159 CDS annotated by NCBI PGAP using locus tag prefix `TK37_RS*`. This assembly was produced by the same group that performed the transcriptomic experiments in Biller et al. (2016, 2018).

2. **GCF_901457835.2** (ICBM, 2022): A complete PacBio genome consisting of a chromosome (LR881168.1, 4.7 Mb) and plasmid (LR881169.1, 203 kb), with 4,028 annotated genes. The GCA annotation uses `MIT1002_NNNNN` (5-digit) locus tags; NCBI PGAP assigned `ALT831_RS*` locus tags to the RefSeq version. This is the canonical assembly used in the knowledge graph.

The Biller et al. publications used a third, independent annotation layer: RAST (Rapid Annotations using Subsystems Technology) applied to the draft genome. RAST assigned `MIT1002_NNNN` (4-digit) identifiers with internal RAST `fig|226.6.peg.N` accessions. These 4-digit identifiers are unrelated to the 5-digit `MIT1002_NNNNN` locus tags from the current assembly — they originate from independent gene-calling runs on different assemblies and cannot be interconverted by zero-padding. Validation of the zero-padding hypothesis showed only 3.8% product description match between the padded pairs, confirming they map to different genes.

A supplementary conversion table (`MIT1002_systematicnames_conversiontable.csv`) provided by Biller et al. (2018) maps RAST fig| IDs to 4-digit `MIT1002_NNNN` genbank IDs and assembly coordinates on the JXRW draft contigs, but contains no cross-references to the current assembly's identifiers.

The proposed bridging strategy uses a two-step chain: (1) match RAST coordinate positions against NCBI PGAP gene positions in the old draft GFF (`GCF_001077695.1`) to obtain `TK37_RS*` locus tags, then (2) use the 3,892 shared WP_ RefSeq protein accessions between the draft and current assemblies to map `TK37_RS*` → `MIT1002_NNNNN`. Alternatively, Coe et al. (2024) supplementary data provides a direct `cds-TK37_RS*` → `MIT1002_NNNNN` mapping for 3,876 genes, which can serve as the bridge for step 2. Author contact has been initiated to obtain a direct RAST-to-current mapping or the RAST protein FASTA for sequence-based verification.

## Diagnostic Output and Quality Monitoring

After building `gene_id_mapping.json`, a diagnostic report (`gene_id_mapping_report.json`) was written for each strain. This report enumerates, for each identifier type, the fraction of identifiers that resolved uniquely (in `specific_lookup`), matched multiple genes (in `multi_lookup`), or produced Tier 1 conflicts. Reclassification warnings were emitted when an identifier type declared as Tier 1 exhibited a multi-match rate exceeding 10%, suggesting the type should be moved to Tier 2.

Regression testing used a snapshot-based approach: `Affects_expression_of` edge counts were recorded per publication before and after each pipeline change. Any publication showing a reduction in edge count triggered investigation of the corresponding `_resolved_report.txt` and re-examination of the paperconfig `id_columns` configuration.

## Phased Validation

The redesigned pipeline was validated strain by strain, beginning with *Prochlorococcus* MIT9312. For each strain, the pipeline proceeded as: (1) rebuild `gene_id_mapping.json` from all annotation sources and paperconfigs; (2) re-resolve all publication CSVs; (3) compare expression edge counts in the knowledge graph against the pre-rebuild snapshot; (4) investigate any publication showing reduced edge counts. This phased approach ensured that improvements in identifier resolution were verified before proceeding to the next strain.
