# Computational Tools for Enhancing Knowledge Graph-Based Exoenzyme Discovery in Marine Bacteria

## Overview

This section evaluates computational tools for improving knowledge graph (KG)-based identification of extracellular enzymes (exoenzymes) in marine bacteria, with emphasis on *Alteromonas* and *Prochlorococcus*. We assess each tool's community adoption, methodological robustness, practical deployment on Linux, and relevance to exoenzyme characterization. Tools are organized by the annotation layer they address: signal peptide prediction, subcellular localization, transmembrane topology, protease classification, protein domain architecture, and secretion system detection.

## 1. Signal Peptide Prediction

### 1.1 SignalP 6.0

**Purpose.** Predicts N-terminal signal peptides that target proteins for export across the cytoplasmic membrane. Signal peptides are the primary determinant distinguishing secreted exoenzymes from cytoplasmic housekeeping proteases.

**Method.** SignalP 6.0 (Teufel et al., 2022) employs a protein language model (fine-tuned ESM-1b transformer) with a conditional random field (CRF) decoder. It is the first tool to predict all five known signal peptide types in a single model: Sec/SPI (standard), Sec/SPII (lipoprotein), Tat/SPI, Tat/SPII, and Sec/SPIII (pilin-like).

**Adoption and robustness.** SignalP 6.0 has accumulated approximately 2,000 citations since its 2022 publication in *Nature Biotechnology* and is the de facto standard for signal peptide prediction. It is integrated into NCBI's Prokaryotic Genome Annotation Pipeline (PGAP), InterProScan (as an optional module), and numerous genome annotation workflows. Benchmark comparisons demonstrate improved accuracy over SignalP 5.0, DeepSig, and Phobius, particularly for lipoprotein (SPII) and Tat signal peptides. The protein language model backbone provides strong generalization to divergent sequences.

**Practical deployment.** Available for local installation on Linux under a free academic license from DTU Health Tech. Distributed as a pip-installable Python package (~4 GB with model weights). Requires Python 3.8--3.11 and PyTorch. A fast single-model mode processes a typical bacterial proteome (~4,000 proteins) in 20--60 minutes on CPU; GPU acceleration reduces this to 1--5 minutes. The tool requires specifying the organism group (gram-negative for all organisms in this study).

**Output.** Tab-delimited prediction results reporting the predicted signal peptide type (SP/LIPO/TAT/TATLIPO/PILIN/OTHER), probability scores for each type, and cleavage site position. GFF3 output and mature sequence FASTA are also produced.

**Relevance.** In a pilot analysis of *Alteromonas* sp. HOT1A3 (4,031 genes), SignalP 6.0 identified 508 proteins with standard signal peptides (SP), 179 lipoproteins (LIPO), 19 Tat-exported proteins, and 16 pilin-like signal peptides. Among annotated proteases, 14 were predicted as secreted, including S8 family serine peptidases with SP signals and M13/M3 family metallopeptidases with lipoprotein signals. This compares to only 140 genes (3.5%) with signal peptide annotation derived from UniProt in the current KG, demonstrating the substantial coverage improvement from direct prediction.

**Reference.** Teufel, F. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat. Biotechnol.* 40, 1023--1025 (2022).

### 1.2 SecretomeP 2.0

**Purpose.** Predicts non-classical (leaderless) protein secretion, i.e., export of proteins lacking canonical signal peptides. Relevant because some bacterial exoenzymes are secreted via non-classical pathways including outer membrane vesicles and specialized secretion systems.

**Method.** SecretomeP 2.0 (Bendtsen et al., 2005) uses a neural network trained on sequence-derived features (amino acid composition, predicted secondary structure, disordered regions). Proteins scoring above 0.5 are predicted as non-classically secreted.

**Adoption and robustness.** The original publication has over 1,500 citations and SecretomeP remains the most commonly cited dedicated tool for non-classical secretion prediction. However, the tool has significant limitations. The models have not been updated since 2005, the training set was small by modern standards, and the false positive rate is notably high. It does not distinguish between different non-classical secretion mechanisms (Type VI effectors, membrane vesicle cargo, holins, etc.). Independent benchmarks report a Matthews correlation coefficient of 0.40--0.50, indicating moderate discriminative power.

**Practical deployment.** Available as a web server at DTU Health Tech and as a standalone download for Linux under academic license. The standalone version has dependencies on SignalP and TMHMM for pre-screening. No updates to the software have been released since the original publication.

**Recommendation.** SecretomeP should be used only as one filter in a multi-tool pipeline (SignalP-negative, TMHMM-negative, PSORTb cytoplasmic/unknown, SecretomeP-positive) rather than as a standalone predictor. Its value is limited by the outdated models and high false positive rate. No robust modern alternative exists specifically for bacterial non-classical secretion prediction, though effector-specific tools (Bastion3, Bastion6) address subsets of this problem.

**Reference.** Bendtsen, J. D. et al. Feature-based prediction of non-classical and leaderless protein secretion. *Protein Eng. Des. Sel.* 18, 249--256 (2005).

## 2. Subcellular Localization Prediction

### 2.1 PSORTb 3.0

**Purpose.** Predicts the subcellular localization of bacterial proteins into five compartments (Gram-negative): cytoplasm, cytoplasmic membrane, periplasm, outer membrane, and extracellular. Complements signal peptide prediction by determining the final destination of exported proteins.

**Method.** PSORTb 3.0 (Yu et al., 2010) uses a multi-module approach combining SVM classifiers, motif-based analysis, profile HMM matching, and homology searches (SCL-BLAST). A Bayesian integration layer combines module outputs into a final prediction. The tool was deliberately designed for high precision: when confidence is low, proteins are classified as "Unknown" rather than receiving a forced prediction.

**Adoption and robustness.** With over 1,500 citations, PSORTb is considered the gold standard for bacterial subcellular localization prediction. Reported precision is approximately 97% across all localization categories, with recall of approximately 92% for Gram-negative bacteria. The conservative "Unknown" classification (assigned to 10--20% of proteins) contributes to this high precision. PSORTb is used in the MicroScope genome annotation platform and is recommended in numerous genome annotation guidelines.

**Practical deployment.** Available via Docker (`brinkmanlab/psortb_commandline`), which is the recommended installation method. Native Linux installation requires legacy BLAST (blastall, not BLAST+) and HMMER 2.x, making it cumbersome on modern systems. Docker usage is straightforward:

```
docker run --rm -v $(pwd):/data brinkmanlab/psortb_commandline:1.0.2 \
  -i /data/proteins.fasta -r /data/results -o terse -n
```

Processing a bacterial proteome takes minutes. The `-n` flag specifies Gram-negative organisms.

**Limitations.** PSORTb does not distinguish thylakoid membranes from cytoplasmic membranes, which is relevant for cyanobacteria such as *Prochlorococcus* that possess extensive thylakoid membrane systems. It does not predict dual-localized proteins. The underlying models have not been substantially updated since 2010.

**Alternatives.** DeepLocPro (2024, DTU Health Tech) is a recent deep learning-based alternative using protein language model embeddings (ESM). It was specifically designed for prokaryotic proteins and claims improved accuracy over PSORTb, particularly for underrepresented localization classes. As a newer tool, it has less community validation but represents the state of the art in methodology. BUSCA (Bologna Unified Subcellular Component Annotator) provides a meta-predictor approach combining multiple tools.

**Reference.** Yu, N. Y. et al. PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. *Bioinformatics* 26, 1608--1615 (2010).

## 3. Transmembrane Topology Prediction

### 3.1 DeepTMHMM

**Purpose.** Predicts transmembrane protein topology, distinguishing alpha-helical transmembrane proteins (inner membrane), beta-barrel transmembrane proteins (outer membrane), signal peptides, and soluble (globular) proteins.

**Method.** DeepTMHMM (Hallgren et al., 2022) is the successor to TMHMM 2.0, one of the most cited bioinformatics tools. It uses protein language model embeddings with a conditional random field (CRF) decoder, providing per-residue topology classification into inside (i), outside (o), transmembrane helix (M), beta-strand in membrane (B), and signal peptide (S) states.

**Key advance over TMHMM 2.0.** DeepTMHMM is the first tool from the DTU group to predict both alpha-helical and beta-barrel transmembrane proteins in a unified model. This is particularly valuable for Gram-negative bacteria where distinguishing inner membrane (alpha-helical) from outer membrane (beta-barrel) proteins is essential for exoenzyme localization analysis. TMHMM 2.0 could only predict alpha-helical transmembrane domains.

**Adoption and robustness.** Published in *Nucleic Acids Research* (2022), DeepTMHMM has been rapidly adopted as the recommended replacement for TMHMM 2.0. Alpha-helical transmembrane protein detection accuracy is approximately 98%, comparable to TMHMM 2.0. Beta-barrel detection accuracy is comparable to specialized tools (PRED-TMBB2, BOCTOPUS2). The tool also provides signal peptide prediction competitive with SignalP 5.0, though for dedicated signal peptide analysis SignalP 6.0 remains more comprehensive (predicting all five SP types).

**Practical deployment.** Available via the BioLib platform:

```
pip install biolib
biolib run DTU/DeepTMHMM --fasta input.fasta
```

The BioLib client handles containerized execution. An internet connection is required for the first run to pull the application container; subsequent runs use cached containers. GPU acceleration is supported. A bacterial proteome processes in 5--15 minutes.

**Relevance.** The current KG stores transmembrane region counts derived from UniProt (TMHMM 2.0-based), covering 3,419 of 35,409 genes (9.7%). DeepTMHMM would provide 100% coverage and add the critical beta-barrel distinction. Outer membrane beta-barrel proteins in *Alteromonas* (porins, outer membrane proteases such as M48 family) are directly relevant to extracellular enzymatic activity.

**Reference.** Hallgren, J. et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. *Nucleic Acids Res.* 50, W633--W640 (2022).

## 4. Protease Classification

### 4.1 MEROPS Database

**Purpose.** Provides systematic, hierarchical classification of proteolytic enzymes (peptidases) and their inhibitors. Essential for distinguishing degradative exoproteases from processing peptidases and for inferring substrate specificity.

**Classification system.** MEROPS (Rawlings et al., 2018) organizes peptidases into a three-level hierarchy:

- **Catalytic type** (9 types): serine (S), cysteine (C), aspartic (A), metallo (M), threonine (T), glutamic (G), asparagine (N), mixed (P), unknown (U).
- **Clan** (superfamily level): groups families sharing evolutionary origin based on 3D structure and catalytic mechanism.
- **Family** (sequence homology level): proteins within a family share statistically significant sequence similarity. Each family has a type peptidase.

A MEROPS classification provides: catalytic residue identities, substrate specificity preferences (Schechter-Berger P4--P4' positions), experimentally verified substrates, known inhibitors, metal ion requirements, and species distribution.

**Adoption and robustness.** MEROPS is the authoritative reference for protease classification, with publications collectively cited tens of thousands of times. The MEROPS family/clan nomenclature (e.g., "S8 family serine peptidase," "clan PA") is standard terminology used across the protease biology literature, InterPro, UniProt, KEGG, and NCBI CDD.

**Practical deployment.** MEROPS provides a downloadable sequence library (`merops_scan.lib`) for local BLAST or Diamond searches:

```
diamond makedb --in merops_scan.lib --db merops
diamond blastp --query proteins.faa --db merops --outfmt 6 --evalue 1e-10 --out merops_hits.tsv
```

The database is hosted at EMBL-EBI (https://www.ebi.ac.uk/merops/). Maintenance has slowed since the database was transferred from the Sanger Institute to EMBL-EBI (~2018), but the classification system is stable and bacterial protease coverage is comprehensive.

**Relevance.** The current KG identifies proteases through text matching of product annotations (e.g., "S8 family peptidase"), which relies on the quality of upstream annotation. Direct MEROPS BLAST would provide systematic family and clan assignments independent of annotation quality, and would enable inference of substrate specificity -- a key gap identified in our HOT1A3 exoenzyme analysis, where the KG encodes enzyme family but not cleavage site preferences.

**Reference.** Rawlings, N. D. et al. The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database. *Nucleic Acids Res.* 46, D624--D632 (2018).

## 5. Protein Domain and Family Annotation

### 5.1 InterProScan

**Purpose.** Provides comprehensive protein domain, family, and functional site annotation by integrating 14+ member databases into a single analysis pipeline. Complements orthology-based annotation (eggNOG-mapper) with precise domain architecture characterization.

**Method.** InterProScan 5 (Jones et al., 2014) runs individual member database algorithms (HMM profiles, patterns, fingerprints) against each query protein and consolidates results through InterPro cross-references. Each member database contributes a distinct annotation perspective.

**Member databases (selected).** Pfam (protein families, HMM-based), PANTHER (families/subfamilies with GO mappings), CDD (NCBI Conserved Domain Database), SMART (signaling/extracellular domains), SUPERFAMILY (SCOP structural domains), Gene3D (CATH structural domains), PRINTS (protein fingerprints), PROSITE (functional sites and domains), HAMAP (microbial protein families), TIGRFAMs/NCBIfam (microbial protein families), SFLD (structure-function linkage), MobiDB-lite (intrinsically disordered regions), and Coils (coiled-coil regions). Optional integration with SignalP, TMHMM, and Phobius is supported but requires separate licenses.

**Adoption and robustness.** InterProScan is the gold standard for integrated protein annotation, with over 10,000 citations. It underpins UniProt's functional annotation pipeline and is used by Ensembl, JGI IMG, and most major genome annotation platforms. The tool is actively maintained by EMBL-EBI with releases approximately every 2 months tracking InterPro database updates.

**Comparison with eggNOG-mapper.** The two tools are complementary rather than redundant. eggNOG-mapper provides orthology-based functional transfer (COG categories, KEGG pathways/modules, fast GO/EC assignment), while InterProScan provides direct domain detection with precise boundary delineation, multi-domain architecture characterization, and structural domain classification. For exoenzyme analysis, InterProScan's domain architecture output reveals domain organization (e.g., signal peptide + protease domain + PPC domain) that is not captured by orthology-based methods.

**Practical deployment.** Distributed as a standalone Java application for Linux (requires Java 11+). Member database files require 80--120 GB of disk space. An optional pre-calculated lookup service for UniProt proteins adds 100--150 GB but dramatically accelerates analysis of known proteins. Processing a bacterial proteome (~4,000 proteins) with all analyses takes 1--4 hours on a multi-core workstation; disabling the slowest analysis (PANTHER) reduces this to 15--45 minutes. For the 13 organisms in this study (~45,000 proteins total), a complete run would require 4--12 hours.

**Output.** TSV, GFF3, JSON, and XML formats. TSV output includes: protein accession, analysis source database, signature accession, signature description, start/stop positions, e-value, InterPro accession, and GO terms.

**Relevance.** InterProScan would add precise domain architectures to gene/protein nodes, enable structural domain classification (SCOP/CATH superfamilies), and provide additional curated GO term assignments (InterPro2GO). For exoenzyme analysis, this would distinguish multi-domain secreted proteases from single-domain cytoplasmic peptidases and enable structural classification of protease domains.

**Reference.** Jones, P. et al. InterProScan 5: genome-scale protein function classification. *Bioinformatics* 30, 1236--1240 (2014).

## 6. Secretion System Detection

### 6.1 MacSyFinder / TXSScan

**Purpose.** Detects complete bacterial secretion systems and related macromolecular complexes in genome sequences. Identifies the cellular machinery responsible for exporting exoenzymes to the extracellular space.

**Method.** MacSyFinder (Abby et al., 2014; N\'eron et al., 2023) uses a two-step detection approach: (1) individual protein components are identified using HMMER profile HMM searches against curated component models, and (2) complete systems are validated by checking whether a sufficient quorum of mandatory and accessory components co-occur within a defined genomic distance (co-localization constraint). This quorum-based logic dramatically reduces false positives compared to single-gene annotation approaches, as isolated HMM hits (which could represent paralogs or distant homologs) are not called as complete systems.

**Systems detected by TXSScan.** The TXSScan model package detects: Type I Secretion System (T1SS, ABC transporter-based), Type II Secretion System (T2SS/General Secretory Pathway), Type III Secretion System (T3SS, injectisome), Type IV Secretion System (T4SS, multiple subtypes), Type V Secretion System (T5SS, autotransporters, subtypes a--d), Type VI Secretion System (T6SS, phage-like contractile injection), Type IX Secretion System (T9SS, Bacteroidetes), Type IV Pilus (T4P), Tight Adherence pilus (Tad), Mannose-Sensitive Hemagglutinin pilus (MSH), competence pseudopilus (Com), and bacterial flagellum.

**Adoption and robustness.** The MacSyFinder paper (2014) has over 300 citations and TXSScan (2016) over 200 citations. The tool is integrated into the MicroScope genome annotation platform at Genoscope/LABGeM and is considered the standard reference tool for systematic secretion system detection. The quorum/co-localization approach has been validated against manually curated secretion system annotations across diverse bacterial genomes with high concordance.

**Practical deployment.** MacSyFinder v2 is pip-installable and available via bioconda:

```
pip install macsyfinder
macsydata install TXSScan
macsyfinder --db-type ordered_replicon --sequence-db proteins.faa --models TXSScan all
```

Requires HMMER (hmmsearch) on the system PATH. The `macsydata` companion tool manages model package installation. Actively maintained by Eduardo Rocha's group at Institut Pasteur.

**Limitations.** Fragmented genome assemblies can split secretion system gene clusters across contigs, causing false negatives when the co-localization constraint fails. Highly divergent systems in deeply-branching phyla may be missed if HMM profiles were trained on limited taxonomic diversity. TXSScan does not cover the Sec and Tat general export pathways (which are detected by SignalP).

**Relevance.** For *Alteromonas* species, which are known to possess T2SS and T6SS, MacSyFinder/TXSScan would identify the complete secretion machinery available for exoenzyme export. T2SS is the primary pathway for extracellular protease secretion in Gammaproteobacteria. T6SS may mediate contact-dependent delivery of effector proteins relevant to *Prochlorococcus*--*Alteromonas* coculture interactions already tracked in the KG. For *Prochlorococcus* strains with minimal genomes, the analysis would reveal which secretion systems have been retained or lost.

**References.** Abby, S. S. et al. MacSyFinder: a program to mine genomes for molecular system components with an application to CRISPR-Cas systems. *PLOS ONE* 9, e110726 (2014). N\'eron, B. et al. MacSyFinder v2: improved modelling, new search engine, and large-scale comparative genomics analysis of macromolecular systems. *Peer Community J.* 3, e28 (2023).

## 7. Current Knowledge Graph Gaps

Analysis of the deployed Neo4j knowledge graph revealed several annotation gaps relevant to exoenzyme characterization:

**Signal peptide coverage.** Only 1,872 of 26,582 Protein nodes (7.0%) and 720 of 35,409 Gene nodes (2.0%) carry signal peptide annotations, derived exclusively from UniProt. *Alteromonas* HOT1A3 has zero genes with signal peptide annotations in the graph despite SignalP 6.0 predicting 706 secreted proteins (508 SP + 179 LIPO + 19 TAT), indicating that the UniProt-derived annotations do not reach HOT1A3 gene nodes due to incomplete protein-gene linkage.

**Subcellular localization.** Zero of 26,582 Protein nodes have the `subcellular_location` property populated, despite the field being defined in the graph schema and downloaded from the UniProt API. This indicates a pipeline defect where the data is fetched but lost during graph construction.

**Transmembrane topology.** Transmembrane region annotations (from UniProt/TMHMM 2.0) are present for 3,419 Gene nodes (9.7%) and 5,825 Protein nodes (21.9%). These are limited to alpha-helical transmembrane domain counts and do not distinguish outer membrane beta-barrel proteins.

**Protease classification.** Proteases are identified by text matching of product annotations rather than systematic classification. No MEROPS family, clan, or substrate specificity data is stored. Several S8 family peptidases in HOT1A3 lack KEGG KO assignments, preventing pathway-level inference.

**Secretion systems.** No secretion system annotations are present in the graph. The presence and completeness of T2SS, T6SS, and other export pathways in each organism is unknown.

## 8. Recommended Integration Strategy

Based on the assessments above, we propose a tiered integration strategy:

**Tier 1 (high impact, proven tools, straightforward automation):**
1. **SignalP 6.0** for all strains -- predicts signal peptide type and cleavage site. Store as `signalp_prediction` (SP/LIPO/TAT/TATLIPO/PILIN/OTHER), `signalp_score`, and `signalp_cs_position` on Gene nodes.
2. **PSORTb 3.0** for all strains -- predicts final subcellular localization. Store as `psortb_localization` and `psortb_score` on Gene nodes.
3. **MEROPS BLAST** for all strains -- systematic protease family/clan classification. Store as `merops_family`, `merops_clan`, and `merops_catalytic_type` on Gene nodes.

**Tier 2 (high value, heavier infrastructure):**
4. **InterProScan** for all strains -- comprehensive domain architecture. Store domain hits as properties or as dedicated Domain nodes with edges to Gene/Protein nodes.
5. **MacSyFinder/TXSScan** for all strains -- secretion system detection. Store as SecretionSystem nodes with component edges to Gene nodes.

**Tier 3 (incremental improvement):**
6. **DeepTMHMM** for all strains -- beta-barrel prediction and improved transmembrane topology. Store as `deeptmhmm_topology` (TM/Beta/SP+TM/Globular) on Gene nodes.
7. **Fix existing pipeline bugs** -- resolve the `subcellular_location` data loss and improve UniProt-to-gene linkage for *Alteromonas* strains.

This strategy would transform the KG from a gene-expression-centric resource into one capable of systematic exoenzyme discovery, enabling queries such as: "Find all genes encoding secreted serine proteases (SignalP SP, PSORTb extracellular, MEROPS clan SB) in organisms possessing a complete T2SS, with differential expression under coculture conditions."
