# Synechococcus + New Paper Integration Design

**Date:** 2026-04-03
**Scope:** Integrate 9 new Synechococcus coculture papers, update Prochlorococcus paperconfigs for newly added data, add ~12 new organisms to the KG (both cyanobacteria and heterotrophic coculture partners), and add `fold_change_type` support to the omics adapter.

## Decisions

- **Both sides of cocultures:** Load expression data for both the cyanobacterium and the heterotrophic partner in each coculture experiment
- **Pipeline structure:** Keep papers under `data/Synechococcus/papers_and_supp/`; extend pipeline to read multiple paperconfig list files
- **Fold-change conversion:** Add `fold_change_type` field to paperconfig (`log2` default, `linear` converts at load time); adapter handles conversion
- **Clustering:** Create `gene_clusters` paperconfig entries now; skip extraction JSONs (adapter handles missing JSON gracefully). Extraction is a separate future task
- **Genome CSV:** Add organisms with DE data and known accessions to `cyanobacteria_genomes.csv`. Organisms that only appear as coculture partners (no DE data, no accession) go into `treatment_organisms.csv` only
- **Cyanorak:** Use downloaded Cyanorak organism tables to match marine Synechococcus strains. Freshwater/thermophilic strains (PCC 7002, PCC 7942, T. elongatus) won't have Cyanorak entries
- **eggNOG:** Run as a separate task later (slow, hours per strain). Gene annotations will initially lack eggNOG fields
- **Iterative paperconfig review:** Stop after each paper's paperconfig for user review before proceeding to the next

## Phase 1: Prochlorococcus Paperconfig Updates

Create/update paperconfigs for papers with strains already in the KG (or being added).

**Workflow:** For each paper:
1. Explore the data (read CSVs, legends, paper if available)
2. Create or update the paperconfig (or decide to skip)
3. Create a README.md in the paper directory summarizing: what data files exist, what analyses are possible, what decision was made and why
4. Present both the paperconfig and README for user approval before moving to the next paper

### Complete Prochlorococcus paper inventory

Every paper directory is listed below. Papers touched by the recent commit are marked with *.

#### Papers needing new paperconfigs

| Paper | Data type | Strains | Action |
|---|---|---|---|
| * Capovilla 2023 | DE (RNA-seq), significant-only tables | MIT9303, MIT9313 | New paperconfig. MIT9303 not in KG -- add genome (see Phase 3). MIT9313 already in KG |
| * zinser 2009 | Time-series expression (50 timepoints) + diel clusters | MED4 | New paperconfig: `gene_clusters` entry for diel clusters. Expression data is raw values (not FC) -- skip for DE edges |
| * alonso 2023 | Soft clustering membership (5 clusters, probability scores) | MED4 | New paperconfig: `gene_clusters` entry only (no DE data) |
| * Wang 2014 | Gene classification (VEG/HEG/MEG) + operon predictions | MED4 | Skip -- not DE or clustering data suitable for KG |
| * Waldbauer 2012 | Supplementary figures + PDF tables only | MED4 | Skip -- no machine-readable expression data |

#### Papers with existing paperconfigs that may need updates

| Paper | What changed in commit | Action |
|---|---|---|
| * Biller 2018 | Added tables S4A/B (periodicity Y/N), S5 (darkness response) | Review -- periodicity markers not standard DE or clusters. Skip new tables for now |
| * coe 2024 | Added Supplemental Table 3 (cluster assignments) | Add `gene_clusters` entry to existing paperconfig |
| * lindell 2007 | Updated supp table 3 | Check if existing paperconfig references this table; update if needed |
| * Capovilla 2023 | New MIT9303 + MIT9313 DE CSVs | Covered above (new paperconfig) |
| * Hennon 2017 | Added two docx files (supp info) | No CSV change -- no paperconfig update needed |
| * tetu 2019 | Added two PDF files (supp info) | No CSV change -- no paperconfig update needed |

#### Papers with existing paperconfigs, no changes needed

| Paper | Status |
|---|---|
| Aharonovich 2016 | Has paperconfig, listed in paperconfig_files.txt |
| Al-Hosani 2015 | Has paperconfig, listed |
| Anjur 2025 | Has paperconfig, listed |
| bagby 2015 | Has paperconfig, listed |
| barreto 2022 | Has paperconfig, listed |
| biller 2016 | Has paperconfig, listed |
| Fang 2019 | Has paperconfig, listed |
| he 2022 | Has paperconfig, listed |
| Lin 2015 | Has paperconfig, listed |
| martiny 2006 | Has paperconfig, listed |
| MIT9313_resources | Has paperconfig, listed |
| Read 2017 | Has paperconfig, listed |
| steglich 2006 | Has paperconfig, listed |
| thompson 2011 | Has paperconfig, listed |
| Thompson 2016 | Has paperconfig, listed |
| tolonen 2006 | Has paperconfig, listed |
| Weissberg_2025 | Has paperconfig, listed |
| ziegler 2025 | Has paperconfig, listed |

#### Papers that are not actionable

| Paper | Reason |
|---|---|
| QPCR | Empty directory, no data |
| Labban 2022 | Commented out -- waiting for author GFF (contacted 2026-03-03) |

### zinser 2009 detail

Table_S1.csv has a unique structure requiring special handling:
- **Expression columns:** 50 timepoints of normalized expression values (not fold-change). These represent absolute expression across a diel cycle
- **Fourier analysis:** Score + FDR for diel periodicity detection
- **Clustering:** Cluster ID + membership score for diel response clusters
- **Action:** Create `gene_clusters` entry for the cluster assignments. The expression time-series data is raw values, not differential -- skip for DE edges unless we compute fold-changes from the time series (out of scope for this project)

### Capovilla 2023 detail

The CSVs list only significant DE genes (upregulated) with log2 fold-change but no p-values. The tables have separate sections for up-regulated and down-regulated genes within the same CSV (need to check exact structure). Gene IDs use `ncbi_cds_locus_tag` and `ncbi_gene_locus_tag` columns. MIT9303 needs a new genome entry.

### alonso 2023 detail

TABLE S5 has soft cluster membership with probability scores per cluster (A-E). Gene IDs are BV-BRC (PATRIC) format (`fig|167546.4.peg.NNN`) plus RefSeq locus tags. The RefSeq locus tags should map directly. `score_col` can capture the max cluster probability.

## Phase 2: Synechococcus Paperconfig Exploration + Creation

### Infrastructure changes

1. **Create `data/Synechococcus/papers_and_supp/paperconfig_files.txt`** -- list of paperconfig paths for Synechococcus papers
2. **Extend `create_knowledge_graph.py`** to read paperconfig list files from both `data/Prochlorococcus/.../paperconfig_files.txt` and `data/Synechococcus/.../paperconfig_files.txt`. Minimal change: read a list of list-file paths (or hardcode both paths)
3. **Add `fold_change_type` field to omics adapter:**
   - New optional field in `statistical_analyses`: `fold_change_type: linear` (default: `log2`)
   - When `linear`: adapter converts via `math.log2(fc)` before storing as `log2_fold_change`
   - Handle edge cases: fc <= 0 (skip row with warning), fc = 1.0 maps to log2FC = 0
   - Update paperconfig validation to accept the new field
4. **Add `fold_change_type` to paperconfig validation vocabulary** in `validate_paperconfig.py`

### Complete Synechococcus paper inventory

**Workflow:** Same as Phase 1 -- for each paper:
1. Explore the data (read CSVs, legends, paper PDF if available)
2. Create a paperconfig (or decide to skip)
3. Create a README.md summarizing data files, analyses, and decisions
4. Present both for user approval before proceeding

All 9 papers are new (no existing paperconfigs). Each is listed below with data assessment.

#### Beliaev 2014 -- Synechococcus PCC 7002 + Shewanella W3-18-1

- **8 CSV tables** from a single XLS, two organisms alternating:
  - Tables s1, s3, s5, s7: Synechococcus PCC 7002 genes
  - Tables s2, s4, s6, s8: Shewanella W3-18-1 genes
- **Columns:** Locus tag, Gene abbreviation, RPKM values (axenic, +Shew lactate, +Shew HCO3-), Fold change columns, annotations
- **Fold change type:** Linear (values like 0.46, 3.32, 16.667)
- **P-values:** None
- **Gene IDs:** `SYNPCC7002_A####` (Synechococcus), `SputW3181_####` (Shewanella) -- NCBI locus tags
- **Experiments:** Multiple comparisons per organism: lactate vs axenic, HCO3- vs axenic, HCO3- vs lactate. Tables s1/s2 = full gene list for one comparison set, s3/s4 = DE genes filtered differently, s5/s6 = another filter, s7/s8 = another. Need to read paper to understand exact table scoping
- **Organisms needed:** Synechococcus sp. PCC 7002, Shewanella putrefaciens W3-18-1

#### Bernstein 2017 -- T. elongatus BP-1 + Meiothermus ruber

- **Data set s1:** Raw transcript abundances (not fold-change) for T. elongatus across 6 conditions (3 light levels x axenic/coculture, plus 6 oxygen conditions). Locus tags: `Tlr####`/`Tll####` format
- **Data set s2 (light):** Clustered average expression per light condition + cluster assignments. Gene IDs: RAST (`fig.6666666.52381.peg.##`) and NCBI (`SY28_RS#####`). Columns: ClustID_light, Clust name_light
- **Data set s2 (oxygen):** Same structure, clustered by oxygen condition. ClustID_ox, Clust name_ox
- **Action:** Create `gene_clusters` entries for light and oxygen clusterings. Skip DE edges (raw values, not fold-change). Coculture partner: Meiothermus ruber
- **Organisms needed:** Thermosynechococcus elongatus BP-1 (note: some references call this strain BP-1, NCBI may list as PCC 6803 variant -- verify). Meiothermus ruber DSM 1279

#### Kratzl 2024 -- S. elongatus PCC 7942 + P. putida KT2440

- **4 CSV tables** (2 RNA-seq + 2 proteomics):
  - Sheet 1: coculture vs P. putida (RNA-seq) -- `gene-M8001_RS#####` IDs, has log2_foldChange + adjusted p_value
  - Sheet 6: coculture vs S. elongatus (RNA-seq) -- `gene-H6G84_RS#####` IDs, has log2_foldChange + adjusted p_value
  - Sheet 11: coculture vs P. putida (proteomics) -- protein names + UniProt IDs, log2_FC + adjusted p_value
  - Sheet 13: coculture vs S. elongatus (proteomics) -- protein names + UniProt IDs, log2-fold change + adjusted p_value
- **Fold change type:** log2 (all tables)
- **Gene IDs:** NCBI locus tag format with `gene-` prefix (RNA-seq), UniProt accessions (proteomics). The `gene-` prefix will need stripping
- **Organisms needed:** Synechococcus elongatus PCC 7942, Pseudomonas putida KT2440
- **Note:** Proteomics tables use protein-level IDs. Will need protein_id or uniprot_accession resolution to gene locus tags

#### Ma 2022 -- S. elongatus PCC 7942 + E. coli

- **3 CSV tables:**
  - S1: Up-regulated transcripts -- Gene ID (`M744_#####`), Description, log2FoldChange. No p-value
  - S2: Down-regulated transcripts -- same structure
  - S4: Differentially expressed proteins -- Protein ID (`M744_#####`), Description, Mean Ratio (linear, not log2)
- **Gene IDs:** `M744_#####` -- need to identify which assembly this corresponds to. Likely S. elongatus UTEX 2973 or PCC 7942 variant
- **Fold change type:** log2 (RNA-seq S1/S2), linear ratio (proteomics S4)
- **Organisms needed:** S. elongatus (verify exact strain from paper), E. coli (strain TBD from paper)
- **Note:** Tables S1 and S2 are split by direction -- will need to handle as two analyses pointing to different files, or merge into one CSV

#### Oleza 2015 -- Synechococcus WH7803 + R. pomeroyi DSS-3

- **Exoproteome (proteomics):**
  - Table s2a: Syn WH7803 proteins -- RefSeq protein IDs (`YP_001224###`), fold change axenic vs co-culture, p-value
  - Table s2b: R. pomeroyi DSS-3 proteins -- NCBI protein IDs, no fold-change (abundance only: NSAF, SC)
- **Fold change type:** Linear (s2a). s2b has no fold-change -- skip for DE edges
- **Gene IDs:** RefSeq protein accessions (YP_ format) for Syn, RefSeq protein accessions for Rpom
- **Organisms needed:** Synechococcus sp. WH7803 (Cyanorak: Syn_WH7803, RefSeq NC_009481.1), Ruegeria pomeroyi DSS-3

#### Oleza 2017 -- Synechococcus WH7803 + R. pomeroyi DSS-3

- **5 CSV tables:**
  - Table s2a: Syn WH7803 proteins (early coculture) -- fold change axenic vs co-culture + p-value. Linear FC
  - Table s2b: R. pomeroyi proteins (early) -- abundance only, no FC. Skip for DE
  - Table s5a: Syn WH7803 proteins (late coculture) -- fold change + p-value (ANOVA). Linear FC
  - Table s5b: R. pomeroyi proteins (late) -- abundance only. Skip for DE
  - Table s7: R. pomeroyi long-term proteins -- has fold change + q-value + significance flag. Linear FC
- **Organisms needed:** Same as Oleza 2015 (WH7803, R. pomeroyi DSS-3)

#### Tal 2009 -- Synechococcus WH8102 + Vibrio parahaemolyticus

- **2 CSV tables:**
  - Table s1: Upregulated genes -- Gene Name (descriptive), Gene ID (`SYNW####`), Score, Log2 fold change
  - Table s2: Downregulated genes -- same structure
- **Fold change type:** log2
- **P-values:** None (has "Score" which may be a statistical metric but not a p-value)
- **Gene IDs:** `SYNW####` -- standard WH8102 locus tags, already in KG
- **Organisms needed:** WH8102 already in KG. Vibrio parahaemolyticus as treatment organism (strain TBD from paper)

#### kaur 2018 -- Synechococcus WH7803 + R. pomeroyi DSS-3

- **4 CSV tables (exoproteome):**
  - Table s4a: Syn WH7803 in seawater (SW) -- RefSeq protein IDs (`YP_001223###`), fold change (multiple timepoints: TP3 vs 1, TP7 vs 1, etc.), q-values
  - Table s4b: Syn WH7803 in artificial seawater (ASW) -- same structure
  - Table s5a: R. pomeroyi in SW -- NCBI protein IDs (`AAV93###`), fold change, q-values
  - Table s5b: R. pomeroyi in ASW -- same structure
- **Fold change type:** Needs verification from paper -- values include negative numbers (-1.6) which could indicate log2. Check paper methods section during paperconfig creation
- **Gene IDs:** RefSeq protein accessions. Will need protein_id resolution
- **Organisms needed:** Same as Oleza papers (WH7803, R. pomeroyi DSS-3)

#### Zhang 2021 -- Synechococcus WH7803

- Doc supplement only -- no machine-readable expression data. **Skip.**

### New organisms inventory (produced by this phase)

#### Cyanobacteria (add to cyanobacteria_genomes.csv)

| Organism | Cyanorak name | RefSeq Acc | Papers |
|---|---|---|---|
| Synechococcus sp. WH7803 | Syn_WH7803 | NC_009481.1 | kaur 2018, Oleza 2015, Oleza 2017 |
| Synechococcus sp. PCC 7002 | not in Cyanorak | TBD (look up on NCBI) | Beliaev 2014 |
| Thermosynechococcus elongatus BP-1 | not in Cyanorak | TBD | Bernstein 2017 |
| Synechococcus elongatus PCC 7942 | not in Cyanorak | TBD | Kratzl 2024, Ma 2022 |
| Prochlorococcus MIT9303 | Pro_MIT9303 | CP000554 (GenBank, no RefSeq) | Capovilla 2023 |

**Strain verification needed during Phase 2 paperconfig creation:**
- Verify exact NCBI taxonomy and assembly for each organism by checking the paper methods section
- Kratzl 2024 gene IDs (`gene-M8001_RS#####` for P. putida, `gene-H6G84_RS#####` for S. elongatus) and Ma 2022 gene IDs (`M744_#####`) suggest different assemblies or locus tag prefixes for what may be the same species -- confirm whether these are the same strain (PCC 7942) or variants (e.g., UTEX 2973)
- T. elongatus BP-1 (Bernstein 2017): verify current NCBI classification -- genus was reclassified; confirm correct accession
- WH7803 Cyanorak RefSeq is NC_009481.1 (chromosome-level) -- need the GCF assembly accession for `cyanobacteria_genomes.csv`

#### Heterotrophs with DE data (add to cyanobacteria_genomes.csv -- full genome)

| Organism | RefSeq Acc | Papers | DE data |
|---|---|---|---|
| Shewanella putrefaciens W3-18-1 | TBD | Beliaev 2014 | Yes (tables s2, s4, s6, s8) |
| Pseudomonas putida KT2440 | TBD | Kratzl 2024 | Yes (RNA-seq + proteomics) |
| Ruegeria pomeroyi DSS-3 | TBD | kaur 2018, Oleza 2015, Oleza 2017 | Yes (exoproteome + long-term) |

#### Heterotrophs without DE data (treatment_organisms.csv only -- no genome)

| Organism | Papers | Why treatment-only |
|---|---|---|
| Meiothermus ruber DSM 1279 | Bernstein 2017 | Only clustering data, no DE for heterotroph side |
| E. coli (strain TBD) | Ma 2022 | DE data is for S. elongatus side only |
| Vibrio parahaemolyticus (strain TBD) | Tal 2009 | DE data is for Synechococcus side only |

**Note:** Organisms with DE data go into `cyanobacteria_genomes.csv` for full genome loading (genes, proteins, annotations). Their organism node comes from the genome pipeline. Organisms without DE data go into `treatment_organisms.csv` only -- they get an OrganismTaxon node for `Tests_coculture_with` edges but no genes/proteins.

## Phase 3: Genome Setup + Data Download

1. Download Cyanorak organism list (already done -- `data/Cyanorak Organism Table synechococcus.csv` and `prochlorococcus.csv`)
2. For each new organism: find RefSeq assembly accession (paper first, then NCBI). Use GenBank accession if no RefSeq exists (e.g., MIT9303)
3. Add rows to `cyanobacteria_genomes.csv`
4. Add heterotroph coculture partners to `treatment_organisms.csv`
5. Run `prepare_data.sh --steps 0` for new strains (NCBI genome + Cyanorak + UniProt download)
6. Run `prepare_data.sh --steps 1 2` (build protein annotations, build gene annotations -- initially without eggNOG)
7. eggNOG: run separately as a later task (`/eggnog-run` for all new strains, then re-run step 2 to merge)
8. Expect iteration: some papers may reference assemblies that don't match current RefSeq, requiring `id_translation` or `annotation_gff` entries

## Phase 4: Gene ID Resolution

1. Run `prepare_data.sh --steps 3 4` for all strains (rebuild gene_id_mapping.json and resolve paper IDs -- new paperconfigs may add id_columns or id_translation entries that affect existing strains too)
2. Run `/check-gene-ids` for all papers with expression data (not just new papers -- existing papers may be affected by mapping changes)
3. Fix mismatches as needed:
   - Add `id_translation` entries for papers using non-standard IDs
   - Add `annotation_gff` entries where GFF bridges are needed
   - For protein-level IDs (kaur 2018, Oleza papers, Kratzl proteomics): may need `id_columns` with `protein_id` type
4. Re-run steps 3-4 after fixes
5. Iterate until match rates are acceptable

### Known ID challenges

| Paper | Gene ID format | Expected difficulty |
|---|---|---|
| Beliaev 2014 | NCBI locus tags (SYNPCC7002_A####, SputW3181_####) | Low -- standard format |
| Bernstein 2017 | RAST (fig.6666666.52381.peg.##) + NCBI (SY28_RS#####) | Medium -- RAST IDs need id_translation, NCBI RS format should resolve |
| Kratzl 2024 RNA | gene-M8001_RS##### / gene-H6G84_RS##### | Low -- strip `gene-` prefix |
| Kratzl 2024 protein | UniProt accessions | Medium -- protein-level resolution |
| Ma 2022 | M744_##### | Medium -- need to identify assembly |
| Oleza 2015/2017 | RefSeq protein accessions (YP_######) | Medium -- protein-level resolution |
| kaur 2018 | RefSeq protein accessions (YP_######, AAV#####) | Medium -- protein-level resolution |
| Tal 2009 | SYNW#### | Low -- standard WH8102 locus tags already in KG |
| Capovilla 2023 | ncbi_cds_locus_tag, ncbi_gene_locus_tag | Low -- standard NCBI format |

## Phase 5: Build and Validate

1. Take edge snapshot with `/omics-edge-snapshot`
2. Build KG: `uv run python create_knowledge_graph.py`
3. Docker deploy: `docker compose up -d`
4. Run post-import scripts
5. Compare edge snapshot -- verify new edges, no regressions on existing papers
6. Run KG validity tests: `pytest -m kg -v`
7. Update test expectations:
   - Organism count (currently 15 -- will increase by ~12)
   - Gene count (currently ~35K)
   - Expression edge count (currently ~188K)
   - Experiment count (currently 76)

## Cluster Extraction (Separate Future Task)

Papers with clustering data created now, extraction later:

| Paper | Cluster type | Organism | Clusters |
|---|---|---|---|
| zinser 2009 | Diel response clusters | MED4 | Multiple (from Cluster column) |
| alonso 2023 | Soft clusters (A-E) | MED4 | 5 clusters with probability scores |
| coe 2024 | Dark-tolerant clusters | MED4 | Multiple (dark-tolerant_cluster column) |
| Bernstein 2017 | Light response + oxygen response | T. elongatus BP-1 | Multiple per condition |

The cluster adapter handles missing extraction JSONs gracefully -- nodes are created with minimal properties (no descriptions, no temporal metadata). Extraction JSONs will be generated in a separate task using the two-stage extraction design.

## Code Changes Summary

| Change | Phase | Files affected |
|---|---|---|
| New Prochlorococcus paperconfigs | 1 | 3-4 new paperconfig.yaml files |
| Synechococcus paperconfig_files.txt | 2 | New file |
| Multiple paperconfig list file support | 2 | `create_knowledge_graph.py` |
| `fold_change_type` adapter support | 2 | `omics_adapter.py` |
| `fold_change_type` validation | 2 | `validate_paperconfig.py` |
| New Synechococcus paperconfigs | 2 | 7-8 new paperconfig.yaml files |
| New genome rows | 3 | `cyanobacteria_genomes.csv` |
| New treatment organisms | 3 | `treatment_organisms.csv` |
| id_translation entries | 4 | Various paperconfig.yaml files |
| KG test expectations | 5 | `tests/kg_validity/` |

## Out of Scope

- Renaming `cyanobacteria_genomes.csv`
- Cluster extraction JSONs (separate task)
- Bernstein 2017 data set s1 DE analysis (raw expression values, not fold-change)
- Biller 2018 S4A/B periodicity data, S5 darkness categories
- Wang 2014 gene classification / operon data
- Waldbauer 2012 (PDF tables)
- Zhang 2021 (doc supplement only)
