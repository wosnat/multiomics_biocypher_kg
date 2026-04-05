# Kratzl 2024

**Title:** Pseudomonas putida as saviour for troubled Synechococcus elongatus in a synthetic co-culture -- interaction studies based on a multi-OMICs approach

**DOI:** https://doi.org/10.1038/s42003-024-06098-5

**Journal:** Communications Biology (2024) 7:452

**Authors:** Franziska Kratzl, Marlene Urban, Jagroop Pandhal, Mengxun Shi, Chen Meng, Karin Kleigrewe, Andreas Kremling, Katharina Pflueger-Grau

## Organisms

- **Synechococcus elongatus PCC 7942** (NCBI Taxid: 1140) -- cscB strain (sucrose-secreting)
  - Assembly: GCF_000012525.1
  - Gene ID prefix in RNA-seq data: `H6G84_RS#####` (with `gene-` prefix that needs stripping)
  - Locus tags in sheet 6: `Synpcc7942_####` style (old locus tags)
- **Pseudomonas putida KT2440** (NCBI Taxid: 160488) -- EM178 att::Tn7 cscRABY strain
  - Assembly: GCF_000007565.2
  - Gene ID prefix in RNA-seq data: `M8001_RS#####` (with `gene-` prefix that needs stripping)
  - Locus tags in sheet 1: `PP_####` style (old locus tags)

Both organisms are **new to the KG** and will need genome data download and annotation before integration.

## Experimental Design

Synthetic co-culture of S. elongatus cscB (phototrophic, sucrose-secreting) and P. putida cscRABY (heterotrophic, sucrose-consuming). S. elongatus cscB was induced with 0.1 mM IPTG. The co-culture was grown in a 9-fold parallel photobioreactor (CellDEG HDC 9.100) with BG11+ medium supplemented with 150 mM NaCl. Exponential light profile: 120 umol photons m-2 s-2 constant for 24h, then exponential rising with td=52h. Samples for multi-OMICs taken at ~60h.

The co-culture was compared to axenic cultures of each organism under matched conditions.

**Key finding:** P. putida cscRABY promotes S. elongatus cscB growth (up to 80% increase), while S. elongatus has a neutral effect on P. putida.

## Data Files

| File | Sheet | Type | Organism | Rows | Description |
|------|-------|------|----------|------|-------------|
| `data set s1 1_...csv` | Sheet 1 | RNA-seq | P. putida | 5218 | All genes, coculture vs P. putida axenic |
| `data set s1 6_...csv` | Sheet 6 | RNA-seq | S. elongatus | 2679 | All genes, coculture vs S. elongatus axenic |
| `data set s1 11_...csv` | Sheet 11 | Proteomics | P. putida | 232 | All proteins, coculture vs P. putida axenic |
| `data set s1 13_...csv` | Sheet 13 | Proteomics | S. elongatus | 567 | All proteins, coculture vs S. elongatus axenic |

### Column Details

**RNA-seq tables (sheets 1 and 6):**
- `gene ID`: NCBI locus tag with `gene-` prefix (e.g., `gene-M8001_RS21020`)
- `log2_foldChange`: log2 fold change
- `p_value`: raw p-value
- `adjusted p_value`: FDR-corrected p-value
- `Protein ID` / `rotein ID` (typo in sheet 1): RefSeq WP_ accession
- `locus tag` / `Locus Tag`: old-style locus tags (PP_#### or Synpcc7942_####)
- `Annotation`: gene product description

**Proteomics table (sheet 11, P. putida):**
- `Gene Name`: gene name or PP_#### locus tag (sometimes combined, e.g., "lptD PP_0404")
- `Protein ID`: UniProt accession
- `log2_FC`: log2 fold change
- `adjusted p_value`: FDR-corrected p-value (no raw p-value column)
- `Protein names`: protein description

**Proteomics table (sheet 13, S. elongatus):**
- `Gene Name`: gene name or Synpcc7942_#### locus tag
- `Protein ID`: UniProt accession
- `log2-fold change`: log2 fold change (different column name from sheet 11)
- `p_value`: raw p-value
- `adjusted p_value`: FDR-corrected p-value
- `Protein.names`: protein description

## Statistical Thresholds

- DEG threshold: |log2(FC)| > 1.0 AND adjusted p-value (FDR) < 0.05
- DAP threshold: |log2(FC)| > 1.0 AND adjusted p-value (FDR) < 0.05
- Sheets 1, 6, 11, 13 contain ALL detected genes/proteins (not filtered)
- Sheets 2-5, 7-10, 12, 14 contain filtered DEGs/DAPs (not used here)

## Integration Notes

1. **Both organisms are new to the KG.** Need to add genome entries to `cyanobacteria_genomes.csv` and run `prepare_data.sh` before integration.
2. **RNA-seq gene IDs have `gene-` prefix** that needs to be stripped to get standard NCBI locus tags (e.g., `gene-M8001_RS21020` -> `M8001_RS21020`). Handle in Phase 4 (id resolution).
3. **Proteomics tables use UniProt accessions** -- will need protein-level resolution to map to gene nodes.
4. **Sheet 11 is missing p_value column** -- only has `adjusted p_value`.
5. **Sheet 1 has typo** in column header: `rotein ID` instead of `Protein ID`.
6. **Sheet 13 uses different column name** for fold change: `log2-fold change` (with hyphen) instead of `log2_FC`.
7. **P. putida KT2440** strain used is actually EM178 att::Tn7 cscRABY (engineered for sucrose transport/metabolism).
8. **S. elongatus PCC 7942** strain used is cscB (engineered for sucrose secretion).
