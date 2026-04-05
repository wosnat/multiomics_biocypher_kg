# Oleza 2015 — Synechococcus WH7803 + R. pomeroyi DSS-3 exoproteome

**Citation:** Christie-Oleza JA, Armengaud J, Bertrand P, Laez B, Boez B, Mariadassou S, Pelletier E (2015) Functional assessments of microbial interactions in the ocean using exoproteomics. *Proteomics* 15:3315-3325.

## Organism

- *Synechococcus* sp. WH7803 (cyanobacterium)
- *Ruegeria pomeroyi* DSS-3 (heterotrophic coculture partner)

## Data files

| File | Description |
|------|-------------|
| `pmic8102-sup-0001-tables1.xlsx` | Table S1 |
| `pmic8102-sup-0002-tables2.xlsx` | Table S2: Synechococcus WH7803 exoproteins (NSAF abundance) |
| `pmic8102-sup-0003-tables3.xlsx` | Table S3 |
| `pmic8102-sup-0004-tables4.xlsx` | Table S4 |
| `pmic8102-sup-0005-tables5.xlsx` | Table S5 |
| `Synechococcus coculture exoproteome Proteomics 2015.pdf` | Paper PDF |

## KG Integration Status

**Skipped** — no paperconfig created.

**Reason:** The proteomics data reports protein abundance (NSAF, spectral counts) rather than fold-change values. There are no DE comparisons with fold-change or p-values suitable for `Changes_expression_of` edges. The KG requires log2 fold-change for expression edges.

## Notes

- All supplementary data is in xlsx format only (no CSVs extracted)
- Table S2 has Synechococcus exoproteins with fold-change axenic vs co-culture and p-value for some entries, but the primary quantification is NSAF abundance
- Table S2b (R. pomeroyi proteins) has abundance only, no fold-change
- The companion paper Oleza 2017 extends this work with proper fold-change comparisons and is integrated separately
