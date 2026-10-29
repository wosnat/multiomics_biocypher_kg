# `_modified.csv` builders

One script per omics paper that needs a derived view of its supplementary CSVs (sign-corrected log2FC, merged up/down halves, normalized column names, extracted locus-tag columns, etc.) before paperconfig ingestion.

## Convention

- One script per paper, named `build_<paperslug>_modified_csv.py`.
- Read the original supp CSV from `data/<Organism>/papers_and_supp/<paper>/`.
- **Never edit the original** — write `<original_name>_modified.csv` in the same directory.
- The paperconfig's `filename:` points at the `*_modified.csv`. Add a YAML comment near each `supp_table_*` block referencing this script for traceability.
- Re-running the script must be idempotent.

## Why a separate script (vs. inline)

The transformation is reproducible from a single script, so a fresh checkout of the repo can regenerate every `_modified.csv` deterministically. Keeps provenance: the original supp file matches what the paper deposited; the derivation lives in code.

## Existing builders

| Paper | Builder | Notes |
|---|---|---|
| Domínguez 2017 | `build_dominguez2017_modified_csv.py` | merges S3 up + S3 down → signed `log2_fold_change` |
| Fuszard 2012 | `build_fuszard2012_modified_csv.py` | adds `log2_fold_change = log2(Ratio of means)` per strain CSV |
| Moreno 2023 (cyano S2 + Marinobacter S4) | `build_moreno2023_modified_csv.py` | adds `log2FC.*` columns; per-strain override for SS120 `FC GmM/L` typo |
| Moreno 2023 (Alteromonas S3) | `build_moreno2023_s3_modified_csv.py` | extracts `DEH24_*` from `Description`, normalizes per-strain column quirks |
| Manck 2020 | `build_manck2020_modified_csv.py` | flattens Data Set S1 (Fe / C sheets) to CSVs; builds the 66-TBDT inventory from Data Set S2 + main-text Table 1 substrates |
| Rubin 2015 | `build_rubin2015_essentiality.py` | sd03 -> rubin2015_essentiality.csv (snake_case essentiality category) |
| Oleza 2012 (not registered) | `build_oleza2012_table_s2.py` | Table S2 panels A-J -> oleza2012_table_s2_wide.csv (signed linear fold -> log2, AC-test p) |
| Koch 2020 | `build_koch2020_unique_genes.py` | Table S7 strain-unique genes -> unique_genes/<strain>_unique_genes.csv (MIT1002, 83-1, AD45, ATCC27126) |
| Barber 2026 | `build_barber2026_rpom_proteomics.py` | R. pomeroyi Perseus Imputed sheet -> barber2026_rpom_proteomics_de.csv; log2fc = -Difference for vs-glucose contrasts |
| Ferrer 2020 | `build_ferrer2020_dss3_rnaseq.py` | DSS3 DESeq sheet -> ferrer2020_dss3_coculture_vs_glucose.csv; log2 fold negated (control minus coculture) |
| Durham 2015 | `build_durham2015_dss3_coculture.py` | ST01 docx -> durham2015_st01_dss3_coculture.csv; linear signed FC -> log2 |
| Tetu 2009 | `build_tetu2009_sam_tables.py` | splits Supp. Table 1 xls (4 sheets) into per-contrast CSVs; SAM numerator -> log2_fold_change, q-value % -> fraction |
| Ostrowski 2010 | `build_ostrowski2010_ptra_tables.py` | ptrA vs WT >2-fold genes: main-text Tables 1/2 from the PDF (pdftotext) + Supp. Table 3 .doc (antiword) |
| Koch 2019 | `build_koch2019_de_nsaf.py` | Table S1 -> koch2019_table_s1.csv; paper "X vs Y" log2FC negated to treatment-vs-control; NSAF "nd" -> blank |
| Oleza 2015b | `build_oleza2015b_exoproteomes.py` | Table S4 -> per-strain NSAF + detection_status CSVs (8 strains, identical-sequence paralogs collapsed); Table S6 -> 4 WH8102 contrast CSVs (signed linear TFold FC -> log2); per-strain protein_seqs_<strain>.faa for diamond id_translations |

## Run

```bash
uv run python scripts/build_modified_csv/build_moreno2023_modified_csv.py
```

Pointers: [memory feedback note on the convention](../../../.claude/projects/-home-osnat-github-multiomics-biocypher-kg/memory/feedback_csv_modified_suffix.md).
