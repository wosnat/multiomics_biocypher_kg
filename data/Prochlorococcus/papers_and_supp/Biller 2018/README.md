# Biller 2018

**Citation:** Biller SJ, Coe A, Roggensack SE, Chisholm SW (2018). Heterotroph interactions alter *Prochlorococcus* transcriptome dynamics during extended periods of darkness. *mSystems* 3(3):e00182-18.
**DOI:** 10.1128/mSystems.00182-18
**Organism(s):** *Prochlorococcus* NATL2A (axenic and coculture); *Alteromonas macleodii* MIT1002 (coculture partner)
**Topic:** RNA-seq time course comparing cultures under extended darkness versus a 13:11 light:dark diel cycle control. Conditions include axenic NATL2A, NATL2A co-cultured with MIT1002, and MIT1002 in coculture. Cultures grown in Pro99 natural seawater medium at 24C, bubbled with air. Additional tables report 24h transcript periodicity (Y/N) and late-darkness transcript survival categories.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `biller-et-al-2018-...pdf` | PDF | Main paper | reference | — |
| `table S3 DE genes natl2a axenic extended darkness sys003182233st3.csv` | CSV | NATL2A DE genes, axenic + coculture, 3 timepoints each; asterisk-encoded significance (DESeq2 P<0.1) | already in | — |
| `table S6B DE genes in extended darkness MIT1002 sys003182233st6.csv` | CSV | MIT1002 DE genes, consistently DE at 1-5h extended darkness (prefiltered) | already in | — |
| `table s4A sys003182233st4.csv` | CSV | NATL2A 24h periodicity Y/N across 4 conditions (axenic/coculture x L:D/darkness) | already in | — |
| `table s4B sys003182233st4.csv` | CSV | MIT1002 24h periodicity Y/N across 2 conditions (coculture L:D, coculture darkness) | already in | — |
| `table S5 sys003182233st5.csv` | CSV | NATL2A darkness survival composite class (3 categorical values) at 72-144h | already in | — |
| `MIT1002_systematicnames_conversiontable.csv` | CSV | RAST → GenBank/assembly ID bridge for MIT1002 | already in | — |
| `MIT1002_RAST_annotation/rast_to_mit1002_id_translation.csv` | CSV | Diamond protein-match RAST fig\| → NCBI locus tag | already in | — |
| `MIT1002_RAST_annotation/226.6.faa` | FASTA | RAST protein FASTA used as diamond source | already in | — |
| `MIT1002_systematicnames_conversiontable.xlsx` | XLSX | XLSX copy of the RAST names conversion table | skip | Duplicate of CSV already in |
| `sys003182233st1.xlsx` | XLSX | Table S1 (sample metadata / summary stats) | skip | Sample-level metadata, not gene evidence |
| `sys003182233st2.xlsx` | XLSX | Table S2 (summary; not per-gene DE) | skip | No per-gene evidence |
| `sys003182233st3.xlsx` | XLSX | Original XLSX of Table S3 | skip | Duplicate — CSV already in |
| `sys003182233st4.xlsx` | XLSX | Original XLSX of Table S4 | skip | Duplicate — CSVs already in |
| `sys003182233st5.xlsx` | XLSX | Original XLSX of Table S5 | skip | Duplicate — CSV already in |
| `DE genes natl2a axenic extended darkness sys003182233st3.xlsx` | XLSX | Alternate XLSX of Table S3 | skip | Duplicate |
| `DE genes in extended darkness MIT1002 sys003182233st6.xlsx` | XLSX | XLSX of Table S6 | skip | Duplicate — S6B CSV already in |
| `periodic genes diel cycle sys003182233st4.xlsx` | XLSX | Alternate XLSX of Table S4 | skip | Duplicate |
| `genes present in extended darkness sys003182233st5.xlsx` | XLSX | Alternate XLSX of Table S5 | skip | Duplicate |
| `41396_2016_bfismej201670_moesm47_esm.docx` | DOCX | Supplementary material from related Biller 2016 paper | reference | Out of scope for this paper |
| `legend for supp tables.txt` | TXT | Legends for S3-S6 | reference | — |
| `cluster_extractions/` | dir | Cluster extraction JSONs (empty / unused — S4/S5 migrated to DerivedMetric) | reference | — |
| `paperconfig.yaml` | YAML | Active paperconfig | reference | — |
| `paperconfig_orig.yaml` | YAML | Earlier paperconfig version | reference | — |

## Current paperconfig summary

- Experiments defined: 3 — `darkness_extended_darkness_natl2a_rnaseq_axenic`, `darkness_extended_darkness_natl2a_rnaseq_coculture`, `darkness_extended_darkness_mit1002_rnaseq`
- Statistical analyses (DE edges): 8 total — 6 in Table S3 (NATL2A axenic x 3tp + coculture x 3tp, asterisk-encoded P<0.1) + 2 in Table S6B (MIT1002 1h, 5h; prefiltered, all significant)
- Supplementary materials entry types: `csv` (S3, S6B), `derived_metrics_table` (S4A split into `s4a_natl2a_axenic` + `s4a_natl2a_coculture`; S4B as `s4b_mit1002`; S5 as `s5_natl2a_survival`), `id_translation` (RAST diamond + systematic names), `annotation_gff` (old draft MIT1002 GFF)
- Organisms covered: Prochlorococcus NATL2A, Alteromonas macleodii MIT1002
- Table scope(s): `significant_any_timepoint` (NATL2A), `significant_only` (MIT1002)
- Non-DE evidence: 4 DerivedMetric entries total (6 boolean periodicity flags across S4A/S4B + 1 categorical class from S5)
- ID resolution: NATL2A via `locus_tag_ncbi` (NATL2_) + `old_locus_tag` (PMN2A_RS*) + `locus_tag_cyanorak` (NCBI ID); MIT1002 via RAST `fig|` IDs bridged through `226.6.faa` diamond match plus 4-digit GenBank IDs via the systematic names table.

## Recommended actions

1. **No action** — Tables S3, S4A, S4B, S5, S6B are all integrated under the correct node types. S4A/S4B/S5 use `derived_metrics_table` (boolean + categorical), which matches the decision rule for Y/N markers and fixed-vocabulary composite labels.
2. **Skip** — all duplicate XLSX copies (S1-S6) carry no gene-level evidence beyond what is already in the ingested CSVs.
3. **Reference** — main PDF and legend TXT are reference-only.

## Notes

- Gene ID formats: NATL2A CSVs carry three parallel ID columns (`NCBI ID` = Cyanorak, `NCBI ID_2` = NCBI, `NCBI ID_3` = old locus tag). MIT1002 S6B uses `RAST_region_ID` bridged via `alternative_locus_tag` through the conversion table; S4B uses 4-digit `MIT1002_0000` GenBank IDs bridged through the same table's `genbank_id` column. The RAST diamond-match file is auto-generated by `prepare_data.sh` step 3 (`generate: diamond_protein_match`).
- `_resolved.csv` / `_resolved_report.txt` files auto-generated by `prepare_data.sh` step 4 — omitted from the inventory above.
- Strain coverage: both organisms are deployed in the KG.
- The S4A split (axenic vs coculture DerivedMetric blocks, same source CSV) is intentional — one DerivedMetric per (Experiment x metric_type) pair. See comments in `paperconfig.yaml` lines 179-189 and `docs/kg-changes/non-de-evidence-extension.md`.
- S5's `darkness_cluster` composite column encodes cross-condition information in the category names (`darkness_axenic+darkness_coculture`, etc.); attributed to the axenic NATL2A experiment as canonical parent.
