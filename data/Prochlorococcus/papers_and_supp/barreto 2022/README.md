# Barreto Filho 2022

**Citation:** Barreto Filho MM, Lu Z, Walker M, Morris JJ. "Community context and pCO2 impact the transcriptome of the 'helper' bacterium *Alteromonas* in co-culture with picocyanobacteria." *ISME Communications* 2, 113 (2022).
**DOI:** 10.1038/s43705-022-00197-2
**Organism(s):** *Prochlorococcus* MIT9312, *Synechococcus* CC9311 and WH8102, *Alteromonas macleodii* EZ55 (3 phototrophs × EZ55 heterotroph × 400 vs 800 ppm pCO2, 14:10 L:D)
**Topic:** Factorial RNA-seq (edgeR) examining how elevated pCO2 (800 vs 400 ppm) and choice of phototroph partner alter EZ55's transcriptome, and how pCO2 reshapes each cyanobacterium's response. Reveals pCO2-linked changes in EZ55 carbohydrate metabolism, stress response, chemotaxis, and transporters; links reduced oxidative-stress help to Prochlorococcus growth loss at high pCO2.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `43705_2022_article_197.pdf` | PDF | Main paper | reference | — |
| `43705_2022_197_moesm1_esm.pdf` | PDF | Supplementary methods + table legends + Figs. S1-S12 | reference | — |
| `pro_9312_anot.csv` | CSV | MIT9312 annotation bridge (GID → old_locus_tag, gene_name, uniprot_acc) | already in (id_translation) | — |
| `pro_9312_anot.xlsx` | XLSX | Excel copy of the above | skip | Duplicate |
| `syn_9311.csv` | CSV | CC9311 annotation bridge (GID, genename, uniprot_acc, definition) | already in (id_translation) | — |
| `syn_9311.xlsx` | XLSX | Excel copy | skip | Duplicate |
| `syn_8102_anot.csv` | CSV | WH8102 annotation bridge | already in (id_translation) | — |
| `syn_8102_anot.xlsx` | XLSX | Excel copy | skip | Duplicate |
| `EZ55.exon.fixed2.gtf` | GTF | EZ55 annotation (bridge for Alteromonas mapping) | already in (annotation_gff) | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S1_MIT9312.csv` | CSV | Per-strain pCO2 DE table (MIT9312) | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S1_CC9311.csv` | CSV | Per-strain pCO2 DE table (CC9311) | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S1_WH8102.csv` | CSV | Per-strain pCO2 DE table (WH8102) | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S2.csv` | CSV | EZ55 pCO2 DE (main effect) | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S9.csv` | CSV | EZ55 coculture DE vs MIT9312 at 400 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S10.csv` | CSV | EZ55 coculture DE vs CC9311 at 400 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S11.csv` | CSV | EZ55 coculture DE vs WH8102 at 400 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S12.csv` | CSV | EZ55 coculture DE vs MIT9312 at 800 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S13.csv` | CSV | EZ55 coculture DE vs CC9311 at 800 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S14.csv` | CSV | EZ55 coculture DE vs WH8102 at 800 ppm | already in | — |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S3.csv` | CSV | EZ55 DE between coculture conditions (partner-vs-partner contrasts): `species,symbol,product,logFC,PValue` | add | Add as `csv` DE entries — one statistical_analysis per partner-pair contrast (needs a treatment/control split if table encodes multiple contrasts) |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S4.csv` | CSV | EZ55 pCO2 DE affected by coculture (pCO2 × coculture interaction effect): `species,symbol,product,logFC,PValue` | add | Add as `csv` interaction-effect statistical_analysis; attach to EZ55 experiment with explanatory `experimental_context` |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S5.csv` | CSV | Filtered subset of EZ55 stress-related genes from S9-S14: `partner,treatment,symbol,product,logFC,logFC_rel_axenic,Pvalue` | skip | Subset of S9-S14 already integrated; no new per-gene evidence |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S6.csv` | CSV | Filtered subset — transporter-related genes | skip | Subset of S9-S14 |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S7.csv` | CSV | Filtered subset — chemotaxis-related genes | skip | Subset of S9-S14 |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S8.csv` | CSV | Filtered subset — metabolism-related genes | skip | Subset of S9-S14 |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S1.csv` | CSV | Combined cyanobacteria pCO2 DE table (merge of S1_MIT9312/S1_CC9311/S1_WH8102) | skip | Redundant with per-strain tables already integrated |
| `transcriptome_supplemental_tables/Supplemental_tables/table_S*_with_locus_tag.csv` | CSV | Pre-mapped versions with added locus_tag column | skip | Superseded by pipeline step 4 `_resolved.csv` outputs |
| `43705_2022_197_moesm2_esm.xlsx`, `PRO DE genes...moesm2_esm.csv` | XLSX/CSV | Xlsx + CSV of Prochlorococcus pCO2 DE (already covered by S1_MIT9312) | skip | Duplicate |
| `43705_2022_197_moesm3_esm.xlsx`, `ALT DE genes elevated pCO2 vs ambient...moesm3_esm.csv/xlsx` | XLSX/CSV | EZ55 pCO2 DE (= table_S2) | skip | Duplicate of S2 already integrated |
| `43705_2022_197_moesm4_esm.xlsx`, `ALT DE genes CC vs AX...moesm4_esm.csv` | XLSX/CSV | EZ55 CC vs AX summary (= S9-S14 combined) | skip | Redundant with S9-S14 |
| `43705_2022_197_moesm5_esm.xlsx` | XLSX | Supplementary Excel (likely = one of S3/S4) | skip | Duplicate; keep S3/S4 CSVs as canonical source |
| `43705_2022_197_moesm6_esm.xlsx`, `ALT DE genes CC vs axenic...moesm6_esm.xlsx` | XLSX | Per-setting EZ55 coculture-vs-axenic DE | skip | Redundant with S9-S14 |
| `43705_2022_197_moesm7_esm.xlsx` – `moesm15_esm.xlsx` | XLSX | Supplementary tables S5-S13 in Excel form | skip | Excel versions of already-inventoried CSVs |
| `paperconfig.yaml` / `paperconfig_orig.yaml` | YAML | Integration config + original | reference | — |
| `mail.md` | Markdown | Author correspondence / deposit notes | reference | — |
| `Alteromonas_analysis.zip` + `Alteromonas_analysis/` | ZIP + DIR | Raw R analysis: fastq→counts, edgeR RData, KO xlsx, ORA/GSEA txt, `EZ55.fasta`, custom `org.AspEZ55.eg.db` SQLite, targets.xlsx | skip | Raw sample-level counts + analysis artefacts; not per-gene evidence beyond what's in S9-S14 |
| `CC9311_analysis.zip` + `CC9311_analysis/` | ZIP + DIR | As above for CC9311 | skip | Sample-level raw counts |
| `MIT9312_analysis.zip` + `MIT9312_analysis/` | ZIP + DIR | As above for MIT9312 | skip | Sample-level raw counts |
| `WH8102_analysis.zip` + `WH8102_analysis/` | ZIP + DIR | As above for WH8102 | skip | Sample-level raw counts |
| `transcriptome_supplemental_tables.zip` | ZIP | Zipped copy of the `Supplemental_tables/` directory | skip | Redundant with extracted directory |
| `transcriptome_supplemental_tables/__MACOSX/` | DIR | macOS resource forks from ZIP extraction | skip | Archive artefact |
| `*_resolved.csv` / `*_resolved_report.txt` | CSV/TXT | Auto-generated by step 4 | reference | — |

## Current paperconfig summary

- Experiments defined: 10 — 4 cyanobacteria-vs-EZ55 pCO2 experiments (MIT9312/CC9311/WH8102/EZ55) + 6 EZ55 coculture-vs-axenic experiments spanning 3 partners × 2 pCO2 levels
- Statistical analyses (DE edges): 10 entries across tables S1 × 3 strains, S2, S9, S10, S11, S12, S13, S14
- Supplementary materials entry types: `id_translation` (×3: pro_9312, syn_9311, syn_8102), `annotation_gff` (EZ55 GTF), `csv` (×10)
- Organisms covered: Prochlorococcus MIT9312, Synechococcus CC9311, Synechococcus WH8102, Alteromonas EZ55
- Table scope(s): `significant_only` on all 10 entries
- Non-DE evidence: none today; tables S3 and S4 are candidate additional DE entries (partner-vs-partner and interaction effects)
- ID resolution: MIT9312 uses `GID` = `old_locus_tag` (P9301_*/P9312_*-family); CC9311/WH8102 use `GID` = `locus_tag` directly; EZ55 via GTF bridge. Annotation tables also feed gene_name and uniprot_accession via v2 multi_lookup (Tier 2/3).

## Recommended actions

1. **Add** — Integrate `table_S3.csv` as one or more `csv` DE entries (partner-pair contrasts in EZ55, e.g. CC9311-coculture vs MIT9312-coculture). Needs explicit `treatment_condition`/`control_condition` per contrast because the source rows tag only the overall test; may require the author's S3 legend to split by contrast axis.
2. **Add** — Integrate `table_S4.csv` as a single `csv` entry representing the pCO2 × coculture interaction effect in EZ55, attached to a new dedicated experiment node (e.g. `interaction_pco2_coculture_ez55_rnaseq`) with clear `experimental_context` noting this is an interaction effect, not a main-effect contrast.
3. **Skip** — S5-S8 are filtered subsets of S9-S14 already integrated; xlsx versions of integrated CSVs are duplicates; `*_with_locus_tag.csv` are superseded by the pipeline's `_resolved.csv`; the four `*_analysis.zip` archives and their extractions are sample-level raw counts + analysis artefacts.
4. **No action** — No periodicity/clustering/classifier tables to migrate; main PDF + supplementary PDF are reference-only.

## Notes

- Gene ID formats: MIT9312 uses GID (old_locus_tag, old P9301-style) via `pro_9312_anot.csv`; CC9311 and WH8102 use GID that maps directly to `locus_tag`; EZ55 genes use `EZ55_#####` locus tags bridged by the GTF.
- Strains in/out of KG: MIT9312, CC9311, WH8102, EZ55 all deployed.
- Anomalies: many `*_with_locus_tag.csv` files are legacy pre-mapped copies from before pipeline step 4; they should not be re-integrated. The `__MACOSX` folder is an extraction artefact.
- Tables S3 and S4 are structured differently from S1/S2/S9-S14 (no per-strain split; `species` column marks the organism per row) — the adapter's single-organism assumption means each contrast must be split to one CSV per organism before wiring up.
