# Read 2017

**Citation:** Read RW, Berube PM, Biller SJ, Neveux I, Cubillos-Ruiz A, Chisholm SW, Grzymski JJ. "Nitrogen cost minimization is promoted by structural changes in the transcriptome of N-deprived *Prochlorococcus* cells." *The ISME Journal* 11(10), 2267–2278 (2017).
**DOI:** 10.1038/ismej.2017.88
**Organism(s):** *Prochlorococcus* MED4 (axenic)
**Topic:** RNA-seq (Rockhopper-analyzed) of axenic MED4 cultures 3 h, 12 h, 24 h after resuspension in N-depleted vs N-replete Pro99. Focus is on transcription-start-site shifts under N-deprivation (shortened transcripts encoding lower-N proteins); supplementary DE tables cover the top 50% of expressed genes.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `41396_2017_article_bfismej201788.pdf` | PDF | Main paper | reference | — |
| `DE_genes_3h_post_starvation.csv` | CSV | 3 h DE: locus_tag, log2_fold_change, p_value, definition | already in | — |
| `DE_genes_12h_post_starvation.csv` | CSV | 12 h DE | already in | — |
| `DE_genes_24h_post_starvation.csv` | CSV | 24 h DE | already in | — |
| `DE_genes_3h_post_starvation_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `DE_genes_12h_post_starvation_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `DE_genes_24h_post_starvation_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `*_resolved_report.txt` (×3) | text | Resolution reports | — | — (auto-generated; skip) |
| `DE genes top 50 percent 3 hours post starvation 41396_2017_bfismej201788_moesm47_esm.pdf` | PDF | PDF version of the 3 h DE table | skip | Skip (PDF duplicate of `DE_genes_3h_post_starvation.csv`). |
| `41396_2017_bfismej201788_moesm38_esm.docx` | DOCX | Supplementary methods / gene-list narrative | reference | — (text-only supplementary material; no per-gene evidence beyond what is in the CSVs) |
| `41396_2017_bfismej201788_moesm48_esm.pdf` | PDF | Supplementary TSS table (48) | skip | Skip — transcription-start-site offsets per gene; not currently modeled in the KG. Could be added later as a `derived_metrics_table` with numeric TSS-shift metric if of interest. |
| `41396_2017_bfismej201788_moesm49_esm.pdf` | PDF | Supplementary table (49), internal-TSS genes list | skip | Skip — same reason as above. Could become a boolean `derived_metrics_table` ("has_internal_TSS_under_N_dep") if prioritized. |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Superseded draft | — | — (skip / historical) |

## Current paperconfig summary

- Experiments defined: 1 — `nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq`.
- Statistical analyses (DE edges): 3 (3 h, 12 h, 24 h). Each CSV feeds one analysis.
- Supplementary materials entry types: 3× `csv`.
- Organisms covered: *Prochlorococcus* MED4 only.
- Table scope: `filtered_subset` — "top 50% of genes by expression level" (not a significance filter; each timepoint CSV is an independent subset).
- Non-DE evidence: none integrated. TSS-shift / internal-TSS evidence is in PDF-only supplementary tables moesm48/49.
- ID resolution: `locus_tag` column → tier-1 `locus_tag` (direct match on MED4 PMM####). Clean, no bridging needed.

## Recommended actions

1. **No action** on DE — three-timepoint integration is fine.
2. **Add (optional)** — extract moesm48/moesm49 PDF tables into CSV and integrate as `derived_metrics_table` entries: moesm48 = numeric "primary_TSS_shift_bp" per gene (`value_kind: numeric`, `rankable: true`); moesm49 = boolean "has_internal_TSS_under_Ndep" (`value_kind: boolean`). Both attach to the existing experiment. These are the paper's headline findings and are currently invisible in the KG.
3. **Skip** — moesm47 PDF and docx are redundant or narrative.
4. **Verify** — `p_value` column is labeled as `adjusted_p_value_col` in the paperconfig; confirm that Rockhopper's p_value is FDR-corrected (if not, rename semantics).

## Notes

- Gene IDs: native `locus_tag` (MED4 PMM####); no bridging required.
- "Top 50% of genes by expression level" is a pre-filter, not a significance filter — the resulting DE edges include rows with high p-values. `expression_status` will correctly classify based on `p_value_threshold: 0.05`.
- The paper's central biological claim (TSS shifts / shorter transcripts → lower-N proteins) is entirely in moesm48/49 and is currently NOT represented in the KG.
