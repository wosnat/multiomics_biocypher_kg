# Lindell 2007

**Citation:** Lindell D, Jaffe JD, Coleman ML, Futschik ME, Axmann IM, Rector T, Kettler G, Sullivan MB, Steen R, Hess WR, Church GM, Chisholm SW. "Genome-wide expression dynamics of a marine virus and host reveal features of co-evolution." *Nature* 449, 83–86 (2007).
**DOI:** 10.1038/nature06130
**Organism(s):** *Prochlorococcus* MED4 (host) infected with T7-like cyanophage P-SSP7
**Topic:** Whole-genome transcriptome time course (0–8 h post-infection) of MED4 cells lytically infected by cyanophage P-SSP7, using Affymetrix microarrays. The supplementary table reports fold-change for upregulated host ORFs only, grouped into two transcription clusters (early transient vs. late sustained induction). Phage genes are also transcriptionally profiled but are not loaded into the KG.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `nature06130 (1).pdf` | PDF | Main paper | reference | — |
| `41586_2007_BFnature06130_MOESM283_ESM.pdf` | PDF | Supplementary methods + figures (includes full supp-table listings) | reference | — |
| `supp table legend.txt` | text | Legend for supp table 3 | reference | — |
| `supp table 3.csv` | CSV | 41 upregulated MED4 ORFs × 9 timepoints (0–8 h), log2FC with asterisk-p-value encoding, TRANSCRIPTION GROUP column (1 or 2) | already in | — (used by both `supp_table_3` DE entry and `med4_phage_transcription_groups` gene_clusters entry) |
| `supp table 3_resolved.csv` | CSV | Auto-generated locus_tag-resolved copy (prepare_data step 4) | — | — (auto-generated; skip) |
| `supp table 3_resolved_report.txt` | text | Resolution report | — | — (auto-generated; skip) |
| `cluster_extractions/med4_phage_transcription_groups.json` | JSON | Extracted per-cluster descriptions for the 2 transcription groups | already in | — (feeds GeneCluster node descriptions) |
| `cluster_extractions/med4_phage_transcription_groups.md` | MD | Markdown companion to the JSON | reference | — |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Superseded draft | — | — (skip / historical) |

## Current paperconfig summary

- Experiments defined: 1 — `viral_phage_infected_med4_microarray` (MED4 + Phage P-SSP7, MICROARRAY, treatment_type=[viral], coculture partner taxid 10239).
- Statistical analyses (DE edges): 8 time-point analyses (1h–8h) on a single source CSV (`supp table 3.csv`, 41 rows) → ~328 expression edges into one Experiment. `pvalue_asterisk_in_logfc: true` because significance is encoded as asterisks on log2FC.
- Supplementary materials entry types: 1× `csv` (DE) + 1× `gene_clusters` (2 transcription groups) — both reading from the same `supp table 3.csv`.
- Organisms covered: *Prochlorococcus* MED4 (host); phage P-SSP7 linked via `Tests_coculture_with`.
- Table scope: `significant_only` — upregulated-only table, phage-induced gene set.
- Non-DE evidence: 1 `gene_clusters` analysis (2 transcription groups, time_course; descriptions in `cluster_extractions/`).
- ID resolution: `ORF ` column → `old_locus_tag` (PMM#### style); resolves directly for MED4 (GCF_000011465.1).

## Recommended actions

1. **No action** — integration appears complete. DE edges and gene-cluster nodes both emit from `supp table 3.csv`.
2. **Verify** — confirm match rate for the 41 ORFs via `/check-gene-ids` (trailing space in `ORF ` column is fragile; ensure adapter still strips it).
3. **Consider** — the asterisk-encoded pseudo-p-values give no continuous adjusted_p_value; `expression_status` will only ever be `significant_up` / `not_significant` here.

## Notes

- Gene IDs: `PMM####` old locus tags (with trailing space in column header and trailing space in values — already accommodated by the paperconfig).
- Only upregulated genes are reported; downregulated host genes in the paper's Supplementary Table (1,716 genes) are NOT in this CSV.
- Downregulated host gene lists and all phage gene transcripts live only in `41586_2007_BFnature06130_MOESM283_ESM.pdf` — not machine-readable from this directory.
- Time-course clustering is appropriately modeled as `gene_clusters` here (two named patterns with temporal_pattern descriptions) rather than `derived_metrics_table` — this is correct soft-clustering semantics.
