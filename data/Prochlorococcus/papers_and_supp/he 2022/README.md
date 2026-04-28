# He 2022

**Citation:** He et al. 2022, "Transcriptomic response of *Prochlorococcus* NATL1A and MED4 to low-salinity acclimation", Frontiers in Microbiology (fmicb-13-1038136).
**DOI:** (not listed in paperconfig; Frontiers article ID 1038136)
**Organism(s):** *Prochlorococcus* NATL1A, MED4
**Topic:** RNA-seq of axenic cultures acclimated through five rounds of transfer in Pro99 medium with reduced salinity (28 psu) vs normal seawater salinity (34 psu); 21C, continuous light (10 uE/m2/s). Tables S1/S2 report only highly DE genes (p < 0.05 and |log2FC| > 1) per strain.

## Classification

**Bucket D - already integrated, nothing to add**

Both strain-specific DE tables (S1 NATL1A, S2 MED4) are wired as `Changes_expression_of` edges with `prefiltered: true`. Both organisms are deployed. The supplementary `Data Sheet 1.docx` is a narrative file with no machine-readable per-gene table. No metabolomics, no clustering. The only conditional follow-up would be if the authors later release a full unfiltered gene x logFC matrix -- not currently available.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `fmicb-13-1038136.pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Older config kept for provenance | reference | — |
| `supp table legends.csv` | CSV | Legends text | reference | — |
| `Data Sheet 1.docx` | DOCX | Supplementary methods / results narrative | reference | — |
| `table s1 Highly differentially expressed genes ... NATL1A .csv` | CSV | S1 — NATL1A DE genes (significant subset) | already in | — |
| `table s2 Highly differentially expressed genes ... MED4 .csv` | CSV | S2 — MED4 DE genes (significant subset) | already in | — |

Auto-generated `*_resolved.csv` / `*_resolved_report.txt` are produced by prepare_data step 4.

## Current paperconfig summary

- Experiments defined: **2** (`salt_low_salinity_acclimation_28_natl1a_rnaseq`, `salt_low_salinity_acclimation_28_med4_rnaseq`)
- Statistical analyses (DE edges): **2** (one per strain; single "acclimated" timepoint)
- Supplementary materials entry types: `csv` × 2
- Organisms covered: NATL1A, MED4
- Table scope(s): `significant_only` (paper's p < 0.05 AND |logFC| > 1 pre-filter); `prefiltered: true`
- Non-DE evidence: none
- ID resolution: `Gene Name` column used for both `old_locus_tag` and `gene_name` (mixed contents: locus tags like `PMN2A_XXXX` or symbolic names)

## Recommended actions

1. **No action** — both strain tables are integrated as significant-only DE, correctly flagged with `prefiltered: true`.
2. **No action** — Data Sheet 1 is a narrative docx without gene-level evidence.
3. **Consider** — if authors publish a full (unfiltered) gene × logFC table later, switch to `table_scope: all_detected_genes` to recover background distribution.

## Notes

- Both tables are pre-filtered to "highly DE" genes (significant-only scope).
- `logfc_col` differs slightly: NATL1A uses `Log2FC`, MED4 uses `logFC` — both mapped correctly in their respective statistical_analyses entry.
- `Gene Name` column holds mixed ID types; declaring it as both `old_locus_tag` and `gene_name` lets resolution fall through Tier 1 → Tier 3 singletons.
- `adjusted_p_value_col: p-value` — paper labels this as padj-equivalent in the legend.
