# Fuszard 2012

**Citation:** Fuszard, Wright & Biggs 2012, "Comparative quantitative proteomics of the marine cyanobacterium *Prochlorococcus* MED4 [three ecotypes] in response to phosphate availability", Aquatic Biosystems 8:7.
**DOI:** 10.1186/2046-9063-8-7
**Organism(s):** *Prochlorococcus* MIT9312 (HL-II), NATL2A (LL-I), SS120 (LL-II/III — reference proteome deployment)
**Topic:** iTRAQ 4-plex quantitative proteomics of three *Prochlorococcus* ecotypes grown in Pro99-based medium, Pi replete (50 uM) vs Pi deplete (10 uM), n=3 biological replicates. Significance: fold-change ratio > 1.6 (up) / < 0.6 (down); paper does not report per-protein p-values. Table S1 lists only significantly changing proteins per strain.

## Classification

**Bucket D - already integrated, nothing to add**

All three strain-specific iTRAQ tables (MIT9312, NATL2A, SS120) are wired as `Changes_expression_of` edges with precomputed signed `log2_fold_change`. SS120 was added as a reference-proteome-match organism in 2026-04-15 specifically to support this paper plus Dominguez 2017. `adjusted_p_value` is null (paper reports only fold-change cutoffs); `logfc_threshold = 0.678` documents the paper's |log2(1.6)| rule. No remaining gene-level evidence; no metabolomics; no clustering.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `2046-9063-8-7.pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | reference | — |
| `12999_2011_8_MOESM2_ESM.DOC` | DOC | Supplementary methods/results text | reference | — |
| `12999_2011_8_MOESM3_ESM.DOC` | DOC | Supp figure/legend | reference | — |
| `12999_2011_8_MOESM3_ESM.pdf` | PDF | Same as above (rendered) | reference | — |
| `12999_2011_8_MOESM1_ESM.XLS` | XLS | Original multi-sheet S1 workbook (MIT9312 + NATL2A + SS120 quantitation) | skip | — (split into per-strain CSVs) |
| `table s1 MIT9312 Quantitation data.csv` | CSV | S1 MIT9312 sheet — raw quant | skip | — (replaced by `_modified.csv`) |
| `table s1 MIT9312 Quantitation data_modified.csv` | CSV | MIT9312 quant + precomputed `log2_fold_change` from `Ratio of means` | already in | — |
| `table s1 NATL2A Quantitation data.csv` | CSV | S1 NATL2A sheet — raw quant | skip | — (replaced by `_modified.csv`) |
| `table s1 NATL2A Quantitation data_modified.csv` | CSV | NATL2A quant + `log2_fold_change` | already in | — |
| `table s1 SS120 Quantitation data.csv` | CSV | S1 SS120 sheet — raw quant | skip | — (replaced by `_modified.csv`) |
| `table s1 SS120 Quantitation data_modified.csv` | CSV | SS120 quant + `log2_fold_change` | already in | — |

Auto-generated `*_resolved.csv` / `*_resolved_report.txt` are produced by prepare_data step 4.

## Current paperconfig summary

- Experiments defined: **3** (`pi_limitation_mit9312_itraq`, `pi_limitation_natl2a_itraq`, `pi_limitation_ss120_itraq`)
- Statistical analyses (DE edges): **3** (one per strain)
- Supplementary materials entry types: `csv` × 3
- Organisms covered: MIT9312, NATL2A, SS120
- Table scope(s): `significant_only` (paper's fold-change cutoff)
- Non-DE evidence: none
- ID resolution: `ORF` → `old_locus_tag` (MIT9312: `PMT9312_XXXX`) or `locus_tag` (NATL2A, SS120); `Gene` → `gene_name` fallback. `log2_fold_change` precomputed by `scripts/build_modified_csv/build_fuszard2012_modified_csv.py` from `Ratio of means`.

## Recommended actions

1. **No action** — all three strain-specific tables are integrated; SS120 is deployed (via reference proteome) as of 2026-04-15.
2. **No action** — the original S1 Excel and raw CSVs are preserved for provenance. Only the `_modified.csv` variants are pipeline-facing.
3. **Consider** — paper reports no per-protein p-values; `adjusted_p_value` is null on all 3 edges. `logfc_threshold = 0.678` (log2(1.6)) applies instead.

## Notes

- Paper's significance rule: fold-change > 1.6 OR < 0.6, based on n=3 replicates; no multiple-testing correction reported.
- Log2 transform is applied in the modified CSV: positive `log2_fold_change` = deplete > replete.
- SS120 deployment as reference-proteome-match organism (see `docs/kg-changes/reference-proteome-match-organisms.md`) — SS120 added 2026-04-15 to support this paper + Domínguez 2017.
- MIT9312 uses `old_locus_tag` (PMT9312_); NATL2A and SS120 use `locus_tag` directly.
