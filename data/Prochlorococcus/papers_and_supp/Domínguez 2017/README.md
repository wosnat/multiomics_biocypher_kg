# Domínguez-Martín 2017

**Citation:** Domínguez-Martín et al. 2017, "Quantitative proteomics shows extensive remodeling induced by nitrogen limitation in Prochlorococcus marinus SS120", mSystems 2:e00008-17.
**DOI:** 10.1128/mSystems.00008-17
**Organism(s):** *Prochlorococcus marinus* SS120 (CCMP1375)
**Topic:** Label-free LC-MS/MS proteomics of SS120 exposed to azaserine (10 uM, a ferredoxin-GOGAT inhibitor used as an acute N-limitation proxy) vs untreated control, cells harvested 24 h after addition; Progenesis QI quantitation with n=3 biological replicates per condition. The paper reports only significantly changing proteins (ANOVA p < 0.05).

## Classification

**Bucket D - already integrated, nothing to add**

The combined merged S3 (`table s3 Combined_modified.csv`) is wired as a single `Changes_expression_of` analysis (`ss120_azaserine_vs_control`) with signed `log2_fold_change` derived from `Max fold change` and ANOVA p-values. SS120 is deployed as a reference-proteome-only strain. S1 (identification list), S2 (BLAST of uncharacterised ORFs), and S4 (carbohydrate subset of S3) carry no quantitative DE evidence beyond what's already in. No metabolomics, no clustering, no remaining gene-level data to add.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `domínguez-martín-et-al-2017-...in.pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | reference | — |
| `legends.txt` | TXT | Supplementary table legends | reference | — |
| `table s3 Combined_modified.csv` | CSV | Merged up+down S3 with `log2_fold_change` signed from `Max fold change` × direction; `extracted_gn` parsed from UniProt `GN=` field | already in | — |
| `table s3 Upregulated prots rel quant sys003172107st8.csv` | CSV | S3 (upregulated) relative quant — source for the combined file | skip | — (duplicated by `_modified.csv`) |
| `table s3 Downregulated prots rel quant sys003172107st8.csv` | CSV | S3 (downregulated) relative quant — source for the combined file | skip | — (duplicated by `_modified.csv`) |
| `sys003172107st8.xlsx` | XLSX | S3 Excel original (up+down sheets) | skip | — (duplicate of split CSVs) |
| `sys003172107st6.pdf` | PDF | S1 — proteins identified in both conditions | skip | — (identification list, no quant/DE) |
| `sys003172107st7.xlsx` | XLSX | S2 — uncharacterised protein BLAST results | skip | — (not gene-level DE) |
| `sys003172107st9.xlsx` | XLSX | S4 — ribosomal + carbohydrate metabolism protein quantification | skip | — (subset duplicate of S3 quant) |
| `table s4 Carbohydrate metabolism prots sys003172107st9.csv` | CSV | S4 carbohydrate subset | skip | — (subset of S3) |
| `sys003172107sf1.pdf` | PDF | Supplementary figure 1 | reference | — |
| `sys003172107sf2.pdf` | PDF | Supplementary figure 2 | reference | — |

Auto-generated `*_resolved.csv` / `*_resolved_report.txt` are produced by prepare_data step 4 and are not inventoried individually.

## Current paperconfig summary

- Experiments defined: **1** (`azaserine_ss120_proteomics`)
- Statistical analyses (DE edges): **1** (`ss120_azaserine_vs_control`) — only significant proteins reported
- Supplementary materials entry types: `csv` × 1 (the merged `_modified` table)
- Organisms covered: SS120
- Table scope(s): `significant_only`
- Non-DE evidence: none
- ID resolution: `Accession` → UniProt entry name (Tier 1); `extracted_gn` → gene_name (Tier 3 singleton fallback). See extensive comment in paperconfig explaining why `Description` can NOT be declared as locus_tag (tokenisation would collapse edges).

## Recommended actions

1. **No action** — integration is working; leave `_modified.csv` as the single canonical CSV and the split up/down CSVs as raw provenance.
2. **No action** — S4 carbohydrate subset is a figure-adjacent slice of S3; no new gene-level evidence beyond S3.
3. **No action** — S1 (identified proteins list) and S2 (BLAST of uncharacterised ORFs) do not contain quantitative DE evidence.

## Notes

- SS120 is deployed as a reference-proteome-only strain; gene IDs in S3 resolve via UniProt entry name (`*_PROMA`) or via the `GN=` token from the Description field.
- No per-protein adjusted p-value column beyond Anova(p); `pvalue_threshold: 0.05` is applied as-is.
- Do not re-declare `Description` as `locus_tag` — the paperconfig comment documents a known edge-count regression (~400 → 1 edges) if you do.
- See `docs/kg-changes/reference-proteome-match-organisms.md` for SS120 deployment details.
