# Fang 2019

**Citation:** Fang et al. 2019, "Transcriptomic responses of *Prochlorococcus* MIT9313 to viral lysis products (vDOM)", Environmental Microbiology (supp file prefix `emi14513`).
**DOI:** (not listed in paperconfig; filename suggests Env. Microbiol. 2019)
**Organism(s):** *Prochlorococcus* MIT9313
**Topic:** RNA-seq time course of axenic MIT9313 in Pro99 medium amended with viral DOM (vDOM, released from phage-lysed MED4) vs Pro99-only control; samples at 0.5, 1, 2, 4, 8, 12, 24, 48, 72 h. Significance is carried in the source table as asterisks on fold-change cells (adjusted p < 0.1 and |log2FC| > 1).

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `Fang 2015 Prochlorococcus mit9313 ... viral lysis products.pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Older / pre-refactor config kept for provenance | reference | — |
| `legend of supp tables.txt` | TXT | S2 legend: asterisk = adjusted-p < 10% and |log2FC| > 1 | reference | — |
| `de genes emi14513-sup-0002-tables2.csv` | CSV | S2 — MIT9313 gene × timepoint log2FC matrix (9 timepoints, asterisk-significance encoding); 5 header rows | already in | — |
| `de genes emi14513-sup-0002-tables2.xlsx` | XLSX | Excel version of S2 | skip | — (duplicate of CSV) |
| `read counts emi14513-sup-0005-tables5.xlsx` | XLSX | S5 — raw/normalised read-count matrix | skip | — (sample-level counts; out of scope) |

Auto-generated `*_resolved.csv` / `*_resolved_report.txt` are produced by prepare_data step 4.

## Current paperconfig summary

- Experiments defined: **1** (`viral_viral_dom_vdom_addition_mit9313_rnaseq`) — viral treatment category, Phage treatment organism (taxid 10239)
- Statistical analyses (DE edges): **9** (one per timepoint: 0.5h, 1h, 2h, 4h, 8h, 12h, 24h, 48h, 72h)
- Supplementary materials entry types: `csv` × 1
- Organisms covered: MIT9313
- Table scope(s): `significant_any_timepoint`
- Non-DE evidence: none
- ID resolution: `Gene ID` → `old_locus_tag` (PMT format), `Gene` → `gene_name`; uses `pvalue_asterisk_in_logfc: true` so significance is parsed from `*` suffix on fold-change cells

## Recommended actions

1. **No action** — full time-course DE is already integrated across 9 timepoints.
2. **Skip** — S5 read-count matrix is sample-level and not suitable for KG integration.
3. **No action** — `paperconfig_orig.yaml` is historical provenance; keep it but do not reference it from the pipeline.

## Notes

- Significance encoding is unusual: the original Excel highlights `*` inside fold-change cells. The omics adapter parses this via `pvalue_asterisk_in_logfc: true` and extracts `logfc_threshold: 1.0` + implicit `adjusted_p_value < 0.1` semantics.
- No dedicated adjusted-p column is present in S2; `adjusted_p_value` on edges will be null (permitted — see CLAUDE.md).
- Treatment is "viral" (vDOM addition, not live phage challenge) but `treatment_organism: Phage` (taxid 10239) is carried on the Experiment for downstream coculture/viral linkage.
