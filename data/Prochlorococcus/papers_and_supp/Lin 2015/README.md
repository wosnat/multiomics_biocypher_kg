# Lin 2015

**Citation:** Lin et al. 2015, "Transcriptomic response during phage infection of a marine cyanobacterium under P limitation", Environmental Microbiology (supp prefix `emi13104`).
**DOI:** (not listed in paperconfig)
**Organism(s):** *Prochlorococcus* NATL2A; phage (cyanophage infection at 47 h)
**Topic:** RNA-seq time-course of NATL2A under P-limited vs P-replete, with and without phage infection; time points span 4, 24, 46, 50 (P added), 59 h uninfected and 47, 48, 49, 51, 55, 55 (P added), 59, 59 (P added) infected. Tables S4A/S4B carry paired log2FC + padj columns per timepoint.

## Classification

**Bucket D — defer / nothing to do**

The two RNA-seq DE tables (S4A uninfected, S4B infected) covering all reported timepoints are fully integrated as `Changes_expression_of` edges. The remaining unintegrated artefacts are out of scope or inactionable: `Table S2` is sample x gene raw read counts (per-sample abundance, no DE), `Table S1` is sample-level metadata, `Supplementary Figures.doc` and `Table S5.doc` are PDF-style narratives, and the SI zip is a duplicate bundle. `Table S3.xls` was earlier flagged as a possible per-gene annotation table to inspect; if it ever turns out to carry per-gene categorical evidence it would graduate to bucket B, but with no clear actionable per-gene content visible today it stays in D.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `Environmental Microbiology - 2015 - Lin - Transcriptomic response during phage infection ...pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Older config kept for provenance | reference | — |
| `legend of supp tables.txt` | TXT | S4 legend (padj < 0.1 red; |log2FC| > 1 bold) | reference | — |
| `Supplementary Figures.doc` | DOC | Supp figures | reference | — |
| `Table S1.docx` | DOCX | Sample/experiment metadata table | skip | — (sample-level metadata, no gene evidence) |
| `emi13104-sup-0001-si.zip` | ZIP | Bundle of original supplementary files | skip | — (redundant; individual CSVs extracted into this directory) |
| `read counts Table S2.xlsx` | XLSX | S2 — raw/normalised read counts per sample | skip | — (sample × gene counts; out of scope) |
| `Table S3.xls` | XLS | S3 — likely annotation/summary table | add | Inspect: if per-gene categorical/numeric annotation (e.g. KEGG pathway assignment), add as `derived_metrics_table`; if sample-level or figure-only, reclassify to skip. |
| `DE genes P Table S4.xls` | XLS | Original S4 Excel workbook (A=uninfected, B=infected) | skip | — (split into CSVs already) |
| `table S4A.csv` | CSV | S4A — NATL2A uninfected: paired log2FC + padj across 6 timepoints (4/24/46/59 h -P; 50/59 h P-re-added) | already in | — |
| `table S4B.csv` | CSV | S4B — NATL2A infected: paired log2FC + padj across 8 timepoints (47/48/49/51/55/59 h -P; 55/59 h P-re-added) | already in | — |
| `Table S5.doc` | DOC | Supp table 5 narrative (likely gene-list discussion or pathway summary) | reference | — (inspect if it contains a tractable per-gene table — currently treated as reference) |

Auto-generated `*_resolved.csv` / `*_resolved_report.txt` are produced by prepare_data step 4.

## Current paperconfig summary

- Experiments defined: **2** (`phosphorus_plimited_natl2a_rnaseq_uninfected`, `phosphorus_plimited_natl2a_rnaseq_infected`)
- Statistical analyses (DE edges): **14** — 6 uninfected (4h, 24h, 46h, 59h -P; 50h, 59h P-re-added) + 8 infected (47h, 48h, 49h, 51h, 55h, 59h -P; 55h, 59h P-re-added)
- Supplementary materials entry types: `csv` × 2 (S4A, S4B)
- Organisms covered: NATL2A; infected experiment carries `viral` background factor (phage lysis at 47 h)
- Table scope(s): `significant_only` (on the Experiment); per-analysis `pvalue_threshold` defaults to the global (0.05) — inspect before trusting counts
- Non-DE evidence: none
- ID resolution: `ID` → `old_locus_tag` (PMN2A_XXXX format), `Gene Name` → `gene_name`

## Recommended actions

1. **Add** — Inspect `Table S3.xls`. If it is a per-gene annotation/category table, wire it as a `derived_metrics_table` bound to one of the two Experiments (probably the uninfected baseline). If it is a sample-level or figure-only table, mark `skip` in a follow-up.
2. **Reference / inspect** — `Table S5.doc` may summarise pathway-level observations. If it embeds a per-gene table, extract it and evaluate for integration. Otherwise leave as reference.
3. **Skip** — `read counts Table S2.xlsx`, `Table S1.docx`, and `emi13104-sup-0001-si.zip` are all out of scope (sample counts, sample metadata, zipped duplicate bundle).
4. **No action** — S4A + S4B paired-column time courses are fully integrated (14 DE analyses).

## Notes

- Infected experiment uses `background_factors: [viral, diel]`; the phage is the background condition (infection at 47 h). It is not a `coculture` experiment in the KG's treatment-category sense.
- "P added" timepoints (50/55/59 h) use `growth_phase: recovery` and share the same Experiment node as the P-limitation analyses, but are distinct analysis edges.
- S4A and S4B pair `N h (-P/control)` fold-change columns with matching `N h (-P/control) padj` columns — the omics adapter resolves these per-timepoint independently.
- No DOI is carried in paperconfig; filename suggests Env. Microbiol. 2015, so confirm before any citation-sensitive downstream use.
