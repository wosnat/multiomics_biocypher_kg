# Capovilla 2023 — Chitin utilization by marine picocyanobacteria

**Paper**: Capovilla, Braakman, Fournier et al. (2023). "Chitin utilization by marine picocyanobacteria and the evolution of a planktonic lifestyle." *PNAS* 120(20): e2213271120.

**DOI**: https://doi.org/10.1073/pnas.2213271120

## Summary

The paper investigates chitin utilization by marine picocyanobacteria (Prochlorococcus and Synechococcus). It shows that certain lineages can attach to and use particulate chitin, proposes the "chitin raft hypothesis" for the evolutionary origin of the planktonic lifestyle, and characterizes gene expression responses to chitosan addition.

## Data files

| File | Description |
|---|---|
| `DE genes RNA-Seq_MIT9313 pnas.2213271120.sd02.csv` | Significantly DE genes for MIT9313 (83 genes: 49 up, 34 down) |
| `DE genes RNA-Seq_MIT9303 pnas.2213271120.sd02.csv` | Significantly DE genes for MIT9303 (36 genes: 18 up, 18 down) |
| `DE genes pnas.2213271120.sd02.xlsx` | Original Excel with both sheets |
| `read counts pnas.2213271120.sd04.csv` | Normalized read counts for MIT9313 (all genes, all samples) |
| `read counts pnas.2213271120.sd05.csv` | Normalized read counts for MIT9303 (all genes, all samples) |
| `pnas.2213271120.sd01.xlsx` | Chitin utilization gene annotations across genomes |
| `metabolites pnas.2213271120.sd03.xlsx` | Metabolomics data (not used in KG) |

## Experimental design

- **Strains**: Prochlorococcus MIT9303 (primary degrader, full chitin pathway) and MIT9313 (secondary degrader, partial pathway — lacks chitinase)
- **Treatment**: Chitosan addition (56 ug/mL) to Pro99 AMP1 medium
- **Control**: Pro99 AMP1 medium without chitosan
- **Conditions**: Continuous light at 11 umol photons m-2 s-1, 24C
- **Statistical test**: DESeq2, adjusted p-value < 0.1
- **Timepoints**: Samples collected on days 1 and 3 after inoculation; DE analysis compares chitosan vs control across these timepoints (not a time-course design in the output)

## CSV structure notes

Both DE CSVs have been reformatted to flat format with a strain name on row 1 (skipped via `skip_rows: 1` in the paperconfig) and merged up/down-regulated sections. Original xlsx had separate sections.

**Columns**: `log2 fold change (+ means upregulated in chitosan)`, `gene` (Cyanorak-style ID), `ncbi_cds_locus_tag`, `ncbi_gene_locus_tag`, `ncbi_gene_old_locus_tag`, `geneFunction`, `NCBI alignment`, `NCBI query cover`, `NCBI percent identity`

**No adjusted p-values in the CSV** — the paper states genes were considered significant at adjusted p < 0.1 (DESeq2 with Benjamini-Hochberg correction), but the actual p-values are not included in the supplementary table. Only significantly DE genes are listed (table_scope: significant_only).

Some `ncbi_gene_locus_tag` values are `nan` (genes not found in current NCBI annotation). These will fail to map.

## Gene ID mapping

- **MIT9313**: `ncbi_gene_locus_tag` uses `AKG35_RS*` format (NCBI RS-style locus tags). The KG uses `PMT*` (Cyanorak) as primary locus_tag, but `AKG35_RS*` is stored as `locus_tag_ncbi` in gene_mapping.csv and should be resolvable via gene_id_mapping.json (tier 1 as locus_tag).
- **MIT9303**: `ncbi_gene_locus_tag` uses `P9303_RS*` format. MIT9303 is NOT in the KG yet (not in cyanobacteria_genomes.csv), so these genes cannot be mapped until MIT9303 is added as a genome.
- **`ncbi_gene_old_locus_tag`** for MIT9313 contains comma-separated old locus tags (e.g., `PMT0148,PMT_0148,RG24_RS00740`) which are alternative resolution paths.

## Decisions and concerns

1. **MIT9303 is NOT in the KG** — the organism is not listed in `cyanobacteria_genomes.csv`. The MIT9303 experiment is included in the paperconfig for completeness but will not produce edges until MIT9303 is deployed as a genome. This requires finding/adding the NCBI GCF accession, taxid, and running the genome data pipeline.

2. **No adjusted p-values** — `adjusted_p_value_col` is omitted from the paperconfig. Expression edges will have null p-values. The `pvalue_threshold: 0.1` documents the paper's significance criterion but won't be used for filtering (all rows are already significant).

3. **CSV reformatted** — Original xlsx had separate up/down-regulated sections with repeated headers. CSVs have been reformatted to flat format with strain name on row 1 (skipped via `skip_rows: 1`).

4. **Treatment type = carbon** — Chitosan is a partially deacetylated form of chitin, used as a supplemental carbon and energy source. Categorized as `carbon` treatment_type. Could also arguably be `chemical` but carbon is more biologically meaningful since the paper frames it as a mixotrophic carbon source.

5. **`nan` values in gene ID columns** — Several genes (7 in MIT9303, 6 in MIT9313) have `nan` for `ncbi_gene_locus_tag`, indicating Cyanorak genes not found in current NCBI annotation. These will be unresolvable.

## TODO before integration

- [x] Pre-process CSVs into flat single-header format (merged, skip_rows: 1)
- [ ] Deploy MIT9303 genome to the KG (add to cyanobacteria_genomes.csv, run prepare_data.sh)
- [ ] Verify MIT9313 `AKG35_RS*` gene IDs resolve correctly via `/check-gene-ids`
- [ ] Add paperconfig path to `paperconfig_files.txt`
