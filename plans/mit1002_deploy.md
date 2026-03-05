# MIT1002 Deploy Plan — Cross-Assembly Gene ID Mapping

## Status: COMPLETE (2026-03-05)

## Problem Summary

Alteromonas macleodii MIT1002 has **two independent genome assemblies** from different research groups, and the Biller 2016/2018 papers use gene IDs from the older draft assembly that cannot be resolved to the current assembly without a verified cross-assembly bridge.

## Two Assemblies

| Property | GCF_001077695.1 (old draft) | GCF_901457835.2 (current) |
|---|---|---|
| Submitter | MIT (Chisholm Lab) | ICBM (Ahrens group, Oldenburg) |
| Year | 2015 | 2022 |
| Assembly level | Contig (draft) | Complete Genome |
| Technology | Illumina | PacBio |
| Contigs | NZ_JXRW01000001–NZ_JXRW01000002 | LR881168.1 (chr) + LR881169.1 (plasmid) |
| Locus tag prefix | TK37_RS* (PGAP) | ALT831_RS* (PGAP), MIT1002_NNNNN (GCA old_locus_tag) |
| Gene count | ~4,159 CDS | 4,028 annotated genes |
| BioProject | PRJNA277277 | PRJEB32626 |

The KG uses `GCF_901457835.2` as the canonical assembly. Gene nodes use `MIT1002_NNNNN` (5-digit) locus tags.

## Papers and Their ID Systems

### Biller 2016 (supp_table_3) — Alteromonas DE in co-culture
- **Gene ID column**: `MIT1002_NNNN` (4-digit, RAST annotation on draft genome)
- **Also has**: `Genomic Region ID` = RAST coordinate format `contigNNNNN_start_stop`
- **766 DE genes, 2 timepoints (24h vs 12h, 48h vs 12h)**
- **Current resolution**: 93.6% via heuristic zero-padding (MIT1002_NNNN → MIT1002_0NNNN) — **WRONG** (3.8% product match between mapped pairs)
- **Needs**: verified RAST → TK37_RS → WP_ → MIT1002_NNNNN bridge

### Biller 2018 (supp_table_s6b) — Alteromonas DE in extended darkness
- **Gene ID column**: `RAST_region_ID` = coordinate format `contig00001_17088_18476`
- **180 DE genes, 2 timepoints (1h, 5h)**
- **Current resolution**: 0% (coordinate format not in any lookup)
- **Needs**: same bridge as Biller 2016 (coordinates are on JXRW draft contigs)

### Coe 2024 (supp_table_4) — Alteromonas DE in diel cycle
- **Gene ID column**: `NCBI ID` = `MIT1002_NNNNN` (5-digit, current canonical)
- **Also has**: `Gene ID` = `cds-TK37_RS*` (old draft assembly)
- **3876 genes, 99.5% resolved** — no issues
- **Bonus**: provides TK37_RS → MIT1002_NNNNN mapping for 3876 genes

### Biller 2018 (supp_table_s3) — Prochlorococcus NATL2A DE
- **Uses NATL2A gene IDs** — not MIT1002 genes. 99.9% resolved. No action needed.

## Why Zero-Padding is Wrong

The 4-digit `MIT1002_NNNN` IDs (from RAST annotation on the draft genome) and the 5-digit `MIT1002_NNNNN` IDs (from IMG/NCBI on the complete genome) come from **independent gene-calling pipelines on different assemblies**. Zero-padding `MIT1002_0001` → `MIT1002_00001` maps to a completely different gene:

- Gene length comparison: wildly mismatched (e.g., 122bp draft vs 1592bp current)
- Product description match: only 3.8% (15/391 compared)
- The old `_with_locus_tag.csv` workaround in the repository used the same wrong zero-padding logic

## Available Resources

### Conversion table (from Biller 2018 supplement)
`MIT1002_systematicnames_conversiontable.csv` (4136 rows):
- `original_rast_id`: fig|226.6.peg.N (RAST internal)
- `genbank_id`: MIT1002_NNNN (4-digit, same as Biller 2016 Gene ID)
- `assembly_coordinates`: contigNNNNN_start_stop (same coordinate system as RAST_region_ID in Biller 2018)

This table bridges RAST IDs ↔ coordinates but does NOT contain any TK37_RS or MIT1002_NNNNN IDs.

### Old assembly GFF (downloaded)
`cache/data/Alteromonas/genomes/MIT1002/genomic_old_draft_GCF_001077695.1.gff`
- Has TK37_RS* locus tags + WP_* protein IDs on NZ_JXRW01000001-2 contigs

### Coe 2024 mapping
The Coe 2024 supplementary table maps `cds-TK37_RS*` → `MIT1002_NNNNN` for 3876 genes.

### Shared WP_ proteins
3,892 WP_ RefSeq protein accessions are shared between the two assemblies, providing a natural bridge.

## Proposed Mapping Chain

```
RAST MIT1002_NNNN (4-digit)
  → (conversion table: coordinate match) →
TK37_RS* (old draft PGAP locus tags)
  → (WP_ protein bridge, 3892 shared) →
ALT831_RS* / MIT1002_NNNNN (current canonical)
```

### Step 1: RAST → TK37_RS (coordinate matching)
Match conversion table coordinates (contigNNNNN_start_stop) against old GFF gene positions. Requires contig name correspondence: `contig00001` ↔ `NZ_JXRW01000001.1` etc.

### Step 2: TK37_RS → MIT1002_NNNNN (WP_ protein bridge)
Use shared WP_ protein accessions between old GFF and current GFF. OR use the Coe 2024 table directly (has cds-TK37_RS → MIT1002_NNNNN for 3876 genes).

## What We Requested from Authors

Priority order:
1. **Direct mapping table**: `MIT1002_NNNN` (4-digit RAST) → `MIT1002_NNNNN` (5-digit canonical)
2. **RAST protein FASTA**: with MIT1002_NNNN headers — enables direct sequence matching against current protein.faa
3. **Confirmation**: that the conversion table `assembly_coordinates` are on the JXRW draft contigs

## Heuristic Zero-Padding — RESOLVED

The heuristic zero-padding issue is now moot: with the correct mapping in `specific_lookup`, `MIT1002_NNNN` (4-digit) IDs are found directly at Tier 1 before the heuristic runs. The heuristic remains available for other strains where it's useful.

## Solution (2026-03-05)

Author (Steve Biller) provided RAST protein FASTA (`226.6.faa`, 4214 proteins with `fig|226.6.peg.N` headers) and RAST GBK/GFF files.

### Approach: Diamond protein matching + transitive closure

1. **Diamond matching**: `scripts/map_img_to_ncbi_proteins.py` mapped RAST proteins (fig| headers) against NCBI MIT1002 proteins:
   - Phase 1 (exact): 3347 (79.4%)
   - Phase 2 (subsequence): +155 (83.1%)
   - Phase 3 (Diamond blastp): +389 (92.3%, 3891/4214)
   - Fragment dedup: 41 shorter fragments discarded
   - Output: `rast_to_mit1002_id_translation.csv` (rast_fig_id → locus_tag)

2. **Two id_translation entries in Biller 2018 paperconfig**:
   - Diamond output: links fig| IDs → locus_tags (old_locus_tag Tier 1)
   - Conversion table: links fig| IDs → MIT1002_NNNN + coordinates (alternative_locus_tag Tier 1)

3. **Transitive closure**: `build_gene_id_mapping.py` converged in 3 passes:
   - Pass 1: fig| → locus_tag (from diamond)
   - Pass 2-3: MIT1002_NNNN and coordinates linked to same locus_tag via shared fig| IDs
   - specific_lookup: 9,637 → 38,029 entries

4. **Resolution results**:
   - Biller 2016: **96.6%** (740/766) — up from ~9% wrong
   - Biller 2018 S6B: **97.2%** (175/180) — up from 0%
   - 26+5 unresolved = RAST-specific ORFs with no NCBI counterpart

### Files added/modified

- `data/.../Biller 2018/MIT1002_RAST_annotation/226.6.faa` — RAST protein FASTA from author
- `data/.../Biller 2018/MIT1002_RAST_annotation/226.6.gbk` — RAST GenBank from author
- `data/.../Biller 2018/MIT1002_RAST_annotation/Alteromonas_A1A_concat.gff` — RAST GFF from author
- `data/.../Biller 2018/MIT1002_RAST_annotation/rast_to_mit1002_id_translation.csv` — diamond output
- `data/.../Biller 2018/paperconfig.yaml` — added 2 id_translation entries, changed S6B id_type
- `data/.../biller 2016/paperconfig.yaml` — added coordinate column to id_columns
- `scripts/map_img_to_ncbi_proteins.py` — generalized (--source-id-col, temp prefix)
- `multiomics_kg/download/build_gene_id_mapping.py` — BOM fix (encoding='utf-8-sig')
- `cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json` — rebuilt with all RAST IDs
