# MIT1002 Deploy Plan — Cross-Assembly Gene ID Mapping

## Status: BLOCKED — waiting for author response

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

## Critical Bug: Heuristic Zero-Padding

The `_heuristic_candidates()` function in `gene_id_utils.py` auto-pads `MIT1002_NNNN` → `MIT1002_0NNNN` and finds a hit in `specific_lookup`. This produces silently wrong mappings. Options:
1. **Disable heuristic for MIT1002**: Add guard in `_heuristic_candidates` to skip when prefix length mismatch exceeds 1 digit
2. **Remove Biller 2016/2018 MIT1002 tables from resolution** until bridge is ready
3. **Mark as known-wrong** and wait for author bridge

Current choice: option 3 — document the issue, leave the heuristic as-is (it's useful for other strains), and plan to replace with a correct bridge once author responds.

## Files Modified in This Investigation

- `data/Prochlorococcus/papers_and_supp/biller 2016/paperconfig.yaml` — reverted from `_with_locus_tag.csv` workaround to original CSV with `name_col: Gene ID`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml` — removed incorrect `id_translation_mit1002` entry
- Deleted: `data/Prochlorococcus/papers_and_supp/Biller 2018/MIT1002_id_bridge.csv` (was wrong)
- Saved: `cache/data/Alteromonas/genomes/MIT1002/genomic_old_draft_GCF_001077695.1.gff`
- Rebuilt: `gene_id_mapping.json` (without incorrect bridge)

## Next Steps (When Author Responds)

1. Build verified RAST → TK37_RS coordinate bridge (or use direct mapping if provided)
2. Chain through WP_ proteins to get RAST → MIT1002_NNNNN
3. Add as `id_translation` to Biller 2016 and Biller 2018 paperconfigs
4. Consider adding old assembly GFF as `annotation_gff` entry to auto-include TK37_RS → WP_ mappings
5. Rebuild gene_id_mapping.json, re-resolve, verify product match rates
6. Rebuild KG + snapshot comparison
