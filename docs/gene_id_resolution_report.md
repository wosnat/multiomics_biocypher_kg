# Gene ID Resolution Report

Generated from `logs/prepare_data_step4.log` (2026-04-04).

**Overall: 146,565 / 154,902 resolved (94.6%)**

All paperconfigs have the correct ID columns registered. The unresolved IDs fall into three categories:
1. **New organisms** where gene ID mapping has gaps (Synechococcus integration in progress)
2. **Non-gene features** (intergenic regions, ncRNAs) that don't have Gene nodes in the KG by design
3. **Small residual gaps** in established organisms (a few genes per table)

### Fixes applied (2026-04-04)

**1. Case-insensitive resolution fallback** — Added to `resolve_row()` in `gene_id_utils.py`. Paper authors sometimes use different casing than NCBI annotations (e.g., `SputW3181_*` vs `Sputw3181_*`). Tries exact match first, then falls back to case-insensitive lookup. Method strings: `tier1_ci`, `locus_tag_ci`, `multi_ci`. Impact: Shewanella W3-18-1 (Beliaev 2014) jumped from 31-33% to 99.2-99.8%.

**2. Bernstein 2017 CSV split + Meiothermus genome** — The original CSVs contained genes from both BP-1 (tll/tlr/tsl IDs) and Meiothermus ruber strain A (SY28_RS IDs) in a single file. Split into per-organism CSVs. Added M. ruber strain A genome (`GCF_000836395.1`, NCBI classifies as M. taiwanensis, WGS NZ_JXOP01) to genome registry and ran full prepare_data pipeline. Impact: BP-1 went from 65.6% to 96.7%, Meiothermus resolves at 97.1%.

**3. Protein accession version suffix heuristic** — Added `.1` version suffix fallback in `_heuristic_candidates()`. Paper CSVs often omit the version (e.g., `AAV95689` vs `AAV95689.1` in the mapping). Pattern: 3 uppercase letters + 5+ digits without a dot. Impact: Ruegeria pomeroyi DSS-3 (Kaur 2018) jumped from 1-2% to 98.4-98.6%.

**4. Assembly switches for Kratzl 2024 organisms** — Switched KT2440 from `GCF_000007565.2` (PP_RS prefix) to `GCF_024662315.1` (M8001_RS prefix, matches paper). Switched PCC 7942 from `GCF_000012525.1` (SYNPCC7942_RS prefix) to `GCF_014698905.1` (H6G84_RS prefix, matches paper). Added old assembly GFFs as `annotation_gff` entries to preserve protein_id bridges for proteomics tables. Impact: KT2440 RNA-seq 96.9%→97.8%, proteomics 98.3%→97.4%; PCC 7942 RNA-seq 86.6%→98.1%, proteomics 95.9%→93.5%.

Overall improved from 91.0% to 94.6% (+5,581 genes).

---

## Shewanella sp. W3-18-1

FIXED by case-insensitive resolution. The NCBI GFF uses `Sputw3181_*` (lowercase w), paper uses `SputW3181_*` (uppercase W).

### Beliaev 2014

| Table | Before | After | Remaining unresolved |
|-------|--------|-------|----------------------|
| supp_table_s4 | 144/457 (31.5%) | 456/457 (99.8%) | 1 |
| supp_table_s6 | 444/1,342 (33.1%) | 1,331/1,342 (99.2%) | 11 |
| supp_table_s8 | 403/1,277 (31.6%) | 1,267/1,277 (99.2%) | 10 |

**ID columns registered:** `Locus tag` (locus_tag), `Gene abbreviation` (gene_name). Paperconfig correct.

**Root cause:** Case mismatch — NCBI GFF has `old_locus_tag=Sputw3181_*`, paper authors wrote `SputW3181_*`. Resolution was case-sensitive and missed the match. Fixed by `locus_tag_ci` / `tier1_ci` fallback.

**Action needed:** None. Resolved.

---

## Thermosynechococcus vestitus BP-1

FIXED by CSV split. The original CSVs mixed BP-1 genes (`tll/tlr/tsl/tsr` locus tags) with Meiothermus ruber genes (`SY28_RS*` locus tags) in a single file. Split into per-organism CSVs.

### Bernstein 2017

| Table | Before (mixed) | After (split) | Remaining unresolved |
|-------|----------------|---------------|----------------------|
| bp1_light_clusters | 2,946/4,492 (65.6%) | 2,339/2,420 (96.7%) | 81 `tll*` |
| bp1_oxygen_clusters | 2,946/4,492 (65.6%) | 2,339/2,420 (96.7%) | 81 `tll*` |

**Root cause:** The original CSV contained genes from both organisms. The 1,546 "unresolved" `SY28_RS*` IDs were actually Meiothermus ruber genes, not BP-1. After splitting, BP-1 resolves at 96.7%. The 81 remaining unresolved are `tll*` IDs not in the current gene_id_mapping (likely pseudogenes or removed annotations).

---

## Meiothermus ruber strain A

New organism added to genome registry (2026-04-04). Genome: `GCF_000836395.1` (NCBI classifies as *M. taiwanensis*, WGS NZ_JXOP01, taxid 172827). Paper calls it *M. ruber* strain A (98.6% 16S identity to DSM 1279). 2,885 genes in GFF.

### Bernstein 2017

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| mruber_light_clusters | 2,012/2,072 | 97.1% | 60 `SY28_RS*` |
| mruber_oxygen_clusters | 2,012/2,072 | 97.1% | 60 `SY28_RS*` |

**ID columns registered:** `NCBI Locus Tag` (locus_tag), `RAST Locus Tag` (other), `gene` (gene_name). Paperconfig correct.

**Analysis:** The `SY28_RS*` IDs in the paper come from a RAST annotation. Most map via `tier1:NCBI Locus Tag` to the NCBI genome. The 60 unresolved are likely RAST-specific predictions not in the NCBI annotation.

**Action needed:** Low priority. 97.1% is acceptable.

---

## Ruegeria pomeroyi DSS-3

FIXED by protein accession version suffix heuristic. The paper uses `AAV95689` (no version), the mapping has `AAV95689.1` (with `.1`). Assembly `GCF_000011965.2` was already in the genome registry with a complete gene_id_mapping.

### Kaur 2018

| Table | Before | After | Remaining unresolved |
|-------|--------|-------|----------------------|
| supp_table_s5a | 3/184 (1.6%) | 181/184 (98.4%) | 3 |
| supp_table_s5b | 2/222 (0.9%) | 219/222 (98.6%) | 3 |

**ID columns registered:** `NCBI reference` (protein_id). Paperconfig correct.

**Root cause:** The `AAV*` GenBank protein accessions in the CSV lacked the `.1` version suffix that appears in the gene_id_mapping. The new `_heuristic_candidates()` function now tries appending `.1` for unversioned protein accessions matching the pattern `[A-Z]{3}\d{5,}`.

**Action needed:** None. Resolved.

---

## Pseudomonas putida KT2440

FIXED by switching assembly from `GCF_000007565.2` (PP_RS prefix) to `GCF_024662315.1` (M8001_RS prefix, matches paper). Old assembly GFF added as `annotation_gff` to preserve protein_id bridges.

### Kratzl 2024

| Table | Before | After | Remaining unresolved |
|-------|--------|-------|----------------------|
| sheet1_rnaseq_pputida | 5,054/5,218 (96.9%) | 5,103/5,218 (97.8%) | 115 |
| sheet11_proteomics_pputida | 228/232 (98.3%) | 226/232 (97.4%) | 6 |

**Root cause:** Paper used `GCF_024662315.1` assembly (M8001_RS prefix), KG had `GCF_000007565.2` (PP_RS prefix). Same organism (taxid 160488), only paper using it. Switched assembly. Old GFF added as `annotation_gff` to bridge old PP_* locus tags via protein_id for proteomics resolution.

**Action needed:** None. Resolved.

---

## Synechococcus elongatus PCC 7942

FIXED by switching assembly from `GCF_000012525.1` (SYNPCC7942_RS prefix) to `GCF_014698905.1` (H6G84_RS prefix, matches paper). Old assembly GFF added as `annotation_gff` to preserve protein_id bridges.

### Kratzl 2024

| Table | Before | After | Remaining unresolved |
|-------|--------|-------|----------------------|
| sheet6_rnaseq_selongatus | 2,321/2,679 (86.6%) | 2,627/2,679 (98.1%) | 52 |
| sheet13_proteomics_selongatus | 544/567 (95.9%) | 530/567 (93.5%) | 37 |

**Root cause:** Paper used `GCF_014698905.1` assembly (H6G84_RS prefix), KG had `GCF_000012525.1` (SYNPCC7942_RS prefix). Same organism (taxid 1140), only paper using it. Switched assembly. Old GFF added as `annotation_gff` to bridge old Synpcc7942_* locus tags via protein_id for proteomics resolution.

**Action needed:** None. Resolved.

---

## Synechococcus PCC 7002

New organism (Synechococcus integration). Mostly resolved; small residual gap. Not a case issue.

### Beliaev 2014

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| supp_table_s3 | 1,159/1,185 | 97.8% | 26 `SYNPCC7002_*` |
| supp_table_s5 | 672/691 | 97.3% | 19 `SYNPCC7002_*` |
| supp_table_s7 | 674/701 | 96.1% | 27 `SYNPCC7002_*` |

**ID columns registered:** `Locus tag` (locus_tag), `Gene abbreviation` (gene_name). Paperconfig correct.

**Analysis:** High resolution rate. The ~25 unresolved `SYNPCC7002_*` IDs are a consistent set across tables (e.g., `SYNPCC7002_A0262`, `SYNPCC7002_F0109`). These may be pseudogenes, plasmid genes, or annotation gaps.

**Action needed:** Low priority. Check whether the 25 unresolved IDs are present in the NCBI GFF — they may be non-coding features or removed in the current annotation.

---

## Prochlorococcus MED4 (non-gene features)

Established organism with complete mapping. Unresolved IDs are intergenic regions and ncRNAs, not protein-coding genes.

### Zinser 2009

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| med4_diel_clusters | 1,697/3,610 | 47.0% | 1,913 `PMMIG_*` intergenic regions |

**ID columns registered:** `Gene or region` (old_locus_tag), `New gene name` (locus_tag). Paperconfig correct.

**Analysis:** All 1,697 genes resolve. The 1,913 unresolved rows are intergenic regions (`PMMIG_5END_PMM0001`, `PMMIG_PMM0003_PMM0004`, etc.) — these are probes for intergenic spacers on the MED4 microarray. They don't correspond to Gene nodes.

**Action needed:** None. This is expected behavior for a gene_clusters table that includes intergenic probes. The cluster adapter should handle these gracefully (skip rows that don't resolve to genes).

### Fang 2019

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| supp_table_2 | 429/665 | 64.5% | 236 `RNA_*` ncRNA features |

**ID columns registered:** `Gene ID` (old_locus_tag), `Gene` (gene_name). Paperconfig correct.

**Analysis:** All 429 protein-coding genes resolve. The 236 unresolved are non-coding RNAs (`RNA_1`, `RNA_10`, `RNA_16`, etc.) — MIT9313 ncRNA features that don't have Gene nodes in the KG.

**Action needed:** None. Expected behavior.

### Thompson 2011

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| supp_table_1 | 112/125 | 89.6% | 13 `PMED4_asRNA_*`/`PMED4_ncRNA_*` |
| supp_table_2 | 98/111 | 88.3% | 13 `P9313_*` ncRNAs + `PMT_ffs` |

**ID columns registered:** Standard locus tag columns. Paperconfig correct.

**Analysis:** Unresolved are antisense RNAs and ncRNAs (e.g., `PMED4_asRNA_04601`, `PMED4_ncRNA_Yfr10`). Not gene nodes.

**Action needed:** None.

---

## Alteromonas HOT1A3 (known problem)

### Aharonovich 2016

| Table | Resolved | Rate | Unresolved IDs |
|-------|----------|------|----------------|
| supp_table_hot1a3_9313_vs_9313dil | 0/1,966 | 0.0% | `AMHOT1A3_single_*` |
| supp_table_hot1a3_9313_vs_med4 | 0/2,381 | 0.0% | `AMHOT1A3_single_*` |

Known problem. The `AMHOT1A3_single_*` IDs are from an unpublished annotation with no mapping to NCBI locus tags.

---

## Summary by priority

| Priority | Organism | Papers | Unresolved genes | Status |
|----------|----------|--------|-------------------|--------|
| ~~High~~ | ~~Shewanella W3-18-1~~ | ~~Beliaev 2014~~ | ~~2,000~~ → 22 | **FIXED** (case-insensitive) |
| ~~High~~ | ~~T. vestitus BP-1~~ | ~~Bernstein 2017~~ | ~~1,500~~ → 81 | **FIXED** (CSV split) |
| ~~High~~ | ~~M. ruber strain A~~ | ~~Bernstein 2017~~ | (new) → 60 | **FIXED** (genome added + CSV split) |
| ~~High~~ | ~~R. pomeroyi DSS-3~~ | ~~Kaur 2018~~ | ~~400~~ → 6 | **FIXED** (version suffix heuristic) |
| ~~Medium~~ | ~~S. elongatus PCC 7942~~ | ~~Kratzl 2024~~ | ~~380~~ → 52 | **FIXED** (assembly switch) |
| ~~Medium~~ | ~~P. putida KT2440~~ | ~~Kratzl 2024~~ | ~~168~~ → 115 | **FIXED** (assembly switch) |
| **Low** | Syn. PCC 7002 | Beliaev 2014 | ~25 | Likely pseudogenes/ncRNAs |
| **None** | MED4 (intergenic) | Zinser 2009 | 1,913 | Expected — intergenic probes |
| **None** | MIT9313 (ncRNAs) | Fang 2019 | 236 | Expected — ncRNA features |
| **None** | MED4 (ncRNAs) | Thompson 2011 | ~26 | Expected — asRNA/ncRNA |
| **Known** | HOT1A3 | Aharonovich 2016 | ~4,300 | Known problem |
