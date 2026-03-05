# Gene ID Mapping v2 — Implementation Status

## Plan Reference
Full approved plan: `/home/osnat/.claude/plans/sleepy-splashing-kettle.md`

## What Was Implemented

### Core files (all complete)

| File | Status | Notes |
|------|--------|-------|
| `multiomics_kg/download/gene_id_graph.py` | ✅ DONE | GeneIdGraph class, three-tier classification, iterative convergence, normalize_id, build_diagnostic_report, to_json_structure |
| `multiomics_kg/download/build_gene_id_mapping.py` | ✅ DONE | Rewritten to use GeneIdGraph; collects all source rows, process_all_rows(), writes v2 JSON + diagnostic report; NO gene_mapping_supp.csv |
| `multiomics_kg/utils/gene_id_utils.py` | ✅ DONE | Added MappingData, load_mapping_v2(), expand_list(), _heuristic_candidates(), resolve_row() |
| `multiomics_kg/download/resolve_paper_ids.py` | ✅ DONE | Uses load_mapping_v2() + resolve_row(), multi-column fallback, writes _resolved_report.txt |
| `tests/test_gene_id_graph.py` | ✅ DONE | 40 unit tests, all passing (897 total tests pass) |
| `.claude/skills/fix-gene-ids/fix_gene_ids.py` | ✅ DONE | Updated to load_mapping_v2() + resolve_row(); fallback to legacy |
| `.claude/skills/fix-gene-ids/SKILL.md` | ✅ DONE | Updated description |
| `.claude/skills/check-gene-ids/check_gene_ids.py` | ✅ DONE | Added check_gene_id_mapping_v2(), new USE_STEP4 / AMBIGUOUS_IN_MAPPING fix strategies |
| `docs/methods_gene_id_mapping.md` | ✅ DONE | ~500-word scientific methods section |
| `CLAUDE.md` | ✅ DONE | Gene ID Mapping section updated with tier table, v2 schema, method strings |

## What Remains

### Phase 1 deployment: MIT9312 ✅ DONE (updated 2026-03-03)
Papers: barreto 2022 (100%), tetu 2019 (100%), Hennon 2017 (**103/103 100%** — was deferred, now resolved).

**Hennon 2017 fix**: Sonya provided `Prochlorococcus_marinus_MIT9312_annot.txt` mapping `lcl|CP000111.1_cds_*` CDS accession IDs → `PMT9312_NNNN` old locus tags. Added as `id_translation_mit9312_annot` to paperconfig with `cds_fna_id` + `old_locus_tag` columns. Also fixed resolver to not skip tables where `name_col == "locus_tag"` but id_columns have non-standard id_types (the CSV column is literally named "locus_tag" but contains CDS accession IDs). Resolver now renames the original column to `original_id` when this collision occurs.

- Mapping stats: 1980 genes, 11593 specific_lookup, 10091 multi_lookup, 2 conflicts, 3 passes

### Phase 2 deployment: MIT9301 ✅ DONE (2026-03-03)

Paper: **Anjur-Dietrich 2025** — `de_genes_freeliving_vs_biofilms.csv` (242 rows, JGI IDs)
- Resolution: **235/242 (97.1%)** via `tier1:ID`
- 7 unresolved: genes in annotation table with null `uniprot_gene_name` AND null `uniprot_entry_name` — no locus tag bridge possible
- Snapshot: `before_MIT9301` (110,028 edges total)
- KG rebuild: deferred (docker skipped)

**Key fix applied**: `annotation_genome_9301` `uniprot_gene_name` column changed from `id_type: gene_name` → `id_type: old_locus_tag` in `data/Prochlorococcus/papers_and_supp/Anjur 2025/paperconfig.yaml`. Values are old locus tags (`P9301_XXXXX`) — bare single tokens need Tier 1 Phase 1 lookup; compound `"gene_name P9301_XXXXX"` values work via Phase 2 whitespace-split.

### Phase 3 deployment: NATL1A ✅ DONE (2026-03-03)

Paper: **He 2022** — `table s1 ... NATL1A .csv` (82 rows, mixed locus tags + gene names)
- Resolution: **81/82 (98.8%)** — 1 empty row
- Methods: `locus_tag:Gene Name`=60, `tier1:Gene Name`=21
- Previous (v1 fix-gene-ids): 75/82 (91.5%) — 6 gene names (`wza`, `cyoA`, `cyoE`, `afuA`, `tal`, `dedA`) unresolvable without v2 mapping

**Key fix applied**: Reverted paperconfig to use original CSV (without `_with_locus_tag` suffix) and `name_col: "Gene Name"`. Removed `id_columns` referencing non-existent `locus_tag` column. The v2 resolver handles the mixed locus-tag/gene-name `Gene Name` column natively — locus tags resolve via `locus_tag:<col>` pass, gene names via `tier1:<col>` (specific_lookup).

**Also rebuilt MED4 v2 mapping** (He 2022 supp_table_s2 covers MED4): 27/31 (87.1%), 4 pre-unresolved (3 ncRNAs + 1 empty).

- KG rebuild: deferred (docker skipped)

### Phase 4 deployment: NATL2A ✅ DONE (2026-03-03)

Papers: **Biller 2016**, **Biller 2018**, **Coe 2024**, **Lin 2015**, **Tetu 2019** — 5 papers, ~37 RNASEQ analyses

Resolution rates (all using standard locus tags — no paperconfig changes needed):
| Paper | Table | Rate | Method |
|-------|-------|------|--------|
| Biller 2016 | supp_table_2 | **353/353 (100%)** | `locus_tag:Original_NCBI_ID` |
| Biller 2018 | supp_table_s3 | **1955/1957 (99.9%)** | `locus_tag:NCBI ID_3` + `tier1:NCBI ID_3` |
| Coe 2024 | supp_table_2 | **2029/2030 (100%)** | `locus_tag:NCBI ID_3` + `tier1:NCBI ID_3` |
| Lin 2015 | supp_table_s4a | **34/34 (100%)** | `locus_tag:ID` |
| Lin 2015 | supp_table_s4b | **34/34 (100%)** | `locus_tag:ID` |
| Tetu 2019 | supp_dataset_4 | **171/171 (100%)** | `locus_tag:NATL2A original locus tag` |

- 3 unresolved total: PMN2A_1906 (not in NCBI annotations — 2 papers) + "(unannotated)" empty row
- Mapping stats: 2214 genes, 10430 specific_lookup, 12503 multi_lookup, 20 conflicts (all `alternative_locus_tag`/`locus_tag_cyanorak` gene family collisions)
- Snapshot: `before_NATL2A` (110,028 edges total)
- KG rebuild: deferred (docker skipped)
- No paperconfig fixes needed — straightforward deployment

### Phase 5 deployment: MIT9313 ✅ DONE (2026-03-03)

Papers: **Aharonovich 2016**, **Fang 2019**, **Martiny 2006**, **Thompson 2011**, **Tolonen 2006** — 7 supplementary tables

MIT9313 was heavily reannotated (4 annotation generations: PMT0001 → PMT_0001 → RG24_RS00005 → AKG35_RS00005). Required three code/data fixes beyond standard deployment.

**Resolution rates:**

| Paper | Table | Rate | Notes |
|-------|-------|------|-------|
| Aharonovich 2016 | supp_table_2 | **2268/2269 (100.0%)** | 1 unresolved: PMT1676 |
| Aharonovich 2016 | supp_table_3 | **2267/2268 (100.0%)** | 1 unresolved: PMT1676 |
| Tolonen 2006 | supp_table_mit9313 | **2241/2250 (99.6%)** | 9 unresolved: PMT_or* + a few PMT0### not in annotations |
| Martiny 2006 | supp_table_2 | **177/181 (97.8%)** | 4 unresolved: PMT1210 + table footnotes |
| Thompson 2011 | supp_table_2 | **95/111 (85.6%)** | See below |
| Fang 2019 | supp_table_2 | **429/665 (64.5%)** | 236 unresolved = `RNA_*` non-coding IDs (expected) |

**Mapping stats:** 2958 genes, 18188 specific_lookup, 17452 multi_lookup, 44 conflicts, 3 passes

**Fix 1 — Position-based merge in gene_mapping.csv** (`build_gene_mapping.py`):
MIT9313 GCF GFF `old_locus_tag` sometimes omits the PMT0### (Cyanorak) form, causing 360 genes to split into separate Cyanorak-only and NCBI-only rows in gene_mapping.csv. A `_position_fallback_merge()` function was added that matches unmatched Cyanorak/NCBI entries by genomic coordinate overlap (≥90% reciprocal overlap, ≤10bp boundary diff). Merged 105 genes. Also added preservation of the consumed Cyanorak locus_tag in `old_locus_tags` so it remains discoverable by gene_id_mapping.

**Fix 2 — GCA GFF `Alternative locus ID` parsing** (`build_gene_id_mapping.py`):
GCA GFF `Note` attribute contains `Alternative locus ID: P9313_NNNNN` — ProPortal-era locus tags used by Thompson 2011 microarrays (GEO platform GPL11412). These are a DIFFERENT numbering scheme from the MIT9313_genbank.tsv P9313_ column (which uses steps of 10). Added regex extraction from the Note field as `alternative_locus_tag` (Tier 1). Also benefits MED4, MIT9312, and NATL2A (all have Alternative locus IDs in their GCA GFF: PMED4_, P9312_, NATL2_ prefixes).

**Fix 3 — Thompson 2011 CSV data cleanup** (`supp table 2.csv`):
47 entries had `P9313 NNNNN` (space instead of underscore) — formatting error from the original supplementary table. Fixed to `P9313_NNNNN`. Also added `+` to the heuristic footnote-stripping in `gene_id_utils.py` (`_heuristic_candidates` now strips trailing `*` and `+`).

**Thompson 2011 unresolved (16):**
- 7 tRNA/ncRNA entries (expected — non-coding, not in gene annotations)
- 5 P9313_ IDs mapping to PMT_ locus_tags absent from gene_mapping.csv (GCA-only genes with no NCBI/Cyanorak entry)
- 2 reversed format (`17651 P9313`, `18001 P9313`)
- 1 `PMT ffs` (sRNA, not a gene)
- 1 `ncRNA Yfr2-5 1`

**MIT9313_resources id_translation** (`MIT9313_genbank.tsv`):
Translation table mapping AKG35_RS → PMT_ → P9313_ (GenBank numbering) → UniProt. Added as `data/Prochlorococcus/papers_and_supp/MIT9313_resources/paperconfig.yaml`. Does NOT have the PMT0### (Cyanorak, no underscore) or the ProPortal P9313_ numbering — these come from GCF GFF old_locus_tag and GCA GFF Note field respectively.

- KG rebuild: deferred (docker skipped)

### Phase 6 deployment: WH8102 ✅ DONE (2026-03-03)

Paper: **Barreto 2022** — `table_S1_WH8102.csv` (13 rows, mixed SYNW#### locus tags + gene names)

Resolution rates:
| Paper | Table | Rate | Notes |
|-------|-------|------|-------|
| Barreto 2022 | supp_table_S1_WH8102 | **12/13 (92.3%)** | 1 unresolved: RNA_15 (tRNA-Met, expected) |

- Mapping stats: 3004 genes, 11248 specific_lookup, 13866 multi_lookup, 26 conflicts, 2 passes
- Methods: `locus_tag:symbol`=6, `tier1:symbol`=6
- No paperconfig fixes needed — straightforward deployment
- `id_translation_syn_8102` from Barreto 2022 provides SYNW→gene_name bridges
- GCF accession verified: **GCF_000195975.1** (Parasynechococcus marenigrum WH 8102, taxid 84588)
- KG rebuild: deferred (docker skipped)

### Phase 7 deployment: CC9311 ✅ DONE (2026-03-03)

Paper: **Barreto 2022** — `table_S1_CC9311.csv` (30 rows, `sync_NNNN` locus tags in `symbol` column)

Resolution rates:
| Paper | Table | Rate | Notes |
|-------|-------|------|-------|
| Barreto 2022 | supp_table_S1_CC9311 | **27/30 (90.0%)** | 3 unresolved: sync_2538 (tRNA-Ile_1), sync_0554 (tRNA-Ile_2), sync_2844 (tRNA-Val) — all tRNAs, expected |

- Mapping stats: 3082 genes, 11766 specific_lookup, 14877 multi_lookup, 0 conflicts, 2 passes
- Methods: `locus_tag:symbol`=24, `tier1:symbol`=3
- No paperconfig fixes needed — straightforward deployment
- `id_translation_syn_9311` from Barreto 2022 provides sync_→gene_name→uniprot bridges
- GCF accession verified: **GCF_000014585.1** (Synechococcus sp. CC9311, taxid 64471)
- KG rebuild: deferred (docker skipped)

### Phase 8 deployment: RSP50 — BLOCKED (2026-03-03)

Paper: **Labban 2022** — `supp table s1.csv` (267 rows × 4 temperature comparisons, `BSR22_NNNNN` IDs)

**Problem**: Paper uses custom annotation with sequential BSR22_ gene IDs (BSR22_00001, BSR22_00002...) that don't match NCBI PGAP annotation (BSR22_00005, BSR22_00010... stepping by 5). Both GCF and GCA GFF have the PGAP numbering. The 46 IDs that numerically coincide are **false matches** — annotations don't correspond (e.g., paper's BSR22_00120="Flavodoxin" vs NCBI BSR22_00120="MauE/DoxX family protein").

**Root cause**: Authors used featureCounts with a custom/earlier annotation GFF against chromosome CP018344. The original submitter annotation was replaced by NCBI PGAP. No mapping file or original GFF available.

**Action taken**:
- Commented out Labban 2022 in `paperconfig_files.txt` (prevents false edges)
- Authors contacted 2026-03-03 for their annotation GFF
- v2 `gene_id_mapping.json` rebuilt for RSP50 (1870 genes, 5415 specific, 7261 multi, 0 conflicts)
- RSP50 has no other papers — strain has zero expression edges until Labban resolved

**TODO**: When authors provide their GFF, add it as `annotation_gff` entry in Labban 2022 paperconfig, uncomment in paperconfig_files.txt, rebuild mapping + resolve.

### Phase 9 deployment: EZ55 ✅ DONE (2026-03-04)

Papers: **Barreto 2022** (7 EZ55 tables), **Hennon 2017** (2 EZ55 tables)

**Barreto 2022**: Uses `EZ55_NNNNN` (5-digit) from GCA_901457815.2 annotation. Resolves at **~95%** via `locus_tag:symbol` pass. No fixes needed.

**Hennon 2017**: Uses `AEZ55_NNNN` (4-digit) from researcher's own annotation of IMG draft genome 2785510739 (JCVI draft, 3 scaffolds, 5016 genes). Author (Gwenn Hennon) provided protein FASTA (`EZ55_annotation/ez55_aa.fasta`, 4930 proteins) + GenBank flat files (`ez55-{1,2,3}.gbf`).

**Cross-assembly protein sequence bridging** via `scripts/map_img_to_ncbi_proteins.py`:

| Phase | Method | Matches | Cumulative |
|-------|--------|---------|------------|
| 1 | Exact protein sequence match | 2,971 | 2,971 (60.3%) |
| 2 | Subsequence match (≥95% overlap) | 136 | 3,107 (63.0%) |
| 3 | Diamond blastp (≥80% id, ≥60% qcov) | 946 | 4,053 (82.2%) |

Fragment deduplication: 496 shorter fragments discarded (draft frameshifts split canonical genes into multiple ORFs; longest fragment kept for expression edges).

**Resolution rates:**

| Paper | Table | Rate | Notes |
|-------|-------|------|-------|
| Hennon 2017 | supp_table_3 | **345/422 (81.8%)** | `tier1:locus tag` via `old_locus_tag` bridge |
| Hennon 2017 | supp_table_4 | **98/125 (78.4%)** | `tier1:locus tag` via `old_locus_tag` bridge |

Unresolved (~20%): 47 draft-specific ORFs (no canonical counterpart), 49 discarded fragments, remainder below Diamond thresholds.

**Mapping stats**: specific_lookup grew from 10,036 → 22,322. 0 conflicts. Converged in 2 passes.

**Key lesson**: `id_translation` entries must declare BOTH the anchor column (`locus_tag` with `id_type: locus_tag`) AND the new mapping column (`aez55_id` with `id_type: old_locus_tag`). Without the anchor, the convergence graph has nothing to attach mappings to → 0% resolution.

- Snapshot: `before_EZ55` (110,333 edges total)
- KG rebuild: deferred (docker skipped)
- See `plans/ez55_deploy.md` for full details

### Phase 10-N: Remaining strains (NOT STARTED)

Next strain sequence:
MIT1002 → HOT1A3

Each strain: snapshot → `build_gene_id_mapping --strains X --force` → `resolve_paper_ids --force` → `check-gene-ids` → rebuild KG → compare snapshot.

## Key Design Decisions (for reference)

- **Three-tier classification**: Tier 1 (gene-unique) → `specific_lookup` 1:1; Tier 2 (protein-level, paralogs OK) → `multi_lookup`; Tier 3 (generic names) → `multi_lookup`
- **Iterative convergence** (not Union-Find): all sources processed together in passes until stable; typically 2-3 passes; order-independent
- **Singleton rule for Tier 2+3**: `multi_lookup` entries used for resolution only when they map to exactly 1 gene in the organism
- **Never silent skip**: every row gets explicit `resolution_method` string in `_resolved.csv`
- **List expansion**: cells like `"PMM0001, PMM0002"` split on `,` and `;`; each part tried; full value tried first

## Files NOT changed (intentionally)

- `multiomics_kg/adapters/omics_adapter.py` — reads `_resolved.csv` unchanged
- `gene_mapping_supp.csv` — not yet deleted (some strains may still need legacy fallback until v2 rebuilt for all)
- `gene_id_mapping.json` files in cache — not yet regenerated (step 3 pending per strain)
