# Known ID Formats Per Strain

Append-only history of which paper uses which ID format, where the ID comes
from, and the resolution outcome. Read this before deploying a new paper to
see whether a similar pattern has already been solved.

| Strain | Paper | ID format | Source | Resolution |
|---|---|---|---|---|
| MIT9312 | Tetu 2019 | `PMT9312_RS00005` | GCF GFF | Tier 1 `locus_tag_ncbi` (100%) |
| MIT9312 | Barreto 2022 | `PMT9312_1938` (old 4-digit) | GCA GFF locus_tag | Tier 1 via `locus_tag` set check (100%) |
| MIT9312 | Barreto 2022 | gene symbols (`rpoA`, etc.) | Annotation CSVs | Tier 3 `gene_name` singleton |
| MIT9312 | Hennon 2017 | `lcl\|CP000111.1_cds_ABS83097.1_354` | Intermediate RefSeq annotation (~2012–2017); **no longer in NCBI** | Unresolved — defer |
| MIT9301 | Anjur 2026 | JGI IDs (`2626311821`) | `annotation_genome_9301` via `uniprot_gene_name` (old locus tags) | 97.1% — change `id_type: gene_name` → `id_type: old_locus_tag` for bare locus-tag columns |
| NATL1A | He 2022 | `NATL1_NNNNN` locus tags + gene names (`wza`, `cyoA`, etc.) | Mixed `Gene Name` column | 98.8% — revert paperconfig to original CSV with `name_col: "Gene Name"`, v2 resolves both natively |
| EZ55 | Barreto 2022 | `EZ55_NNNNN` (5-digit) | GCA_901457815.2 (ENA/NCBI) canonical locus tags | **95%** — already working via `locus_tag:symbol` |
| EZ55 | Hennon 2017 | `AEZ55_NNNN` (4-digit) | Researcher's own gene-calling on IMG draft genome 2785510739 | **84%** — cross-assembly protein bridging via `map_img_to_ncbi_proteins.py` (3-phase: exact+subsequence+Diamond) |
| NATL2A | Tetu 2019 | Standard RS locus tags | GCF | 100% |
| MIT9313 | Aharonovich 2016 | `PMT####` (no underscore, Cyanorak) | GCF old_locus_tag + Cyanorak GBK | 100.0% — standard Tier 1 |
| MIT9313 | Tolonen 2006 | `PMT####` + `PMT_or####` | Same | 99.6% |
| MIT9313 | Thompson 2011 | `P9313_NNNNN` (ProPortal) + `P9313_NNNNN (PMT####)` composite | GCA GFF `Alternative locus ID` in Note field; GEO platform GPL11412 | 85.6% — required CSV space→underscore fix + Note field parsing |
| MIT9313 | Fang 2019 | `PMT####` + `RNA_*` | Same as Aharonovich | 64.5% — RNA_* non-coding expected unresolved |
| CC9311 | Barreto 2022 | `sync_NNNN` | Locus tags (primary) + annotation CSV | 90.0% — 3 tRNA unresolved (expected) |
| WH8102 | Barreto 2022 | `SYNWNNNN` / gene symbols | Locus tags + annotation CSV | 92.3% — 1 RNA_15 tRNA unresolved (expected) |
| MIT1002 | Coe 2024 | `MIT1002_NNNNN` (5-digit) | Canonical locus tags from GCF_901457835.2 | **99.5%** — resolves natively |
| MIT1002 | Biller 2016 | `MIT1002_NNNN` (4-digit RAST) | RAST annotation on draft genome GCF_001077695.1 | **96.6%** — cross-assembly protein bridging via `map_img_to_ncbi_proteins.py` + transitive closure through conversion table |
| MIT1002 | Biller 2018 | `RAST_region_ID` (coordinate format `contig00001_start_stop`) | RAST annotation on draft genome | **97.2%** — resolved via conversion table `assembly_coordinates` → `alternative_locus_tag` (Tier 1) through transitive closure |
| WH8109 | Doron 2016 | `Syncc8109_####_at` array probes | Strip `_at` → NATIVE locus tags of GCF_000161795.2 (canonical form IS `Syncc8109_####`) | **93.9%** — 168 unresolved are Moore/JCVI-annotation features dropped by RefSeq PGAP re-annotation of the same assembly (2,757 array features vs 2,707 current genes); no bridge work possible or needed |
| 83-1 | Koch 2019, Koch 2020 | `Alt831_NNNNN` (submitter / IMG annotation) | Matches the GCF_901446755.2 `old_locus_tag` `ALT831_NNNNN` case-insensitively (`locus_tag_ci`). **Trap:** this assembly's own locus tags are `MIT1002_RS*`, and our real MIT1002 (GCF_901457835.2) uses `ALT831_RS*` — the submitter swapped prefixes; sequence identity confirms the strains | **95.1%** (Koch 2019 S1), ~92% (Koch 2020 S7) — misses are tRNA/rRNA/hypotheticals/IS elements the RefSeq re-annotation dropped |
| MIT1002 | Koch 2020 | `MIT1002_NNNN` (4-digit RAST) + `contigNNNNN_start_stop` | Same RAST draft names Biller 2018 registered as `alternative_locus_tag`; 724/725 agree with a stop-codon coordinate check against the JXRW draft GFF. Do NOT zero-pad to `MIT1002_0NNNN` — that is the complete genome's numbering, a different gene | S3 not integrated (duplicate of Biller 2016) |
| WH8102, CC9311, BL107, RS9917 | Oleza 2015b | 2012 `NP_` / `YP_` / `ZP_` protein accessions | `id_translation` with `generate: diamond_protein_match` from the sequences in Table S4/S6 (IPG cannot map WH8102 NP_ to the GCF; `ZP_` not recognised as legacy) | ≥ 99% |

Add new rows above as strains are deployed. Keep columns consistent so the
patterns are pattern-greppable by future deployments.

## Already-deployed strains (no action needed)

MIT9312, MIT9301, NATL1A, MED4, NATL2A, MIT9313, WH8102, CC9311, EZ55, MIT1002, AS9601, RSP50, HOT1A3, SS120, BL107, HP15 (Marinobacter MarRef v6), Alt_MarRef, WH8109, 83-1.

## Backlog

- Pseudohoeflea, Thalassospira (ziegler 2025)
- Metabolomics paper integration strains — see `docs/superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md`
