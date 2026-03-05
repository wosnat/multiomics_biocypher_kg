---
name: deploy-strain
description: Deploy v2 gene ID mapping for a new strain. Snapshot KG, rebuild gene_id_mapping.json, re-resolve paper CSVs, verify match rates, rebuild KG, compare snapshot.
argument-hint: <strain-name>
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *), Bash(docker *)
---

# Deploy Strain (v2 Gene ID Mapping)

End-to-end checklist for deploying v2 gene ID mapping for one strain. Run these steps in order; each one is a checkpoint.

## Step-by-step

```bash
STRAIN=MIT9312   # change this

# 0 — snapshot current KG (Neo4j must be running)
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save before_${STRAIN}

# 1 — rebuild gene_id_mapping.json (+ gca gff + cds_from_genomic.fna downloaded automatically)
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains $STRAIN --force

# 2 — re-resolve paper CSVs
uv run python -m multiomics_kg.download.resolve_paper_ids --force

# 3 — verify match rates
uv run python .claude/skills/check-gene-ids/check_gene_ids.py

# 3b — if any paperconfig fix was needed: validate paperconfig + add unit tests
# uv run pytest tests/test_paperconfig_validation.py -v   # all 22 paperconfigs must pass
# uv run pytest tests/test_gene_id_graph.py -q            # add tests for new fix pattern

# 4 — rebuild KG
docker compose up -d --build

# 5 — compare snapshot
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before_${STRAIN}
```

Review the diagnostic report after step 1:
```
cache/data/<Organism>/genomes/<Strain>/gene_id_mapping_report.json
```

## What build_gene_id_mapping auto-includes (no paperconfig needed)

Beyond standard gene_annotations_merged.json, these are automatically picked up:

| Source | File | What it adds |
|--------|------|-------------|
| GCF cds_from_genomic.fna | `genomic_gca.gff` in cache dir | `cds_fna_id` Tier 1 — full `lcl|<chr>_cds_<protein>_<n>` IDs; also `locus_tag_ncbi`, `old_locus_tag`, `protein_id_refseq`, `gene_name` from FNA headers |
| GCA GFF | `genomic_gca.gff` in cache dir | Old locus_tag format (e.g. `PMT9312_1938`), old protein IDs (e.g. `ABS83190.1`), Cyanorak locus tag cross-refs, `Alternative locus ID` from Note field (ProPortal-era IDs like `P9313_NNNNN`, `PMED4_NNNNN`) |

Both files are downloaded by `scripts/prepare_data.sh` (step 0 NCBI sub-steps). No manual action needed.

## Adding non-standard ID sources (paperconfig)

When a paper uses IDs not in the NCBI/Cyanorak annotations, add entries to the paper's `paperconfig.yaml` **before** running step 1.

### annotation CSV (`id_translation`)

For files like `pro_9312_anot.csv` that map gene symbols / JGI IDs / UniProt names:

```yaml
id_translation_pro_9312:
  type: id_translation
  filename: "data/Prochlorococcus/papers_and_supp/<Author Year>/pro_9312_anot.csv"
  organism: "Prochlorococcus MIT9312"
  id_columns:
    - column: "GID"
      id_type: old_locus_tag      # Tier 1 — e.g. PMT9312_1938
    - column: "SYMBOL"
      id_type: gene_name          # Tier 3
    - column: "uniprot_acc"
      id_type: uniprot_accession  # Tier 2
  product_columns:
    - column: "PRODUCT"
```

### GFF annotation file (`annotation_gff`)

For a supplementary GFF that bridges protein_ids or alternate locus tags:

```yaml
annotation_gff_mit9312:
  type: annotation_gff
  filename: "data/Prochlorococcus/papers_and_supp/<Author Year>/annotation.gff"
  organism: "Prochlorococcus MIT9312"
```

The GFF `locus_tag` attribute → `locus_tag_ncbi` (Tier 1), `old_locus_tag` → `old_locus_tag` (Tier 1), `protein_id` → `protein_id_refseq` (Tier 2).

## Known ID formats per strain (learnt during migration)

| Strain | Paper | ID format | Source | Resolution |
|--------|-------|-----------|--------|------------|
| MIT9312 | Tetu 2019 | `PMT9312_RS00005` | GCF GFF | Tier 1 `locus_tag_ncbi` (100%) |
| MIT9312 | Barreto 2022 | `PMT9312_1938` (old 4-digit) | GCA GFF locus_tag | Tier 1 via `locus_tag` set check (100%) |
| MIT9312 | Barreto 2022 | gene symbols (`rpoA`, etc.) | Annotation CSVs | Tier 3 `gene_name` singleton |
| MIT9312 | Hennon 2017 | `lcl\|CP000111.1_cds_ABS83097.1_354` | Intermediate RefSeq annotation (~2012–2017); **no longer in NCBI** | Unresolved — defer |
| MIT9301 | Anjur 2025 | JGI IDs (`2626311821`) | `annotation_genome_9301` via `uniprot_gene_name` (old locus tags) | 97.1% — change `id_type: gene_name` → `id_type: old_locus_tag` for bare locus-tag columns |
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
| MIT1002 | Biller 2016 | `MIT1002_NNNN` (4-digit RAST) | RAST annotation on draft genome GCF_001077695.1 | **93.6% WRONG** — heuristic zero-pads to 5-digit, maps to wrong genes (3.8% product match). Needs author bridge |
| MIT1002 | Biller 2018 | `RAST_region_ID` (coordinate format `contig00001_start_stop`) | RAST annotation on draft genome | **0%** — coordinate format, needs coordinate→TK37_RS→WP_→MIT1002_NNNNN bridge |

## Removing `_with_locus_tag.csv` workarounds

Earlier strains were fixed using the `/fix-gene-ids` skill, which created `_with_locus_tag.csv` copies with an added `locus_tag` column and changed the paperconfig to `name_col: "locus_tag"`. With v2 mapping, this workaround is unnecessary — the resolver handles mixed ID columns natively.

**Before deploying a strain**, check if its paperconfig was previously modified by fix-gene-ids:

1. Check if `filename` points to a `_with_locus_tag.csv` file
2. Check if `name_col` was changed to `"locus_tag"` (from the original column name)
3. Check if `id_columns` references a `locus_tag` column that doesn't exist in the original CSV

**To revert**: use `git show <commit>:<paperconfig path>` to find the original paperconfig, then:
- Change `filename` back to the original CSV (without `_with_locus_tag` suffix)
- Change `name_col` back to the original column (e.g. `"Gene Name"`)
- Remove `id_columns` entries for columns that don't exist in the original CSV
- Keep valid additions like `prefiltered`, `pvalue_threshold`, `product_columns`, `environmental_conditions`

The v2 resolver creates `_resolved.csv` files with `locus_tag` and `resolution_method` columns — replacing the manual `_with_locus_tag.csv` approach cleanly.

## Diagnosing resolution failures

### IDs that are themselves locus_tags but not in specific_lookup

This was a bug (now fixed in `resolve_row`): canonical locus_tags from old GCA annotations (e.g. `PMT9312_1938`) are in the `genes` dict but were not in `specific_lookup` (which maps alt_ids → canonical locus_tags). The fix adds a `locus_tag:<col>` pass that checks `mapping_data.locus_tags` directly.

If you see high `unresolved` rates for IDs that look like valid locus_tags, check whether they appear in `gene_id_mapping.json` `genes` dict:

```python
import json
with open('cache/data/Prochlorococcus/genomes/MIT9312/gene_id_mapping.json') as f:
    m = json.load(f)
print('in genes:', 'PMT9312_1938' in m['genes'])
print('in specific:', 'PMT9312_1938' in m['specific_lookup'])
```

If `in genes=True` and `in specific=False`, the ID is a canonical locus_tag and should now resolve via `locus_tag:<col>`.

### IDs with intermediate/historical protein accessions (ABS*, ABB*)

NCBI now returns only WP_* proteins for GCF and ABB_* for GCA. Intermediate accessions (e.g. ABS83097.1 used in Hennon 2017) are no longer available via any current NCBI API. Options:
- **Defer**: skip these papers in the current strain cycle
- **Manual mapping**: fetch individual accessions via `efetch` and build a hand-crafted `id_translation` CSV

### Annotation table column contains old locus tags but declared as gene_name

Symptom: id_translation resolution is low even though the annotation table has locus-tag-like values in a column. Classic example: MIT9301 Anjur 2025 `9301_annotated_genome.csv`, column `uniprot_gene_name`.

Root cause: `_find_anchor` has three phases:
1. Phase 1 — Tier 1 normalized match in `specific_lookup` (skips Tier 3)
2. Phase 2 — whitespace-split only fires when the value **has a space** (skips bare single-token values)
3. Phase 3 — Tier 2 singleton in `multi_lookup` (skips Tier 1 and Tier 3)

If a column contains bare old locus tags like `P9301_06681` (no space, Tier 3 `gene_name`), all three phases skip it. But compound values like `"dnaA P9301_00001"` work via Phase 2.

**Fix**: change the column's `id_type` from `gene_name` to `old_locus_tag`:
- Bare values (e.g. `P9301_06681`) → Phase 1 hits `specific_lookup['P9301_06681']` directly
- Compound values (e.g. `"dnaN P9301_00001"`) → Phase 2 whitespace-split still works
- Multi-word compounds (e.g. `"hisI hisIE P9301_06041"`) → Phase 2 split tries all tokens

Diagnostic — check if values in the column look like locus tags:
```python
import pandas as pd, re
ann = pd.read_csv('path/to/annotation_table.csv')
col = ann['uniprot_gene_name'].dropna()
locus_pat = re.compile(r'[A-Z0-9]+_\d{4,5}')
print('rows with locus-tag-like tokens:', col.str.contains(locus_pat).sum())
print('sample:', col[col.str.contains(locus_pat)].head(5).tolist())
```

If most values end in a `STRAIN_XXXXX` token, the column should be `id_type: old_locus_tag`.

### Compound / comma-separated IDs

The resolver automatically splits on `,` and `;` (via `expand_list()`). Check if IDs like `"csoS1, ccmk1"` resolve after split — the first matching part is used.

### Multiple annotation generations (MIT9313 pattern)

Some strains were reannotated multiple times, creating multiple locus tag styles for the same gene. MIT9313 has four generations:
1. `PMT0001` — original Cyanorak annotation (no underscore)
2. `PMT_0001` — second NCBI annotation (with underscore, different numbering)
3. `RG24_RS00005` — intermediate RS format
4. `AKG35_RS00005` — current NCBI locus_tag

**IMPORTANT**: `PMT_0003 != PMT0003` — these are distinct IDs from different annotation rounds, NOT string variants. Do not try to match them by string manipulation.

The GCF GFF `old_locus_tag` field lists all old IDs (URL-encoded comma-separated), but for ~360 MIT9313 genes the PMT0### form is missing. A position-based fallback merge in `build_gene_mapping.py` catches 105 of these. The remaining 569 are Cyanorak-only genes with no NCBI match.

### ProPortal-era Alternative locus IDs (GCA GFF Note field)

GCA GFF files for MED4, MIT9312, MIT9313, and NATL2A contain `Note=Alternative locus ID:P9313_NNNNN` — ProPortal-era locus tags that differ from the GenBank TSV numbering. These are used by older microarray studies (e.g. Thompson 2011, GEO platform GPL11412).

The build_gene_id_mapping GCA GFF parser extracts these from the Note field as `alternative_locus_tag` (Tier 1). If a paper uses IDs that look like `P9313_NNNNN` but don't match the id_translation TSV P9313_ column, check the GCA GFF Note field — different numbering schemes exist.

### Footnote artifacts in CSV values

`_heuristic_candidates()` in `gene_id_utils.py` strips trailing `*` and `+` (footnote markers from supplementary tables). If a paper uses other footnote suffixes (e.g. `†`, `‡`), add them to the `rstrip()` call.

### CSV formatting errors (space vs underscore)

Some supplementary tables have spaces where underscores should be (e.g. `P9313 01731` instead of `P9313_01731`). This is a data extraction artifact, not a different ID format. Fix the source CSV directly — the resolver does not normalize spaces to underscores.

### Diagnosing Alteromonas strains (MIT1002, EZ55, HOT1A3)

**Core challenge**: Unlike Prochlorococcus (well-curated model organisms with stable NCBI/UniProt entries), the Alteromonas strains were sequenced by individual research groups and lack proper long-term NCBI/UniProt curation. This means:
- **Multiple genome versions**: Different papers may use different assemblies/annotations of the same strain (e.g., EZ55 has a JCVI draft via IMG AND a later ENA complete genome)
- **Accession verification needed**: The NCBI accession in `cyanobacteria_genomes.csv` must be verified against what each paper actually used. Some accessions have been superseded or replaced (e.g., MIT1002 `GCF_001077695.1` was replaced by `GCF_901457835.2`)
- **GFF provenance matters**: Check whether each paper's gene IDs come from NCBI PGAP, IMG, RAST, or a custom annotation pipeline — each produces different locus tags for the same genome
- **No Cyanorak**: No Cyanorak entries exist for Alteromonas, so there are no Cyanorak locus tags, GBK files, or cluster assignments. All genomic data is NCBI-only

Additional differences from Prochlorococcus/Synechococcus:

**Shared taxid**: All three strains share NCBI taxid `28108`. UniProt data is downloaded once to `cache/data/Alteromonas/uniprot/28108/`. This means WP_ protein accessions can be shared across strains — a protein matching a WP_ ID might belong to MIT1002, EZ55, or HOT1A3. The `gene_mapping.csv` per strain disambiguates via per-genome protein_id assignments.

**Locus tag prefixes**:

| Strain | Primary locus_tag | NCBI RS locus_tag | Old locus_tag | Notes |
|--------|-------------------|-------------------|---------------|-------|
| MIT1002 | `MIT1002_NNNNN` | `ALT831_RS*` | `MIT1002_NNNNN` | Step-by-5 numbering |
| EZ55 | `EZ55_NNNNN` | `ALTBGP6_RS*` | `EZ55_NNNNN` | Step-by-5 numbering |
| HOT1A3 | `ACZ81_NNNNN` | `ACZ81_RS*` | `ACZ81_NNNNN` | Step-by-5 numbering; locus_tag = RS prefix |

**Organism name matching**: `ORGANISM_TO_GENOME_DIR` in `gene_id_utils.py` has entries for both `"alteromonas macleodii <strain>"` and `"alteromonas <strain>"` (without "macleodii"). The matching uses substring containment (`key in norm or norm in key`), so both forms work. However, bare `"Alteromonas"` (without strain name, as used in ziegler 2025) does NOT match any specific strain — it matches the genus-level treatment organism in `treatment_organisms.csv` (taxid 28108). This is correct behavior: ziegler 2025's Alteromonas data refers to a heterotroph community, not a specific strain.

**Papers per strain**:

| Strain | Papers | Role |
|--------|--------|------|
| MIT1002 | Biller 2018 (DE, dark stress), Coe 2021 (DE, diel), biller 2016 (coculture treatment) | Both organism (DE tables) and treatment source (coculture with Pro) |
| EZ55 | Barreto 2022 (DE, pCO2 + coculture), Hennon 2017 (DE, elevated CO2) | Organism |
| HOT1A3 | Aharonovich 2016 (coculture treatment only) | Treatment source only — no DE tables targeting HOT1A3 genes |

**HOT1A3 special case**: HOT1A3 only appears as `treatment_organism` (edge source) in coculture experiments. It has no DE tables where HOT1A3 genes are the target. Therefore, gene ID resolution is only needed for the HOT1A3→gene edges in Prochlorococcus coculture papers (the Pro gene IDs need resolving, not HOT1A3's). HOT1A3's own gene_id_mapping.json exists but may not be exercised by any paper CSV.

**MIT1002 dual-assembly problem (2026-03-04)**: MIT1002 has TWO assemblies from different groups:
- `GCF_001077695.1` (MIT/Chisholm Lab, 2015): draft Illumina, `TK37_RS*` locus tags
- `GCF_901457835.2` (ICBM, 2022): complete PacBio, `ALT831_RS*` / `MIT1002_NNNNN` (5-digit) — used in KG

Biller 2016/2018 used **RAST annotation on the draft genome**, producing `MIT1002_NNNN` (4-digit) IDs that are unrelated to the current 5-digit locus tags. Zero-padding 4→5 digits maps to WRONG genes (validated: 3.8% product match). The old `_with_locus_tag.csv` workaround was also wrong — reverted.

**Current status**:
- Biller 2016 supp_table_3: 93.6% "resolved" via heuristic zero-padding — **INCORRECT mappings**
- Biller 2018 supp_table_s6b: 0% resolved (RAST coordinate format)
- Coe 2024 supp_table_4: 99.5% resolved correctly (5-digit MIT1002_NNNNN)

**Mapping chain needed**: RAST MIT1002_NNNN → (coordinate match on JXRW contigs via conversion table) → TK37_RS* → (WP_ bridge or Coe 2024 mapping) → MIT1002_NNNNN. Author contacted for direct mapping table or RAST protein FASTA.

**Old draft GFF saved**: `cache/data/Alteromonas/genomes/MIT1002/genomic_old_draft_GCF_001077695.1.gff`

See `plans/mit1002_deploy.md` for full details.

**Coe 2024 MIT1002**: Uses `NCBI ID` = 5-digit MIT1002_NNNNN (canonical locus tags). Resolves at 99.5%. Also contains `Gene ID` = `cds-TK37_RS*` mapping from old draft → current assembly for 3876 genes.

### Cross-assembly protein sequence bridging

When a paper uses gene IDs from a **different genome assembly or annotation pipeline** than the canonical one in the KG, no string transformation can bridge the IDs. Instead, use protein sequence matching to build a verified cross-assembly translation table.

**When to use**: The paper's gene IDs come from a draft assembly, RAST annotation, IMG annotation, or any independent gene-calling run whose locus tags have no overlap with the canonical NCBI locus tags. Typical signs: zero-padding to match canonical IDs produces wrong genes (verify by product description match <10%).

**Script**: `scripts/map_img_to_ncbi_proteins.py`

**Requires**: `diamond` (v2.1.9+) for Phase 3. Install via `apt-get install diamond-aligner` or download from github.

**Inputs**:
- Draft/old protein FASTA with the paper's gene IDs as headers (from author, IMG, or RAST)
- NCBI canonical protein FASTA: `cache/data/<Organism>/genomes/<Strain>/protein.faa`
- Gene mapping: `cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv` (WP_ → locus_tag)
- Optional: DE CSVs for coverage reporting

**Three-phase matching**:

| Phase | Method | Typical yield |
|-------|--------|---------------|
| 1 | Exact protein sequence match | ~60% |
| 2 | Subsequence match (≥95% overlap) | +3% |
| 3 | Diamond blastp (≥80% id, ≥60% query coverage, no subject coverage filter) | +19% |

**Fragment deduplication**: Draft genome frameshifts split single canonical genes into multiple shorter ORFs. When multiple draft IDs hit the same canonical locus tag, keep only the **longest fragment** (it captured the most RNA-seq reads). Discarded fragments are logged. No subject-coverage filter is applied because fragments have high identity but low subject coverage by design.

**Usage**:
```bash
uv run python scripts/map_img_to_ncbi_proteins.py \
  --img-faa "path/to/draft_proteins.fasta" \
  --ncbi-faa cache/data/<Organism>/genomes/<Strain>/protein.faa \
  --gene-mapping cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv \
  --output "path/to/id_translation.csv" \
  --de-csvs "path/to/de_table1.csv" "path/to/de_table2.csv"
```

**After mapping**, add to paperconfig as `id_translation`:
```yaml
id_translation_draft_author:
  type: id_translation
  filename: "path/to/id_translation.csv"
  organism: "<Organism Strain>"
  id_columns:
    - column: "locus_tag"
      id_type: locus_tag        # ANCHOR — must be declared
    - column: "draft_id"
      id_type: old_locus_tag    # NEW mapping column
```

**CRITICAL**: Both the anchor column (`locus_tag`, containing canonical locus tags) AND the new mapping column must be declared in `id_columns`. If only the new column is declared, the convergence graph has no anchor to attach the mappings to and resolution will be 0%.

**Expected coverage**: ~80-85% of DE genes. Unresolved genes are typically:
- Draft-specific ORFs with no canonical counterpart (draft predicts more genes than finished genome)
- Discarded fragments (shorter sibling of same canonical gene already mapped)

### Completed cross-assembly deployments

#### EZ55 — Hennon 2017 (DONE)

Hennon 2017 uses `AEZ55_NNNN` (4-digit) gene IDs from the researcher's own annotation of IMG draft genome 2785510739. Author (Gwenn Hennon) provided protein FASTA (`EZ55_annotation/ez55_aa.fasta`, 4930 proteins).

Results: 4053/4930 matched (82.2%). DE coverage: 345/422 (81.8%) table 3, 98/125 (78.4%) table 4. 496 fragments discarded (longest-wins dedup). Output: `aez55_to_ez55_id_translation.csv`.

See `plans/ez55_deploy.md` for full details.

#### MIT1002 — Biller 2016/2018 (BLOCKED — waiting for author)

Biller 2016/2018 use RAST `MIT1002_NNNN` (4-digit) IDs and coordinate-format `RAST_region_ID` from a draft genome (GCF_001077695.1). Zero-padding to 5-digit is WRONG (3.8% product match). Need author's RAST protein FASTA or direct mapping table.

Proposed mapping chain: RAST MIT1002_NNNN → (coordinate match via conversion table) → TK37_RS* → (WP_ bridge or Coe 2024 mapping) → MIT1002_NNNNN.

See `plans/mit1002_deploy.md` for full details.

## After all steps pass

If any new fix pattern was needed (paperconfig change, id_type correction, etc.), add unit tests documenting the failure mode and fix before marking done:

```bash
uv run pytest tests/test_paperconfig_validation.py -v   # validate all paperconfigs
# Add tests to tests/test_gene_id_graph.py covering:
# - the failure mode (what broke and why)
# - the fix (what the correct id_type / config is)
uv run pytest tests/test_gene_id_graph.py -q
```

Then update `plans/gene_id_mapping_v2_status.md` and the "Known ID formats per strain" table in this skill. Proceed to the next strain:
`AS9601 → RSP50 → MIT1002 → HOT1A3`

Already deployed: MIT9312, MIT9301, NATL1A, MED4, NATL2A, MIT9313, WH8102, CC9311, EZ55
Blocked: MIT1002 (waiting for author RAST protein FASTA — see `plans/mit1002_deploy.md`)
