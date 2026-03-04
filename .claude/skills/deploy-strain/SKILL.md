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
| EZ55 | Hennon 2017 | `AEZ55_0520` | Old 4-digit locus tag | Needs GCA GFF for EZ55 |
| NATL2A | Tetu 2019 | Standard RS locus tags | GCF | 100% |
| MIT9313 | Aharonovich 2016 | `PMT####` (no underscore, Cyanorak) | GCF old_locus_tag + Cyanorak GBK | 100.0% — standard Tier 1 |
| MIT9313 | Tolonen 2006 | `PMT####` + `PMT_or####` | Same | 99.6% |
| MIT9313 | Thompson 2011 | `P9313_NNNNN` (ProPortal) + `P9313_NNNNN (PMT####)` composite | GCA GFF `Alternative locus ID` in Note field; GEO platform GPL11412 | 85.6% — required CSV space→underscore fix + Note field parsing |
| MIT9313 | Fang 2019 | `PMT####` + `RNA_*` | Same as Aharonovich | 64.5% — RNA_* non-coding expected unresolved |
| CC9311 | Barreto 2022 | `sync_NNNN` | Locus tags (primary) + annotation CSV | 90.0% — 3 tRNA unresolved (expected) |
| WH8102 | Barreto 2022 | `SYNWNNNN` / gene symbols | Locus tags + annotation CSV | 92.3% — 1 RNA_15 tRNA unresolved (expected) |

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
`AS9601 → RSP50 → MIT1002 → EZ55 → HOT1A3`

Already deployed: MIT9312, MIT9301, NATL1A, MED4, NATL2A, MIT9313, WH8102, CC9311
