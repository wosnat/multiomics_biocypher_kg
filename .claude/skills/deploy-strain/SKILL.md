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
| GCA GFF | `genomic_gca.gff` in cache dir | Old locus_tag format (e.g. `PMT9312_1938`), old protein IDs (e.g. `ABS83190.1`), Cyanorak locus tag cross-refs |

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
| EZ55 | Hennon 2017 | `AEZ55_0520` | Old 4-digit locus tag | Needs GCA GFF for EZ55 |
| NATL2A | Tetu 2019 | Standard RS locus tags | GCF | 100% |
| CC9311 | Barreto 2022 | `sync_NNNN` | Cyanorak locus tag | Tier 1 `locus_tag_cyanorak` |
| WH8102 | Barreto 2022 | `SYNWNNNN` / gene symbols | Locus tags + annotation CSV | ~95% |

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

### Compound / comma-separated IDs

The resolver automatically splits on `,` and `;` (via `expand_list()`). Check if IDs like `"csoS1, ccmk1"` resolve after split — the first matching part is used.

## After all steps pass

Update `plans/gene_id_mapping_v2_status.md` to mark the strain done, then proceed to the next strain in the sequence:
`MED4 → MIT9313 → NATL2A → MIT9301 → AS9601 → NATL1A → RSP50 → CC9311 → WH8102 → MIT1002 → EZ55 → HOT1A3`
