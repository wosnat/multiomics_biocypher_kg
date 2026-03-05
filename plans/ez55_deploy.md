# EZ55 Deployment Plan — v2 Gene ID Mapping

## Background

Alteromonas macleodii EZ55 has two papers with expression data:

| Paper | ID format | Source annotation | Current resolution |
|-------|-----------|-------------------|-------------------|
| Barreto 2022 | `EZ55_NNNNN` (5-digit) | GCA_901457815.2 (ENA/NCBI) | **95%** via `locus_tag:symbol` — already working |
| Hennon 2017 | `AEZ55_NNNN` (4-digit) | Researcher's own annotation of IMG genome 2785510739 | **0%** — completely different annotation |

The `AEZ55_` IDs are from the researcher's independent gene-calling (`genbank_gmh_script`) on a JCVI-sequenced draft genome (IMG genome 2785510739, 3 scaffolds). The draft was never deposited as a genome assembly in NCBI (only raw reads at SRX022631, BioProject PRJNA50031). The complete genome (GCF_901457815.2) was later submitted to ENA with a fresh annotation using `EZ55_NNNNN` locus tags.

No string transformation (padding, prefix change) can bridge `AEZ55_` to `EZ55_` — the gene-calling and numbering are independent.

## Key Findings (2026-03-04)

IMG files received from JGI tape archive. **Critical discovery**: the IMG download contains TWO independent annotations of the same draft genome:

| Annotation | ID format | Source | Genes | Scaffolds |
|---|---|---|---|---|
| IMG (img_core_v400) | `altEZ55_NNNNN` (5-digit) | IMG GFF + protein FASTA | 4,237 proteins / 4,347 CDS | 2: `cAltEZ55` (4,713 kb), `pAltEZ55` (204 kb) |
| Researcher (genbank_gmh_script) | `AEZ55_NNNN` (4-digit) | `ez55.gff` from Sonya Hennon | 5,016 genes | 3: `scf7180000008833` (64 kb), `scf7180000008834` (2,076 kb), `scf7180000008835` (2,780 kb) |

**`AEZ55` ≠ `altEZ55`**: Despite similar names, these are completely different gene-callings. Same ordinal number = different gene (e.g., `AEZ55_0001` at scf...33:1-467(-) vs `altEZ55_00001` at cAltEZ55:1-1593(+) = different coordinates, different strand, different product).

**The paper DE tables use `AEZ55_NNNN`** (from the researcher's annotation, `ez55.gff`).

**IMG protein FASTA headers**: `>2785712097 altEZ55_00001 product...` — numeric OID first, then `altEZ55_NNNNN` as second token.

**Scaffold correspondence**: The 3 researcher scaffolds (~4,920 kb total) were merged into 2 IMG scaffolds (~4,917 kb total). Exact scaffold-to-scaffold mapping requires DNA alignment since names differ and sizes don't directly correspond.

## Prerequisites

- [x] IMG files requested from JGI tape archive (IMG genome 2785510739)
- [x] IMG files received — at `data/Prochlorococcus/papers_and_supp/Hennon 2017/IMG_2785510739/`
- [x] Researcher's GFF received — `data/Prochlorococcus/papers_and_supp/Hennon 2017/ez55.gff`
- [x] **Author's protein FASTA + GBF received (2026-03-04)** — `data/Prochlorococcus/papers_and_supp/Hennon 2017/EZ55_annotation/` (ez55_aa.fasta: 4930 proteins; ez55-1/2/3.gbf: 3 contigs)

### Available IMG files

In `IMG_2785510739/IMG Assembled Data/`:
- `2785510739.genes.faa` — IMG protein FASTA, `altEZ55_NNNNN` IDs (4,237 proteins)
- `2785510739.genes.fna` — IMG gene nucleotide FASTA
- `2785510739.gff` — IMG GFF with `locus_tag=altEZ55_NNNNN`
- `2785510739.fna` — scaffold FASTA (cAltEZ55 + pAltEZ55, ~4.9 MB)
- `2785510739.cog.tab.txt`, `.ko.tab.txt`, `.pfam.tab.txt`, `.tigrfam.tab.txt`, etc. — functional annotations

In `IMG_2785510739/IMG Data/`:
- `2785510739.tar.gz` — archive of the same assembled data + README
- `genbank_file_191530_106661.gb` — GenBank flat file with `altEZ55_NNNNN` annotation

## Step 0: Build AEZ55 → EZ55 protein mapping — DONE

**Script**: `scripts/map_img_to_ncbi_proteins.py`
**Status**: COMPLETE — author (Gwenn Hennon) provided protein FASTA + GBF files (2026-03-04)

### Method: 3-phase protein sequence matching

Uses researcher's protein FASTA (`EZ55_annotation/ez55_aa.fasta`, 4930 proteins with `AEZ55_NNNN` headers):

| Phase | Method | Matches | Cumulative |
|-------|--------|---------|------------|
| 1 | Exact sequence match | 2,971 | 2,971 (60.3%) |
| 2 | Subsequence match (≥95% overlap) | 136 | 3,107 (63.0%) |
| 3 | Diamond blastp (≥80% id, ≥60% qcov) | 946 | 4,053 (82.2%) |

**Fragment deduplication**: The researcher's draft genome annotation splits many NCBI genes into multiple smaller ORFs (due to frameshifts/assembly gaps in the draft). When multiple AEZ55 fragments hit the same NCBI gene, only the **longest fragment** is kept — it captured the most RNA-seq reads and best represents the gene's expression. 496 shorter fragments were discarded to avoid duplicate/conflicting expression edges. Only query coverage is required (not subject coverage) because a valid fragment matches well across its own length but may cover only part of the larger NCBI protein.

**DE gene coverage**: 345/422 (81.8%) for table 3, 98/125 (78.4%) for table 4.

**Remaining unmapped** (~20%):
- 47 DE genes: draft-specific ORFs with no NCBI counterpart (researcher called 4930 proteins vs NCBI's 4096)
- 49 DE genes: discarded fragments (shorter sibling of same NCBI gene already mapped)
- 381 total unmatched proteins (not DE)

### Input/Output files

- **Researcher FASTA**: `data/Prochlorococcus/papers_and_supp/Hennon 2017/EZ55_annotation/ez55_aa.fasta`
- **Researcher GBF**: `data/Prochlorococcus/papers_and_supp/Hennon 2017/EZ55_annotation/ez55-{1,2,3}.gbf`
- NCBI protein FASTA: `cache/data/Alteromonas/genomes/EZ55/protein.faa` (4096 proteins)
- Gene mapping: `cache/data/Alteromonas/genomes/EZ55/gene_mapping.csv`
- **Output**: `data/Prochlorococcus/papers_and_supp/Hennon 2017/aez55_to_ez55_id_translation.csv` (4053 rows)

### Reproducing

Requires `diamond` (v2.1.9+):
```bash
uv run python scripts/map_img_to_ncbi_proteins.py \
  --img-faa "data/Prochlorococcus/papers_and_supp/Hennon 2017/EZ55_annotation/ez55_aa.fasta" \
  --ncbi-faa cache/data/Alteromonas/genomes/EZ55/protein.faa \
  --gene-mapping cache/data/Alteromonas/genomes/EZ55/gene_mapping.csv \
  --output "data/Prochlorococcus/papers_and_supp/Hennon 2017/aez55_to_ez55_id_translation.csv" \
  --de-csvs "data/Prochlorococcus/papers_and_supp/Hennon 2017/ALT DE genes 41396_2018_bfismej2017189_moesm29_esm.csv" \
            "data/Prochlorococcus/papers_and_supp/Hennon 2017/ALT axenic DE genes 41396_2018_bfismej2017189_moesm30_esm.csv"
```

## Step 1: Add id_translation to Hennon 2017 paperconfig — DONE

Added to `data/Prochlorococcus/papers_and_supp/Hennon 2017/paperconfig.yaml`:

```yaml
id_translation_ez55_author:
  type: id_translation
  filename: "data/Prochlorococcus/papers_and_supp/Hennon 2017/aez55_to_ez55_id_translation.csv"
  organism: "Alteromonas EZ55"
  id_columns:
    - column: "locus_tag"
      id_type: locus_tag
    - column: "aez55_id"
      id_type: old_locus_tag
```

**Important**: Both columns must be declared — the `locus_tag` column provides the anchor (known canonical locus tags) and `aez55_id` provides the new mapping.

Also updated `supp_table_3` and `supp_table_4` `id_columns` to `id_type: old_locus_tag`.

## Step 2: Snapshot current KG — DONE

```
before_EZ55: 110,333 total edges, 20 publications
```

## Step 3: Rebuild gene_id_mapping.json — DONE

```bash
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains EZ55 --force
```

Results: `specific_lookup` grew from 10,036 → 22,322 entries. 0 conflicts. Converged in 2 passes.

## Step 4: Re-resolve paper CSVs — DONE

```bash
uv run python -m multiomics_kg.download.resolve_paper_ids --force --papers "Hennon 2017"
```

Results: 84.0% overall (546/650). Table 3: 345/422 (81.8%). Table 4: 98/125 (78.4%). All via `tier1:locus tag`.

## Step 5: Verify match rates — DONE

check-gene-ids: All 3 analyses MATCH. 546 matched, 104 pre-unresolved (no edges), 0 mismatches.

## Step 6: Rebuild KG and compare snapshot — PENDING

```bash
docker compose up -d --build
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before_EZ55
```

## Step 7: Update status — PENDING

- Update `plans/gene_id_mapping_v2_status.md` Phase 9 (EZ55)
- Update MEMORY.md — DONE

## Resolved Questions

1. **IMG protein FASTA header format**: Headers use numeric OIDs as first token, `altEZ55_NNNNN` as second token (e.g., `>2785712097 altEZ55_00001 product...`). The `altEZ55` IDs are from IMG's own annotation (img_core_v400), NOT the researcher's `AEZ55` IDs. The `remap_img_headers` function needs to extract the second token.

2. **AEZ55 vs altEZ55**: These are TWO DIFFERENT gene-callings on the same draft genome. `AEZ55_NNNN` (researcher's, 5,016 genes) ≠ `altEZ55_NNNNN` (IMG's, 4,237 proteins). Same ordinal = different gene. A coordinate-based mapping step is needed.

3. **Scaffold FASTA**: Available at `IMG_2785510739/IMG Assembled Data/2785510739.fna` (2 scaffolds: cAltEZ55 4,713 kb + pAltEZ55 204 kb). Researcher's GFF has 3 scaffolds (~64 + 2,076 + 2,780 kb). Scaffolds were merged/rearranged between the two annotations.

## Resolved Questions (2026-03-04)

1. **Researcher response**: Gwenn Hennon provided protein FASTA + GBF files. The preferred (single-step) approach was used.

2. **High-numbered AEZ55 genes**: Researcher's FASTA has 4930 proteins up to AEZ55_5024, covering the high-ID DE genes. No longer an issue.

3. **Barreto GTF annotation_gff**: The `annotation_gtf_ez55` collected 0 GFF rows (confirmed in build_gene_id_mapping output). It adds no new bridges — can be removed from paperconfig if desired.

## Open Questions

1. **Draft gene fragmentation**: The researcher's draft annotation splits ~496 NCBI genes into 2+ smaller ORFs (frameshifts in draft). Longest-fragment-wins dedup handles this for expression edges, but the underlying issue is real assembly differences between the draft and finished genome.

## Cleanup

After successful deployment:
- Delete `*_with_locus_tag.csv` files from Hennon 2017 (empty locus_tag columns, useless)
- Delete `fix_gene_ids_report.md` from Hennon 2017 (outdated, 0% mapping)
- Consider deleting `*_with_locus_tag.csv` from Barreto 2022 (identity mapping, superseded by resolver)
