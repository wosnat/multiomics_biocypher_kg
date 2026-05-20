# Cross-assembly Protein Sequence Bridging

When a paper uses gene IDs from a **different genome assembly or annotation
pipeline** than the canonical one in the KG, no string transformation can
bridge the IDs. Instead, use protein sequence matching to build a verified
cross-assembly translation table.

## When to use

The paper's gene IDs come from a draft assembly, RAST annotation, IMG
annotation, or any independent gene-calling run whose locus tags have no
overlap with the canonical NCBI locus tags. Typical signs:

- Zero-padding to match canonical IDs produces wrong genes
- Verify by product description match < 10% on a sample of mapped IDs

## Script

`scripts/map_img_to_ncbi_proteins.py`

Requires `diamond` (v2.1.9+) — install via `apt-get install diamond-aligner`
or download from GitHub.

## Inputs

- Draft/old protein FASTA with the paper's gene IDs as headers (from author, IMG, or RAST)
- NCBI canonical protein FASTA: `cache/data/<Organism>/genomes/<Strain>/protein.faa`
- Gene mapping: `cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv` (WP_ → locus_tag)
- Optional: DE CSVs for coverage reporting

## Three-phase matching

| Phase | Method | Typical yield |
|---|---|---|
| 1 | Exact protein sequence match | ~60% |
| 2 | Subsequence match (≥95% overlap) | +3% |
| 3 | Diamond blastp (≥80% id, ≥60% query coverage, no subject coverage filter) | +19% |

**Fragment deduplication**: Draft genome frameshifts split single canonical
genes into multiple shorter ORFs. When multiple draft IDs hit the same
canonical locus tag, keep only the **longest fragment** (it captured the most
RNA-seq reads). Discarded fragments are logged. No subject-coverage filter is
applied because fragments have high identity but low subject coverage by
design.

## Paperconfig entry

Add `id_translation` with a `generate` block so `build_gene_id_mapping.py`
(step 3) automatically runs the diamond script when the output is missing or
`--force` is given:

```yaml
id_translation_draft_author:
  type: id_translation
  filename: "path/to/id_translation.csv"
  organism: "<Organism Strain>"
  generate:
    method: diamond_protein_match
    source_fasta: "path/to/draft_proteins.fasta"
    source_id_col: draft_id
    img_gff: "path/to/draft.gff"   # optional, for header remapping
  id_columns:
    - column: "locus_tag"
      id_type: locus_tag        # ANCHOR — must be declared
    - column: "draft_id"
      id_type: old_locus_tag    # NEW mapping column
```

`generate` block fields:
- `method`: `diamond_protein_match` (runs `scripts/map_img_to_ncbi_proteins.py`)
- `source_fasta`: path to draft/old protein FASTA (relative to project root)
- `source_id_col`: column name for source IDs in output (maps to `--source-id-col`)
- `img_gff`: optional GFF for header remapping (maps to `--img-gff`)

`--ncbi-faa` and `--gene-mapping` are derived automatically from the
organism's genome_dir. `--output` comes from `filename`.

## Manual usage

```bash
uv run python scripts/map_img_to_ncbi_proteins.py \
  --img-faa "path/to/draft_proteins.fasta" \
  --ncbi-faa cache/data/<Organism>/genomes/<Strain>/protein.faa \
  --gene-mapping cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv \
  --output "path/to/id_translation.csv" \
  --source-id-col draft_id \
  --de-csvs "path/to/de_table1.csv" "path/to/de_table2.csv"
```

**CRITICAL**: Both the anchor column (`locus_tag`, containing canonical locus
tags) AND the new mapping column must be declared in `id_columns`. If only the
new column is declared, the convergence graph has no anchor to attach the
mappings to and resolution will be 0%.

## Expected coverage

~80–85% of DE genes. Unresolved genes are typically:

- Draft-specific ORFs with no canonical counterpart (draft predicts more genes than finished genome)
- Discarded fragments (shorter sibling of same canonical gene already mapped)

## Completed deployments

### EZ55 — Hennon 2017

Hennon 2017 uses `AEZ55_NNNN` (4-digit) gene IDs from the researcher's own
annotation of IMG draft genome 2785510739. Author (Gwenn Hennon) provided
protein FASTA (`EZ55_annotation/ez55_aa.fasta`, 4930 proteins).

**Results**: 4053/4930 matched (82.2%). DE coverage: 345/422 (81.8%) table 3,
98/125 (78.4%) table 4. 496 fragments discarded (longest-wins dedup). Output:
`aez55_to_ez55_id_translation.csv`.

See `plans/ez55_deploy.md` for full details.

### MIT1002 — Biller 2016/2018

Biller 2016/2018 use RAST `MIT1002_NNNN` (4-digit) IDs and coordinate-format
`RAST_region_ID` from a draft genome (GCF_001077695.1). Author (Steve Biller)
provided RAST protein FASTA (`226.6.faa`, 4214 proteins with `fig|226.6.peg.N`
headers).

Diamond-based matching: 3891/4214 (92.3%) mapped to canonical locus tags. The
output CSV (`rast_fig_id`, `locus_tag`) was registered as an id_translation.
The existing conversion table (`MIT1002_systematicnames_conversiontable.csv`)
was registered as a second id_translation, linking `fig|` IDs to 4-digit
`MIT1002_NNNN` genbank IDs and `contigNNNNN_start_stop` coordinates (both as
`alternative_locus_tag`, Tier 1). Transitive closure in
`build_gene_id_mapping.py` linked all three ID formats to the correct
locus_tags in 3 passes.

**Results**: Biller 2016 96.6% (740/766), Biller 2018 S6B 97.2% (175/180).
26+5 unresolved = RAST-specific ORFs with no Diamond match above threshold.

See `plans/mit1002_deploy.md` for full details.
