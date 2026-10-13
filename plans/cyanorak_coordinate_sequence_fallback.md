# Cyanorak coordinate + sequence fallback for genes NCBI never called

**Status:** code + caches done 2026-10-13; Docker rebuild + `pytest -m kg` pending (user's call — staging holds :7687)
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P1 item 1,
extended on 2026-10-13 to carry the amino-acid sequence as well.

## Problem

3,837 genes across the 21 Cyanorak-annotated strains (MED4 96, MIT9313 517, MIT9303 504,
WH8109 362, CC9311 339, …) exist in the KG with a product, roles and clusters but no
`contig` / `start` / `end` / `strand`, no `protein_id` and no `sequence`. Every one is a
Cyanorak-only call: PGAP never annotated it, and the Gene coordinate properties are read
from the NCBI columns alone. Cyanorak does place and translate all of them —
`gene_mapping.csv` already carries `start_cyanorak` / `end_cyanorak` / `strand_cyanorak`
for every one (3,837 / 3,837), and the Cyanorak GFF carries a `translation=` attribute
on every one (1,150 / 1,150 over the four strains checked).

Consequences today: `gene_neighbors` cannot place them, the composite
`Gene(organism_name, contig, start)` index skips them, neighbourhood rules in downstream
analyses must report them as "unplaceable", and `gene_aa_sequence` returns nothing.

## Facts the design rests on (measured 2026-10-13)

- **Closed genomes (19 of 21):** one NCBI contig; Cyanorak coordinates equal NCBI
  coordinates (median |Δstart| = 0 on the genes both callers agree on).
- **Drafts (PAC1 20 contigs, SB 4):** Cyanorak's single record is the NCBI contigs
  concatenated in accession order with a **50 bp spacer** between each pair. Verified by
  tiling: the per-contig offsets `_contig_offsets` already derives, plus the contig lengths
  from the NCBI GFF `##sequence-region` lines, line up with a 50 bp gap at every boundary,
  and PAC1's last tile ends at 1,842,113 = the Cyanorak record length exactly. SB tiles the
  same way. So contig recovery is a deterministic inversion of the layout.
- **Translations:** the attribute is URL-encoded and ends in one of three stop glyphs
  (`#` 1,386 / `+` 289 / `*` 285 on MED4 — one per stop codon). After stripping the glyph,
  over MED4 + MIT9313 + PAC1 + CC9311: 5,845 identical to the NCBI protein, 1,376 differ
  only in the first residue (Cyanorak writes the literal start-codon residue V/L/I where
  NCBI writes M), 1,603 are in-frame start re-calls, 5 differ otherwise. No non-standard
  letters inside.

## Scope

**Changing**
1. `multiomics_kg/download/build_gene_mapping.py` (step 0 sub-step 5):
   - parse NCBI contig lengths from the GFF `##sequence-region` directives;
   - build Cyanorak-record tiles from `_contig_offsets` + lengths;
   - for every Cyanorak-only row, locate the tile containing both ends, set
     `seqid` / `start` / `end` / `strand` (contig-relative) and `coordinate_source = cyanorak`;
     a row whose interval crosses a tile boundary, or falls outside every tile, is left
     blank and warned about (count reported in the merge log line);
   - every NCBI-sourced row (and every row in NCBI-only mode) gets `coordinate_source = ncbi`.
2. `config/gene_annotations_config.yaml`: `coordinate_source` passthrough field.
3. `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`:
   - load `{data_dir}/cyanorak/*.gff` translations keyed by Cyanorak id (line parser,
     mirrors `_load_protein_sequences`): URL-decode, strip one trailing `*`/`#`/`+`,
     normalise the first residue to `M`;
   - when the RefSeq join yields no sequence and the gene has a `locus_tag_cyanorak` with a
     translation, emit `sequence` from it; emit `sequence_source ∈ {ncbi, cyanorak}` beside
     every `sequence` (sparse: absent when there is no sequence).
4. `config/schema_config.yaml`: `coordinate_source: str`, `sequence_source: str` on gene.
5. `config/controlled_vocabularies.yaml`: `Gene.coordinate_source` (closed, dense),
   `Gene.sequence_source` (closed, sparse).
6. Tests: `tests/test_build_gene_mapping.py` (tiling inversion incl. the boundary case and
   closed-genome identity), `tests/test_cyanorak_ncbi_adapter.py` (translation parser +
   fallback precedence), `tests/kg_validity/test_structure.py` (no Cyanorak-tagged gene
   lacks a contig; `coordinate_source` present on every gene).
7. Regenerate committed caches for the 21 Cyanorak strains: step 0 sub-step 5
   (`gene_mapping.csv`) and step 2 (`gene_annotations_merged.json`).
8. Docs: CLAUDE.md Key graph facts (Gene `sequence` bullet), CHANGELOG `[Unreleased]`,
   backlog item ticked, `docs/kg-changes/cyanorak-coordinate-fallback.md`.

**Not changing**
- No new `Protein` nodes or `Gene_encodes_protein` edges — these genes still have no
  RefSeq accession; the sequence rides on the Gene as the retrieval-only payload it is.
- No change to the position-fallback merge itself, to `gene_id_mapping` (step 3), or to
  post-import Cypher (the composite index picks the genes up unchanged).
- No Docker rebuild in this task: the staging container serves :7687 (see memory); the
  rebuild + `pytest -m kg` + snapshot run is the user's call.
- Where NCBI and Cyanorak both called a gene, NCBI's coordinates and protein stay
  authoritative — the fallback only fills blanks.

## Acceptance criteria

- `gene_mapping.csv` for every Cyanorak strain: 0 rows with `locus_tag_cyanorak` set and
  `start` blank, except rows the build warned about by name (expected 0).
- PAC1 / SB: every filled row's `seqid` is a real NCBI contig and `1 ≤ start < end ≤ length`.
- `coordinate_source` is dense on every gene in every strain's `gene_annotations_merged.json`.
- Adapter dry run on MED4: PMM0573 (cyabrB2), PMM0815 (hli), PMM1718 (petL) emit
  `contig = NC_005072.1`, integer `start`/`end`, `coordinate_source = cyanorak`, a sequence
  starting with `M`, `sequence_source = cyanorak`; PMM0001 keeps its NCBI values and
  `sequence_source = ncbi`.
- `pytest -m "not slow and not kg"` green.
