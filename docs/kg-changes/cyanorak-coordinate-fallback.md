# KG Change: coordinates + sequence for genes NCBI never called (Cyanorak fallback)

**Date:** 2026-10-13
**Plan:** `plans/cyanorak_coordinate_sequence_fallback.md`
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P1 item 1
(friction log 2026-10-07 "91 MED4 genes have no coordinates").

## What changed

3,837 genes across the 21 Cyanorak-annotated strains had a product, roles and clusters but
**no `contig` / `start` / `end` / `strand`, no `protein_id` and no `sequence`**. Every one is
a gene Cyanorak called and PGAP did not: the Gene coordinate properties were read from the
NCBI columns only, and `sequence` was joined from `protein.faa` via the RefSeq accession.

| Strain | Genes | Cyanorak-only (now placed) |
|---|---|---|
| MIT9313 | 2,906 | 517 |
| MIT9303 | 3,091 | 504 |
| WH8109 | 2,681 | 362 |
| CC9311 | 3,038 | 339 |
| PAC1 (20 contigs) | 2,317 | 198 |
| WH8102 | 2,830 | 197 |
| MIT0801 | 2,330 | 185 |
| NATL2A | 2,210 | 152 |
| BL107 | 2,567 | 144 |
| NATL1A | 2,221 | 138 |
| WH7803 | 2,590 | 134 |
| MIT9202 | 2,005 | 124 |
| SB (4 contigs) | 1,950 | 111 |
| MIT0604 | 2,111 | 108 |
| SS120 | 1,959 | 100 |
| AS9601 | 1,945 | 99 |
| MED4 | 1,973 | 96 |
| MIT9312 | 1,970 | 94 |
| MIT9215 | 2,022 | 92 |
| MIT9301 | 1,926 | 76 |
| MIT9515 | 1,945 | 67 |
| **Total** | | **3,837 placed, 0 unplaced** |

### Coordinates (step 0 sub-step 5, `build_gene_mapping.py`)

Cyanorak places every one of these genes — `gene_mapping.csv` already carried
`start_cyanorak` / `end_cyanorak` / `strand_cyanorak` — but in **Cyanorak's coordinate
space**, which has no contig: Cyanorak stores each strain as one record.

- **Closed genomes (19 of 21):** one contig, Cyanorak coordinates = NCBI coordinates
  (median |Δstart| = 0 on the genes both callers agree on). The fallback is a copy.
- **Drafts (PAC1, SB):** Cyanorak's record is the NCBI contigs concatenated in accession
  order with a **50 bp spacer** between neighbours. Verified by tiling the per-contig offsets
  `_contig_offsets` already derives (for the stop-codon position merge) against the
  `##sequence-region` lengths from the NCBI GFF: a 50 bp gap at every boundary and PAC1's last
  tile ends at 1,842,113 = the Cyanorak record length. The contig is therefore recovered by
  **inverting the tiling**: find the tile containing both ends of the gene, its accession is
  the contig, subtract the tile offset for contig-relative coordinates. A gene whose interval
  crosses a tile boundary would span the spacer (an assembly artefact) and is left unplaced
  with a warning — none exists.

New helpers: `_ncbi_contig_lengths`, `_cyanorak_tiles`, `_fill_cyanorak_only_coordinates`.
Every row now carries **`coordinate_source ∈ {ncbi, cyanorak}`** (dense), passed through to
the Gene node. Where both callers annotated a gene, NCBI's coordinates stay authoritative —
the fallback only fills blanks.

### Sequence (build time, `cyanorak_ncbi_adapter.py`)

The Cyanorak GFF carries a URL-encoded `translation=` attribute on every CDS (1,150 / 1,150
Cyanorak-only genes over the four strains checked). The adapter loads it keyed by Cyanorak
ORF id (`_load_cyanorak_translations`) and uses it **only when the RefSeq join yields no
sequence**. Two normalisations make it follow `protein.faa`'s convention:

- the trailing **stop glyph** is stripped — Cyanorak writes one glyph per stop codon
  (`#` 1,386 / `+` 289 / `*` 285 on MED4);
- the **first residue is normalised to `M`** — Cyanorak writes the literal start-codon
  residue (V/L/I) where RefSeq writes Met.

Agreement check on genes both callers annotated (MED4 + MIT9313 + PAC1 + CC9311, after the
two normalisations): 5,845 identical to the RefSeq protein, 1,376 would have differed only
in the first residue, 1,603 are in-frame start re-calls, 5 differ otherwise.

Every `sequence` now carries a sibling **`sequence_source ∈ {ncbi, cyanorak}`** (sparse —
present exactly when `sequence` is).

## What did NOT change

- No `Protein` nodes or `Gene_encodes_protein` edges for these genes — they have no RefSeq
  accession. The sequence rides on the Gene as the retrieval-only payload it already was
  (not indexed, not in any routing signal).
- The position-fallback merge, `gene_id_mapping` (step 3) and post-import Cypher are
  untouched. The composite `Gene(organism_name, contig, start)` index and the explorer's
  `gene_neighbors` seek pick the placed genes up with no change.
- NCBI-only strains gain only the dense `coordinate_source = ncbi` column.

## Vocabulary

`Gene.coordinate_source` (closed, dense) and `Gene.sequence_source` (closed, sparse) in
`config/controlled_vocabularies.yaml`; both are R5 two-state strings naming both states.

## Verification

- Unit: `tests/test_build_gene_mapping.py::TestCyanorakOnlyCoordinates` (directive parser,
  tile builder, closed-genome identity, draft inversion, boundary guard, NCBI-only mode);
  `tests/test_cyanorak_ncbi_adapter.py::TestLoadCyanorakTranslations` /
  `TestCyanorakSequenceFallback` (glyph stripping, M-normalisation, precedence of the RefSeq
  sequence, absent-directory path).
- Cache hygiene (2026-10-13): over MED4, PAC1, CC9311, KT2440 the regenerated
  `gene_mapping.csv` differs from HEAD only in the new column and the four coordinate cells
  of Cyanorak-only rows; every placed draft row sits inside a real contig with
  `1 ≤ start < end ≤ length`.
- KG validity (after the next rebuild): `test_cyanorak_annotated_genes_have_coordinates`,
  `test_coordinate_source_dense_and_closed`, `test_sequence_source_present_iff_sequence` in
  `tests/kg_validity/test_structure.py`. They fail against a graph built before this change.

## Consumer notes

- Downstream neighbourhood rules no longer need an "unplaceable genes" bucket for
  Cyanorak strains.
- A `coordinate_source = cyanorak` gene's start codon is Cyanorak's call with no PGAP
  second opinion; on genes both callers annotate, starts differ by up to ~1 kb in rare
  cases, so treat the 5′ end as one caller's opinion.
- `sequence_source = cyanorak` sequences have no RefSeq accession to cite; cite the Cyanorak
  ORF id (`all_identifiers` carries it).
