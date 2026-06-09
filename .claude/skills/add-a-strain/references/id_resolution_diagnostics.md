# ID Resolution Diagnostics

Open this when `/check-gene-ids` returns < 80% on a paper. Each section
describes a specific failure mode, the diagnostic to confirm it, and the fix.

## IDs that are themselves locus_tags but not in `specific_lookup`

This was a v2 bug (now fixed in `resolve_row`): canonical locus_tags from old
GCA annotations (e.g. `PMT9312_1938`) are in the `genes` dict but were not in
`specific_lookup` (which maps alt_ids → canonical locus_tags). The fix adds a
`locus_tag:<col>` pass that checks `mapping_data.locus_tags` directly.

If you see high `unresolved` rates for IDs that look like valid locus_tags,
check whether they appear in `gene_id_mapping.json` `genes` dict:

```python
import json
with open('cache/data/Prochlorococcus/genomes/MIT9312/gene_id_mapping.json') as f:
    m = json.load(f)
print('in genes:', 'PMT9312_1938' in m['genes'])
print('in specific:', 'PMT9312_1938' in m['specific_lookup'])
```

If `in genes=True` and `in specific=False`, the ID is a canonical locus_tag
and should now resolve via `locus_tag:<col>` after rebuilding the mapping.

## IDs with intermediate/historical protein accessions (ABS*, ABB*)

NCBI now returns only WP_* proteins for GCF and ABB_* for GCA. Intermediate
accessions (e.g. ABS83097.1 used in Hennon 2017) are no longer available via
any current NCBI API. Options:

- **Defer**: skip these papers in the current strain cycle
- **Manual mapping**: fetch individual accessions via `efetch` and build a
  hand-crafted `id_translation` CSV

## Annotation table column contains old locus tags but declared as `gene_name`

**Symptom**: id_translation resolution is low even though the annotation table
has locus-tag-like values in a column. Classic example: MIT9301 Anjur 2025
`9301_annotated_genome.csv`, column `uniprot_gene_name`.

**Root cause**: `_find_anchor` has three phases:

1. Phase 1 — Tier 1 normalized match in `specific_lookup` (skips Tier 3)
2. Phase 2 — whitespace-split only fires when the value **has a space**
   (skips bare single-token values)
3. Phase 3 — Tier 2 singleton in `multi_lookup` (skips Tier 1 and Tier 3)

If a column contains bare old locus tags like `P9301_06681` (no space, Tier 3
`gene_name`), all three phases skip it. But compound values like
`"dnaA P9301_00001"` work via Phase 2.

**Fix**: change the column's `id_type` from `gene_name` to `old_locus_tag`:

- Bare values (e.g. `P9301_06681`) → Phase 1 hits `specific_lookup['P9301_06681']` directly
- Compound values (e.g. `"dnaN P9301_00001"`) → Phase 2 whitespace-split still works
- Multi-word compounds (e.g. `"hisI hisIE P9301_06041"`) → Phase 2 split tries all tokens

**Diagnostic** — check if values in the column look like locus tags:

```python
import pandas as pd, re
ann = pd.read_csv('path/to/annotation_table.csv')
col = ann['uniprot_gene_name'].dropna()
locus_pat = re.compile(r'[A-Z0-9]+_\d{4,5}')
print('rows with locus-tag-like tokens:', col.str.contains(locus_pat).sum())
print('sample:', col[col.str.contains(locus_pat)].head(5).tolist())
```

If most values end in a `STRAIN_XXXXX` token, the column should be
`id_type: old_locus_tag`.

## Compound / comma-separated IDs

The resolver automatically splits on `,` and `;` (via `expand_list()`). Check
if IDs like `"csoS1, ccmk1"` resolve after split — the first matching part is
used.

## Multiple annotation generations (MIT9313 pattern)

Some strains were reannotated multiple times, creating multiple locus tag
styles for the same gene. MIT9313 has four generations:

1. `PMT0001` — original Cyanorak annotation (no underscore)
2. `PMT_0001` — second NCBI annotation (with underscore, different numbering)
3. `RG24_RS00005` — intermediate RS format
4. `AKG35_RS00005` — current NCBI locus_tag

**IMPORTANT**: `PMT_0003 ≠ PMT0003` — these are distinct IDs from different
annotation rounds, NOT string variants. Do not try to match them by string
manipulation.

The GCF GFF `old_locus_tag` field lists all old IDs (URL-encoded
comma-separated), but for ~360 MIT9313 genes the `PMT0###` form is missing. A
position-based fallback merge in `build_gene_mapping.py` catches 105 of these.
The remaining 569 are Cyanorak-only genes with no NCBI match.

## ProPortal-era Alternative locus IDs (GCA GFF Note field)

GCA GFF files for MED4, MIT9312, MIT9313, and NATL2A contain
`Note=Alternative locus ID:P9313_NNNNN` — ProPortal-era locus tags that differ
from the GenBank TSV numbering. These are used by older microarray studies
(e.g. Thompson 2011, GEO platform GPL11412).

The `build_gene_id_mapping` GCA GFF parser extracts these from the Note field
as `alternative_locus_tag` (Tier 1). If a paper uses IDs that look like
`P9313_NNNNN` but don't match the id_translation TSV `P9313_` column, check
the GCA GFF Note field — different numbering schemes exist.

## Footnote artifacts in CSV values

`_heuristic_candidates()` in `gene_id_utils.py` strips trailing `*` and `+`
(footnote markers from supplementary tables). If a paper uses other footnote
suffixes (e.g. `†`, `‡`), add them to the `rstrip()` call.

## CSV formatting errors (space vs underscore)

Some supplementary tables have spaces where underscores should be (e.g.
`P9313 01731` instead of `P9313_01731`). This is a data extraction artifact,
not a different ID format. Fix the source CSV directly — the resolver does not
normalize spaces to underscores.

Per the project convention, write `<stem>_modified.csv` rather than editing
the source CSV in place, and point the paperconfig at the modified file.

## Diagnosing Alteromonas strains (MIT1002, EZ55, HOT1A3)

**Core challenge**: Unlike Prochlorococcus (well-curated model organisms with
stable NCBI/UniProt entries), the Alteromonas strains were sequenced by
individual research groups and lack proper long-term NCBI/UniProt curation.
This means:

- **Multiple genome versions**: Different papers may use different
  assemblies/annotations of the same strain (e.g., EZ55 has a JCVI draft via
  IMG AND a later ENA complete genome)
- **Accession verification needed**: The NCBI accession in
  `cyanobacteria_genomes.csv` must be verified against what each paper
  actually used. Some accessions have been superseded or replaced (e.g.,
  MIT1002 `GCF_001077695.1` was replaced by `GCF_901457835.2`)
- **GFF provenance matters**: Check whether each paper's gene IDs come from
  NCBI PGAP, IMG, RAST, or a custom annotation pipeline — each produces
  different locus tags for the same genome
- **No Cyanorak**: No Cyanorak entries exist for Alteromonas, so there are no
  Cyanorak locus tags, GBK files, or cluster assignments. All genomic data is
  NCBI-only

**Shared taxid**: All three strains share NCBI taxid `28108`. UniProt data is
downloaded once to `cache/data/Alteromonas/uniprot/28108/`. This means WP_
protein accessions can be shared across strains — a protein matching a WP_ ID
might belong to MIT1002, EZ55, or HOT1A3. The `gene_mapping.csv` per strain
disambiguates via per-genome protein_id assignments.

**Locus tag prefixes**:

| Strain | Primary locus_tag | NCBI RS locus_tag | Old locus_tag | Notes |
|---|---|---|---|---|
| MIT1002 | `MIT1002_NNNNN` | `ALT831_RS*` | `MIT1002_NNNNN` | Step-by-5 numbering |
| EZ55 | `EZ55_NNNNN` | `ALTBGP6_RS*` | `EZ55_NNNNN` | Step-by-5 numbering |
| HOT1A3 | `ACZ81_NNNNN` | `ACZ81_RS*` | `ACZ81_NNNNN` | Step-by-5 numbering; locus_tag = RS prefix |

**Organism name matching**: `ORGANISM_TO_GENOME_DIR` in `gene_id_utils.py` has
entries for both `"alteromonas macleodii <strain>"` and `"alteromonas
<strain>"` (without "macleodii"). The matching uses substring containment
(`key in norm or norm in key`), so both forms work. However, bare
`"Alteromonas"` (without strain name, as used in ziegler 2025) does NOT match
any specific strain — it matches the genus-level treatment organism in
`treatment_organisms.csv` (taxid 28108). This is correct behavior.

**MIT1002 dual-assembly trap (2026-03-04)**:

- `GCF_001077695.1` (MIT/Chisholm Lab, 2015): draft Illumina, `TK37_RS*` locus tags
- `GCF_901457835.2` (ICBM, 2022): complete PacBio, `ALT831_RS*` / `MIT1002_NNNNN` (5-digit) — used in KG

Biller 2016/2018 used **RAST annotation on the draft genome**, producing
`MIT1002_NNNN` (4-digit) IDs that are unrelated to the current 5-digit locus
tags. Zero-padding 4→5 digits maps to WRONG genes (validated: 3.8% product
match). The old `_with_locus_tag.csv` workaround was also wrong — reverted.

**Resolution (2026-03-05)**: Author (Steve Biller) provided RAST protein FASTA
(`226.6.faa`, 4214 proteins). Diamond-based protein matching mapped 3891/4214
(92.3%) `fig|` IDs to canonical locus tags. Conversion table registered as a
second `id_translation` for transitive closure. See `plans/mit1002_deploy.md`.

**Coe 2024 MIT1002**: Uses `NCBI ID` = 5-digit `MIT1002_NNNNN` (canonical locus
tags). Resolves at 99.5% natively.
