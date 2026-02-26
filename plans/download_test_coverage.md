# Test Coverage Improvement Plan: `multiomics_kg/download/`

Generated: 2026-02-26

## Overview

The download module has highly uneven test coverage. `build_gene_annotations.py` and most of
`download_genome_data.py` are well tested. Three files are almost entirely untested despite
containing complex, data-transforming logic:

- `build_gene_mapping.py` — zero tests; highest-risk file
- `build_protein_annotations.py` — zero tests
- `download_uniprot.py` — preprocessing helpers always mocked away

---

## Source File Inventory

### `download_genome_data.py`

| Function | Tested? |
|---|---|
| `_read_genomes_csv` | YES |
| `_get_org_group` | YES |
| `_ncbi_download_genome` | YES |
| `_cyanorak_download_file` | YES (but has a bug — see Priority 4) |
| `step1_ncbi` | YES |
| `step2_cyanorak` | YES |
| `step3_uniprot` | YES (force flag not tested) |
| `step4_eggnog` | YES |
| `step5_gene_mapping` | **NO — zero coverage** |
| `main` | NO (low priority) |

### `download_uniprot.py`

| Function | Tested? |
|---|---|
| `UniprotNodeField` enum methods | NO |
| `_extract_go_id` | **NO** |
| `_split_field` | **NO** |
| `fetch_raw_uniprot` | NO (always mocked) |
| `preprocess_uniprot_data` | **NO (always mocked)** |
| `download_uniprot` | YES (fetch+preprocess always mocked) |

### `build_gene_mapping.py`

| Function | Tested? |
|---|---|
| `load_gff` | **NO** |
| `ncbi_merge_cds_and_gene_entries` | **NO** |
| `_get_cynaorak_ID` | **NO** |
| `_get_cyanorak_id_map_from_gbk` | **NO** |
| `_get_ec_numbers_from_gbff` | **NO** |
| `load_gff_from_ncbi_only` | **NO** |
| `load_gff_from_ncbi_and_cyanorak` | **NO** |
| `build_gene_mapping` | **NO** |

### `build_gene_annotations.py`

| Function | Tested? |
|---|---|
| `_nonempty` | YES |
| `_split` | YES |
| `_coerce_to_tokens` | YES |
| All 6 transform functions | YES |
| `infer_organism_group` | YES |
| `load_gene_mapping` | YES |
| `load_eggnog` | YES |
| `load_uniprot` | YES |
| `AnnotationBuilder.build_wide` | YES |
| `AnnotationBuilder.build_merged` | YES (missing `passthrough_list`) |
| `process_strain` | YES |
| `main` | NO (low priority) |

### `build_protein_annotations.py`

| Function | Tested? |
|---|---|
| `_nonempty`, `_split`, `_coerce_to_tokens` | **NO** |
| `_tx_add_go_prefix`, `_tx_strip_function_prefix` | **NO** |
| `load_uniprot_columnar` | **NO** |
| `ProteinAnnotationBuilder.build_merged` | **NO** |
| `process_taxid` | **NO** |
| `main` | NO (low priority) |

---

## Prioritized Implementation Plan

### Priority 1 — New file: `tests/test_build_gene_mapping.py`

All functions are untested. This is the highest-risk file: bugs silently lose genes or
produce wrong locus-tag mappings with no error output.

#### 1a. `ncbi_merge_cds_and_gene_entries`

The most complex function in the module. Merges gene/CDS rows from a GFF3 DataFrame,
URL-decodes old_locus_tag values, splits multi-tag comma-separated values, applies a
column rename map, and handles missing columns gracefully.

Tests to write:

- Given a DataFrame with matching gene and CDS rows (joined on ID/Parent), verify the
  merge produces one row per gene with correctly renamed columns
- When `old_locus_tag_gene` column is absent (e.g. Alteromonas genomes), `locus_tag`
  falls back to `locus_tag_ncbi` without error
- When `old_locus_tag` contains URL-encoded characters (e.g. `PMT0003%2CPMT_0003`),
  the value is decoded and stored correctly in `old_locus_tags`
- When a gene has multiple old_locus_tags (comma-separated after decoding), each tag
  appears as a separate row in the exploded output before dedup

#### 1b. `load_gff_from_ncbi_and_cyanorak`

The outer merge + deduplication strategy is the trickiest logic in the module. Genes
with multiple old_locus_tags are exploded, then deduplicated preferring Cyanorak-matched
rows. Cyanorak-only genes must not be dropped.

Tests to write:

- NCBI gene matched to Cyanorak via locus_tag: result row carries data from both sources
- NCBI gene with no Cyanorak match: appears in output with Cyanorak columns as NaN
- Cyanorak gene with no NCBI match (e.g. plasmid ORF): appears in output under its
  Cyanorak locus_tag (not dropped by the outer join)
- NCBI gene with two old_locus_tags: only one output row after dedup (by `locus_tag_ncbi`)
- Final column rename via `_FINAL_MERGE_RENAME` is applied to the result

#### 1c. `build_gene_mapping` (entry point)

Routing logic and EC-number merge are both untested.

Tests to write:

- With only `ncbi_gff_file` (no Cyanorak files), delegates to `load_gff_from_ncbi_only`
- With Cyanorak files provided, delegates to `load_gff_from_ncbi_and_cyanorak`
- When `ncbi_gbff_file` exists, EC numbers are merged into `ec_numbers` column
- When `ncbi_gbff_file` is None or missing, `ec_numbers` column is absent (no crash)

#### 1d. `_get_ec_numbers_from_gbff`

Tests to write:

- Minimal synthetic GBFF with one CDS having one EC_number qualifier → correct dict
  `{"locus_tag": ["1.2.3.4"]}`
- CDS with no EC_number qualifier → not included in output
- Same locus_tag on multiple records → EC numbers accumulated (extend, not overwrite)

#### 1e. `_get_cyanorak_id_map_from_gbk` and `_get_cynaorak_ID`

Tests to write:

- Minimal GBK with one CDS having `cyanorak ORF Id: CYAN001` in qualifiers →
  result maps `"CYAN001"` to the locus_tag
- CDS with no cyanorak ORF qualifier → not included in output
- CDS with two cyanorak IDs → warning printed, first ID used

#### 1f. `step5_gene_mapping`

This orchestration step has zero coverage and contains meaningful skip/force/exception
logic separate from the lower-level functions.

Tests to write:

- Cache hit (`gene_mapping.csv` exists, `force=False`): `build_gene_mapping` not called,
  result row has `SKIP_EXISTS`
- Missing GFF (`step1_ncbi` not yet run): result row has `SKIP_NO_GFF`, no crash
- GBFF absent: `ncbi_gbff=None` passed to `build_gene_mapping`
- Cyanorak files absent: `cyan_gff=None, cyan_gbk=None` passed
- `build_gene_mapping` raises exception: result row has `FAILED`, exception not propagated
- `force=True`: existing `gene_mapping.csv` is regenerated

---

### Priority 2 — New file: `tests/test_build_protein_annotations.py`

Mirror the structure from `tests/test_build_gene_annotations.py`.

#### 2a. `ProteinAnnotationBuilder.build_merged`

The `bool` type resolver has a subtle `False`-not-`None` distinction:
`if val is not None and (_nonempty(val) or val is False)`. A regression here would
silently drop `is_reviewed=False` from all unreviewed proteins.

Tests to write:

- `reviewed="Reviewed"` → `is_reviewed=True` in output
- `reviewed="Unreviewed"` → `is_reviewed=False` present in output (must not be filtered
  out by the `_nonempty` check)
- `reviewed=None` → `is_reviewed` key absent from output
- `reviewed="unknown"` → `is_reviewed` key absent from output
- `passthrough_list` field with `add_go_prefix` transform: each token gets prefix applied
- `integer` field: string `"1234.0"` → `1234`
- `integer` non-parseable value → field absent from output

#### 2b. `load_uniprot_columnar`

Tests to write:

- Column-oriented JSON (field → {uniprot_acc: value}) → row-oriented dict
  ({uniprot_acc: {field: value}})
- All UniProt accessions in the source appear as keys in the result
- Fields missing for a given accession → absent from that accession's dict (no KeyError)

#### 2c. `process_taxid`

Tests to write:

- Missing `uniprot_preprocess_data.json` → prints skip message, no error, no output file
- Cache hit (`force=False`, output file exists) → skip message, file unchanged
- `force=True` → overwrites existing output file
- Normal run → creates output file with expected per-protein keys

---

### Priority 3 — Extend `tests/test_download_uniprot.py`

#### 3a. `preprocess_uniprot_data`

This function modifies the raw dict in-place with six distinct code paths and is always
mocked in existing tests. Bugs here would silently corrupt cached JSON files.

Tests to write:

- Integer field conversion: `LENGTH="1234"` → `1234`; value with comma `"1,234"` → `1234`
- String sanitisation: value containing `|` → replaced with `,`; `'` → `^`
- GO ID extraction: `go_c` field with `"term [GO:0005737]"` → `go_c_id` has `"0005737"`
- KEGG colon split: `"pro:PMM0001"` → `"PMM0001"`
- ENTREZ take-first: list `["1234567", "7654321"]` → `"1234567"`
- Single-element list unwrapping: `["value"]` → `"value"` for non-split fields

#### 3b. `_extract_go_id`

Tests to write:

- `"aspartate-semialdehyde dehydrogenase activity [GO:0004073]"` → `"0004073"`
- `"no GO term here"` → `None`
- Empty string → `None`

#### 3c. `_split_field`

Tests to write:

- `_split_field("xref_proteomes", "UP000000625,UP000002311")` → list of two accessions
- `_split_field("gene_names", "dnaN repA")` → list split on space
- `_split_field("xref_kegg", "pro:PMM0001")` → `"PMM0001"` (take part after colon)
- `_split_field("xref_geneid", "1234567;7654321")` → `"1234567"` (take first only)
- `_split_field("xref_proteomes", None)` → `None`
- Single-element result is unwrapped to scalar for most field types

---

### Priority 4 — Fixes and gaps in existing tests

#### 4a. Fix `test_creates_cyanorak_subdir_if_missing` in `tests/test_download_genome_data.py`

**Bug:** the test patches `_CURL_CLS` instead of `_REQUESTS_GET`, so `requests.get`
attempts a real network call. The directory assertion passes coincidentally because
`os.makedirs` runs before `requests.get`. This gives a false green in environments
with network access and fails unexpectedly when offline.

Fix: replace the broken patch target with the correct one for `requests.get`.
Also add an assertion that the URL passed to `requests.get` has the expected format.

#### 4b. `step3_uniprot` — force flag not tested

Add a test: `step3_uniprot(genomes, force=True)` passes `force=True` to the
`download_uniprot` call (currently no test checks this argument).

#### 4c. `passthrough_list` resolver in `AnnotationBuilder`

`_resolve_passthrough_list` is never invoked because `MINIMAL_CONFIG` in the test
module does not include a `passthrough_list`-type field.

Add: a `passthrough_list` entry to `MINIMAL_CONFIG` and verify a delimited string
produces a proper list; verify a pre-existing list passes through unchanged.

---

## Conventions to Follow

These patterns are used consistently in the existing test suite and must be followed
in all new test files.

### File and class structure

- One file per source module: `test_build_gene_mapping.py`, `test_build_protein_annotations.py`
- One `class Test<FunctionName>:` per function under test
- Descriptive method names: `def test_<scenario>(self, tmp_path):`
- `setup_method` to construct shared objects (e.g. builder instances)

### Fixture data

- Define synthetic CSV/TSV/GFF/GBK content as **module-level string constants**
  (e.g. `GFF_CONTENT = "##gff-version 3\n..."`)
- Write them to `tmp_path` inside fixtures or test methods using `Path.write_text()`
- For DataFrames passed directly: construct with `pd.DataFrame({...})` inline

### Mocking

- Use `unittest.mock.patch` as a context manager, not a decorator
- Patch at the **source module** for lazy imports inside function bodies:
  - `from Bio import SeqIO` inside a function → patch `"Bio.SeqIO"`
  - `from Bio import SeqIO` at module top → patch `"multiomics_kg.download.build_gene_mapping.SeqIO"`
- Patch module-level constants on the calling module:
  - `patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", str(tmp_path))`

### Path patching

When a function uses `PROJECT_ROOT` or other module-level path constants to build
file paths, patch the constant to `str(tmp_path)` so tests are fully self-contained:

```python
with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", str(tmp_path)):
    result = process_strain(row, config, force=False)
```

### Assertion style

- Prefer `assert result == expected` over `assertEqual` (no unittest base class)
- For DataFrame results: `assert set(df.columns) == {...}` and `assert len(df) == N`
- For file output: `assert Path(tmp_path / "gene_mapping.csv").exists()`
- For mock call verification: `mock_fn.assert_called_once_with(...)`
