# Plan: Gene Annotation Preprocessing Improvements

## Context

Companion to `protein_annotation_improvements.md` (already implemented).
The protein annotation pipeline now outputs proper arrays for `catalytic_activities`,
`transmembrane_regions`, `cofactor_names`, `pathways`, `keywords`, `signal_peptide`.
Gene annotations need similar cleanup, plus gene-specific issues.

---

## What changed in protein_annotations.json (commit e517525)

The following fields in `protein_annotations.json` changed format. All cache files are updated.

| Field | Old type | New type | Notes |
|---|---|---|---|
| `catalytic_activities` | `str` (raw blob) | `list[str]` | split on `CATALYTIC ACTIVITY:`, cleaned |
| `transmembrane_regions` | `str` (raw blob) | `list[str]` | split into `'start..end'` range strings |
| `signal_peptide` | `str` or absent | `str` (`'1..26'`) | extract_signal_range transform applied |
| `cofactor_names` | `str` | `list[str]` | split on `COFACTOR:`, name extracted |
| `pathways` | `str` (`pathway_description`) | `list[str]` | renamed + split on `PATHWAY:` |
| `keywords` | `str` | `list[str]` | split on `;` |
| `pfam_ids` | absent / `str` | `list[str]` | from UniProt `xref_pfam` only (clean PF* IDs) |

**Field renames in protein_annotations_config.yaml** (also reflected in protein_annotations.json keys):
- `catalytic_activity` → `catalytic_activities`
- `cofactors` → `cofactor_names`
- `pathway_description` → `pathways`

These renames do NOT affect gene_annotations_config.yaml because it references the same field
names under `source: uniprot`. Verify that the `field:` values in gene config still match the
protein_annotations.json keys. Quick check:

```
python3 -c "
import json; data=json.load(open('cache/data/Alteromonas/uniprot/28108/protein_annotations.json'))
entry=list(data.values())[1]
print([k for k in entry if 'catalyt' in k or 'transmem' in k or 'pfam' in k])
"
# should print: ['pfam_ids', 'catalytic_activities', 'transmembrane_regions']
```

---

## Issues identified (from inspection of `gene_annotations_merged.json`)

### A. Fields inherited from protein_annotations.json — RESOLVED ✓

`catalytic_activities` and `transmembrane_regions` are now `list[str]` in protein_annotations.json.
The `gene_annotations_config.yaml` uses `type: passthrough` for both, which passes the raw value
through as-is (including lists) — no stringification occurs. Verified in EZ55 merged output:

```
catalytic_activities  list  ['Reaction=ATP-dependent breakage...']
transmembrane_regions list  ['12..31', '37..57', '69..90', ...]
```

**Config change (cosmetic, optional):** changing to `passthrough_list` would make the semantic
intent clearer but has no behavioral effect since the data is already a list and there is no
`transform` configured. Skip unless doing a broader cleanup pass.

**IMPORTANT: Verify field names in gene config.** The protein config renamed some fields:
- `catalytic_activity` → `catalytic_activities` (gene config uses `field: catalytic_activities` — OK ✓)
- `cofactors` → `cofactor_names` (gene config does NOT pull `cofactor_names` — not a problem, gene config doesn't have this)
- `pathway_description` → `pathways` (gene config does NOT pull `pathways` — not a problem)

Note: `subcellular_location` remains `type: passthrough` / `str` in both configs — no change needed.

### B. `cog_category` → `list[str]` *(pending)*

Currently a plain string. Multi-letter values exist (e.g. `'LU'`, `'ET'`, `'EGP'`)
meaning the gene belongs to multiple COG categories. Should be a list of single chars.

- Source: eggnog `COG_category` column (plain string like `"LU"`)
- Fix: add new transform `split_cog_category` that iterates characters → `['L', 'U']`;
  keep `type: passthrough` (since `passthrough_list` splits on a delimiter, not per-character)
- Schema: change `cog_category: str` → `cog_category: str[]` in `config/schema_config.yaml`
- Adapter: `COG_CATEGORY` already in `GeneNodeField` enum; the adapter's `isinstance(value, list)`
  branch (line 308) passes the list through transparently — **no adapter code change needed**

### C. `pfam_ids` → split into `pfam_ids` + `pfam_names` *(pending)*

EggNOG `PFAMs` column mixes PF* IDs (`PF00712`) and human-readable shortnames
(`DNA_pol3_beta`, `HATPase_c`). From EZ55 data:
```
pfam_ids list ['PF00712', 'PF02767', 'PF02768', 'DNA_pol3_beta', 'DNA_pol3_beta_2', 'DNA_pol3_beta_3']
```
Note: UniProt source (`protein_annotations.json.pfam_ids`) is now **clean** — only PF* IDs from
`xref_pfam`. The shortname contamination comes only from the eggnog `PFAMs` source.

Fix: add a second field `pfam_names` using a new transform `extract_pfam_names`
(inverse of existing `extract_pfam_ids`). Config change + one new transform.
The existing `pfam_ids` union entry for eggnog should also apply `extract_pfam_ids` to filter.

- Schema: add `pfam_names: str[]` to Gene node in `config/schema_config.yaml`
- Adapter: add `PFAM_NAMES = 'pfam_names'` to `GeneNodeField` enum in
  `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` (under Pfam section after `PFAM_DESCRIPTIONS`)

### D. `kegg_pathway` — deduplicate `ko` / `map` mirrors *(pending)*

EggNOG `KEGG_Pathway` emits both `ko03030` and `map03030` for the same pathway.
From EZ55 data: both `ko*` and `map*` entries present in same list (100% duplication).

Fix: add a `filter` key to the `kegg_pathway` union (like the existing `filter: "^GO:"` on
`go_terms`) to keep only `ko*` entries:
```yaml
kegg_pathway:
  type: union
  filter: "^ko"   # drop map* duplicates
  sources:
    - source: eggnog
      field: KEGG_Pathway
      delimiter: ","
```
This requires no new transform — the existing union `filter` mechanism handles it.

### E. `gene_synonyms` → add `gene_name_synonyms` + `alternative_locus_tags` *(pending)*

**Decision: Option 2** — add two typed fields while keeping `gene_synonyms`.

**Root cause of the current bug:** `gene_synonyms` config has `delimiter: " "` at the *field*
level, but `_resolve_union` reads `delimiter` from the *source* config. Result: tokens are split
on `,` (default), not space. So `"dnaN PMM0001"` becomes a single opaque token instead of two.

**Confirmed data patterns (MED4 + EZ55):**
- `gene_names` in gene_mapping.csv: space-separated, mixes gene names + locus tag aliases
  - e.g., `"dnaN PMM0001"`, `"purF PMM0004"`, `"TX50_RS00025"`, `"ALTBGP6_RS00025"`
- `gene_names_cyanorak`: mostly empty for Alteromonas; may mix name + locus tag for Pro
- `uniprot.gene_names` (protein_annotations.json): already `list[str]`, locus-tag-like entries
  - e.g., `['AVL55_09440']`, `['ALFOR1_30888']`

**Locus tag vs gene name classification:**
Examined all tokens across MED4 + EZ55. Two patterns uniquely identify locus tags:
1. All-uppercase + 4+ trailing digits: `PMM0001`, `SYNW1234` (old format, no underscore)
2. Contains underscore: `ALTBGP6_RS00025`, `TX50_RS00025`, `BFV95_0002`

Gene names that are all-uppercase also exist (`ALDH`, `ALB3`, `PTOX`, `APE1`) but have
at most 1–3 digit suffix → not caught by rule 1. The regex works cleanly on the full dataset.

**Classification regex (locus tag pattern):**
```
^([A-Z][A-Z0-9]*\d{4,}|.+_.+)$
```
- `PMM0001`, `SYNW1234` → matches ✓ | `ALDH`, `ALB3`, `APE1` → does NOT match ✓
- `ALTBGP6_RS00025`, `BFV95_0002` → matches ✓ | `dnaN`, `purF` → does NOT match ✓

**Implementation:**

**1. Fix `gene_synonyms`** (keep it — removing would break gene_id_utils.py, cyanorak adapter,
schema, snapshot tests): move `delimiter: " "` from field level → per-source:
```yaml
gene_synonyms:
  type: union
  sources:
    - source: gene_mapping
      field: gene_names_cyanorak
      delimiter: " "
    - source: gene_mapping
      field: gene_names
      delimiter: " "
    - source: uniprot
      field: gene_names   # already a list, delimiter ignored
```
Result: `["dnaN", "PMM0001"]` instead of `["dnaN PMM0001"]`.

**2. Add `gene_name_synonyms`** (gene-name-like tokens only):
```yaml
gene_name_synonyms:
  type: union
  filter_not: "^([A-Z][A-Z0-9]*\\d{4,}|.+_.+)$"
  sources:
    - source: gene_mapping
      field: gene_names_cyanorak
      delimiter: " "
    - source: gene_mapping
      field: gene_names
      delimiter: " "
    - source: uniprot
      field: gene_names
```

**3. Add `alternative_locus_tags`** (locus-tag-like tokens only):
```yaml
alternative_locus_tags:
  type: union
  filter: "^([A-Z][A-Z0-9]*\\d{4,}|.+_.+)$"
  sources:
    - source: gene_mapping
      field: gene_names_cyanorak
      delimiter: " "
    - source: gene_mapping
      field: gene_names
      delimiter: " "
    - source: uniprot
      field: gene_names
```

**4. Update `build_gene_annotations.py`** — generalize the hardcoded `gene_synonyms` exclusion
block to also remove canonical `gene_name` from `gene_name_synonyms`:
```python
# (existing block for gene_synonyms stays)
# Add analogous block for gene_name_synonyms:
if gene_name and "gene_name_synonyms" in result:
    filtered = [s for s in result["gene_name_synonyms"] if s != gene_name]
    if filtered:
        result["gene_name_synonyms"] = filtered
    else:
        del result["gene_name_synonyms"]
```

**5. Update `config/schema_config.yaml`** — add to Gene node properties:
```yaml
gene_name_synonyms: str[]
alternative_locus_tags: str[]
```

**6. Update `multiomics_kg/utils/gene_id_utils.py`** — add `"gene_name_synonyms"` and
`"alternative_locus_tags"` to the lookup field list for better ID matching coverage.

**7. Update `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`** — add two new enum members
to `GeneNodeField` (after `GENE_SYNONYMS`):
```python
GENE_NAME_SYNONYMS = 'gene_name_synonyms'
ALTERNATIVE_LOCUS_TAGS = 'alternative_locus_tags'
```
The adapter's `isinstance(value, list)` branch (line 308) already handles list values correctly.
`GENE_SYNONYMS` itself stays unchanged.

**8. Add `filter_not` support to `build_gene_annotations.py`** — implement in `_resolve_union`:
after applying the optional `filter` regex (keep matches), also apply `filter_not` (drop matches).
This is needed for `gene_name_synonyms` which uses `filter_not` to exclude locus-tag-like tokens.

---

## Files to change

| File | Changes |
|---|---|
| `config/gene_annotations_config.yaml` | D: filter kegg_pathway; B: cog_category transform; C: pfam_names field + fix eggnog pfam_ids entry; E: fix gene_synonyms delimiters + add 2 new fields |
| `multiomics_kg/download/utils/annotation_transforms.py` | B: `split_cog_category`; C: `extract_pfam_names` |
| `multiomics_kg/download/build_gene_annotations.py` | E: add `filter_not` to `_resolve_union`; generalize gene_name exclusion for gene_name_synonyms |
| `config/schema_config.yaml` | B: `cog_category` str→str[]; C: add `pfam_names: str[]`; E: add `gene_name_synonyms: str[]` + `alternative_locus_tags: str[]` to Gene node |
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | C: add `PFAM_NAMES`; E: add `GENE_NAME_SYNONYMS` + `ALTERNATIVE_LOCUS_TAGS` to `GeneNodeField` enum |
| `multiomics_kg/utils/gene_id_utils.py` | E: add new fields to lookup list |
| `tests/test_build_gene_annotations.py` | New test classes for each change |

---

## Implementation order (suggested)

1. ~~Fix A (type mismatch)~~ — **already resolved** by new protein_annotations.json format
2. Fix D (`kegg_pathway` dedup) — config-only, one line, zero risk
3. Fix B (`cog_category` → list) — single new transform, config line change
4. Fix C (`pfam_names`) — single new transform
5. Fix E (`gene_synonyms` + new typed fields) — several files, but each step is small

---

## Prompt for next chat

```
We are working on a BioCypher knowledge graph for marine cyanobacteria
(Prochlorococcus / Alteromonas). The preprocessing pipeline builds gene
annotation JSON files from three sources: gene_mapping.csv (NCBI+Cyanorak),
eggnog emapper annotations, and protein_annotations.json (UniProt, already
improved in a previous session).

The plan is in plans/gene_annotation_improvements.md. Read it carefully before
starting — it contains the full design for each issue including exact YAML snippets
and regex patterns.

Key files:
  config/gene_annotations_config.yaml              — merge rules (edit this most)
  config/schema_config.yaml                         — Gene node property definitions
  multiomics_kg/download/build_gene_annotations.py
  multiomics_kg/download/utils/annotation_transforms.py
  multiomics_kg/adapters/cyanorak_ncbi_adapter.py  — GeneNodeField enum (add new members)
  multiomics_kg/utils/gene_id_utils.py             — ID lookup field list
  tests/test_build_gene_annotations.py

Issue A (catalytic_activities/transmembrane_regions type mismatch) is already
resolved — passthrough passes lists through correctly. Skip A.

Please implement the remaining improvements in the order given in the plan:
  D → B → C → E

For each issue:
  1. Make the code/config changes described in the plan
  2. Add/update tests in tests/test_build_gene_annotations.py
  3. Run: pytest tests/test_build_gene_annotations.py -v

After all changes rebuild cache for EZ55:
  uv run python multiomics_kg/download/build_gene_annotations.py --strains EZ55 --force

Verify output:
  python3 -c "
  import json
  with open('cache/data/Alteromonas/genomes/EZ55/gene_annotations_merged.json') as f:
      data = json.load(f)
  entry = list(data.values())[0]
  for k, v in entry.items():
      print(k, type(v).__name__, repr(v)[:80])
  "

Run full test suite: pytest -m 'not slow and not kg' -v

IMPORTANT — adapter changes for new fields:
The GeneNodeField enum in cyanorak_ncbi_adapter.py is hardcoded. Fields absent from
the enum are silently dropped from Gene nodes even if present in the JSON. For each
new field, add an enum member:
  - PFAM_NAMES = 'pfam_names'               (Fix C — under Pfam section)
  - GENE_NAME_SYNONYMS = 'gene_name_synonyms' (Fix E — under Gene naming section)
  - ALTERNATIVE_LOCUS_TAGS = 'alternative_locus_tags' (Fix E — under Gene naming section)
cog_category (Fix B) already has COG_CATEGORY in the enum; no change needed there.
```
