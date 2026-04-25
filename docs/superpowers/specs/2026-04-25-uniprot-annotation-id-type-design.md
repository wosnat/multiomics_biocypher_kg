# `uniprot_annotation_string` id_type — design

**Date:** 2026-04-25
**Status:** Approved
**Driver paper:** Biller 2022 (vesicle proteomics, MIT9312 + MIT9313)
**Replaces workaround:** per-paper `_modified.csv` builders in `scripts/build_modified_csv/` whose only job is to extract `GN=<token>` into a clean column.

## Problem

Many proteomics supplementary tables include a column whose values are the canonical UniProt FASTA-header annotation:

```
Q31DF2_PROM9 Type II secretion system protein-like protein OS=Prochlorococcus marinus (strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1
RL33_PROM9 50S ribosomal protein L33 OS=Prochlorococcus marinus (strain MIT 9312) GN=rpmG PE=3 SV=2
```

The string embeds two resolvable IDs:
- the leading `<entry>_<ORG>` token — a UniProt **entry name** (Tier 1, 1:1 with a protein)
- the `GN=<token>` capture — a **gene name** or **locus_tag** (Tier 3 fallback; gets resolved via Tier 1 specific_lookup if the token is itself a real locus_tag)

Today, declaring such a column in `id_columns` with any existing id_type either fails (the long string never matches anything) or, if declared as `id_type: locus_tag`, pollutes the strain's `specific_lookup` with thousands of junk alt_ids — see the `FREE_TEXT_COLUMN_NAMES` anti-pattern guard in `scripts/validate_paperconfig.py`. The current escape hatch is a per-paper modified-CSV builder; this is friction that a recurring paper pattern shouldn't require.

## Solution at a glance

A new id_type `uniprot_annotation_string`. When `build_gene_id_mapping` (step 3) and `resolve_paper_ids` (step 4) encounter a column declared with this id_type, they parse the cell value with two regexes and feed the captured tokens through the existing v2 mapping machinery as `uniprot_entry_name` (Tier 1) and `gene_name` (Tier 3) respectively.

No regex in YAML; no new resolution path; no changes to `Cyanorak`/GFF/JGI handling. The narrow scope was confirmed during brainstorming.

## Paperconfig surface

```yaml
# biller 2022 S4 — after the change
s4_mit9312_enrichment:
  type: derived_metrics_table
  filename: data/.../emi15834-sup-0005-tables4_MIT9312_modified.csv
  organism: Prochlorococcus MIT9312
  experiment: vesicle_proteomics_mit9312
  name_col: "Protein ID"
  id_columns:
    - column: "Protein ID"
      id_type: uniprot_accession
    - column: "UniProt Annotation"
      id_type: uniprot_annotation_string      # NEW
  product_columns:
    - column: "UniProt Annotation"             # column may be in BOTH lists
    - column: "Prokka Annotation"
```

A column declared with `id_type: uniprot_annotation_string` is *expected* to hold free-text narrative; the existing `FREE_TEXT_COLUMN_NAMES` warning (which fires only on `id_type: locus_tag`) is unaffected.

## Extractor

Lives in `multiomics_kg/utils/gene_id_utils.py` (single source of truth — used by both step 3 and step 4):

```python
_UNIPROT_ANNOT_ENTRY_RE = re.compile(r"^([A-Z0-9]+_[A-Z0-9]+)\b")
_UNIPROT_ANNOT_GN_RE    = re.compile(r"\bGN=(\S+)")

def extract_uniprot_annotation_tokens(value: str) -> list[tuple[str, str]]:
    """Parse a UniProt FASTA-header style annotation string into (token, id_type) pairs.

    Returns [] when neither pattern matches — the cell is then a no-op for ID purposes.
    """
    out: list[tuple[str, str]] = []
    if not isinstance(value, str):
        return out
    s = value.strip()
    if not s:
        return out
    if (m := _UNIPROT_ANNOT_ENTRY_RE.match(s)):
        out.append((m.group(1), "uniprot_entry_name"))   # Tier 1
    if (m := _UNIPROT_ANNOT_GN_RE.search(s)):
        out.append((m.group(1), "gene_name"))             # Tier 3
    return out
```

Examples (verifiable against the biller 2022 S4 CSVs):

| Input fragment | Output |
|---|---|
| `Q31DF2_PROM9 ... GN=PMT9312_0032 PE=4 SV=1` | `[("Q31DF2_PROM9","uniprot_entry_name"), ("PMT9312_0032","gene_name")]` |
| `RL33_PROM9 ... GN=rpmG PE=3 SV=2` | `[("RL33_PROM9","uniprot_entry_name"), ("rpmG","gene_name")]` |
| `Q7V5C8_PROMM ... GN=PMT_1636 PE=4 SV=1` | `[("Q7V5C8_PROMM","uniprot_entry_name"), ("PMT_1636","gene_name")]` |
| `hypothetical protein` | `[]` |

GN= captures may be locus_tags (e.g. `PMT9312_0032`). They are still emitted with id_type `gene_name`. At resolution time those tokens hit `specific_lookup` via the existing locus_tag path before falling through to Tier 3 — no extra work needed in `resolve_row`. Treating GN= captures as Tier 3 is intentional: locus_tag-shaped tokens are 1:1 already (caught by Tier 1), while bare gene symbols (`rpmG`) are correctly handled as Tier 3 multi_lookup singletons.

## Integration points

### Step 3 — `multiomics_kg/download/build_gene_id_mapping.py`

In each of the four row-extraction functions, when `id_type == "uniprot_annotation_string"`, replace the single `(val, id_type)` append with iteration over `extract_uniprot_annotation_tokens(val)`. Empty list → cell contributes no pairs (skip).

Affected functions:
- `extract_rows_from_id_translation` (~L78–120) — id_columns list comprehension at L109–116
- `extract_rows_from_csv` (~L320 onward) — name_col block (L326–337) and id_columns block (L342–345)
- `extract_rows_from_gene_clusters` (~L385–410)
- `extract_rows_from_derived_metrics_table` (~L445–475) — name_col block and id_columns block

A small helper `_emit_pairs(val, id_type) -> list[tuple[str, str]]` keeps the four call sites symmetric: returns `extract_uniprot_annotation_tokens(val)` when `id_type == "uniprot_annotation_string"`, otherwise `[(val, id_type)]`.

### Step 4 — `multiomics_kg/utils/gene_id_utils.py::resolve_row`

`resolve_row` currently iterates `id_columns` (list of dicts) and uses only the `column` field — the `id_type` is ignored at lookup time. Add a small per-call dict `id_type_by_col = {c["column"]: c.get("id_type","other") for c in id_columns}`, and at the top of each per-column pass, if `id_type_by_col.get(col) == "uniprot_annotation_string"`, replace the cell value with the extracted tokens and run the existing Pass 1/1b/2/3 logic against each token in turn (preserving the first match's method string, e.g. `tier1:UniProt Annotation` when the entry name resolves).

The `name_col` is treated identically — if the paperconfig declares `name_col` to be a column whose `id_columns` entry has id_type `uniprot_annotation_string`, the same extraction happens.

### Validator — `scripts/validate_paperconfig.py:92`

Add `"uniprot_annotation_string"` to `VALID_ID_TYPES`. No other validator changes — the free-text-column guard at L342 only fires for `id_type: locus_tag`, so the new id_type is intentionally exempt (it's *meant* for narrative columns).

### Skill / docs

- `.claude/skills/paperconfig/SKILL.md` — add `uniprot_annotation_string` row to the **ID Types** table with a one-line description; mention briefly that columns of this type can be in `id_columns` and `product_columns` simultaneously.
- `CLAUDE.md` — no change required (the relevant table is in SKILL.md).

## Tests

### Unit tests (`tests/test_gene_id_utils.py` or new `tests/test_uniprot_annotation_extraction.py`)

`extract_uniprot_annotation_tokens` cases:
1. Full input — entry + GN= → 2 tokens, correct types
2. Entry only, no GN= → 1 token (entry only)
3. GN= only, no entry-shaped prefix → 1 token (gene_name only)
4. Neither — `"hypothetical protein"` → `[]`
5. Empty / whitespace / `None` / non-string → `[]`
6. Trailing punctuation: `GN=PMT_1636.` should capture `PMT_1636.`? Decision: capture verbatim; downstream Tier 1/3 lookup is the right place to handle this (matches existing `_heuristic_candidates` behavior of stripping `*+`).
7. Multiple GN= occurrences — pick the first (regex `search` semantics).

### Integration test

Load biller 2022 S4 MIT9312 CSV through `resolve_paper_ids` after step 3 has rebuilt `gene_id_mapping.json` for MIT9312 with this paperconfig included. Assert match rate ≥ 95% on the 207 rows. The expected resolution path:
- `Protein ID` column (`Q31DF2`) — already resolves via existing `uniprot_accession` path for known proteins.
- `UniProt Annotation` column — provides a backup via `uniprot_entry_name` (Tier 1) and GN= → locus_tag (Tier 1 via specific_lookup hit on the GN= token).

The cross-check: turning the new id_columns entry off should not regress match rate below the existing baseline (the new id_type is purely additive).

### Manual smoke check

After implementation, run `/check-gene-ids` against biller 2022 with the new id_type wired in; confirm match rate is ≥ 95% for both strains.

## Acceptance criteria

1. `extract_uniprot_annotation_tokens` returns the documented `(token, id_type)` pairs on the 6 test inputs above.
2. `build_gene_id_mapping.py` step 3 with biller 2022 in the registry produces a `gene_id_mapping.json` for MIT9312 and MIT9313 that contains:
   - `Q31DF2_PROM9` (and similar entry names) in `specific_lookup` mapping to a locus_tag
   - GN=-derived tokens registered as gene_name (multi_lookup)
   - No pollution: `specific_lookup` size for MIT9312 grows by O(rows) entries, not by O(words-in-narrative).
3. `resolve_paper_ids` step 4 against biller 2022 S4 produces `_resolved.csv` files where `locus_tag` is non-null on ≥ 95% of rows for both strains.
4. `pytest -m "not slow and not kg"` passes including the new unit tests.
5. Biller 2022 paperconfig validator passes with the updated `id_columns`.

## Out of scope

- No generalisation to other embedded-token patterns (`cds-`, `lcl|`, paren-alt, etc.). User-confirmed narrow scope (option A in brainstorming).
- No deprecation of `scripts/build_modified_csv/` builders — those handle distinct cases (column derivation beyond regex) and remain useful.
- No changes to `omics_adapter.py` — it consumes `_resolved.csv`, which remains the same shape.
- No new tests for the legacy `build_id_lookup` API — the new path is v2-only (`resolve_row` + `MappingData`).
- No retroactive reprocessing of older paperconfigs to use the new id_type. Existing modified-CSV workflows continue to work; new papers may opt in.
