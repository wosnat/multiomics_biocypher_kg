# `uniprot_annotation_string` id_type — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new `uniprot_annotation_string` id_type that auto-extracts the UniProt entry name (Tier 1) and `GN=` token (Tier 3) from FASTA-header-style annotation columns, replacing the current per-paper `_modified.csv` workaround.

**Architecture:** A single shared extractor (`extract_uniprot_annotation_tokens`) lives in `multiomics_kg/utils/gene_id_utils.py`. Step 3 (`build_gene_id_mapping.py`) calls it via a new `_emit_pairs` helper at row-extraction time so the mapping is enriched with the right tokens at the right tier. Step 4 (`resolve_row`) calls the same extractor at lookup time so paper rows resolve via the existing Pass 1/1b/2/3 logic against the extracted tokens. Validator + paperconfig skill table get the new id_type.

**Tech Stack:** Python 3, pandas, pytest. Project uses `uv` for dependency management; tests gated by `pytest -m "not slow and not kg"`.

**Spec:** [docs/superpowers/specs/2026-04-25-uniprot-annotation-id-type-design.md](docs/superpowers/specs/2026-04-25-uniprot-annotation-id-type-design.md)

---

## File Structure

| File | Responsibility | Action |
|------|----------------|--------|
| `multiomics_kg/utils/gene_id_utils.py` | Shared ID utilities — extractor + `resolve_row` | Modify (add extractor; thread id_type into `resolve_row` per-column passes) |
| `multiomics_kg/download/build_gene_id_mapping.py` | Step 3 — build per-strain `gene_id_mapping.json` | Modify (introduce `_emit_pairs`; replace 6 append sites across 4 extraction functions) |
| `scripts/validate_paperconfig.py` | Paperconfig validator | Modify (add `uniprot_annotation_string` to `VALID_ID_TYPES`) |
| `tests/test_gene_id_utils.py` | Unit tests for shared utilities | Modify (extractor tests + `resolve_row` extractor-aware tests) |
| `tests/test_build_gene_id_mapping.py` | Tests for step-3 extraction | Modify or create (one test asserting the new id_type expands at row-extraction time) |
| `.claude/skills/paperconfig/SKILL.md` | Paperconfig authoring guide | Modify (ID Types table — new row) |
| `data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml` | Driver paper config | Modify (add `UniProt Annotation` to `id_columns` of both `s4_*_enrichment` entries) |

Each task below is self-contained: failing test → minimal implementation → passing test → commit.

---

## Task 1: Extractor function `extract_uniprot_annotation_tokens`

**Files:**
- Modify: `multiomics_kg/utils/gene_id_utils.py` (add a function near top of file, after the `expand_list` helper at ~L286)
- Test: `tests/test_gene_id_utils.py` (add a new test class `TestExtractUniprotAnnotationTokens`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_gene_id_utils.py` (after the existing tests, before any module-level code if present at end of file):

```python
class TestExtractUniprotAnnotationTokens:
    """Unit tests for extract_uniprot_annotation_tokens (Biller 2022 driver case)."""

    def test_full_input_entry_and_gn_locus_tag(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = ("Q31DF2_PROM9 Type II secretion system protein-like protein "
             "OS=Prochlorococcus marinus (strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")
        assert extract_uniprot_annotation_tokens(s) == [
            ("Q31DF2_PROM9", "uniprot_entry_name"),
            ("PMT9312_0032", "gene_name"),
        ]

    def test_full_input_entry_and_gn_gene_symbol(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = ("RL33_PROM9 50S ribosomal protein L33 "
             "OS=Prochlorococcus marinus (strain MIT 9312) GN=rpmG PE=3 SV=2")
        assert extract_uniprot_annotation_tokens(s) == [
            ("RL33_PROM9", "uniprot_entry_name"),
            ("rpmG", "gene_name"),
        ]

    def test_entry_only_no_gn(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = "Q7V5C8_PROMM Possible serine protease OS=Prochlorococcus marinus"
        assert extract_uniprot_annotation_tokens(s) == [
            ("Q7V5C8_PROMM", "uniprot_entry_name"),
        ]

    def test_gn_only_no_entry_prefix(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        # No leading <entry>_<ORG> shape; only GN= token
        s = "putative serine protease GN=PMT_1636 PE=4 SV=1"
        assert extract_uniprot_annotation_tokens(s) == [
            ("PMT_1636", "gene_name"),
        ]

    def test_neither_pattern_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        assert extract_uniprot_annotation_tokens("hypothetical protein") == []

    def test_empty_and_non_string_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        assert extract_uniprot_annotation_tokens("") == []
        assert extract_uniprot_annotation_tokens("   ") == []
        assert extract_uniprot_annotation_tokens(None) == []  # type: ignore[arg-type]
        assert extract_uniprot_annotation_tokens(123) == []   # type: ignore[arg-type]

    def test_first_gn_wins_when_multiple(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = "X_Y abcd GN=first PE=1 GN=second"
        result = extract_uniprot_annotation_tokens(s)
        # First entry-name match + first GN= match
        assert ("first", "gene_name") in result
        assert ("second", "gene_name") not in result
```

- [ ] **Step 2: Run the failing tests**

```bash
uv run pytest tests/test_gene_id_utils.py::TestExtractUniprotAnnotationTokens -v
```

Expected: all 7 tests fail with `ImportError: cannot import name 'extract_uniprot_annotation_tokens'`.

- [ ] **Step 3: Implement the extractor**

In `multiomics_kg/utils/gene_id_utils.py`, locate the existing `expand_list` function (around line 286). Add the following block immediately after `expand_list` ends (just before `_heuristic_candidates`):

```python
# ─── UniProt FASTA-header style annotation parsing ────────────────────────────

# Matches the leading "<entry>_<ORG>" token that begins UniProt-style annotation
# strings (e.g. "Q31DF2_PROM9", "RL33_PROM9", "Q7V5C8_PROMM").
_UNIPROT_ANNOT_ENTRY_RE = re.compile(r"^([A-Z0-9]+_[A-Z0-9]+)\b")

# Matches the GN=<token> capture (gene name OR locus tag, depending on entry).
_UNIPROT_ANNOT_GN_RE = re.compile(r"\bGN=(\S+)")


def extract_uniprot_annotation_tokens(value) -> list[tuple[str, str]]:
    """Parse a UniProt FASTA-header style annotation string into (token, id_type) pairs.

    The leading "<entry>_<ORG>" token is emitted as `uniprot_entry_name` (Tier 1).
    The GN=<token> capture is emitted as `gene_name` (Tier 3); when GN= happens to
    hold a real locus_tag, the existing Tier 1 specific_lookup catches it first.

    Returns an empty list when neither pattern matches (e.g. plain product
    descriptions, blank cells, non-string inputs); the cell is then a no-op for
    ID purposes.
    """
    out: list[tuple[str, str]] = []
    if not isinstance(value, str):
        return out
    s = value.strip()
    if not s:
        return out
    if (m := _UNIPROT_ANNOT_ENTRY_RE.match(s)):
        out.append((m.group(1), "uniprot_entry_name"))
    if (m := _UNIPROT_ANNOT_GN_RE.search(s)):
        out.append((m.group(1), "gene_name"))
    return out
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest tests/test_gene_id_utils.py::TestExtractUniprotAnnotationTokens -v
```

Expected: all 7 tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/gene_id_utils.py tests/test_gene_id_utils.py
git commit -m "$(cat <<'EOF'
gene_id_utils: extract_uniprot_annotation_tokens for FASTA-header columns

Parses "<entry>_<ORG> ... GN=<token> ..." into (uniprot_entry_name, Tier 1)
and (gene_name, Tier 3) pairs. Driver: Biller 2022 vesicle proteomics S4.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Whitelist new id_type in validator

**Files:**
- Modify: `scripts/validate_paperconfig.py:92-98` — `VALID_ID_TYPES` set

- [ ] **Step 1: Write the failing test**

Add to `tests/test_gene_id_utils.py` (a separate test class — keeps the validator concern testable from inside the existing test suite):

```python
class TestUniprotAnnotationStringValidatorAccepts:
    """Validator must accept the new id_type; the free-text guard must NOT fire."""

    def test_id_type_in_valid_set(self):
        from scripts.validate_paperconfig import VALID_ID_TYPES
        assert "uniprot_annotation_string" in VALID_ID_TYPES
```

- [ ] **Step 2: Run the failing test**

```bash
uv run pytest tests/test_gene_id_utils.py::TestUniprotAnnotationStringValidatorAccepts -v
```

Expected: FAIL — `'uniprot_annotation_string' not in VALID_ID_TYPES`.

- [ ] **Step 3: Add the new id_type**

In `scripts/validate_paperconfig.py`, line 92-98, edit `VALID_ID_TYPES`:

```python
VALID_ID_TYPES = {
    "locus_tag", "locus_tag_ncbi", "locus_tag_cyanorak",
    "old_locus_tag", "alternative_locus_tag",
    "gene_name", "gene_synonym",
    "protein_id_refseq", "uniprot_accession", "uniprot_entry_name",
    "uniprot_annotation_string",  # NEW — extracts entry_name + GN= via gene_id_utils
    "jgi_id", "probeset", "rast_id", "annotation_specific", "other",
}
```

- [ ] **Step 4: Run the test to confirm it passes**

```bash
uv run pytest tests/test_gene_id_utils.py::TestUniprotAnnotationStringValidatorAccepts -v
```

Expected: PASS.

- [ ] **Step 5: Verify the existing free-text-column guard still does NOT fire on the new id_type**

The guard at `scripts/validate_paperconfig.py:342` reads:

```python
if id_type == "locus_tag" and col_name in FREE_TEXT_COLUMN_NAMES:
```

Confirm by running the full validator against `data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml` after adding the new id_type to its `s4_*` entries (Task 6); for now, just `grep` confirm:

```bash
uv run python -c "from scripts.validate_paperconfig import VALID_ID_TYPES, FREE_TEXT_COLUMN_NAMES; print('uniprot_annotation_string' in VALID_ID_TYPES); print(FREE_TEXT_COLUMN_NAMES)"
```

Expected: `True` printed, then the existing free-text column set (unchanged).

- [ ] **Step 6: Commit**

```bash
git add scripts/validate_paperconfig.py tests/test_gene_id_utils.py
git commit -m "$(cat <<'EOF'
validate_paperconfig: accept uniprot_annotation_string id_type

Free-text-column guard at L342 only fires for id_type=locus_tag, so the
new id_type is intentionally exempt — it is meant for narrative columns.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Step 3 — wire extractor into `build_gene_id_mapping.py`

**Files:**
- Modify: `multiomics_kg/download/build_gene_id_mapping.py` — introduce `_emit_pairs` helper; replace 6 `row_pairs.append((val, id_type))` sites (lines ~116, 337, 345, 399, 407, 468, 475 — exact line numbers shift as imports/helpers are added; locate by surrounding context)
- Test: `tests/test_build_gene_id_mapping.py` (create new file or extend existing — search first)

- [ ] **Step 1: Locate or create the test file**

Run:

```bash
uv run pytest --collect-only tests/ 2>&1 | grep -i "build_gene_id_mapping" | head -5
```

If a base test file `tests/test_build_gene_id_mapping.py` does not exist, create it with the standard header. Existing file `tests/test_build_gene_id_mapping_derived_metrics.py` covers derived_metrics extraction; we add a new file for the cross-cutting `_emit_pairs` behavior.

- [ ] **Step 2: Write failing tests**

Create `tests/test_build_gene_id_mapping_uniprot_annotation.py`:

```python
"""Tests for uniprot_annotation_string id_type expansion in step-3 row extraction.

The extractor lives in gene_id_utils; this file verifies that
build_gene_id_mapping calls it at the right place for each of the four
extraction functions (id_translation, csv, gene_clusters, derived_metrics).
"""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import pytest

from multiomics_kg.download.build_gene_id_mapping import (
    extract_rows_from_csv_table,
    extract_rows_from_derived_metrics_table,
    extract_rows_from_gene_clusters_table,
    extract_rows_from_id_translation,
)


def _write_csv(path: Path, rows: list[dict], header: list[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)


@pytest.fixture
def annot_string():
    return ("Q31DF2_PROM9 Type II secretion system protein-like protein "
            "OS=Prochlorococcus marinus (strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")


def _flatten_pairs(rows):
    """Helper: flatten (row_pairs, source_label) tuples into a single list of pairs."""
    return [pair for row_pairs, _ in rows for pair in row_pairs]


def test_id_translation_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "id_translation.csv"
    _write_csv(csv_path, [{"Locus": "PMT9312_0032", "Annot": annot_string}],
               ["Locus", "Annot"])
    entry = {
        "type": "id_translation",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "id_columns": [
            {"column": "Locus", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    # build_gene_id_mapping reads paths relative to PROJECT_ROOT; patch via cwd
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_id_translation(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs
    # The raw long string must NOT be emitted with id_type=uniprot_annotation_string
    assert not any(t == "uniprot_annotation_string" for _, t in pairs)


def test_id_translation_no_match_emits_no_extra_pairs(tmp_path):
    csv_path = tmp_path / "id_translation.csv"
    _write_csv(csv_path, [{"Locus": "PMT9312_0032", "Annot": "hypothetical protein"}],
               ["Locus", "Annot"])
    entry = {
        "type": "id_translation",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "id_columns": [
            {"column": "Locus", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_id_translation(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    # Locus pair still emitted; Annot column contributes nothing
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert all(t in {"locus_tag"} for _, t in pairs)


def test_derived_metrics_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "dm.csv"
    _write_csv(csv_path,
               [{"Protein ID": "Q31DF2", "UniProt Annotation": annot_string,
                 "value_col": "1.5"}],
               ["Protein ID", "UniProt Annotation", "value_col"])
    entry = {
        "type": "derived_metrics_table",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "name_col": "Protein ID",
        "id_columns": [
            {"column": "Protein ID", "id_type": "uniprot_accession"},
            {"column": "UniProt Annotation", "id_type": "uniprot_annotation_string"},
        ],
        "metrics": [],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_derived_metrics_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("Q31DF2", "uniprot_accession") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs


def test_csv_table_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "csv_table.csv"
    _write_csv(csv_path,
               [{"Protein ID": "Q31DF2", "UniProt Annotation": annot_string,
                 "logfc": "1.5"}],
               ["Protein ID", "UniProt Annotation", "logfc"])
    entry = {
        "type": "csv",
        "filename": str(csv_path.relative_to(tmp_path)),
        "id_columns": [
            {"column": "Protein ID", "id_type": "uniprot_accession"},
            {"column": "UniProt Annotation", "id_type": "uniprot_annotation_string"},
        ],
        "statistical_analyses": [
            {"id": "test_de", "experiment": "exp", "name_col": "Protein ID",
             "logfc_col": "logfc"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_csv_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs


def test_gene_clusters_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "clusters.csv"
    _write_csv(csv_path,
               [{"gene_id": "PMT9312_0032", "Annot": annot_string, "cluster": "C1"}],
               ["gene_id", "Annot", "cluster"])
    entry = {
        "type": "gene_clusters",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "gene_id_col": "gene_id",
        "cluster_col": "cluster",
        "id_columns": [
            {"column": "gene_id", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_gene_clusters_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs
```

- [ ] **Step 3: Run failing tests**

```bash
uv run pytest tests/test_build_gene_id_mapping_uniprot_annotation.py -v
```

Expected: 5 tests fail — extractor not wired into the four extraction functions yet, so the long annotation string is appended verbatim with id_type `uniprot_annotation_string`.

- [ ] **Step 4: Add `_emit_pairs` helper at top of `build_gene_id_mapping.py`**

In `multiomics_kg/download/build_gene_id_mapping.py`, immediately above the `_safe_str` helper (currently line 71), add:

```python
from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens


def _emit_pairs(val: str, id_type: str) -> list[tuple[str, str]]:
    """Convert a single (cell_value, declared_id_type) pair into a list of
    (token, id_type) pairs for the gene_id mapping graph.

    For most id_types this is a 1-element list. For id_types that wrap
    multiple embedded tokens (currently only `uniprot_annotation_string`),
    returns the list of extracted (token, canonical_id_type) pairs; an
    empty list when nothing parseable is found, so the cell is skipped.
    """
    if id_type == "uniprot_annotation_string":
        return extract_uniprot_annotation_tokens(val)
    return [(val, id_type)]
```

- [ ] **Step 5: Replace 6 append sites with `_emit_pairs`**

In `multiomics_kg/download/build_gene_id_mapping.py`, replace each of the six existing `row_pairs.append((val, ...))` calls in the four extraction functions:

**Site 1** — `extract_rows_from_id_translation`, around line 116:

```python
# BEFORE:
val = _safe_str(row.get(col, ""))
if val:
    row_pairs.append((val, id_type))

# AFTER:
val = _safe_str(row.get(col, ""))
if val:
    row_pairs.extend(_emit_pairs(val, id_type))
```

**Site 2** — `extract_rows_from_csv_table`, name_col block (around line 337):

```python
# BEFORE:
row_pairs.append((val, nc_id_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, nc_id_type))
```

**Site 3** — `extract_rows_from_csv_table`, id_columns block (around line 345):

```python
# BEFORE:
row_pairs.append((val, id_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, id_type))
```

**Site 4** — `extract_rows_from_gene_clusters_table`, gene_id_col block (around line 399):

```python
# BEFORE:
row_pairs.append((val, gid_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, gid_type))
```

**Site 5** — `extract_rows_from_gene_clusters_table`, id_columns block (around line 407):

```python
# BEFORE:
row_pairs.append((val, id_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, id_type))
```

**Site 6** — `extract_rows_from_derived_metrics_table`, name_col block (around line 468):

```python
# BEFORE:
row_pairs.append((val, name_col_id_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, name_col_id_type))
```

**Site 7** — `extract_rows_from_derived_metrics_table`, id_columns block (around line 475):

```python
# BEFORE:
row_pairs.append((val, id_type))

# AFTER:
row_pairs.extend(_emit_pairs(val, id_type))
```

(Seven sites total — five in the body and the two name_col/gene_id_col sites in csv + gene_clusters + derived_metrics. Search for `row_pairs.append((val,` to confirm you've caught all of them; there should be zero remaining matches after this edit.)

Verify with:

```bash
uv run python -c "
import re, pathlib
src = pathlib.Path('multiomics_kg/download/build_gene_id_mapping.py').read_text()
appends = re.findall(r'row_pairs\.append\(\(val,', src)
extends = re.findall(r'row_pairs\.extend\(_emit_pairs\(val,', src)
print(f'appends={len(appends)}  extends={len(extends)}')
"
```

Expected: `appends=0  extends=7`. (Note the GFF function `extract_rows_from_annotation_gff` appends differently — it uses static GFF_ATTR_TYPES values, not paperconfig-declared id_type, so it's intentionally not changed.)

- [ ] **Step 6: Run the unit tests to confirm they pass**

```bash
uv run pytest tests/test_build_gene_id_mapping_uniprot_annotation.py -v
```

Expected: all 5 tests pass.

- [ ] **Step 7: Run the broader unit suite for regressions**

```bash
uv run pytest -m "not slow and not kg" tests/test_gene_id_utils.py tests/test_build_gene_id_mapping_uniprot_annotation.py tests/test_build_gene_id_mapping_derived_metrics.py -v
```

Expected: all pass; no regressions in existing extraction tests.

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/download/build_gene_id_mapping.py tests/test_build_gene_id_mapping_uniprot_annotation.py
git commit -m "$(cat <<'EOF'
build_gene_id_mapping: expand uniprot_annotation_string at row extraction

New _emit_pairs helper routes uniprot_annotation_string through the
extractor; all 4 extraction functions (id_translation, csv, gene_clusters,
derived_metrics) call it via row_pairs.extend(_emit_pairs(...)).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Step 4 — wire extractor into `resolve_row`

**Files:**
- Modify: `multiomics_kg/utils/gene_id_utils.py` — `resolve_row` function (lines 347–477) — thread id_type per-column lookup through each pass
- Test: `tests/test_gene_id_utils.py` — add a new test class `TestResolveRowUniprotAnnotation`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_gene_id_utils.py`:

```python
class TestResolveRowUniprotAnnotation:
    """resolve_row must extract uniprot_annotation_string columns at lookup time."""

    @pytest.fixture
    def mapping_data(self):
        from multiomics_kg.utils.gene_id_utils import MappingData
        # PMT9312_0032 is a real locus tag; specific_lookup also has the
        # uniprot_entry_name token Q31DF2_PROM9.
        return MappingData(
            locus_tags={"PMT9312_0032"},
            specific_lookup={
                "Q31DF2_PROM9": "PMT9312_0032",
                "PMT9312_0032": "PMT9312_0032",
            },
            multi_lookup={"rpmG": ["PMT9312_0107"]},
            conflicts={},
        )

    def test_resolves_via_entry_name(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        row = {"Annot": ("Q31DF2_PROM9 protein OS=... "
                        "GN=NOT_A_KNOWN_TAG PE=4 SV=1")}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"
        assert "Annot" in method  # method string mentions the column

    def test_resolves_via_gn_locus_tag_when_entry_missing(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        # No leading entry name → only GN= token; that token is itself a real
        # locus_tag (Tier 1 specific_lookup hit).
        row = {"Annot": "putative protein GN=PMT9312_0032 PE=3 SV=1"}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, _ = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"

    def test_unresolvable_when_no_tokens_match(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        row = {"Annot": "ZZZZZ_NOTREAL protein OS=Foo GN=NOT_A_GENE PE=4 SV=1"}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt is None
        assert method == "unresolved"

    def test_skips_cell_with_no_extractable_tokens(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        # When the cell has no entry/GN= tokens, that column contributes
        # nothing — but other columns can still resolve.
        row = {
            "Annot": "hypothetical protein",
            "Other": "PMT9312_0032",
        }
        id_columns = [
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
            {"column": "Other", "id_type": "locus_tag"},
        ]
        lt, _ = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"

    def test_does_not_match_long_annotation_as_literal(self, mapping_data):
        """Regression guard: the raw long annotation string must NOT be looked up
        verbatim — that would either miss (best case) or pollute the fallback
        path (if the literal happens to coincide with a multi_lookup key)."""
        from multiomics_kg.utils.gene_id_utils import resolve_row
        long_annot = ("Q31DF2_PROM9 protein OS=Prochlorococcus marinus "
                      "(strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")
        row = {"Annot": long_annot}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        # Resolves via the extracted entry_name, not the literal long string
        assert lt == "PMT9312_0032"
```

- [ ] **Step 2: Run the failing tests**

```bash
uv run pytest tests/test_gene_id_utils.py::TestResolveRowUniprotAnnotation -v
```

Expected: tests fail because `resolve_row` currently calls `expand_list(str(raw))` on the literal long string, which produces no usable tokens.

- [ ] **Step 3: Add an extractor-aware helper inside `resolve_row`**

In `multiomics_kg/utils/gene_id_utils.py`, locate `resolve_row` (line 347). Before the existing `# ── Pass 1` block, add:

```python
    # Build per-column id_type index so passes can apply id_type-specific
    # extraction (currently only `uniprot_annotation_string`).
    id_type_by_col: dict[str, str] = {}
    for c in id_columns:
        col_name = c.get("column", "")
        if col_name:
            id_type_by_col[col_name] = c.get("id_type", "other")

    def _candidate_values(col: str, raw) -> list[str]:
        """Return candidate string values for a column, applying any
        id_type-specific extraction. Default: expand_list on the raw cell."""
        if isinstance(raw, float) and pd.isna(raw):
            return []
        if raw is None:
            return []
        s = str(raw)
        if id_type_by_col.get(col) == "uniprot_annotation_string":
            return [tok for tok, _ in extract_uniprot_annotation_tokens(s)]
        return expand_list(s)
```

Then replace each occurrence of `for val in expand_list(str(raw)):` inside the four passes (Pass 1, 1b, 2, 3, 3b — five sites) with `for val in _candidate_values(col, raw):` and remove the now-redundant `if raw is None or (isinstance(raw, float) and pd.isna(raw)):` skip block (the helper handles it). Concretely, each pass changes from:

```python
# BEFORE (one of the five passes):
for col in all_cols:
    raw = row.get(col)
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        continue
    for val in expand_list(str(raw)):
        if not val:
            continue
        # ... existing lookup logic ...

# AFTER:
for col in all_cols:
    raw = row.get(col)
    for val in _candidate_values(col, raw):
        if not val:
            continue
        # ... existing lookup logic UNCHANGED ...
```

The five passes are at:
- Pass 1: ~line 387–401 (`# ── Pass 1: specific_lookup with list expansion`)
- Pass 1b: ~line 402–415 (`# ── Pass 1b: case-insensitive fallback`)
- Pass 2: ~line 417–427 (`# ── Pass 2: heuristics`)
- Pass 3: ~line 429–442 (`# ── Pass 3: multi_lookup`)
- Pass 3b: ~line 444–457 (`# ── Pass 3b: case-insensitive fallback on multi_lookup`)

The trailing diagnostic block (lines 459–477) does its own row.get(name_col)/expand_list call — leave it unchanged; it only runs after all passes have failed and uses `name_col` directly, so the literal long annotation in the failure path is harmless (it would already have produced no extracted tokens via the new path, leaving `parts` from `expand_list` benign).

- [ ] **Step 4: Run the tests to confirm they pass**

```bash
uv run pytest tests/test_gene_id_utils.py::TestResolveRowUniprotAnnotation -v
```

Expected: all 5 tests pass.

- [ ] **Step 5: Run the full gene_id_utils suite to check for regressions**

```bash
uv run pytest tests/test_gene_id_utils.py -v
```

Expected: all tests pass — pre-existing tests should be unaffected since the helper falls back to `expand_list` for any column whose declared id_type is not `uniprot_annotation_string`.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/gene_id_utils.py tests/test_gene_id_utils.py
git commit -m "$(cat <<'EOF'
resolve_row: id_type-aware candidate expansion for uniprot_annotation_string

Per-column dispatch via _candidate_values: when id_type is
uniprot_annotation_string, replace expand_list with the extractor;
all other id_types are unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Document new id_type in paperconfig SKILL.md

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md` — `ID Types` table

- [ ] **Step 1: Locate the ID Types table**

Run:

```bash
grep -n "id_type.*Description.*Example" .claude/skills/paperconfig/SKILL.md | head -3
```

Find the `| `id_type` | Description | Example |` table (around line 477 in the current revision; locate by content search above).

- [ ] **Step 2: Add the new row**

Add a new row for `uniprot_annotation_string` immediately after the `uniprot_entry_name` row, e.g.:

```markdown
| `uniprot_entry_name` | UniProt entry name (build script strips `_ORGANISM` suffix) | `DNAA_PROM0` |
| `uniprot_annotation_string` | UniProt FASTA-header style annotation; build pipeline auto-extracts the `<entry>_<ORG>` token (Tier 1) and `GN=<token>` (Tier 3). Use on narrative columns named "UniProt Annotation", "Description", etc. that embed these markers. A column declared this way may also appear in `product_columns` for description harvest. | `Q31DF2_PROM9 ... GN=PMT9312_0032 PE=4 SV=1` |
| `jgi_id` | JGI IMG gene catalog ID (integer string) | `2626311743` |
```

(Insert between `uniprot_entry_name` and the existing `jgi_id` row — exact existing text to match in the diff is the `uniprot_entry_name` line plus the next line.)

- [ ] **Step 3: Visual verify the table renders**

```bash
grep -A1 "uniprot_annotation_string" .claude/skills/paperconfig/SKILL.md | head -2
```

Expected: the new row prints; no markdown table breakage.

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md
git commit -m "$(cat <<'EOF'
paperconfig SKILL: document uniprot_annotation_string id_type

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Migrate Biller 2022 paperconfig to new id_type

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml` — `s4_mit9312_enrichment` and `s4_mit9313_enrichment` entries

- [ ] **Step 1: Add `UniProt Annotation` to `id_columns` in both S4 entries**

In `data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml`, locate `s4_mit9312_enrichment.id_columns` and add a second entry:

```yaml
    s4_mit9312_enrichment:
      type: derived_metrics_table
      filename: data/Prochlorococcus/papers_and_supp/biller 2022/emi15834-sup-0005-tables4_MIT9312_modified.csv
      organism: Prochlorococcus MIT9312
      experiment: vesicle_proteomics_mit9312
      name_col: Protein ID
      id_columns:
        - column: Protein ID
          id_type: uniprot_accession
        - column: UniProt Annotation        # NEW
          id_type: uniprot_annotation_string
      product_columns:
        - column: UniProt Annotation         # unchanged — still in product_columns
        - column: Prokka Annotation
      ...
```

Apply the same change to `s4_mit9313_enrichment.id_columns` (same column name, same id_type).

- [ ] **Step 2: Re-run validator**

```bash
uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml"
```

Expected: `VALIDATION PASSED` (the same two warnings about no DE analyses are still expected — they pre-date this change).

- [ ] **Step 3: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/biller 2022/paperconfig.yaml"
git commit -m "$(cat <<'EOF'
biller 2022: declare UniProt Annotation as uniprot_annotation_string

Removes the need for a per-paper modified-CSV builder. The column is
declared in both id_columns (for token extraction) and product_columns
(for description harvest) — they are processed independently.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: End-to-end verification — rebuild MIT9312/MIT9313 ID mappings + check match rates

**Files:**
- No file changes — pipeline run + diagnostic check.

- [ ] **Step 1: Rebuild the gene_id mapping for MIT9312 and MIT9313**

```bash
bash scripts/prepare_data.sh --steps 3 4 --strains MIT9312 MIT9313 --force
```

Expected: writes `cache/data/Prochlorococcus/genomes/MIT9312/gene_id_mapping.json` and `MIT9313/gene_id_mapping.json`; logs in `logs/prepare_data_step3.log` and `logs/prepare_data_step4.log`. Pipeline should finish without errors.

- [ ] **Step 2: Sanity-check that the new tokens landed in the mapping**

```bash
uv run python <<'PY'
import json
for strain in ("MIT9312", "MIT9313"):
    p = f"cache/data/Prochlorococcus/genomes/{strain}/gene_id_mapping.json"
    m = json.load(open(p))
    sl = m["specific_lookup"]
    ml = m["multi_lookup"]
    # Sample tokens we expect to see for MIT9312
    if strain == "MIT9312":
        print(f"{strain}: Q31DF2_PROM9 in specific_lookup ->",
              sl.get("Q31DF2_PROM9"))
        print(f"{strain}: rpmG in multi_lookup ->", ml.get("rpmG"))
    else:
        print(f"{strain}: Q7V5C8_PROMM in specific_lookup ->",
              sl.get("Q7V5C8_PROMM"))
    print(f"{strain}: specific_lookup size = {len(sl)}; multi_lookup size = {len(ml)}")
PY
```

Expected: `Q31DF2_PROM9` resolves to a real locus_tag (e.g. `PMT9312_0032`); `Q7V5C8_PROMM` resolves to `PMT_1636`; sizes are within ~1% of pre-change baseline (no pollution). If `specific_lookup` size has grown by thousands, that signals the regex captured something unintended — investigate before continuing.

- [ ] **Step 3: Verify resolution rates on the S4 modified CSVs**

```bash
uv run python <<'PY'
import pandas as pd
for strain in ("MIT9312", "MIT9313"):
    p = (f"data/Prochlorococcus/papers_and_supp/biller 2022/"
         f"emi15834-sup-0005-tables4_{strain}_modified_resolved.csv")
    df = pd.read_csv(p)
    total = len(df)
    matched = df["locus_tag"].notna().sum()
    rate = matched / total * 100 if total else 0
    print(f"{strain}: {matched}/{total} resolved = {rate:.1f}%")
    # Per resolution method
    print(df["resolution_method"].value_counts().to_dict())
PY
```

Expected: ≥ 95% resolution for both strains. Methods should include `tier1:Protein ID` (UniProt accession hits) and `tier1:UniProt Annotation` (entry name hits) as the dominant paths.

- [ ] **Step 4: Run check-gene-ids skill for a final cross-reference**

```bash
# This invokes the diagnostic skill — manual inspection of output below.
uv run python -c "
from multiomics_kg.utils.gene_id_utils import load_mapping_v2, resolve_row
import pandas as pd
md = load_mapping_v2('cache/data/Prochlorococcus/genomes/MIT9312')
df = pd.read_csv('data/Prochlorococcus/papers_and_supp/biller 2022/emi15834-sup-0005-tables4_MIT9312_modified.csv')
unresolved = []
id_columns = [
    {'column': 'Protein ID', 'id_type': 'uniprot_accession'},
    {'column': 'UniProt Annotation', 'id_type': 'uniprot_annotation_string'},
]
for _, row in df.iterrows():
    lt, method = resolve_row(row.to_dict(), 'Protein ID', id_columns, md)
    if lt is None:
        unresolved.append((row.get('Protein ID'), row.get('UniProt Annotation', '')[:80]))
print(f'unresolved: {len(unresolved)} / {len(df)}')
for u in unresolved[:5]:
    print(' ', u)
"
```

Expected: small number unresolved (<5%); printed examples are real exceptions (proteins not in the strain's gene_annotations_merged.json).

- [ ] **Step 5: No commit — verification only.** If steps 1–4 all show expected results, the implementation is complete. If anything is off (lower match rate, mapping pollution, pipeline errors), reopen the relevant earlier task to debug before declaring done.

---

## Self-Review

**Spec coverage (each spec section → task):**

| Spec section | Implemented in |
|---|---|
| 1. Paperconfig surface | Task 6 |
| 2. Extractor function | Task 1 |
| 3. Step 3 integration | Task 3 |
| 3. Step 4 integration | Task 4 |
| 3. Validator | Task 2 |
| 3. Skill / docs | Task 5 |
| 4. Tests (unit + integration) | Tasks 1, 3, 4 (unit), Task 7 (e2e) |
| 5. Acceptance criteria | Task 7 verifies criteria 1–5 |

**Placeholder scan:** none — every step has exact paths, full code, exact commands, and expected output.

**Type/name consistency:**
- `extract_uniprot_annotation_tokens` — defined Task 1, used in Tasks 3 (`_emit_pairs`) and 4 (`_candidate_values`). ✓
- `_emit_pairs(val, id_type) -> list[tuple[str,str]]` — defined Task 3, used 7 sites. ✓
- `id_type_by_col`, `_candidate_values(col, raw)` — defined Task 4. ✓
- `VALID_ID_TYPES` — modified Task 2; matches the import in Task 2 test. ✓
- Method-string convention `tier1:<col>` — preserved in Task 4 because the existing pass logic is unchanged after token substitution. ✓
