# Paperconfig DOI Override — Design Spec

**Date:** 2026-04-16
**Status:** Draft

## Problem

Publication node IDs and Experiment node IDs derive from the DOI extracted by `pdf_publication_extraction.py`. When the PDF extractor fails to find a DOI, the Publication ID falls back to `pub_<title_first_20_chars>`, producing ugly IDs and a failing KG validity test (`test_publication_doi_or_pmid_present`).

Concrete case: Fuszard 2012 extracts to `pub_Comparative quantitative prote` instead of `10.1186/2046-9063-8-7`.

## Solution

Add an optional `doi:` field to the `publication:` block of paperconfig.yaml. When present, it overrides the PDF-extracted DOI for:

- Publication node `doi` property
- Publication node ID (`add_prefix_to_id(prefix="doi", identifier=doi)`)
- Experiment node IDs (`{doi}_{experiment_key}`)

## Decisions

| Question | Decision |
|---|---|
| Precedence when both config and PDF have a DOI | Config wins; log warning if they disagree |
| Add `pmid:` override too? | No — YAGNI |
| Back-fill title/authors from config? | No — different PR |
| Behavior when PDF extraction fails entirely | Hard error (raise), replacing today's silent skip |
| Resource-only configs (no `publication:` block) | Unaffected — `download_data()` returns early before PDF check |

## Paperconfig Schema Change

```yaml
publication:
  papername: "Fuszard 2012"
  doi: "10.1186/2046-9063-8-7"   # optional, overrides PDF-extracted DOI
  papermainpdf: "data/.../paper.pdf"
```

`doi:` must match `^10\.\d{4,9}/\S+$` when present.

## Adapter Changes (`omics_adapter.py`)

### New helper: `_get_override_doi()`

Reads `config_data["publication"]["doi"]`, returns stripped string or `None`.

### `get_publication_id()` — override precedence

1. Read `_get_override_doi()`.
2. If override and PDF-extracted DOI both present and differ → `logger.warning(...)`.
3. If override present → return it.
4. Otherwise, fall through to existing logic (PDF publication_id → doi → pubmed_id → title fallback).

### `get_publication_nodes()` — `doi` property

```python
"doi": self._get_override_doi() or pub.get("doi"),
```

All other Publication properties (title, authors, journal, abstract, etc.) remain sourced from PDF extraction.

### `download_data()` — hard error on PDF extraction failure

Replace the existing `logger.warning(...) Skipping publication nodes` with:

- `FileNotFoundError` if `papermainpdf` path is missing or doesn't exist.
- `RuntimeError` if `pdf_extractor.extract_from_pdf()` returns empty or has no `"publication"` key.

Resource-only configs (no `publication:` block) are unaffected — `download_data()` returns early at line 254.

### Experiment IDs — no change needed

`f"{pub_id}_{exp_key}"` already calls `get_publication_id()`.

## Validator Changes (`scripts/validate_paperconfig.py`)

Add after the `papername` check:

- Module-level: `DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")`
- If `pub.get("doi")` is not None: validate it matches `DOI_RE`. Error if invalid, print if valid.

## Tests

### Validator test (`tests/test_paperconfig_validation.py`)

- Valid DOI passes.
- Missing `doi:` key passes (optional).
- Invalid DOIs (`"not-a-doi"`, `"https://doi.org/10.1186/..."`) produce errors.

### Adapter test (new or existing omics adapter test file)

- Paperconfig with `doi: "10.9999/testdoi"` + one experiment.
- Mock `pdf_extractor.extract_from_pdf` returning a publication block with a different DOI.
- Assert: Publication node ID is `doi:10.9999/testdoi`, `doi` property is `"10.9999/testdoi"`, Experiment IDs start with `"10.9999/testdoi_"`.
- Assert: warning logged about disagreement.
- Second case: PDF extraction returns empty → `download_data()` raises `RuntimeError`.

### Existing tests

All existing tests pass unchanged (no paperconfig has `doi:` yet).

## Data Fix

Add `doi: "10.1186/2046-9063-8-7"` to Fuszard 2012 paperconfig. After KG rebuild + Docker redeploy, `test_publication_doi_or_pmid_present` passes.

## Risks

1. **DOIs contain `/`.** Experiment IDs become e.g. `10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq`. Neo4j handles `/` in string properties fine. Spot-check that MCP tools and Cypher fixtures don't parse IDs on `/`.

2. **Hard error on PDF extraction failure is a behavior change.** Verified: all 36 active paperconfigs with `publication:` blocks have valid PDF paths. Resource-only configs (like MIT9313_resources) are unaffected.

## Out of Scope

- `pmid:` override field.
- Title/authors/journal/abstract overrides from config.
- Changes to `pdf_publication_extraction.py` or `pdf_extraction_cache.json`.
- Migration of existing graph data — users rebuild KG.
