# Paperconfig DOI Override — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow paperconfig authors to set an optional `doi:` field under `publication:` that overrides the PDF-extracted DOI for Publication and Experiment node IDs.

**Architecture:** A `_get_override_doi()` helper on `OMICSAdapter` reads the config field and feeds it into `get_publication_id()` (with a warning on disagreement) and `get_publication_nodes()` (for the `doi` property). The validator gets a DOI regex check. `download_data()` becomes strict: hard-error when a publication block exists but PDF extraction fails.

**Tech Stack:** Python, pytest, pyyaml, regex

**Spec:** `docs/superpowers/specs/2026-04-16-paperconfig-doi-override-design.md`

---

### Task 1: Create feature branch

**Files:** None

- [ ] **Step 1: Create and switch to feature branch**

```bash
git checkout -b feat/paperconfig-doi-override
```

- [ ] **Step 2: Verify clean state**

```bash
git status
```

Expected: `On branch feat/paperconfig-doi-override`, nothing to commit.

---

### Task 2: Validator — DOI format validation

**Files:**
- Modify: `scripts/validate_paperconfig.py:28` (add `import re`), `:846-860` (add DOI check after papername)
- Test: `tests/test_paperconfig_validation.py`

- [ ] **Step 1: Write the failing tests**

Add these tests to `tests/test_paperconfig_validation.py`, inside a new class at the end of the file:

```python
class TestDoiOverride:
    """Validator accepts/rejects optional publication.doi field."""

    def test_valid_doi_passes(self, tmp_path, monkeypatch):
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "10.1186/2046-9063-8-7"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is True

    def test_missing_doi_passes(self, tmp_path, monkeypatch):
        """doi is optional — omitting it is fine."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        # No doi key at all
        assert "doi" not in config["publication"]
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is True

    def test_invalid_doi_rejected(self, tmp_path, monkeypatch):
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "not-a-doi"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is False

    def test_doi_url_rejected(self, tmp_path, monkeypatch):
        """Full URL form should be rejected — we want the bare DOI."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "https://doi.org/10.1186/2046-9063-8-7"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is False
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_paperconfig_validation.py::TestDoiOverride -v
```

Expected: `test_valid_doi_passes` passes (no DOI check exists = no error), `test_missing_doi_passes` passes, but `test_invalid_doi_rejected` and `test_doi_url_rejected` FAIL (validate returns True for bad DOIs since no check exists yet).

- [ ] **Step 3: Add DOI regex and validation to the validator**

In `scripts/validate_paperconfig.py`:

1. Add `import re` at the top (after `import json`, around line 31).

2. Add a module-level constant (after the existing constant blocks, e.g. after line 56):

```python
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
```

3. After the `papername` check (after line 851, before the `# --- PDF ---` comment at line 853), add:

```python
        # --- DOI override ---
        doi_override = pub.get("doi")
        if doi_override is not None:
            if not isinstance(doi_override, str) or not DOI_RE.match(doi_override.strip()):
                errors.append(f"publication.doi '{doi_override}' is not a valid DOI (expected pattern 10.NNNN/...)")
            else:
                print(f"  doi (override): {doi_override.strip()}")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_paperconfig_validation.py::TestDoiOverride -v
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Run the full validator test suite**

```bash
pytest tests/test_paperconfig_validation.py -v
```

Expected: all tests PASS (no regressions).

- [ ] **Step 6: Commit**

```bash
git add scripts/validate_paperconfig.py tests/test_paperconfig_validation.py
git commit -m "feat(validator): add optional doi field validation to paperconfig

Validates that publication.doi, when present, matches the DOI pattern
^10.\d{4,9}/\S+$. The field is optional — omitting it is fine."
```

---

### Task 3: Adapter — `_get_override_doi()` helper and `get_publication_id()` override

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py:842-851` (`get_publication_id`), add `_get_override_doi` nearby
- Test: `tests/test_omics_adapter.py`

- [ ] **Step 1: Write the failing tests**

Add these tests to the end of `tests/test_omics_adapter.py`:

```python
class TestDoiOverride:
    """Config doi overrides PDF-extracted doi for pub and experiment IDs."""

    def _make_config_with_doi(self, tmp_path, doi_override=None):
        """Build a minimal paperconfig with optional doi override."""
        de_csv = tmp_path / "de.csv"
        pd.DataFrame({
            "gene": ["PMM0001"],
            "log2fc": [1.5],
        }).to_csv(de_csv, index=False)

        paperconfig = {
            "publication": {
                "papername": "Test",
                "papermainpdf": str(tmp_path / "dummy.pdf"),
                "experiments": {
                    "exp1": {
                        "name": "Test experiment",
                        "organism": "Prochlorococcus MED4",
                        "treatment_condition": "test",
                        "control_condition": "control",
                        "omics_type": "RNASEQ",
                        "test_type": "DESeq2",
                        "treatment_type": ["nitrogen"],
                    },
                },
                "supplementary_materials": {
                    "tbl": {
                        "type": "csv",
                        "filename": str(de_csv),
                        "statistical_analyses": [{
                            "id": "DE_test",
                            "experiment": "exp1",
                            "name_col": "gene",
                            "logfc_col": "log2fc",
                        }],
                    },
                },
            },
        }
        if doi_override:
            paperconfig["publication"]["doi"] = doi_override

        config_file = tmp_path / "paperconfig.yaml"
        config_file.write_text(yaml.dump(paperconfig))
        return str(config_file)

    def test_config_doi_overrides_publication_id(self, tmp_path):
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            "publication": {
                "publication_id": "pub_SomethingElse",
                "doi": "10.0000/other",
                "title": "Some Title",
            },
        }

        pub_id = adapter.get_publication_id()
        assert pub_id == "10.9999/testdoi"

    def test_config_doi_overrides_experiment_ids(self, tmp_path):
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            "publication": {
                "publication_id": "pub_SomethingElse",
                "doi": "10.0000/other",
                "title": "Some Title",
            },
        }

        nodes = adapter.get_nodes()
        experiment_nodes = [n for n in nodes if n[1] == "experiment"]
        assert len(experiment_nodes) == 1
        exp_id = experiment_nodes[0][0]
        assert exp_id.startswith("10.9999/testdoi_"), f"Experiment id should start with override DOI, got {exp_id}"

    def test_config_doi_flows_to_publication_node_property(self, tmp_path):
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            "publication": {
                "publication_id": "pub_SomethingElse",
                "doi": "10.0000/other",
                "title": "Some Title",
            },
        }

        pub_nodes = adapter.get_publication_nodes()
        assert len(pub_nodes) == 1
        props = pub_nodes[0][2]
        assert props["doi"] == "10.9999/testdoi"

    def test_no_config_doi_uses_pdf_extracted(self, tmp_path):
        config_file = self._make_config_with_doi(tmp_path, doi_override=None)
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            "publication": {
                "publication_id": "10.1234/from.pdf",
                "doi": "10.1234/from.pdf",
                "title": "PDF Title",
            },
        }

        pub_id = adapter.get_publication_id()
        assert pub_id == "10.1234/from.pdf"

    def test_config_doi_warns_on_disagreement(self, tmp_path, caplog):
        import logging
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            "publication": {
                "doi": "10.0000/different",
                "title": "Title",
            },
        }

        with caplog.at_level(logging.WARNING):
            pub_id = adapter.get_publication_id()

        assert pub_id == "10.9999/testdoi"
        assert "disagrees" in caplog.text.lower() or "Config doi" in caplog.text
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_omics_adapter.py::TestDoiOverride -v
```

Expected: all 5 tests FAIL (no `_get_override_doi` method yet; `get_publication_id` ignores config doi).

- [ ] **Step 3: Implement `_get_override_doi()` and update `get_publication_id()`**

In `multiomics_kg/adapters/omics_adapter.py`, add the helper just before `get_publication_id()` (before line 842):

```python
    def _get_override_doi(self) -> str | None:
        """Read optional doi from paperconfig publication block."""
        pub = self.config_data.get("publication", {}) or {}
        val = pub.get("doi")
        if not isinstance(val, str):
            return None
        val = val.strip()
        return val or None
```

Replace `get_publication_id()` (lines 842-851) with:

```python
    def get_publication_id(self) -> str:
        """Get the publication ID from config doi override, PDF extraction, or fallback."""
        override = self._get_override_doi()
        if hasattr(self, 'extracted_data') and "publication" in self.extracted_data:
            pub = self.extracted_data["publication"]
            extracted_doi = pub.get("doi")
            if override:
                if extracted_doi and extracted_doi != override:
                    logger.warning(
                        f"Config doi '{override}' disagrees with PDF-extracted doi "
                        f"'{extracted_doi}'; using config value."
                    )
                return override
            pub_id = pub.get("publication_id") or extracted_doi or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:20]}"
            return str(pub_id)
        if override:
            return override
        pub = self.config_data.get('publication', {})
        pub_id = pub.get("pubmed_id") or pub.get("papername") or "unknown"
        return str(pub_id)
```

- [ ] **Step 4: Update `get_publication_nodes()` doi property**

In `get_publication_nodes()`, change line 285 from:

```python
                "doi": pub.get("doi"),
```

to:

```python
                "doi": self._get_override_doi() or pub.get("doi"),
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/test_omics_adapter.py::TestDoiOverride -v
```

Expected: all 5 tests PASS.

- [ ] **Step 6: Run full omics adapter tests**

```bash
pytest tests/test_omics_adapter.py tests/test_omics_adapter_organism_gene.py tests/test_omics_adapter_fold_change_type.py -v
```

Expected: all existing tests PASS (no regressions).

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py tests/test_omics_adapter.py
git commit -m "feat(adapter): add doi override from paperconfig

Config publication.doi takes precedence over PDF-extracted DOI for
Publication node ID, doi property, and Experiment IDs. Warns when
config and PDF DOIs disagree."
```

---

### Task 4: Adapter — hard-error on PDF extraction failure

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py:240-262` (`download_data`)
- Test: `tests/test_omics_adapter.py`

- [ ] **Step 1: Write the failing tests**

Add to the `TestDoiOverride` class in `tests/test_omics_adapter.py`:

```python
    def test_download_data_raises_on_missing_pdf(self, tmp_path):
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        # papermainpdf points to tmp_path/dummy.pdf which does not exist
        with pytest.raises(FileNotFoundError):
            adapter.download_data()

    def test_download_data_raises_on_empty_extraction(self, tmp_path, monkeypatch):
        # Create the dummy PDF so the path check passes
        dummy_pdf = tmp_path / "dummy.pdf"
        dummy_pdf.write_text("fake pdf")
        config_file = self._make_config_with_doi(tmp_path, doi_override="10.9999/testdoi")
        adapter = OMICSAdapter(config_file=config_file)
        # Mock extractor to return empty dict
        monkeypatch.setattr(adapter.pdf_extractor, "extract_from_pdf", lambda path: {})
        with pytest.raises(RuntimeError, match="no publication block"):
            adapter.download_data()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_omics_adapter.py::TestDoiOverride::test_download_data_raises_on_missing_pdf tests/test_omics_adapter.py::TestDoiOverride::test_download_data_raises_on_empty_extraction -v
```

Expected: both FAIL (currently logs a warning instead of raising).

- [ ] **Step 3: Update `download_data()` to raise on failure**

Replace lines 257-262 of `multiomics_kg/adapters/omics_adapter.py` (from `# extract publication metadata from pdf` through the `else: logger.warning(...)` block) with:

```python
        # extract publication metadata from pdf
        pdf_path = publication.get('papermainpdf', None)
        if not pdf_path or not os.path.exists(pdf_path):
            raise FileNotFoundError(
                f"PDF path missing or not found for "
                f"{publication.get('papername', '?')}: {pdf_path!r}"
            )
        self.extracted_data = self.pdf_extractor.extract_from_pdf(pdf_path)
        if not self.extracted_data or "publication" not in self.extracted_data:
            raise RuntimeError(
                f"PDF extraction produced no publication block for {pdf_path}"
            )
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_omics_adapter.py::TestDoiOverride -v
```

Expected: all 7 tests PASS.

- [ ] **Step 5: Run full test suite (non-slow, non-kg)**

```bash
pytest -m "not slow and not kg" -v --tb=short 2>&1 | tail -20
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py tests/test_omics_adapter.py
git commit -m "fix(adapter): hard-error when PDF extraction fails

Replace silent logger.warning + skip with FileNotFoundError (missing
PDF) or RuntimeError (empty extraction). Resource-only configs without
a publication block are unaffected."
```

---

### Task 5: Fuszard 2012 data fix

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Fuszard 2012/paperconfig.yaml:1-3`

- [ ] **Step 1: Add doi to Fuszard 2012 paperconfig**

In `data/Prochlorococcus/papers_and_supp/Fuszard 2012/paperconfig.yaml`, change line 2 from:

```yaml
  papername: "Fuszard 2012"
```

to:

```yaml
  papername: "Fuszard 2012"
  doi: "10.1186/2046-9063-8-7"
```

(Insert `doi:` as a new line 3, between `papername:` and the comment block.)

- [ ] **Step 2: Run the validator on Fuszard 2012**

```bash
uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Fuszard 2012/paperconfig.yaml"
```

Expected: passes, prints `doi (override): 10.1186/2046-9063-8-7`.

- [ ] **Step 3: Run the full parametrized validator tests**

```bash
pytest tests/test_paperconfig_validation.py::test_paperconfig_validates -v --tb=short 2>&1 | tail -40
```

Expected: all paperconfigs pass (including Fuszard 2012).

- [ ] **Step 4: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Fuszard 2012/paperconfig.yaml"
git commit -m "data: add doi override to Fuszard 2012 paperconfig

Fixes Publication node ID from 'pub_Comparative quantitative prote'
to '10.1186/2046-9063-8-7'. Resolves kg-validity test failure in
test_publication_doi_or_pmid_present."
```

---

### Task 6: Spot-check DOI-with-slash in Experiment IDs

**Files:** None (read-only investigation)

- [ ] **Step 1: Grep for ID-parsing on `/` in MCP and Cypher**

```bash
grep -rn 'split.*/' multiomics_kg/mcp/ scripts/post-import.sh scripts/post-import.cypher tests/kg_validity/ 2>/dev/null | grep -v '.pyc' | grep -v __pycache__ || echo "No slash-splitting found"
```

Expected: no hits that parse node IDs on `/`. DOIs with `/` are safe.

- [ ] **Step 2: Document result**

If clean, no action needed. If any parsing on `/` is found, document it and add a follow-up fix to this branch before merging.

---

### Task 7: Update CLAUDE.md and skill docs

**Files:**
- Modify: `CLAUDE.md` (paperconfig example block)
- Modify: `.claude/skills/paperconfig/SKILL.md` (if it has a template — add `doi:` as optional field)

- [ ] **Step 1: Update CLAUDE.md paperconfig example**

In the "Adding Omics Data from Publications" section of `CLAUDE.md`, update the example `publication:` block to show the optional `doi:` field:

```yaml
publication:
  papername: "Author Year"
  doi: "10.NNNN/xxxxx"  # optional — overrides PDF-extracted DOI
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Author Year/paper.pdf"
```

- [ ] **Step 2: Update SKILL.md if it has a template**

Check `.claude/skills/paperconfig/SKILL.md` for a template block. If it shows a `publication:` example, add `doi:` as an optional field with a comment.

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md .claude/skills/paperconfig/SKILL.md
git commit -m "docs: document optional doi field in paperconfig publication block"
```

---

### Task 8: Final verification

**Files:** None

- [ ] **Step 1: Run the full non-slow test suite**

```bash
pytest -m "not slow and not kg" -v --tb=short 2>&1 | tail -30
```

Expected: all tests PASS.

- [ ] **Step 2: Review the branch diff**

```bash
git log --oneline main..HEAD
git diff main --stat
```

Expected: 4-5 commits, changes in `omics_adapter.py`, `validate_paperconfig.py`, `test_omics_adapter.py`, `test_paperconfig_validation.py`, Fuszard paperconfig, `CLAUDE.md`, and optionally `SKILL.md`.
