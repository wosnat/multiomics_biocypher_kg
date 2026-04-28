# Phase 1.1A — Metabolite foundation: download + skeletons

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land sub-step 6 (raw downloads of MNX/TCDB/CAZy reference data) and stub-out the rest of Phase 1.1, so we can audit the actual downloaded files before designing the resolver/accessor implementations.

**Architecture:** sub-step 6 is real (HTTP downloads + caching). Sub-step 7 (resolver builder), the three `utils/` accessor modules, and the gene-annotation transforms are stubbed (importable, raise `NotImplementedError`). After this plan ships, we run the download once, inspect actual file shapes, and update Spec 1.1 inline before writing the next plan (Phase 1.1B).

**Tech Stack:** Python 3.10+, `requests` for HTTP, pytest with `monkeypatch` for HTTP mocking, the existing `prepare_data.sh` / `download_genome_data.py` orchestrator.

**Spec:** [`docs/superpowers/specs/2026-04-28-metabolite-foundation-design.md`](../specs/2026-04-28-metabolite-foundation-design.md)

---

## File structure

**Create (new files):**

| Path | Role |
|---|---|
| `multiomics_kg/utils/metabolite_utils.py` | Skeleton: resolver accessor API. All functions raise `NotImplementedError` until 1.1B. |
| `multiomics_kg/utils/tcdb_utils.py` | Skeleton: TCDB hierarchy accessor API. |
| `multiomics_kg/utils/cazy_utils.py` | Skeleton: CAZy hierarchy accessor API. |
| `multiomics_kg/download/build_metabolite_resolver.py` | Skeleton: sub-step 7 module. `main()` raises `NotImplementedError`. |
| `multiomics_kg/download/download_metabolism_reference.py` | **Real**: sub-step 6 implementation. HTTP downloads → `cache/data/{mnx,tcdb,cazy}/`. |
| `tests/test_download_metabolism_reference.py` | **Real**: mocked HTTP + file-write tests for sub-step 6. |
| `tests/test_metabolism_skeletons.py` | **Real**: smoke tests that the four skeleton modules import and that `NotImplementedError` is raised by their stubs. |

**Modify (existing files):**

| Path | Lines | Why |
|---|---|---|
| `multiomics_kg/download/download_genome_data.py` | argparse + dispatch (~lines 460–509) | Extend `--steps choices` to `[1..7]`; add `step6_metabolism_reference` + `step7_metabolite_resolver` functions and dispatch. |
| `scripts/prepare_data.sh` | DOWNLOAD_SUBSTEPS lines (~127–131) | Append `6 7` to both DOWNLOAD_SUBSTEPS variants. |
| `CLAUDE.md` | "Genome Data Download Pipeline" section | Document sub-steps 6 + 7 and the granular re-run command. |
| `docs/superpowers/specs/2026-04-28-metabolite-foundation-design.md` | inline | Post-audit: replace tentative schemas with confirmed ones. (Task 11 only.) |

---

## Tasks

### Task 1: Skeleton — `utils/metabolite_utils.py`

**Files:**
- Create: `multiomics_kg/utils/metabolite_utils.py`
- Test: `tests/test_metabolism_skeletons.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_metabolism_skeletons.py
"""Smoke tests for Phase 1.1A skeleton modules.

These verify the modules import and expose the expected API surface.
Real behaviour ships in Phase 1.1B; until then every stub raises
NotImplementedError.
"""
import sqlite3
import pytest


def test_metabolite_utils_imports_and_raises():
    from multiomics_kg.utils import metabolite_utils

    # API surface is present
    assert hasattr(metabolite_utils, "open_resolver")
    assert hasattr(metabolite_utils, "resolve_metabolite")
    assert hasattr(metabolite_utils, "resolve_reaction")

    # Stubs raise NotImplementedError
    with pytest.raises(NotImplementedError):
        metabolite_utils.open_resolver()
    fake_conn = sqlite3.connect(":memory:")
    with pytest.raises(NotImplementedError):
        metabolite_utils.resolve_metabolite("C00031", fake_conn)
    with pytest.raises(NotImplementedError):
        metabolite_utils.resolve_reaction("R00200", fake_conn)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_metabolism_skeletons.py::test_metabolite_utils_imports_and_raises -v`
Expected: FAIL — `ModuleNotFoundError: multiomics_kg.utils.metabolite_utils`.

- [ ] **Step 3: Write the skeleton**

```python
# multiomics_kg/utils/metabolite_utils.py
"""Resolver accessor API.

Phase 1.1A skeleton: every function raises NotImplementedError. Real
implementations land in Phase 1.1B once the actual MNX file shapes are
confirmed via the audit.

API contract (consumed by download/build_gene_annotations.py via the
transforms framework, by the Spec 1.2 scaffold builder, and by Phase 2
paper-measurement extraction):

    open_resolver(path: Path | None = None) -> sqlite3.Connection
    resolve_metabolite(value: str, conn) -> tuple[str | None, str]
    resolve_reaction (value: str, conn) -> tuple[str | None, str]
"""
from __future__ import annotations

import sqlite3
from pathlib import Path


def open_resolver(path: Path | None = None) -> sqlite3.Connection:
    raise NotImplementedError("Phase 1.1B")


def resolve_metabolite(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    raise NotImplementedError("Phase 1.1B")


def resolve_reaction(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    raise NotImplementedError("Phase 1.1B")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_metabolism_skeletons.py::test_metabolite_utils_imports_and_raises -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/metabolite_utils.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: skeleton metabolite_utils accessor (Phase 1.1A)"
```

---

### Task 2: Skeleton — `utils/tcdb_utils.py`

**Files:**
- Create: `multiomics_kg/utils/tcdb_utils.py`
- Test: `tests/test_metabolism_skeletons.py` (append to existing file)

- [ ] **Step 1: Append the failing test**

```python
# tests/test_metabolism_skeletons.py — append below the metabolite_utils test
def test_tcdb_utils_imports_and_raises():
    from multiomics_kg.utils import tcdb_utils

    assert hasattr(tcdb_utils, "load_tcdb")
    assert hasattr(tcdb_utils, "tcdb_ancestors")
    assert hasattr(tcdb_utils, "is_valid_tcdb")

    with pytest.raises(NotImplementedError):
        tcdb_utils.load_tcdb()
    with pytest.raises(NotImplementedError):
        tcdb_utils.tcdb_ancestors("1.A.1.1.1")
    with pytest.raises(NotImplementedError):
        tcdb_utils.is_valid_tcdb("1.A.1.1.1")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_metabolism_skeletons.py::test_tcdb_utils_imports_and_raises -v`
Expected: FAIL — `ModuleNotFoundError: multiomics_kg.utils.tcdb_utils`.

- [ ] **Step 3: Write the skeleton**

```python
# multiomics_kg/utils/tcdb_utils.py
"""TCDB hierarchy accessor API.

Phase 1.1A skeleton — every function raises NotImplementedError until 1.1B.

API:
    load_tcdb()                 -> dict[str, dict]      # cached at module level
    tcdb_ancestors(tc_id: str)  -> list[str]            # ["1", "1.A", "1.A.1", "1.A.1.1"] for "1.A.1.1.1"
    is_valid_tcdb(tc_id: str)   -> bool
"""
from __future__ import annotations


def load_tcdb() -> dict[str, dict]:
    raise NotImplementedError("Phase 1.1B")


def tcdb_ancestors(tc_id: str) -> list[str]:
    raise NotImplementedError("Phase 1.1B")


def is_valid_tcdb(tc_id: str) -> bool:
    raise NotImplementedError("Phase 1.1B")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_metabolism_skeletons.py::test_tcdb_utils_imports_and_raises -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_utils.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: skeleton tcdb_utils accessor (Phase 1.1A)"
```

---

### Task 3: Skeleton — `utils/cazy_utils.py`

**Files:**
- Create: `multiomics_kg/utils/cazy_utils.py`
- Test: `tests/test_metabolism_skeletons.py` (append)

- [ ] **Step 1: Append the failing test**

```python
def test_cazy_utils_imports_and_raises():
    from multiomics_kg.utils import cazy_utils

    assert hasattr(cazy_utils, "load_cazy")
    assert hasattr(cazy_utils, "cazy_ancestors")
    assert hasattr(cazy_utils, "is_valid_cazy")

    with pytest.raises(NotImplementedError):
        cazy_utils.load_cazy()
    with pytest.raises(NotImplementedError):
        cazy_utils.cazy_ancestors("GH13_1")
    with pytest.raises(NotImplementedError):
        cazy_utils.is_valid_cazy("GH13_1")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_metabolism_skeletons.py::test_cazy_utils_imports_and_raises -v`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 3: Write the skeleton**

```python
# multiomics_kg/utils/cazy_utils.py
"""CAZy hierarchy accessor API.

Phase 1.1A skeleton — every function raises NotImplementedError until 1.1B.

API:
    load_cazy()                  -> dict[str, dict]     # cached at module level
    cazy_ancestors(cazy_id: str) -> list[str]
    is_valid_cazy(cazy_id: str)  -> bool
"""
from __future__ import annotations


def load_cazy() -> dict[str, dict]:
    raise NotImplementedError("Phase 1.1B")


def cazy_ancestors(cazy_id: str) -> list[str]:
    raise NotImplementedError("Phase 1.1B")


def is_valid_cazy(cazy_id: str) -> bool:
    raise NotImplementedError("Phase 1.1B")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_metabolism_skeletons.py::test_cazy_utils_imports_and_raises -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/cazy_utils.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: skeleton cazy_utils accessor (Phase 1.1A)"
```

---

### Task 4: Skeleton — `download/build_metabolite_resolver.py`

**Files:**
- Create: `multiomics_kg/download/build_metabolite_resolver.py`
- Test: `tests/test_metabolism_skeletons.py` (append)

- [ ] **Step 1: Append the failing test**

```python
def test_build_metabolite_resolver_imports_and_raises():
    from multiomics_kg.download import build_metabolite_resolver

    assert hasattr(build_metabolite_resolver, "main")

    with pytest.raises(NotImplementedError):
        build_metabolite_resolver.main(force=False)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_metabolism_skeletons.py::test_build_metabolite_resolver_imports_and_raises -v`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 3: Write the skeleton**

```python
# multiomics_kg/download/build_metabolite_resolver.py
"""Step 0 sub-step 7 — Build resolver + hierarchy caches.

Phase 1.1A skeleton — main() raises NotImplementedError. Real implementation
lands in Phase 1.1B once the actual MNX/TCDB/CAZy file shapes are confirmed.

Outputs (when implemented):
- cache/data/mnx/metabolite_resolver.db        (SQLite)
- cache/data/tcdb/tcdb_hierarchy.json
- cache/data/cazy/cazy_hierarchy.json
- cache/data/mnx/metabolite_id_mapping_report.json
"""
from __future__ import annotations

import argparse


def main(force: bool = False) -> None:
    raise NotImplementedError("Phase 1.1B")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild caches even if they exist.")
    args = parser.parse_args()
    main(force=args.force)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_metabolism_skeletons.py -v`
Expected: all four tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: skeleton build_metabolite_resolver (Phase 1.1A)"
```

---

### Task 5: Implement `download/download_metabolism_reference.py`

**Files:**
- Create: `multiomics_kg/download/download_metabolism_reference.py`
- Test: `tests/test_download_metabolism_reference.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_download_metabolism_reference.py
"""Tests for sub-step 6 download module.

HTTP is mocked via monkeypatch; tests assert that the right URLs are
requested and that bytes are written to the right cache paths.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from multiomics_kg.download import download_metabolism_reference as dmr


class _FakeResponse:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code
        self._chunks = [content]

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int = 8192):
        return iter(self._chunks)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _patch_requests_get(monkeypatch, response_map: dict[str, _FakeResponse]):
    """Patch requests.get to look up responses by URL."""
    def fake_get(url, *args, **kwargs):
        if url not in response_map:
            raise AssertionError(f"unexpected URL: {url}")
        return response_map[url]
    monkeypatch.setattr(dmr.requests, "get", fake_get)


def test_sources_table_has_six_entries():
    """Six expected sources: 4 MNX TSVs + 1 TCDB + 1 CAZy."""
    assert len(dmr.SOURCES) == 6
    assert {"mnx_chem_prop", "mnx_chem_xref",
            "mnx_reac_prop", "mnx_reac_xref",
            "tcdb_families", "cazy_families"} == set(dmr.SOURCES)


def test_download_writes_all_files(monkeypatch, tmp_path):
    """One download per source; each writes to the configured relative path."""
    fake_responses = {
        url: _FakeResponse(f"BODY-{key}".encode())
        for key, (url, _rel) in dmr.SOURCES.items()
    }
    _patch_requests_get(monkeypatch, fake_responses)

    dmr.download_all(cache_root=tmp_path, force=True)

    for key, (_url, rel) in dmr.SOURCES.items():
        out = tmp_path / rel
        assert out.exists(), f"missing output for {key}: {out}"
        assert out.read_bytes() == f"BODY-{key}".encode()


def test_skip_when_cached_unless_force(monkeypatch, tmp_path):
    """Existing files are skipped when force=False; re-downloaded when force=True."""
    # Pre-create the MNX chem_prop cache file with stale content
    target = tmp_path / dmr.SOURCES["mnx_chem_prop"][1]
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b"STALE")

    # Mock all URLs to return fresh content
    fake_responses = {
        url: _FakeResponse(b"FRESH")
        for url, _rel in dmr.SOURCES.values()
    }
    _patch_requests_get(monkeypatch, fake_responses)

    # Without --force: stale content preserved
    dmr.download_all(cache_root=tmp_path, force=False)
    assert target.read_bytes() == b"STALE"

    # With --force: stale content overwritten
    dmr.download_all(cache_root=tmp_path, force=True)
    assert target.read_bytes() == b"FRESH"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_download_metabolism_reference.py -v`
Expected: FAIL — `ModuleNotFoundError: multiomics_kg.download.download_metabolism_reference`.

- [ ] **Step 3: Implement the download module**

```python
# multiomics_kg/download/download_metabolism_reference.py
"""Step 0 sub-step 6 — Download MNX, TCDB, CAZy reference data.

Six small downloads (~50–100 MB total, dominated by MNX). Files are cached
under cache/data/{mnx,tcdb,cazy}/ and skipped on re-run unless --force.

URLs are best-effort guesses for Phase 1.1A and may need adjustment after the
first real run; see the audit task in plan 1.1A.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import requests

log = logging.getLogger(__name__)

# (URL, cache-relative path)
SOURCES: dict[str, tuple[str, str]] = {
    "mnx_chem_prop":  ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",  "mnx/chem_prop.tsv"),
    "mnx_chem_xref":  ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",  "mnx/chem_xref.tsv"),
    "mnx_reac_prop":  ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",  "mnx/reac_prop.tsv"),
    "mnx_reac_xref":  ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",  "mnx/reac_xref.tsv"),
    "tcdb_families":  ("https://www.tcdb.org/cgi-bin/projectv/public/families.py",      "tcdb/families.tsv"),
    "cazy_families":  ("http://www.cazy.org/IMG/cazy_data/family_classification.txt",   "cazy/families.tsv"),
}

DEFAULT_CACHE_ROOT = Path("cache/data")


def download_one(url: str, dest: Path, force: bool = False) -> bool:
    """Download URL → dest. Returns True if downloaded, False if skipped."""
    if dest.exists() and not force:
        log.info(f"  skip (cached): {dest}")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"  GET {url} → {dest}")
    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=64 * 1024):
                if chunk:
                    f.write(chunk)
    return True


def download_all(cache_root: Path = DEFAULT_CACHE_ROOT, force: bool = False) -> None:
    """Download every source in SOURCES into cache_root."""
    log.info(f"download_metabolism_reference: cache_root={cache_root} force={force}")
    n_downloaded = 0
    for key, (url, rel) in SOURCES.items():
        if download_one(url, cache_root / rel, force=force):
            n_downloaded += 1
    log.info(f"  done — {n_downloaded}/{len(SOURCES)} downloaded, "
             f"{len(SOURCES) - n_downloaded} cached.")


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    download_all(force=force)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if cache files exist.")
    args = parser.parse_args()
    main(force=args.force)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_download_metabolism_reference.py -v`
Expected: all three tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/download_metabolism_reference.py tests/test_download_metabolism_reference.py
git commit -m "metabolism: implement sub-step 6 download module (Phase 1.1A)"
```

---

### Task 6: Wire sub-steps 6 + 7 into `download_genome_data.py`

**Files:**
- Modify: `multiomics_kg/download/download_genome_data.py` (argparse choices + dispatch + two new step functions)

- [ ] **Step 1: Read the existing dispatch block**

Read [`multiomics_kg/download/download_genome_data.py:441-511`](../../multiomics_kg/download/download_genome_data.py#L441-L511) so you know exactly what to replace. Note the existing pattern: `step1_ncbi`, …, `step5_gene_mapping` functions defined at module level, dispatched by `if N in steps:` lines in `main()`.

- [ ] **Step 2: Add the two new step functions above `def main()`**

Insert immediately before `def main() -> None:` (~line 441):

```python
def step6_metabolism_reference(force: bool) -> None:
    """Sub-step 6: download MNX/TCDB/CAZy reference data."""
    from multiomics_kg.download.download_metabolism_reference import download_all
    log.info("─── Step 6: Download MNX/TCDB/CAZy reference data ───")
    download_all(force=force)


def step7_metabolite_resolver(force: bool) -> None:
    """Sub-step 7: build metabolite resolver + hierarchy caches."""
    from multiomics_kg.download.build_metabolite_resolver import main as build_main
    log.info("─── Step 7: Build metabolite resolver + hierarchy caches ───")
    build_main(force=force)
```

- [ ] **Step 3: Extend argparse `--steps` choices and help text**

In `def main()`, replace the `parser.add_argument("--steps", ...)` block with:

```python
    parser.add_argument(
        "--steps", nargs="+", type=int, choices=[1, 2, 3, 4, 5, 6, 7],
        default=[1, 2, 3, 4, 5, 6, 7],
        help=("Steps to run (default: all). 1=NCBI 2=Cyanorak 3=UniProt "
              "4=eggNOG 5=gene_mapping 6=metabolism_reference 7=metabolite_resolver"),
    )
```

Also extend the `epilog` string (~line 446) to list 6 + 7:

```python
        epilog="""
Steps:
  1  NCBI genome (GFF + protein FASTA + GBFF)
  2  Cyanorak GFF/GBK (strains with cyanorak_organism only)
  3  UniProt (one download per unique taxid)
  4  eggNOG-mapper (requires EGGNOG_DATA_DIR in .env)
  5  gene_mapping.csv (requires steps 1+2)
  6  Download MNX/TCDB/CAZy reference data
  7  Build metabolite resolver + hierarchy caches (requires step 6)

Examples:
  uv run python multiomics_kg/download/download_genome_data.py
  uv run python multiomics_kg/download/download_genome_data.py --steps 1 2 3
  uv run python multiomics_kg/download/download_genome_data.py --steps 6 7 --force
  uv run python multiomics_kg/download/download_genome_data.py --strains MED4 --force
        """,
```

- [ ] **Step 4: Add the dispatch lines**

After the `if 5 in steps: step5_gene_mapping(...)` line in `main()`, append:

```python
    if 6 in steps:
        step6_metabolism_reference(force=args.force)
    if 7 in steps:
        step7_metabolite_resolver(force=args.force)
```

- [ ] **Step 5: Smoke-test the wiring**

Run: `uv run python multiomics_kg/download/download_genome_data.py --steps 7 --force 2>&1 | tail -5`
Expected: `NotImplementedError: Phase 1.1B` (sub-step 7 stub raises). This proves dispatch reaches the new function.

Run: `uv run python multiomics_kg/download/download_genome_data.py --help`
Expected: help output lists steps 6 and 7 with their descriptions.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/download_genome_data.py
git commit -m "metabolism: wire sub-steps 6 + 7 into download_genome_data dispatch"
```

---

### Task 7: Wire sub-steps 6 + 7 into `prepare_data.sh`

**Files:**
- Modify: `scripts/prepare_data.sh:127-131`

- [ ] **Step 1: Read the existing block**

```bash
grep -n DOWNLOAD_SUBSTEPS scripts/prepare_data.sh
```

You should see two assignments controlled by `$SKIP_CYANORAK`.

- [ ] **Step 2: Update both DOWNLOAD_SUBSTEPS assignments**

Replace the two-arm `if` with the extended sub-step lists (append `6 7` to both):

```bash
            if [[ $SKIP_CYANORAK -eq 1 ]]; then
                DOWNLOAD_SUBSTEPS="1 3 5 6 7"
                STEP0_LABEL="Download genome data (NCBI + UniProt + gene_mapping + metabolism reference + resolver; Cyanorak SKIPPED)"
            else
                DOWNLOAD_SUBSTEPS="1 2 3 5 6 7"
                STEP0_LABEL="Download genome data (NCBI + Cyanorak + UniProt + gene_mapping + metabolism reference + resolver)"
            fi
```

- [ ] **Step 3: Smoke-test by running step 0 with `--steps 0` (dry pass)**

Verify that the `STEPS` parser still accepts `--steps 0`. Look at the script: `--steps 0` runs step 0 which now calls `download_genome_data.py --steps 1 2 3 5 6 7`. Sub-step 7 is still a stub and will raise `NotImplementedError` — that's expected and isolates failure to one place.

```bash
bash scripts/prepare_data.sh --steps 0 --strains MED4 2>&1 | tail -10
```

Expected: most sub-steps run, sub-step 7 raises `NotImplementedError: Phase 1.1B`. The script will exit non-zero. This is the expected state at end of Phase 1.1A.

- [ ] **Step 4: Commit**

```bash
git add scripts/prepare_data.sh
git commit -m "metabolism: extend prepare_data.sh DOWNLOAD_SUBSTEPS to include 6 + 7"
```

---

### Task 8: Update `CLAUDE.md`

**Files:**
- Modify: `CLAUDE.md` ("Genome Data Download Pipeline" section)

- [ ] **Step 1: Find the Step 0 sub-steps section**

```bash
grep -n "Step 0.*download_genome_data" CLAUDE.md
```

You should land on a line like `**Step 0** (\`multiomics_kg/download/download_genome_data.py\`) — sub-steps:` followed by 5 bullets.

- [ ] **Step 2: Append two bullets for sub-steps 6 + 7**

After the existing `- 5: Build gene_mapping.csv` line, insert:

```markdown
- 6: Download MNX/TCDB/CAZy reference data → `cache/data/{mnx,tcdb,cazy}/`
- 7: Build metabolite resolver + hierarchy caches (requires sub-step 6). Phase 1.1A: skeleton only — full build lands in Phase 1.1B.
```

- [ ] **Step 3: Add a "Granular metabolism re-run" note**

Locate the `bash scripts/prepare_data.sh ...` examples block (a few `# rebuild ... only` comments). Append this example:

```bash
# Refresh metabolism caches only (bypasses prepare_data.sh; avoids re-running NCBI/UniProt):
uv run python multiomics_kg/download/download_genome_data.py --steps 6 7 --force
# Or just rebuild the resolver from already-downloaded raw files:
uv run python multiomics_kg/download/download_genome_data.py --steps 7 --force
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document Phase 1.1A sub-steps 6 + 7 in CLAUDE.md"
```

---

### Task 9: Run the real download (manual)

**Files:** none modified — produces files under `cache/data/{mnx,tcdb,cazy}/`.

This is **not a test step** — we're populating the cache once with real bytes so we can audit them in Task 10.

- [ ] **Step 1: Run the download**

```bash
uv run python multiomics_kg/download/download_genome_data.py --steps 6 --force
```

Expected: ~30 seconds to a few minutes depending on bandwidth. Six files downloaded.

- [ ] **Step 2: Verify outputs landed**

```bash
ls -lh cache/data/mnx/ cache/data/tcdb/ cache/data/cazy/
```

Expected: six files, sizes in the rough ballpark below. Record actual sizes — they go into the audit notes.

| Expected file | Rough size |
|---|---|
| `cache/data/mnx/chem_prop.tsv` | 30–60 MB |
| `cache/data/mnx/chem_xref.tsv` | 30–80 MB |
| `cache/data/mnx/reac_prop.tsv` | 5–20 MB |
| `cache/data/mnx/reac_xref.tsv` | 5–20 MB |
| `cache/data/tcdb/families.tsv` | < 1 MB |
| `cache/data/cazy/families.tsv` | < 1 MB |

- [ ] **Step 3: If any URL 404s or returns HTML instead of data**

This is expected for TCDB and CAZy — the URLs are guesses. Capture the failure, note which URL failed and what the server returned, then in Task 11 update `SOURCES` in `download_metabolism_reference.py` with the corrected URL and re-run. **Do not commit yet** — finish Tasks 10 + 11 first so the URL fix lands together with the audit-driven schema updates.

- [ ] **Step 4: No commit yet — proceed to Task 10**

The cache files themselves are not committed by this step (gitignore decision is part of Task 11). The Python code already downloaded by previous tasks is committed; this task is a runtime invocation.

---

### Task 10: Schema audit (manual)

**Files:** none modified — produces audit notes that feed into Task 11.

The point of this task: open every downloaded file with `head` / `wc -l` / a small Python REPL session, write down what each looks like, and identify every place the spec's tentative schemas are wrong. The audit notes are scratch — they don't get committed; they're the input to Task 11's spec edits.

- [ ] **Step 1: Audit MNX `chem_prop.tsv`**

```bash
head -20 cache/data/mnx/chem_prop.tsv
wc -l cache/data/mnx/chem_prop.tsv
```

Record:
- Comment-line marker (`#` or `##`?)
- Header line columns (in order)
- Whether each row maps cleanly to the spec's `compounds (mnxm_id, name, formula, inchikey, smiles, charge, mass)` or has different/extra/missing columns
- Total row count

- [ ] **Step 2: Audit MNX `chem_xref.tsv`**

```bash
head -20 cache/data/mnx/chem_xref.tsv
awk -F'\t' 'NR>5 {print $1}' cache/data/mnx/chem_xref.tsv | head -30
awk -F'\t' 'NR>5 {split($1,a,":"); print a[1]}' cache/data/mnx/chem_xref.tsv | sort -u | head -40
```

Record:
- Column structure (combined `source:value` in one column, or separate `source` + `value` columns?)
- Distinct alias-source vocabulary (`kegg.compound`, `chebi`, `bigg`, …) — full list
- Whether MNX uses bioregistry-style prefixes or raw IDs

- [ ] **Step 3: Audit MNX `reac_prop.tsv` + `reac_xref.tsv`**

```bash
head -20 cache/data/mnx/reac_prop.tsv
head -20 cache/data/mnx/reac_xref.tsv
```

Record:
- `reac_prop.tsv` columns (does it carry `direction` info? equation text format?)
- Reaction xref source vocabulary
- Whether the spec's `direction_source ∈ {reversible, directional, unknown}` matches what MNX actually serves

- [ ] **Step 4: Audit TCDB `families.tsv`**

```bash
head -30 cache/data/tcdb/families.tsv
file cache/data/tcdb/families.tsv
wc -l cache/data/tcdb/families.tsv
```

Record:
- Real format (TSV? JSON? HTML? something else?)
- Hierarchy depth (does the file contain all 5 levels — class/subclass/family/subfamily/specificity — or just classes? families?)
- Whether `substrate_classes` info exists in this file or only in the full TCDB DB
- If format is wildly different from TSV, note the actual download URL needed

- [ ] **Step 5: Audit CAZy `families.tsv`**

```bash
head -30 cache/data/cazy/families.tsv
wc -l cache/data/cazy/families.tsv
```

Record:
- Real format
- Whether 3 levels (class/family/subfamily) are present
- Whether class info is implicit (from family prefix `GH13_1` → class `GH`) or carried as a separate column

- [ ] **Step 6: Audit eggNOG output for one strain**

The eggNOG output is **already cached** from prior runs — no new download needed.

```bash
head -5 cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations | tr '\t' '\n' | nl
awk -F'\t' '!/^#/ {print NF; exit}' cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations
```

Record:
- Actual column index of `KEGG_Reaction`, `KEGG_TC`, `CAZy` (verify spec's claim of cols 15, 18, 19)
- Token separator within each column (comma? semicolon? mix?)
- Prefix conventions (`ko:K00001` vs `K00001`; is `KEGG_Reaction` populated as `R00200` or `kegg.reaction:R00200`?)
- Sentinel values for empty fields (`-`? empty string?)

- [ ] **Step 7: Compile audit notes**

Open a temporary file (or just use comments in the next task) with sections matching the spec's "Schemas are tentative" note. For each schema:
- ✓ Confirmed as-spec
- ⚠️ Needs update: <what changed>

You'll feed these into Task 11.

---

### Task 11: Update Spec 1.1 with audit findings

**Files:**
- Modify: `docs/superpowers/specs/2026-04-28-metabolite-foundation-design.md` (inline)
- Modify: `multiomics_kg/download/download_metabolism_reference.py` (only if any URL was wrong)
- Modify: `.gitignore` (only if Task 10 file sizes warrant — see Step 3 below)

- [ ] **Step 1: Apply audit-driven spec edits**

For each `⚠️ Needs update` from Task 10, edit Spec 1.1 inline:

- SQLite schema section ([line 118 onwards in spec](../specs/2026-04-28-metabolite-foundation-design.md)) — adjust `compounds` / `compound_aliases` / `reactions` / `reaction_aliases` columns to match `chem_prop.tsv` / `chem_xref.tsv` reality. If MNX uses combined `source:value` strings rather than separate columns, document the parsing convention here.
- TCDB hierarchy section — replace the speculative shape with what the actual file gives us. If `substrate_classes` isn't available from this file, drop it from the schema (Spec 1.3 will need to find another source).
- CAZy hierarchy section — same.
- eggNOG columns table in "Step 2" section — confirm or correct column indices and note any prefix-stripping the transforms will need to do.

- [ ] **Step 2: Drop or downgrade the "Schemas are tentative" callout**

If the audit confirms most of the spec, change the callout from `> **Schemas are tentative.**` to `> **Schemas confirmed against MNX release <X>, TCDB <date>, CAZy <date>.**` followed by the residual unknowns (if any). If significant rework was needed, leave the callout in but note which specific parts are now confirmed.

- [ ] **Step 3: Decide on `.gitignore` for cache files**

Project convention: most `cache/data/` content **is committed** (genome data, gene_annotations_merged.json, eggNOG output, etc. — ~530 tracked files at the time of this plan). The decision for MNX/TCDB/CAZy:

- **Commit** if total bundle ≤ ~50 MB and files refresh < quarterly (matches existing committed reference data like `cache/data/eggnog/og_descriptions.json`, `cache/data/pfam/pfam_reference.json`).
- **Gitignore** if total bundle > 100 MB or files refresh monthly. Add:

  ```
  # Phase 1.1 metabolism reference data (downloaded by sub-step 6)
  cache/data/mnx/
  cache/data/tcdb/
  cache/data/cazy/
  ```

  And document the re-fetch command (`uv run python multiomics_kg/download/download_genome_data.py --steps 6 --force`) under CLAUDE.md "Data Locations".

Borderline (50–100 MB): default to commit; the fetch URL might change before the next refresh, and committed bytes are reproducible without depending on upstream availability.

- [ ] **Step 4: If any URL was wrong, fix `SOURCES`**

Edit `multiomics_kg/download/download_metabolism_reference.py` `SOURCES` dict with the corrected URL(s) found in Task 9 Step 3. Re-run `uv run python multiomics_kg/download/download_genome_data.py --steps 6 --force` and verify the file lands.

- [ ] **Step 5: Commit the audit-driven changes**

```bash
git add docs/superpowers/specs/2026-04-28-metabolite-foundation-design.md
# include any of these only if you actually edited them:
git add multiomics_kg/download/download_metabolism_reference.py 2>/dev/null
git add .gitignore CLAUDE.md 2>/dev/null
git commit -m "metabolism: confirm Phase 1.1A schemas via audit; update spec inline"
```

---

## Done state

After all 11 tasks:

- ✅ `cache/data/{mnx,tcdb,cazy}/` populated with real reference files.
- ✅ Sub-step 6 (download) runs end-to-end via `prepare_data.sh` step 0 OR via `download_genome_data.py --steps 6`.
- ✅ Sub-step 7 dispatches but raises `NotImplementedError: Phase 1.1B` (intentional — implementation is the next plan).
- ✅ Three accessor modules (`metabolite_utils`, `tcdb_utils`, `cazy_utils`) importable; all functions stubbed.
- ✅ Spec 1.1 internally consistent and confirmed against actual files.
- ✅ All non-slow tests pass: `pytest -m "not slow and not kg" -v`.
- ✅ No KG impact: deployed graph byte-identical.

**Next plan (Phase 1.1B):** implement sub-step 7 (resolver + hierarchies) and the three accessor modules and the gene-annotation transforms, against the now-confirmed schemas. Spec 1.1 is the contract; tests drive each module.
