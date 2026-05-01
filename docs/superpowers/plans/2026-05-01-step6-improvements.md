# Phase 1.2.2 — Step 6 Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Three independent improvements to the Phase 1.2.1 metabolism layer: (1) move MNX download + resolver build out of `prepare_data.sh` into a standalone refresh script (TCDB/CAZy stays in step 0), (2) make `kegg_data.json` git-friendly via `indent=2 + sort_keys=True`, (3) replace per-entity SQLite enrichment in step 6 with batched IN-clause queries (~30-60s → <2s).

**Architecture:**
- **MNX extraction**: New `scripts/refresh_mnx.sh` calls `download_metabolism_reference.py --sources mnx` + `build_metabolite_resolver.py`. The `--sources` flag is added to `download_metabolism_reference.py` (currently downloads all 7 sources unconditionally). Step 0 sub-step 6 is renamed to "Download TCDB reference data" and runs only the 3 TCDB URLs. Step 0 sub-step 7 (resolver build) is removed entirely. The `prepare_data.sh` orchestrator only invokes step 0 sub-steps 1-6 now.
- **JSON formatting**: One-line change in `build_pruned_kegg_data` to `json.dumps(out, indent=2, sort_keys=True)`. Keys are already sorted at creation time, but `sort_keys` adds a guarantee. File grows from 2 MB → ~3-4 MB.
- **Vector enrichment**: Refactor `_enrich_reaction` and `_enrich_compound` from per-entity calls (~16K queries) into a 2-phase batched approach. Phase 1: bulk-resolve all KEGG IDs → MNX IDs in one batched IN query per entity type. Phase 2: bulk-fetch all MNX-side properties (`reactions.classifs/is_balanced/is_transport`, `compounds.formula/mass/...`, `compound_aliases.chebi/hmdb`, `reaction_aliases.rhea`) using batched IN queries on MNX IDs. Build dict lookups, then assemble per-entity output dicts in pure Python.

**Tech Stack:** Python 3.10+, `sqlite3` (stdlib), `pytest` with `tmp_path` fixtures, the existing MNX SQLite resolver schema.

**Spec basis:** Conversation 2026-05-01 (post-Phase-1.2.1 ship).

**Depends on:** Phase 1.2.1 complete (commits `0ad040b3` → `9793a996` on branch `metabolite-scaffold`).

---

## Tasks must run strictly 1 → 2 → 3 → 4 → 5

Tasks 1 and 2 can run in either order (independent), but Task 3 (vector enrichment) and Task 4 (JSON formatting) must precede Task 5 (regenerate `kegg_data.json` + final snapshot diff).

---

## File structure

**Create:**

| Path | Why |
|---|---|
| `scripts/refresh_mnx.sh` | Standalone wrapper for MNX download + resolver build; documents that this only needs to run when MNX releases a new version |

**Modify:**

| Path | Why |
|---|---|
| `multiomics_kg/download/download_metabolism_reference.py` | Add `--sources` flag (mnx, tcdb, all). Rename module/file: keep filename for now (downstream paths) but split SOURCES into MNX_SOURCES + TCDB_SOURCES dicts. |
| `multiomics_kg/download/download_genome_data.py` | Step 0 sub-step 6 calls `download_all(sources=["tcdb"])` instead of all sources. Remove sub-step 7 entirely (resolver build no longer in step 0). Update `--steps` choices, default, and help text. |
| `scripts/prepare_data.sh` | Header docs reflect step-0 sub-step 6 scope change. The wrapping step-6 entry that calls `build_kegg_metabolism_xrefs` already documents an MNX-resolver dependency in its preconditions docstring; just update wording from "step 0 sub-steps 6+7" to "scripts/refresh_mnx.sh". |
| `multiomics_kg/download/build_kegg_metabolism_xrefs.py` | (a) Change `json.dumps(out, separators=(",", ":"))` → `json.dumps(out, indent=2, sort_keys=True)`. (b) Replace `_enrich_reaction` + `_enrich_compound` per-entity loops with batched bulk-fetch helpers. |
| `tests/test_build_kegg_metabolism_xrefs.py` | Update fixtures + add new tests for batched enrichment. Existing 3 tests should still pass against new vectorized code. |
| `cache/data/kegg/kegg_data.json` | Regenerate with new compact-but-pretty JSON shape (Task 5). |
| `CLAUDE.md` | Update "Step 0" sub-steps list (drop 7, change 6 to TCDB-only). Add note about `scripts/refresh_mnx.sh`. Update `kegg_data.json` size mention from "~2 MB" to "~3-4 MB". |

**Delete:** none.

---

## Task 1: Add `--sources` flag to `download_metabolism_reference.py`

**Files:**
- Modify: `multiomics_kg/download/download_metabolism_reference.py`
- Modify: `tests/test_download_genome_data.py` (if it tests metabolism reference download — check first)

The current module downloads all 7 sources unconditionally. We need it to optionally download a subset.

- [ ] **Step 1: Read the current module**

Read `multiomics_kg/download/download_metabolism_reference.py` end-to-end (~75 lines). Note the `SOURCES` dict (lines 23-31) which is a single flat dict.

- [ ] **Step 2: Split SOURCES into MNX_SOURCES + TCDB_SOURCES**

Replace the `SOURCES` dict (lines 23-31) with two dicts:

```python
MNX_SOURCES: dict[str, tuple[str, str]] = {
    "mnx_chem_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",       "mnx/chem_prop.tsv"),
    "mnx_chem_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",       "mnx/chem_xref.tsv"),
    "mnx_reac_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",       "mnx/reac_prop.tsv"),
    "mnx_reac_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",       "mnx/reac_xref.tsv"),
}

TCDB_SOURCES: dict[str, tuple[str, str]] = {
    "tcdb_families":      ("https://www.tcdb.org/cgi-bin/projectv/public/families.py",           "tcdb/families.tsv"),
    "tcdb_substrates":    ("https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py",           "tcdb/substrates.tsv"),
    "tcdb_superfamilies": ("https://www.tcdb.org/cgi-bin/substrates/listSuperfamilies.py",       "tcdb/superfamilies.tsv"),
}

SOURCES_BY_GROUP: dict[str, dict[str, tuple[str, str]]] = {
    "mnx": MNX_SOURCES,
    "tcdb": TCDB_SOURCES,
}
```

- [ ] **Step 3: Add `sources` parameter to `download_all`**

Replace `download_all(cache_root, force)` with:

```python
def download_all(
    cache_root: Path = DEFAULT_CACHE_ROOT,
    force: bool = False,
    sources: list[str] | None = None,
) -> None:
    """Download requested source groups into cache_root.

    Args:
        cache_root: where to write files (default cache/data)
        force: re-download even if cached
        sources: list of source-group keys (e.g. ['mnx', 'tcdb']); None = all groups
    """
    selected_groups = sources or list(SOURCES_BY_GROUP.keys())
    invalid = [s for s in selected_groups if s not in SOURCES_BY_GROUP]
    if invalid:
        raise ValueError(
            f"Unknown source group(s): {invalid}. "
            f"Valid: {sorted(SOURCES_BY_GROUP.keys())}"
        )

    log.info(
        f"download_metabolism_reference: cache_root={cache_root} "
        f"force={force} groups={selected_groups}"
    )
    n_downloaded = 0
    n_total = 0
    for group in selected_groups:
        for key, (url, rel) in SOURCES_BY_GROUP[group].items():
            n_total += 1
            if download_one(url, cache_root / rel, force=force):
                n_downloaded += 1
    log.info(f"  done — {n_downloaded}/{n_total} downloaded, "
             f"{n_total - n_downloaded} cached.")
```

- [ ] **Step 4: Update `main()` and CLI**

Replace the existing `main()` and `__main__` block:

```python
def main(force: bool = False, sources: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    download_all(force=force, sources=sources)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if cache files exist.")
    parser.add_argument(
        "--sources", nargs="+", choices=sorted(SOURCES_BY_GROUP.keys()),
        default=None,
        help=f"Source groups to download (default: all = {sorted(SOURCES_BY_GROUP.keys())})",
    )
    args = parser.parse_args()
    main(force=args.force, sources=args.sources)
```

- [ ] **Step 5: Update module docstring**

Replace the current docstring at the top:

```python
"""Download MNX and/or TCDB reference data.

Default behaviour downloads both groups (4 MNX TSVs + 3 TCDB TSVs). Use
``--sources mnx`` or ``--sources tcdb`` to download a single group. MNX dominates
total size (~1.5 GB unzipped); TCDB tables are <1 MB each. Files are cached
under cache/data/{mnx,tcdb}/ and skipped on re-run unless --force.

Used by:
- prepare_data.sh step 0 sub-step 6 (TCDB only — MNX moved to scripts/refresh_mnx.sh)
- scripts/refresh_mnx.sh (MNX only)

CAZy family hierarchy is NOT downloaded here. It is bootstrapped in Phase 1.1B
from observed eggNOG `CAZy` columns plus mechanical ID parsing — the format
itself encodes the hierarchy (`GH13_1` → family `GH13` → class `GH`), and we
only care about families our genes actually reference.
"""
```

- [ ] **Step 6: Add unit tests**

Read existing tests for `download_metabolism_reference.py` first:
```bash
grep -rln "download_metabolism_reference\|download_all" tests/
```

If no existing test file: create `tests/test_download_metabolism_reference.py` with this content:

```python
"""Unit tests for download_metabolism_reference module."""
from __future__ import annotations
from pathlib import Path
from unittest.mock import patch

import pytest

from multiomics_kg.download import download_metabolism_reference as dmr


def test_sources_by_group_keys():
    assert set(dmr.SOURCES_BY_GROUP.keys()) == {"mnx", "tcdb"}
    assert len(dmr.MNX_SOURCES) == 4
    assert len(dmr.TCDB_SOURCES) == 3


def test_download_all_defaults_to_all_groups(tmp_path, monkeypatch):
    """No --sources arg → downloads both groups."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=None)
    assert len(called) == 7  # 4 MNX + 3 TCDB


def test_download_all_filters_to_mnx(tmp_path, monkeypatch):
    """--sources mnx → only MNX URLs."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=["mnx"])
    assert len(called) == 4
    assert all("metanetx.org" in u for u in called)


def test_download_all_filters_to_tcdb(tmp_path, monkeypatch):
    """--sources tcdb → only TCDB URLs."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=["tcdb"])
    assert len(called) == 3
    assert all("tcdb.org" in u for u in called)


def test_download_all_invalid_source_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown source group"):
        dmr.download_all(cache_root=tmp_path, sources=["bogus"])
```

- [ ] **Step 7: Run tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_download_metabolism_reference.py -v
```
Expected: 5 PASS.

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/download/download_metabolism_reference.py \
        tests/test_download_metabolism_reference.py
git commit -m "metabolism: 1.2.2 — split metabolism reference download into mnx/tcdb groups"
```

---

## Task 2: Standalone MNX refresh script + remove from step 0

**Files:**
- Create: `scripts/refresh_mnx.sh`
- Modify: `multiomics_kg/download/download_genome_data.py`

- [ ] **Step 1: Create `scripts/refresh_mnx.sh`**

Create with this content:

```bash
#!/usr/bin/env bash
# Refresh MNX (MetaNetX) cross-reference cache + SQLite resolver.
#
# Run this ONLY when MNX releases a new version. The resolver is heavy (~2.5 GB
# SQLite, ~30 minutes wall time) and does not depend on which strains are in
# the project — it's a global metabolism cross-reference lookup.
#
# After this completes, run `bash scripts/prepare_data.sh --steps 6 --force`
# to regenerate cache/data/kegg/kegg_data.json with the new MNX cross-refs.
#
# Outputs:
#   cache/data/mnx/{chem_prop,chem_xref,reac_prop,reac_xref}.tsv  (~1.5 GB raw)
#   cache/data/mnx/metabolite_resolver.db                         (~2.5 GB SQLite)
#
# Usage:
#   bash scripts/refresh_mnx.sh           # download + build (skip if cached)
#   bash scripts/refresh_mnx.sh --force   # re-download + rebuild

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$LOG_DIR"

FORCE=""
if [[ "${1:-}" == "--force" ]]; then
    FORCE="--force"
fi

cd "$PROJECT_ROOT"

echo "══════════════════════════════════════════════════════════════════════"
echo "  refresh_mnx.sh — MNX TSV download + SQLite resolver build"
echo "  Log: $LOG_DIR/refresh_mnx.log"
echo "══════════════════════════════════════════════════════════════════════"

{
    echo ""
    echo "── 1. Download MNX TSVs (4 files, ~1.5 GB) ──"
    uv run python -m multiomics_kg.download.download_metabolism_reference \
        --sources mnx $FORCE

    echo ""
    echo "── 2. Build metabolite resolver SQLite (~30 min, ~2.5 GB output) ──"
    uv run python -m multiomics_kg.download.build_metabolite_resolver $FORCE
} 2>&1 | tee "$LOG_DIR/refresh_mnx.log"

echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo "  refresh_mnx.sh complete. Next: bash scripts/prepare_data.sh --steps 6 --force"
echo "══════════════════════════════════════════════════════════════════════"
```

Make executable:
```bash
chmod +x scripts/refresh_mnx.sh
```

- [ ] **Step 2: Update `download_genome_data.py` step 6 (remove MNX) and remove step 7**

Open `multiomics_kg/download/download_genome_data.py`. Three changes:

(A) Replace the body of `step6_metabolism_reference` (around lines 445-449) with:

```python
def step6_tcdb_reference(force: bool) -> None:
    """Sub-step 6: download TCDB reference data (MNX moved to scripts/refresh_mnx.sh)."""
    from multiomics_kg.download.download_metabolism_reference import download_all
    log.info("─── Step 6: Download TCDB reference data ───")
    download_all(force=force, sources=["tcdb"])
```

Rename the function from `step6_metabolism_reference` → `step6_tcdb_reference` and update the call site (search for the old function name in the same file — there's a dispatch around line 540-ish).

(B) DELETE `step7_metabolite_resolver` (around lines 452-456) entirely.

(C) Update the `--steps` argparse choices (around line 484) from `[1, 2, 3, 4, 5, 6, 7]` to `[1, 2, 3, 4, 5, 6]`. Update the default from `[1, 2, 3, 4, 5, 6, 7]` to `[1, 2, 3, 4, 5, 6]`. Update the `--steps` help text from `6=metabolism_reference 7=metabolite_resolver` to `6=tcdb_reference`.

(D) Update the docstring near the top (around lines 21 and 30) to remove mention of step 7 and update step 6 to "Download TCDB reference data".

(E) Update the dispatch loop (search for `if 7 in steps:` and `if 6 in steps:`) — remove the step-7 dispatch entirely, and update the step-6 dispatch to call `step6_tcdb_reference(force=args.force)`.

(F) Update the epilog Examples block (around lines 477-480): remove the `--steps 6 7 --force` example or change it to `--steps 6 --force`.

- [ ] **Step 3: Verify the module imports + parses**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run python -c "from multiomics_kg.download import download_genome_data; print('OK')"
```
Expected: prints `OK`, no traceback.

```bash
uv run python -m multiomics_kg.download.download_genome_data --help 2>&1 | head -30
```
Expected: help text shows steps 1-6 (no 7), and step 6 is described as TCDB.

- [ ] **Step 4: Test the standalone script's syntax**

```bash
bash -n scripts/refresh_mnx.sh && echo "syntax OK"
```
Expected: prints `syntax OK`.

(Don't actually run `refresh_mnx.sh` end-to-end — it would re-download 1.5 GB and take 30 minutes. Manual verification.)

- [ ] **Step 5: Update existing test for `step7_metabolite_resolver` if any**

Check:
```bash
grep -rn "step7_metabolite_resolver\|step6_metabolism_reference" tests/
```
If found, either remove the test (for step7) or update the test (for step6 — rename to `step6_tcdb_reference` and adjust expected source filter).

- [ ] **Step 6: Run the full unit suite to catch any breakage**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest -m "not slow and not kg" -q 2>&1 | tail -5
```
Expected: all PASS (1726+).

- [ ] **Step 7: Commit**

```bash
git add scripts/refresh_mnx.sh multiomics_kg/download/download_genome_data.py
# Plus any tests if updated:
# git add tests/test_download_genome_data.py
git commit -m "metabolism: 1.2.2 — extract MNX download+resolver to scripts/refresh_mnx.sh"
```

---

## Task 3: Vectorize step-6 enrichment (batched IN-clause queries)

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Modify: `tests/test_build_kegg_metabolism_xrefs.py`

The existing `_enrich_reaction(rxn_id, ...)` and `_enrich_compound(cpd_id, ...)` make 3 and 4 SQL queries per entity respectively. For 2,340 reactions + 2,184 compounds = ~16K SQL queries. Each ~1-3 ms with SQLite → 30-60s total runtime.

Replace with a 2-phase batched approach:
- **Phase 1**: bulk-resolve KEGG IDs → MNX IDs (one batched IN query per entity type)
- **Phase 2**: bulk-fetch all MNX-side properties (reaction classifications, balance, transport flags; compound formula/mass/inchikey/smiles; chebi/hmdb aliases; rhea aliases) — one batched IN query per attribute, indexed by MNX ID

Then assemble per-entity dicts in pure Python.

Helper utility: `_chunked_in_query(conn, sql_template, ids, chunk_size=900)` — SQLite has a default limit of 999 host parameters, so chunk at 900 to be safe.

- [ ] **Step 1: Read the current enrichment code**

Read `multiomics_kg/download/build_kegg_metabolism_xrefs.py` lines 96-186 (the `_enrich_reaction` and `_enrich_compound` functions). Note the SQL queries used.

- [ ] **Step 2: Write failing tests for the new batched API**

Open `tests/test_build_kegg_metabolism_xrefs.py`. After the existing 3 tests, add 2 new tests for the new batched helpers. Keep the existing 3 tests as-is (they're end-to-end against `build_pruned_kegg_data` and should still pass after refactor).

```python
def test_bulk_enrich_reactions_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_reactions(conn, kegg_ids, allowed_pathways, raw) returns
    a dict mapping kegg_reaction_id → enrichment dict, with the same fields
    as the per-entity _enrich_reaction.
    """
    conn = _make_resolver(tmp_path)
    raw = {
        "reaction_names": {"R00200": "ATP:pyruvate ..."},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
        "reaction_to_compounds": {"R00200": ["C00031"]},
    }
    allowed = {"ko00010", "ko00710"}

    result = bx._bulk_enrich_reactions(conn, ["R00200", "R99999"], allowed, raw)

    # Every requested ID is in the result (R99999 has no MNX entry but still gets a stub)
    assert set(result.keys()) == {"R00200", "R99999"}

    r = result["R00200"]
    assert r["name"].startswith("ATP:pyruvate")
    assert r["mnxr_id"] == "MNXR101234"
    assert r["ec_numbers"] == ["2.7.1.40"]
    assert r["mass_balance"] == "balanced"
    assert r["reaction_class"] == "chemical"
    assert r["pathways"] == ["ko00010", "ko00710"]
    assert r["compounds"] == ["C00031"]

    r99 = result["R99999"]
    assert r99["mnxr_id"] is None
    assert r99["ec_numbers"] == []
    assert r99["mass_balance"] == "unbalanced"  # default
    assert r99["reaction_class"] == "chemical"  # default


def test_bulk_enrich_compounds_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_compounds returns a dict mapping kegg_compound_id → enrichment."""
    conn = _make_resolver(tmp_path)
    raw = {
        "compound_names": {"C00031": "D-glucose", "C99999": "obscure"},
        "compound_to_pathways": {
            "C00031": ["ko00010", "ko00500"],
            "C99999": [],
        },
    }
    allowed = {"ko00010"}

    result = bx._bulk_enrich_compounds(conn, ["C00031", "C99999"], allowed, raw)

    assert set(result.keys()) == {"C00031", "C99999"}

    glucose = result["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["pathways"] == ["ko00010"]  # ko00500 filtered out by allowed

    obscure = result["C99999"]
    assert obscure["mnxm_id"] is None
    assert obscure["chebi_id"] is None
    assert obscure["formula"] is None
```

Run them:
```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v
```
Expected: 3 EXISTING PASS, 2 NEW FAIL with `AttributeError: module 'build_kegg_metabolism_xrefs' has no attribute '_bulk_enrich_reactions'`.

- [ ] **Step 3: Add the chunked-IN helper**

Insert into `multiomics_kg/download/build_kegg_metabolism_xrefs.py` near the top (after imports, before `_gene_reachable_sets`):

```python
# SQLite default: 999 host parameters per statement. Use 900 for safety.
_BATCH_SIZE = 900


def _batched(seq: list, n: int = _BATCH_SIZE):
    """Yield successive n-sized chunks from seq."""
    for i in range(0, len(seq), n):
        yield seq[i:i + n]
```

- [ ] **Step 4: Add `_bulk_enrich_reactions`**

Insert after the `_batched` helper (replacing the old `_enrich_reaction` location around lines 96-135):

```python
def _bulk_enrich_reactions(
    conn: sqlite3.Connection,
    kegg_ids: list[str],
    allowed_pathways: set[str],
    raw: dict,
) -> dict[str, dict]:
    """Batched enrichment for all KEGG reactions in one go.

    Returns dict mapping kegg_reaction_id → enrichment dict (same shape as the
    per-entity _enrich_reaction output). Reactions with no MNX entry get a stub
    with mnxr_id=None and default fields, so every requested ID appears in the
    output.

    Replaces ~3 queries × N reactions with ~3 batched queries total (per chunk).
    """
    cur = conn.cursor()

    # Initialize stubs for every requested reaction
    rxn_names = raw.get("reaction_names", {})
    rxn_to_pw = raw.get("reaction_to_pathways", {})
    rxn_to_cpds = raw.get("reaction_to_compounds", {})

    out: dict[str, dict] = {}
    for rxn_id in kegg_ids:
        out[rxn_id] = {
            "name": rxn_names.get(rxn_id, ""),
            "pathways": [p for p in rxn_to_pw.get(rxn_id, []) if p in allowed_pathways],
            "compounds": list(rxn_to_cpds.get(rxn_id, [])),
            "ec_numbers": [],
            "mnxr_id": None,
            "rhea_ids": [],
            "mass_balance": "unbalanced",
            "reaction_class": "chemical",
        }

    # Phase 1: bulk-resolve kegg.reaction → mnxr_id
    kegg_to_mnxr: dict[str, str] = {}
    for chunk in _batched(kegg_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxr_id FROM reaction_aliases "
            f"WHERE source='kegg.reaction' AND value IN ({placeholders})",
            chunk,
        )
        for value, mnxr_id in cur.fetchall():
            kegg_to_mnxr[value] = mnxr_id

    # Apply phase-1 results
    for rxn_id, mnxr_id in kegg_to_mnxr.items():
        out[rxn_id]["mnxr_id"] = mnxr_id

    if not kegg_to_mnxr:
        return out

    # Phase 2: bulk-fetch MNX-side properties keyed by mnxr_id
    mnxr_ids = list(set(kegg_to_mnxr.values()))

    # 2a. reactions table → ec_numbers, mass_balance, reaction_class
    rxn_props: dict[str, dict] = {}
    for chunk in _batched(mnxr_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxr_id, classifs, is_balanced, is_transport "
            f"FROM reactions WHERE mnxr_id IN ({placeholders})",
            chunk,
        )
        for mnxr_id, classifs, is_balanced, is_transport in cur.fetchall():
            rxn_props[mnxr_id] = {
                "ec_numbers": (
                    [c.strip() for c in classifs.split(";") if c.strip()]
                    if classifs else []
                ),
                "mass_balance": "balanced" if (is_balanced or "").upper() == "B" else "unbalanced",
                "reaction_class": "transport" if (is_transport or "").upper() == "T" else "chemical",
            }

    # 2b. reaction_aliases → rhea_ids
    rhea_by_mnxr: dict[str, list[str]] = {}
    for chunk in _batched(mnxr_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxr_id, value FROM reaction_aliases "
            f"WHERE source='rhea' AND mnxr_id IN ({placeholders}) ORDER BY value",
            chunk,
        )
        for mnxr_id, value in cur.fetchall():
            rhea_by_mnxr.setdefault(mnxr_id, []).append(value)

    # Stitch phase-2 results into per-reaction output
    for rxn_id, mnxr_id in kegg_to_mnxr.items():
        props = rxn_props.get(mnxr_id, {})
        out[rxn_id]["ec_numbers"] = props.get("ec_numbers", [])
        out[rxn_id]["mass_balance"] = props.get("mass_balance", "unbalanced")
        out[rxn_id]["reaction_class"] = props.get("reaction_class", "chemical")
        out[rxn_id]["rhea_ids"] = rhea_by_mnxr.get(mnxr_id, [])

    return out
```

- [ ] **Step 5: Add `_bulk_enrich_compounds`**

Insert after `_bulk_enrich_reactions`:

```python
def _bulk_enrich_compounds(
    conn: sqlite3.Connection,
    kegg_ids: list[str],
    allowed_pathways: set[str],
    raw: dict,
) -> dict[str, dict]:
    """Batched enrichment for all KEGG compounds.

    Same pattern as _bulk_enrich_reactions: phase 1 resolves kegg.compound →
    mnxm_id, phase 2 bulk-fetches MNX-side properties (formula/mass/inchikey/
    smiles + chebi/hmdb aliases), then assembles per-compound output.
    """
    cur = conn.cursor()

    cpd_names = raw.get("compound_names", {})
    cpd_to_pw = raw.get("compound_to_pathways", {})

    out: dict[str, dict] = {}
    for cpd_id in kegg_ids:
        out[cpd_id] = {
            "name": cpd_names.get(cpd_id, ""),
            "formula": None, "mass": None, "inchikey": None, "smiles": None,
            "mnxm_id": None, "chebi_id": None, "hmdb_id": None,
            "pathways": [p for p in cpd_to_pw.get(cpd_id, []) if p in allowed_pathways],
        }

    # Phase 1: bulk-resolve kegg.compound → mnxm_id
    kegg_to_mnxm: dict[str, str] = {}
    for chunk in _batched(kegg_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxm_id FROM compound_aliases "
            f"WHERE source='kegg.compound' AND value IN ({placeholders})",
            chunk,
        )
        for value, mnxm_id in cur.fetchall():
            kegg_to_mnxm[value] = mnxm_id

    for cpd_id, mnxm_id in kegg_to_mnxm.items():
        out[cpd_id]["mnxm_id"] = mnxm_id

    if not kegg_to_mnxm:
        return out

    mnxm_ids = list(set(kegg_to_mnxm.values()))

    # 2a. compounds table → formula, mass, inchikey, smiles
    cpd_props: dict[str, dict] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, formula, mass, inchikey, smiles "
            f"FROM compounds WHERE mnxm_id IN ({placeholders})",
            chunk,
        )
        for mnxm_id, formula, mass, inchikey, smiles in cur.fetchall():
            cpd_props[mnxm_id] = {
                "formula": formula or None,
                "mass": mass,  # numeric, no `or None` coercion needed
                "inchikey": inchikey or None,
                "smiles": smiles or None,
            }

    # 2b. compound_aliases → chebi_id (one lowest-id per mnxm)
    chebi_by_mnxm: dict[str, str] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='chebi' AND mnxm_id IN ({placeholders}) ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            # Keep first (lowest-sorted) chebi id per mnxm
            chebi_by_mnxm.setdefault(mnxm_id, value)

    # 2c. compound_aliases → hmdb_id (same pattern)
    hmdb_by_mnxm: dict[str, str] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='hmdb' AND mnxm_id IN ({placeholders}) ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            hmdb_by_mnxm.setdefault(mnxm_id, value)

    # Stitch phase-2 results into per-compound output
    for cpd_id, mnxm_id in kegg_to_mnxm.items():
        props = cpd_props.get(mnxm_id, {})
        out[cpd_id]["formula"] = props.get("formula")
        out[cpd_id]["mass"] = props.get("mass")
        out[cpd_id]["inchikey"] = props.get("inchikey")
        out[cpd_id]["smiles"] = props.get("smiles")
        out[cpd_id]["chebi_id"] = chebi_by_mnxm.get(mnxm_id)
        out[cpd_id]["hmdb_id"] = hmdb_by_mnxm.get(mnxm_id)

    return out
```

- [ ] **Step 6: Update `build_pruned_kegg_data` to call the bulk helpers**

Find `build_pruned_kegg_data` in the same file (around lines 191-234). Replace the `reactions` and `compounds` dict-comprehensions:

```python
        "reactions": {
            r: _enrich_reaction(r, raw, conn, pws) for r in sorted(rxns)
        },
        "compounds": {
            c: _enrich_compound(c, raw, conn, pws) for c in sorted(cpds)
        },
```

with:

```python
    rxn_enriched = _bulk_enrich_reactions(conn, sorted(rxns), pws, raw)
    cpd_enriched = _bulk_enrich_compounds(conn, sorted(cpds), pws, raw)
```

Then in the `out = {...}` dict literal, the `"reactions"` and `"compounds"` keys become:

```python
        "reactions": rxn_enriched,
        "compounds": cpd_enriched,
```

(Note: `_bulk_enrich_*` already return dicts keyed by KEGG ID; if `sorted(...)` order matters for downstream, Python 3.7+ preserves insertion order — and `_bulk_enrich_*` insert in the order of the input `kegg_ids` list. We pass `sorted(rxns)` and `sorted(cpds)`, so the resulting dicts are in sorted order.)

- [ ] **Step 7: Delete the old per-entity helpers**

Delete `_enrich_reaction` (around lines 96-135) and `_enrich_compound` (around lines 138-186). They're replaced.

- [ ] **Step 8: Run the unit tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v
```
Expected: 5 PASS (3 existing end-to-end + 2 new bulk-helper tests).

If the existing 3 tests fail because they reference `_enrich_reaction` / `_enrich_compound` symbols that no longer exist, update those test usages too.

- [ ] **Step 9: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        tests/test_build_kegg_metabolism_xrefs.py
git commit -m "metabolism: 1.2.2 — vectorize step-6 enrichment via batched IN queries"
```

---

## Task 4: JSON formatting (indent=2 + sort_keys=True)

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Modify: `tests/test_build_kegg_metabolism_xrefs.py` (assertions on file contents may need adjusting if any check raw JSON bytes)

- [ ] **Step 1: Change the JSON dump line**

In `multiomics_kg/download/build_kegg_metabolism_xrefs.py`, find the line in `build_pruned_kegg_data` (around line 234):

```python
    out_path.write_text(json.dumps(out, separators=(",", ":")))
```

Replace with:

```python
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
```

- [ ] **Step 2: Run the existing tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v
```
Expected: 5 PASS. Tests use `json.loads(out_path.read_text())` which doesn't care about formatting.

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py
git commit -m "metabolism: 1.2.2 — kegg_data.json uses indent=2 + sort_keys=True (git-friendly)"
```

---

## Task 5: Regenerate `kegg_data.json` + verify

**Files:**
- Modify (regenerate): `cache/data/kegg/kegg_data.json`
- Modify: `CLAUDE.md`
- Modify: `scripts/prepare_data.sh` (header comment — drop step-0 sub-step 7 reference)

- [ ] **Step 1: Regenerate `kegg_data.json`**

```bash
cd /home/osnat/github/multiomics_biocypher_kg
bash scripts/prepare_data.sh --steps 6 --force 2>&1 | tail -20
```
Expected: completes successfully. Step-6 runtime should drop from ~30-60s to <5s thanks to the batched queries. New `kegg_data.json` will be ~3-4 MB (was 2 MB) due to the indented format.

- [ ] **Step 2: Verify size + structure**

```bash
ls -lh cache/data/kegg/kegg_data.json
uv run python -c "
import json
data = json.loads(open('cache/data/kegg/kegg_data.json').read())
print('Top-level keys (sorted):', sorted(data.keys()))
print('Counts:')
for k in ('kos','pathways','subcategories','categories','reactions','compounds'):
    print(f'  {k}: {len(data.get(k, {}))}')
# First key in each section should be lexically smallest (sort_keys=True)
print('First KO:', list(data['kos'].keys())[0])
print('First reaction:', list(data['reactions'].keys())[0])
print('First compound:', list(data['compounds'].keys())[0])
"
```
Expected: counts unchanged from Phase 1.2.1 baseline (~4534 KOs, ~377 pathways, ~46 subcategories, ~6 categories, ~2340 reactions, ~2184 compounds). First keys in lexical order (e.g., K00001, R00001, C00001).

- [ ] **Step 3: Eyeball the diff to confirm git-friendliness**

```bash
git diff --stat cache/data/kegg/kegg_data.json
head -50 cache/data/kegg/kegg_data.json
```
Expected: file shows clean indented structure. Git diff shows the file changed. Don't try to read the full diff — it's the entire file that changed.

- [ ] **Step 4: Run the full unit suite**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest -m "not slow and not kg" -q 2>&1 | tail -5
```
Expected: all PASS.

- [ ] **Step 5: Update `CLAUDE.md`**

Find these sections and update:

(A) **Step 0 sub-steps list** in "Genome Data Download Pipeline" section. Currently includes:
> - 6: Download MNX/TCDB/CAZy reference data → `cache/data/{mnx,tcdb,cazy}/`
> - 7: Build metabolite resolver + hierarchy caches (requires sub-step 6). Phase 1.1A: skeleton only — full build lands in Phase 1.1B.

Change to:
> - 6: Download TCDB reference data → `cache/data/tcdb/`. (MNX moved to standalone `scripts/refresh_mnx.sh` as of Phase 1.2.2 — only needs to run when MNX releases a new version.)

(B) **Step 6 description** (the top-level `prepare_data.sh` step 6, in the same Pipeline section). Find the existing description and update the file size from "~2 MB" to "~3-4 MB (indented JSON for git-friendly diffs)". Add a note about the standalone MNX script:

> ... writes a single `cache/data/kegg/kegg_data.json` (~3-4 MB, indented JSON for git-friendly diffs). Both `kegg_annotation_adapter` and `metabolism_adapter` read this file. The 2.6 GB MNX resolver opens here only — it must be built first via `bash scripts/refresh_mnx.sh` (one-time setup; rerun only when MNX releases). Run as module: `uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]`. Requires step 5 + the MNX resolver.

(C) Search for any other references to `step 7` in the context of MNX/metabolism in CLAUDE.md (e.g., in commands or examples). Update or remove as needed.

- [ ] **Step 6: Update `scripts/prepare_data.sh` header comment**

In `scripts/prepare_data.sh` header (the comment block at the top), find the description of step 6 (currently mentions MNX dependency on "step 0 sub-steps 6+7"). Change to "scripts/refresh_mnx.sh".

Look for this line:
```
# Step 6 — Build pruned KEGG data cache (kegg_data.json) with metabolism enrichment
#           ...
#           Requires step 0 sub-steps 6+7 (MNX/TCDB/CAZy reference + resolver) and step 2
```

Change to:
```
# Step 6 — Build pruned KEGG data cache (kegg_data.json) with metabolism enrichment
#           ...
#           Requires step 0 sub-step 6 (TCDB reference) + step 2 + scripts/refresh_mnx.sh
#           (MNX resolver — heavy one-time build, rerun only when MNX releases)
```

- [ ] **Step 7: Commit**

```bash
git add cache/data/kegg/kegg_data.json CLAUDE.md scripts/prepare_data.sh
git commit -m "metabolism: 1.2.2 — regenerate kegg_data.json (indented + vectorized) + docs"
```

---

## Validation gate (manual, after all 5 tasks land)

```bash
# 1. Re-confirm step 6 runtime improvement (compare to before)
time bash scripts/prepare_data.sh --steps 6 --force
# Expected: < 5s wall time (was ~30-60s)

# 2. Confirm CLAUDE.md is consistent — no orphan "step 7" / "MNX/TCDB/CAZy" refs
grep -nE "step 0.*7|MNX/TCDB" CLAUDE.md
# Expected: no hits (or only intentional ones in historical notes)

# 3. Run KG validity tests against the existing deployed graph (no rebuild needed —
#    kegg_data.json size/format change doesn't affect graph content, only the
#    on-disk JSON shape)
uv run pytest -m kg -v
# Expected: same pass count as before Phase 1.2.2

# 4. (Optional) Test the standalone MNX script's `--force` end-to-end on a clean cache
#    -- skip this unless validating the script itself; takes ~30 minutes
# bash scripts/refresh_mnx.sh --force
```

Acceptance:
- `kegg_data.json` is committed with `indent=2 + sort_keys=True` (eyeball: starts with `{\n  "categories": ...`)
- Step 6 runtime drops to <5s
- `scripts/refresh_mnx.sh` exists, is executable, and `bash -n` passes
- `download_genome_data.py --steps` no longer accepts `7` as a value
- All existing unit tests + KG validity tests still pass
- CLAUDE.md reflects the new Step 0 sub-step 6 scope and mentions `scripts/refresh_mnx.sh`

---

## Self-review checklist

1. **Spec coverage:** Each of the 3 user-stated improvements (MNX extraction, JSON git-friendliness, vector enrichment) maps to a task: Task 1+2 (MNX), Task 4 (JSON), Task 3 (vector). ✓
2. **No placeholders:** Every step has actual code or commands. ✓
3. **Type / name consistency:** `_bulk_enrich_reactions` / `_bulk_enrich_compounds` consistent. `step6_tcdb_reference` rename consistent. `SOURCES_BY_GROUP` keys (`mnx`, `tcdb`) used consistently in tests + `--sources` CLI choices. ✓
4. **Order discipline:** Tasks 3 + 4 must precede Task 5 (regen). Tasks 1 + 2 are independent of Tasks 3-5 but it's simpler to land them first to keep the diff focused. Sequential 1→2→3→4→5 works cleanly. ✓
5. **Rollback safety:** Each commit is independent. Task 1 alone is harmless (extends API). Task 2 removes step 7 from `prepare_data.sh`'s default flow but keeps the underlying `build_metabolite_resolver.py` callable. Task 3 is a refactor with same I/O contract (verified by existing 3 end-to-end tests). Task 4 is a one-line format change. Task 5 regenerates the data file using the post-Task-3+4 code. ✓
