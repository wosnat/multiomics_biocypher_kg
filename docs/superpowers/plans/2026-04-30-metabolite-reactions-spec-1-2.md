# Phase 1.2 — Reactions / Metabolites scaffold (KEGG-native)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the chemistry layer of the metabolism scaffold: `Reaction` and `Metabolite` nodes plus the four edges (`Gene_catalyzes_reaction`, `Reaction_has_metabolite`, `Reaction_in_kegg_pathway`, `Organism_has_metabolite`) that connect them to genes, pathways, and organisms.

**Architecture:** New step 6 in the prepare-data pipeline (`build_kegg_metabolism_xrefs.py`) walks every strain's `gene_annotations_merged.json` to identify gene-reachable KEGG R-numbers, derives reachable C-numbers via KEGG `/link/compound/reaction`, and enriches each with KEGG metadata + MNX xrefs into a small JSON cache (`cache/data/kegg/kegg_metabolism_xrefs.json`, ~1 MB). The build-time adapter (`metabolism_adapter.py`) is a pure file reader — no SQLite, no KEGG REST. The MNX 2.6 GB resolver opens at *prepare-data* time only, mirroring the `build_og_descriptions.py` step-5 pattern.

**Tech Stack:** Python 3.10+, `sqlite3` (stdlib), KEGG REST API (`https://rest.kegg.jp`), `requests`, `pytest` with `tmp_path` for synthetic fixtures, BioCypher schema framework, the existing `multiomics_kg/utils/{kegg_utils,metabolite_utils}.py` modules.

**Spec:** [`docs/superpowers/specs/2026-04-30-metabolite-reactions-scaffold-revised.md`](../specs/2026-04-30-metabolite-reactions-scaffold-revised.md) (KEGG-native; supersedes the 2026-04-28 spec).

**Depends on:** Phase 1.1B (resolver + accessors + transforms — DONE on branch `metabolite-scaffold`, 15 commits between `8743e0e` and `0102115`).

---

## File structure

**Create:**

| Path | Role |
|---|---|
| `multiomics_kg/download/build_kegg_metabolism_xrefs.py` | Step 6: prune-then-enrich KEGG R/C IDs gene-reachable from any strain |
| `multiomics_kg/adapters/metabolism_adapter.py` | Build-time adapter emitting Reaction + Metabolite nodes + 3 edge types |
| `tests/test_kegg_utils_metabolism.py` | Unit tests for the 5 new KEGG endpoint parsers |
| `tests/test_metabolite_utils_primary_id.py` | Unit tests for `mnxm_to_primary_id()` / `mnxr_to_primary_id()` |
| `tests/test_build_kegg_metabolism_xrefs.py` | Unit tests for step 6 |
| `tests/test_metabolism_adapter.py` | Unit tests for the metabolism adapter |
| `tests/kg_validity/test_metabolism.py` | KG validity tests against the deployed graph |
| `tests/test_metabolism_smoke.py` | Slow integration test: real cache → adapter → output |

**Modify:**

| Path | Why |
|---|---|
| `multiomics_kg/utils/kegg_utils.py` | Add 5 endpoint constants + 5 parsers + extend `download_kegg_data()` |
| `multiomics_kg/download/utils/annotation_transforms.py` | Drop `_tx_resolve_kegg_reaction_to_mnxr` + registry entry |
| `config/gene_annotations_config.yaml` | Drop `transform: resolve_kegg_reaction_to_mnxr` line on `kegg_reactions` |
| `tests/test_annotation_transforms_metabolism.py` | Drop the KEGG-reaction transform test cases (keep TCDB + CAZy) |
| `multiomics_kg/download/build_metabolite_resolver.py` | Add 2 indexes (`idx_compound_aliases_mnxm`, `idx_reaction_aliases_mnxr`) |
| `multiomics_kg/utils/metabolite_utils.py` | Add `mnxm_to_primary_id()` + `mnxr_to_primary_id()` helpers |
| `scripts/prepare_data.sh` | Add a top-level step 6 invocation |
| `config/schema_config.yaml` | Add `reaction` + `metabolite` nodes and 4 edges |
| `create_knowledge_graph.py` | Wire `MultiMetabolismAdapter` into the pipeline |
| `scripts/post-import.sh` | Add 7 scalar + 2 full-text indexes; add 4 rollup queries; materialize `Organism_has_metabolite` |
| `scripts/post-import.cypher` | Mirror the same Cypher (kept in sync per CLAUDE.md) |
| `CLAUDE.md` | Document new node/edge labels + step 6 + new fields |

---

## Tasks

### Task 1: Add 5 new KEGG endpoint parsers (tests-first)

**Files:**
- Modify: `multiomics_kg/utils/kegg_utils.py`
- Test: `tests/test_kegg_utils_metabolism.py` (new)

The 5 new endpoints are pure list/link queries that follow the existing `_fetch_text` + `_parse_*` idiom in `kegg_utils.py`. Add them as private parser functions, then in Task 2 wire them into `download_kegg_data()`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_kegg_utils_metabolism.py
"""Unit tests for the 5 new KEGG endpoints introduced in Spec 1.2."""
from __future__ import annotations

import textwrap

from multiomics_kg.utils import kegg_utils


REACTION_LIST_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpyruvate kinase reaction
    rn:R00010\tphosphofructokinase reaction
    invalid_line
    rn:R12345\t
""")

COMPOUND_LIST_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tD-glucose; alpha-D-glucopyranose
    cpd:C00002\tATP; adenosine 5'-triphosphate
    bad_prefix:Cxxxx\tjunk
""")

LINK_CR_FIXTURE = textwrap.dedent("""\
    cpd:C00074\trn:R00200
    cpd:C00008\trn:R00200
    cpd:C00031\trn:R00010
""")

LINK_PR_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpath:rn00010
    rn:R00200\tpath:rn00710
    rn:R00010\tpath:rn00010
""")

LINK_PC_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tpath:map00010
    cpd:C00031\tpath:map00500
""")


def test_parse_reaction_names():
    out = kegg_utils._parse_reaction_names(REACTION_LIST_FIXTURE)
    assert out == {
        "R00200": "pyruvate kinase reaction",
        "R00010": "phosphofructokinase reaction",
        "R12345": "",
    }


def test_parse_compound_names():
    out = kegg_utils._parse_compound_names(COMPOUND_LIST_FIXTURE)
    # Multi-name field: take first synonym (before first '; ')
    assert out["C00031"] == "D-glucose"
    assert out["C00002"] == "ATP"
    assert "Cxxxx" not in out


def test_parse_reaction_to_compounds():
    out = kegg_utils._parse_reaction_to_compounds(LINK_CR_FIXTURE)
    assert sorted(out["R00200"]) == ["C00008", "C00074"]
    assert out["R00010"] == ["C00031"]


def test_parse_reaction_to_pathways_strips_rn_prefix():
    """KEGG `/link/pathway/reaction` returns rn-prefixed pathway IDs.
    We normalize them to the ko-prefixed form used elsewhere in the KG.
    """
    out = kegg_utils._parse_reaction_to_pathways(LINK_PR_FIXTURE)
    assert sorted(out["R00200"]) == ["ko00010", "ko00710"]
    assert out["R00010"] == ["ko00010"]


def test_parse_compound_to_pathways_strips_map_prefix():
    out = kegg_utils._parse_compound_to_pathways(LINK_PC_FIXTURE)
    # map-prefixed pathways → ko-prefixed for consistency with KeggTerm node IDs
    assert sorted(out["C00031"]) == ["ko00010", "ko00500"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_kegg_utils_metabolism.py -v`
Expected: 5 FAIL — `AttributeError: module 'kegg_utils' has no attribute '_parse_reaction_names'` (and the other 4).

- [ ] **Step 3: Add the 5 parsers + 5 URL constants**

Open [multiomics_kg/utils/kegg_utils.py](multiomics_kg/utils/kegg_utils.py) and after the line `_PATHWAY_KO_LIST_URL = ...` (around line 26), add:

```python
_REACTION_LIST_URL = f"{_KEGG_BASE}/list/reaction"
_COMPOUND_LIST_URL = f"{_KEGG_BASE}/list/compound"
_LINK_COMPOUND_REACTION_URL = f"{_KEGG_BASE}/link/compound/reaction"
_LINK_PATHWAY_REACTION_URL = f"{_KEGG_BASE}/link/pathway/reaction"
_LINK_PATHWAY_COMPOUND_URL = f"{_KEGG_BASE}/link/pathway/compound"
```

Then **after** `_parse_ko_to_pathways` (around line 100) and **before** `_parse_brite_hierarchy`, add the 5 parsers:

```python
def _parse_reaction_names(text: str) -> dict[str, str]:
    """Parse `/list/reaction` response into {R#####: name_str}."""
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_id, name = parts
        rxn_id = raw_id.removeprefix("rn:")
        if rxn_id.startswith("R") and rxn_id[1:].isdigit():
            result[rxn_id] = name.strip()
    logger.info(f"Parsed {len(result)} reaction names")
    return result


def _parse_compound_names(text: str) -> dict[str, str]:
    """Parse `/list/compound` into {C#####: first_synonym}.

    KEGG returns semicolon-separated synonyms; we keep the first as canonical.
    """
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_id, names = parts
        cpd_id = raw_id.removeprefix("cpd:")
        if cpd_id.startswith("C") and cpd_id[1:].isdigit():
            first = names.split(";", 1)[0].strip()
            result[cpd_id] = first
    logger.info(f"Parsed {len(result)} compound names")
    return result


def _parse_reaction_to_compounds(text: str) -> dict[str, list[str]]:
    """Parse `/link/compound/reaction` into {R#####: [C#####, ...]}.

    Source line format: `cpd:C00074\trn:R00200`.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_cpd, raw_rxn = parts
        cpd_id = raw_cpd.removeprefix("cpd:")
        rxn_id = raw_rxn.removeprefix("rn:")
        if not (cpd_id.startswith("C") and rxn_id.startswith("R")):
            continue
        result.setdefault(rxn_id, []).append(cpd_id)
    logger.info(f"Parsed compound-reaction links for {len(result)} reactions")
    return result


def _parse_reaction_to_pathways(text: str) -> dict[str, list[str]]:
    """Parse `/link/pathway/reaction` into {R#####: [ko#####, ...]}.

    KEGG returns rn-prefixed pathway IDs (e.g. `path:rn00010`); we rewrite them
    to ko-prefixed form so they match existing KeggTerm pathway node IDs.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_rxn, raw_pw = parts
        rxn_id = raw_rxn.removeprefix("rn:")
        pw_id = raw_pw.removeprefix("path:")
        if not rxn_id.startswith("R"):
            continue
        # Normalize rn00010 → ko00010 (the form used by existing KeggTerm nodes)
        if pw_id.startswith("rn"):
            pw_id = "ko" + pw_id[2:]
        if not pw_id.startswith("ko"):
            continue
        result.setdefault(rxn_id, []).append(pw_id)
    logger.info(f"Parsed reaction-pathway links for {len(result)} reactions")
    return result


def _parse_compound_to_pathways(text: str) -> dict[str, list[str]]:
    """Parse `/link/pathway/compound` into {C#####: [ko#####, ...]}.

    KEGG returns map-prefixed pathway IDs; we rewrite to ko-prefixed for
    consistency with the rest of the KG.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_cpd, raw_pw = parts
        cpd_id = raw_cpd.removeprefix("cpd:")
        pw_id = raw_pw.removeprefix("path:")
        if not cpd_id.startswith("C"):
            continue
        if pw_id.startswith("map"):
            pw_id = "ko" + pw_id[3:]
        if not pw_id.startswith("ko"):
            continue
        result.setdefault(cpd_id, []).append(pw_id)
    logger.info(f"Parsed compound-pathway links for {len(result)} compounds")
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_kegg_utils_metabolism.py -v`
Expected: 5 PASS.

- [ ] **Step 5: Run the full kegg_utils test suite to confirm no regression**

Run: `uv run pytest tests/test_kegg_utils.py tests/test_kegg_utils_metabolism.py -v`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/kegg_utils.py tests/test_kegg_utils_metabolism.py
git commit -m "metabolism: 1.2 — 5 new KEGG endpoint parsers (reaction/compound list + 3 link)"
```

---

### Task 2: Wire 5 new endpoints into `download_kegg_data()`

**Files:**
- Modify: `multiomics_kg/utils/kegg_utils.py` (extend `download_kegg_data`)
- Test: `tests/test_kegg_utils_metabolism.py` (extend)

The 5 new keys must be added to the cached `kegg_data.json` so callers can read them without hitting the network on every build.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_kegg_utils_metabolism.py`:

```python
def test_download_kegg_data_includes_metabolism_keys(tmp_path, monkeypatch):
    """download_kegg_data should populate the 5 new metabolism keys."""

    # Stub the network: return synthetic fixtures for every URL we know about.
    def fake_fetch_text(url: str) -> str:
        if url.endswith("/list/ko"):
            return "ko:K02338\tDNA polymerase III\n"
        if url.endswith("/link/pathway/ko"):
            return "ko:K02338\tpath:ko03030\n"
        if url.endswith("/list/pathway/ko"):
            return "path:ko03030\tDNA replication\n"
        if url.endswith("/list/reaction"):
            return REACTION_LIST_FIXTURE
        if url.endswith("/list/compound"):
            return COMPOUND_LIST_FIXTURE
        if url.endswith("/link/compound/reaction"):
            return LINK_CR_FIXTURE
        if url.endswith("/link/pathway/reaction"):
            return LINK_PR_FIXTURE
        if url.endswith("/link/pathway/compound"):
            return LINK_PC_FIXTURE
        raise AssertionError(f"unexpected URL {url}")

    def fake_fetch_json(url: str) -> dict:
        # Minimal BRITE skeleton — children empty, parsers handle that gracefully
        return {"children": []}

    monkeypatch.setattr(kegg_utils, "_fetch_text", fake_fetch_text)
    monkeypatch.setattr(kegg_utils, "_fetch_json", fake_fetch_json)

    data = kegg_utils.download_kegg_data(tmp_path, force=True)

    # Pre-existing keys still present
    for key in ("ko_names", "pathway_names", "ko_to_pathways"):
        assert key in data

    # New Spec 1.2 keys
    assert data["reaction_names"] == {
        "R00200": "pyruvate kinase reaction",
        "R00010": "phosphofructokinase reaction",
        "R12345": "",
    }
    assert data["compound_names"]["C00031"] == "D-glucose"
    assert sorted(data["reaction_to_compounds"]["R00200"]) == ["C00008", "C00074"]
    assert sorted(data["reaction_to_pathways"]["R00200"]) == ["ko00010", "ko00710"]
    assert sorted(data["compound_to_pathways"]["C00031"]) == ["ko00010", "ko00500"]

    # Cache file written and round-trips
    cache_file = tmp_path / "kegg" / "kegg_data.json"
    assert cache_file.exists()
    import json
    reloaded = json.loads(cache_file.read_text())
    assert reloaded["reaction_names"] == data["reaction_names"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_kegg_utils_metabolism.py::test_download_kegg_data_includes_metabolism_keys -v`
Expected: FAIL — `KeyError: 'reaction_names'` (the function doesn't populate the new keys yet).

- [ ] **Step 3: Extend `download_kegg_data()`**

In [multiomics_kg/utils/kegg_utils.py](multiomics_kg/utils/kegg_utils.py), inside `download_kegg_data` (after `pathway_names = {**api_pw_names, **brite_pw_names}`), insert before the `data = {...}` dict:

```python
    # ── Spec 1.2 metabolism endpoints ─────────────────────────────────
    reaction_names = _parse_reaction_names(_fetch_text(_REACTION_LIST_URL))
    compound_names = _parse_compound_names(_fetch_text(_COMPOUND_LIST_URL))
    reaction_to_compounds = _parse_reaction_to_compounds(_fetch_text(_LINK_COMPOUND_REACTION_URL))
    reaction_to_pathways = _parse_reaction_to_pathways(_fetch_text(_LINK_PATHWAY_REACTION_URL))
    compound_to_pathways = _parse_compound_to_pathways(_fetch_text(_LINK_PATHWAY_COMPOUND_URL))
```

Then extend the `data = {...}` dict with the 5 new keys:

```python
    data = {
        "ko_names": ko_names,
        "pathway_names": pathway_names,
        "ko_to_pathways": ko_to_pathways,
        "pathway_to_subcategory": pathway_to_subcat,
        "subcategory_names": subcat_names,
        "subcategory_to_category": subcat_to_cat,
        "category_names": cat_names,
        # Spec 1.2 — metabolism
        "reaction_names": reaction_names,
        "compound_names": compound_names,
        "reaction_to_compounds": reaction_to_compounds,
        "reaction_to_pathways": reaction_to_pathways,
        "compound_to_pathways": compound_to_pathways,
    }
```

Update the docstring of `download_kegg_data` to mention the 5 new keys.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_kegg_utils_metabolism.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/kegg_utils.py tests/test_kegg_utils_metabolism.py
git commit -m "metabolism: 1.2 — wire 5 new endpoints into download_kegg_data()"
```

---

### Task 3: Drop `resolve_kegg_reaction_to_mnxr` transform

**Files:**
- Modify: `multiomics_kg/download/utils/annotation_transforms.py`
- Modify: `config/gene_annotations_config.yaml`
- Modify: `tests/test_annotation_transforms_metabolism.py`

The Spec 1.2 pivot says `gene["kegg_reactions"]` must contain raw KEGG R-numbers, not MNXR. Remove the transform completely (its tests + YAML wire) so the field passes through unmodified.

- [ ] **Step 1: Write the regression test (raw R-numbers survive)**

In `tests/test_annotation_transforms_metabolism.py`, replace any `_tx_resolve_kegg_reaction_to_mnxr` test with a regression that confirms the transform is **gone**:

```python
def test_resolve_kegg_reaction_transform_removed():
    """Spec 1.2 pivot: KEGG reactions stay as raw R-numbers (no MNX resolution)."""
    from multiomics_kg.download.utils.annotation_transforms import _TRANSFORMS
    assert "resolve_kegg_reaction_to_mnxr" not in _TRANSFORMS
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_annotation_transforms_metabolism.py::test_resolve_kegg_reaction_transform_removed -v`
Expected: FAIL — the registry still contains the key.

- [ ] **Step 3: Remove the transform from the registry**

In [multiomics_kg/download/utils/annotation_transforms.py](multiomics_kg/download/utils/annotation_transforms.py):

1. Remove the `_tx_resolve_kegg_reaction_to_mnxr` function definition (~lines 266-269).
2. Remove the `from multiomics_kg.utils.metabolite_utils import open_resolver, resolve_reaction` import line (~line 251) **only if** `_RESOLVER_CONN` and `_get_resolver_conn` are no longer referenced (they aren't — TCDB and CAZy use their own utils).
3. Remove the `_RESOLVER_CONN`/`_get_resolver_conn` block (~lines 255-263).
4. Remove the `"resolve_kegg_reaction_to_mnxr": _tx_resolve_kegg_reaction_to_mnxr,` entry from `_TRANSFORMS` (~line 298).

**Keep**: `_tx_validate_tcdb`, `_tx_validate_cazy`, and their imports / registry entries.

- [ ] **Step 4: Drop the YAML transform line**

In [config/gene_annotations_config.yaml](config/gene_annotations_config.yaml), under the `kegg_reactions` field (~lines 406-412), remove the `transform: resolve_kegg_reaction_to_mnxr` line. The block becomes:

```yaml
  kegg_reactions:
    type: union
    sources:
      - source: eggnog
        field: KEGG_Reaction
        delimiter: ","
```

- [ ] **Step 5: Run the regression test**

Run: `uv run pytest tests/test_annotation_transforms_metabolism.py -v`
Expected: PASS — registry cleansed; existing TCDB/CAZy tests still pass.

- [ ] **Step 6: Run the wider gene-annotations test (no regression on shape)**

Run: `uv run pytest tests/test_build_gene_annotations_metabolism.py -v`
Expected: any test that asserted MNXR shape on `kegg_reactions` must be updated to assert raw `R*` values. If those tests were absent or only checked TCDB/CAZy, they should still PASS unchanged.

If a test fails because it checked for `MNXR` content, update it to assert raw R-number content (e.g. `assert all(v.startswith("R") for v in gene["kegg_reactions"])`).

- [ ] **Step 7: Run the full test suite to confirm**

Run: `uv run pytest -m "not slow and not kg" -q`
Expected: all PASS.

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/download/utils/annotation_transforms.py \
        config/gene_annotations_config.yaml \
        tests/test_annotation_transforms_metabolism.py \
        tests/test_build_gene_annotations_metabolism.py
git commit -m "metabolism: 1.2 — drop resolve_kegg_reaction_to_mnxr transform; keep raw R-numbers"
```

---

### Task 4: Add 2 indexes to MNX resolver DDL

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py`
- Test: `tests/test_build_metabolite_resolver.py` (extend)

The new helpers in Task 5 query the resolver by the MNX-side ID (e.g. "give me all aliases of MNXM41"). The existing `idx_compound_aliases_value` indexes the *external* side (chebi:17234 etc). We need indexes on the MNX side too.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_build_metabolite_resolver.py`:

```python
def test_resolver_has_mnx_side_indexes(tmp_path):
    """Spec 1.2: resolver DB must have indexes on mnxm_id / mnxr_id columns."""
    import sqlite3
    from multiomics_kg.download import build_metabolite_resolver as bmr

    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(
        "#xref\tmnxm_id\tdescription\n"
        "chebi:17234\tMNXM41\tD-glucose\n"
    )
    reac_xref = tmp_path / "reac_xref.tsv"
    reac_xref.write_text(
        "#xref\tmnxr_id\n"
        "kegg.reaction:R00200\tMNXR101234\n"
    )

    db_path = tmp_path / "resolver.db"
    conn = sqlite3.connect(db_path)
    bmr.build_compound_aliases_table(conn, chem_xref)
    bmr.build_reaction_aliases_table(conn, reac_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='index'")
    indexes = {row[0] for row in cur.fetchall()}

    assert "idx_compound_aliases_mnxm" in indexes
    assert "idx_reaction_aliases_mnxr" in indexes
    conn.close()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_resolver_has_mnx_side_indexes -v`
Expected: FAIL — `assert 'idx_compound_aliases_mnxm' in {...}`.

- [ ] **Step 3: Add the indexes to the DDL**

In [multiomics_kg/download/build_metabolite_resolver.py](multiomics_kg/download/build_metabolite_resolver.py):

In `_ALIAS_DDL` (~line 145), append:
```sql
CREATE INDEX IF NOT EXISTS idx_compound_aliases_mnxm
    ON compound_aliases(mnxm_id);
```

In `_REAC_ALIAS_DDL` (~line 265), append:
```sql
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_mnxr
    ON reaction_aliases(mnxr_id);
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_resolver_has_mnx_side_indexes -v`
Expected: PASS.

- [ ] **Step 5: Run the full resolver test suite to confirm no regression**

Run: `uv run pytest tests/test_build_metabolite_resolver.py -v`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py \
        tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.2 — resolver DDL adds mnxm/mnxr-side indexes"
```

---

### Task 5: Add `mnxm_to_primary_id()` / `mnxr_to_primary_id()` helpers

**Files:**
- Modify: `multiomics_kg/utils/metabolite_utils.py`
- Test: `tests/test_metabolite_utils_primary_id.py` (new)

These helpers map an MNX ID to the canonical KG primary ID (KEGG > ChEBI > MNX for compounds; KEGG > Rhea > MNX for reactions). Used by step 6 to assign primary IDs deterministically when a metabolite isn't in KEGG.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_metabolite_utils_primary_id.py
"""Unit tests for mnxm_to_primary_id() / mnxr_to_primary_id()."""
from __future__ import annotations

import sqlite3

from multiomics_kg.utils import metabolite_utils as mu


def _make_db(tmp_path):
    """Create a minimal in-memory resolver DB with synthetic aliases."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        -- MNXM41 has both KEGG and ChEBI: KEGG should win
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');

        -- MNXM999 has only ChEBI: chebi should win
        INSERT INTO compound_aliases VALUES ('chebi', '12345', 'MNXM999');

        -- MNXM_orphan has nothing: fallback to mnx
        -- (no rows inserted)

        -- MNXR101234 has both KEGG and Rhea: KEGG should win
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
        INSERT INTO reaction_aliases VALUES ('rhea', '10828', 'MNXR101234');

        -- MNXR_rhea_only has only Rhea
        INSERT INTO reaction_aliases VALUES ('rhea', '99999', 'MNXR_rhea_only');
    """)
    conn.commit()
    return conn


def test_mnxm_kegg_wins(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM41", conn) == "kegg.compound:C00031"


def test_mnxm_chebi_when_no_kegg(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM999", conn) == "chebi:12345"


def test_mnxm_orphan_falls_back_to_mnx(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM_orphan", conn) == "mnx:MNXM_orphan"


def test_mnxr_kegg_wins(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR101234", conn) == "kegg.reaction:R00200"


def test_mnxr_rhea_when_no_kegg(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR_rhea_only", conn) == "rhea:99999"


def test_mnxr_orphan_falls_back_to_mnx(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR_unknown", conn) == "mnx:MNXR_unknown"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_metabolite_utils_primary_id.py -v`
Expected: 6 FAIL — `AttributeError: module 'metabolite_utils' has no attribute 'mnxm_to_primary_id'`.

- [ ] **Step 3: Implement the helpers**

Append to [multiomics_kg/utils/metabolite_utils.py](multiomics_kg/utils/metabolite_utils.py):

```python
# ── Primary-ID resolution (KG-side canonical IDs) ────────────────────────────
#
# Per Spec 1.2: KEGG-native primary IDs win, with bioregistry-style fallbacks.
# Compound: kegg.compound > chebi > mnx
# Reaction: kegg.reaction > rhea > mnx

_COMPOUND_PRIMARY_PRIORITY = ("kegg.compound", "chebi")
_REACTION_PRIMARY_PRIORITY = ("kegg.reaction", "rhea")


def mnxm_to_primary_id(mnxm_id: str, conn: sqlite3.Connection) -> str:
    """Map an MNXM* ID to its canonical KG primary ID.

    Priority: kegg.compound > chebi > fallback `mnx:<MNXM>`.
    """
    cur = conn.cursor()
    for source in _COMPOUND_PRIMARY_PRIORITY:
        cur.execute(
            "SELECT value FROM compound_aliases WHERE mnxm_id = ? AND source = ? ORDER BY value LIMIT 1",
            (mnxm_id, source),
        )
        row = cur.fetchone()
        if row:
            return f"{source}:{row[0]}"
    return f"mnx:{mnxm_id}"


def mnxr_to_primary_id(mnxr_id: str, conn: sqlite3.Connection) -> str:
    """Map an MNXR* ID to its canonical KG primary ID.

    Priority: kegg.reaction > rhea > fallback `mnx:<MNXR>`.
    """
    cur = conn.cursor()
    for source in _REACTION_PRIMARY_PRIORITY:
        cur.execute(
            "SELECT value FROM reaction_aliases WHERE mnxr_id = ? AND source = ? ORDER BY value LIMIT 1",
            (mnxr_id, source),
        )
        row = cur.fetchone()
        if row:
            return f"{source}:{row[0]}"
    return f"mnx:{mnxr_id}"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_metabolite_utils_primary_id.py -v`
Expected: all 6 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/metabolite_utils.py \
        tests/test_metabolite_utils_primary_id.py
git commit -m "metabolism: 1.2 — mnxm_to_primary_id / mnxr_to_primary_id helpers"
```

---

### Task 6: `build_kegg_metabolism_xrefs.py` — pruning + structure

**Files:**
- Create: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py` (new)

The script walks every strain's `gene_annotations_merged.json`, collects raw KEGG R-numbers from `gene["kegg_reactions"]`, then expands to the C-numbers that participate in those reactions via `kegg_data["reaction_to_compounds"]`. This is the *prune* part (analogous to step 5).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_kegg_metabolism_xrefs.py
"""Unit tests for step 6: build_kegg_metabolism_xrefs."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def test_collect_gene_reachable_kegg_ids(tmp_path, monkeypatch):
    # Two synthetic strains
    s1 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MED4"
    s2 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MIT9301"
    _write_strain_annotations(s1, {
        "PMM0001": {"kegg_reactions": ["R00200", "R00010"]},
        "PMM0002": {"kegg_reactions": ["R00010"]},
        "PMM0003": {},  # no reactions
    })
    _write_strain_annotations(s2, {
        "P9301_001": {"kegg_reactions": ["R12345"]},
    })

    # Stub: load_genome_rows returns rows pointing at the two strains
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [
        {"data_dir": str(s1)},
        {"data_dir": str(s2)},
    ])

    # Stub kegg_data.json with a tiny reaction_to_compounds map
    kegg_data = {
        "reaction_names": {
            "R00200": "pyruvate kinase reaction",
            "R00010": "phosphofructokinase",
            "R12345": "test rxn",
        },
        "compound_names": {
            "C00074": "phosphoenolpyruvate",
            "C00008": "ADP",
            "C00031": "D-glucose",
            "C00002": "ATP",
            "C99999": "test cpd",
        },
        "reaction_to_compounds": {
            "R00200": ["C00074", "C00008"],
            "R00010": ["C00031", "C00002"],
            "R12345": ["C99999"],
            "R_unused": ["C_should_not_appear"],
        },
        "reaction_to_pathways": {
            "R00200": ["ko00010"],
        },
    }

    rxns, cpds = bx._collect_gene_reachable_ids(kegg_data)

    # Only gene-reachable R-numbers
    assert rxns == {"R00200", "R00010", "R12345"}
    # C-numbers reachable from those reactions
    assert cpds == {"C00074", "C00008", "C00031", "C00002", "C99999"}
    # R_unused / C_should_not_appear are pruned
    assert "R_unused" not in rxns
    assert "C_should_not_appear" not in cpds


def test_collect_handles_missing_kegg_reactions_field(tmp_path, monkeypatch):
    s = tmp_path / "strain"
    _write_strain_annotations(s, {
        "g1": {},  # no kegg_reactions key
        "g2": {"kegg_reactions": []},  # empty list
        "g3": {"kegg_reactions": ["R00001"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(s)}])
    kegg_data = {"reaction_to_compounds": {"R00001": ["C00001"]}}
    rxns, cpds = bx._collect_gene_reachable_ids(kegg_data)
    assert rxns == {"R00001"}
    assert cpds == {"C00001"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.download.build_kegg_metabolism_xrefs'`.

- [ ] **Step 3: Create the module skeleton + pruning logic**

```python
# multiomics_kg/download/build_kegg_metabolism_xrefs.py
"""Step 6 — prune-then-enrich KEGG metabolism cache.

Walks every strain's gene_annotations_merged.json to identify gene-reachable
KEGG R-numbers, expands to the C-numbers that participate in those reactions,
and enriches each with KEGG metadata (name, EC, pathways) + MNX cross-refs
(MNXR/MNXM, ChEBI, HMDB, InChIKey, formula, mass).

Output: cache/data/kegg/kegg_metabolism_xrefs.json (~1 MB).

This file is read by metabolism_adapter.py at KG build time. The 2.6 GB MNX
SQLite resolver is opened only here (prepare-data time), not at build time.

Usage:
    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]
"""
from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations
from multiomics_kg.utils import metabolite_utils as mu

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
OUTPUT_FILE = KEGG_CACHE_DIR / "kegg_metabolism_xrefs.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"

log = logging.getLogger(__name__)


def _collect_gene_reachable_ids(kegg_data: dict) -> tuple[set[str], set[str]]:
    """Walk all strains' gene_annotations_merged.json and KEGG link data.

    Returns (R_numbers, C_numbers): the gene-reachable subsets.
    """
    reachable_rxns: set[str] = set()
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for gene in genes.values():
            for rxn in gene.get("kegg_reactions", []) or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    reachable_rxns.add(rxn)

    rxn_to_cpds = kegg_data.get("reaction_to_compounds", {})
    reachable_cpds: set[str] = set()
    for rxn in reachable_rxns:
        for cpd in rxn_to_cpds.get(rxn, []):
            reachable_cpds.add(cpd)

    log.info(f"Pruned to {len(reachable_rxns)} reactions, {len(reachable_cpds)} compounds")
    return reachable_rxns, reachable_cpds
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v`
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        tests/test_build_kegg_metabolism_xrefs.py
git commit -m "metabolism: 1.2 — step 6 module + gene-reachable pruning"
```

---

### Task 7: Step 6 — enrichment + main()

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (add enrichment + main)
- Test: `tests/test_build_kegg_metabolism_xrefs.py` (extend)

For each gene-reachable R-number: pull KEGG name + pathway IDs from `kegg_data`, then for `mnxr_id`/`rhea_ids`/`reaction_class`/`mass_balance` query the resolver. For each C-number: same pattern with `mnxm_id`, `chebi_id`, `hmdb_id`, `inchikey`, `formula`, `mass`, `smiles`.

EC numbers come from KEGG `/get/rn:` *(out of scope here — defer to MNX `reactions.classifs` from the resolver)*.

- [ ] **Step 1: Write failing test for enrichment**

Append to `tests/test_build_kegg_metabolism_xrefs.py`:

```python
def test_enrich_reaction(tmp_path):
    """Reaction enrichment: KEGG name + pathway + MNX cross-refs."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db) if (lambda: db)() else None  # noqa
    import sqlite3
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
            reference TEXT, classifs TEXT, is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        INSERT INTO reactions VALUES
            ('MNXR101234', '1 C00074 + 1 C00008 = 1 C00022 + 1 C00002',
             'kegg.reaction:R00200', '2.7.1.40', 'B', NULL);
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
        INSERT INTO reaction_aliases VALUES ('rhea', '10828', 'MNXR101234');

        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'InChI=1/C6H12O6',
             'WQZGKKKJIJFFOK-GASJEMHNSA-N', 'OC[CH]1OC(O)[CH](O)[CH](O)[CH]1O',
             0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('hmdb', 'HMDB0000122', 'MNXM41');
    """)
    conn.commit()

    kegg_data = {
        "reaction_names": {"R00200": "pyruvate kinase reaction"},
        "compound_names": {"C00031": "D-glucose"},
        "reaction_to_compounds": {"R00200": ["C00074", "C00008"]},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
    }

    rxn = bx._enrich_reaction("R00200", kegg_data, conn)
    assert rxn["name"] == "pyruvate kinase reaction"
    assert rxn["compound_ids"] == ["C00074", "C00008"]
    assert rxn["kegg_pathway_ids"] == ["ko00010", "ko00710"]
    assert rxn["mnxr_id"] == "MNXR101234"
    assert rxn["rhea_ids"] == ["10828"]
    assert rxn["mass_balance"] == "balanced"
    assert rxn["reaction_class"] == "chemical"
    assert rxn["ec_numbers"] == ["2.7.1.40"]


def test_enrich_compound(tmp_path):
    import sqlite3
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'InChI=...',
             'WQZGKKKJIJFFOK-GASJEMHNSA-N', 'OC[CH]1...', 0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('hmdb', 'HMDB0000122', 'MNXM41');
    """)
    conn.commit()

    kegg_data = {"compound_names": {"C00031": "D-glucose"}}
    cpd = bx._enrich_compound("C00031", kegg_data, conn)
    assert cpd["name"] == "D-glucose"
    assert cpd["formula"] == "C6H12O6"
    assert cpd["mass"] == 180.063
    assert cpd["inchikey"] == "WQZGKKKJIJFFOK-GASJEMHNSA-N"
    assert cpd["mnxm_id"] == "MNXM41"
    assert cpd["chebi_id"] == "17234"
    assert cpd["hmdb_id"] == "HMDB0000122"


def test_enrich_compound_orphan_no_mnx_match(tmp_path):
    """C-number with no MNX entry: enrichment returns just KEGG name."""
    import sqlite3
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
    """)
    conn.commit()

    kegg_data = {"compound_names": {"C99999": "obscure compound"}}
    cpd = bx._enrich_compound("C99999", kegg_data, conn)
    assert cpd["name"] == "obscure compound"
    assert cpd.get("mnxm_id") is None
    assert cpd.get("formula") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v`
Expected: 3 new FAIL — `AttributeError: module ... has no attribute '_enrich_reaction'`.

- [ ] **Step 3: Add enrichment + main**

Append to [multiomics_kg/download/build_kegg_metabolism_xrefs.py](multiomics_kg/download/build_kegg_metabolism_xrefs.py):

```python
def _resolve_mnx_for_kegg_reaction(kegg_id: str, conn: sqlite3.Connection) -> str | None:
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxr_id FROM reaction_aliases WHERE source = 'kegg.reaction' AND value = ? LIMIT 1",
        (kegg_id,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _resolve_mnx_for_kegg_compound(kegg_id: str, conn: sqlite3.Connection) -> str | None:
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxm_id FROM compound_aliases WHERE source = 'kegg.compound' AND value = ? LIMIT 1",
        (kegg_id,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _query_aliases(table: str, mnx_col: str, mnx_id: str, source: str,
                   conn: sqlite3.Connection) -> list[str]:
    cur = conn.cursor()
    cur.execute(
        f"SELECT value FROM {table} WHERE {mnx_col} = ? AND source = ? ORDER BY value",
        (mnx_id, source),
    )
    return [r[0] for r in cur.fetchall()]


def _enrich_reaction(rxn_id: str, kegg_data: dict, conn: sqlite3.Connection) -> dict:
    """Build the per-reaction enriched dict written to xrefs JSON."""
    out: dict = {
        "name": kegg_data.get("reaction_names", {}).get(rxn_id, ""),
        "compound_ids": kegg_data.get("reaction_to_compounds", {}).get(rxn_id, []),
        "kegg_pathway_ids": kegg_data.get("reaction_to_pathways", {}).get(rxn_id, []),
        "ec_numbers": [],
        "mnxr_id": None,
        "rhea_ids": [],
        "mass_balance": "unbalanced",
        "reaction_class": "chemical",
    }
    mnxr_id = _resolve_mnx_for_kegg_reaction(rxn_id, conn)
    if mnxr_id is None:
        return out
    out["mnxr_id"] = mnxr_id

    out["rhea_ids"] = _query_aliases("reaction_aliases", "mnxr_id", mnxr_id, "rhea", conn)

    cur = conn.cursor()
    cur.execute(
        "SELECT classifs, is_balanced, is_transport FROM reactions WHERE mnxr_id = ?",
        (mnxr_id,),
    )
    row = cur.fetchone()
    if row:
        classifs, is_balanced, is_transport = row
        if classifs:
            out["ec_numbers"] = [c.strip() for c in classifs.split(";") if c.strip()]
        out["mass_balance"] = "balanced" if (is_balanced or "").upper() == "B" else "unbalanced"
        out["reaction_class"] = "transport" if (is_transport or "").upper() == "T" else "chemical"
    return out


def _enrich_compound(cpd_id: str, kegg_data: dict, conn: sqlite3.Connection) -> dict:
    """Build the per-compound enriched dict written to xrefs JSON."""
    out: dict = {
        "name": kegg_data.get("compound_names", {}).get(cpd_id, ""),
        "mnxm_id": None,
        "formula": None,
        "mass": None,
        "inchikey": None,
        "smiles": None,
        "chebi_id": None,
        "hmdb_id": None,
    }
    mnxm_id = _resolve_mnx_for_kegg_compound(cpd_id, conn)
    if mnxm_id is None:
        return out
    out["mnxm_id"] = mnxm_id

    cur = conn.cursor()
    cur.execute(
        "SELECT formula, mass, inchikey, smiles FROM compounds WHERE mnxm_id = ?",
        (mnxm_id,),
    )
    row = cur.fetchone()
    if row:
        out["formula"] = row[0] or None
        out["mass"] = row[1] if row[1] is not None else None
        out["inchikey"] = row[2] or None
        out["smiles"] = row[3] or None

    chebi = _query_aliases("compound_aliases", "mnxm_id", mnxm_id, "chebi", conn)
    if chebi:
        out["chebi_id"] = chebi[0]
    hmdb = _query_aliases("compound_aliases", "mnxm_id", mnxm_id, "hmdb", conn)
    if hmdb:
        out["hmdb_id"] = hmdb[0]
    return out


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if OUTPUT_FILE.exists() and not force:
        log.info(f"{OUTPUT_FILE} exists; use --force to rebuild.")
        return

    if not KEGG_DATA_FILE.exists():
        raise FileNotFoundError(
            f"{KEGG_DATA_FILE} missing — run KEGG download first (kegg_anno_adapter)."
        )
    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run prepare_data.sh step 0 sub-step 7 first."
        )

    log.info("Loading KEGG data ...")
    kegg_data = json.loads(KEGG_DATA_FILE.read_text())

    log.info("Pruning to gene-reachable subset ...")
    rxn_ids, cpd_ids = _collect_gene_reachable_ids(kegg_data)

    log.info("Opening MNX resolver ...")
    conn = mu.open_resolver(RESOLVER_DB)

    log.info(f"Enriching {len(rxn_ids)} reactions ...")
    reactions = {rxn_id: _enrich_reaction(rxn_id, kegg_data, conn) for rxn_id in sorted(rxn_ids)}

    log.info(f"Enriching {len(cpd_ids)} compounds ...")
    compounds = {cpd_id: _enrich_compound(cpd_id, kegg_data, conn) for cpd_id in sorted(cpd_ids)}

    conn.close()

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps({"reactions": reactions, "compounds": compounds}, indent=2))
    log.info(f"Wrote {OUTPUT_FILE} ({OUTPUT_FILE.stat().st_size:,} bytes, "
             f"{len(reactions)} reactions, {len(compounds)} compounds).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild the output cache even if it exists.")
    args = parser.parse_args()
    main(force=args.force)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v`
Expected: 5 PASS (2 from Task 6 + 3 from Task 7).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        tests/test_build_kegg_metabolism_xrefs.py
git commit -m "metabolism: 1.2 — step 6 enrichment (KEGG name + MNX xrefs) + main()"
```

---

### Task 8: Wire step 6 into `prepare_data.sh`

**Files:**
- Modify: `scripts/prepare_data.sh`

The new step 6 sits between step 5 (`build_og_descriptions`) and the (currently-absent) end of pipeline. Follow the same `with_step_log` shell pattern as steps 1-5.

- [ ] **Step 1: Read the existing pattern**

Read [scripts/prepare_data.sh](scripts/prepare_data.sh) lines 165-185 to see how step 5 is invoked. The step 6 block mirrors it.

- [ ] **Step 2: Add step 6 in the dispatch loop**

In `scripts/prepare_data.sh`, after the step 5 block, add:

```bash
            6)
                echo "=== Step 6: build KEGG metabolism xrefs ==="
                with_step_log \
                    "$LOG_DIR/prepare_data_step6.log" \
                    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs \
                        $FORCE_FLAG
                ;;
```

Update the `STEPS_DEFAULT` (or whichever variable lists the default step set) to include `6`. Update the comment block at the top of the file (the "Logs:" line at ~line 31) from "step5" → "step6".

If the script gates on a `--steps` arg, also update the help text and the validation regex (e.g. `^[0-6]$` instead of `^[0-5]$`).

- [ ] **Step 3: Smoke test the dispatch (no actual run)**

Run: `bash scripts/prepare_data.sh --steps 6 --dry-run` if dry-run is supported. Otherwise just confirm the bash script parses with `bash -n scripts/prepare_data.sh`:

```bash
bash -n scripts/prepare_data.sh
```
Expected: no output, exit 0.

- [ ] **Step 4: Commit**

```bash
git add scripts/prepare_data.sh
git commit -m "metabolism: 1.2 — wire step 6 (build_kegg_metabolism_xrefs) into prepare_data.sh"
```

---

### Task 9: `schema_config.yaml` — 2 new nodes + 4 new edges

**Files:**
- Modify: `config/schema_config.yaml`

Mirror the OrthologGroup / Pfam patterns. **No EnvironmentalCondition or other deprecated forms — purely additive.**

- [ ] **Step 1: Add the `reaction` node block**

After the existing `pfam clan:` block in [config/schema_config.yaml](config/schema_config.yaml) (or in a logical location near other functional-annotation nodes), add:

```yaml
reaction:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.reaction
  label_in_input: reaction
  properties:
    kegg_reaction_id: str        # raw KEGG R-number, e.g. "R00200"
    name: str                    # KEGG reaction name
    ec_numbers: str[]            # e.g. ["2.7.1.40"]
    kegg_pathway_ids: str[]      # ko-prefixed pathway IDs, e.g. ["ko00010"]
    mnxr_id: str                 # nullable; MNX cross-reference
    rhea_ids: str[]
    reaction_class: str          # "transport" | "chemical"
    mass_balance: str            # "balanced" | "unbalanced"
    organisms: str[]             # post-import: distinct organism names with a gene catalyzing
    organism_count: int          # post-import
    gene_count: int              # post-import: distinct genes catalyzing
```

- [ ] **Step 2: Add the `metabolite` node block**

```yaml
metabolite:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.compound
  label_in_input: metabolite
  properties:
    kegg_compound_id: str        # raw KEGG C-number; nullable for ChEBI/MNX-only
    name: str
    formula: str
    mass: float
    inchikey: str
    smiles: str
    mnxm_id: str                 # nullable
    chebi_id: str                # nullable
    hmdb_id: str                 # nullable
    organism_count: int          # post-import
    gene_count: int              # post-import
```

- [ ] **Step 3: Add the 4 edges**

Place near the existing Gene-association edges:

```yaml
gene to reaction association:
  is_a: ReactionToParticipantAssociation
  represented_as: edge
  label_as_edge: gene_catalyzes_reaction
  source: gene
  target: reaction
  label_in_input: gene_catalyzes_reaction

reaction to metabolite association:
  is_a: association
  represented_as: edge
  label_as_edge: reaction_has_metabolite
  source: reaction
  target: metabolite
  label_in_input: reaction_has_metabolite

reaction to kegg pathway association:
  is_a: association
  represented_as: edge
  label_as_edge: reaction_in_kegg_pathway
  source: reaction
  target: kegg term
  label_in_input: reaction_in_kegg_pathway

organism to metabolite association:
  is_a: association
  represented_as: edge
  label_as_edge: organism_has_metabolite
  source: organism taxon
  target: metabolite
  label_in_input: organism_has_metabolite
```

- [ ] **Step 4: Sanity-check the YAML parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output, exit 0.

- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml
git commit -m "metabolism: 1.2 — schema additions (Reaction + Metabolite + 4 edges)"
```

---

### Task 10: `metabolism_adapter.py` — `get_nodes()`

**Files:**
- Create: `multiomics_kg/adapters/metabolism_adapter.py`
- Test: `tests/test_metabolism_adapter.py` (new)

The adapter is a flat (non-strain-scoped) reader of the xrefs JSON. The Multi wrapper isn't strictly needed for nodes (one source file), but we provide one for consistency with the other adapters and to keep `get_edges()` strain-scoped.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_metabolism_adapter.py
"""Unit tests for MetabolismAdapter and MultiMetabolismAdapter."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.adapters import metabolism_adapter as ma


XREFS_FIXTURE = {
    "reactions": {
        "R00200": {
            "name": "pyruvate kinase reaction",
            "compound_ids": ["C00074", "C00008"],
            "kegg_pathway_ids": ["ko00010", "ko00710"],
            "ec_numbers": ["2.7.1.40"],
            "mnxr_id": "MNXR101234",
            "rhea_ids": ["10828"],
            "mass_balance": "balanced",
            "reaction_class": "chemical",
        },
        "R12345": {
            "name": "test reaction",
            "compound_ids": ["C99999"],
            "kegg_pathway_ids": [],
            "ec_numbers": [],
            "mnxr_id": None,
            "rhea_ids": [],
            "mass_balance": "unbalanced",
            "reaction_class": "chemical",
        },
    },
    "compounds": {
        "C00074": {
            "name": "phosphoenolpyruvate",
            "formula": "C3H5O6P",
            "mass": 168.04,
            "inchikey": "DTBNBXWJWCWCIK-UHFFFAOYSA-N",
            "smiles": "OC(=O)C(=C)OP(O)(O)=O",
            "mnxm_id": "MNXM73",
            "chebi_id": "44897",
            "hmdb_id": "HMDB0000263",
        },
        "C00008": {
            "name": "ADP",
            "formula": "C10H15N5O10P2",
            "mass": 427.20,
            "inchikey": "XTWYTFMLZFPYCI-KQYNXXCUSA-N",
            "smiles": "Nc1ncnc2c1ncn2[CH]1O[CH]...",
            "mnxm_id": "MNXM7",
            "chebi_id": "16761",
            "hmdb_id": "HMDB0001341",
        },
        "C99999": {
            "name": "obscure compound",
            "formula": None,
            "mass": None,
            "inchikey": None,
            "smiles": None,
            "mnxm_id": None,
            "chebi_id": None,
            "hmdb_id": None,
        },
    },
}


@pytest.fixture
def xrefs_file(tmp_path):
    p = tmp_path / "kegg_metabolism_xrefs.json"
    p.write_text(json.dumps(XREFS_FIXTURE))
    return p


def test_reaction_nodes_have_kegg_primary_id(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    nodes = list(adapter.get_nodes())

    rxn_ids = [n[0] for n in nodes if n[1] == "reaction"]
    assert "kegg.reaction:R00200" in rxn_ids
    assert "kegg.reaction:R12345" in rxn_ids


def test_reaction_node_properties(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    rxn = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.reaction:R00200" and n[1] == "reaction")
    _, _, props = rxn
    assert props["kegg_reaction_id"] == "R00200"
    assert props["name"] == "pyruvate kinase reaction"
    assert props["ec_numbers"] == ["2.7.1.40"]
    assert props["kegg_pathway_ids"] == ["ko00010", "ko00710"]
    assert props["mnxr_id"] == "MNXR101234"
    assert props["rhea_ids"] == ["10828"]
    assert props["mass_balance"] == "balanced"
    assert props["reaction_class"] == "chemical"


def test_metabolite_nodes_have_kegg_primary_id(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    nodes = list(adapter.get_nodes())
    cpd_ids = [n[0] for n in nodes if n[1] == "metabolite"]
    assert "kegg.compound:C00074" in cpd_ids
    assert "kegg.compound:C99999" in cpd_ids


def test_metabolite_node_properties_drop_null_fields(xrefs_file):
    """Sparse properties: nulls should be omitted from the property dict."""
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    cpd = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.compound:C99999" and n[1] == "metabolite")
    _, _, props = cpd
    assert props["kegg_compound_id"] == "C99999"
    assert props["name"] == "obscure compound"
    # All other fields are null in the fixture; should be absent from the dict
    assert "formula" not in props
    assert "mnxm_id" not in props


def test_metabolite_full_properties_preserved(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    cpd = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.compound:C00074" and n[1] == "metabolite")
    _, _, props = cpd
    assert props["formula"] == "C3H5O6P"
    assert props["mass"] == 168.04
    assert props["mnxm_id"] == "MNXM73"
    assert props["chebi_id"] == "44897"
    assert props["hmdb_id"] == "HMDB0000263"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_metabolism_adapter.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.adapters.metabolism_adapter'`.

- [ ] **Step 3: Create the adapter**

```python
# multiomics_kg/adapters/metabolism_adapter.py
"""Metabolism adapter — emits Reaction + Metabolite nodes and 3 edge types.

Pure file reader: consumes
- cache/data/kegg/kegg_metabolism_xrefs.json (built by step 6)
- per-strain cache/data/<organism>/genomes/<strain>/gene_annotations_merged.json

No SQLite / no KEGG REST at build time.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

log = logging.getLogger(__name__)


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_XREFS = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_metabolism_xrefs.json"


def _clean_str(value):
    """Remove characters that break neo4j-admin import (single quote, pipe)."""
    if not isinstance(value, str):
        return value
    return value.replace("'", "^").replace("|", "")


def _drop_nulls(props: dict) -> dict:
    """Sparse-output convention: drop keys whose value is None or empty list."""
    return {k: v for k, v in props.items() if v is not None and v != []}


class MetabolismAdapter:
    """Single-pass adapter emitting Reaction + Metabolite nodes from the xrefs cache."""

    def __init__(self, xrefs_path: Path | str | None = None, test_mode: bool = False):
        self.xrefs_path = Path(xrefs_path) if xrefs_path else DEFAULT_XREFS
        self.test_mode = test_mode
        self._xrefs: dict | None = None

    def _load(self) -> dict:
        if self._xrefs is None:
            if not self.xrefs_path.exists():
                log.warning(f"{self.xrefs_path} missing; emitting no metabolism nodes")
                self._xrefs = {"reactions": {}, "compounds": {}}
            else:
                self._xrefs = json.loads(self.xrefs_path.read_text())
        return self._xrefs

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        data = self._load()
        n = 0
        for rxn_id, rxn in data.get("reactions", {}).items():
            if self.test_mode and n >= 100:
                break
            node_id = f"kegg.reaction:{rxn_id}"
            props = _drop_nulls({
                "kegg_reaction_id": rxn_id,
                "name": _clean_str(rxn.get("name", "")),
                "ec_numbers": list(rxn.get("ec_numbers", [])),
                "kegg_pathway_ids": list(rxn.get("kegg_pathway_ids", [])),
                "mnxr_id": rxn.get("mnxr_id"),
                "rhea_ids": list(rxn.get("rhea_ids", [])),
                "mass_balance": rxn.get("mass_balance"),
                "reaction_class": rxn.get("reaction_class"),
            })
            yield node_id, "reaction", props
            n += 1

        n = 0
        for cpd_id, cpd in data.get("compounds", {}).items():
            if self.test_mode and n >= 100:
                break
            node_id = f"kegg.compound:{cpd_id}"
            props = _drop_nulls({
                "kegg_compound_id": cpd_id,
                "name": _clean_str(cpd.get("name", "")),
                "formula": _clean_str(cpd.get("formula")),
                "mass": cpd.get("mass"),
                "inchikey": _clean_str(cpd.get("inchikey")),
                "smiles": _clean_str(cpd.get("smiles")),
                "mnxm_id": cpd.get("mnxm_id"),
                "chebi_id": cpd.get("chebi_id"),
                "hmdb_id": cpd.get("hmdb_id"),
            })
            yield node_id, "metabolite", props
            n += 1


class MultiMetabolismAdapter(MetabolismAdapter):
    """Multi-strain wrapper. Nodes come from the global xrefs cache; edges
    come from per-strain gene_annotations_merged.json (see Task 11)."""

    def __init__(self, genome_config_file: str, xrefs_path: Path | str | None = None,
                 test_mode: bool = False):
        super().__init__(xrefs_path=xrefs_path, test_mode=test_mode)
        self.genome_config_file = genome_config_file
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_metabolism_adapter.py -v`
Expected: 5 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolism_adapter.py \
        tests/test_metabolism_adapter.py
git commit -m "metabolism: 1.2 — MetabolismAdapter + MultiMetabolismAdapter (nodes)"
```

---

### Task 11: `metabolism_adapter.py` — `get_edges()`

**Files:**
- Modify: `multiomics_kg/adapters/metabolism_adapter.py` (add `get_edges`)
- Test: `tests/test_metabolism_adapter.py` (extend)

Edges:
- `Gene_catalyzes_reaction`: per-strain, derived from `gene["kegg_reactions"]`
- `Reaction_has_metabolite`: from xrefs `reactions[r].compound_ids`
- `Reaction_in_kegg_pathway`: from xrefs `reactions[r].kegg_pathway_ids`

Edge IDs follow `<source>_<predicate>_<target>` so reruns are idempotent.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_metabolism_adapter.py`:

```python
def test_reaction_metabolite_edges(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    edges = [e for e in adapter._reaction_metabolite_edges()]
    pairs = {(e[1], e[2]) for e in edges}
    assert ("kegg.reaction:R00200", "kegg.compound:C00074") in pairs
    assert ("kegg.reaction:R00200", "kegg.compound:C00008") in pairs
    assert ("kegg.reaction:R12345", "kegg.compound:C99999") in pairs
    # All edges have label "reaction_has_metabolite"
    assert {e[3] for e in edges} == {"reaction_has_metabolite"}


def test_reaction_pathway_edges_target_kegg_term(xrefs_file):
    """Pathway target uses the existing KeggTerm primary-ID convention."""
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    edges = [e for e in adapter._reaction_pathway_edges()]
    targets = {e[2] for e in edges if e[1] == "kegg.reaction:R00200"}
    # KeggTerm IDs use the format "kegg.pathway:ko00010" in the existing graph.
    # Confirm by checking the kegg_anno_adapter or schema (do NOT use kegg.reaction prefix here).
    assert "kegg.pathway:ko00010" in targets
    assert "kegg.pathway:ko00710" in targets


def test_gene_reaction_edges(tmp_path, xrefs_file, monkeypatch):
    # Synthetic strain with 2 genes, both catalyzing R00200
    strain_dir = tmp_path / "strain"
    strain_dir.mkdir()
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM0001": {"kegg_reactions": ["R00200"]},
        "PMM0002": {"kegg_reactions": ["R00200", "R12345"]},
        "PMM0003": {"kegg_reactions": []},
    }))
    monkeypatch.setattr(ma, "load_genome_rows",
                        lambda: [{"data_dir": str(strain_dir)}])

    adapter = ma.MultiMetabolismAdapter(
        genome_config_file="ignored", xrefs_path=xrefs_file,
    )
    edges = [e for e in adapter._gene_reaction_edges()]
    sources = {(e[1], e[2]) for e in edges}
    assert ("ncbigene:PMM0001", "kegg.reaction:R00200") in sources
    assert ("ncbigene:PMM0002", "kegg.reaction:R00200") in sources
    assert ("ncbigene:PMM0002", "kegg.reaction:R12345") in sources
    # PMM0003 has no kegg_reactions → no edges
    assert not any(e[1] == "ncbigene:PMM0003" for e in edges)


def test_gene_reaction_edges_skip_unknown_reactions(tmp_path, xrefs_file, monkeypatch):
    """A gene's R-number not present in the xrefs JSON must not produce a dangling edge."""
    strain_dir = tmp_path / "strain"
    strain_dir.mkdir()
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_dangling": {"kegg_reactions": ["R_not_in_xrefs"]},
    }))
    monkeypatch.setattr(ma, "load_genome_rows",
                        lambda: [{"data_dir": str(strain_dir)}])
    adapter = ma.MultiMetabolismAdapter(
        genome_config_file="ignored", xrefs_path=xrefs_file,
    )
    edges = [e for e in adapter._gene_reaction_edges()]
    assert not any(e[2] == "kegg.reaction:R_not_in_xrefs" for e in edges)


def test_get_edges_yields_all_three_types(tmp_path, xrefs_file, monkeypatch):
    strain_dir = tmp_path / "strain"
    strain_dir.mkdir()
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM0001": {"kegg_reactions": ["R00200"]},
    }))
    monkeypatch.setattr(ma, "load_genome_rows",
                        lambda: [{"data_dir": str(strain_dir)}])

    adapter = ma.MultiMetabolismAdapter(
        genome_config_file="ignored", xrefs_path=xrefs_file,
    )
    labels = {e[3] for e in adapter.get_edges()}
    assert labels == {
        "gene_catalyzes_reaction",
        "reaction_has_metabolite",
        "reaction_in_kegg_pathway",
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_metabolism_adapter.py -v`
Expected: 5 new FAIL.

- [ ] **Step 3: Verify the existing KeggTerm pathway primary-ID format**

Before implementing, confirm the pathway primary-ID convention used by the existing kegg_annotation_adapter. The `Reaction_in_kegg_pathway` edge target must match what `MultiKeggAnnotationAdapter` emits.

Run: `grep -nE "kegg.pathway|pathway:" multiomics_kg/adapters/kegg_annotation_adapter.py | head -10`

If the format is `kegg.pathway:<koXXXXX>`, use that. If it differs (e.g. just `<koXXXXX>` or `kegg:<koXXXXX>`), update the test fixture and the implementation below to match. **Do not invent a new primary-ID format — reuse the existing one.**

- [ ] **Step 4: Implement edge generation**

Add to [multiomics_kg/adapters/metabolism_adapter.py](multiomics_kg/adapters/metabolism_adapter.py):

```python
    def _reaction_metabolite_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        for rxn_id, rxn in data.get("reactions", {}).items():
            for cpd_id in rxn.get("compound_ids", []):
                edge_id = f"r2m:{rxn_id}:{cpd_id}"
                yield (
                    edge_id,
                    f"kegg.reaction:{rxn_id}",
                    f"kegg.compound:{cpd_id}",
                    "reaction_has_metabolite",
                    {},
                )

    def _reaction_pathway_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        for rxn_id, rxn in data.get("reactions", {}).items():
            for pw_id in rxn.get("kegg_pathway_ids", []):
                edge_id = f"r2p:{rxn_id}:{pw_id}"
                yield (
                    edge_id,
                    f"kegg.reaction:{rxn_id}",
                    f"kegg.pathway:{pw_id}",       # see Step 3 — confirm format
                    "reaction_in_kegg_pathway",
                    {},
                )
```

For `_gene_reaction_edges`, add to `MultiMetabolismAdapter`:

```python
    def _gene_reaction_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        known_rxns = set(data.get("reactions", {}).keys())
        for row in load_genome_rows():
            genes = load_gene_annotations(row["data_dir"])
            if not genes:
                continue
            for locus_tag, gene in genes.items():
                for rxn_id in gene.get("kegg_reactions", []) or []:
                    if rxn_id not in known_rxns:
                        continue  # skip dangling: keeps KG free of orphan edges
                    edge_id = f"g2r:{locus_tag}:{rxn_id}"
                    yield (
                        edge_id,
                        f"ncbigene:{locus_tag}",
                        f"kegg.reaction:{rxn_id}",
                        "gene_catalyzes_reaction",
                        {},
                    )

    def get_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        yield from self._gene_reaction_edges()
        yield from self._reaction_metabolite_edges()
        yield from self._reaction_pathway_edges()
```

Also override `MetabolismAdapter.get_edges()` to yield only the two non-strain-scoped edge types (so the unit test `test_reaction_metabolite_edges` works without needing strains):

```python
    def get_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        yield from self._reaction_metabolite_edges()
        yield from self._reaction_pathway_edges()
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_metabolism_adapter.py -v`
Expected: 10 PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/metabolism_adapter.py \
        tests/test_metabolism_adapter.py
git commit -m "metabolism: 1.2 — adapter edge generation (3 edge types)"
```

---

### Task 12: Wire `MultiMetabolismAdapter` into `create_knowledge_graph.py`

**Files:**
- Modify: `create_knowledge_graph.py`

The adapter is *additive* — slot it after the OrthologGroup adapter and before the omics_adapter (so reaction/metabolite nodes exist before any post-import that joins on them).

- [ ] **Step 1: Add the import**

In [create_knowledge_graph.py](create_knowledge_graph.py), after `from multiomics_kg.adapters.brite_adapter import MultiBriteAdapter` (~line 21), add:

```python
from multiomics_kg.adapters.metabolism_adapter import MultiMetabolismAdapter
```

- [ ] **Step 2: Wire the adapter call**

After the OrthologGroup adapter block (~line 67-68) and before the UniProt adapter, add:

```python
    # Metabolism adapter: emits Reaction + Metabolite nodes + 3 edge types.
    # Reads cache/data/kegg/kegg_metabolism_xrefs.json (built by prepare_data step 6).
    # Pure file reader — no SQLite / no REST at build time.
    metabolism_adapter = MultiMetabolismAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    bc.write_nodes(metabolism_adapter.get_nodes())
    bc.write_edges(metabolism_adapter.get_edges())
```

- [ ] **Step 3: Smoke-test the import**

Run: `uv run python -c "from multiomics_kg.adapters.metabolism_adapter import MultiMetabolismAdapter; print(MultiMetabolismAdapter.__module__)"`
Expected: prints the module name; exits 0.

- [ ] **Step 4: Commit**

```bash
git add create_knowledge_graph.py
git commit -m "metabolism: 1.2 — wire MultiMetabolismAdapter into create_knowledge_graph"
```

---

### Task 13: Post-import Cypher (indexes + rollups + Organism_has_metabolite)

**Files:**
- Modify: `scripts/post-import.sh`
- Modify: `scripts/post-import.cypher` (mirror; CLAUDE.md requires byte-identical Cypher)

Add to the existing Group 1 (indexes) and Group 3 (CALL IN TRANSACTIONS / writes) blocks. Mirror the BRITE pattern.

- [ ] **Step 1: Add scalar + full-text indexes (Group 1)**

In [scripts/post-import.sh](scripts/post-import.sh), inside the existing Group 1 cypher-shell heredoc (after the BriteCategory indexes), add:

```cypher
// Reaction
CREATE INDEX reaction_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.id);
CREATE INDEX reaction_kegg_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.kegg_reaction_id);
CREATE INDEX reaction_mnxr_idx IF NOT EXISTS FOR (r:Reaction) ON (r.mnxr_id);
CREATE FULLTEXT INDEX reactionFullText IF NOT EXISTS FOR (r:Reaction) ON EACH [r.name];

// Metabolite
CREATE INDEX metabolite_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.id);
CREATE INDEX metabolite_kegg_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.kegg_compound_id);
CREATE INDEX metabolite_mnxm_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.mnxm_id);
CREATE INDEX metabolite_chebi_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.chebi_id);
CREATE FULLTEXT INDEX metaboliteFullText IF NOT EXISTS FOR (m:Metabolite) ON EACH [m.name];
```

- [ ] **Step 2: Add the rollup queries (Group 3)**

In `scripts/post-import.sh`, inside the third cypher-shell block (writes / CALL IN TRANSACTIONS), add:

```cypher
// ── Metabolism rollups ────────────────────────────────────────────────────

// Reaction.gene_count, organism_count, organisms[]
CALL {
  MATCH (r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
  WITH r, count(DISTINCT g) AS gene_count, collect(DISTINCT g.organism_name) AS organisms
  SET r.gene_count = gene_count,
      r.organism_count = size(organisms),
      r.organisms = organisms
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.gene_count, organism_count
CALL {
  MATCH (m:Metabolite)<-[:Reaction_has_metabolite]-(r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
  WITH m, count(DISTINCT g) AS gene_count, count(DISTINCT g.organism_name) AS organism_count
  SET m.gene_count = gene_count, m.organism_count = organism_count
} IN TRANSACTIONS OF 1000 ROWS;

// Materialize Organism_has_metabolite (2-hop save)
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_catalyzes_reaction]->(:Reaction)
        -[:Reaction_has_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[:Organism_has_metabolite]->(m)
} IN TRANSACTIONS OF 1000 ROWS;

// Organism rollup props
CALL {
  MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
  WITH o, count(DISTINCT m) AS metabolite_count
  SET o.metabolite_count = metabolite_count
} IN TRANSACTIONS OF 100 ROWS;

CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)
        -[:Gene_catalyzes_reaction]->(r:Reaction)
  WITH o, count(DISTINCT r) AS reaction_count
  SET o.reaction_count = reaction_count
} IN TRANSACTIONS OF 100 ROWS;
```

- [ ] **Step 3: Mirror in `post-import.cypher`**

Open [scripts/post-import.cypher](scripts/post-import.cypher) and add the same Cypher (without the `time cypher-shell` shell wrapper). Per CLAUDE.md, the two files must contain *byte-identical* Cypher logic.

- [ ] **Step 4: Smoke-test the post-import script parses**

Run: `bash -n scripts/post-import.sh`
Expected: no output, exit 0.

(The Cypher itself can only be validated against a running graph — that happens in Tasks 14 / 15.)

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "metabolism: 1.2 — post-import indexes + rollups + Organism_has_metabolite"
```

---

### Task 14: KG validity tests + integration smoke + CLAUDE.md + validation gate

**Files:**
- Create: `tests/kg_validity/test_metabolism.py`
- Create: `tests/test_metabolism_smoke.py`
- Modify: `CLAUDE.md`

This single task collapses items 9-12 of the user-listed task buckets because they're all "verify the build works against real data + document". Each subsection is a separate commit.

- [ ] **Step 1: Write KG validity tests**

```python
# tests/kg_validity/test_metabolism.py
"""KG validity tests for the metabolism layer (Reaction + Metabolite + 4 edges).

Requires: a running Neo4j instance with the deployed graph (skipped otherwise).
Mark: @pytest.mark.kg
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.kg]


@pytest.fixture(scope="module")
def neo4j_session():
    from neo4j import GraphDatabase
    try:
        driver = GraphDatabase.driver("bolt://localhost:7687")
        with driver.session() as s:
            s.run("RETURN 1").consume()
        with driver.session() as s:
            yield s
    except Exception as e:
        pytest.skip(f"Neo4j not reachable: {e}")
    finally:
        driver.close()


def test_reaction_count_in_expected_range(neo4j_session):
    """Spec acceptance: ≥ 2,000 reactions (gene-reachable subset of KEGG ~14K)."""
    n = neo4j_session.run("MATCH (r:Reaction) RETURN count(r) AS n").single()["n"]
    assert 2000 <= n <= 6000, f"Reaction count {n} outside expected 2000-6000"


def test_metabolite_count_in_expected_range(neo4j_session):
    """Spec acceptance: ≥ 1,000 metabolites."""
    n = neo4j_session.run("MATCH (m:Metabolite) RETURN count(m) AS n").single()["n"]
    assert 1000 <= n <= 5000, f"Metabolite count {n} outside expected 1000-5000"


def test_all_reactions_have_kegg_primary_id(neo4j_session):
    """Primary IDs all start with kegg.reaction: per Spec 1.2."""
    n_bad = neo4j_session.run(
        "MATCH (r:Reaction) WHERE NOT r.id STARTS WITH 'kegg.reaction:' RETURN count(r) AS n"
    ).single()["n"]
    assert n_bad == 0, f"{n_bad} Reaction nodes lack kegg.reaction: prefix"


def test_metabolite_primary_ids_predominantly_kegg(neo4j_session):
    """Most metabolites should have kegg.compound: prefix; a few may fall back to chebi/mnx."""
    total = neo4j_session.run("MATCH (m:Metabolite) RETURN count(m) AS n").single()["n"]
    kegg_count = neo4j_session.run(
        "MATCH (m:Metabolite) WHERE m.id STARTS WITH 'kegg.compound:' RETURN count(m) AS n"
    ).single()["n"]
    assert kegg_count / total >= 0.95, (
        f"Only {kegg_count}/{total} metabolites have kegg.compound prefix"
    )


def test_kegg_reaction_id_property_matches_id_suffix(neo4j_session):
    """Reaction.kegg_reaction_id should match the suffix of the primary id."""
    n_bad = neo4j_session.run(
        "MATCH (r:Reaction) "
        "WHERE 'kegg.reaction:' + r.kegg_reaction_id <> r.id "
        "RETURN count(r) AS n"
    ).single()["n"]
    assert n_bad == 0


def test_every_reaction_has_at_least_one_metabolite(neo4j_session):
    """Every Reaction should connect to ≥1 Metabolite (some KEGG-only orphans tolerated)."""
    n_bad = neo4j_session.run("""
        MATCH (r:Reaction)
        WHERE NOT EXISTS { MATCH (r)-[:Reaction_has_metabolite]->() }
        RETURN count(r) AS n
    """).single()["n"]
    # Allow up to 5% orphans (compound IDs missing from KEGG link data)
    total = neo4j_session.run("MATCH (r:Reaction) RETURN count(r) AS n").single()["n"]
    assert n_bad / total < 0.05, f"{n_bad}/{total} reactions lack any metabolite edge"


def test_organism_has_metabolite_edges_materialized(neo4j_session):
    """≥ 1K Organism_has_metabolite edges should exist after post-import."""
    n = neo4j_session.run(
        "MATCH ()-[r:Organism_has_metabolite]->() RETURN count(r) AS n"
    ).single()["n"]
    assert n >= 1000, f"Only {n} Organism_has_metabolite edges (expected ≥ 1000)"


def test_pyruvate_kinase_spot_check(neo4j_session):
    """Spec spot-check: kegg.reaction:R00200 (pyruvate kinase) exists with expected props."""
    rec = neo4j_session.run(
        "MATCH (r:Reaction {id: 'kegg.reaction:R00200'}) RETURN r"
    ).single()
    assert rec is not None, "kegg.reaction:R00200 not found"
    rxn = dict(rec["r"])
    assert rxn["kegg_reaction_id"] == "R00200"
    # 4 metabolites: PEP, ADP, pyruvate, ATP
    n_meta = neo4j_session.run(
        "MATCH (:Reaction {id: 'kegg.reaction:R00200'})-[:Reaction_has_metabolite]->(m) "
        "RETURN count(m) AS n"
    ).single()["n"]
    assert n_meta == 4, f"R00200 has {n_meta} metabolites; expected 4"
    # In glycolysis pathway
    has_glycolysis = neo4j_session.run(
        "MATCH (:Reaction {id: 'kegg.reaction:R00200'})-[:Reaction_in_kegg_pathway]->(p) "
        "WHERE p.id CONTAINS 'ko00010' RETURN count(p) AS n"
    ).single()["n"]
    assert has_glycolysis >= 1


def test_reaction_class_values_valid(neo4j_session):
    """reaction_class enum: 'transport' or 'chemical'."""
    bad = neo4j_session.run(
        "MATCH (r:Reaction) WHERE NOT r.reaction_class IN ['transport', 'chemical'] "
        "RETURN count(r) AS n"
    ).single()["n"]
    assert bad == 0


def test_mass_balance_values_valid(neo4j_session):
    """mass_balance enum: 'balanced' or 'unbalanced'."""
    bad = neo4j_session.run(
        "MATCH (r:Reaction) WHERE NOT r.mass_balance IN ['balanced', 'unbalanced'] "
        "RETURN count(r) AS n"
    ).single()["n"]
    assert bad == 0
```

- [ ] **Step 2: Run KG validity locally if Neo4j is up**

Run: `uv run pytest tests/kg_validity/test_metabolism.py -v`
Expected: PASS if Neo4j is running and the graph has been rebuilt with the new adapter; SKIP otherwise.

If Neo4j is not running yet, defer this assertion until the next end-to-end build (covered in Step 5 below).

- [ ] **Step 3: Commit KG validity test**

```bash
git add tests/kg_validity/test_metabolism.py
git commit -m "metabolism: 1.2 — KG validity tests (Reaction + Metabolite + edge counts + spot checks)"
```

- [ ] **Step 4: Write integration smoke test (slow marker)**

```python
# tests/test_metabolism_smoke.py
"""Integration smoke test: run prepare_data step 6 + adapter against the real cache.

Slow (~2 min). Skipped unless invoked with -m slow or "slow or kg".
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parent.parent
KEGG_DATA = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_data.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"
XREFS_OUT = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_metabolism_xrefs.json"


@pytest.mark.slow
def test_step6_smoke_run_against_real_cache():
    if not KEGG_DATA.exists():
        pytest.skip(f"{KEGG_DATA} missing — run KEGG download first")
    if not RESOLVER_DB.exists():
        pytest.skip(f"{RESOLVER_DB} missing — run prepare_data step 0 sub-step 7 first")

    result = subprocess.run(
        [sys.executable, "-m", "multiomics_kg.download.build_kegg_metabolism_xrefs", "--force"],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    assert result.returncode == 0, f"step 6 failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    assert XREFS_OUT.exists()
    # Sanity bounds
    import json
    data = json.loads(XREFS_OUT.read_text())
    assert 2000 <= len(data["reactions"]) <= 6000
    assert 1000 <= len(data["compounds"]) <= 5000
    # All keys are R-numbers / C-numbers
    assert all(k.startswith("R") for k in data["reactions"])
    assert all(k.startswith("C") for k in data["compounds"])


@pytest.mark.slow
def test_metabolism_adapter_against_real_xrefs():
    if not XREFS_OUT.exists():
        pytest.skip(f"{XREFS_OUT} missing — run step 6 first")
    from multiomics_kg.adapters.metabolism_adapter import MultiMetabolismAdapter

    adapter = MultiMetabolismAdapter(
        genome_config_file=str(PROJECT_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"),
    )
    nodes = list(adapter.get_nodes())
    assert any(n[1] == "reaction" for n in nodes)
    assert any(n[1] == "metabolite" for n in nodes)

    edges = list(adapter.get_edges())
    labels = {e[3] for e in edges}
    assert labels == {"gene_catalyzes_reaction", "reaction_has_metabolite", "reaction_in_kegg_pathway"}
```

- [ ] **Step 5: Run smoke tests against real cache**

Run: `uv run pytest tests/test_metabolism_smoke.py -m slow -v`
Expected: PASS (both tests). If `kegg_data.json` is missing, run `uv run python create_knowledge_graph.py --test` first to populate it.

- [ ] **Step 6: Commit smoke tests**

```bash
git add tests/test_metabolism_smoke.py
git commit -m "metabolism: 1.2 — integration smoke (real cache, slow marker)"
```

- [ ] **Step 7: Validation gate against MED4 (manual)**

Rebuild the graph end-to-end:

```bash
# 1. Run prepare-data step 6
bash scripts/prepare_data.sh --steps 6 --force

# 2. Confirm xrefs file size is reasonable (~1 MB)
ls -lh cache/data/kegg/kegg_metabolism_xrefs.json

# 3. Run KG build (test mode for speed; remove --test for full run)
uv run python create_knowledge_graph.py --test

# 4. Bring up Docker stack to import + post-process
docker compose down
docker compose up -d build
docker compose up -d import post-process deploy

# 5. Wait for post-process to complete, then run KG validity
docker compose logs -f post-process | grep -E "Done|ERROR" &
# Wait for "Post-process: Done" then Ctrl-C the log follow
uv run pytest -m kg -v
```

Acceptance criteria from the spec:
- Reaction nodes ≥ 2,000 with primary IDs all `kegg.reaction:R*`
- Metabolite nodes ≥ 1,000 with primary IDs predominantly `kegg.compound:C*`
- `Gene_catalyzes_reaction` edges: at least one per gene with non-empty `kegg_reactions`
- `Reaction_has_metabolite` edges: average ~4 per reaction
- `Reaction_in_kegg_pathway` edges: present and reusing existing `KeggTerm` pathway nodes
- `Organism_has_metabolite` edges: materialized for all 25 genome organisms

Run a focused MED4 spot-check via cypher-shell:

```cypher
MATCH (o:OrganismTaxon {preferred_name: 'Prochlorococcus MED4'})
      -[:Organism_has_metabolite]->(m:Metabolite)
RETURN count(m) AS med4_metabolites;
```
Expected: roughly 800-1500 metabolites.

```cypher
MATCH (o:OrganismTaxon {preferred_name: 'Prochlorococcus MED4'})
      <-[:Gene_belongs_to_organism]-(g:Gene)
      -[:Gene_catalyzes_reaction]->(r:Reaction)
RETURN count(DISTINCT r) AS med4_reactions;
```
Expected: roughly 1500-2500 reactions.

If either count is wildly off (zero, or 10x larger), STOP and inspect the import report + post-import logs before proceeding.

- [ ] **Step 8: Update `CLAUDE.md`**

In [CLAUDE.md](CLAUDE.md):

1. In the **Architecture → Key Adapters** list, add a bullet for `metabolism_adapter.py`.
2. In the **Genome Data Download Pipeline → Step 6** description, replace the Phase 1.1A/B placeholder text with a description of `build_kegg_metabolism_xrefs.py`.
3. In the **Notes → Key graph facts** list (the "BioCypher output: ..." paragraph), append the new node and edge counts: e.g. `~3K Reaction nodes, ~1.5K Metabolite nodes, ~12K Gene_catalyzes_reaction edges, ~12K Reaction_has_metabolite edges, ~3K Reaction_in_kegg_pathway edges, ~XXK Organism_has_metabolite edges.`
4. In the **Actual Neo4j labels** list, add `Reaction`, `Metabolite` to nodes and `Gene_catalyzes_reaction`, `Reaction_has_metabolite`, `Reaction_in_kegg_pathway`, `Organism_has_metabolite` to relationships.
5. In **Post-import indexes**, add the new scalar + full-text indexes.
6. In **Notes**, document the `reaction_class` and `mass_balance` enums.

- [ ] **Step 9: Commit doc updates**

```bash
git add CLAUDE.md
git commit -m "metabolism: 1.2 — CLAUDE.md docs (new nodes, edges, step 6, indexes)"
```

---

## Self-review checklist

Before declaring this plan complete, verify:

1. **Spec coverage:** every section of the revised spec has at least one task above.
   - [x] 5 KEGG endpoints → Tasks 1-2
   - [x] YAML transform drop → Task 3
   - [x] DDL indexes → Task 4
   - [x] Primary-ID helpers → Task 5
   - [x] step 6 builder → Tasks 6-7
   - [x] prepare_data wiring → Task 8
   - [x] schema → Task 9
   - [x] adapter → Tasks 10-11
   - [x] create_knowledge_graph wire → Task 12
   - [x] post-import Cypher → Task 13
   - [x] KG validity tests → Task 14
   - [x] integration smoke + validation gate + CLAUDE.md → Task 14

2. **No placeholders:** every step has actual code or exact commands.

3. **Type / name consistency** across tasks:
   - `mnxm_to_primary_id` / `mnxr_to_primary_id` (Task 5) → consumed by Task 7's enrichment? **No** — Task 7's enrichment uses `_resolve_mnx_for_kegg_*` (forward direction) directly. The primary-ID helpers are for non-KEGG-primary fallback (compounds with no KEGG entry). Currently the enrichment code in Task 7 always emits `kegg.reaction:R*` / `kegg.compound:C*` primary IDs because the pruning starts from KEGG IDs. The helpers will be load-bearing for **Phase 2** paper integration. **Confirmed: no inconsistency.**
   - `mass_balance` / `reaction_class` enums match between adapter (Task 10), spec, and KG validity tests (Task 14). ✅
   - `kegg.pathway:` primary-ID prefix on `Reaction_in_kegg_pathway` target — flagged in Task 11 Step 3 as "verify before implementing".

4. **Commit cadence:** 14 tasks → 14+ commits, each is reviewable in isolation.

5. **Order discipline:** Task 3 (drop transform) sequenced *first* among the surgical changes so Tasks 6+ build on raw R-numbers. Task 4 (DDL indexes) before Task 5 (helpers that use them) before Task 7 (enrichment that calls helpers' SQL patterns). ✅
