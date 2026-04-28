# Phase 1.1B — Metabolite foundation: resolver + accessors + step-2 transforms

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the four Phase 1.1A skeletons with real implementations, plus add the three step-2 transforms + YAML wiring + per-strain metabolism report. End state: `gene_annotations_merged.json` contains MNXR-resolved `kegg_reactions`, TCDB-validated `transporter_classification`, and CAZy-validated `cazy_ids` for every gene.

**Architecture:** sub-step 7 (`build_metabolite_resolver.py`) parses the seven sub-step-6 cache files plus raw eggNOG output (for CAZy bootstrap), producing a SQLite resolver DB + two hierarchy JSONs + a diagnostic report. Three accessor modules in `utils/` wrap the cache files with a clean read API. Three new per-token transforms in `annotation_transforms.py` consume the accessors and slot into the existing YAML transforms framework via three `transforms:` lines in `gene_annotations_config.yaml`.

**Tech Stack:** Python 3.10+, `sqlite3` (stdlib), `pytest` with `tmp_path` for synthetic fixtures, the existing YAML transforms framework in `multiomics_kg/download/utils/annotation_transforms.py`.

**Spec:** [`docs/superpowers/specs/2026-04-28-metabolite-foundation-design.md`](../specs/2026-04-28-metabolite-foundation-design.md) — schemas confirmed against MNXref 4.5 + live TCDB + eggNOG v2.x. Real cache files already on disk in `cache/data/{mnx,tcdb}/` from Phase 1.1A.

---

## File structure

**Create / fully implement (skeleton bodies replaced):**

| Path | Role |
|---|---|
| `multiomics_kg/download/build_metabolite_resolver.py` | sub-step 7 builder (replaces Phase 1.1A skeleton) |
| `multiomics_kg/utils/metabolite_utils.py` | resolver accessor (replaces skeleton) |
| `multiomics_kg/utils/tcdb_utils.py` | TCDB hierarchy accessor (replaces skeleton) |
| `multiomics_kg/utils/cazy_utils.py` | CAZy hierarchy accessor (replaces skeleton) |

**Modify:**

| Path | Why |
|---|---|
| `multiomics_kg/download/utils/annotation_transforms.py` | add three `_tx_*` functions |
| `config/gene_annotations_config.yaml` | rename `kegg_reaction` → `kegg_reactions`; add `transforms:` lines on 3 fields |
| `multiomics_kg/download/build_gene_annotations.py` | invoke per-strain metabolism report at end of build |

**Test files (all new):**

| Path | Role |
|---|---|
| `tests/test_build_metabolite_resolver.py` | sub-step 7 unit tests (synthetic input, in-memory SQLite) |
| `tests/test_metabolite_utils.py` | resolver accessor tests |
| `tests/test_tcdb_utils.py` | TCDB accessor tests |
| `tests/test_cazy_utils.py` | CAZy accessor tests |
| `tests/test_annotation_transforms_metabolism.py` | the 3 new transforms |
| `tests/test_build_gene_annotations_metabolism.py` | end-to-end YAML transform integration |
| `tests/test_prepare_data_step2_metabolism_smoke.py` | integration smoke (`@pytest.mark.slow`) |

**Skeleton tests removed at the end:** `tests/test_metabolism_skeletons.py` is dropped once the four modules have real tests covering them.

---

## Tasks

### Task 1: SQLite `compounds` table from chem_prop.tsv

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (replace skeleton body)
- Test: `tests/test_build_metabolite_resolver.py` (new)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_metabolite_resolver.py
"""Unit tests for sub-step 7 resolver builder."""
from __future__ import annotations

import sqlite3
import textwrap
from pathlib import Path

from multiomics_kg.download import build_metabolite_resolver as bmr


CHEM_PROP_FIXTURE = textwrap.dedent("""\
    ### MetaNetX/MNXref reconciliation ###
    #VERSION:   4.5
    #DATE:      2025/08/13
    #ID\tname\treference\tformula\tcharge\tmass\tInChI\tInChIKey\tSMILES
    MNXM01\tPMF\tmnx:PMF\tH\t1\t1.00794\tInChI=1S/p+1\tGPRLSGONYQIRFK-UHFFFAOYSA-N\t[H+]
    MNXM41\tD-glucose\tchebi:17234\tC6H12O6\t0\t180.063\tInChI=1S/C6H12O6\tWQZGKKKJIJFFOK-GASJEMHNSA-N\tOC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O
    BIOMASS\tBIOMASS\tmnx:BIOMASS\t\t\t\t\t\t
""")


def test_build_compounds_table(tmp_path):
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    db_path = tmp_path / "resolver.db"

    conn = sqlite3.connect(db_path)
    bmr.build_compounds_table(conn, chem_prop)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles FROM compounds ORDER BY mnxm_id")
    rows = cur.fetchall()

    # 3 rows expected (BIOMASS, MNXM01, MNXM41), but BIOMASS has empty fields
    assert len(rows) == 3
    biomass = [r for r in rows if r[0] == "BIOMASS"][0]
    assert biomass[1] == "BIOMASS"
    assert biomass[3] == ""  # empty formula

    pmf = [r for r in rows if r[0] == "MNXM01"][0]
    assert pmf[1] == "PMF"
    assert pmf[2] == "mnx:PMF"
    assert pmf[3] == "H"
    assert pmf[4] == 1
    assert abs(pmf[5] - 1.00794) < 1e-5
    assert pmf[7] == "GPRLSGONYQIRFK-UHFFFAOYSA-N"

    glucose = [r for r in rows if r[0] == "MNXM41"][0]
    assert glucose[1] == "D-glucose"
    assert glucose[3] == "C6H12O6"
    assert abs(glucose[5] - 180.063) < 1e-3

    conn.close()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_compounds_table -v`
Expected: FAIL — `AttributeError: module 'build_metabolite_resolver' has no attribute 'build_compounds_table'`.

- [ ] **Step 3: Replace the skeleton with the real builder header + compounds-table function**

Open `multiomics_kg/download/build_metabolite_resolver.py` and **replace its entire contents** with:

```python
"""Step 0 sub-step 7 — Build resolver + hierarchy caches.

Reads the seven sub-step-6 cache files (4 MNX TSVs + 3 TCDB TSVs) plus raw
eggNOG annotation files from configured strains, and produces:

- cache/data/mnx/metabolite_resolver.db    (SQLite)
- cache/data/tcdb/tcdb_hierarchy.json
- cache/data/cazy/cazy_hierarchy.json
- cache/data/mnx/metabolite_id_mapping_report.json
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sqlite3
from pathlib import Path

log = logging.getLogger(__name__)


# ── compounds ─────────────────────────────────────────────────────────────────

_COMPOUNDS_DDL = """
CREATE TABLE IF NOT EXISTS compounds (
    mnxm_id   TEXT PRIMARY KEY,
    name      TEXT,
    reference TEXT,
    formula   TEXT,
    charge    INTEGER,
    mass      REAL,
    inchi     TEXT,
    inchikey  TEXT,
    smiles    TEXT
);
"""


def _iter_data_rows(path: Path):
    """Yield non-comment, non-blank lines from a tab-delimited MNX file."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def _coerce_int(s: str) -> int | None:
    return int(s) if s and s.lstrip("-").isdigit() else None


def _coerce_float(s: str) -> float | None:
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def build_compounds_table(conn: sqlite3.Connection, chem_prop_path: Path) -> int:
    """Populate the `compounds` table from chem_prop.tsv. Returns row count."""
    cur = conn.cursor()
    cur.executescript(_COMPOUNDS_DDL)
    n = 0
    for fields in _iter_data_rows(chem_prop_path):
        # 9 columns expected; pad with empties if shorter
        fields = (fields + [""] * 9)[:9]
        mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles = fields
        cur.execute(
            "INSERT OR REPLACE INTO compounds "
            "(mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (mnxm_id, name, reference, formula,
             _coerce_int(charge), _coerce_float(mass),
             inchi, inchikey, smiles),
        )
        n += 1
    log.info(f"  compounds: {n} rows")
    return n


def main(force: bool = False) -> None:
    raise NotImplementedError("Phase 1.1B — see plan for follow-up tasks")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild caches even if they exist.")
    args = parser.parse_args()
    main(force=args.force)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_compounds_table -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.1B — compounds table from chem_prop.tsv"
```

---

### Task 2: SQLite alias + name tables from chem_xref.tsv

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (add 3 functions + a normalization map)
- Test: `tests/test_build_metabolite_resolver.py` (append)

- [ ] **Step 1: Append the failing test**

```python
# Append to tests/test_build_metabolite_resolver.py

CHEM_XREF_FIXTURE = textwrap.dedent("""\
    #source\tID\tdescription
    MNXM41\tMNXM41\tD-glucose
    chebi:17234\tMNXM41\tD-glucose||Dextrose
    CHEBI:17234\tMNXM41\tD-glucose
    kegg.compound:C00031\tMNXM41\tD-glucose
    keggC:C00031\tMNXM41\tD-glucose
    bigg.metabolite:glc__D\tMNXM41\tD-glucose
    biggM:glc__D\tMNXM41\tD-glucose
    mnx:PMF\tMNXM01\tPMF
""")


def test_build_compound_aliases_normalizes_dual_prefixes(tmp_path):
    """Dual short/long source forms (CHEBI/chebi, keggC/kegg.compound) collapse
    to a single canonical prefix; values dedup per (source, value, mnxm)."""
    # First populate compounds (alias FK target)
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(CHEM_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    bmr.build_compounds_table(conn, chem_prop)

    bmr.build_compound_aliases_table(conn, chem_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT source, value, mnxm_id FROM compound_aliases ORDER BY source, value")
    rows = cur.fetchall()

    # CHEBI:17234 should be normalized to source='chebi', value='17234', and
    # both `chebi:17234` and `CHEBI:17234` rows should collapse to one
    chebi_rows = [r for r in rows if r[0] == "chebi"]
    assert len(chebi_rows) == 1
    assert chebi_rows[0] == ("chebi", "17234", "MNXM41")

    # kegg.compound:C00031 and keggC:C00031 collapse to one
    kegg_rows = [r for r in rows if r[0] == "kegg.compound"]
    assert kegg_rows == [("kegg.compound", "C00031", "MNXM41")]

    # bigg.metabolite/biggM collapse
    bigg_rows = [r for r in rows if r[0] == "bigg.metabolite"]
    assert bigg_rows == [("bigg.metabolite", "glc__D", "MNXM41")]

    # Self-references (bare MNXM41) skipped
    assert not any(r[0] == "" for r in rows)

    conn.close()


def test_build_compound_names_normalizes(tmp_path):
    """Names from `name` and description (||-split) get normalized + indexed."""
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(CHEM_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    bmr.build_compounds_table(conn, chem_prop)
    bmr.build_compound_names_table(conn, chem_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT DISTINCT name_normalized, mnxm_id FROM compound_names "
                "WHERE mnxm_id = 'MNXM41' ORDER BY name_normalized")
    rows = cur.fetchall()
    names = {r[0] for r in rows}
    assert "d-glucose" in names
    assert "dextrose" in names

    conn.close()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_build_metabolite_resolver.py -v`
Expected: 2 new tests FAIL with `AttributeError: build_compound_aliases_table` / `build_compound_names_table`.

- [ ] **Step 3: Add the source normalization map and the two builder functions**

Insert in `multiomics_kg/download/build_metabolite_resolver.py` between `build_compounds_table` and `def main`:

```python
# Source-prefix normalization: collapse MNX's dual short/long forms.
# Long bioregistry-style spelling wins. Keys: any prefix MNX uses; values:
# canonical form. Prefixes not in the map are passed through unchanged.
_SOURCE_NORM = {
    "CHEBI":           "chebi",
    "chebi":           "chebi",
    "keggC":           "kegg.compound",
    "kegg.compound":   "kegg.compound",
    "keggD":           "kegg.drug",
    "kegg.drug":       "kegg.drug",
    "keggG":           "kegg.glycan",
    "kegg.glycan":     "kegg.glycan",
    "biggM":           "bigg.metabolite",
    "bigg.metabolite": "bigg.metabolite",
    "metacycM":        "metacyc.compound",
    "metacyc.compound":"metacyc.compound",
    "seedM":           "seed.compound",
    "seed.compound":   "seed.compound",
    "vmhM":            "vmh",
    "vmh":             "vmh",
    "lipidmapsM":      "lipidmaps",
    "lipidmaps":       "lipidmaps",
    "SLM":             "slm",
    "slm":             "slm",
    "reactomeM":       "reactome",
    "reactome":        "reactome",
    "envipathM":       "envipath",
    "envipath":        "envipath",
    "sabiorkM":        "sabiork",
    "sabiork":         "sabiork",
    # Reaction-side equivalents
    "keggR":           "kegg.reaction",
    "kegg.reaction":   "kegg.reaction",
    "biggR":           "bigg.reaction",
    "bigg.reaction":   "bigg.reaction",
    "seedR":           "seed.reaction",
    "seed.reaction":   "seed.reaction",
    "vmhreaction":     "vmhreaction",
    "vmhR":            "vmhreaction",
    "metacycR":        "metacyc.reaction",
    "metacyc.reaction":"metacyc.reaction",
    "sabiork.reaction":"sabiork.reaction",
    "rh":              "rhea",
    "rhea":            "rhea",
}


def _split_source_value(combined: str) -> tuple[str, str] | None:
    """Split MNX xref col-1 ('chebi:17234') into (canonical_source, value).

    Returns None for self-references (bare MNXM/MNXR ids without a colon).
    """
    if ":" not in combined:
        return None
    src, _, value = combined.partition(":")
    canonical = _SOURCE_NORM.get(src, src)
    return canonical, value


_ALIAS_DDL = """
CREATE TABLE IF NOT EXISTS compound_aliases (
    source  TEXT,
    value   TEXT,
    mnxm_id TEXT,
    PRIMARY KEY (source, value, mnxm_id)
);
CREATE INDEX IF NOT EXISTS idx_compound_aliases_value
    ON compound_aliases(value);
CREATE INDEX IF NOT EXISTS idx_compound_aliases_source_value
    ON compound_aliases(source, value);
"""


def build_compound_aliases_table(conn: sqlite3.Connection, chem_xref_path: Path) -> int:
    """Populate compound_aliases from chem_xref.tsv with source normalization."""
    cur = conn.cursor()
    cur.executescript(_ALIAS_DDL)
    n = 0
    for fields in _iter_data_rows(chem_xref_path):
        if len(fields) < 2:
            continue
        combined, mnxm_id = fields[0], fields[1]
        parsed = _split_source_value(combined)
        if parsed is None:
            continue  # skip self-references
        source, value = parsed
        cur.execute(
            "INSERT OR IGNORE INTO compound_aliases (source, value, mnxm_id) "
            "VALUES (?, ?, ?)",
            (source, value, mnxm_id),
        )
        n += 1
    log.info(f"  compound_aliases: {n} rows scanned")
    return n


_NAMES_DDL = """
CREATE TABLE IF NOT EXISTS compound_names (
    name_normalized TEXT,
    mnxm_id         TEXT,
    PRIMARY KEY (name_normalized, mnxm_id)
);
CREATE INDEX IF NOT EXISTS idx_compound_names
    ON compound_names(name_normalized);
"""


_NAME_NORM_RE = re.compile(r"[^\w+\-/ ]")


def _normalize_name(s: str) -> str:
    """Lowercase, collapse whitespace, strip punctuation except + - /."""
    s = _NAME_NORM_RE.sub(" ", s.lower())
    return " ".join(s.split())


def build_compound_names_table(conn: sqlite3.Connection, chem_xref_path: Path) -> int:
    """Populate compound_names from the description column of chem_xref.tsv.

    Description uses '||' as alternate-name separator. Each non-empty alt-name
    contributes one (name_normalized, mnxm_id) row.
    """
    cur = conn.cursor()
    cur.executescript(_NAMES_DDL)
    n = 0
    for fields in _iter_data_rows(chem_xref_path):
        if len(fields) < 3:
            continue
        mnxm_id, description = fields[1], fields[2]
        for raw in description.split("||"):
            normalized = _normalize_name(raw)
            if not normalized:
                continue
            cur.execute(
                "INSERT OR IGNORE INTO compound_names (name_normalized, mnxm_id) "
                "VALUES (?, ?)",
                (normalized, mnxm_id),
            )
            n += 1
    log.info(f"  compound_names: {n} rows scanned")
    return n
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_build_metabolite_resolver.py -v`
Expected: 3 tests PASS (the original `test_build_compounds_table` plus the 2 new ones).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.1B — alias + name tables from chem_xref.tsv with source normalization"
```

---

### Task 3: SQLite reactions + reaction_aliases tables

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (add reactions + reaction_aliases builders)
- Test: `tests/test_build_metabolite_resolver.py` (append)

- [ ] **Step 1: Append the failing tests**

```python
# Append to tests/test_build_metabolite_resolver.py

REAC_PROP_FIXTURE = textwrap.dedent("""\
    #ID\tmnx_equation\treference\tclassifs\tis_balanced\tis_transport
    MNXR101234\t1 MNXM3@MNXD1 + 1 MNXM41@MNXD1 = 1 MNXM7@MNXD1 + 1 MNXM58@MNXD1\tmnx:MNXR101234\t2.7.1.1\tB\t
    MNXR02\t1 MNXM1@MNXD1 = 1 MNXM1@MNXD2\tmnx:MNXR02\t\tB\tT
    EMPTY\t = \tmnx:EMPTY\t\tB\t
""")


REAC_XREF_FIXTURE = textwrap.dedent("""\
    #source\tID\tdescription
    MNXR101234\tMNXR101234\thexokinase
    kegg.reaction:R00299\tMNXR101234\thexokinase
    keggR:R00299\tMNXR101234\thexokinase
    rhea:16332\tMNXR101234\thexokinase
    rh:16332\tMNXR101234\thexokinase
""")


def test_build_reactions_table(tmp_path):
    """Reactions table mirrors reac_prop.tsv 6 columns."""
    reac_prop = tmp_path / "reac_prop.tsv"
    reac_prop.write_text(REAC_PROP_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)

    bmr.build_reactions_table(conn, reac_prop)
    conn.commit()

    cur = conn.cursor()
    cur.execute(
        "SELECT mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport "
        "FROM reactions ORDER BY mnxr_id"
    )
    rows = cur.fetchall()
    assert len(rows) == 3

    # Look up MNXR101234
    hk = [r for r in rows if r[0] == "MNXR101234"][0]
    assert "MNXM41" in hk[1]
    assert "MNXM7" in hk[1]
    assert hk[3] == "2.7.1.1"
    assert hk[4] == "B"
    assert hk[5] in (None, "")

    # MNXR02 is a transport reaction
    transport = [r for r in rows if r[0] == "MNXR02"][0]
    assert transport[5] == "T"

    conn.close()


def test_build_reaction_aliases_normalizes(tmp_path):
    """Reaction xrefs canonicalize keggR→kegg.reaction and rh→rhea."""
    reac_xref = tmp_path / "reac_xref.tsv"
    reac_xref.write_text(REAC_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)

    bmr.build_reaction_aliases_table(conn, reac_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT source, value, mnxr_id FROM reaction_aliases ORDER BY source, value")
    rows = cur.fetchall()

    kegg = [r for r in rows if r[0] == "kegg.reaction"]
    assert kegg == [("kegg.reaction", "R00299", "MNXR101234")]

    rhea = [r for r in rows if r[0] == "rhea"]
    assert rhea == [("rhea", "16332", "MNXR101234")]

    # Self-reference dropped
    assert not any(r[0] == "" for r in rows)

    conn.close()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_reactions_table tests/test_build_metabolite_resolver.py::test_build_reaction_aliases_normalizes -v`
Expected: FAIL — `AttributeError: build_reactions_table` / `build_reaction_aliases_table`.

- [ ] **Step 3: Add the two reaction-table builders**

Insert in `multiomics_kg/download/build_metabolite_resolver.py` after `build_compound_names_table`:

```python
# ── reactions ─────────────────────────────────────────────────────────────────

_REACTIONS_DDL = """
CREATE TABLE IF NOT EXISTS reactions (
    mnxr_id      TEXT PRIMARY KEY,
    mnx_equation TEXT,
    reference    TEXT,
    classifs     TEXT,
    is_balanced  TEXT,
    is_transport TEXT
);
"""


def build_reactions_table(conn: sqlite3.Connection, reac_prop_path: Path) -> int:
    """Populate the `reactions` table from reac_prop.tsv (6 columns)."""
    cur = conn.cursor()
    cur.executescript(_REACTIONS_DDL)
    n = 0
    for fields in _iter_data_rows(reac_prop_path):
        fields = (fields + [""] * 6)[:6]
        mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport = fields
        cur.execute(
            "INSERT OR REPLACE INTO reactions "
            "(mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (mnxr_id, mnx_equation, reference, classifs,
             is_balanced or None, is_transport or None),
        )
        n += 1
    log.info(f"  reactions: {n} rows")
    return n


_REAC_ALIAS_DDL = """
CREATE TABLE IF NOT EXISTS reaction_aliases (
    source  TEXT,
    value   TEXT,
    mnxr_id TEXT,
    PRIMARY KEY (source, value, mnxr_id)
);
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_value
    ON reaction_aliases(value);
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_source_value
    ON reaction_aliases(source, value);
"""


def build_reaction_aliases_table(conn: sqlite3.Connection, reac_xref_path: Path) -> int:
    """Populate reaction_aliases from reac_xref.tsv with source normalization."""
    cur = conn.cursor()
    cur.executescript(_REAC_ALIAS_DDL)
    n = 0
    for fields in _iter_data_rows(reac_xref_path):
        if len(fields) < 2:
            continue
        combined, mnxr_id = fields[0], fields[1]
        parsed = _split_source_value(combined)
        if parsed is None:
            continue
        source, value = parsed
        cur.execute(
            "INSERT OR IGNORE INTO reaction_aliases (source, value, mnxr_id) "
            "VALUES (?, ?, ?)",
            (source, value, mnxr_id),
        )
        n += 1
    log.info(f"  reaction_aliases: {n} rows scanned")
    return n
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_build_metabolite_resolver.py -v`
Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.1B — reactions + reaction_aliases tables"
```

---

### Task 4: TCDB hierarchy JSON from 3 TCDB sources

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (add `build_tcdb_hierarchy`)
- Test: `tests/test_build_metabolite_resolver.py` (append)

- [ ] **Step 1: Append the failing test**

```python
# Append to tests/test_build_metabolite_resolver.py

TCDB_FAMILIES_FIXTURE = textwrap.dedent("""\
    1.A.1\tThe Voltage-gated Ion Channel (VIC) Superfamily 
    1.A.10\tThe Glutamate-gated Ion Channel (GIC) Family of Neurotransmitter Receptors
""")


TCDB_SUPERFAMILIES_FIXTURE = textwrap.dedent("""\
    #TCID\tSubfamily\tFamily\tFam_abbreviation\tSuperfamily
    1.A.1.1.1\t1.A.1.1\t1.A.1\tVIC\tVIC Superfamily
    1.A.1.1.2\t1.A.1.1\t1.A.1\tVIC\tVIC Superfamily
    1.A.10.1.1\t1.A.10.1\t1.A.10\tGIC\tGIC Superfamily
""")


TCDB_SUBSTRATES_FIXTURE = textwrap.dedent("""\
    1.A.1.1.1\tCHEBI:29103;potassium(1+)
    1.A.10.1.1\tCHEBI:29987;glutamate(2-)|CHEBI:33709;amino acid
""")


def test_build_tcdb_hierarchy(tmp_path):
    """TCDB hierarchy JSON joins families + superfamilies + substrates."""
    fams = tmp_path / "families.tsv"
    fams.write_text(TCDB_FAMILIES_FIXTURE)
    supers = tmp_path / "superfamilies.tsv"
    supers.write_text(TCDB_SUPERFAMILIES_FIXTURE)
    subs = tmp_path / "substrates.tsv"
    subs.write_text(TCDB_SUBSTRATES_FIXTURE)
    out = tmp_path / "tcdb_hierarchy.json"

    bmr.build_tcdb_hierarchy(out, fams, supers, subs)

    h = json.loads(out.read_text())

    # Class (level 0) synthesized
    assert h["1"]["level"] == 0
    assert h["1"]["level_kind"] == "tc_class"
    assert h["1"]["parent"] is None

    # Subclass (level 1) synthesized
    assert h["1.A"]["level"] == 1
    assert h["1.A"]["parent"] == "1"

    # Family (level 2) — name from families.tsv
    assert h["1.A.1"]["level"] == 2
    assert "Voltage-gated Ion Channel" in h["1.A.1"]["name"]
    assert h["1.A.1"]["parent"] == "1.A"
    assert h["1.A.1"]["abbreviation"] == "VIC"

    # Subfamily (level 3) — derived from col 2 of superfamilies
    assert h["1.A.1.1"]["level"] == 3
    assert h["1.A.1.1"]["parent"] == "1.A.1"

    # Specificity (level 4) — substrate joined
    sp = h["1.A.1.1.1"]
    assert sp["level"] == 4
    assert sp["parent"] == "1.A.1.1"
    assert sp["substrate_classes"] == ["potassium(1+)"]
    assert sp["superfamily"] == "VIC Superfamily"

    # Multi-substrate
    glu = h["1.A.10.1.1"]
    assert set(glu["substrate_classes"]) == {"glutamate(2-)", "amino acid"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_tcdb_hierarchy -v`
Expected: FAIL — `AttributeError: build_tcdb_hierarchy`.

- [ ] **Step 3: Add the TCDB hierarchy builder**

Insert in `build_metabolite_resolver.py` after the reaction builders:

```python
# ── TCDB hierarchy ────────────────────────────────────────────────────────────

_TC_CLASS_NAMES = {
    "1": "Channels and Pores",
    "2": "Electrochemical Potential-driven Transporters",
    "3": "Primary Active Transporters",
    "4": "Group Translocators",
    "5": "Transmembrane Electron Carriers",
    "8": "Auxiliary Transport Proteins",
    "9": "Incompletely Characterized Transport Systems",
}


def _parse_tcdb_families(path: Path) -> dict[str, str]:
    """Parse families.tsv → {tc_family_id: description}."""
    out: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                out[parts[0]] = parts[1].strip()
    return out


def _parse_tcdb_superfamilies(path: Path) -> list[tuple[str, str, str, str, str]]:
    """Parse superfamilies.tsv → [(tcid, subfamily, family, abbreviation, superfamily), ...]."""
    rows: list[tuple[str, str, str, str, str]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                rows.append((parts[0], parts[1], parts[2], parts[3], parts[4]))
    return rows


def _parse_tcdb_substrates(path: Path) -> dict[str, list[str]]:
    """Parse substrates.tsv → {tcid_specificity: [substrate_name, ...]}.

    Substrate column is pipe-separated 'CHEBI:N;name|CHEBI:N;name'. We keep
    just the name (post-`;`) for substrate_classes; full CHEBI mapping is a
    Phase 1.3 enhancement.
    """
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            tcid, raw = parts
            names = []
            for entry in raw.split("|"):
                if ";" in entry:
                    _, _, name = entry.partition(";")
                    if name.strip():
                        names.append(name.strip())
            if names:
                out[tcid] = names
    return out


def build_tcdb_hierarchy(
    out_path: Path,
    families_path: Path,
    superfamilies_path: Path,
    substrates_path: Path,
) -> int:
    """Build tcdb_hierarchy.json by joining the 3 TCDB sources."""
    fam_descs = _parse_tcdb_families(families_path)
    super_rows = _parse_tcdb_superfamilies(superfamilies_path)
    substrates = _parse_tcdb_substrates(substrates_path)

    h: dict[str, dict] = {}

    # Synthesize class + subclass entries from any TCID prefix we observe
    seen_classes: set[str] = set()
    seen_subclasses: set[str] = set()

    for tcid, subfam, fam, abbr, superfam in super_rows:
        parts = tcid.split(".")
        # Need at least class.subclass.family.subfam.specificity = 5 parts
        if len(parts) < 5:
            continue
        cls = parts[0]
        subcls = ".".join(parts[:2])

        # Class (level 0)
        if cls not in seen_classes:
            h[cls] = {
                "name": _TC_CLASS_NAMES.get(cls, ""),
                "level": 0,
                "level_kind": "tc_class",
                "parent": None,
            }
            seen_classes.add(cls)

        # Subclass (level 1)
        if subcls not in seen_subclasses:
            h[subcls] = {
                "name": "",
                "level": 1,
                "level_kind": "tc_subclass",
                "parent": cls,
            }
            seen_subclasses.add(subcls)

        # Family (level 2)
        if fam not in h:
            h[fam] = {
                "name": fam_descs.get(fam, ""),
                "level": 2,
                "level_kind": "tc_family",
                "parent": subcls,
                "abbreviation": abbr,
            }

        # Subfamily (level 3)
        if subfam not in h:
            h[subfam] = {
                "name": "",
                "level": 3,
                "level_kind": "tc_subfamily",
                "parent": fam,
            }

        # Specificity (level 4)
        node: dict = {
            "name": "",
            "level": 4,
            "level_kind": "tc_specificity",
            "parent": subfam,
        }
        if tcid in substrates:
            node["substrate_classes"] = substrates[tcid]
        if superfam:
            node["superfamily"] = superfam
        h[tcid] = node

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(h, indent=2, sort_keys=True))
    log.info(f"  tcdb_hierarchy.json: {len(h)} entries")
    return len(h)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_tcdb_hierarchy -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.1B — TCDB hierarchy from families+superfamilies+substrates"
```

---

### Task 5: CAZy hierarchy bootstrap from raw eggNOG annotations

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (add `build_cazy_hierarchy`)
- Test: `tests/test_build_metabolite_resolver.py` (append)

- [ ] **Step 1: Append the failing test**

```python
# Append to tests/test_build_metabolite_resolver.py

EGGNOG_FIXTURE = textwrap.dedent("""\
    ##   eggnog-mapper                                                                                                                                                                                                                  
    #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
    WP_001.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGT19\t-\t-
    WP_002.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGH13_1,CBM48\t-\t-
    WP_003.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-
    WP_004.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGH32\t-\t-
""")


def test_build_cazy_hierarchy(tmp_path):
    """CAZy hierarchy bootstrapped from raw eggNOG col 19, all 6 classes hardcoded."""
    eggnog = tmp_path / "MED4.emapper.annotations"
    eggnog.write_text(EGGNOG_FIXTURE)
    out = tmp_path / "cazy_hierarchy.json"

    bmr.build_cazy_hierarchy(out, [eggnog])

    h = json.loads(out.read_text())

    # All 6 classes present even if not all observed
    for cls in ("GH", "GT", "PL", "CE", "AA", "CBM"):
        assert h[cls]["level"] == 0
        assert h[cls]["level_kind"] == "cazy_class"
        assert h[cls]["parent"] is None
        assert h[cls]["class"] == cls

    # Observed families derived
    assert h["GT19"]["level"] == 1
    assert h["GT19"]["parent"] == "GT"
    assert h["GT19"]["class"] == "GT"

    assert h["GH13"]["level"] == 1
    assert h["GH13"]["parent"] == "GH"

    # Observed subfamily derived (and its parent family also present)
    assert h["GH13_1"]["level"] == 2
    assert h["GH13_1"]["level_kind"] == "cazy_subfamily"
    assert h["GH13_1"]["parent"] == "GH13"
    assert h["GH13_1"]["class"] == "GH"

    assert h["CBM48"]["parent"] == "CBM"
    assert h["GH32"]["parent"] == "GH"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_cazy_hierarchy -v`
Expected: FAIL — `AttributeError: build_cazy_hierarchy`.

- [ ] **Step 3: Add the CAZy bootstrap function**

Insert in `build_metabolite_resolver.py` after `build_tcdb_hierarchy`:

```python
# ── CAZy hierarchy (bootstrapped from eggNOG observations) ────────────────────

_CAZY_CLASSES = {
    "GH":  "Glycoside Hydrolases",
    "GT":  "GlycosylTransferases",
    "PL":  "Polysaccharide Lyases",
    "CE":  "Carbohydrate Esterases",
    "AA":  "Auxiliary Activities",
    "CBM": "Carbohydrate-Binding Modules",
}

# eggNOG column index (1-based) for the CAZy field; 0-based here.
_EGGNOG_CAZY_COL = 18


# A CAZy ID is `<class><digits>` optionally followed by `_<digits>` subfamily
# Examples: GH13, GH13_1, CBM48, AA10
_CAZY_FAMILY_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)(?:_(\d+))?$")


def _parse_cazy_id(token: str) -> tuple[str, str | None] | None:
    """Return (family_id, subfamily_id_or_None) or None if malformed.

    'GH13'    → ('GH13', None)
    'GH13_1'  → ('GH13', 'GH13_1')
    'CBM48'   → ('CBM48', None)
    'invalid' → None
    """
    m = _CAZY_FAMILY_RE.match(token.strip())
    if not m:
        return None
    cls, fam_num, sub_num = m.groups()
    family = f"{cls}{fam_num}"
    subfamily = f"{family}_{sub_num}" if sub_num else None
    return family, subfamily


def _collect_cazy_ids(eggnog_paths: list[Path]) -> set[str]:
    """Iterate over all eggNOG annotation files; collect distinct CAZy tokens."""
    ids: set[str] = set()
    for path in eggnog_paths:
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= _EGGNOG_CAZY_COL:
                    continue
                cell = cols[_EGGNOG_CAZY_COL]
                if cell in ("", "-"):
                    continue
                for token in cell.split(","):
                    token = token.strip()
                    if token:
                        ids.add(token)
    return ids


def build_cazy_hierarchy(out_path: Path, eggnog_paths: list[Path]) -> int:
    """Build cazy_hierarchy.json from CAZy IDs observed in eggNOG output.

    Always emits the 6 top-level classes (level 0) so consumers have a stable
    set of root nodes. Family + subfamily entries are derived per observation.
    """
    h: dict[str, dict] = {
        cls: {
            "name": display_name,
            "level": 0,
            "level_kind": "cazy_class",
            "parent": None,
            "class": cls,
        }
        for cls, display_name in _CAZY_CLASSES.items()
    }

    observed_ids = _collect_cazy_ids(eggnog_paths)
    for token in sorted(observed_ids):
        parsed = _parse_cazy_id(token)
        if parsed is None:
            log.warning(f"  cazy: skipping malformed token: {token!r}")
            continue
        family, subfamily = parsed
        cls = _CAZY_FAMILY_RE.match(family).group(1)
        if family not in h:
            h[family] = {
                "name": "",
                "level": 1,
                "level_kind": "cazy_family",
                "parent": cls,
                "class": cls,
            }
        if subfamily and subfamily not in h:
            h[subfamily] = {
                "name": "",
                "level": 2,
                "level_kind": "cazy_subfamily",
                "parent": family,
                "class": cls,
            }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(h, indent=2, sort_keys=True))
    log.info(f"  cazy_hierarchy.json: {len(h)} entries ({len(observed_ids)} eggNOG observations)")
    return len(h)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_cazy_hierarchy -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py
git commit -m "metabolism: 1.1B — CAZy hierarchy bootstrap from eggNOG observations"
```

---

### Task 6: Wire `main()` end-to-end + diagnostic report

**Files:**
- Modify: `multiomics_kg/download/build_metabolite_resolver.py` (replace stub `main()`)
- Test: `tests/test_build_metabolite_resolver.py` (append integration test)

- [ ] **Step 1: Append the failing test**

```python
# Append to tests/test_build_metabolite_resolver.py

def test_build_main_end_to_end(tmp_path, monkeypatch):
    """main() wires all builders + writes diagnostic report. Synthetic cache."""
    cache = tmp_path / "cache" / "data"
    (cache / "mnx").mkdir(parents=True)
    (cache / "tcdb").mkdir(parents=True)

    # Synthesize all the input files
    (cache / "mnx" / "chem_prop.tsv").write_text(CHEM_PROP_FIXTURE)
    (cache / "mnx" / "chem_xref.tsv").write_text(CHEM_XREF_FIXTURE)
    (cache / "mnx" / "reac_prop.tsv").write_text(REAC_PROP_FIXTURE)
    (cache / "mnx" / "reac_xref.tsv").write_text(REAC_XREF_FIXTURE)
    (cache / "tcdb" / "families.tsv").write_text(TCDB_FAMILIES_FIXTURE)
    (cache / "tcdb" / "superfamilies.tsv").write_text(TCDB_SUPERFAMILIES_FIXTURE)
    (cache / "tcdb" / "substrates.tsv").write_text(TCDB_SUBSTRATES_FIXTURE)

    # Synthesize a tiny eggNOG annotation for one strain
    strain_dir = cache / "Prochlorococcus" / "genomes" / "MED4" / "eggnog"
    strain_dir.mkdir(parents=True)
    (strain_dir / "MED4.emapper.annotations").write_text(EGGNOG_FIXTURE)

    monkeypatch.chdir(tmp_path)
    bmr.main(force=True)

    # All four outputs exist
    assert (cache / "mnx" / "metabolite_resolver.db").exists()
    assert (cache / "tcdb" / "tcdb_hierarchy.json").exists()
    assert (cache / "cazy" / "cazy_hierarchy.json").exists()
    assert (cache / "mnx" / "metabolite_id_mapping_report.json").exists()

    # Diagnostic report has expected keys
    report = json.loads((cache / "mnx" / "metabolite_id_mapping_report.json").read_text())
    assert report["mnx_release"] == "MNXref 4.5 (2025-08-13)"
    assert report["compound_count"] == 3
    assert report["reaction_count"] == 3
    assert report["tcdb_hierarchy_entry_count"] >= 7
    assert report["cazy_hierarchy_entry_count"] >= 6  # at least the 6 classes
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_metabolite_resolver.py::test_build_main_end_to_end -v`
Expected: FAIL — currently `main()` raises `NotImplementedError`.

- [ ] **Step 3: Replace `main()` with the wired-up implementation**

In `multiomics_kg/download/build_metabolite_resolver.py`, **replace** the existing `main` function (and remove the stub `NotImplementedError`) with:

```python
# ── orchestration ─────────────────────────────────────────────────────────────

CACHE_ROOT = Path("cache/data")


def _resolve_paths(cache_root: Path) -> dict[str, Path]:
    return {
        "chem_prop":     cache_root / "mnx" / "chem_prop.tsv",
        "chem_xref":     cache_root / "mnx" / "chem_xref.tsv",
        "reac_prop":     cache_root / "mnx" / "reac_prop.tsv",
        "reac_xref":     cache_root / "mnx" / "reac_xref.tsv",
        "families":      cache_root / "tcdb" / "families.tsv",
        "superfamilies": cache_root / "tcdb" / "superfamilies.tsv",
        "substrates":    cache_root / "tcdb" / "substrates.tsv",
        "resolver_db":   cache_root / "mnx" / "metabolite_resolver.db",
        "tcdb_json":     cache_root / "tcdb" / "tcdb_hierarchy.json",
        "cazy_json":     cache_root / "cazy" / "cazy_hierarchy.json",
        "report":        cache_root / "mnx" / "metabolite_id_mapping_report.json",
    }


def _read_mnx_release(chem_prop_path: Path) -> str:
    """Extract '#VERSION' and '#DATE' from chem_prop.tsv comment header."""
    version = ""
    date = ""
    with open(chem_prop_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if "VERSION" in line:
                version = line.split(":", 1)[1].strip().rstrip("\t")
            elif "DATE" in line:
                date = line.split(":", 1)[1].strip().rstrip("\t")
            if version and date:
                break
    return f"MNXref {version} ({date})" if version and date else "MNXref unknown"


def _find_eggnog_annotations(cache_root: Path) -> list[Path]:
    """Discover all <strain>.emapper.annotations files under cache_root."""
    return sorted(cache_root.glob("*/genomes/*/eggnog/*.emapper.annotations"))


def _alias_count_by_source(conn: sqlite3.Connection) -> dict[str, int]:
    cur = conn.cursor()
    cur.execute("SELECT source, COUNT(*) FROM compound_aliases GROUP BY source ORDER BY 2 DESC")
    return dict(cur.fetchall())


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    paths = _resolve_paths(CACHE_ROOT)

    if not force and paths["resolver_db"].exists() and paths["tcdb_json"].exists() and paths["cazy_json"].exists():
        log.info("All caches present; skip (use --force to rebuild).")
        return

    log.info("Building metabolite resolver + hierarchy caches")

    # SQLite resolver
    paths["resolver_db"].parent.mkdir(parents=True, exist_ok=True)
    if paths["resolver_db"].exists():
        paths["resolver_db"].unlink()
    conn = sqlite3.connect(paths["resolver_db"])
    n_compounds = build_compounds_table(conn, paths["chem_prop"])
    build_compound_aliases_table(conn, paths["chem_xref"])
    build_compound_names_table(conn, paths["chem_xref"])
    n_reactions = build_reactions_table(conn, paths["reac_prop"])
    build_reaction_aliases_table(conn, paths["reac_xref"])
    conn.commit()

    alias_counts = _alias_count_by_source(conn)
    conn.close()

    # TCDB hierarchy
    n_tcdb = build_tcdb_hierarchy(
        paths["tcdb_json"], paths["families"], paths["superfamilies"], paths["substrates"]
    )

    # CAZy hierarchy (bootstrapped from eggNOG observations)
    eggnog_paths = _find_eggnog_annotations(CACHE_ROOT)
    n_cazy = build_cazy_hierarchy(paths["cazy_json"], eggnog_paths)

    # Diagnostic report
    report = {
        "mnx_release":              _read_mnx_release(paths["chem_prop"]),
        "compound_count":           n_compounds,
        "reaction_count":           n_reactions,
        "alias_counts_by_source":   alias_counts,
        "tcdb_hierarchy_entry_count": n_tcdb,
        "cazy_hierarchy_entry_count": n_cazy,
        "eggnog_annotation_files_scanned": len(eggnog_paths),
    }
    paths["report"].write_text(json.dumps(report, indent=2))
    log.info(f"  metabolite_id_mapping_report.json written ({len(report)} fields)")
    log.info("Done.")
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `uv run pytest tests/test_build_metabolite_resolver.py -v`
Expected: 6 tests PASS.

Also re-run the existing skeleton tests to make sure none of them broke:
Run: `uv run pytest tests/test_metabolism_skeletons.py::test_build_metabolite_resolver_imports_and_raises -v`
Expected: this test will now FAIL because `main()` no longer raises `NotImplementedError`. **Delete this test** — it's superseded by the real one.

```python
# Edit tests/test_metabolism_skeletons.py — remove the entire
# test_build_metabolite_resolver_imports_and_raises function (and its helpers
# if exclusively used by it). Leave the other 3 skeleton tests for now;
# they'll be deleted as their respective accessors get real impls in tasks 8-10.
```

- [ ] **Step 5: Run all tests once more to confirm clean**

Run: `uv run pytest tests/test_metabolism_skeletons.py tests/test_build_metabolite_resolver.py -v`
Expected: 9 tests PASS (3 remaining skeleton + 6 build).

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_metabolite_resolver.py tests/test_build_metabolite_resolver.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: 1.1B — wire main() + diagnostic report (sub-step 7 done)"
```

---

### Task 7: Implement `metabolite_utils` accessor

**Files:**
- Modify: `multiomics_kg/utils/metabolite_utils.py` (replace skeleton)
- Test: `tests/test_metabolite_utils.py` (new)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_metabolite_utils.py
"""Tests for utils/metabolite_utils.py — resolver accessors."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.utils import metabolite_utils as mu


@pytest.fixture
def tiny_resolver(tmp_path) -> Path:
    """Build a tiny resolver DB inline: 2 compounds, 4 aliases, 3 names; 2 reactions, 3 aliases."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY (source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY (name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY (source, value, mnxr_id));
    """)
    cur.execute("INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', '', '', 0, 0, '', '', '')")
    cur.execute("INSERT INTO compounds VALUES ('MNXM999', 'ambig-name', '', '', 0, 0, '', '', '')")
    cur.executemany("INSERT INTO compound_aliases VALUES (?, ?, ?)", [
        ("chebi", "17234", "MNXM41"),
        ("kegg.compound", "C00031", "MNXM41"),
        ("ambiguous-x", "X", "MNXM41"),
        ("ambiguous-x", "X", "MNXM999"),
    ])
    cur.executemany("INSERT INTO compound_names VALUES (?, ?)", [
        ("d-glucose", "MNXM41"),
        ("dextrose", "MNXM41"),
        ("ambig-name", "MNXM999"),
    ])
    cur.execute("INSERT INTO reactions VALUES ('MNXR101234', '1 MNXM3 = 1 MNXM41', '', '', 'B', NULL)")
    cur.executemany("INSERT INTO reaction_aliases VALUES (?, ?, ?)", [
        ("kegg.reaction", "R00299", "MNXR101234"),
        ("rhea", "16332", "MNXR101234"),
        ("rhea", "DUP", "MNXR101234"),
    ])
    conn.commit()
    conn.close()
    return db


def test_open_resolver_reads_db(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM compounds")
    assert cur.fetchone()[0] == 2
    conn.close()


def test_resolve_metabolite_direct_mnxm(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("MNXM41", conn) == ("MNXM41", "xref:exact")


def test_resolve_metabolite_alias(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("17234", conn) == ("MNXM41", "xref:exact")
    assert mu.resolve_metabolite("C00031", conn) == ("MNXM41", "xref:exact")


def test_resolve_metabolite_alias_ambiguous(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    mnxm, method = mu.resolve_metabolite("X", conn)
    assert mnxm in ("MNXM41", "MNXM999")
    assert method == "xref:ambiguous"


def test_resolve_metabolite_name(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("D-Glucose", conn) == ("MNXM41", "name:normalized")
    # Whitespace and case variation also resolves
    assert mu.resolve_metabolite("  Dextrose  ", conn) == ("MNXM41", "name:normalized")


def test_resolve_metabolite_unresolved(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("does-not-exist", conn) == (None, "unresolved")


def test_resolve_reaction_direct(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("MNXR101234", conn) == ("MNXR101234", "xref:exact")


def test_resolve_reaction_alias(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("R00299", conn) == ("MNXR101234", "xref:exact")
    assert mu.resolve_reaction("16332", conn) == ("MNXR101234", "xref:exact")


def test_resolve_reaction_unresolved(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("R99999", conn) == (None, "unresolved")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_metabolite_utils.py -v`
Expected: FAIL — `NotImplementedError: Phase 1.1B`.

- [ ] **Step 3: Replace `metabolite_utils.py` skeleton with real impl**

Replace the entire contents of `multiomics_kg/utils/metabolite_utils.py` with:

```python
"""Resolver accessor API.

Reads the SQLite resolver DB built by sub-step 7. Used by step 2 transforms,
the Spec 1.2 scaffold builder, and Phase 2 paper-measurement extraction.
"""
from __future__ import annotations

import re
import sqlite3
from pathlib import Path

DEFAULT_DB_PATH = Path("cache/data/mnx/metabolite_resolver.db")

_NAME_NORM_RE = re.compile(r"[^\w+\-/ ]")


def _normalize_name(s: str) -> str:
    """Lowercase, collapse whitespace, strip punctuation except + - /."""
    s = _NAME_NORM_RE.sub(" ", s.lower())
    return " ".join(s.split())


def open_resolver(path: Path | None = None) -> sqlite3.Connection:
    """Open a read-only connection to the resolver DB."""
    p = Path(path) if path else DEFAULT_DB_PATH
    return sqlite3.connect(f"file:{p}?mode=ro", uri=True)


def resolve_metabolite(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """Resolve a value to an MNXM ID. See spec for method enum."""
    value = value.strip()
    cur = conn.cursor()

    # 1. Direct MNXM
    if value.startswith("MNXM"):
        cur.execute("SELECT 1 FROM compounds WHERE mnxm_id = ?", (value,))
        if cur.fetchone():
            return value, "xref:exact"

    # 2. Alias match (across all sources)
    cur.execute("SELECT mnxm_id FROM compound_aliases WHERE value = ? ORDER BY mnxm_id", (value,))
    rows = cur.fetchall()
    if len(rows) == 1:
        return rows[0][0], "xref:exact"
    if len(rows) > 1:
        return rows[0][0], "xref:ambiguous"

    # 3. Name match (after normalization)
    normalized = _normalize_name(value)
    if normalized:
        cur.execute("SELECT mnxm_id FROM compound_names WHERE name_normalized = ? ORDER BY mnxm_id", (normalized,))
        rows = cur.fetchall()
        if len(rows) == 1:
            return rows[0][0], "name:normalized"
        if len(rows) > 1:
            return rows[0][0], "ambiguous"

    return None, "unresolved"


def resolve_reaction(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """Resolve a value to an MNXR ID. Alias-only (no name normalization)."""
    value = value.strip()
    cur = conn.cursor()

    if value.startswith("MNXR"):
        cur.execute("SELECT 1 FROM reactions WHERE mnxr_id = ?", (value,))
        if cur.fetchone():
            return value, "xref:exact"

    cur.execute("SELECT mnxr_id FROM reaction_aliases WHERE value = ? ORDER BY mnxr_id", (value,))
    rows = cur.fetchall()
    if len(rows) == 1:
        return rows[0][0], "xref:exact"
    if len(rows) > 1:
        return rows[0][0], "xref:ambiguous"

    return None, "unresolved"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_metabolite_utils.py -v`
Expected: 9 PASS.

- [ ] **Step 5: Remove the now-obsolete skeleton test**

Delete `test_metabolite_utils_imports_and_raises` from `tests/test_metabolism_skeletons.py` (it now fails because the module no longer raises `NotImplementedError`).

Run: `uv run pytest tests/test_metabolism_skeletons.py -v`
Expected: 2 tests remain (`tcdb_utils` + `cazy_utils` skeletons), both PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/metabolite_utils.py tests/test_metabolite_utils.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: 1.1B — metabolite_utils real accessor impl"
```

---

### Task 8: Implement `tcdb_utils` accessor

**Files:**
- Modify: `multiomics_kg/utils/tcdb_utils.py`
- Test: `tests/test_tcdb_utils.py` (new)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_tcdb_utils.py
"""Tests for utils/tcdb_utils.py — TCDB hierarchy accessor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.utils import tcdb_utils as tu


TINY_TCDB = {
    "1":         {"name": "Channels and Pores", "level": 0, "level_kind": "tc_class", "parent": None},
    "1.A":       {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
    "1.A.1":     {"name": "VIC Family", "level": 2, "level_kind": "tc_family", "parent": "1.A", "abbreviation": "VIC"},
    "1.A.1.1":   {"name": "", "level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.1.1": {"name": "", "level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.1",
                  "substrate_classes": ["potassium(1+)"]},
}


@pytest.fixture(autouse=True)
def patch_tcdb(monkeypatch, tmp_path):
    """Point load_tcdb at a tiny synthetic JSON; reset the module cache."""
    p = tmp_path / "tcdb_hierarchy.json"
    p.write_text(json.dumps(TINY_TCDB))
    monkeypatch.setattr(tu, "DEFAULT_PATH", p)
    monkeypatch.setattr(tu, "_CACHE", None)
    yield


def test_load_tcdb_returns_dict():
    h = tu.load_tcdb()
    assert isinstance(h, dict)
    assert h["1"]["level"] == 0


def test_load_tcdb_caches():
    """Second call returns the same object without re-reading the file."""
    a = tu.load_tcdb()
    b = tu.load_tcdb()
    assert a is b


def test_is_valid_tcdb_present():
    assert tu.is_valid_tcdb("1.A.1.1.1") is True
    assert tu.is_valid_tcdb("1.A") is True
    assert tu.is_valid_tcdb("1") is True


def test_is_valid_tcdb_absent():
    assert tu.is_valid_tcdb("99.X.99") is False
    assert tu.is_valid_tcdb("") is False


def test_tcdb_ancestors_full_chain():
    """Specificity-level TCID → list of all ancestors (root → parent)."""
    assert tu.tcdb_ancestors("1.A.1.1.1") == ["1", "1.A", "1.A.1", "1.A.1.1"]


def test_tcdb_ancestors_partial():
    assert tu.tcdb_ancestors("1.A.1") == ["1", "1.A"]
    assert tu.tcdb_ancestors("1") == []


def test_tcdb_ancestors_unknown_returns_empty():
    assert tu.tcdb_ancestors("99.X.99") == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_tcdb_utils.py -v`
Expected: FAIL — `NotImplementedError: Phase 1.1B`.

- [ ] **Step 3: Replace `tcdb_utils.py` skeleton with real impl**

Replace the entire contents of `multiomics_kg/utils/tcdb_utils.py` with:

```python
"""TCDB hierarchy accessor API.

Reads cache/data/tcdb/tcdb_hierarchy.json (built by sub-step 7).
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT_PATH = Path("cache/data/tcdb/tcdb_hierarchy.json")
_CACHE: dict[str, dict] | None = None


def load_tcdb() -> dict[str, dict]:
    """Load the TCDB hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_tcdb(tc_id: str) -> bool:
    """Return True if `tc_id` exists as a key in the hierarchy."""
    if not tc_id:
        return False
    return tc_id in load_tcdb()


def tcdb_ancestors(tc_id: str) -> list[str]:
    """Return root-to-parent ancestor chain for `tc_id`. Empty for unknown / root."""
    h = load_tcdb()
    if tc_id not in h:
        return []
    chain: list[str] = []
    parent = h[tc_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_tcdb_utils.py -v`
Expected: 7 tests PASS.

- [ ] **Step 5: Remove the now-obsolete skeleton test**

Delete `test_tcdb_utils_imports_and_raises` from `tests/test_metabolism_skeletons.py`.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/tcdb_utils.py tests/test_tcdb_utils.py tests/test_metabolism_skeletons.py
git commit -m "metabolism: 1.1B — tcdb_utils real accessor impl"
```

---

### Task 9: Implement `cazy_utils` accessor

**Files:**
- Modify: `multiomics_kg/utils/cazy_utils.py`
- Test: `tests/test_cazy_utils.py` (new)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cazy_utils.py
"""Tests for utils/cazy_utils.py — CAZy hierarchy accessor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.utils import cazy_utils as cu


TINY_CAZY = {
    "GH":     {"name": "Glycoside Hydrolases", "level": 0, "level_kind": "cazy_class",     "parent": None,    "class": "GH"},
    "GT":     {"name": "GlycosylTransferases", "level": 0, "level_kind": "cazy_class",     "parent": None,    "class": "GT"},
    "GH13":   {"name": "",                     "level": 1, "level_kind": "cazy_family",    "parent": "GH",    "class": "GH"},
    "GH13_1": {"name": "",                     "level": 2, "level_kind": "cazy_subfamily", "parent": "GH13",  "class": "GH"},
    "GT19":   {"name": "",                     "level": 1, "level_kind": "cazy_family",    "parent": "GT",    "class": "GT"},
}


@pytest.fixture(autouse=True)
def patch_cazy(monkeypatch, tmp_path):
    p = tmp_path / "cazy_hierarchy.json"
    p.write_text(json.dumps(TINY_CAZY))
    monkeypatch.setattr(cu, "DEFAULT_PATH", p)
    monkeypatch.setattr(cu, "_CACHE", None)
    yield


def test_load_cazy_returns_dict():
    h = cu.load_cazy()
    assert isinstance(h, dict)
    assert h["GH"]["class"] == "GH"


def test_is_valid_cazy_present():
    assert cu.is_valid_cazy("GH") is True
    assert cu.is_valid_cazy("GH13") is True
    assert cu.is_valid_cazy("GH13_1") is True


def test_is_valid_cazy_absent():
    assert cu.is_valid_cazy("XX99") is False
    assert cu.is_valid_cazy("GH99999") is False  # plausible format but not in hierarchy
    assert cu.is_valid_cazy("") is False


def test_cazy_ancestors_subfamily():
    assert cu.cazy_ancestors("GH13_1") == ["GH", "GH13"]


def test_cazy_ancestors_family():
    assert cu.cazy_ancestors("GH13") == ["GH"]


def test_cazy_ancestors_class():
    assert cu.cazy_ancestors("GH") == []


def test_cazy_ancestors_unknown():
    assert cu.cazy_ancestors("XX99") == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cazy_utils.py -v`
Expected: FAIL — `NotImplementedError: Phase 1.1B`.

- [ ] **Step 3: Replace `cazy_utils.py` skeleton with real impl**

Replace the entire contents of `multiomics_kg/utils/cazy_utils.py` with:

```python
"""CAZy hierarchy accessor API.

Reads cache/data/cazy/cazy_hierarchy.json (built by sub-step 7 from observed
eggNOG `CAZy` columns).
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT_PATH = Path("cache/data/cazy/cazy_hierarchy.json")
_CACHE: dict[str, dict] | None = None


def load_cazy() -> dict[str, dict]:
    """Load the CAZy hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_cazy(cazy_id: str) -> bool:
    if not cazy_id:
        return False
    return cazy_id in load_cazy()


def cazy_ancestors(cazy_id: str) -> list[str]:
    h = load_cazy()
    if cazy_id not in h:
        return []
    chain: list[str] = []
    parent = h[cazy_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_cazy_utils.py -v`
Expected: 7 tests PASS.

- [ ] **Step 5: Remove the obsolete skeleton test**

Delete `test_cazy_utils_imports_and_raises` from `tests/test_metabolism_skeletons.py`. The file should now contain **zero** tests (or only `import` boilerplate). If empty, delete the file:

```bash
test -s tests/test_metabolism_skeletons.py || rm tests/test_metabolism_skeletons.py
```

Actually, simpler: after the deletion, the file has only the module docstring + imports. Just delete the whole file with `git rm`:

```bash
git rm tests/test_metabolism_skeletons.py
```

- [ ] **Step 6: Run full test sweep**

Run: `uv run pytest tests/test_metabolite_utils.py tests/test_tcdb_utils.py tests/test_cazy_utils.py tests/test_build_metabolite_resolver.py -v`
Expected: 9 + 7 + 7 + 6 = 29 tests PASS.

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/utils/cazy_utils.py tests/test_cazy_utils.py
git rm tests/test_metabolism_skeletons.py
git commit -m "metabolism: 1.1B — cazy_utils real accessor impl + drop skeleton tests"
```

---

### Task 10: Three transforms in `annotation_transforms.py`

**Files:**
- Modify: `multiomics_kg/download/utils/annotation_transforms.py` (add 3 `_tx_*` functions + accessor cache)
- Test: `tests/test_annotation_transforms_metabolism.py` (new)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_annotation_transforms_metabolism.py
"""Tests for the three Phase 1.1B metabolism transforms."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download.utils import annotation_transforms as at


@pytest.fixture(autouse=True)
def patch_caches(monkeypatch, tmp_path):
    """Set up tiny resolver DB + tiny TCDB/CAZy hierarchies, redirect module caches."""
    # Resolver DB
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY(source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY(name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY(source, value, mnxr_id));
    """)
    conn.execute("INSERT INTO reactions VALUES ('MNXR101234', '', '', '', 'B', NULL)")
    conn.execute("INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00299', 'MNXR101234')")
    conn.commit()
    conn.close()

    from multiomics_kg.utils import metabolite_utils as mu
    monkeypatch.setattr(mu, "DEFAULT_DB_PATH", db)
    monkeypatch.setattr(at, "_RESOLVER_CONN", None)

    # TCDB hierarchy
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1": {}, "3.A.1.1": {}, "3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    # CAZy hierarchy
    from multiomics_kg.utils import cazy_utils as cu
    cazy_path = tmp_path / "cazy_hierarchy.json"
    cazy_path.write_text(json.dumps({"GH": {}, "GH13": {}, "GH13_1": {}, "GT": {}, "GT19": {}}))
    monkeypatch.setattr(cu, "DEFAULT_PATH", cazy_path)
    monkeypatch.setattr(cu, "_CACHE", None)

    yield


def test_resolve_kegg_reaction_resolves_known():
    assert at._tx_resolve_kegg_reaction_to_mnxr("R00299") == "MNXR101234"


def test_resolve_kegg_reaction_drops_unknown():
    assert at._tx_resolve_kegg_reaction_to_mnxr("R99999") is None


def test_validate_tcdb_keeps_known():
    assert at._tx_validate_tcdb("3.A.1.1.1") == "3.A.1.1.1"
    assert at._tx_validate_tcdb("3.A.1") == "3.A.1"


def test_validate_tcdb_drops_unknown():
    assert at._tx_validate_tcdb("99.X.99") is None


def test_validate_cazy_keeps_known():
    assert at._tx_validate_cazy("GH13_1") == "GH13_1"
    assert at._tx_validate_cazy("GT19") == "GT19"


def test_validate_cazy_drops_unknown():
    assert at._tx_validate_cazy("XX99") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_annotation_transforms_metabolism.py -v`
Expected: FAIL — `AttributeError: _tx_resolve_kegg_reaction_to_mnxr` etc.

- [ ] **Step 3: Add the three transforms + cache machinery**

Append to `multiomics_kg/download/utils/annotation_transforms.py`:

```python
# ── Phase 1.1B metabolism transforms ──────────────────────────────────────────

from multiomics_kg.utils.metabolite_utils import open_resolver, resolve_reaction
from multiomics_kg.utils.tcdb_utils import is_valid_tcdb
from multiomics_kg.utils.cazy_utils import is_valid_cazy

_RESOLVER_CONN = None


def _get_resolver_conn():
    """Lazy open + cache module-level resolver connection."""
    global _RESOLVER_CONN
    if _RESOLVER_CONN is None:
        _RESOLVER_CONN = open_resolver()
    return _RESOLVER_CONN


def _tx_resolve_kegg_reaction_to_mnxr(value: str) -> str | None:
    """KEGG R-number → MNXR ID. Returns None if unresolved (drops the token)."""
    mnxr, _method = resolve_reaction(value, _get_resolver_conn())
    return mnxr


def _tx_validate_tcdb(value: str) -> str | None:
    """Drop TC IDs not present in the TCDB hierarchy."""
    return value if is_valid_tcdb(value) else None


def _tx_validate_cazy(value: str) -> str | None:
    """Drop CAZy IDs not present in the CAZy hierarchy."""
    return value if is_valid_cazy(value) else None
```

Now register the three transforms in the framework's transform-name → callable lookup. Find the existing dispatch table — likely a dict at the bottom of `annotation_transforms.py` or a `TRANSFORMS = {...}` constant. Read the surrounding code:

```bash
grep -n "TRANSFORMS\|first_token_space\|add_go_prefix" multiomics_kg/download/utils/annotation_transforms.py
```

Find the dict that maps existing transform names (`first_token_space`, `add_go_prefix`, `strip_function_prefix`, etc.) to functions, and add three new entries:

```python
    "resolve_kegg_reaction_to_mnxr": _tx_resolve_kegg_reaction_to_mnxr,
    "validate_tcdb":                  _tx_validate_tcdb,
    "validate_cazy":                  _tx_validate_cazy,
```

(The exact key style follows the existing convention in the file — use the same naming pattern as `first_token_space`, `add_go_prefix`, etc.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_annotation_transforms_metabolism.py -v`
Expected: 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/utils/annotation_transforms.py tests/test_annotation_transforms_metabolism.py
git commit -m "metabolism: 1.1B — three transforms (resolve_kegg_reaction, validate_tcdb, validate_cazy)"
```

---

### Task 11: Wire YAML config + per-strain metabolism report

**Files:**
- Modify: `config/gene_annotations_config.yaml` (rename + add `transforms:` lines on 3 fields)
- Modify: `multiomics_kg/download/build_gene_annotations.py` (invoke per-strain report at end of build)
- Test: `tests/test_build_gene_annotations_metabolism.py` (new — end-to-end YAML transform integration)

- [ ] **Step 1: Read the existing YAML config block**

```bash
grep -n "^  kegg_reaction:\|^  transporter_classification:\|^  cazy_ids:" config/gene_annotations_config.yaml
```

Note the line numbers of the three fields. They each have shape:

```yaml
field_name:
  type: union
  sources:
    - source: eggnog
      field: <UPPER_CASE_NAME>
      delimiter: ","
```

- [ ] **Step 2: Edit the YAML — rename + add transforms**

For `kegg_reaction` (rename + transform): change

```yaml
  kegg_reaction:
    type: union
    sources:
      - source: eggnog
        field: KEGG_Reaction
        delimiter: ","
```

to

```yaml
  kegg_reactions:
    type: union
    transforms: [resolve_kegg_reaction_to_mnxr]
    sources:
      - source: eggnog
        field: KEGG_Reaction
        delimiter: ","
```

For `transporter_classification` (just add transform):

```yaml
  transporter_classification:
    type: union
    transforms: [validate_tcdb]
    sources:
      - source: eggnog
        field: KEGG_TC
        delimiter: ","
```

For `cazy_ids` (just add transform):

```yaml
  cazy_ids:
    type: union
    transforms: [validate_cazy]
    sources:
      - source: eggnog
        field: CAZy
        delimiter: ","
```

- [ ] **Step 3: Write the failing integration test**

```python
# tests/test_build_gene_annotations_metabolism.py
"""End-to-end test: YAML transform pipeline produces resolved fields."""
from __future__ import annotations

import json
import sqlite3
import textwrap
from pathlib import Path

import pytest


EGGNOG_LINES = textwrap.dedent("""\
    ##
    #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
    WP_001.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\tR00299\t-\t-\t3.A.1.1.1,99.X.99\tGH13_1,XX99\t-\t-
""")


@pytest.fixture(autouse=True)
def patch_metabolism_caches(monkeypatch, tmp_path):
    """Tiny resolver/tcdb/cazy fixtures matching the eggNOG test data above."""
    # Resolver: R00299 → MNXR101234
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY(source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY(name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY(source, value, mnxr_id));
    """)
    conn.execute("INSERT INTO reactions VALUES ('MNXR101234', '', '', '', 'B', NULL)")
    conn.execute("INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00299', 'MNXR101234')")
    conn.commit()
    conn.close()

    from multiomics_kg.utils import metabolite_utils as mu
    monkeypatch.setattr(mu, "DEFAULT_DB_PATH", db)

    from multiomics_kg.download.utils import annotation_transforms as at
    monkeypatch.setattr(at, "_RESOLVER_CONN", None)

    # TCDB: 3.A.1.1.1 valid, 99.X.99 invalid
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    # CAZy: GH13_1 valid (and GH13 / GH for closure), XX99 invalid
    from multiomics_kg.utils import cazy_utils as cu
    cazy_path = tmp_path / "cazy_hierarchy.json"
    cazy_path.write_text(json.dumps({"GH": {}, "GH13": {}, "GH13_1": {}}))
    monkeypatch.setattr(cu, "DEFAULT_PATH", cazy_path)
    monkeypatch.setattr(cu, "_CACHE", None)

    yield


def test_yaml_transforms_produce_resolved_fields(tmp_path, monkeypatch):
    """Run the YAML pipeline against a tiny eggNOG file and assert merged fields.

    Calls the real process_strain(row, config) entry point.
    """
    import yaml
    from multiomics_kg.download import build_gene_annotations as bga

    # Synthesize the data_dir layout that process_strain expects:
    #   <data_dir>/eggnog/<strain>.emapper.annotations
    #   <data_dir>/gene_mapping.csv
    data_dir = tmp_path / "MED4_data"
    eggnog_dir = data_dir / "eggnog"
    eggnog_dir.mkdir(parents=True)
    (eggnog_dir / "MED4.emapper.annotations").write_text(EGGNOG_LINES)
    (data_dir / "gene_mapping.csv").write_text("locus_tag,protein_id\nPMM0001,WP_001.1\n")

    # Load the real YAML config from the project tree
    project_root = Path(__file__).parent.parent
    config = yaml.safe_load((project_root / "config/gene_annotations_config.yaml").read_text())

    row = {
        "strain_name":     "MED4",
        "preferred_name":  "Prochlorococcus MED4",
        "data_dir":        str(data_dir),
        "ncbi_taxon_id":   "59919",
    }
    bga.process_strain(row, config, force=True, pfam_data=None)

    merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
    gene = merged["PMM0001"]

    assert gene["kegg_reactions"] == ["MNXR101234"]
    assert gene["transporter_classification"] == ["3.A.1.1.1"]
    assert gene["cazy_ids"] == ["GH13_1"]
```

- [ ] **Step 4: Run the test**

Run: `uv run pytest tests/test_build_gene_annotations_metabolism.py -v`
Expected: PASS if `merge_for_strain()` exists and the framework is wired correctly. If `merge_for_strain` does not exist by exactly that name, the test will FAIL with `AttributeError`. **Read** `multiomics_kg/download/build_gene_annotations.py` to find the correct entry point — it might be `build_strain` or a different name. Adjust the test's call accordingly. (The function is the strain-level merge that returns the dict mapping locus_tag → fields.)

- [ ] **Step 5: Run the existing skeleton + sub-step-7 tests to confirm no regression**

Run: `uv run pytest tests/test_metabolite_utils.py tests/test_tcdb_utils.py tests/test_cazy_utils.py tests/test_annotation_transforms_metabolism.py tests/test_build_gene_annotations_metabolism.py tests/test_build_metabolite_resolver.py tests/test_download_metabolism_reference.py -v`
Expected: 9 + 7 + 7 + 6 + 1 + 6 + 3 = 39 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add config/gene_annotations_config.yaml tests/test_build_gene_annotations_metabolism.py
git commit -m "metabolism: 1.1B — wire 3 transforms into YAML config + integration test"
```

---

### Task 12: Per-strain metabolism report

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py` (add post-merge sweep that writes the report)
- Test: `tests/test_build_gene_annotations_metabolism.py` (append report test)

- [ ] **Step 1: Append the failing test**

```python
# Append to tests/test_build_gene_annotations_metabolism.py

def test_per_strain_metabolism_report_written(tmp_path, monkeypatch):
    """After process_strain runs, step2_metabolism_report.json sits next to gene_annotations_merged.json."""
    import yaml
    from multiomics_kg.download import build_gene_annotations as bga

    data_dir = tmp_path / "MED4_data"
    eggnog_dir = data_dir / "eggnog"
    eggnog_dir.mkdir(parents=True)
    (eggnog_dir / "MED4.emapper.annotations").write_text(EGGNOG_LINES)
    (data_dir / "gene_mapping.csv").write_text("locus_tag,protein_id\nPMM0001,WP_001.1\n")

    project_root = Path(__file__).parent.parent
    config = yaml.safe_load((project_root / "config/gene_annotations_config.yaml").read_text())

    row = {
        "strain_name":     "MED4",
        "preferred_name":  "Prochlorococcus MED4",
        "data_dir":        str(data_dir),
        "ncbi_taxon_id":   "59919",
    }
    bga.process_strain(row, config, force=True, pfam_data=None)

    report_path = data_dir / "step2_metabolism_report.json"
    assert report_path.exists()
    report = json.loads(report_path.read_text())

    assert report["strain"] == "MED4"
    assert report["gene_count"] == 1

    kr = report["kegg_reactions"]
    assert kr["raw_total"] == 1
    assert kr["resolved_total"] == 1
    assert kr["resolved_unique_mnxr"] == 1
    assert kr["unresolved_unique"] == 0

    tc = report["transporter_classification"]
    assert tc["raw_total"] == 2
    assert tc["validated_total"] == 1
    assert "99.X.99" in tc["invalid_examples"]

    cz = report["cazy_ids"]
    assert cz["raw_total"] == 2
    assert cz["validated_total"] == 1
    assert "XX99" in cz["invalid_examples"]
```

- [ ] **Step 2: Run test**

Run: `uv run pytest tests/test_build_gene_annotations_metabolism.py::test_per_strain_metabolism_report_written -v`
Expected: FAIL (the report writer function doesn't exist yet).

- [ ] **Step 3: Implement the report writer**

In `multiomics_kg/download/build_gene_annotations.py`, the strain-level entry point is `process_strain(row, config, force, pfam_data)` (around line 773). It already writes `gene_annotations_merged.json` near the end. Add a call to a new helper just after the merged-JSON write:

```python
# At the end of process_strain, after the merged JSON is written:
_write_metabolism_report(
    strain_name=strain_name,
    data_dir=Path(data_dir),
    raw_eggnog_path=Path(data_dir) / "eggnog" / f"{strain_name}.emapper.annotations",
    merged=result,  # the dict that was just serialized
)
```

Helper definition (add to the module, near the other private `_*` helpers):

```python
def _write_metabolism_report(strain_name: str, data_dir: Path,
                              raw_eggnog_path: Path,
                              merged: dict) -> None:
    """Generate step2_metabolism_report.json for one strain.

    Counts raw values (from eggNOG TSV cols 15/18/19) vs resolved/validated
    values (from the merged JSON's kegg_reactions / transporter_classification /
    cazy_ids fields).
    """
    import json
    from collections import Counter

    raw_kegg: list[str] = []
    raw_tc: list[str] = []
    raw_cazy: list[str] = []

    if raw_eggnog_path.exists():
        with open(raw_eggnog_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= 18:
                    continue
                for col_idx, target in ((14, raw_kegg), (17, raw_tc), (18, raw_cazy)):
                    cell = cols[col_idx] if col_idx < len(cols) else ""
                    if cell in ("", "-"):
                        continue
                    for tok in cell.split(","):
                        tok = tok.strip()
                        if tok:
                            target.append(tok)

    resolved_mnxr: set[str] = set()
    validated_tc: set[str] = set()
    validated_cazy: set[str] = set()
    for gene in merged.values():
        resolved_mnxr.update(gene.get("kegg_reactions", []))
        validated_tc.update(gene.get("transporter_classification", []))
        validated_cazy.update(gene.get("cazy_ids", []))

    invalid_tc = sorted(set(raw_tc) - validated_tc)[:10]
    invalid_cazy = sorted(set(raw_cazy) - validated_cazy)[:10]

    report = {
        "strain": strain_name,
        "gene_count": len(merged),
        "kegg_reactions": {
            "raw_total":             len(raw_kegg),
            "raw_unique":            len(set(raw_kegg)),
            "resolved_total":        sum(len(g.get("kegg_reactions", [])) for g in merged.values()),
            "resolved_unique_mnxr":  len(resolved_mnxr),
            "unresolved_unique":     len(set(raw_kegg)) - len(resolved_mnxr),
        },
        "transporter_classification": {
            "raw_total":         len(raw_tc),
            "raw_unique":        len(set(raw_tc)),
            "validated_total":   sum(len(g.get("transporter_classification", [])) for g in merged.values()),
            "validated_unique":  len(validated_tc),
            "invalid_examples":  invalid_tc,
        },
        "cazy_ids": {
            "raw_total":         len(raw_cazy),
            "raw_unique":        len(set(raw_cazy)),
            "validated_total":   sum(len(g.get("cazy_ids", [])) for g in merged.values()),
            "validated_unique":  len(validated_cazy),
            "invalid_examples":  invalid_cazy,
        },
    }

    out_path = data_dir / "step2_metabolism_report.json"
    out_path.write_text(json.dumps(report, indent=2))
```

Note: `merged` here is the `dict[locus_tag → fields]` that the caller already has in scope from the merge logic. Look at the existing code's variable name (likely `result` or `merged_data`) and pass it through. The helper does not re-read `gene_annotations_merged.json`.

- [ ] **Step 4: Run test**

Run: `uv run pytest tests/test_build_gene_annotations_metabolism.py -v`
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_gene_annotations.py tests/test_build_gene_annotations_metabolism.py
git commit -m "metabolism: 1.1B — per-strain metabolism report"
```

---

### Task 13: Integration smoke test

**Files:**
- Test: `tests/test_prepare_data_step2_metabolism_smoke.py` (new)

- [ ] **Step 1: Write the smoke test**

```python
# tests/test_prepare_data_step2_metabolism_smoke.py
"""End-to-end smoke test for Phase 1.1B against real cache files.

Assumes:
- cache/data/mnx/, cache/data/tcdb/ are populated (sub-step 6 has run)
- cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations exists
- cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv exists

If any of those is missing, the test is skipped.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


REQUIRED_FIXTURES = [
    "cache/data/mnx/chem_prop.tsv",
    "cache/data/mnx/chem_xref.tsv",
    "cache/data/mnx/reac_prop.tsv",
    "cache/data/mnx/reac_xref.tsv",
    "cache/data/tcdb/families.tsv",
    "cache/data/tcdb/superfamilies.tsv",
    "cache/data/tcdb/substrates.tsv",
    "cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations",
    "cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv",
]


@pytest.mark.slow
def test_phase_1_1b_full_pipeline_med4():
    """Run sub-step 7 + step 2 against real cache. Assert metabolism report sane."""
    project_root = Path(__file__).parent.parent
    for rel in REQUIRED_FIXTURES:
        if not (project_root / rel).exists():
            pytest.skip(f"Missing fixture: {rel}")

    # 1. Build resolver + hierarchies (sub-step 7)
    subprocess.run(
        ["uv", "run", "python", "-m",
         "multiomics_kg.download.download_genome_data", "--steps", "7", "--force"],
        cwd=project_root, check=True,
    )

    assert (project_root / "cache/data/mnx/metabolite_resolver.db").exists()
    assert (project_root / "cache/data/tcdb/tcdb_hierarchy.json").exists()
    assert (project_root / "cache/data/cazy/cazy_hierarchy.json").exists()
    report = json.loads((project_root / "cache/data/mnx/metabolite_id_mapping_report.json").read_text())
    assert report["compound_count"] >= 100_000

    # 2. Run the gene annotation merge for MED4
    subprocess.run(
        ["uv", "run", "python", "-m",
         "multiomics_kg.download.build_gene_annotations", "--strains", "MED4", "--force"],
        cwd=project_root, check=True,
    )

    # 3. Assert metabolism report shows expected sanity values for MED4
    strain_dir = project_root / "cache/data/Prochlorococcus/genomes/MED4"
    report_path = strain_dir / "step2_metabolism_report.json"
    assert report_path.exists()
    metabolism = json.loads(report_path.read_text())

    assert metabolism["strain"] == "MED4"
    assert metabolism["gene_count"] >= 1500  # MED4 has 1838 genes
    # 526 KEGG_Reaction-annotated genes; resolution rate is the unknown.
    # ≥ 400 unique resolved is the spec's lower bound.
    assert metabolism["kegg_reactions"]["resolved_unique_mnxr"] >= 400
    # 121 KEGG_TC-annotated genes; spec lower bound 50.
    assert metabolism["transporter_classification"]["validated_unique"] >= 50
    # CAZy coverage on MED4 is sparse (1.3%); spec lower bound 5.
    assert metabolism["cazy_ids"]["validated_unique"] >= 5
```

- [ ] **Step 2: Run the smoke test**

Run: `uv run pytest tests/test_prepare_data_step2_metabolism_smoke.py -v -m slow`
Expected: PASS (with all fixtures present). If any fixture is missing, the test SKIPs; investigate before continuing.

- [ ] **Step 3: Run the full unit-test sweep to confirm no regression**

Run: `uv run pytest -m "not slow and not kg" -v 2>&1 | tail -20`
Expected: all metabolism unit tests pass; nothing else broken.

- [ ] **Step 4: Commit**

```bash
git add tests/test_prepare_data_step2_metabolism_smoke.py
git commit -m "metabolism: 1.1B — integration smoke test (slow marker)"
```

---

## Done state

After all 13 tasks:

- ✅ `cache/data/mnx/metabolite_resolver.db` populated (SQLite with compounds + aliases + names + reactions + reaction_aliases tables, ≥ 100K compounds).
- ✅ `cache/data/tcdb/tcdb_hierarchy.json` populated (5-level hierarchy + superfamilies + substrate_classes).
- ✅ `cache/data/cazy/cazy_hierarchy.json` populated (3-level hierarchy bootstrapped from observed eggNOG `CAZy` columns).
- ✅ `cache/data/mnx/metabolite_id_mapping_report.json` written with diagnostic counts.
- ✅ Three accessor modules (`metabolite_utils`, `tcdb_utils`, `cazy_utils`) have real implementations + tests.
- ✅ Three transforms (`_tx_resolve_kegg_reaction_to_mnxr`, `_tx_validate_tcdb`, `_tx_validate_cazy`) added to `annotation_transforms.py`.
- ✅ YAML config wired: `kegg_reactions` (renamed from singular), `transporter_classification`, `cazy_ids` carry the right transforms.
- ✅ Per-strain `step2_metabolism_report.json` written next to `gene_annotations_merged.json` for every strain.
- ✅ All non-slow tests pass: `pytest -m "not slow and not kg" -v`.
- ✅ Slow integration test passes for MED4: `pytest tests/test_prepare_data_step2_metabolism_smoke.py -m slow`.
- ✅ No KG impact: `schema_config.yaml` unchanged, no adapter change, no post-import change. The deployed graph remains byte-identical.

**Next plan (Phase 1.2):** the reactions scaffold (`Metabolite` + `Reaction` nodes; `Gene_catalyzes_reaction`; `Organism_has_metabolite` post-import edge; driver query). Spec already exists at `docs/superpowers/specs/2026-04-28-metabolite-reactions-scaffold-design.md`.
