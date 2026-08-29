# TigrRole hierarchy + NcbifamFamily→TigrRole bridge + equivalog gene roles — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `TigrRole` a two-level ontology, bridge `NcbifamFamily` TIGR* nodes to it from JCVI's frozen TIGRFAMs 15.0 role archive, emit equivalog-gated `Gene_has_tigr_role` edges for every strain, and use the same inference to fill `gene_category = Unknown` and add `[tigr_role_inferred]` description lines.

**Architecture:** A new committed reference `cache/data/ncbifam/tigr_roles.json` (built by prepare_data step 9) feeds three consumers: the step-2 merge (`build_gene_annotations.py`: category fill + description lines), `functional_annotation_adapter.py` (TigrRole nodes at two levels, `is_a` edges, inferred gene edges merged with curated ones), and `ncbifam_adapter.py` (the `Ncbifam_family_has_tigr_role` bridge). One pure helper module `multiomics_kg/utils/tigr_roles.py` owns the equivalog gate and slug rules so all three consumers agree. Post-import Cypher switches the TigrRole rollup to the CyanorakRole subtree form.

**Tech Stack:** Python 3.12 (uv), BioCypher adapters, pytest, Neo4j Cypher (`scripts/post-import.sh` + `.cypher` kept identical), Docker rebuild in the second clone.

**Spec:** `docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md`

## Global Constraints

- String sanitization: every string property yielded by an adapter passes through the adapter's local `_clean_str` (`'`→`^`, `|` removed). Never use `|` or `'` as separators in computed strings.
- Vocabulary rules R1–R5 (`docs/kg-changes/vocabulary-contract.md`): values the KG mints are lowercase `snake_case`; every `sources` value must be a `DataSource` id (`data_source:<value>`); no native `bool`; new closed value sets are declared in `config/controlled_vocabularies.yaml`.
- Equivalog gate: a `Gene_has_tigr_role` edge, a `gene_category` fill, or a `[tigr_role_inferred]` line is produced ONLY from `ncbifam_ids` accessions whose reference `family_type == "equivalog"` (exact string) AND which have a role in `tigr_roles.json["family_role"]`.
- Edge id for gene→TigrRole is `{locus_tag}-tigrrole-{code}` for BOTH curated and inferred edges; the adapter merges them before yielding (one edge per (gene, role)).
- Node ids: subroles `tigr.role:<numeric code>` (unchanged); mainroles `tigr.role:<slug>` where slug = lowercase, runs of non-`[a-z0-9]` → `_`, stripped of leading/trailing `_`. One prefix only — never `tigr.mainrole:`.
- Role `719` (unnamed in the archive) never appears anywhere.
- `scripts/post-import.sh` and `scripts/post-import.cypher` must carry identical Cypher logic.
- Docker build happens in the OTHER clone (`multiomics_biocypher_kg`); commit + push from this repo, the user pulls and rebuilds there.
- Run unit tests with `uv run pytest -m "not slow and not kg" -q`; KG tests with `uv run pytest -m kg -v` against a running graph.

---

### Task 1: Archive parsers in `multiomics_kg/utils/ncbifam.py`

**Files:**
- Modify: `multiomics_kg/utils/ncbifam.py` (append after `parse_hmm_pgap_rows`)
- Test: `tests/test_ncbifam.py` (append)

**Interfaces:**
- Produces: `parse_tigr_role_names(lines: Iterable[str]) -> dict[str, dict]` → `{role_id: {"mainrole": str, "sub1role": str}}`, only roles with a non-empty mainrole.
- Produces: `parse_tigr_role_link(lines: Iterable[str], roles: dict[str, dict]) -> dict[str, str]` → `{"TIGR00001": "158"}`, dropping links whose role is not in `roles`.
- Both raise `ValueError` when given non-empty input that yields zero entries (fail-loud, `_require_parsed` precedent in `kegg_utils`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ncbifam.py`:

```python
import pytest

from multiomics_kg.utils.ncbifam import parse_tigr_role_link, parse_tigr_role_names


ROLE_NAMES = """\
role_id:\t100\tmainrole:\tCentral intermediary metabolism
role_id:\t100\tsub1role:\tAmino sugars
role_id:\t132\tmainrole:\tDNA metabolism
role_id:\t132\tsub1role:\tDNA replication, recombination, and repair
role_id:\t719\tsub1role:\tOrphan without mainrole
""".splitlines()

ROLE_LINK = """\
TIGR00001\t132
TIGR00002\t100
TIGR00003\t719
TIGR00004\t999
""".splitlines()


def test_parse_tigr_role_names_pairs_main_and_sub():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert roles["100"] == {"mainrole": "Central intermediary metabolism", "sub1role": "Amino sugars"}
    assert roles["132"]["sub1role"] == "DNA replication, recombination, and repair"


def test_parse_tigr_role_names_drops_roles_without_mainrole():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert "719" not in roles


def test_parse_tigr_role_link_keeps_only_named_roles():
    roles = parse_tigr_role_names(ROLE_NAMES)
    link = parse_tigr_role_link(ROLE_LINK, roles)
    assert link == {"TIGR00001": "132", "TIGR00002": "100"}


def test_parse_tigr_role_link_strips_version_suffix():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert parse_tigr_role_link(["TIGR00005.1\t132"], roles) == {"TIGR00005": "132"}


def test_parsers_fail_loud_on_nonempty_garbage():
    with pytest.raises(ValueError):
        parse_tigr_role_names(["this is not a role line", "neither is this"])
    with pytest.raises(ValueError):
        parse_tigr_role_link(["garbage"], {"132": {"mainrole": "x", "sub1role": "y"}})


def test_parsers_accept_empty_input():
    assert parse_tigr_role_names([]) == {}
    assert parse_tigr_role_link([], {}) == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ncbifam.py -q -k tigr_role`
Expected: FAIL with `ImportError: cannot import name 'parse_tigr_role_link'`

- [ ] **Step 3: Implement the parsers**

Append to `multiomics_kg/utils/ncbifam.py`:

```python
# ── TIGRFAMs 15.0 role archive (frozen 2018) ───────────────────────────────
#
# NCBI FTP: https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/
#   TIGR_ROLE_NAMES   -- "role_id:\t<id>\tmainrole:\t<text>" / "...\tsub1role:\t<text>"
#   TIGRFAMS_ROLE_LINK -- "<TIGRxxxxx>\t<role_id>"
# NCBIfam kept the TIGR accessions but dropped the role column, so this archive
# is the ONLY source of family -> JCVI role. Role ids are the same id space the
# Cyanorak GFF `tIGR_Role` field uses (110/110 shared codes name-identical,
# measured 2026-08-28).


def parse_tigr_role_names(lines: Iterable[str]) -> dict[str, dict]:
    """Parse ``TIGR_ROLE_NAMES`` → ``{role_id: {"mainrole", "sub1role"}}``.

    Roles with no ``mainrole`` line (e.g. ``719`` in release 15.0) are dropped —
    an unnamed role cannot become a node. Raises ``ValueError`` when non-empty
    input yields nothing (format drift, never a real result).
    """
    roles: dict[str, dict] = {}
    n_lines = 0
    for line in lines:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        n_lines += 1
        parts = line.split("\t")
        if len(parts) < 4 or parts[0].strip() != "role_id:":
            continue
        role_id = parts[1].strip()
        kind = parts[2].strip().rstrip(":")
        text = parts[3].strip()
        if kind not in ("mainrole", "sub1role") or not role_id or not text:
            continue
        roles.setdefault(role_id, {})[kind] = text
    if n_lines and not roles:
        raise ValueError("TIGR_ROLE_NAMES: non-empty input parsed to zero roles — format drift?")
    return {
        rid: {"mainrole": r["mainrole"], "sub1role": r.get("sub1role", "")}
        for rid, r in roles.items()
        if r.get("mainrole")
    }


def parse_tigr_role_link(lines: Iterable[str], roles: dict[str, dict]) -> dict[str, str]:
    """Parse ``TIGRFAMS_ROLE_LINK`` → ``{unversioned_TIGR_acc: role_id}``.

    Links to roles absent from *roles* (unnamed or unknown) are dropped so the
    result is closed over the named-role set. Raises ``ValueError`` when
    non-empty input yields no parseable pair.
    """
    out: dict[str, str] = {}
    n_lines = n_parsed = 0
    for line in lines:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        n_lines += 1
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip().split(".", 1)[0]
        role_id = parts[1].strip()
        if not acc.startswith("TIGR") or not role_id:
            continue
        n_parsed += 1
        if role_id in roles:
            out[acc] = role_id
    if n_lines and not n_parsed:
        raise ValueError("TIGRFAMS_ROLE_LINK: non-empty input parsed to zero links — format drift?")
    return out
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_ncbifam.py -q`
Expected: all PASS (existing + 6 new)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/ncbifam.py tests/test_ncbifam.py
git commit -m "feat(ncbifam): parse TIGRFAMs 15.0 role archive (TIGR_ROLE_NAMES + TIGRFAMS_ROLE_LINK)"
```

---

### Task 2: Step 9 writes `cache/data/ncbifam/tigr_roles.json`

**Files:**
- Modify: `multiomics_kg/download/build_ncbifam_reference.py`
- Modify: `scripts/prepare_data.sh:59-71` (step 9 comment block)
- Modify: `.gitignore:187-189` comment (raw archive files are also gitignored under `cache/data/ncbifam/raw/` — already covered by the existing pattern; only the comment changes)
- Create (generated, committed): `cache/data/ncbifam/tigr_roles.json`
- Test: `tests/test_build_ncbifam_reference.py` (create)

**Interfaces:**
- Produces: `build_tigr_roles(force: bool = False, refetch_raw: bool = False) -> dict` writing `TIGR_ROLES_JSON = CACHE_DIR / "tigr_roles.json"` with shape `{"release": "TIGRFAMs 15.0 (frozen 2018)", "roles": {...}, "family_role": {...}}`.
- Produces: module constants `TIGR_ROLE_NAMES_URL`, `TIGRFAMS_ROLE_LINK_URL`, `TIGR_ROLE_NAMES_RAW`, `TIGRFAMS_ROLE_LINK_RAW`, `TIGR_ROLES_JSON`.
- `build()` (existing) now also calls `build_tigr_roles(force, refetch_raw)` so `prepare_data.sh --steps 9` needs no new command.

- [ ] **Step 1: Write the failing test**

Create `tests/test_build_ncbifam_reference.py`:

```python
"""Tests for the tigr_roles.json sub-builder of step 9 (offline: raw files are
written into a tmp cache dir; no network)."""
import json

import pytest

from multiomics_kg.download import build_ncbifam_reference as b


ROLE_NAMES = (
    "role_id:\t132\tmainrole:\tDNA metabolism\n"
    "role_id:\t132\tsub1role:\tDNA replication, recombination, and repair\n"
    "role_id:\t719\tsub1role:\tOrphan\n"
)
ROLE_LINK = "TIGR00001\t132\nTIGR00003\t719\n"


@pytest.fixture
def cache(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "TIGR_ROLE_NAMES").write_text(ROLE_NAMES)
    (raw / "TIGRFAMS_ROLE_LINK").write_text(ROLE_LINK)
    monkeypatch.setattr(b, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(b, "RAW_DIR", raw)
    monkeypatch.setattr(b, "TIGR_ROLE_NAMES_RAW", raw / "TIGR_ROLE_NAMES")
    monkeypatch.setattr(b, "TIGRFAMS_ROLE_LINK_RAW", raw / "TIGRFAMS_ROLE_LINK")
    monkeypatch.setattr(b, "TIGR_ROLES_JSON", tmp_path / "tigr_roles.json")
    return tmp_path


def test_build_tigr_roles_writes_expected_shape(cache):
    out = b.build_tigr_roles(force=True)
    assert out["release"].startswith("TIGRFAMs 15.0")
    assert out["roles"] == {"132": {"mainrole": "DNA metabolism",
                                    "sub1role": "DNA replication, recombination, and repair"}}
    assert out["family_role"] == {"TIGR00001": "132"}
    on_disk = json.loads((cache / "tigr_roles.json").read_text())
    assert on_disk == out


def test_build_tigr_roles_reuses_existing_without_force(cache):
    (cache / "tigr_roles.json").write_text(json.dumps({"release": "x", "roles": {}, "family_role": {}}))
    assert b.build_tigr_roles(force=False) == {"release": "x", "roles": {}, "family_role": {}}


def test_build_tigr_roles_download_failure_reuses_committed(cache, monkeypatch):
    (cache / "tigr_roles.json").write_text(json.dumps({"release": "kept", "roles": {}, "family_role": {}}))
    def boom(url, dest):
        raise RuntimeError("ftp down")
    monkeypatch.setattr(b, "_download", boom)
    out = b.build_tigr_roles(force=True, refetch_raw=True)
    assert out["release"] == "kept"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_build_ncbifam_reference.py -q`
Expected: FAIL with `AttributeError: module ... has no attribute 'TIGR_ROLE_NAMES_RAW'`

- [ ] **Step 3: Implement the sub-builder**

In `multiomics_kg/download/build_ncbifam_reference.py`:

Change the import line to:
```python
from multiomics_kg.utils.ncbifam import (
    HYPOTH_FAMILY_TYPES,
    parse_hmm_pgap_rows,
    parse_tigr_role_link,
    parse_tigr_role_names,
)
```

After `REFERENCE_JSON = CACHE_DIR / "ncbifam_reference.json"` add:
```python
# TIGRFAMs 15.0 role archive (frozen 2018) — the only surviving family→role map.
TIGRFAMS_15_BASE = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0"
TIGR_ROLE_NAMES_URL = f"{TIGRFAMS_15_BASE}/TIGR_ROLE_NAMES"
TIGRFAMS_ROLE_LINK_URL = f"{TIGRFAMS_15_BASE}/TIGRFAMS_ROLE_LINK"
TIGR_ROLE_NAMES_RAW = RAW_DIR / "TIGR_ROLE_NAMES"
TIGRFAMS_ROLE_LINK_RAW = RAW_DIR / "TIGRFAMS_ROLE_LINK"
TIGR_ROLES_JSON = CACHE_DIR / "tigr_roles.json"
TIGR_ROLES_RELEASE = "TIGRFAMs 15.0 (frozen 2018)"
```

After `_ensure_raw` add:
```python
def _load_tigr_roles_json() -> dict | None:
    if TIGR_ROLES_JSON.exists():
        with open(TIGR_ROLES_JSON, encoding="utf-8") as fh:
            return json.load(fh)
    return None


def build_tigr_roles(force: bool = False, refetch_raw: bool = False) -> dict:
    """Build (and cache) ``tigr_roles.json`` from the frozen TIGRFAMs 15.0 archive.

    ``{"release", "roles": {role_id: {mainrole, sub1role}}, "family_role":
    {TIGR_acc: role_id}}``. Unnamed roles (``719``) are excluded from both maps.
    On a download failure the committed file is reused with a warning (TCDB
    outage precedent); only a missing file is fatal.
    """
    existing = _load_tigr_roles_json()
    if existing is not None and not force and not refetch_raw:
        logger.info("TIGR roles cache exists: %s (use --force to rebuild)", TIGR_ROLES_JSON)
        return existing

    try:
        for url, dest in ((TIGR_ROLE_NAMES_URL, TIGR_ROLE_NAMES_RAW),
                          (TIGRFAMS_ROLE_LINK_URL, TIGRFAMS_ROLE_LINK_RAW)):
            if refetch_raw or not dest.exists():
                _download(url, dest)
    except Exception as exc:  # network / FTP outage
        if existing is not None:
            logger.warning("TIGRFAMs archive download failed (%s); reusing committed %s",
                           exc, TIGR_ROLES_JSON)
            return existing
        raise

    with open(TIGR_ROLE_NAMES_RAW, encoding="utf-8") as fh:
        roles = parse_tigr_role_names(fh)
    with open(TIGRFAMS_ROLE_LINK_RAW, encoding="utf-8") as fh:
        family_role = parse_tigr_role_link(fh, roles)

    out = {"release": TIGR_ROLES_RELEASE, "roles": roles, "family_role": family_role}
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(TIGR_ROLES_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    logger.info("Wrote %s: %d named roles, %d family→role links",
                TIGR_ROLES_JSON, len(roles), len(family_role))
    return out
```

In `build()`, immediately before the final `return ref` (after the "Wrote %s" log), add:
```python
    build_tigr_roles(force=force, refetch_raw=refetch_raw)
```

Also in `build()`'s early-return branch (the `if REFERENCE_JSON.exists() and not force and not refetch_raw:` block), add `build_tigr_roles(force=False, refetch_raw=False)` before the `return json.load(fh)` so a first run on a checkout that has `ncbifam_reference.json` but no `tigr_roles.json` still creates it. Concretely replace that block with:
```python
    if REFERENCE_JSON.exists() and not force and not refetch_raw:
        logger.info("NCBIfam reference cache exists: %s (use --force to rebuild)", REFERENCE_JSON)
        build_tigr_roles(force=False, refetch_raw=False)
        with open(REFERENCE_JSON, encoding="utf-8") as fh:
            return json.load(fh)
```

Update the module docstring's "One source file:" list to add:
```
- ``TIGR_ROLE_NAMES`` + ``TIGRFAMS_ROLE_LINK`` (NCBI FTP,
  ``https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/``) — the frozen
  2018 JCVI role archive; written to ``cache/data/ncbifam/tigr_roles.json``
  (separate file so ``ncbifam_reference.json`` keeps its flat ``{acc: …}`` shape).
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_build_ncbifam_reference.py tests/test_ncbifam.py -q`
Expected: PASS

- [ ] **Step 5: Generate the committed artefact and check it**

Run: `uv run python -m multiomics_kg.download.build_ncbifam_reference --force`
(Needs network for the two small archive files; `hmm_PGAP.tsv` is reused from `cache/data/ncbifam/raw/` if present, else downloaded — ~5 MB.)

Then verify:
```bash
uv run python - <<'EOF'
import json
d = json.load(open("cache/data/ncbifam/tigr_roles.json"))
print(len(d["roles"]), len(d["family_role"]))
assert "719" not in d["roles"] and "719" not in set(d["family_role"].values())
assert d["roles"]["132"]["mainrole"] == "DNA metabolism"
print(sorted({r["mainrole"] for r in d["roles"].values()}))
EOF
```
Expected: `116 2920` (116 named roles; 2,963 links minus the 43 that pointed at `719`), and 21 distinct mainrole strings, every one present as a key in `TIGR_TO_CATEGORY` in `multiomics_kg/download/build_gene_annotations.py:135` (check by eye; Task 4 asserts it in code).

- [ ] **Step 6: Update the prepare_data.sh step-9 comment**

In `scripts/prepare_data.sh` step-9 comment block (lines 59–71) append one line after the `cache/data/ncbifam/ncbifam_reference.json` line:
```
#           Also writes cache/data/ncbifam/tigr_roles.json (TIGRFAMs 15.0 frozen role
#           archive: roles + TIGR family→role), consumed by step 2 (gene_category fill +
#           [tigr_role_inferred] lines) and by functional_annotation/ncbifam adapters.
```

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/download/build_ncbifam_reference.py tests/test_build_ncbifam_reference.py scripts/prepare_data.sh cache/data/ncbifam/tigr_roles.json
git commit -m "feat(step9): build committed tigr_roles.json from the TIGRFAMs 15.0 role archive"
```

---

### Task 3: Shared inference helper `multiomics_kg/utils/tigr_roles.py`

**Files:**
- Create: `multiomics_kg/utils/tigr_roles.py`
- Test: `tests/test_tigr_roles.py` (create)

**Interfaces:**
- Produces: `EQUIVALOG_TYPES: frozenset[str] = frozenset({"equivalog"})`
- Produces: `load_tigr_roles(cache_root: Path) -> dict | None` — reads `<cache_root>/ncbifam/tigr_roles.json`, `None` when absent.
- Produces: `mainrole_slug(mainrole: str) -> str`.
- Produces: `role_name(role_id: str, tigr_roles: dict) -> str` → `"<mainrole> / <sub1role>"` (or bare mainrole when sub1role empty).
- Produces: `split_role_name(name: str) -> tuple[str, str | None]` → `(mainrole, sub1role|None)` splitting on the first `" / "`.
- Produces: `inferred_roles(ncbifam_ids: list[str] | None, tigr_roles: dict, ncbifam_ref: dict) -> dict[str, list[str]]` → `{role_id: [supporting TIGR accessions]}`, equivalog-gated.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_tigr_roles.py`:

```python
import json

from multiomics_kg.utils.tigr_roles import (
    EQUIVALOG_TYPES,
    inferred_roles,
    load_tigr_roles,
    mainrole_slug,
    role_name,
    split_role_name,
)

TIGR_ROLES = {
    "release": "TIGRFAMs 15.0 (frozen 2018)",
    "roles": {
        "132": {"mainrole": "DNA metabolism", "sub1role": "DNA replication, recombination, and repair"},
        "120": {"mainrole": "Energy metabolism", "sub1role": "TCA cycle"},
        "156": {"mainrole": "Hypothetical proteins", "sub1role": "Conserved"},
    },
    "family_role": {"TIGR00001": "132", "TIGR00002": "120", "TIGR00003": "156", "TIGR00004": "120"},
}
NCBIFAM_REF = {
    "TIGR00001": {"name": "a", "family_type": "equivalog"},
    "TIGR00002": {"name": "b", "family_type": "subfamily"},
    "TIGR00003": {"name": "c", "family_type": "equivalog"},
    "TIGR00004": {"name": "d", "family_type": "equivalog"},
}


def test_equivalog_types_is_exactly_equivalog():
    assert EQUIVALOG_TYPES == frozenset({"equivalog"})


def test_mainrole_slug():
    assert mainrole_slug("Energy metabolism") == "energy_metabolism"
    assert mainrole_slug("Purines, pyrimidines, nucleosides, and nucleotides") == \
        "purines_pyrimidines_nucleosides_and_nucleotides"
    assert mainrole_slug("Biosynthesis of cofactors, prosthetic groups, and carriers") == \
        "biosynthesis_of_cofactors_prosthetic_groups_and_carriers"


def test_role_name_and_split_round_trip():
    assert role_name("120", TIGR_ROLES) == "Energy metabolism / TCA cycle"
    assert split_role_name("Energy metabolism / TCA cycle") == ("Energy metabolism", "TCA cycle")
    assert split_role_name("Not Found") == ("Not Found", None)


def test_inferred_roles_equivalog_gate():
    out = inferred_roles(["TIGR00001", "TIGR00002", "NF000001"], TIGR_ROLES, NCBIFAM_REF)
    assert out == {"132": ["TIGR00001"]}          # TIGR00002 is subfamily → dropped


def test_inferred_roles_groups_supporting_accessions():
    out = inferred_roles(["TIGR00004", "TIGR00002"], TIGR_ROLES, NCBIFAM_REF)
    assert out == {"120": ["TIGR00004"]}


def test_inferred_roles_keeps_junk_roles():
    assert inferred_roles(["TIGR00003"], TIGR_ROLES, NCBIFAM_REF) == {"156": ["TIGR00003"]}


def test_inferred_roles_empty_inputs():
    assert inferred_roles(None, TIGR_ROLES, NCBIFAM_REF) == {}
    assert inferred_roles(["TIGR00001"], TIGR_ROLES, {}) == {}   # no ref → no family_type → no edge


def test_load_tigr_roles(tmp_path):
    assert load_tigr_roles(tmp_path) is None
    (tmp_path / "ncbifam").mkdir()
    (tmp_path / "ncbifam" / "tigr_roles.json").write_text(json.dumps(TIGR_ROLES))
    assert load_tigr_roles(tmp_path) == TIGR_ROLES
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_tigr_roles.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'multiomics_kg.utils.tigr_roles'`

- [ ] **Step 3: Implement the helper**

Create `multiomics_kg/utils/tigr_roles.py`:

```python
"""Shared rules for inferring JCVI TIGR roles from NCBIfam TIGR* family hits.

Three consumers must agree on the gate and the naming, so the rules live here:
- step-2 merge (``build_gene_annotations``): ``gene_category`` fill-only rule +
  ``[tigr_role_inferred]`` description lines;
- ``functional_annotation_adapter``: inferred ``Gene_has_tigr_role`` edges +
  two-level ``TigrRole`` nodes;
- ``ncbifam_adapter``: ``Ncbifam_family_has_tigr_role`` bridge (ungated —
  the bridge is ontology→ontology and records every archive link).

Design: docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md
"""

from __future__ import annotations

import json
import re
from pathlib import Path

# Only equivalog families ("same function in every member") transfer a role to
# a gene. subfamily / equivalog_domain / hypoth_equivalog / domain … never do —
# measured 2026-08-29: the gate cuts multi-role genes 375→15 and
# Cyanorak contradictions 7%→5.4% at the cost of ~1/3 of heterotroph reach.
EQUIVALOG_TYPES: frozenset[str] = frozenset({"equivalog"})

_SLUG_RE = re.compile(r"[^a-z0-9]+")
ROLE_SEP = " / "


def load_tigr_roles(cache_root: Path) -> dict | None:
    """Load ``<cache_root>/ncbifam/tigr_roles.json`` (step 9); ``None`` if absent."""
    path = Path(cache_root) / "ncbifam" / "tigr_roles.json"
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def mainrole_slug(mainrole: str) -> str:
    """``"Energy metabolism"`` → ``"energy_metabolism"`` (mainrole node local id)."""
    return _SLUG_RE.sub("_", mainrole.strip().lower()).strip("_")


def role_name(role_id: str, tigr_roles: dict) -> str:
    """Compound display name ``"<mainrole> / <sub1role>"`` (Cyanorak convention)."""
    r = tigr_roles["roles"][role_id]
    return f"{r['mainrole']}{ROLE_SEP}{r['sub1role']}" if r.get("sub1role") else r["mainrole"]


def split_role_name(name: str) -> tuple[str, str | None]:
    """Inverse of :func:`role_name` for compound names carried by Cyanorak."""
    if ROLE_SEP in name:
        main, sub = name.split(ROLE_SEP, 1)
        return main.strip(), sub.strip()
    return name.strip(), None


def inferred_roles(
    ncbifam_ids: list[str] | None,
    tigr_roles: dict,
    ncbifam_ref: dict,
) -> dict[str, list[str]]:
    """Equivalog-gated ``{role_id: [supporting TIGR accessions]}`` for one gene.

    An accession contributes iff ``ncbifam_ref[acc]["family_type"]`` is in
    :data:`EQUIVALOG_TYPES` AND ``tigr_roles["family_role"]`` maps it. Junk
    roles (156/157/…) are returned like any other — the nodes carry
    ``is_uninformative`` downstream; hiding them here would hide that a
    family is "hypothetical".
    """
    out: dict[str, list[str]] = {}
    if not ncbifam_ids or not tigr_roles or not ncbifam_ref:
        return out
    family_role = tigr_roles.get("family_role") or {}
    for acc in ncbifam_ids:
        if not acc or acc not in family_role:
            continue
        if (ncbifam_ref.get(acc) or {}).get("family_type") not in EQUIVALOG_TYPES:
            continue
        out.setdefault(family_role[acc], []).append(acc)
    return {k: sorted(v) for k, v in sorted(out.items())}
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_tigr_roles.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tigr_roles.py tests/test_tigr_roles.py
git commit -m "feat(tigr-roles): shared equivalog-gated role inference helper"
```

---

### Task 4: Step-2 merge — `gene_category` fill + `[tigr_role_inferred]` lines

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py` (new function after `enrich_pfam_fields`; call in `process_strain` after `enrich_interpro_fields`; `main()` loads `tigr_roles`; `process_strain` gains a `tigr_roles` kwarg)
- Test: `tests/test_build_gene_annotations.py` (append a class)

**Interfaces:**
- Consumes: `inferred_roles`, `role_name`, `load_tigr_roles` from Task 3; `TIGR_TO_CATEGORY` (existing, line 135).
- Produces: `apply_tigr_role_inference(gene: dict, tigr_roles: dict | None, ncbifam_ref: dict | None) -> None` (in-place).
- Produces: `process_strain(..., tigr_roles: dict | None = None)`.
- Produces: `_check_tigr_roles_mapped(tigr_roles: dict) -> None` — raises `AssertionError` listing any archive mainrole missing from `TIGR_TO_CATEGORY`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_build_gene_annotations.py`:

```python
# ─── TIGR role inference from equivalog NCBIfam hits (step-2, spec §3.5/§3.6) ──

from multiomics_kg.download.build_gene_annotations import (
    apply_tigr_role_inference,
    _check_tigr_roles_mapped,
)

_TIGR_ROLES = {
    "release": "t",
    "roles": {
        "120": {"mainrole": "Energy metabolism", "sub1role": "TCA cycle"},
        "132": {"mainrole": "DNA metabolism", "sub1role": "DNA replication, recombination, and repair"},
        "157": {"mainrole": "Unknown function", "sub1role": "General"},
    },
    "family_role": {"TIGR00001": "120", "TIGR00002": "132", "TIGR00003": "157", "TIGR00005": "120"},
}
_NCBIFAM_REF = {
    "TIGR00001": {"name": "a", "family_type": "equivalog"},
    "TIGR00002": {"name": "b", "family_type": "equivalog"},
    "TIGR00003": {"name": "c", "family_type": "equivalog"},
    "TIGR00004": {"name": "d", "family_type": "subfamily"},
    "TIGR00005": {"name": "e", "family_type": "equivalog"},
}


class TestApplyTigrRoleInference:
    def test_fills_unknown_category(self):
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00001"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Energy production"

    def test_never_overwrites_existing_category(self):
        g = {"gene_category": "Translation", "ncbifam_ids": ["TIGR00001"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Translation"

    def test_most_frequent_category_wins_then_alphabetical(self):
        # 120 twice (Energy production) vs 132 once (Replication and repair)
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00001", "TIGR00005", "TIGR00002"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Energy production"
        # tie: Energy production vs Replication and repair → alphabetical
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00001", "TIGR00002"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Energy production"

    def test_junk_role_does_not_fill_category_but_adds_description(self):
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00003"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Unknown"
        assert g["alternate_functional_descriptions"] == ["[tigr_role_inferred] Unknown function / General"]

    def test_non_equivalog_hit_is_ignored(self):
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00004"]}
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["gene_category"] == "Unknown"
        assert "alternate_functional_descriptions" not in g

    def test_description_lines_deduped_against_curated(self):
        g = {
            "gene_category": "Energy production",
            "ncbifam_ids": ["TIGR00001", "TIGR00002"],
            "alternate_functional_descriptions": ["[tigr_role] Energy metabolism / TCA cycle"],
        }
        apply_tigr_role_inference(g, _TIGR_ROLES, _NCBIFAM_REF)
        assert g["alternate_functional_descriptions"] == [
            "[tigr_role] Energy metabolism / TCA cycle",
            "[tigr_role_inferred] DNA metabolism / DNA replication, recombination, and repair",
        ]

    def test_noop_when_reference_missing(self):
        g = {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00001"]}
        apply_tigr_role_inference(g, None, _NCBIFAM_REF)
        apply_tigr_role_inference(g, _TIGR_ROLES, None)
        assert g == {"gene_category": "Unknown", "ncbifam_ids": ["TIGR00001"]}


def test_check_tigr_roles_mapped_passes_for_known_mainroles():
    _check_tigr_roles_mapped(_TIGR_ROLES)


def test_check_tigr_roles_mapped_raises_on_unknown_mainrole():
    bad = {"roles": {"1": {"mainrole": "Quantum metabolism", "sub1role": "x"}}, "family_role": {}}
    with pytest.raises(AssertionError, match="Quantum metabolism"):
        _check_tigr_roles_mapped(bad)


def test_committed_tigr_roles_all_mapped():
    """Every mainrole in the committed archive JSON has a TIGR_TO_CATEGORY entry."""
    import json
    from pathlib import Path
    p = Path("cache/data/ncbifam/tigr_roles.json")
    if not p.exists():
        pytest.skip("tigr_roles.json not built")
    _check_tigr_roles_mapped(json.loads(p.read_text()))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_build_gene_annotations.py -q -k "TigrRoleInference or tigr_roles_mapped"`
Expected: FAIL with `ImportError: cannot import name 'apply_tigr_role_inference'`

- [ ] **Step 3: Implement**

In `multiomics_kg/download/build_gene_annotations.py`:

Add to the imports (after `from multiomics_kg.utils.pfam_utils import PfamData, load_pfam_data`):
```python
from multiomics_kg.utils.tigr_roles import inferred_roles, load_tigr_roles, role_name
```

After the `enrich_pfam_fields` function add:
```python
# ─── TIGR role inference from equivalog NCBIfam hits (spec §3.5 / §3.6) ──────

_TIGR_INFERRED_LABEL = "[tigr_role_inferred]"
_TIGR_CURATED_LABEL = "[tigr_role]"


def _check_tigr_roles_mapped(tigr_roles: dict) -> None:
    """Fail loud if the archive names a mainrole TIGR_TO_CATEGORY cannot map
    (an archive refresh must never silently produce 'Unknown')."""
    missing = sorted({r["mainrole"] for r in (tigr_roles.get("roles") or {}).values()}
                     - set(TIGR_TO_CATEGORY))
    assert not missing, f"tigr_roles.json mainroles missing from TIGR_TO_CATEGORY: {missing}"


def apply_tigr_role_inference(gene: dict, tigr_roles: dict | None,
                              ncbifam_ref: dict | None) -> None:
    """Post-merge: fill ``gene_category`` (only when 'Unknown') and append
    ``[tigr_role_inferred] <Main> / <Sub>`` description lines from the gene's
    equivalog TIGR* hits. Priority 4 in the category chain (after COG):
    Cyanorak role → Cyanorak TIGR role → COG → NCBIfam-bridged TIGR role.
    """
    if not tigr_roles or not ncbifam_ref:
        return
    roles = inferred_roles(gene.get("ncbifam_ids"), tigr_roles, ncbifam_ref)
    if not roles:
        return

    # 1. gene_category fill-only: most frequent non-Unknown category, ties → alphabetical.
    if gene.get("gene_category", "Unknown") == "Unknown":
        votes: dict[str, int] = {}
        for role_id, accs in roles.items():
            cat = TIGR_TO_CATEGORY.get(tigr_roles["roles"][role_id]["mainrole"])
            if cat and cat != "Unknown":
                votes[cat] = votes.get(cat, 0) + len(accs)
        if votes:
            gene["gene_category"] = sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]

    # 2. [tigr_role_inferred] lines, deduped against the curated [tigr_role] text.
    afd = list(gene.get("alternate_functional_descriptions") or [])
    curated_texts = {s[len(_TIGR_CURATED_LABEL):].strip() for s in afd
                     if s.startswith(_TIGR_CURATED_LABEL + " ")}
    present = set(afd)
    for role_id in roles:
        text = role_name(role_id, tigr_roles)
        line = f"{_TIGR_INFERRED_LABEL} {text}"
        if text in curated_texts or line in present:
            continue
        afd.append(line)
        present.add(line)
    if afd:
        gene["alternate_functional_descriptions"] = afd
```

In `process_strain`, change the signature to add `tigr_roles: dict | None = None,` after `ncbifam_ref`, and right after the `enrich_interpro_fields(...)` call (inside the `if interpro_ref is not None:` block's sibling position — i.e. after that `if` block, before `merged_out[locus_tag] = merged`) add:
```python
        apply_tigr_role_inference(merged, tigr_roles, ncbifam_ref)
```

In `main()`, after `print(f"NCBIfam reference: {len(ncbifam_ref)} entries")` add:
```python
    tigr_roles = load_tigr_roles(cache_root)
    if tigr_roles is None:
        print("WARNING: cache/data/ncbifam/tigr_roles.json missing — TIGR-role inference "
              "(gene_category fill, [tigr_role_inferred]) disabled. Run prepare_data.sh --steps 9.")
    else:
        _check_tigr_roles_mapped(tigr_roles)
        print(f"TIGR roles archive: {len(tigr_roles['roles'])} roles, "
              f"{len(tigr_roles['family_role'])} family links")
```
and pass `tigr_roles=tigr_roles` in the `process_strain(...)` call.

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_build_gene_annotations.py -q`
Expected: PASS (all existing + new)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_gene_annotations.py tests/test_build_gene_annotations.py
git commit -m "feat(merge): equivalog TIGR-role inference — gene_category fill-only + [tigr_role_inferred] lines"
```

---

### Task 5: Two-level `TigrRole` nodes + `Tigr_role_is_a_tigr_role` (adapter + schema)

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py` (`MultiCogRoleAnnotationAdapter.__init__`, `_all_tigr_codes`, `get_nodes`, `get_edges`; new `tigr_role_node_ids`)
- Modify: `config/schema_config.yaml:813-822` (add `level_kind`), and add the hierarchy edge type after `cyanorak role hierarchical association` (line ~1771)
- Test: `tests/test_cog_role_annotation_adapter.py` (append)

**Interfaces:**
- Consumes: `mainrole_slug`, `split_role_name` from Task 3.
- Produces: `MultiCogRoleAnnotationAdapter(genome_config_file, role_tree_file, test_mode=False, extra_tigr_roles: dict[str, str] | None = None, tigr_roles: dict | None = None, ncbifam_ref: dict | None = None)` — `tigr_roles`/`ncbifam_ref` are used by Task 6; accept them now so the signature is final.
- Produces: `_tigr_mainrole_node_id(mainrole: str) -> str` = `f"tigr.role:{mainrole_slug(mainrole)}"`.
- Produces: `MultiCogRoleAnnotationAdapter.tigr_role_node_ids() -> set[str]` — every TigrRole node id `get_nodes()` emits (subroles + mainroles).
- Node props: subrole `{"code", "name", "level": 1, "level_kind": "tigr_subrole"}`; mainrole `{"code": slug, "name": mainrole, "level": 0, "level_kind": "tigr_mainrole"}`; a code whose name has no `" / "` → `{"code", "name", "level": 0, "level_kind": "tigr_mainrole"}` and no parent edge.
- Edge: `(f"{code}-tigr_is_a-{slug}", tigr.role:<code>, tigr.role:<slug>, "tigr_role_is_a_tigr_role", {})`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_cog_role_annotation_adapter.py` (uses the module's existing `multi_adapter` fixture and `MINI_COG_GENE_DATA`, whose tIGR descriptions are `"Some tIGR role"` / `"Another tIGR role"` / `"Third tIGR role"` — no `" / "`; extend the fixture data through `extra_tigr_roles` instead):

```python
from multiomics_kg.adapters.functional_annotation_adapter import _tigr_mainrole_node_id


class TestTigrRoleHierarchy:
    @pytest.fixture
    def adapter_with_extra(self, multi_adapter_factory):
        return multi_adapter_factory(extra_tigr_roles={
            "120": "Energy metabolism / TCA cycle",
            "108": "Energy metabolism / Aerobic",
            "132": "DNA metabolism / DNA replication, recombination, and repair",
        })

    def test_mainrole_node_id(self):
        assert _tigr_mainrole_node_id("Energy metabolism") == "tigr.role:energy_metabolism"

    def test_subroles_level_one_and_mainroles_level_zero(self, adapter_with_extra):
        nodes = {n[0]: n[2] for n in adapter_with_extra.get_nodes() if n[1] == "tigr role"}
        assert nodes["tigr.role:120"] == {"code": "120", "name": "Energy metabolism / TCA cycle",
                                          "level": 1, "level_kind": "tigr_subrole"}
        assert nodes["tigr.role:energy_metabolism"] == {"code": "energy_metabolism",
                                                        "name": "Energy metabolism",
                                                        "level": 0, "level_kind": "tigr_mainrole"}
        # mainrole emitted ONCE although two subroles share it
        assert sum(1 for k in nodes if k == "tigr.role:energy_metabolism") == 1
        assert "tigr.role:dna_metabolism" in nodes

    def test_code_without_separator_is_a_root(self, adapter_with_extra):
        nodes = {n[0]: n[2] for n in adapter_with_extra.get_nodes() if n[1] == "tigr role"}
        # "12345" → "Some tIGR role" (fixture) has no " / "
        assert nodes["tigr.role:12345"]["level"] == 0
        assert nodes["tigr.role:12345"]["level_kind"] == "tigr_mainrole"

    def test_is_a_edges(self, adapter_with_extra):
        edges = [e for e in adapter_with_extra.get_edges() if e[3] == "tigr_role_is_a_tigr_role"]
        pairs = {(e[1], e[2]) for e in edges}
        assert ("tigr.role:120", "tigr.role:energy_metabolism") in pairs
        assert ("tigr.role:108", "tigr.role:energy_metabolism") in pairs
        assert ("tigr.role:132", "tigr.role:dna_metabolism") in pairs
        assert not any(src == "tigr.role:12345" for src, _ in pairs)
        assert all(e[4] == {} for e in edges)

    def test_tigr_role_node_ids_matches_emitted_nodes(self, adapter_with_extra):
        emitted = {n[0] for n in adapter_with_extra.get_nodes() if n[1] == "tigr role"}
        assert adapter_with_extra.tigr_role_node_ids() == emitted
```

The existing `multi_adapter` fixture (line 134) is:
```python
@pytest.fixture
def multi_adapter(genome_config_csv, roles_csv_file):
    return MultiCogRoleAnnotationAdapter(
        genome_config_file=genome_config_csv,
        role_tree_file=roles_csv_file,
    )
```
Add this sibling fixture directly after it:
```python
@pytest.fixture
def multi_adapter_factory(genome_config_csv, roles_csv_file):
    """Same inputs as `multi_adapter`, forwarding kwargs (extra_tigr_roles, tigr_roles, …)."""
    def _make(**kwargs):
        return MultiCogRoleAnnotationAdapter(
            genome_config_file=genome_config_csv,
            role_tree_file=roles_csv_file,
            **kwargs,
        )
    return _make
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cog_role_annotation_adapter.py -q -k TigrRoleHierarchy`
Expected: FAIL with `ImportError: cannot import name '_tigr_mainrole_node_id'`

- [ ] **Step 3: Implement**

In `multiomics_kg/adapters/functional_annotation_adapter.py`:

Add import near the other utils imports:
```python
from multiomics_kg.utils.tigr_roles import mainrole_slug, split_role_name
```

After `_tigr_role_node_id` add:
```python
def _tigr_mainrole_node_id(mainrole: str) -> str:
    """Node ID for a TIGR mainrole — same house-minted prefix as subroles
    (one prefix per ontology, CyanorakRole precedent); local id is a slug,
    never numeric, so it cannot collide with subrole codes."""
    return f"tigr.role:{mainrole_slug(mainrole)}"
```

Change `MultiCogRoleAnnotationAdapter.__init__` to:
```python
    def __init__(
        self,
        genome_config_file: str,
        role_tree_file: Path,
        test_mode: bool = False,
        extra_tigr_roles: dict[str, str] | None = None,
        tigr_roles: dict | None = None,
        ncbifam_ref: dict | None = None,
    ) -> None:
        self.test_mode = test_mode
        self.role_tree = parse_cyanorak_role_tree(Path(role_tree_file))
        logger.info(f"MultiCogRoleAnnotationAdapter: loaded {len(self.role_tree)} Cyanorak role nodes")
        # {code: "Main / Sub"} for TIGR roles reached only by NCBIfam inference
        # (bridge or equivalog gene edge) and absent from every strain's Cyanorak data.
        self.extra_tigr_roles: dict[str, str] = dict(extra_tigr_roles or {})
        # Task 6: archive + NCBIfam reference for inferred gene edges.
        self.tigr_roles = tigr_roles
        self.ncbifam_ref = ncbifam_ref
        self._strain_adapters: list[CogRoleAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)
```

Replace `_all_tigr_codes` with:
```python
    def _all_tigr_codes(self) -> dict[str, str]:
        """All (code → compound description) pairs: Cyanorak-observed ∪ extra_tigr_roles.
        Cyanorak's text wins when both name a code (they are identical in practice —
        110/110 measured)."""
        codes: dict[str, str] = {}
        for adapter in self._strain_adapters:
            for code, desc in adapter.get_all_tigr_codes():
                if code not in codes:
                    codes[code] = desc
        for code, desc in self.extra_tigr_roles.items():
            codes.setdefault(code, desc)
        return codes

    def _tigr_nodes(self) -> tuple[list[tuple[str, dict]], list[tuple[str, str]]]:
        """Compute (nodes, is_a pairs) for the TigrRole ontology.

        nodes: [(node_id, props)] — subroles (level 1) then mainroles (level 0, deduped).
        pairs: [(subrole_id, mainrole_id)].
        """
        nodes: list[tuple[str, dict]] = []
        pairs: list[tuple[str, str]] = []
        mainroles: dict[str, str] = {}  # node_id -> mainrole text
        for code, desc in sorted(self._all_tigr_codes().items()):
            main, sub = split_role_name(desc)
            if sub is None or not main:
                nodes.append((_tigr_role_node_id(code),
                              {"code": code, "name": _clean_str(desc),
                               "level": 0, "level_kind": "tigr_mainrole"}))
                continue
            nodes.append((_tigr_role_node_id(code),
                          {"code": code, "name": _clean_str(desc),
                           "level": 1, "level_kind": "tigr_subrole"}))
            mid = _tigr_mainrole_node_id(main)
            mainroles.setdefault(mid, main)
            pairs.append((_tigr_role_node_id(code), mid))
        for mid, main in sorted(mainroles.items()):
            nodes.append((mid, {"code": mainrole_slug(main), "name": _clean_str(main),
                                "level": 0, "level_kind": "tigr_mainrole"}))
        return nodes, pairs

    def tigr_role_node_ids(self) -> set[str]:
        """Every TigrRole node id get_nodes() emits (for bridge dangling guards)."""
        return {nid for nid, _ in self._tigr_nodes()[0]}
```

In `get_nodes()`, replace the "3. TigrRole nodes" block with:
```python
        # 3. TigrRole nodes — Cyanorak-observed ∪ inference-reached codes, two levels
        tigr_nodes, _ = self._tigr_nodes()
        for node_id, props in tigr_nodes:
            yield (node_id, "tigr role", props)
        logger.info(f"MultiCogRoleAnnotationAdapter.get_nodes: {len(tigr_nodes)} TigrRole nodes "
                    f"({sum(1 for _, p in tigr_nodes if p['level'] == 1)} subroles, "
                    f"{sum(1 for _, p in tigr_nodes if p['level'] == 0)} mainroles/roots)")
```

In `get_edges()`, after the CyanorakRole hierarchy block (end of the method), add:
```python
        # 5. TigrRole hierarchy: subrole → mainrole
        _, pairs = self._tigr_nodes()
        for sub_id, main_id in pairs:
            yield (
                f"{sub_id.split(':', 1)[1]}-tigr_is_a-{main_id.split(':', 1)[1]}",
                sub_id,
                main_id,
                "tigr_role_is_a_tigr_role",
                {},
            )
        logger.info(f"MultiCogRoleAnnotationAdapter.get_edges: {len(pairs)} TigrRole is_a edges")
```
and update the `get_edges` docstring list to include `5. TigrRole→mainrole`.

In `config/schema_config.yaml`, in `tigr role:` properties add `level_kind: str` after `level: int`. After the `cyanorak role hierarchical association:` block add:
```yaml
tigr role hierarchical association:
  is_a: association
  represented_as: edge
  label_as_edge: tigr_role_is_a_tigr_role
  source: tigr role
  target: tigr role
  label_in_input: tigr_role_is_a_tigr_role
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_cog_role_annotation_adapter.py -q`
Expected: PASS. (If an existing test asserts every TigrRole node has `level == 0`, update it to assert `level in (0, 1)`.)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py config/schema_config.yaml tests/test_cog_role_annotation_adapter.py
git commit -m "feat(tigrrole): two-level TigrRole nodes + Tigr_role_is_a_tigr_role hierarchy"
```

---

### Task 6: Equivalog-gated inferred `Gene_has_tigr_role` edges, merged with curated

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py` (`CogRoleAnnotationAdapter.__init__`/`get_edges`, `MultiCogRoleAnnotationAdapter._build_strain_adapters`, new `get_inferred_tigr_codes`)
- Test: `tests/test_cog_role_annotation_adapter.py` (append)

**Interfaces:**
- Consumes: `inferred_roles` (Task 3); `self.tigr_roles`, `self.ncbifam_ref` from Task 5's constructor.
- Produces: `CogRoleAnnotationAdapter(genome_dir, test_mode=False, tigr_roles=None, ncbifam_ref=None)`.
- Produces: `CogRoleAnnotationAdapter.get_inferred_tigr_codes() -> set[str]` — role ids any gene of this strain reaches via the gate (used by `create_knowledge_graph` in Task 8 to build `extra_tigr_roles`).
- Produces: module constant `_INFERRED_TIGR_EDGE_PROPS = {"sources": ["interproscan"], "evidence": "family_inferred"}`.
- Edge merge rule: for one gene, curated roles (from `tIGR_Role`) and inferred roles (gate) are merged per code: `sources` = sorted union, `evidence` = `"curated"` if curated present else `"family_inferred"`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_cog_role_annotation_adapter.py`:

```python
import json as _json

_TR = {
    "release": "t",
    "roles": {
        "120": {"mainrole": "Energy metabolism", "sub1role": "TCA cycle"},
        "132": {"mainrole": "DNA metabolism", "sub1role": "DNA replication, recombination, and repair"},
        "156": {"mainrole": "Hypothetical proteins", "sub1role": "Conserved"},
    },
    "family_role": {"TIGR00001": "120", "TIGR00002": "132", "TIGR00003": "156", "TIGR00004": "120"},
}
_NR = {
    "TIGR00001": {"name": "a", "family_type": "equivalog"},
    "TIGR00002": {"name": "b", "family_type": "equivalog"},
    "TIGR00003": {"name": "c", "family_type": "equivalog"},
    "TIGR00004": {"name": "d", "family_type": "subfamily"},
}


def _strain_dir(tmp_path, genes: dict):
    d = tmp_path / "STRAIN"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(_json.dumps(genes))
    return d


def _tigr_edges(adapter, locus_tag):
    return {e[2]: e[4] for e in adapter.get_edges()
            if e[3] == "gene_has_tigr_role" and e[1].endswith(locus_tag)}


class TestInferredTigrRoleEdges:
    def test_inferred_edge_props(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "ncbifam_ids": ["TIGR00001"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        edges = _tigr_edges(a, "G1")
        assert edges == {"tigr.role:120": {"sources": ["interproscan"], "evidence": "family_inferred"}}

    def test_subfamily_hit_gives_no_edge(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "ncbifam_ids": ["TIGR00004"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        assert _tigr_edges(a, "G1") == {}

    def test_curated_and_inferred_same_role_merge_into_one_edge(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "tIGR_Role": ["120"],
                   "tIGR_Role_description": ["Energy metabolism / TCA cycle"],
                   "ncbifam_ids": ["TIGR00001"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        all_edges = [e for e in a.get_edges() if e[3] == "gene_has_tigr_role"]
        assert len(all_edges) == 1
        eid, src, tgt, _, props = all_edges[0]
        assert eid == "G1-tigrrole-120"
        assert props == {"sources": ["cyanorak", "interproscan"], "evidence": "curated"}

    def test_disagreement_yields_two_edges(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "tIGR_Role": ["120"],
                   "tIGR_Role_description": ["Energy metabolism / TCA cycle"],
                   "ncbifam_ids": ["TIGR00002"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        edges = _tigr_edges(a, "G1")
        assert edges == {
            "tigr.role:120": {"sources": ["cyanorak"], "evidence": "curated"},
            "tigr.role:132": {"sources": ["interproscan"], "evidence": "family_inferred"},
        }

    def test_junk_role_still_emitted(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "ncbifam_ids": ["TIGR00003"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        assert "tigr.role:156" in _tigr_edges(a, "G1")

    def test_no_reference_means_curated_only(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "tIGR_Role": ["120"],
                   "tIGR_Role_description": ["Energy metabolism / TCA cycle"],
                   "ncbifam_ids": ["TIGR00001"]}}))
        assert _tigr_edges(a, "G1") == {"tigr.role:120": {"sources": ["cyanorak"], "evidence": "curated"}}

    def test_get_inferred_tigr_codes(self, tmp_path):
        a = CogRoleAnnotationAdapter(_strain_dir(tmp_path, {
            "G1": {"locus_tag": "G1", "ncbifam_ids": ["TIGR00001", "TIGR00004"]},
            "G2": {"locus_tag": "G2", "ncbifam_ids": ["TIGR00003"]}}), tigr_roles=_TR, ncbifam_ref=_NR)
        assert a.get_inferred_tigr_codes() == {"120", "156"}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cog_role_annotation_adapter.py -q -k InferredTigrRoleEdges`
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'tigr_roles'`

- [ ] **Step 3: Implement**

In `multiomics_kg/adapters/functional_annotation_adapter.py`:

Extend the import from Task 5 to:
```python
from multiomics_kg.utils.tigr_roles import inferred_roles, mainrole_slug, split_role_name
```

After `_CYANORAK_EDGE_PROPS = ...` add:
```python
# Inferred TIGR role: an equivalog NCBIfam family (InterProScan NCBIFAM facet)
# whose JCVI role is known from the frozen TIGRFAMs 15.0 archive. Merged with
# the curated Cyanorak edge for the same (gene, role): sources = union,
# evidence = 'curated' when Cyanorak agrees (curated outranks family_inferred
# on the ladder; TCDB eggNOG+diamond precedent).
_INFERRED_TIGR_EDGE_PROPS = {"sources": ["interproscan"], "evidence": "family_inferred"}
```

Change `CogRoleAnnotationAdapter.__init__` to:
```python
    def __init__(self, genome_dir: Path, test_mode: bool = False,
                 tigr_roles: dict | None = None, ncbifam_ref: dict | None = None) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self.tigr_roles = tigr_roles
        self.ncbifam_ref = ncbifam_ref
        self._genes: dict = {}
        self._load()
```

After `get_all_tigr_codes` add:
```python
    def _gene_tigr_roles(self, gene: dict) -> dict[str, dict]:
        """{role_code: edge_props} for one gene — curated (tIGR_Role) merged with
        equivalog-inferred (ncbifam_ids); one entry per role."""
        roles: dict[str, dict] = {}
        for code in gene.get("tIGR_Role") or []:
            if code:
                roles[code] = dict(_CYANORAK_EDGE_PROPS)
        if self.tigr_roles and self.ncbifam_ref:
            for code in inferred_roles(gene.get("ncbifam_ids"), self.tigr_roles, self.ncbifam_ref):
                if code in roles:
                    roles[code] = {
                        "sources": sorted(set(roles[code]["sources"]) | set(_INFERRED_TIGR_EDGE_PROPS["sources"])),
                        "evidence": "curated",
                    }
                else:
                    roles[code] = dict(_INFERRED_TIGR_EDGE_PROPS)
        return roles

    def get_inferred_tigr_codes(self) -> set[str]:
        """Role codes reached by any gene of this strain through the equivalog gate."""
        if not (self.tigr_roles and self.ncbifam_ref):
            return set()
        out: set[str] = set()
        for gene in self._genes.values():
            out.update(inferred_roles(gene.get("ncbifam_ids"), self.tigr_roles, self.ncbifam_ref))
        return out
```

In `CogRoleAnnotationAdapter.get_edges`, replace the `# gene → tIGR role` loop with:
```python
            # gene → tIGR role (curated ∪ equivalog-inferred, merged per role)
            for code, props in sorted(self._gene_tigr_roles(gene).items()):
                yield (
                    f"{locus_tag}-tigrrole-{code}",
                    _gene_node_id(locus_tag),
                    _tigr_role_node_id(code),
                    "gene_has_tigr_role",
                    props,
                )
                tigr_count += 1
                if self.test_mode and tigr_count >= 100:
                    logger.debug(
                        f"CogRoleAnnotationAdapter({self.genome_dir.name}): test_mode stop (TigrRole)"
                    )
                    return
```
(Note `props` is already a fresh dict per call — no `dict(...)` copy needed.)

In `MultiCogRoleAnnotationAdapter._build_strain_adapters`, pass the references through:
```python
            self._strain_adapters.append(
                CogRoleAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode,
                                         tigr_roles=self.tigr_roles, ncbifam_ref=self.ncbifam_ref)
            )
```
(`self.tigr_roles`/`self.ncbifam_ref` are set in `__init__` before `_build_strain_adapters` is called — Task 5 placed them there.)

Add to `MultiCogRoleAnnotationAdapter`:
```python
    def inferred_tigr_codes(self) -> set[str]:
        """Union of every strain's equivalog-reached role codes."""
        out: set[str] = set()
        for adapter in self._strain_adapters:
            out |= adapter.get_inferred_tigr_codes()
        return out
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_cog_role_annotation_adapter.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py tests/test_cog_role_annotation_adapter.py
git commit -m "feat(tigrrole): equivalog-gated inferred Gene_has_tigr_role edges merged with curated"
```

---

### Task 7: `Ncbifam_family_has_tigr_role` bridge (adapter + schema)

**Files:**
- Modify: `multiomics_kg/adapters/ncbifam_adapter.py` (`MultiNcbifamAdapter.__init__`, `download_data`, `get_edges`)
- Modify: `config/schema_config.yaml` (after the `ncbifam family to interpro entry association` block, ~line 1719)
- Test: `tests/test_ncbifam_adapter.py` (append)

**Interfaces:**
- Consumes: `load_tigr_roles` (Task 3).
- Produces: `MultiNcbifamAdapter(..., tigr_role_node_ids: set[str] | None = None)`; `self._tigr_roles: dict | None` loaded in `download_data()`.
- Produces: `MultiNcbifamAdapter.observed_ids() -> set[str]` — public alias of the existing private `_observed_ids()` (add `def observed_ids(self) -> set[str]: return self._observed_ids()` right after it), used by `create_knowledge_graph` in Task 8.
- Produces: `_tigr_role_node_id(code) -> f"tigr.role:{code}"` (local copy; must equal `functional_annotation_adapter._tigr_role_node_id`).
- Edge: `(f"{acc}-has_tigr_role-{code}", ncbifam:<acc>, tigr.role:<code>, "ncbifam_family_has_tigr_role", {})`, only for observed `TIGR*` accessions whose archive role node id is in `tigr_role_node_ids`; `None` → no bridge edges.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ncbifam_adapter.py` (reuse the file's `_write_strain` helper and the pattern used by its existing Multi-adapter tests to build a genome-config CSV — copy the CSV-writing lines from the existing `MultiNcbifamAdapter` test in that file):

```python
def _multi(tmp_path, genes, calls, **kwargs):
    genome_dir = _write_strain(tmp_path, genes, calls)
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("strain_name,data_dir\nTESTSTRAIN," + str(genome_dir) + "\n")
    cache_root = tmp_path / "cache"
    (cache_root / "ncbifam").mkdir(parents=True)
    (cache_root / "ncbifam" / "ncbifam_reference.json").write_text(json.dumps({
        "TIGR00001": {"name": "a", "family_type": "equivalog"},
        "TIGR00002": {"name": "b", "family_type": "subfamily"},
        "NF000001": {"name": "n", "family_type": "equivalog"},
    }))
    (cache_root / "ncbifam" / "tigr_roles.json").write_text(json.dumps({
        "release": "t",
        "roles": {"120": {"mainrole": "Energy metabolism", "sub1role": "TCA cycle"},
                  "132": {"mainrole": "DNA metabolism", "sub1role": "x"}},
        "family_role": {"TIGR00001": "120", "TIGR00002": "132"},
    }))
    a = MultiNcbifamAdapter(genome_config_file=str(cfg), cache_root=cache_root, **kwargs)
    a.download_data()
    return a


_GENES = {"LT001": {"protein_id": "WP_1", "ncbifam_ids": ["TIGR00001", "TIGR00002", "NF000001"]}}
_CALLS = {"WP_1": {"libraries": {"NCBIFAM": [
    {"accession": "TIGR00001", "name": "a", "start": 1, "end": 9, "evalue": 1e-9, "score": 10.0},
    {"accession": "TIGR00002", "name": "b", "start": 1, "end": 9, "evalue": 1e-9, "score": 10.0},
    {"accession": "NF000001", "name": "n", "start": 1, "end": 9, "evalue": 1e-9, "score": 10.0},
]}}}


def _bridge(adapter):
    return {(e[1], e[2]): e[4] for e in adapter.get_edges() if e[3] == "ncbifam_family_has_tigr_role"}


def test_tigr_role_bridge_ungated_and_dangling_proof(tmp_path):
    a = _multi(tmp_path, _GENES, _CALLS, tigr_role_node_ids={"tigr.role:120", "tigr.role:132"})
    assert _bridge(a) == {
        ("ncbifam:TIGR00001", "tigr.role:120"): {},
        ("ncbifam:TIGR00002", "tigr.role:132"): {},   # subfamily still bridges (ontology-level)
    }


def test_tigr_role_bridge_skips_missing_target_node(tmp_path):
    a = _multi(tmp_path, _GENES, _CALLS, tigr_role_node_ids={"tigr.role:120"})
    assert set(_bridge(a)) == {("ncbifam:TIGR00001", "tigr.role:120")}


def test_tigr_role_bridge_none_means_no_edges(tmp_path):
    a = _multi(tmp_path, _GENES, _CALLS)
    assert _bridge(a) == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ncbifam_adapter.py -q -k tigr_role_bridge`
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'tigr_role_node_ids'`

- [ ] **Step 3: Implement**

In `multiomics_kg/adapters/ncbifam_adapter.py`:

Add import:
```python
from multiomics_kg.utils.tigr_roles import load_tigr_roles
```

After `_ncbifam_node_id` add:
```python
def _tigr_role_node_id(code: str) -> str:
    """Must match functional_annotation_adapter._tigr_role_node_id (bridge target)."""
    return f"tigr.role:{code}"
```

Right after `_observed_ids` add:
```python
    def observed_ids(self) -> set[str]:
        """Public: distinct NCBIfam accessions observed across all strains."""
        return self._observed_ids()
```

In `MultiNcbifamAdapter.__init__` add parameter `tigr_role_node_ids: set[str] | None = None,` after `interpro_kept_ids`, and in the body after `self.interpro_kept_ids = interpro_kept_ids`:
```python
        # None = TigrRole node-id set not provided → emit NO TigrRole bridge edges.
        self.tigr_role_node_ids = tigr_role_node_ids
        self._tigr_roles: dict | None = None
```

In `download_data`, after the reference load log line:
```python
        self._tigr_roles = load_tigr_roles(self.cache_root)
        if self._tigr_roles is None:
            logger.warning("tigr_roles.json missing under %s — no Ncbifam_family_has_tigr_role edges",
                           self.cache_root / "ncbifam")
```

In `get_edges`, after the InterPro bridge loop (before `# 2. Gene → NcbifamFamily edges`), add:
```python
        # 1b. NcbifamFamily → TigrRole bridge (TIGRFAMs 15.0 archive; TIGR* only,
        # UNGATED by family_type — ontology→ontology record of JCVI's assignment;
        # the equivalog gate applies only to gene-level inference elsewhere).
        role_bridge = skipped = 0
        family_role = (self._tigr_roles or {}).get("family_role") or {}
        if self.tigr_role_node_ids is not None and family_role:
            for acc in sorted(observed):
                code = family_role.get(acc)
                if not code:
                    continue
                target = _tigr_role_node_id(code)
                if target not in self.tigr_role_node_ids:
                    skipped += 1
                    continue
                yield (
                    f"{acc}-has_tigr_role-{code}",
                    _ncbifam_node_id(acc),
                    target,
                    "ncbifam_family_has_tigr_role",
                    {},
                )
                role_bridge += 1
```
and extend the final log line to `f"MultiNcbifamAdapter.get_edges: {bridge} interpro-bridge, {role_bridge} tigr-role-bridge ({skipped} skipped: no TigrRole node), {gene} gene edges"`.

In `config/schema_config.yaml`, after the `ncbifam family to interpro entry association` block:
```yaml
ncbifam family to tigr role association:
  # TIGRFAMs 15.0 frozen archive (JCVI role assignment per family). No
  # properties: single frozen source, provenance lives on the edge type.
  # Router semantics — never assign a gene role from it (see
  # docs/kg-changes/tigr-role-bridge.md).
  is_a: association
  represented_as: edge
  label_as_edge: Ncbifam_family_has_tigr_role
  source: ncbifam family
  target: tigr role
  label_in_input: ncbifam_family_has_tigr_role
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_ncbifam_adapter.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/ncbifam_adapter.py config/schema_config.yaml tests/test_ncbifam_adapter.py
git commit -m "feat(ncbifam): Ncbifam_family_has_tigr_role bridge from the TIGRFAMs 15.0 archive"
```

---

### Task 8: Orchestration, vocabularies, uninformative flags, smoke build

**Files:**
- Modify: `create_knowledge_graph.py:255-262` and `:344-351`
- Modify: `config/controlled_vocabularies.yaml` (`Gene_has_tigr_role.sources`, `Gene_has_tigr_role.evidence`; new `TigrRole.level_kind` next to `CazyFamily.level_kind`)
- Modify: `config/uninformative_terms.yaml` (`tigr_role.ids`)
- Test: `tests/test_controlled_vocab.py` (existing static gates) + smoke build

**Interfaces:**
- Consumes: `MultiCogRoleAnnotationAdapter(extra_tigr_roles=, tigr_roles=, ncbifam_ref=)`, `.inferred_tigr_codes()`, `.tigr_role_node_ids()` (Tasks 5–6); `MultiNcbifamAdapter(tigr_role_node_ids=)` (Task 7); `load_tigr_roles`, `role_name` (Task 3).

- [ ] **Step 1: Wire `create_knowledge_graph.py`**

Replace the COG/role block (lines 255–262) with:
```python
    # COG functional categories + Cyanorak roles + tIGR roles.
    # TIGR roles: Cyanorak-curated (Pro/Syn) ∪ equivalog-inferred from NCBIfam
    # TIGR* hits via the frozen TIGRFAMs 15.0 archive (all strains). Roles that
    # only inference reaches (bridge or gene edge) get nodes through
    # extra_tigr_roles so nothing dangles. See
    # docs/kg-changes/tigr-role-bridge.md.
    from multiomics_kg.utils.tigr_roles import load_tigr_roles, role_name
    tigr_roles = load_tigr_roles(Path("cache/data"))
    ncbifam_ref_path = Path("cache/data/ncbifam/ncbifam_reference.json")
    ncbifam_ref = json.loads(ncbifam_ref_path.read_text(encoding="utf-8")) if ncbifam_ref_path.exists() else {}
    if tigr_roles is None:
        logger.warning("cache/data/ncbifam/tigr_roles.json missing — TigrRole inference + bridge disabled "
                       "(run prepare_data.sh --steps 9)")
    # Roles reachable by the ontology bridge = every observed TIGR* accession's archive role.
    from multiomics_kg.adapters.ncbifam_adapter import MultiNcbifamAdapter
    ncbifam_adapter = MultiNcbifamAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_root="cache/data",
        interpro_kept_ids=None,  # set below once interpro_adapter exists
        test_mode=TEST_MODE,
    )
    bridge_codes: set[str] = set()
    if tigr_roles is not None:
        family_role = tigr_roles["family_role"]
        bridge_codes = {family_role[a] for a in ncbifam_adapter.observed_ids() if a in family_role}
    cog_role_adapter = MultiCogRoleAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        role_tree_file=Path("data/cyanorak_roles.csv"),
        test_mode=TEST_MODE,
        tigr_roles=tigr_roles,
        ncbifam_ref=ncbifam_ref,
    )
    if tigr_roles is not None:
        reached = bridge_codes | cog_role_adapter.inferred_tigr_codes()
        cog_role_adapter.extra_tigr_roles = {c: role_name(c, tigr_roles) for c in reached}
    bc.write_nodes(cog_role_adapter.get_nodes())
    bc.write_edges(cog_role_adapter.get_edges())
    tigr_role_node_ids = cog_role_adapter.tigr_role_node_ids()
```
(`create_knowledge_graph.py` does NOT import `json` today — add `import json` to its stdlib import block at the top.)

Then in the NCBIfam block (lines ~344–351) delete the `from ... import MultiNcbifamAdapter` line and the `MultiNcbifamAdapter(...)` construction, replacing them with:
```python
    ncbifam_adapter.interpro_kept_ids = interpro_adapter.kept_node_accessions()
    ncbifam_adapter.tigr_role_node_ids = tigr_role_node_ids
```
keeping the existing `ncbifam_adapter.download_data(cache=CACHE)` and the materialize/guard lines unchanged. (`_observed_ids()` reads only each strain's merged JSON; constructing the adapter early costs one JSON load per strain that the adapter caches — no double work.)

- [ ] **Step 2: Vocabularies**

In `config/controlled_vocabularies.yaml`:

`Gene_has_tigr_role.sources` → `values: [cyanorak, interproscan]` and description:
```
    Who asserted this annotation (uniform provenance, KG-SYNC-005 / ONT-007).
    Each value joins a DataSource node via id = 'data_source:' + value.
    'cyanorak' = curated TIGR role carried by the Cyanorak GFF (Pro/Syn);
    'interproscan' = equivalog NCBIfam TIGR* hit → JCVI role via the frozen
    TIGRFAMs 15.0 archive (all strains, 2026-08-29). Both on one edge = the
    two agree. Prefer `'cyanorak' IN r.sources` over list equality.
```
`Gene_has_tigr_role.evidence` → `values: [curated, family_inferred]` and append to its description: `family_inferred = equivalog-family role transfer (no Cyanorak curation); curated wins when both sources name the role.`

Add after `CazyFamily.level_kind`:
```yaml
TigrRole.level_kind:
  applies_to: TigrRole
  applies_to_kind: node
  property: level_kind
  value_type: string
  closed: true
  values: [tigr_mainrole, tigr_subrole]
  description: >
    JCVI two-level role scheme (2026-08-29): tigr_mainrole = level 0 (node id
    tigr.role:<slug>, or a Cyanorak-only code whose name has no ' / ');
    tigr_subrole = level 1 (tigr.role:<numeric code>, compound name
    '<Main> / <Sub>'). Hierarchy edge Tigr_role_is_a_tigr_role (sub → main).
```

In `config/uninformative_terms.yaml`, under `tigr_role.ids` add:
```yaml
    - tigr.role:hypothetical_proteins   # mainrole (level 0) — 2026-08-29 hierarchy
    - tigr.role:unknown_function        # mainrole
    - tigr.role:unclassified            # mainrole
```
and extend the trailing comment: `# 856 "Not Found" is itself a root (no ' / '), already listed above.`

- [ ] **Step 3: Run the static gates**

Run: `uv run pytest tests/test_controlled_vocab.py tests/test_annotation_quality_buckets.py -q`
Expected: PASS. If `test_evidence_ladder_placements` or `test_every_gene_ontology_edge_declares_sources_and_evidence` fails, read the failing assertion and align the YAML (the values above are the intended truth).

- [ ] **Step 4: Smoke build in test mode**

Run: `uv run python create_knowledge_graph.py --test --output-dir ./biocypher-log/tigr_smoke/ 2>&1 | tee logs/tigr_smoke.log | grep -E "TigrRole|tigr-role-bridge|is_a edges|Traceback|Error"`
Expected: log lines `… TigrRole nodes (… subroles, … mainroles/roots)`, `… TigrRole is_a edges`, `… tigr-role-bridge (…)`, no Traceback. Then check the CSV headers exist:
```bash
ls biocypher-log/tigr_smoke/*/ | grep -i -E "TigrRole|Ncbifam_family_has_tigr_role|Tigr_role_is_a"
head -2 biocypher-log/tigr_smoke/*/TigrRole-header.csv
```
Expected: `TigrRole-header.csv` contains `level_kind`; both edge header files exist.

- [ ] **Step 5: Commit**

```bash
git add create_knowledge_graph.py config/controlled_vocabularies.yaml config/uninformative_terms.yaml
git commit -m "feat(kg): wire TigrRole inference + bridge into the build; register vocabularies + junk mainrole flags"
```

---

### Task 9: Post-import Cypher (subtree rollups, `ncbifam_family_count`, indexes, flags)

**Files:**
- Modify: `scripts/post-import.cypher` (index block ~line 68; F1.1 flags ~line 613; TigrRole rollup ~line 1281)
- Modify: `scripts/post-import.sh` (same three places: ~87, ~262, ~1283) — identical Cypher

**Interfaces:**
- Produces on `TigrRole`: `gene_count` (subtree), `direct_gene_count`, `organism_count`, `ncbifam_family_count` (subtree count of incoming `Ncbifam_family_has_tigr_role`).
- Indexes: `tigr_role_level_idx`, `tigr_role_level_kind_idx`.

- [ ] **Step 1: Indexes**

In both files, next to the `cazy_family_level_idx` / `cazy_family_level_kind_idx` lines add:
```cypher
CREATE INDEX tigr_role_level_idx IF NOT EXISTS FOR (t:TigrRole) ON (t.level);
CREATE INDEX tigr_role_level_kind_idx IF NOT EXISTS FOR (t:TigrRole) ON (t.level_kind);
```

- [ ] **Step 2: Flags**

In both files replace the TigrRole F1.1 statement with:
```cypher
MATCH (t:TigrRole)
WHERE t.id IN ['tigr.role:156','tigr.role:704','tigr.role:856',
               'tigr.role:185','tigr.role:157',
               'tigr.role:hypothetical_proteins','tigr.role:unknown_function',
               'tigr.role:unclassified']
SET t.is_uninformative = 'true';
```

- [ ] **Step 3: Rollups**

In both files replace the `MATCH (n:TigrRole) CALL { … }` rollup with the CyanorakRole subtree form plus the bridge count:
```cypher
// TigrRole — two-level since 2026-08-29: subtree gene_count over
// Tigr_role_is_a_tigr_role (CyanorakRole pattern) + direct_gene_count +
// ncbifam_family_count (incoming Ncbifam_family_has_tigr_role, subtree).
MATCH (n:TigrRole)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Tigr_role_is_a_tigr_role*0..]-(desc:TigrRole)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_has_tigr_role]-(g:Gene)
  WITH n, descs, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  UNWIND descs AS d2
  OPTIONAL MATCH (d2)<-[:Ncbifam_family_has_tigr_role]-(f:NcbifamFamily)
  WITH n, gc, dgc, orgs, count(DISTINCT f) AS fc
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL]),
      n.ncbifam_family_count = fc
} IN TRANSACTIONS OF 100 ROWS;
```

- [ ] **Step 4: Verify the two files stay in sync**

Run: `diff <(grep -A20 "MATCH (n:TigrRole)" scripts/post-import.cypher) <(grep -A20 "MATCH (n:TigrRole)" scripts/post-import.sh)` and the same for the flag and index lines.
Expected: no output (identical).

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "feat(post-import): TigrRole subtree rollups, ncbifam_family_count, level indexes, mainrole junk flags"
```

---

### Task 10: KG validity tests

**Files:**
- Modify: `tests/kg_validity/test_ontology_level.py:36-43,60-67`
- Modify: `tests/kg_validity/test_annotation_trust.py:19-26,44-52`
- Modify: `tests/kg_validity/test_functional_annotation.py` (append after line ~665)
- Modify: `tests/kg_validity/test_structure.py` — no change needed (`TigrRole` already listed); verify the relationship-type list there includes the two new edge labels if such a list exists (`grep -n "Ncbifam_family_in_interpro_entry" tests/kg_validity/test_structure.py`) and add `Tigr_role_is_a_tigr_role`, `Ncbifam_family_has_tigr_role` alongside it.

- [ ] **Step 1: `test_ontology_level.py`**

Remove `"TigrRole"` from `NON_GO_LABELS` if that list drives a "flat" assertion, and replace `test_flat_ontologies_all_level_zero`'s TigrRole block with:
```python
    rows = run_query(
        "MATCH (t:TigrRole) RETURN count(t) AS n, min(t.level) AS mn, max(t.level) AS mx, "
        "count(CASE WHEN t.level = 1 THEN 1 END) AS subs, count(CASE WHEN t.level = 0 THEN 1 END) AS mains"
    )[0]
    assert rows["mn"] == 0 and rows["mx"] == 1, rows
    assert rows["subs"] >= 110, f"TigrRole subroles {rows['subs']} < 110"
    assert 15 <= rows["mains"] <= 40, f"TigrRole mainroles/roots {rows['mains']} out of range"
```
Keep the CogFunctionalCategory block as is. Update the docstring to "CogFunctionalCategory is flat; TigrRole is two-level since 2026-08-29".

- [ ] **Step 2: `test_annotation_trust.py`**

Move `"TigrRole"` from `FLAT` into `HIERARCHICAL`. In `test_evidence_rungs_by_edge` change `"Gene_has_tigr_role": {"curated"}` to `{"curated", "family_inferred"}`. Check what the HIERARCHICAL loop asserts (it expects `direct_gene_count` — Task 9 sets it).

- [ ] **Step 3: New assertions in `test_functional_annotation.py`**

Append:
```python
# ---------------------------------------------------------------------------
# TigrRole hierarchy + NCBIfam bridge + inferred gene roles (2026-08-29)
# ---------------------------------------------------------------------------

def test_tigr_role_every_subrole_has_exactly_one_parent(run_query):
    bad = run_query("""
        MATCH (s:TigrRole {level: 1})
        OPTIONAL MATCH (s)-[r:Tigr_role_is_a_tigr_role]->(m:TigrRole {level: 0})
        WITH s, count(r) AS k WHERE k <> 1 RETURN count(s) AS bad
    """)[0]["bad"]
    assert bad == 0


def test_tigr_role_mainrole_gene_count_is_subtree(run_query):
    rows = run_query("""
        MATCH (m:TigrRole {level: 0})<-[:Tigr_role_is_a_tigr_role]-(s:TigrRole)
        WITH m, collect(s) AS subs
        MATCH (g:Gene)-[:Gene_has_tigr_role]->(x:TigrRole) WHERE x = m OR x IN subs
        WITH m, count(DISTINCT g) AS expected
        WHERE m.gene_count <> expected RETURN count(m) AS bad
    """)
    assert rows[0]["bad"] == 0


def test_ncbifam_tigr_role_bridge(run_query):
    row = run_query("""
        MATCH (f:NcbifamFamily)-[r:Ncbifam_family_has_tigr_role]->(t:TigrRole)
        RETURN count(r) AS n, count(CASE WHEN f.ncbifam_id STARTS WITH 'TIGR' THEN 1 END) AS tigr,
               count(DISTINCT f) AS fams, size(keys(r)) AS props
    """)[0]
    assert row["n"] >= 1600, row
    assert row["tigr"] == row["n"], "bridge sources must all be TIGR*"
    assert row["fams"] == row["n"], "one role per family"


def test_inferred_tigr_role_edges_span_all_organisms(run_query):
    row = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tigr_role]->()
        RETURN count(DISTINCT g.organism_name) AS orgs,
               count(CASE WHEN 'interproscan' IN r.sources THEN 1 END) AS inferred,
               count(CASE WHEN r.evidence = 'family_inferred' AND 'cyanorak' IN r.sources THEN 1 END) AS bad_merge,
               count(CASE WHEN NOT all(s IN r.sources WHERE s IN ['cyanorak','interproscan']) THEN 1 END) AS bad_src
    """)[0]
    assert row["orgs"] >= 40, row
    assert row["inferred"] >= 13_000, row
    assert row["bad_merge"] == 0 and row["bad_src"] == 0, row


def test_tigr_role_ncbifam_family_count_present(run_query):
    row = run_query("MATCH (t:TigrRole) RETURN count(t) AS n, count(t.ncbifam_family_count) AS c")[0]
    assert row["c"] == row["n"]
```

- [ ] **Step 4: Run the unit suite (KG tests auto-skip without Neo4j)**

Run: `uv run pytest -m "not slow and not kg" -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/kg_validity/
git commit -m "test(kg): TigrRole hierarchy, NCBIfam bridge and inferred-role invariants"
```

---

### Task 11: Documentation

**Files:**
- Create: `docs/kg-changes/tigr-role-bridge.md`
- Modify: `CLAUDE.md` (Key Adapters `ncbifam_adapter.py` bullet; Actual Neo4j labels relationships list; Key graph facts — new TigrRole bullet + NcbifamFamily bullet; Step 9 description; Data Locations; `informative_annotation_types` note unchanged)
- Modify: `CHANGELOG.md` under `## [Unreleased]`
- Modify: `plans/backlog.md:27-44` and `:65-76`
- Modify: `docs/kg-changes/vocabulary-contract.md` (new entries)
- Modify: `plans/interpro_redesign_backlog.md:37-40` (mark done)

- [ ] **Step 1: Write `docs/kg-changes/tigr-role-bridge.md`**

Content (write it fully; sections): **What changed** (two-level TigrRole with ids/level_kind; `Tigr_role_is_a_tigr_role`; `Ncbifam_family_has_tigr_role` from TIGRFAMs 15.0 with the "NCBI kept accessions, dropped roles" history; equivalog-gated `Gene_has_tigr_role` with `sources`/`evidence` and the merge rule; `gene_category` priority chain now Cyanorak → Cyanorak TIGR → COG → NCBIfam-bridged TIGR, fill-only; `[tigr_role_inferred]` lines; `TigrRole.ncbifam_family_count`; new flags; indexes), **Semantics — read this before querying** (bridge is a router; equivalog gate rationale with the table from spec §2; 557 curated/inferred contradictions are facet choices, kept as two edges; coverage bias 18% vs ~90% — ORA needs a per-organism background of genes with ≥1 role edge; `annotation_quality` unchanged by measurement), **Query cookbook** (heterotroph genes by role via `Gene_has_tigr_role` with `'interproscan' IN r.sources`; disagreements: `MATCH (g)-[a]->(x),(g)-[b]->(y) WHERE a.evidence='curated' AND b.evidence='family_inferred'`; role→families via the bridge), **Refresh procedure** (`prepare_data.sh --steps 9 --force` then `--steps 2 --force`, Docker rebuild), **Numbers measured 2026-08-29** (copy spec §2 tables).

- [ ] **Step 2: CLAUDE.md edits**

- Relationships list: add `Tigr_role_is_a_tigr_role`, `Ncbifam_family_has_tigr_role` after `Ncbifam_family_in_interpro_entry`.
- Add a Key-graph-facts bullet **TigrRole nodes (two-level since 2026-08-29)** summarising: ~140 nodes (level 0 mainroles `tigr.role:<slug>` + Cyanorak-only roots; level 1 subroles `tigr.role:<code>` with compound names), `level_kind` values, `Tigr_role_is_a_tigr_role`, computed `gene_count` (subtree) / `direct_gene_count` / `organism_count` / `ncbifam_family_count`; `Gene_has_tigr_role` on all 43 organisms — `sources` ∈ {cyanorak, interproscan}, `evidence` curated | family_inferred, equivalog gate, merge rule; coverage caveat; pointer to `docs/kg-changes/tigr-role-bridge.md`.
- NcbifamFamily bullet: append one sentence about `Ncbifam_family_has_tigr_role` (TIGR* only, ~1,720, no props, router).
- Gene `gene_category` mention (schema comment line 559 says "from Cyanorak/TIGR/COG"): add "→ NCBIfam-bridged TIGR (fill-only)" wherever the priority chain is described (grep `Cyanorak Role → TIGR Role → COG` in CLAUDE.md and `build_gene_annotations.py` docstrings and update).
- Step 9 paragraph + Data Locations: add `cache/data/ncbifam/tigr_roles.json` (built by step 9; roles + family→role; frozen TIGRFAMs 15.0; consumed by step 2 + two adapters).
- `alternate_functional_descriptions` labels: mention `[tigr_role_inferred]` next to wherever `[tigr_role]` is documented (grep).

- [ ] **Step 3: CHANGELOG**

Under `## [Unreleased]`:
- `### Breaking`: "`TigrRole` is now two-level (114 → ~140 nodes; `level` 0/1, new `level_kind`); `gene_count` on new mainrole nodes is a subtree count. `Gene_has_tigr_role` now spans all 43 organisms (~13.6K → ~27K edges); `sources` may be `[cyanorak, interproscan]` — filter with `'cyanorak' IN r.sources`, not list equality; `evidence` gains `family_inferred`. `gene_category` changes `Unknown` → a category on ~720 genes."
- `### Added`: hierarchy edge; `Ncbifam_family_has_tigr_role` bridge (TIGRFAMs 15.0 archive, ~1,720 edges); equivalog-gated inferred gene roles; `[tigr_role_inferred]` description lines; `TigrRole.ncbifam_family_count`; `tigr_roles.json` step-9 artefact; indexes; 3 mainrole junk flags.
- `### Highlights`: "Heterotroph genes (Alteromonas, Shewanella, Pseudomonas, …) can now be found by JCVI functional role — the same role vocabulary Prochlorococcus/Synechococcus already carried from Cyanorak — via NCBIfam equivalog families."

- [ ] **Step 4: Backlog + vocabulary contract**

- `plans/backlog.md` lines 27–44: retitle to "**TIGRFAM→TigrRole — SHIPPED 2026-08-29 as (a) ontology bridge + (b) equivalog-gated gene edges**", keep the 2026-08-18 measurement paragraph, add: "Gene-level rejection superseded after gating on `family_type = equivalog`: multi-role genes 375→15, contradictions 7%→5.4% (facet choices, not errors), `annotation_quality` movement 0 — see spec + `docs/kg-changes/tigr-role-bridge.md`." Lines 65–76 (hierarchy normalization): mark `[x]` done 2026-08-29. Add a new entry: "**Multi-source `Gene_has_ncbifam_family` — measured 2026-08-28, deferred.** Cyanorak `protein_domains` TIGR* tokens are 93% already in InterProScan hits (1,244 extra ids / 903 genes, Pro/Syn only); PGAP `inference` HMM ids are 84% NF* `domain`/`PfamEq` models InterPro's NCBIfam member DB excludes by design (0/1,860 PfamEq, 108/14,234 domain observed across 43 genomes; not a version gap — ids interleave). Would add corroboration only (~480 TIGR + ~205 NF PGAP confirmations); gate on `for_naming`/family-grade types if ever done."
- `plans/interpro_redesign_backlog.md:37-40`: mark `[x]` with pointer.
- `docs/kg-changes/vocabulary-contract.md`: add rows for `TigrRole.level_kind`, widened `Gene_has_tigr_role.sources` / `.evidence`, and note `Ncbifam_family_has_tigr_role` carries no properties (R3/R5).

- [ ] **Step 5: Commit**

```bash
git add docs/kg-changes/tigr-role-bridge.md CLAUDE.md CHANGELOG.md plans/backlog.md plans/interpro_redesign_backlog.md docs/kg-changes/vocabulary-contract.md
git commit -m "docs: TigrRole hierarchy + NCBIfam bridge + inferred roles (kg-changes note, CLAUDE.md, CHANGELOG, backlog)"
```

---

### Task 12: Regenerate merged annotations, full build, Docker rebuild, verification

**Files:**
- Regenerated (committed): `cache/data/*/genomes/*/gene_annotations_merged.json` (+ `gene_annotations_wide.json` if step 2 rewrites it)
- Regenerated: `tests/kg_validity/snapshot_data.json`, `tests/kg_validity/annotation_state_baseline.json`

- [ ] **Step 1: Baselines before the rebuild (against the currently deployed graph)**

```bash
uv run python tests/kg_validity/capture_annotation_state.py --save
bash scripts/post-import-validate.sh > logs/tigr_baseline_postimport.txt
```
Also run `/omics-edge-snapshot` (save mode) per its SKILL.md.

- [ ] **Step 2: Re-run step 2 for all strains (uses the new tigr_roles.json)**

Run: `bash scripts/prepare_data.sh --steps 2 --force` (takes a while; log `logs/prepare_data_step2.log`). Confirm the header line `TIGR roles archive: 116 roles, 2920 family links` appears and no `WARNING: … tigr_roles.json missing`.

Then measure the merge effect:
```bash
uv run python - <<'EOF'
import json, glob, collections
c = collections.Counter()
for f in glob.glob("cache/data/*/genomes/*/gene_annotations_merged.json"):
    genes = json.load(open(f)); genes = genes.get("genes", genes)
    for g in (genes.values() if isinstance(genes, dict) else genes):
        c["genes"] += 1
        if any(s.startswith("[tigr_role_inferred]") for s in g.get("alternate_functional_descriptions") or []):
            c["with_inferred_line"] += 1
print(c)
EOF
git diff --stat -- cache/data | tail -1
```
Expected: `with_inferred_line` ≈ 24K–25K (13.7K heterotroph + ~10–11K Cyanorak-strain genes where the inferred role differs from or adds to the curated line). Category shift: compare `cat_Unknown` totals in `logs/prepare_data_step2.log` against the previous run's log (expect ≈ −720).

- [ ] **Step 3: Full local build sanity (no Docker)**

Run: `uv run python create_knowledge_graph.py --output-dir ./biocypher-log/tigr_full/ 2>&1 | tee logs/tigr_full_build.log | grep -E "TigrRole|tigr-role-bridge|is_a edges|Traceback"`
Expected: `~140 TigrRole nodes`, `~110 TigrRole is_a edges`, `~1,720 tigr-role-bridge`, `Gene_has_tigr_role` total ≈ 27K in the adapter log, no Traceback.

- [ ] **Step 4: Commit regenerated caches, push, Docker rebuild in the other clone**

```bash
git add cache/data
git commit -m "data: regenerate gene_annotations_merged.json with TIGR-role inference (gene_category fill + [tigr_role_inferred])"
git push -u origin feat/tigrrole-hierarchy-ncbifam-bridge
```
Then in `multiomics_biocypher_kg` (the Docker clone — user's workflow): `git fetch && git checkout feat/tigrrole-hierarchy-ncbifam-bridge && git pull`, `bash scripts/prepare_data.sh --steps 9` (creates `tigr_roles.json` raw-free from the committed file — it already exists, so this is a no-op check), `docker compose down deploy app && docker compose up -d` and wait for `post-process` to exit 0 (check `[timing]` lines in `docker compose logs post-process`).

- [ ] **Step 5: Verify against the live graph**

```bash
uv run pytest -m kg -v
uv run python tests/kg_validity/capture_annotation_state.py --compare
bash scripts/post-import-validate.sh > logs/tigr_after_postimport.txt; diff logs/tigr_baseline_postimport.txt logs/tigr_after_postimport.txt | head -40
```
Expected: all KG tests pass (including Task 10's); `--compare` shows 0 (or ≤ a handful) `annotation_state` moves; the post-import diff shows ONLY TigrRole/Gene_has_tigr_role/`gene_category`-related lines. Run `/omics-edge-snapshot` compare — no lost expression edges.

Spot checks via `/cypher-queries` or `run_cypher`:
```cypher
MATCH (t:TigrRole) RETURN t.level_kind, count(*) ORDER BY t.level_kind;
MATCH (t:TigrRole {id:'tigr.role:energy_metabolism'}) RETURN t.name, t.gene_count, t.direct_gene_count, t.ncbifam_family_count;
MATCH (g:Gene {organism_name:'Alteromonas macleodii HOT1A3'})-[r:Gene_has_tigr_role]->(t) RETURN count(DISTINCT g), count(r);
MATCH (g:Gene)-[a:Gene_has_tigr_role]->(x),(g)-[b:Gene_has_tigr_role]->(y) WHERE a.evidence='curated' AND b.evidence='family_inferred' RETURN count(DISTINCT g);
```
Expected: 2 level_kind values; energy_metabolism gene_count > direct_gene_count; HOT1A3 ≈ 700 genes; disagreement genes ≈ 550.

- [ ] **Step 6: Regenerate the snapshot fixture and commit**

```bash
uv run python tests/kg_validity/generate_snapshot.py
git add tests/kg_validity/snapshot_data.json tests/kg_validity/annotation_state_baseline.json
git commit -m "test(kg): refresh snapshot + annotation-state baseline after TigrRole rebuild"
git push
```

- [ ] **Step 7: Finish the branch**

Use the `superpowers:finishing-a-development-branch` skill (merge to `main` per the repo's usual flow; the user rebuilds Docker from `main` in the other clone).
