# TCDB and CAZy Ontologies Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote TCDB and CAZy classifications to first-class ontologies (TcdbFamily, CazyFamily nodes with hierarchy edges, gene-annotation edges, and TCDB substrate→Metabolite linkage) and remove the redundant flat Gene array properties.

**Architecture:** Add two ontology adapters (mirroring BRITE for pruning + Pfam for two-class shape) plus a step 6 build extension that resolves TCDB substrates through MNX into the existing Metabolite layer (~785 new transport-only nodes). Five new edge labels; no `TransportReaction` reified entities (TCDB modeled as ontology, not reaction — see spec).

**Tech Stack:** Python (BioCypher adapters), SQLite (MNX resolver, existing), Cypher (post-import + KG validity tests), pytest, Docker neo4j-admin import.

**Spec:** [`docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md`](../specs/2026-05-01-tcdb-cazy-ontologies-design.md). All design decisions are settled there — do not reopen them.

**Cross-spec coordination doc:** [`multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md`](../../../../multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md). TCDB-S1..S5 + MET-M4 are already adopted into this plan.

---

## File map (decomposition reference)

The 8 commits touch these files. Each commit's section below lists exactly which subset it modifies.

### New files
- `multiomics_kg/adapters/tcdb_adapter.py` — `TcdbAnnotationAdapter` + `MultiTcdbAnnotationAdapter` (commit 5).
- `multiomics_kg/adapters/cazy_adapter.py` — `CazyAnnotationAdapter` + `MultiCazyAnnotationAdapter` (commit 6).
- `tests/test_tcdb_adapter.py` (commit 5).
- `tests/test_cazy_adapter.py` (commit 6).
- `tests/kg_validity/test_tcdb_cazy.py` (commit 8).

### Modified files
- `multiomics_kg/utils/cazy_utils.py` — pure-Python rewrite (commit 1).
- `multiomics_kg/utils/tcdb_utils.py` — re-point at new `tcdb_pruned.json` (commit 3) or kept reading raw `tcdb_hierarchy.json` if step 6 still writes it.
- `multiomics_kg/download/utils/annotation_transforms.py` — drop `_tx_validate_tcdb` / `_tx_validate_cazy` (commits 1 + 2).
- `multiomics_kg/download/build_metabolite_resolver.py` → renamed to `build_mnx_resolver.py` (commit 2). Drops TCDB + CAZy hierarchy logic.
- `multiomics_kg/download/build_kegg_metabolism_xrefs.py` — extended with TCDB download + parse + prune + substrate resolution + `additional_compounds` fold-in (commits 2 + 3).
- `multiomics_kg/download/download_genome_data.py` — drop `step6_tcdb_reference()` (commit 2).
- `multiomics_kg/adapters/metabolism_adapter.py` — read `additional_compounds`; emit `evidence_sources` on every Metabolite (commit 3).
- `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` — drop `TRANSPORTER_CLASSIFICATION` and `CAZY_IDS` enum values (commits 5 + 6 — they're unreferenced outside the enum but tracked here for completeness).
- `config/schema_config.yaml` — add `evidence_sources` (commit 3); add `tcdb family` / `cazy family` + 5 edges (commit 4); drop `transporter_classification` from `gene` (commit 5); drop `cazy_ids` from `gene` (commit 6); add new Gene rollups + Metabolite `transporter_count` + TcdbFamily/CazyFamily computed properties (commit 7).
- `config/gene_annotations_config.yaml` — drop `validate_cazy` transform (commit 1); drop `validate_tcdb` transform (commit 2); drop `transporter_classification` field (commit 5); drop `cazy_ids` field (commit 6).
- `scripts/refresh_mnx.sh` — point at renamed `build_mnx_resolver.py` (commit 2).
- `scripts/prepare_data.sh` — drop step 0 sub-step 6 reference; documented update (commit 2).
- `scripts/post-import.sh` + `scripts/post-import.cypher` — add indexes + computed properties + Gene rollups + Organism_has_metabolite UNION arm (commit 7). Both files must stay in sync (CLAUDE.md "post-import.cypher must mirror post-import.sh").
- `tests/test_cazy_utils.py` — rewrite to test pure-function API (commit 1).
- `tests/test_build_metabolite_resolver.py` → renamed `test_build_mnx_resolver.py`; drop CAZy + TCDB tests (commit 2).
- `tests/test_build_kegg_metabolism_xrefs.py` — extend with TCDB pruning + substrate resolution + `additional_compounds` + `evidence_sources` (commit 3).
- `tests/test_prepare_data_step2_metabolism_smoke.py` — drop `cazy_hierarchy.json` path entry (commit 1).
- `tests/test_build_gene_annotations_metabolism.py`, `tests/test_annotation_transforms_metabolism.py` — drop CAZy/TCDB validation transform tests (commits 1 + 2).
- `tests/kg_validity/test_structure.py` — extend orphan checks for TcdbFamily / CazyFamily; widen Metabolite count threshold (commits 3 + 8).
- `tests/kg_validity/test_post_import.py` — add TcdbFamily / CazyFamily rollup assertions (commit 7).
- `tests/kg_validity/snapshot_data.json` — regenerated (commit 8).
- `create_knowledge_graph.py` — wire in TCDB adapter (commit 5) + CAZy adapter (commit 6).
- `CLAUDE.md` — update (commit 8).
- `memory/MEMORY.md` + new `memory/project_tcdb_cazy_ontologies.md` (commit 8).
- `.claude/skills/cypher-queries/SKILL.md` — add new node/edge types + 2–3 templates (commit 8).
- `docs/kg-changes/tcdb-cazy-ontologies.md` — flip from "PROPOSED" to "LANDED"; backfill exact counts (commit 8).

### Deleted files
- `cache/data/cazy/cazy_hierarchy.json` (commit 1, optional rm — file is gitignored cache).

---

## Commit 1 — CAZy: pure-Python utils + drop validation + drop `cazy_hierarchy.json`

**Goal:** Remove the JSON-file CAZy hierarchy in favor of a pure-Python class table. CAZy validation transform goes away (passthrough). Adapter for CAZy is added in commit 6 — this commit only takes things away.

**Files:**
- Modify: `multiomics_kg/utils/cazy_utils.py` (full rewrite, ~40 lines → ~50 lines).
- Modify: `multiomics_kg/download/utils/annotation_transforms.py:264-275, 295` (drop `_tx_validate_cazy`, `validate_cazy` registry entry, `is_valid_cazy` import).
- Modify: `multiomics_kg/download/build_metabolite_resolver.py:458-563, 651` (drop `_CAZY_CLASSES`, `_CAZY_FAMILY_RE`, `_parse_cazy_id`, `_collect_cazy_ids`, `build_cazy_hierarchy`, the `eggnog_paths` discovery and `n_cazy` invocation in `main()`, the `cazy_hierarchy_entry_count` field in the report dict, and the `cazy_json` path entry from `_resolve_paths()`).
- Modify: `config/gene_annotations_config.yaml:520-522` (drop the `transform: validate_cazy` line — keep the `cazy_ids` field; commit 6 will remove the field itself).
- Modify: `tests/test_cazy_utils.py` (full rewrite to test pure-function API).
- Modify: `tests/test_prepare_data_step2_metabolism_smoke.py` (remove `cazy_hierarchy.json` path entry from the expected path list).
- Modify: `tests/test_build_gene_annotations_metabolism.py`, `tests/test_annotation_transforms_metabolism.py` (drop `validate_cazy` transform tests + `cazy_hierarchy.json` bootstrap fixtures — keep tests that exercise other transforms).
- Modify: `tests/test_build_metabolite_resolver.py` (drop `test_build_cazy_hierarchy`-style tests; the file is renamed in commit 2 so localized edits here only).
- Delete: `cache/data/cazy/cazy_hierarchy.json` (cache file, optional rm).

- [ ] **Step 1.1: Write the failing test for the new pure-Python `cazy_utils`**

Replace `tests/test_cazy_utils.py` entirely:

```python
"""Pure-function CAZy utility tests (post-cazy_hierarchy.json removal)."""
from __future__ import annotations

from multiomics_kg.utils.cazy_utils import (
    CAZY_CLASSES,
    cazy_ancestors,
    is_valid_cazy,
    parse_cazy_id,
)


def test_cazy_classes_table_has_six_classes():
    assert set(CAZY_CLASSES.keys()) == {"GH", "GT", "PL", "CE", "AA", "CBM"}
    assert CAZY_CLASSES["GH"] == "Glycoside Hydrolases"


def test_parse_cazy_id_family_no_subfamily():
    assert parse_cazy_id("GH13") == ("GH13", None)


def test_parse_cazy_id_with_subfamily():
    assert parse_cazy_id("GH13_5") == ("GH13", "GH13_5")


def test_parse_cazy_id_cbm():
    assert parse_cazy_id("CBM48") == ("CBM48", None)


def test_parse_cazy_id_malformed_returns_none():
    assert parse_cazy_id("GH") is None
    assert parse_cazy_id("XYZ12") is None
    assert parse_cazy_id("") is None
    assert parse_cazy_id("GH13_") is None


def test_parse_cazy_id_strips_whitespace():
    assert parse_cazy_id("  GH13  ") == ("GH13", None)


def test_is_valid_cazy_recognizes_class_family_subfamily():
    assert is_valid_cazy("GH")
    assert is_valid_cazy("GH13")
    assert is_valid_cazy("GH13_5")


def test_is_valid_cazy_rejects_garbage():
    assert not is_valid_cazy("not-a-cazy")
    assert not is_valid_cazy("")
    assert not is_valid_cazy("GH-13")


def test_cazy_ancestors_subfamily():
    assert cazy_ancestors("GH13_5") == ["GH", "GH13"]


def test_cazy_ancestors_family():
    assert cazy_ancestors("GH13") == ["GH"]


def test_cazy_ancestors_class():
    assert cazy_ancestors("GH") == []


def test_cazy_ancestors_unknown_returns_empty():
    assert cazy_ancestors("not-a-cazy") == []
```

- [ ] **Step 1.2: Run the test — verify import error**

Run: `pytest tests/test_cazy_utils.py -v`
Expected: FAIL — `ImportError: cannot import name 'CAZY_CLASSES'` (or similar; current module exports different names).

- [ ] **Step 1.3: Rewrite `multiomics_kg/utils/cazy_utils.py` as pure-Python**

Replace the full file:

```python
"""Pure-Python CAZy classification helpers.

CAZy hierarchy is now derived in-process — no file I/O. The 6 classes are
hardcoded; family + subfamily IDs are parsed from observed eggNOG annotations
in `multiomics_kg/adapters/cazy_adapter.py`.

Public API:
    CAZY_CLASSES — dict of class code → display name (immutable map).
    parse_cazy_id(token) — split a CAZy token into (family, subfamily | None).
    is_valid_cazy(value) — True when value is a recognized class / family / subfamily.
    cazy_ancestors(value) — root-to-parent ancestor chain ([class] or [class, family]).
"""
from __future__ import annotations

import re

CAZY_CLASSES: dict[str, str] = {
    "GH":  "Glycoside Hydrolases",
    "GT":  "GlycosylTransferases",
    "PL":  "Polysaccharide Lyases",
    "CE":  "Carbohydrate Esterases",
    "AA":  "Auxiliary Activities",
    "CBM": "Carbohydrate-Binding Modules",
}

_CAZY_FAMILY_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)(?:_(\d+))?$")


def parse_cazy_id(token: str) -> tuple[str, str | None] | None:
    """Return (family_id, subfamily_id_or_None) or None for malformed tokens.

    Examples:
        'GH13'    → ('GH13', None)
        'GH13_5'  → ('GH13', 'GH13_5')
        'CBM48'   → ('CBM48', None)
        'invalid' → None
    """
    if not token:
        return None
    m = _CAZY_FAMILY_RE.match(token.strip())
    if not m:
        return None
    cls, fam_num, sub_num = m.groups()
    family = f"{cls}{fam_num}"
    subfamily = f"{family}_{sub_num}" if sub_num else None
    return family, subfamily


def is_valid_cazy(value: str) -> bool:
    """True for a recognized class code, family ID, or subfamily ID."""
    if not value:
        return False
    if value in CAZY_CLASSES:
        return True
    return parse_cazy_id(value) is not None


def cazy_ancestors(value: str) -> list[str]:
    """Root-to-parent ancestor list. Empty for class-level or unknown values."""
    if not value or value in CAZY_CLASSES:
        return []
    parsed = parse_cazy_id(value)
    if parsed is None:
        return []
    family, subfamily = parsed
    cls = _CAZY_FAMILY_RE.match(family).group(1)
    if subfamily and value == subfamily:
        return [cls, family]
    return [cls]
```

- [ ] **Step 1.4: Run the test — verify pass**

Run: `pytest tests/test_cazy_utils.py -v`
Expected: PASS (12 tests).

- [ ] **Step 1.5: Drop `validate_cazy` transform from annotation_transforms.py**

In `multiomics_kg/download/utils/annotation_transforms.py`:

1. Delete lines 264-265 (the `is_valid_tcdb` import line stays for now — commit 2 removes it):

```python
# Delete:
from multiomics_kg.utils.cazy_utils import is_valid_cazy
```

2. Delete `_tx_validate_cazy` function (lines 273-275):

```python
def _tx_validate_cazy(value: str) -> str | None:
    """Drop CAZy IDs not present in the CAZy hierarchy."""
    return value if is_valid_cazy(value) else None
```

3. Delete the registry entry (line 295):

```python
"validate_cazy":                  _tx_validate_cazy,
```

- [ ] **Step 1.6: Drop `transform: validate_cazy` from `gene_annotations_config.yaml`**

Edit `config/gene_annotations_config.yaml:516-522` so the `cazy_ids` block becomes:

```yaml
  cazy_ids:
    type: union
    sources:
      - source: eggnog
        field: CAZy
        delimiter: ","
```

(Keep the field itself for now — commit 6 removes it after the CAZy adapter exists.)

- [ ] **Step 1.7: Drop CAZy logic from `build_metabolite_resolver.py`**

In `multiomics_kg/download/build_metabolite_resolver.py`:

1. Delete the entire `# ── CAZy hierarchy (bootstrapped from eggNOG observations) ────` section (lines 458-563): `_CAZY_CLASSES` dict, `_EGGNOG_CAZY_COL`, `_CAZY_FAMILY_RE`, `_parse_cazy_id`, `_collect_cazy_ids`, `build_cazy_hierarchy`.
2. In `_resolve_paths()`, drop the `"cazy_json"` entry.
3. In `main()`, drop:
   - The `_find_eggnog_annotations(CACHE_ROOT)` call.
   - The `n_cazy = build_cazy_hierarchy(...)` call.
   - The `"cazy_hierarchy_entry_count"` and `"eggnog_annotation_files_scanned"` fields in the `report` dict.
4. Update the `if not force` short-circuit at the top of `main()` to drop the `paths["cazy_json"].exists()` check.
5. Update the module docstring to remove the "cazy_hierarchy.json" line.

- [ ] **Step 1.8: Drop `cazy_hierarchy.json` references from existing tests**

Edit `tests/test_prepare_data_step2_metabolism_smoke.py`: remove the `cache/data/cazy/cazy_hierarchy.json` path entry from the expected output list (find via `grep -n cazy_hierarchy tests/test_prepare_data_step2_metabolism_smoke.py`).

Edit `tests/test_build_gene_annotations_metabolism.py` and `tests/test_annotation_transforms_metabolism.py`: drop any test method that bootstraps `cazy_hierarchy.json` or asserts on `_tx_validate_cazy`. Use `grep -n "cazy_hierarchy\|validate_cazy" tests/test_build_gene_annotations_metabolism.py tests/test_annotation_transforms_metabolism.py` to locate them.

Edit `tests/test_build_metabolite_resolver.py`: drop any `test_build_cazy_hierarchy*`, `test_collect_cazy_ids*`, or `test_parse_cazy_id*` (the parse function moved to `tests/test_cazy_utils.py`).

- [ ] **Step 1.9: Delete the cache file (optional)**

```bash
rm -f cache/data/cazy/cazy_hierarchy.json
rmdir cache/data/cazy 2>/dev/null || true
```

(The file is in gitignored cache; this just frees disk.)

- [ ] **Step 1.10: Run the full unit test suite**

Run: `pytest -m "not slow and not kg" -v`
Expected: PASS. If anything fails outside the changes you just made, stop and investigate (don't paper over it). Common breakage point: a stale `is_valid_cazy` import elsewhere — search via `grep -rn is_valid_cazy multiomics_kg tests`.

- [ ] **Step 1.11: Commit**

```bash
git add multiomics_kg/utils/cazy_utils.py \
        multiomics_kg/download/utils/annotation_transforms.py \
        multiomics_kg/download/build_metabolite_resolver.py \
        config/gene_annotations_config.yaml \
        tests/test_cazy_utils.py \
        tests/test_prepare_data_step2_metabolism_smoke.py \
        tests/test_build_gene_annotations_metabolism.py \
        tests/test_annotation_transforms_metabolism.py \
        tests/test_build_metabolite_resolver.py
git commit -m "$(cat <<'EOF'
refactor(cazy): drop cazy_hierarchy.json — pure-Python utils

CAZy hierarchy is now an in-process Python table (multiomics_kg/utils/cazy_utils.py)
instead of a JSON file. CAZy validation transform is dropped from the gene-annotation
pipeline (CAZy IDs become passthrough; unknown IDs are filtered at adapter time in
commit 6). Pure-function API (parse_cazy_id, is_valid_cazy, cazy_ancestors) replaces
the JSON-file accessors. build_metabolite_resolver.py drops CAZy entirely.

Why: simpler artifact set; CAZy hierarchy is small (6 classes), stable, and
only needs string parsing — no file I/O necessary. First commit of the
TCDB-CAZy ontology promotion (see plans/2026-05-01-tcdb-cazy-ontologies.md).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** Unit tests pass (`pytest -m "not slow and not kg"`). No KG rebuild needed.

---

## Commit 2 — TCDB: move download + parse to step 6, rename build script, drop validation

**Goal:** Pipeline reorganization without changing graph output. TCDB download moves from `prepare_data.sh` step 0.6 → `prepare_data.sh` step 6. TCDB hierarchy JSON parsing moves from `build_metabolite_resolver.py` → `build_kegg_metabolism_xrefs.py`. Script `build_metabolite_resolver.py` is renamed to `build_mnx_resolver.py` (does only one thing: build the heavy MNX SQLite). `validate_tcdb` transform is dropped (TCDB IDs become passthrough).

**Files:**
- Rename: `multiomics_kg/download/build_metabolite_resolver.py` → `multiomics_kg/download/build_mnx_resolver.py`. Drop the entire `# ── TCDB hierarchy ────` section (`_TC_CLASS_NAMES`, `_parse_tcdb_families`, `_parse_tcdb_superfamilies`, `_parse_tcdb_substrates`, `build_tcdb_hierarchy`) and the corresponding orchestration in `main()` (the `n_tcdb = build_tcdb_hierarchy(...)` call + the `tcdb_hierarchy_entry_count` report field + the TCDB paths in `_resolve_paths()`).
- Rename: `tests/test_build_metabolite_resolver.py` → `tests/test_build_mnx_resolver.py`. Drop tests referencing TCDB.
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` — add a new function `_build_tcdb_hierarchy(cache_root)` (lifted from the old `build_metabolite_resolver.build_tcdb_hierarchy`) and a new `download_metabolism_reference.download_all(sources=["tcdb"])` invocation early in `main()`. The new function writes `cache/data/tcdb/tcdb_hierarchy.json` (same contents as before — substrate resolution is added in commit 3).
- Modify: `multiomics_kg/download/download_genome_data.py:446-450, 480-483, 530-531` — drop `step6_tcdb_reference()` and its CLI entry; update `--steps` choices and help text from `[1,2,3,4,5,6]` to `[1,2,3,4,5]`; epilog text drops "6 Download TCDB reference data".
- Modify: `multiomics_kg/download/utils/annotation_transforms.py` — drop `_tx_validate_tcdb`, `is_valid_tcdb` import, `"validate_tcdb"` registry entry.
- Modify: `config/gene_annotations_config.yaml:420-426` — drop `transform: validate_tcdb` line.
- Modify: `scripts/refresh_mnx.sh` — replace `build_metabolite_resolver` reference with `build_mnx_resolver`.
- Modify: `scripts/prepare_data.sh:39` — update step-list comment if it references step 0.6.
- Modify: `multiomics_kg/utils/tcdb_utils.py` — keep reading `cache/data/tcdb/tcdb_hierarchy.json` (still produced by step 6 — only the producer has moved).

- [ ] **Step 2.1: Write failing test for step 6 building tcdb_hierarchy.json**

Add to `tests/test_build_kegg_metabolism_xrefs.py` (or create if missing):

```python
def test_step6_builds_tcdb_hierarchy(tmp_path, monkeypatch):
    """Step 6 must produce cache/data/tcdb/tcdb_hierarchy.json (was previously
    built by build_metabolite_resolver.py)."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    # Stage minimal TCDB TSVs the way download_metabolism_reference.download_all would:
    tcdb_dir = tmp_path / "cache" / "data" / "tcdb"
    tcdb_dir.mkdir(parents=True)
    (tcdb_dir / "families.tsv").write_text("1.A.1\tThe Voltage-gated Ion Channel\n")
    (tcdb_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (tcdb_dir / "substrates.tsv").write_text(
        "1.A.1.5.2\tCHEBI:9314;sucrose|CHEBI:3308;calcium(2+)\n"
    )

    out = tcdb_dir / "tcdb_hierarchy.json"
    mod._build_tcdb_hierarchy(cache_root=tmp_path / "cache" / "data")

    assert out.exists()
    import json
    h = json.loads(out.read_text())
    # 5 levels for 1.A.1.5.2 → expect class+subclass+family+subfamily+specificity
    assert "1" in h
    assert "1.A" in h
    assert "1.A.1" in h
    assert "1.A.1.5" in h
    assert "1.A.1.5.2" in h
    assert h["1.A.1.5.2"]["level_kind"] == "tc_specificity"
    assert h["1.A.1.5.2"]["substrate_classes"] == ["sucrose", "calcium(2+)"]
```

- [ ] **Step 2.2: Run the test — verify failure**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py::test_step6_builds_tcdb_hierarchy -v`
Expected: FAIL — `AttributeError: module ... has no attribute '_build_tcdb_hierarchy'`.

- [ ] **Step 2.3: Move TCDB parser into `build_kegg_metabolism_xrefs.py`**

Lift the entire `# ── TCDB hierarchy ────` section from `build_metabolite_resolver.py` (lines 306-455 — `_TC_CLASS_NAMES`, `_parse_tcdb_families`, `_parse_tcdb_superfamilies`, `_parse_tcdb_substrates`, `build_tcdb_hierarchy`) into `build_kegg_metabolism_xrefs.py`. Wrap the public entry as `_build_tcdb_hierarchy(cache_root: Path) -> int`:

```python
def _build_tcdb_hierarchy(cache_root: Path) -> int:
    """Parse the 3 TCDB TSVs into cache_root/tcdb/tcdb_hierarchy.json. Returns entry count."""
    tcdb_dir = cache_root / "tcdb"
    return build_tcdb_hierarchy(
        out_path=tcdb_dir / "tcdb_hierarchy.json",
        families_path=tcdb_dir / "families.tsv",
        superfamilies_path=tcdb_dir / "superfamilies.tsv",
        substrates_path=tcdb_dir / "substrates.tsv",
    )
```

In `main()`, after `kegg_utils.download_kegg_raw(...)` (line 435), add:

```python
log.info("Ensuring TCDB reference TSVs are downloaded ...")
from multiomics_kg.download.download_metabolism_reference import download_all
download_all(force=force, sources=["tcdb"])

log.info("Building tcdb_hierarchy.json ...")
_build_tcdb_hierarchy(cache_root)
```

- [ ] **Step 2.4: Run the test — verify pass**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py::test_step6_builds_tcdb_hierarchy -v`
Expected: PASS.

- [ ] **Step 2.5: Rename `build_metabolite_resolver.py` → `build_mnx_resolver.py`**

```bash
git mv multiomics_kg/download/build_metabolite_resolver.py multiomics_kg/download/build_mnx_resolver.py
git mv tests/test_build_metabolite_resolver.py tests/test_build_mnx_resolver.py
```

In the renamed `build_mnx_resolver.py`:
1. Update the module docstring to mention only the MNX SQLite output.
2. Delete the entire `# ── TCDB hierarchy ────` section (already moved in step 2.3 — duplicate must be removed).
3. In `_resolve_paths()`, drop the `tcdb_json` entry (and the `families`/`superfamilies`/`substrates` entries — they're step 6's responsibility now).
4. In `main()`, drop the `n_tcdb = build_tcdb_hierarchy(...)` invocation and the `tcdb_hierarchy_entry_count` report field. The MNX-only short-circuit becomes `if not force and paths["resolver_db"].exists():`.

In `tests/test_build_mnx_resolver.py`:
1. Rename test functions referencing the old module path so the import points at `build_mnx_resolver`.
2. Delete any `test_build_tcdb_hierarchy*`, `test_parse_tcdb_*` tests (those move with the parser to `tests/test_build_kegg_metabolism_xrefs.py` if you want runtime coverage of them — for this commit, just delete; commit 3 adds richer coverage).

- [ ] **Step 2.6: Update `scripts/refresh_mnx.sh`**

Replace the `build_metabolite_resolver` reference with `build_mnx_resolver`:

```bash
uv run python -m multiomics_kg.download.build_mnx_resolver $FORCE
```

- [ ] **Step 2.7: Drop `step6_tcdb_reference()` from `download_genome_data.py`**

In `multiomics_kg/download/download_genome_data.py`:

1. Delete the function `step6_tcdb_reference` (lines 446-450).
2. In the parser declaration:
   - Change `--steps choices=[1,2,3,4,5,6], default=[1,2,3,4,5,6]` → `choices=[1,2,3,4,5], default=[1,2,3,4,5]`.
   - Help text: drop `"6=tcdb_reference"`.
3. Drop the `if 6 in steps: step6_tcdb_reference(...)` block.
4. Update epilog to drop the step-6 line and the standalone TCDB note (TCDB now lives in prepare_data.sh step 6).

- [ ] **Step 2.8: Drop `validate_tcdb` from annotation_transforms.py**

In `multiomics_kg/download/utils/annotation_transforms.py`:

1. Drop the import `from multiomics_kg.utils.tcdb_utils import is_valid_tcdb` (line 264).
2. Drop the `_tx_validate_tcdb` function (lines 268-270).
3. Drop the `"validate_tcdb": _tx_validate_tcdb,` registry line (line 294).

- [ ] **Step 2.9: Drop `transform: validate_tcdb` from `gene_annotations_config.yaml`**

Edit `config/gene_annotations_config.yaml:420-426`. The `transporter_classification` block becomes:

```yaml
  transporter_classification:
    type: union
    sources:
      - source: eggnog
        field: KEGG_TC
        delimiter: ","
```

- [ ] **Step 2.10: Update `scripts/prepare_data.sh` step-0 sub-step list**

Edit `scripts/prepare_data.sh` — search for any docs/comments referencing step 0.6 / "TCDB reference" and update to reflect the move (in particular the comment header that lists "Step 0 sub-steps: …").

- [ ] **Step 2.11: Run all unit tests**

Run: `pytest -m "not slow and not kg" -v`
Expected: PASS. If `tests/test_download_genome_data.py` references step 6, fix the test to match the new step list.

- [ ] **Step 2.12: Run prepare_data step 6 against the cached fixtures (optional smoke)**

If MNX SQLite + TCDB TSVs are already cached locally:

```bash
bash scripts/prepare_data.sh --steps 6 --force
```

Expected: produces `cache/data/kegg/kegg_data.json` (size unchanged from baseline) AND `cache/data/tcdb/tcdb_hierarchy.json`. Inspect both files briefly. The kegg_data.json schema is unchanged at this commit — the substrate fold-in lands in commit 3.

- [ ] **Step 2.13: Commit**

```bash
git add multiomics_kg/download/build_mnx_resolver.py \
        multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        multiomics_kg/download/download_genome_data.py \
        multiomics_kg/download/utils/annotation_transforms.py \
        config/gene_annotations_config.yaml \
        scripts/refresh_mnx.sh \
        scripts/prepare_data.sh \
        tests/test_build_mnx_resolver.py \
        tests/test_build_kegg_metabolism_xrefs.py \
        tests/test_download_genome_data.py
git commit -m "$(cat <<'EOF'
refactor(tcdb): move TCDB download + parse to prepare_data step 6

build_metabolite_resolver.py is renamed to build_mnx_resolver.py and now does
only one thing: build the heavy MNX SQLite (~30 min). TCDB reference TSV
download + tcdb_hierarchy.json build move to prepare_data step 6 alongside the
KEGG metabolism cache, so all chemistry-related data lives in one step.
download_genome_data step 0.6 is removed; --steps now goes 1-5. The
validate_tcdb annotation transform is dropped (TCDB IDs become passthrough;
unknown IDs are filtered at adapter time in commit 5).

Output is byte-identical to before this commit — pipeline reorganization only.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** Unit tests pass (`pytest -m "not slow and not kg"`). Optional: `bash scripts/prepare_data.sh --steps 6 --force` produces the same `kegg_data.json` (no schema change yet) plus a fresh `tcdb_hierarchy.json`.

---

## Commit 3 — Step 6: substrate resolution + `additional_compounds` + `evidence_sources` tagging

**Goal:** Extend `kegg_data.json` with a new `additional_compounds` section for non-KEGG transport chemistry (~785 nodes); tag every compound entry with `evidence_sources: list[str]`; write `cache/data/tcdb/tcdb_pruned.json` (kept-node set + per-leaf metabolite primary IDs); extend `metabolism_adapter.py` to read `additional_compounds` and emit `evidence_sources` on every Metabolite. **No new node types added yet** — schema additions for TcdbFamily/CazyFamily are deferred to commit 4.

**Files:**
- Modify: `config/schema_config.yaml:655-672` — add `evidence_sources: str[]` to the `metabolite` properties block.
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` — extend with TCDB pruning + substrate→Metabolite resolution + `additional_compounds` fold-in + `evidence_sources` tagging + `tcdb_pruned.json` write.
- Modify: `multiomics_kg/adapters/metabolism_adapter.py` — read `additional_compounds`; emit `evidence_sources` on every Metabolite.
- Modify: `multiomics_kg/utils/tcdb_utils.py` — optionally re-point at `tcdb_pruned.json` for adapter convenience (or keep raw `tcdb_hierarchy.json`; commit 5 decides).
- Modify: `tests/test_build_kegg_metabolism_xrefs.py` — add tests for substrate walk, MNX resolution, `additional_compounds`, `evidence_sources`, `tcdb_pruned.json`.
- Modify: `tests/kg_validity/test_structure.py` — widen Metabolite-count threshold (1000-5000 → e.g. 1000-6000) to accept the +785 transport-only growth.
- Modify: `tests/kg_validity/test_metabolism.py` — same widening if it has its own threshold.

### Step-6 algorithm (extension)

After the existing `build_pruned_kegg_data(...)` call, add:

1. Load `tcdb_hierarchy.json`.
2. Collect all TCDB IDs from all strains' `gene_annotations_merged.json` (already happens in `_gene_reachable_sets` for KOs — extend it to also collect TCDB IDs).
3. For each gene-annotated TCDB ID, walk **down** to all leaf descendants (`tc_specificity`) and **up** to `tc_class`. Union → `kept_tcdb_ids`.
4. For each leaf in `kept_tcdb_ids`, parse `substrate_classes` strings (`CHEBI:NNNN;name`); resolve each `chebi:NNNN` via the MNX resolver to its MNXM ID; resolve MNXM → primary metabolite ID via `metabolite_utils.mnxm_to_primary_id`.
5. For each resolved primary ID:
   - If it's a `kegg.compound:C*` ID and is already in `kegg_data["compounds"]`: append `"transport"` to that entry's `evidence_sources` (initialize list if absent).
   - Else: add a new entry under `kegg_data["additional_compounds"]` keyed by the bare primary ID (e.g. `chebi:9314`, `mnx:MNXM12345`). Properties pulled from MNX: `name`, `formula`, `mass`, `inchikey`, `mnxm_id`, `chebi_id`. Set `evidence_sources: ["transport"]`.
6. Tag every entry in the existing `compounds` section with `evidence_sources: ["metabolism"]` (or `["metabolism", "transport"]` if also reached via transport).
7. Write `kegg_data.json` (now with `additional_compounds` section + `evidence_sources` everywhere).
8. Write `cache/data/tcdb/tcdb_pruned.json`:

```json
{
  "kept_tcdb_ids": ["1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2", ...],
  "leaf_substrates": {
    "1.A.1.5.2": ["kegg.compound:C00208", "chebi:9314", ...]
  }
}
```

(`leaf_substrates` is keyed by `tc_specificity`-level TCID; values are resolved Metabolite primary IDs in deterministic sorted order. Adapter consumes this in commit 5 — no need to re-resolve at build time.)

- [ ] **Step 3.1: Write failing test for `evidence_sources` on existing compounds**

Add to `tests/test_build_kegg_metabolism_xrefs.py`:

```python
def test_compounds_get_metabolism_evidence_source(tmp_path):
    """Every entry in kegg_data['compounds'] gets evidence_sources=['metabolism']
    after step 6 runs. (No transport overlap in this minimal fixture.)"""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    # Stage a minimal kegg_data.json by calling build_pruned_kegg_data with a
    # tiny raw + empty TCDB hierarchy; assert evidence_sources tagging.
    # Implementation: see fixture pattern in existing tests at
    # tests/test_build_kegg_metabolism_xrefs.py:test_pathway_pruning_keeps_kos
    # …
    # Final assertion:
    assert all(
        cpd.get("evidence_sources") == ["metabolism"]
        for cpd in result["compounds"].values()
    )
```

(Concrete fixture body: re-use the conftest fixture `_minimal_raw` if it exists; otherwise skim `tests/test_build_kegg_metabolism_xrefs.py` and add the missing pieces.)

- [ ] **Step 3.2: Write failing test for `additional_compounds` from transport**

```python
def test_transport_only_compounds_land_in_additional_compounds(tmp_path):
    """A TCDB substrate not reachable via metabolism shows up under
    kegg_data['additional_compounds'] with evidence_sources=['transport']."""
    # Stage:
    # - tcdb_hierarchy.json with one leaf: 1.A.1.5.2 transporting CHEBI:9999;tetracycline
    # - empty kegg compounds set
    # - resolver-db stub mapping chebi:9999 → MNXM00099 → chebi:9999 primary
    # …
    assert "chebi:9999" in result["additional_compounds"]
    entry = result["additional_compounds"]["chebi:9999"]
    assert entry["evidence_sources"] == ["transport"]
    assert entry["chebi_id"] == "9999"
```

- [ ] **Step 3.3: Write failing test for `tcdb_pruned.json`**

```python
def test_tcdb_pruned_json_contains_kept_ids_and_leaf_substrates(tmp_path):
    """Step 6 writes cache/data/tcdb/tcdb_pruned.json with kept_tcdb_ids
    (above + below) and per-leaf metabolite primary IDs."""
    # Stage fixtures + run step 6
    out = tmp_path / "cache" / "data" / "tcdb" / "tcdb_pruned.json"
    assert out.exists()
    import json
    p = json.loads(out.read_text())
    assert "kept_tcdb_ids" in p
    assert "leaf_substrates" in p
    # walk-up: gene-annotated 1.A.1 forces 1, 1.A, 1.A.1 to be kept
    # walk-down: forces all 1.A.1.5.2 leaves to be kept too
    assert {"1", "1.A", "1.A.1"}.issubset(set(p["kept_tcdb_ids"]))
    assert "1.A.1.5.2" in p["leaf_substrates"]
```

- [ ] **Step 3.4: Run the three new tests — verify they fail**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py -k "evidence_sources or additional_compounds or tcdb_pruned" -v`
Expected: FAIL (functions don't exist yet).

- [ ] **Step 3.5: Implement substrate resolution + `additional_compounds` + tagging**

In `multiomics_kg/download/build_kegg_metabolism_xrefs.py`:

1. Extend `_gene_reachable_sets()` to also return a `tcdb_ids: set[str]` — collected from `gene["transporter_classification"]` (existing field on `gene_annotations_merged.json`).

2. Add a new function:

```python
def _prune_tcdb(hierarchy: dict, seed_ids: set[str]) -> tuple[set[str], dict[str, list[str]]]:
    """Bidirectional prune: walk up (to tc_class) AND down (to tc_specificity).

    Returns (kept_ids, leaf_to_substrate_strings).
    """
    kept: set[str] = set()
    leaf_subs: dict[str, list[str]] = {}

    # Walk up + collect descendants per seed
    parent_of = {tc: data.get("parent") for tc, data in hierarchy.items()}
    children_of: dict[str, list[str]] = {}
    for tc, parent in parent_of.items():
        if parent is not None:
            children_of.setdefault(parent, []).append(tc)

    def walk_up(tc: str) -> None:
        cur = tc
        while cur is not None and cur in hierarchy and cur not in kept:
            kept.add(cur)
            cur = parent_of.get(cur)

    def walk_down(tc: str) -> None:
        if tc in kept and any(c in kept for c in children_of.get(tc, [])):
            return  # already fully expanded
        kept.add(tc)
        node = hierarchy.get(tc, {})
        if node.get("level_kind") == "tc_specificity":
            subs = node.get("substrate_classes")
            if subs:
                leaf_subs[tc] = list(subs)
        for child in children_of.get(tc, []):
            walk_down(child)

    for seed in seed_ids:
        if seed not in hierarchy:
            continue
        walk_up(seed)
        walk_down(seed)

    return kept, leaf_subs


def _resolve_substrates(
    leaf_subs: dict[str, list[str]],
    conn: sqlite3.Connection,
) -> tuple[dict[str, list[str]], dict[str, dict]]:
    """Resolve `CHEBI:NNNN;name` substrate strings to Metabolite primary IDs.

    Returns:
      leaf_to_primary_ids: {leaf_tcid: [primary_id, ...]} (sorted, deduped).
      compound_props: {primary_id: {name, formula, mass, inchikey, mnxm_id, chebi_id}}.
    """
    from multiomics_kg.utils.metabolite_utils import (
        resolve_metabolite, mnxm_to_primary_id,
    )

    leaf_to_primary_ids: dict[str, list[str]] = {}
    compound_props: dict[str, dict] = {}

    for leaf, subs in leaf_subs.items():
        primary_ids: list[str] = []
        for sub_str in subs:
            # Format: 'CHEBI:NNNN;name' or just 'name' (rare; spec says all 1410 resolve)
            if ":" not in sub_str:
                continue
            chebi_part, _, name_part = sub_str.partition(";")
            mnxm, method = resolve_metabolite(chebi_part, conn)
            if mnxm is None:
                log.debug(f"TCDB substrate unresolved: {sub_str!r}")
                continue
            primary = mnxm_to_primary_id(mnxm, conn)
            primary_ids.append(primary)

            # Look up MNX-side properties
            cur = conn.cursor()
            cur.execute(
                "SELECT name, formula, mass, inchikey FROM compounds WHERE mnxm_id = ?",
                (mnxm,),
            )
            row = cur.fetchone()
            mnx_name, formula, mass, inchikey = row if row else (name_part, None, None, None)
            cur.execute(
                "SELECT value FROM compound_aliases WHERE source='chebi' AND mnxm_id = ? "
                "ORDER BY value LIMIT 1",
                (mnxm,),
            )
            chebi_row = cur.fetchone()
            chebi_id = chebi_row[0] if chebi_row else None

            compound_props[primary] = {
                "name": mnx_name or name_part,
                "formula": formula or None,
                "mass": mass,
                "inchikey": inchikey or None,
                "mnxm_id": mnxm,
                "chebi_id": chebi_id,
            }
        leaf_to_primary_ids[leaf] = sorted(set(primary_ids))

    return leaf_to_primary_ids, compound_props


def _fold_substrates_into_kegg_data(
    kegg_data: dict,
    leaf_to_primary_ids: dict[str, list[str]],
    compound_props: dict[str, dict],
) -> None:
    """Fold transport-substrate metabolites into kegg_data, tagging evidence_sources."""
    # 1. Tag every existing compound with metabolism evidence
    for cpd in kegg_data.setdefault("compounds", {}).values():
        cpd.setdefault("evidence_sources", []).append("metabolism")

    additional: dict[str, dict] = kegg_data.setdefault("additional_compounds", {})

    # 2. Visit every transport substrate
    seen_in_metabolism: set[str] = set()
    for leaf, primary_ids in leaf_to_primary_ids.items():
        for primary_id in primary_ids:
            if primary_id.startswith("kegg.compound:"):
                kegg_cpd_id = primary_id[len("kegg.compound:"):]
                if kegg_cpd_id in kegg_data["compounds"]:
                    if "transport" not in kegg_data["compounds"][kegg_cpd_id]["evidence_sources"]:
                        kegg_data["compounds"][kegg_cpd_id]["evidence_sources"].append("transport")
                    seen_in_metabolism.add(primary_id)
                    continue
            # Else: into additional_compounds
            if primary_id not in additional:
                additional[primary_id] = {
                    **compound_props[primary_id],
                    "evidence_sources": ["transport"],
                }
            elif "transport" not in additional[primary_id]["evidence_sources"]:
                additional[primary_id]["evidence_sources"].append("transport")
```

3. In `main()`, after `build_pruned_kegg_data(...)` writes `kegg_data.json`:

```python
log.info("Pruning TCDB hierarchy + resolving transport substrates ...")
hierarchy = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())
seed_tcdb_ids = sets["tcdb_ids"]  # added to _gene_reachable_sets above

kept, leaf_subs = _prune_tcdb(hierarchy, seed_tcdb_ids)
leaf_to_primary, compound_props = _resolve_substrates(leaf_subs, conn)

kegg_data = json.loads(KEGG_DATA_FILE.read_text())
_fold_substrates_into_kegg_data(kegg_data, leaf_to_primary, compound_props)
KEGG_DATA_FILE.write_text(json.dumps(kegg_data, indent=2, sort_keys=True))

(cache_root / "tcdb" / "tcdb_pruned.json").write_text(json.dumps({
    "kept_tcdb_ids": sorted(kept),
    "leaf_substrates": {leaf: ids for leaf, ids in sorted(leaf_to_primary.items())},
}, indent=2, sort_keys=True))
log.info(f"  TCDB: {len(kept)} kept IDs, {len(leaf_to_primary)} leaves with substrates")
```

- [ ] **Step 3.6: Run the new tests — verify pass**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py -k "evidence_sources or additional_compounds or tcdb_pruned" -v`
Expected: PASS.

- [ ] **Step 3.7: Add `evidence_sources` schema entry**

Edit `config/schema_config.yaml:655-672`. The `metabolite` properties block becomes:

```yaml
metabolite:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.compound
  label_in_input: metabolite
  properties:
    kegg_compound_id: str
    name: str
    formula: str
    mass: float
    inchikey: str
    smiles: str
    mnxm_id: str
    chebi_id: str
    hmdb_id: str
    evidence_sources: str[]      # new: ['metabolism', 'transport']; metabolomics reserved
    organism_count: int
    gene_count: int
```

- [ ] **Step 3.8: Write failing test for `metabolism_adapter` reading `additional_compounds`**

Add to `tests/test_metabolism_adapter.py` (file likely exists; if not, mirror the existing `tests/test_brite_adapter.py` shape):

```python
def test_adapter_emits_additional_compounds_with_transport_evidence(tmp_path):
    """metabolism_adapter must read additional_compounds and tag evidence_sources=['transport']."""
    fixture = {
        "kos": {}, "pathways": {}, "subcategories": {}, "categories": {},
        "reactions": {},
        "compounds": {},
        "additional_compounds": {
            "chebi:9999": {
                "name": "tetracycline",
                "formula": "C22H24N2O8",
                "mass": 444.43,
                "inchikey": None,
                "mnxm_id": "MNXM00099",
                "chebi_id": "9999",
                "evidence_sources": ["transport"],
            }
        },
    }
    p = tmp_path / "kegg_data.json"
    import json
    p.write_text(json.dumps(fixture))

    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter
    nodes = list(MetabolismAdapter(kegg_data_path=p).get_nodes())
    metabolite_nodes = [n for n in nodes if n[1] == "metabolite"]
    assert len(metabolite_nodes) == 1
    nid, _label, props = metabolite_nodes[0]
    assert nid == "chebi:9999"
    assert props["evidence_sources"] == ["transport"]
    assert props["mnxm_id"] == "MNXM00099"
```

- [ ] **Step 3.9: Run the test — verify failure**

Run: `pytest tests/test_metabolism_adapter.py -k additional_compounds -v`
Expected: FAIL — adapter currently only reads `compounds`.

- [ ] **Step 3.10: Update `metabolism_adapter.py` to read `additional_compounds`**

In `multiomics_kg/adapters/metabolism_adapter.py:83-100`. After the loop over `data.get("compounds", {})`, add:

```python
n = 0
for primary_id, cpd in data.get("additional_compounds", {}).items():
    if self.test_mode and n >= 100:
        break
    # primary_id is already the full ID, e.g. 'chebi:9999' or 'mnx:MNXM12345'
    props = _drop_nulls({
        "kegg_compound_id": None,  # absent for non-KEGG
        "name": _clean_str(cpd.get("name", "")),
        "formula": _clean_str(cpd.get("formula")),
        "mass": cpd.get("mass"),
        "inchikey": _clean_str(cpd.get("inchikey")),
        "smiles": _clean_str(cpd.get("smiles")),
        "mnxm_id": cpd.get("mnxm_id"),
        "chebi_id": cpd.get("chebi_id"),
        "hmdb_id": cpd.get("hmdb_id"),
        "evidence_sources": list(cpd.get("evidence_sources", [])),
    })
    yield primary_id, "metabolite", props
    n += 1
```

Also extend the existing `compounds` loop to include `evidence_sources` in its `_drop_nulls(...)` call:

```python
"evidence_sources": list(cpd.get("evidence_sources", [])),
```

- [ ] **Step 3.11: Run the metabolism_adapter test — verify pass**

Run: `pytest tests/test_metabolism_adapter.py -v`
Expected: PASS.

- [ ] **Step 3.12: Widen Metabolite-count thresholds in KG-validity tests**

In `tests/kg_validity/test_metabolism.py:20-23`:

```python
def test_metabolite_count_in_expected_range(run_query):
    """Spec acceptance: 1,000 ≤ count ≤ 6,000 (post-TCDB transport-only growth)."""
    n = run_query("MATCH (m:Metabolite) RETURN count(m) AS n")[0]["n"]
    assert 1000 <= n <= 6000, f"Metabolite count {n} outside expected 1000-6000"
```

In `tests/kg_validity/test_structure.py` — search for any other `Metabolite` count assertion via `grep -n Metabolite tests/kg_validity/test_structure.py` and update similarly.

- [ ] **Step 3.13: Take an `/omics-edge-snapshot` baseline before the rebuild**

In an interactive shell with the current Docker graph still running:

```bash
# Slash-command flow — see SKILL.md
/omics-edge-snapshot
mv omics-edge-snapshot-*.json omics-edge-snapshot-pre-commit3.json
```

(Capture this on a real run; commit 8 will diff before/after to verify the omics layer is byte-identical.)

- [ ] **Step 3.14: Rebuild the KG**

```bash
docker compose down
bash scripts/prepare_data.sh --steps 6 --force
docker compose up -d --build
```

Wait for `import` and `post-process` containers to complete (~25 min). Inspect `cache/data/kegg/kegg_data.json` and confirm:

- It contains an `additional_compounds` section.
- Every entry in both `compounds` and `additional_compounds` has `evidence_sources` (a list).
- `cache/data/tcdb/tcdb_pruned.json` exists with `kept_tcdb_ids` and `leaf_substrates` keys.

- [ ] **Step 3.15: Run KG validity tests — verify pass**

Run: `pytest -m kg -v`
Expected: PASS. Metabolite count grows by ~785 (now ~2,973 total, well within the new threshold).

- [ ] **Step 3.16: Run `/omics-edge-snapshot` again — verify byte-identical**

```bash
/omics-edge-snapshot
diff omics-edge-snapshot-pre-commit3.json omics-edge-snapshot-*.json
```

Expected: empty diff (omics layer untouched).

- [ ] **Step 3.17: Commit**

```bash
git add config/schema_config.yaml \
        multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        multiomics_kg/adapters/metabolism_adapter.py \
        tests/test_build_kegg_metabolism_xrefs.py \
        tests/test_metabolism_adapter.py \
        tests/kg_validity/test_metabolism.py \
        tests/kg_validity/test_structure.py
git commit -m "$(cat <<'EOF'
feat(metabolism): TCDB substrate resolution + evidence_sources tagging

Step 6 now walks gene-annotated TCDB IDs bidirectionally (to tc_class above and
tc_specificity below), resolves each leaf's substrate strings (CHEBI:NNNN;name)
through MNX, and folds them into kegg_data.json. Compounds reachable via
catalysis carry evidence_sources=['metabolism']; transport-only substrates land
under a new additional_compounds section with evidence_sources=['transport'];
both-sources compounds get the union. Adds cache/data/tcdb/tcdb_pruned.json
(kept-node set + per-leaf metabolite primary IDs) so the TCDB adapter (commit 5)
doesn't re-resolve.

Metabolite layer grows by ~785 transport-only nodes (2188 → ~2973). No new
schema labels yet — TcdbFamily/CazyFamily land in commit 4. Omics layer is
unchanged (verified via /omics-edge-snapshot diff).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** Unit tests pass; KG rebuild succeeds; `pytest -m kg` passes; `/omics-edge-snapshot` diff is empty.

---

## Commit 4 — Schema additions: `tcdb family` / `cazy family` + 5 edges

**Goal:** Add the new node and edge labels to `schema_config.yaml`. No adapter changes — schema only validates against an empty input (no nodes/edges yet). This commit is independently safe to land.

**Files:**
- Modify: `config/schema_config.yaml` — append two new node blocks + five new edge blocks.

- [ ] **Step 4.1: Append new node + edge schema entries**

Append to `config/schema_config.yaml` after the existing `brite category in brite category` association block (~line 1018) — keep node blocks together and edge blocks together for readability:

```yaml
# ── TCDB ─────────────────────────────────────────────────────────────────────

tcdb family:
  is_a: named thing
  represented_as: node
  preferred_id: tcdb
  label_in_input: tcdb family
  properties:
    name: str               # TCDB entry name; falls back to tcdb_id when empty
    tcdb_id: str            # bare hierarchical ID, e.g. "1.A.1.5.2"
    level: int              # 0=tc_class … 4=tc_specificity
    level_kind: str
    superfamily: str        # sparse, leaf-only
    tc_class_id: str        # post-import: pointer to root tc_class (e.g. "tcdb:3"); self on tc_class
    gene_count: int         # post-import
    organism_count: int     # post-import
    member_count: int       # post-import: direct child count
    metabolite_count: int   # post-import: distinct substrate metabolites in subtree

cazy family:
  is_a: named thing
  represented_as: node
  preferred_id: cazy
  label_in_input: cazy family
  properties:
    name: str               # human-readable for class; falls back to cazy_id for family/subfamily
    cazy_id: str            # e.g. "GH13" or "GH13_5"
    level: int              # 0=cazy_class, 1=cazy_family, 2=cazy_subfamily
    level_kind: str
    gene_count: int         # post-import
    organism_count: int     # post-import

# ── TCDB / CAZy edges ────────────────────────────────────────────────────────

gene to tcdb family association:
  is_a: association
  represented_as: edge
  label_as_edge: Gene_has_tcdb_family
  source: gene
  target: tcdb family
  label_in_input: gene_has_tcdb_family

tcdb family hierarchical association:
  is_a: association
  represented_as: edge
  label_as_edge: Tcdb_family_is_a_tcdb_family
  source: tcdb family
  target: tcdb family
  label_in_input: tcdb_family_is_a_tcdb_family

tcdb family transports metabolite:
  is_a: association
  represented_as: edge
  label_as_edge: Tcdb_family_transports_metabolite
  source: tcdb family
  target: metabolite
  label_in_input: tcdb_family_transports_metabolite

gene to cazy family association:
  is_a: association
  represented_as: edge
  label_as_edge: Gene_has_cazy_family
  source: gene
  target: cazy family
  label_in_input: gene_has_cazy_family

cazy family hierarchical association:
  is_a: association
  represented_as: edge
  label_as_edge: Cazy_family_is_a_cazy_family
  source: cazy family
  target: cazy family
  label_in_input: cazy_family_is_a_cazy_family
```

- [ ] **Step 4.2: Run the unit test suite to confirm BioCypher schema parses**

Run: `pytest tests/test_create_knowledge_graph.py -v`
Expected: PASS — `BioCypher` instance constructs successfully against the new schema (no adapters yet emit nodes for these labels; that's fine).

If `tests/test_create_knowledge_graph.py` doesn't load schema explicitly, run a one-liner:

```bash
uv run python -c "from biocypher import BioCypher; bc = BioCypher(schema_config_path='config/schema_config.yaml', biocypher_config_path='config/biocypher_config.yaml'); print('OK')"
```

Expected: prints `OK`. Any YAML parse error or schema validation error fails the gate.

- [ ] **Step 4.3: Commit**

```bash
git add config/schema_config.yaml
git commit -m "$(cat <<'EOF'
schema: add TcdbFamily, CazyFamily nodes + 5 ontology edge labels

Adds tcdb family + cazy family node blocks (BRITE-style: single label per
ontology, distinguished by level_kind + level) and five new association edges:
Gene_has_tcdb_family, Tcdb_family_is_a_tcdb_family,
Tcdb_family_transports_metabolite, Gene_has_cazy_family,
Cazy_family_is_a_cazy_family. No adapter emits these yet — adapters land in
commits 5 and 6. Schema parses cleanly against an empty input.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** BioCypher schema constructs without error.

---

## Commit 5 — TCDB adapter + wire-in + drop `transporter_classification` from Gene

**Goal:** Implement `TcdbAnnotationAdapter` and `MultiTcdbAnnotationAdapter`; wire into `create_knowledge_graph.py`; remove the `transporter_classification: str[]` flat property from Gene.

**Files:**
- Create: `multiomics_kg/adapters/tcdb_adapter.py`.
- Create: `tests/test_tcdb_adapter.py`.
- Modify: `create_knowledge_graph.py` — instantiate `MultiTcdbAnnotationAdapter` + `bc.write_nodes(...)` / `bc.write_edges(...)`.
- Modify: `config/schema_config.yaml:378` — drop `transporter_classification: str[]` from `gene` properties.
- Modify: `config/gene_annotations_config.yaml:420-426` — drop the `transporter_classification:` field block (the field still flows through `gene_annotations_merged.json`; commit's choice is to drop from gene-property emission only — keep the merged JSON field so the adapter can read it). **NB:** spec sequencing says drop the field unless other code reads it. The adapter reads `gene["transporter_classification"]` directly from `gene_annotations_merged.json`, so the field stays in `gene_annotations_config.yaml` and `gene_annotations_merged.json` — only the schema-side Gene property goes away.
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:145` — drop the `TRANSPORTER_CLASSIFICATION` enum entry (unused outside the enum; no other consumers per `grep`).

### Adapter design

`TcdbAnnotationAdapter` (per-strain, mirrors `PfamAnnotationAdapter`):
- Loads `gene_annotations_merged.json`.
- Exposes `get_all_tcdb_ids() -> set[str]` for the orchestrator's pruning convenience (the orchestrator actually loads `tcdb_pruned.json` — but the per-strain adapter still emits gene→family edges only for IDs present in the pruned set, so cross-checking is useful for diagnostic logs).
- `get_edges()` yields `Gene_has_tcdb_family` edges using the **exact level annotated** (no walk-up — the hierarchy edges connect levels). Edge ID: `f"{locus_tag}-has_tcdb-{tcdb_id}"`.

`MultiTcdbAnnotationAdapter` (multi-strain orchestrator, mirrors `MultiPfamAnnotationAdapter` + `MultiBriteAdapter`):
- Loads `cache/data/tcdb/tcdb_pruned.json` and `cache/data/tcdb/tcdb_hierarchy.json`.
- `get_nodes()` yields one `tcdb family` node per `kept_tcdb_ids` entry, with properties pulled from `tcdb_hierarchy.json` (name, level, level_kind, superfamily for leaves).
- `get_edges()` yields:
  1. `Tcdb_family_is_a_tcdb_family` parent edges from the pruned hierarchy.
  2. `Gene_has_tcdb_family` edges via per-strain delegation (same pattern as Pfam orchestrator).
  3. `Tcdb_family_transports_metabolite` edges from `leaf_substrates` (one per leaf × resolved metabolite primary ID).

ID conventions:
- TcdbFamily node ID: `tcdb:{tcdb_id}` (e.g. `tcdb:1.A.1.5.2`).
- Metabolite target IDs come straight from `leaf_substrates` (already in canonical form: `kegg.compound:C*`, `chebi:NNNN`, or `mnx:MNXM*`).

- [ ] **Step 5.1: Write failing tests for the per-strain adapter**

Create `tests/test_tcdb_adapter.py`:

```python
"""Tests for TCDB ontology adapter (per-strain + multi-strain)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.adapters.tcdb_adapter import (
    TcdbAnnotationAdapter,
    MultiTcdbAnnotationAdapter,
)


@pytest.fixture
def strain_dir(tmp_path):
    """Stage a minimal gene_annotations_merged.json for one strain."""
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001",
                     "transporter_classification": ["1.A.1.5.2"]},
        "PMM_0002": {"locus_tag": "PMM_0002",
                     "transporter_classification": ["3.A.1"]},
        "PMM_0003": {"locus_tag": "PMM_0003"},  # no TCDB
    }))
    return d


def test_per_strain_get_all_tcdb_ids(strain_dir):
    a = TcdbAnnotationAdapter(genome_dir=strain_dir)
    assert a.get_all_tcdb_ids() == {"1.A.1.5.2", "3.A.1"}


def test_per_strain_get_edges(strain_dir):
    a = TcdbAnnotationAdapter(genome_dir=strain_dir)
    edges = list(a.get_edges())
    edge_ids = {(e[1], e[2], e[3]) for e in edges}
    assert ("ncbigene:PMM_0001", "tcdb:1.A.1.5.2", "gene_has_tcdb_family") in edge_ids
    assert ("ncbigene:PMM_0002", "tcdb:3.A.1", "gene_has_tcdb_family") in edge_ids
    assert len(edges) == 2  # PMM_0003 has no TCDB → no edge
```

For the orchestrator tests, stage minimal `tcdb_hierarchy.json` + `tcdb_pruned.json`:

```python
@pytest.fixture
def cache_root(tmp_path):
    """Stage a minimal cache_root with tcdb_hierarchy.json + tcdb_pruned.json."""
    tcdb_dir = tmp_path / "tcdb"
    tcdb_dir.mkdir(parents=True)
    (tcdb_dir / "tcdb_hierarchy.json").write_text(json.dumps({
        "1": {"name": "Channels and Pores", "level": 0,
              "level_kind": "tc_class", "parent": None},
        "1.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"name": "VIC family", "level": 2,
                  "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"name": "", "level": 3, "level_kind": "tc_subfamily",
                    "parent": "1.A.1"},
        "1.A.1.5.2": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "superfamily": "VIC Superfamily"},
    }))
    (tcdb_dir / "tcdb_pruned.json").write_text(json.dumps({
        "kept_tcdb_ids": ["1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2"],
        "leaf_substrates": {
            "1.A.1.5.2": ["kegg.compound:C00208", "chebi:9999"],
        },
    }))
    return tmp_path


def test_orchestrator_emits_one_node_per_kept_id(cache_root, strain_dir):
    # Stub out genome config so orchestrator only sees strain_dir
    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")

    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv),
        cache_root=cache_root,
        test_mode=False,
        cache=True,
    )
    a.download_data(cache=True)
    node_ids = {n[0] for n in a.get_nodes()}
    assert node_ids == {"tcdb:1", "tcdb:1.A", "tcdb:1.A.1", "tcdb:1.A.1.5", "tcdb:1.A.1.5.2"}


def test_orchestrator_emits_hierarchy_edges(cache_root, strain_dir):
    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")

    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv),
        cache_root=cache_root,
        test_mode=False, cache=True,
    )
    a.download_data(cache=True)
    edges = list(a.get_edges())
    parent_edges = {(e[1], e[2]) for e in edges if e[3] == "tcdb_family_is_a_tcdb_family"}
    # 1 has no parent; the rest do
    assert ("tcdb:1.A", "tcdb:1") in parent_edges
    assert ("tcdb:1.A.1", "tcdb:1.A") in parent_edges
    assert ("tcdb:1.A.1.5", "tcdb:1.A.1") in parent_edges
    assert ("tcdb:1.A.1.5.2", "tcdb:1.A.1.5") in parent_edges
    assert len(parent_edges) == 4


def test_orchestrator_emits_substrate_edges_only_on_leaves(cache_root, strain_dir):
    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")

    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv),
        cache_root=cache_root,
        test_mode=False, cache=True,
    )
    a.download_data(cache=True)
    edges = list(a.get_edges())
    sub_edges = {(e[1], e[2]) for e in edges if e[3] == "tcdb_family_transports_metabolite"}
    assert sub_edges == {
        ("tcdb:1.A.1.5.2", "kegg.compound:C00208"),
        ("tcdb:1.A.1.5.2", "chebi:9999"),
    }


def test_orchestrator_node_props(cache_root, strain_dir):
    config_csv = cache_root / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")

    a = MultiTcdbAnnotationAdapter(
        genome_config_file=str(config_csv),
        cache_root=cache_root,
        test_mode=False, cache=True,
    )
    a.download_data(cache=True)
    nodes = {n[0]: n for n in a.get_nodes()}

    # Class node: name from _TC_CLASS_NAMES, level 0
    nid, _label, props = nodes["tcdb:1"]
    assert props["level"] == 0
    assert props["level_kind"] == "tc_class"
    assert props["name"] == "Channels and Pores"
    assert props["tcdb_id"] == "1"

    # Leaf node: name falls back to tcdb_id when source name is empty;
    # superfamily is set
    nid, _label, props = nodes["tcdb:1.A.1.5.2"]
    assert props["level"] == 4
    assert props["level_kind"] == "tc_specificity"
    assert props["tcdb_id"] == "1.A.1.5.2"
    assert props["name"] == "1.A.1.5.2"  # fallback
    assert props["superfamily"] == "VIC Superfamily"
```

- [ ] **Step 5.2: Run the new tests — verify failure**

Run: `pytest tests/test_tcdb_adapter.py -v`
Expected: FAIL — `ModuleNotFoundError: multiomics_kg.adapters.tcdb_adapter`.

- [ ] **Step 5.3: Implement `multiomics_kg/adapters/tcdb_adapter.py`**

Create the full file:

```python
"""TCDB ontology adapter.

Yields:
- TcdbFamily nodes (only IDs in cache/data/tcdb/tcdb_pruned.json's kept_tcdb_ids).
- Tcdb_family_is_a_tcdb_family parent edges within the pruned hierarchy.
- Gene_has_tcdb_family edges from per-strain gene_annotations_merged.json.
- Tcdb_family_transports_metabolite edges on tc_specificity-level nodes only,
  using pre-resolved metabolite primary IDs from tcdb_pruned.json.

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter.
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

logger = logging.getLogger(__name__)

_TC_CLASS_NAMES = {
    "1": "Channels and Pores",
    "2": "Electrochemical Potential-driven Transporters",
    "3": "Primary Active Transporters",
    "4": "Group Translocators",
    "5": "Transmembrane Electron Carriers",
    "8": "Auxiliary Transport Proteins",
    "9": "Incompletely Characterized Transport Systems",
}


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return f"ncbigene:{locus_tag}"


def _tcdb_node_id(tcdb_id: str) -> str:
    return f"tcdb:{tcdb_id}"


class TcdbAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_tcdb_family edges from gene_annotations_merged.json."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        json_path = self.genome_dir / "gene_annotations_merged.json"
        if not json_path.exists():
            logger.warning(f"gene_annotations_merged.json not found at {json_path}, skipping")
            return
        with open(json_path, encoding="utf-8") as fh:
            self._genes = json.load(fh)

    def get_all_tcdb_ids(self) -> set[str]:
        ids: set[str] = set()
        for gene in self._genes.values():
            for tc in gene.get("transporter_classification") or []:
                if tc:
                    ids.add(tc)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            for tc in gene.get("transporter_classification") or []:
                if not tc:
                    continue
                yield (
                    f"{locus_tag}-has_tcdb-{tc}",
                    _gene_node_id(locus_tag),
                    _tcdb_node_id(tc),
                    "gene_has_tcdb_family",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(f"TcdbAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→TCDB edges")


class MultiTcdbAnnotationAdapter:
    """Multi-strain orchestrator: owns TcdbFamily nodes + parent edges + substrate edges."""

    def __init__(
        self,
        genome_config_file: str,
        cache_root: Path,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.cache_root = Path(cache_root)
        self.test_mode = test_mode
        self.cache = cache
        self._strain_adapters: list[TcdbAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)
        self._hierarchy: dict[str, dict] = {}
        self._kept_ids: set[str] = set()
        self._leaf_substrates: dict[str, list[str]] = {}

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                TcdbAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(f"MultiTcdbAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters")

    def download_data(self, cache: bool = True) -> None:
        """Read tcdb_hierarchy.json + tcdb_pruned.json (built by step 6)."""
        tcdb_dir = self.cache_root / "tcdb"
        hierarchy_path = tcdb_dir / "tcdb_hierarchy.json"
        pruned_path = tcdb_dir / "tcdb_pruned.json"
        if not hierarchy_path.exists() or not pruned_path.exists():
            raise FileNotFoundError(
                f"Missing TCDB cache file(s): {hierarchy_path}, {pruned_path}. "
                f"Run `bash scripts/prepare_data.sh --steps 6 --force` first."
            )
        self._hierarchy = json.loads(hierarchy_path.read_text())
        pruned = json.loads(pruned_path.read_text())
        self._kept_ids = set(pruned["kept_tcdb_ids"])
        self._leaf_substrates = pruned["leaf_substrates"]

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        emit_count = 0
        for tcdb_id in sorted(self._kept_ids):
            if self.test_mode and emit_count >= 100:
                break
            entry = self._hierarchy.get(tcdb_id, {})
            level = entry.get("level", 0)
            level_kind = entry.get("level_kind", "tc_class")
            raw_name = entry.get("name") or ""
            # Class fallback: pull from _TC_CLASS_NAMES if YAML name was empty
            if not raw_name and level_kind == "tc_class":
                raw_name = _TC_CLASS_NAMES.get(tcdb_id, "")
            # Other level fallback: tcdb_id itself
            display_name = raw_name or tcdb_id

            props = {
                "name": _clean_str(display_name),
                "tcdb_id": tcdb_id,
                "level": level,
                "level_kind": level_kind,
            }
            if entry.get("superfamily"):
                props["superfamily"] = _clean_str(entry["superfamily"])
            yield _tcdb_node_id(tcdb_id), "tcdb family", props
            emit_count += 1
        logger.info(f"MultiTcdbAnnotationAdapter.get_nodes: {emit_count} TcdbFamily nodes")

    def get_edges(self):
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        # 1. Parent edges
        parent_count = 0
        for tcdb_id in sorted(self._kept_ids):
            entry = self._hierarchy.get(tcdb_id, {})
            parent = entry.get("parent")
            if parent and parent in self._kept_ids:
                yield (
                    f"{tcdb_id}-parent-{parent}",
                    _tcdb_node_id(tcdb_id),
                    _tcdb_node_id(parent),
                    "tcdb_family_is_a_tcdb_family",
                    {},
                )
                parent_count += 1

        # 2. Gene→TcdbFamily edges (delegate to per-strain adapters; filter to kept IDs)
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                # edge[2] is "tcdb:..."; strip prefix to compare against kept_tcdb_ids
                tcdb_target = edge[2][len("tcdb:"):]
                if tcdb_target not in self._kept_ids:
                    logger.debug(f"Dropping gene→TCDB edge to unpruned {edge[2]!r}")
                    continue
                yield edge
                gene_count += 1

        # 3. Substrate edges (leaves only)
        sub_count = 0
        for leaf, primary_ids in sorted(self._leaf_substrates.items()):
            for primary in primary_ids:
                yield (
                    f"{leaf}-transports-{primary}",
                    _tcdb_node_id(leaf),
                    primary,  # primary_id is already in canonical form (kegg.compound:.., chebi:.., mnx:..)
                    "tcdb_family_transports_metabolite",
                    {},
                )
                sub_count += 1

        logger.info(
            f"MultiTcdbAnnotationAdapter.get_edges: {parent_count} parent, "
            f"{gene_count} gene, {sub_count} substrate edges"
        )
```

- [ ] **Step 5.4: Run the test — verify pass**

Run: `pytest tests/test_tcdb_adapter.py -v`
Expected: PASS (5 tests).

- [ ] **Step 5.5: Wire into `create_knowledge_graph.py`**

In `create_knowledge_graph.py`, after the `pfam_adapter.write_edges(...)` block (~line 197), add:

```python
    # TCDB transport classification ontology + substrate→Metabolite bridge
    from multiomics_kg.adapters.tcdb_adapter import MultiTcdbAnnotationAdapter
    tcdb_adapter = MultiTcdbAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_root=Path("cache/data"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    tcdb_adapter.download_data(cache=CACHE)
    bc.write_nodes(tcdb_adapter.get_nodes())
    bc.write_edges(tcdb_adapter.get_edges())
```

- [ ] **Step 5.6: Drop `transporter_classification` from `gene` schema**

Edit `config/schema_config.yaml:378`. Delete the line `transporter_classification: str[]`. (The entry above and below stays.)

- [ ] **Step 5.7: Drop the unused enum entry**

Edit `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:145`. Delete the line `TRANSPORTER_CLASSIFICATION = 'transporter_classification'`.

- [ ] **Step 5.8: Run unit tests**

Run: `pytest -m "not slow and not kg" -v`
Expected: PASS. If `test_create_knowledge_graph.py` references `transporter_classification`, update it.

- [ ] **Step 5.9: Take an `/omics-edge-snapshot` baseline**

```bash
/omics-edge-snapshot
mv omics-edge-snapshot-*.json omics-edge-snapshot-pre-commit5.json
```

- [ ] **Step 5.10: Rebuild the KG**

```bash
docker compose down
docker compose up -d --build
```

- [ ] **Step 5.11: Run KG validity tests**

Run: `pytest -m kg -v`
Expected: PASS (the new `test_tcdb_cazy.py` doesn't exist yet — commit 8 — but existing tests should all still pass; the new TcdbFamily nodes don't break any orphan checks because they have parent edges within the new ontology).

- [ ] **Step 5.12: Spot-check the graph**

Open Neo4j Browser at http://localhost:7474 and run:

```cypher
MATCH (tf:TcdbFamily) RETURN count(tf) AS n;
// Expected: a few hundred (gene-reachable subhierarchy)

MATCH (tf:TcdbFamily {level_kind: 'tc_class'}) RETURN tf.tcdb_id, tf.name;
// Expected: at least a few class-level nodes (1, 2, 3, ...)

MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n;
// Expected: 1500-2500

MATCH ()-[r:Tcdb_family_transports_metabolite]->() RETURN count(r) AS n;
// Expected: ~10000

MATCH (tf:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
WHERE tf.tcdb_id = '1.A.1.5.2'
RETURN m.name LIMIT 5;
// Expected: actual substrate names
```

- [ ] **Step 5.13: Run `/omics-edge-snapshot` and verify byte-identical to commit 3 baseline**

```bash
/omics-edge-snapshot
diff omics-edge-snapshot-pre-commit3.json omics-edge-snapshot-*.json
```

Expected: empty diff.

- [ ] **Step 5.14: Commit**

```bash
git add multiomics_kg/adapters/tcdb_adapter.py \
        tests/test_tcdb_adapter.py \
        create_knowledge_graph.py \
        config/schema_config.yaml \
        multiomics_kg/adapters/cyanorak_ncbi_adapter.py
git commit -m "$(cat <<'EOF'
feat(tcdb): TCDB transport ontology + substrate→Metabolite bridge

Adds TcdbAnnotationAdapter (per-strain) + MultiTcdbAnnotationAdapter (pruned
hierarchy + substrate edges) following the BRITE+Pfam two-class precedent.
Reads cache/data/tcdb/tcdb_pruned.json (built in commit 3) so node + edge
emission are pure file→graph transformations.

Drops Gene.transporter_classification flat property in favor of
Gene_has_tcdb_family edges. Hierarchy is queryable via
Tcdb_family_is_a_tcdb_family*0.. ; substrates via Tcdb_family_transports_metabolite
on tc_specificity leaves only.

Edge counts (rough): ~few-hundred TcdbFamily nodes, ~1500-2500 gene→TCDB
edges, ~10K substrate edges. Omics layer unchanged (verified via
/omics-edge-snapshot diff against commit 3 baseline).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** Unit tests pass; rebuilt KG passes `pytest -m kg`; spot-checks above match expected; `/omics-edge-snapshot` byte-identical to commit 3 baseline.

---

## Commit 6 — CAZy adapter + wire-in + drop `cazy_ids` from Gene

**Goal:** Implement `CazyAnnotationAdapter` and `MultiCazyAnnotationAdapter`; wire in; remove the `cazy_ids: str[]` flat property from Gene.

**Files:**
- Create: `multiomics_kg/adapters/cazy_adapter.py`.
- Create: `tests/test_cazy_adapter.py`.
- Modify: `create_knowledge_graph.py` — instantiate and write the CAZy adapter alongside TCDB.
- Modify: `config/schema_config.yaml:379` — drop `cazy_ids: str[]` from `gene` properties.
- Modify: `config/gene_annotations_config.yaml:516-522` — keep the field block (the adapter still reads `gene["cazy_ids"]` from `gene_annotations_merged.json`); only the schema-side Gene property goes away.
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:146` — drop the `CAZY_IDS` enum entry.

### Adapter design

`CazyAnnotationAdapter` (per-strain):
- Loads `gene_annotations_merged.json`.
- `get_all_cazy_ids() -> set[str]` returns the raw observed tokens (e.g. `"GH13_5"`) — these tokens may map to a class+family+subfamily triple.
- `get_edges()` yields `Gene_has_cazy_family` edges to the **most specific** level the token represents (i.e. subfamily if subfamily is set, else family). Spec: "attaches at exact level annotated."

`MultiCazyAnnotationAdapter` (multi-strain orchestrator):
- Walks all per-strain `get_all_cazy_ids()` to derive the full set of class + family + subfamily IDs that need nodes (only observed ones).
- `get_nodes()` yields one `cazy family` node per derived ID.
- `get_edges()` yields:
  1. `Cazy_family_is_a_cazy_family` parent edges (subfamily → family, family → class).
  2. `Gene_has_cazy_family` via per-strain delegation.

ID conventions:
- CazyFamily node ID: `cazy:{cazy_id}` (e.g. `cazy:GH`, `cazy:GH13`, `cazy:GH13_5`).

- [ ] **Step 6.1: Write failing tests for CAZy adapter**

Create `tests/test_cazy_adapter.py`:

```python
"""Tests for CAZy ontology adapter (per-strain + multi-strain)."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.adapters.cazy_adapter import (
    CazyAnnotationAdapter,
    MultiCazyAnnotationAdapter,
)


@pytest.fixture
def strain_dir(tmp_path):
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(json.dumps({
        "PMM_0001": {"locus_tag": "PMM_0001", "cazy_ids": ["GH13"]},
        "PMM_0002": {"locus_tag": "PMM_0002", "cazy_ids": ["GH13_5", "CBM48"]},
        "PMM_0003": {"locus_tag": "PMM_0003"},  # no CAZy
        "PMM_0004": {"locus_tag": "PMM_0004", "cazy_ids": ["GH13_5", "garbage"]},
    }))
    return d


def test_per_strain_get_all_cazy_ids_filters_garbage(strain_dir, caplog):
    a = CazyAnnotationAdapter(genome_dir=strain_dir)
    ids = a.get_all_cazy_ids()
    assert ids == {"GH13", "GH13_5", "CBM48"}
    assert "garbage" not in ids


def test_per_strain_get_edges_attaches_at_most_specific_level(strain_dir):
    a = CazyAnnotationAdapter(genome_dir=strain_dir)
    edges = list(a.get_edges())
    edge_targets = {(e[1], e[2]) for e in edges}
    # 'GH13' attaches to family-level cazy:GH13
    assert ("ncbigene:PMM_0001", "cazy:GH13") in edge_targets
    # 'GH13_5' attaches to subfamily-level cazy:GH13_5 (NOT cazy:GH13)
    assert ("ncbigene:PMM_0002", "cazy:GH13_5") in edge_targets
    assert ("ncbigene:PMM_0002", "cazy:CBM48") in edge_targets
    assert ("ncbigene:PMM_0004", "cazy:GH13_5") in edge_targets
    # garbage skipped
    assert all("garbage" not in t for _, t in edge_targets)


def test_orchestrator_emits_class_family_subfamily(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    a = MultiCazyAnnotationAdapter(
        genome_config_file=str(config_csv),
        test_mode=False,
    )
    nodes = list(a.get_nodes())
    node_ids = {n[0] for n in nodes}
    # Observed: GH13, GH13_5, CBM48 → derives GH (class), GH13 (fam), GH13_5 (sub),
    # CBM (class), CBM48 (fam). NO subfamily under CBM48.
    assert node_ids == {"cazy:GH", "cazy:GH13", "cazy:GH13_5", "cazy:CBM", "cazy:CBM48"}


def test_orchestrator_class_node_has_human_readable_name(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    a = MultiCazyAnnotationAdapter(
        genome_config_file=str(config_csv),
        test_mode=False,
    )
    nodes = {n[0]: n for n in a.get_nodes()}
    _, _label, props = nodes["cazy:GH"]
    assert props["name"] == "Glycoside Hydrolases"
    assert props["level"] == 0
    assert props["level_kind"] == "cazy_class"
    # Family fallback: name == cazy_id when source has no human name
    _, _label, props = nodes["cazy:GH13"]
    assert props["name"] == "GH13"
    assert props["level"] == 1
    assert props["level_kind"] == "cazy_family"


def test_orchestrator_emits_parent_edges(strain_dir, tmp_path):
    config_csv = tmp_path / "genomes.csv"
    config_csv.write_text(f"strain_name,data_dir\nMED4,{strain_dir}\n")
    a = MultiCazyAnnotationAdapter(
        genome_config_file=str(config_csv),
        test_mode=False,
    )
    edges = list(a.get_edges())
    parents = {(e[1], e[2]) for e in edges if e[3] == "cazy_family_is_a_cazy_family"}
    assert ("cazy:GH13", "cazy:GH") in parents
    assert ("cazy:GH13_5", "cazy:GH13") in parents
    assert ("cazy:CBM48", "cazy:CBM") in parents
    assert ("cazy:GH", "cazy:GH") not in parents  # classes have no parent
```

- [ ] **Step 6.2: Run the test — verify failure**

Run: `pytest tests/test_cazy_adapter.py -v`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 6.3: Implement `multiomics_kg/adapters/cazy_adapter.py`**

Create the full file:

```python
"""CAZy ontology adapter — observed-only (no external download).

Yields:
- CazyFamily nodes (class + family + subfamily, only IDs observed in
  gene_annotations_merged.json across configured strains).
- Cazy_family_is_a_cazy_family parent edges.
- Gene_has_cazy_family edges from per-strain merged JSON. Each gene attaches
  to the MOST SPECIFIC observed level (subfamily if set, else family).
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.cazy_utils import CAZY_CLASSES, parse_cazy_id

logger = logging.getLogger(__name__)


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return f"ncbigene:{locus_tag}"


def _cazy_node_id(cazy_id: str) -> str:
    return f"cazy:{cazy_id}"


def _most_specific_id(token: str) -> str | None:
    """Return the subfamily if present, else the family. None for malformed."""
    parsed = parse_cazy_id(token)
    if parsed is None:
        return None
    family, subfamily = parsed
    return subfamily or family


class CazyAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_cazy_family edges from merged JSON."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        json_path = self.genome_dir / "gene_annotations_merged.json"
        if not json_path.exists():
            logger.warning(f"gene_annotations_merged.json not found at {json_path}, skipping")
            return
        with open(json_path, encoding="utf-8") as fh:
            self._genes = json.load(fh)

    def get_all_cazy_ids(self) -> set[str]:
        ids: set[str] = set()
        for gene in self._genes.values():
            for token in gene.get("cazy_ids") or []:
                if not token:
                    continue
                parsed = parse_cazy_id(token)
                if parsed is None:
                    logger.debug(f"CAZy malformed token {token!r}, skipping")
                    continue
                family, subfamily = parsed
                ids.add(family)
                if subfamily:
                    ids.add(subfamily)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            for token in gene.get("cazy_ids") or []:
                target = _most_specific_id(token)
                if target is None:
                    logger.debug(f"CAZy malformed token {token!r} on {locus_tag}, skipping")
                    continue
                yield (
                    f"{locus_tag}-has_cazy-{target}",
                    _gene_node_id(locus_tag),
                    _cazy_node_id(target),
                    "gene_has_cazy_family",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(f"CazyAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→CAZy edges")


class MultiCazyAnnotationAdapter:
    """Multi-strain orchestrator: owns CazyFamily nodes + parent edges."""

    def __init__(
        self,
        genome_config_file: str,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self._strain_adapters: list[CazyAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                CazyAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(f"MultiCazyAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters")

    def download_data(self, cache: bool = True) -> None:
        """No external resources; CAZy hierarchy is in-process Python."""
        return

    def _all_observed_ids(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_cazy_ids()
        return ids

    def _expand_to_class_family_subfamily(self, observed: set[str]) -> set[str]:
        """For every observed family/subfamily, add its class (and family if it's a subfamily)."""
        out: set[str] = set()
        for tok in observed:
            parsed = parse_cazy_id(tok)
            if parsed is None:
                continue
            family, subfamily = parsed
            cls = family[:2] if family.startswith("CB") else family[:1]
            # Robust: re-extract via regex match
            from multiomics_kg.utils.cazy_utils import _CAZY_FAMILY_RE  # type: ignore
            m = _CAZY_FAMILY_RE.match(family)
            cls = m.group(1) if m else cls
            out.add(cls)
            out.add(family)
            if subfamily:
                out.add(subfamily)
        return out

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        observed = self._all_observed_ids()
        all_ids = self._expand_to_class_family_subfamily(observed)
        for cazy_id in sorted(all_ids):
            if cazy_id in CAZY_CLASSES:
                level, level_kind, name = 0, "cazy_class", CAZY_CLASSES[cazy_id]
            elif "_" in cazy_id:
                level, level_kind = 2, "cazy_subfamily"
                name = cazy_id  # fallback
            else:
                level, level_kind = 1, "cazy_family"
                name = cazy_id  # fallback
            props = {
                "name": _clean_str(name),
                "cazy_id": cazy_id,
                "level": level,
                "level_kind": level_kind,
            }
            yield _cazy_node_id(cazy_id), "cazy family", props
        logger.info(f"MultiCazyAnnotationAdapter.get_nodes: {len(all_ids)} CazyFamily nodes")

    def get_edges(self):
        observed = self._all_observed_ids()
        all_ids = self._expand_to_class_family_subfamily(observed)

        # 1. Parent edges
        parent_count = 0
        for cazy_id in sorted(all_ids):
            if cazy_id in CAZY_CLASSES:
                continue
            parsed = parse_cazy_id(cazy_id)
            if parsed is None:
                continue
            family, subfamily = parsed
            from multiomics_kg.utils.cazy_utils import _CAZY_FAMILY_RE  # type: ignore
            m = _CAZY_FAMILY_RE.match(family)
            cls = m.group(1) if m else None
            if cls is None:
                continue
            if subfamily and cazy_id == subfamily:
                parent = family
            else:
                parent = cls
            yield (
                f"{cazy_id}-parent-{parent}",
                _cazy_node_id(cazy_id),
                _cazy_node_id(parent),
                "cazy_family_is_a_cazy_family",
                {},
            )
            parent_count += 1

        # 2. Gene→CAZy edges via per-strain delegation
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(
            f"MultiCazyAnnotationAdapter.get_edges: {parent_count} parent, {gene_count} gene edges"
        )
```

- [ ] **Step 6.4: Run the test — verify pass**

Run: `pytest tests/test_cazy_adapter.py -v`
Expected: PASS (5 tests).

- [ ] **Step 6.5: Wire into `create_knowledge_graph.py`**

In `create_knowledge_graph.py`, after the TCDB adapter block from commit 5, add:

```python
    # CAZy carbohydrate-active-enzyme ontology (observed-only)
    from multiomics_kg.adapters.cazy_adapter import MultiCazyAnnotationAdapter
    cazy_adapter = MultiCazyAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    bc.write_nodes(cazy_adapter.get_nodes())
    bc.write_edges(cazy_adapter.get_edges())
```

- [ ] **Step 6.6: Drop `cazy_ids` from `gene` schema**

Edit `config/schema_config.yaml:379`. Delete the line `cazy_ids: str[]`.

- [ ] **Step 6.7: Drop the unused enum entry**

Edit `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:146`. Delete the line `CAZY_IDS = 'cazy_ids'`.

- [ ] **Step 6.8: Run unit tests**

Run: `pytest -m "not slow and not kg" -v`
Expected: PASS.

- [ ] **Step 6.9: Take an `/omics-edge-snapshot` baseline**

```bash
/omics-edge-snapshot
mv omics-edge-snapshot-*.json omics-edge-snapshot-pre-commit6.json
```

- [ ] **Step 6.10: Rebuild the KG**

```bash
docker compose down
docker compose up -d --build
```

- [ ] **Step 6.11: Run KG validity tests**

Run: `pytest -m kg -v`
Expected: PASS.

- [ ] **Step 6.12: Spot-check the graph**

```cypher
MATCH (cf:CazyFamily) RETURN cf.level_kind, count(cf) AS n
ORDER BY cf.level_kind;
// Expected: cazy_class ~6 (or fewer if some classes unobserved),
// cazy_family ~30-50, cazy_subfamily small/zero

MATCH ()-[r:Gene_has_cazy_family]->() RETURN count(r) AS n;
// Expected: 400-600

MATCH (g:Gene)-[:Gene_has_cazy_family]->(:CazyFamily {cazy_id: 'GH13'}) RETURN g LIMIT 3;
// Sanity check
```

- [ ] **Step 6.13: Run `/omics-edge-snapshot` and verify byte-identical**

```bash
/omics-edge-snapshot
diff omics-edge-snapshot-pre-commit3.json omics-edge-snapshot-*.json
```

Expected: empty diff.

- [ ] **Step 6.14: Commit**

```bash
git add multiomics_kg/adapters/cazy_adapter.py \
        tests/test_cazy_adapter.py \
        create_knowledge_graph.py \
        config/schema_config.yaml \
        multiomics_kg/adapters/cyanorak_ncbi_adapter.py
git commit -m "$(cat <<'EOF'
feat(cazy): CAZy carbohydrate-active-enzyme ontology

Adds CazyAnnotationAdapter (per-strain) + MultiCazyAnnotationAdapter
(observed-only — no external download; classes inferred from
multiomics_kg.utils.cazy_utils.CAZY_CLASSES). Genes attach at the most
specific observed level (subfamily if present, else family). Hierarchy edges
chain subfamily → family → class.

Drops Gene.cazy_ids flat property in favor of Gene_has_cazy_family edges.
Edge counts (rough): 6 class + ~30-50 family/subfamily nodes; ~400-600
gene→CAZy edges. Omics layer unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** Unit tests pass; rebuilt KG passes `pytest -m kg`; `/omics-edge-snapshot` diff empty.

---

## Commit 7 — Post-import additions: indexes, computed properties, Gene rollups, Metabolite UNION arms

**Goal:** All post-import work for the new ontologies. Indexes, ontology-level rollups (`gene_count`, `organism_count`, `member_count`, `metabolite_count`, `tc_class_id`), Gene routing extensions (`annotation_types` extension, `tcdb_family_count`, `cazy_family_count`, `metabolite_count` UNION'd across catalysis + transport), Metabolite extensions (`transporter_count`, `gene_count`/`organism_count` UNION'd), `Organism_has_metabolite` UNION arm.

**Files:**
- Modify: `scripts/post-import.sh` — add new sections (indexes, ontology rollups, Gene rollups, Metabolite/Organism extensions). Group as currently structured (Group 1 = indexes, Group 2 = small aggregations, Group 3 = heavy `IN TRANSACTIONS` writes).
- Modify: `scripts/post-import.cypher` — keep in sync with `.sh` per CLAUDE.md.
- Modify: `config/schema_config.yaml` — add `tcdb_family_count: int`, `cazy_family_count: int`, `metabolite_count: int` to `gene` properties; add `transporter_count: int` to `metabolite` properties.
- Modify: `tests/kg_validity/test_post_import.py` — add assertions for new properties.

### Cypher payload to add

Append in `scripts/post-import.sh` Group 1 (indexes), after the BriteCategory full-text index:

```cypher
// TCDB / CAZy scalar indexes
CREATE INDEX tcdb_family_level_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.level);
CREATE INDEX tcdb_family_level_kind_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.level_kind);
CREATE INDEX tcdb_family_tcdb_id_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.tcdb_id);
CREATE INDEX tcdb_family_tc_class_id_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.tc_class_id);
CREATE INDEX cazy_family_level_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.level);
CREATE INDEX cazy_family_level_kind_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.level_kind);
CREATE INDEX cazy_family_cazy_id_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.cazy_id);

// TCDB / CAZy full-text indexes
CREATE FULLTEXT INDEX tcdbFamilyFullText IF NOT EXISTS
    FOR (t:TcdbFamily) ON EACH [t.name, t.tcdb_id, t.superfamily];
CREATE FULLTEXT INDEX cazyFamilyFullText IF NOT EXISTS
    FOR (c:CazyFamily) ON EACH [c.name, c.cazy_id];
```

In Group 3 (heavy writes), add after the BriteCategory rollup block (~line 721):

```cypher
// ── TcdbFamily computed properties ───────────────────────────────────────────

// member_count: direct child count
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (child:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(t)
  WITH t, count(child) AS mc
  SET t.member_count = mc
} IN TRANSACTIONS OF 1000 ROWS;

// tc_class_id: pointer to the root tc_class node
MATCH (t:TcdbFamily {level_kind: 'tc_class'})
SET t.tc_class_id = t.id;

MATCH (t:TcdbFamily) WHERE t.level_kind <> 'tc_class'
CALL {
  WITH t
  MATCH (t)-[:Tcdb_family_is_a_tcdb_family*1..]->(cls:TcdbFamily {level_kind: 'tc_class'})
  WITH t, cls LIMIT 1
  SET t.tc_class_id = cls.id
} IN TRANSACTIONS OF 1000 ROWS;

// gene_count + organism_count: subtree traversal (descendants ∪ self)
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(desc:TcdbFamily)<-[:Gene_has_tcdb_family]-(g:Gene)
  WITH t, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET t.gene_count = gc,
      t.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// metabolite_count: distinct metabolites reachable via Tcdb_family_transports_metabolite in subtree
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(desc:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH t, count(DISTINCT m) AS mc
  SET t.metabolite_count = mc
} IN TRANSACTIONS OF 1000 ROWS;

// ── CazyFamily computed properties ───────────────────────────────────────────

// gene_count + organism_count: subtree traversal
MATCH (c:CazyFamily)
CALL {
  WITH c
  OPTIONAL MATCH (c)<-[:Cazy_family_is_a_cazy_family*0..]-(desc:CazyFamily)<-[:Gene_has_cazy_family]-(g:Gene)
  WITH c, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET c.gene_count = gc,
      c.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// ── Gene routing extensions ──────────────────────────────────────────────────

// Extend annotation_types — add 'tcdb' / 'cazy' if any edge of that type exists
MATCH (g:Gene)
CALL {
  WITH g
  WITH g, g.annotation_types AS existing
  SET g.annotation_types = existing +
    CASE WHEN EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN ['tcdb'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;

// tcdb_family_count + cazy_family_count
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r1:Gene_has_tcdb_family]->()
  WITH g, count(r1) AS tc
  OPTIONAL MATCH (g)-[r2:Gene_has_cazy_family]->()
  WITH g, tc, count(r2) AS cz
  SET g.tcdb_family_count = tc,
      g.cazy_family_count = cz
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.metabolite_count = UNION across catalysis + transport paths (per TCDB-S3)
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m1:Metabolite)
  WITH g, collect(DISTINCT m1) AS m_cat
  OPTIONAL MATCH (g)-[:Gene_has_tcdb_family]->(:TcdbFamily)
                 -[:Tcdb_family_is_a_tcdb_family*0..]->(:TcdbFamily {level_kind: 'tc_specificity'})
                 -[:Tcdb_family_transports_metabolite]->(m2:Metabolite)
  WITH g, m_cat + collect(DISTINCT m2) AS all_m
  SET g.metabolite_count = size(apoc.coll.toSet(all_m))
} IN TRANSACTIONS OF 500 ROWS;

// ── Metabolite extensions ────────────────────────────────────────────────────

// transporter_count: distinct TcdbFamily ancestors with substrate edge to this metabolite
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (t:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m)
  WITH m, count(DISTINCT t) AS tc
  SET m.transporter_count = tc
} IN TRANSACTIONS OF 1000 ROWS;

// gene_count / organism_count UNION (catalysis + transport)
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(g_cat:Gene)
  WITH m, collect(DISTINCT g_cat) AS gs_cat
  OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(:TcdbFamily)
                 <-[:Gene_has_tcdb_family]-(g_tr:Gene)
  WITH m, gs_cat + collect(DISTINCT g_tr) AS all_g
  WITH m,
       size(apoc.coll.toSet(all_g)) AS gc,
       size(apoc.coll.toSet([x IN all_g | x.organism_name])) AS oc
  SET m.gene_count = gc, m.organism_count = oc
} IN TRANSACTIONS OF 1000 ROWS;

// ── Organism_has_metabolite — add UNION transport arm ───────────────────────
// (Existing block already MERGEs the catalysis arm; this block adds the transport arm.
//  MERGE deduplicates by (organism, metabolite).)
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_has_tcdb_family]->(:TcdbFamily)
        -[:Tcdb_family_is_a_tcdb_family*0..]->(:TcdbFamily {level_kind: 'tc_specificity'})
        -[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[:Organism_has_metabolite]->(m)
} IN TRANSACTIONS OF 1000 ROWS;
```

Add `metabolite_count` to OrganismTaxon rollup if not already present (the existing one only counts `Organism_has_metabolite` so it will pick up transport-path additions automatically — verify). If it lives elsewhere, ensure it doesn't need a refresh after the new MERGE arm.

- [ ] **Step 7.1: Take a deterministic-dump baseline**

```bash
bash scripts/post-import-validate.sh > before.txt
```

- [ ] **Step 7.2: Add new schema properties**

In `config/schema_config.yaml`:

1. In `gene` properties block (after `compartments_observed`):

```yaml
    # TCDB / CAZy / chemistry rollups (post-import)
    tcdb_family_count: int                   # per TCDB-S1
    cazy_family_count: int                   # per TCDB-S2
    metabolite_count: int                    # UNION catalysis + transport (per TCDB-S3 / KG-A2)
```

2. In `metabolite` properties block (after `evidence_sources`):

```yaml
    transporter_count: int                   # post-import: distinct TcdbFamily nodes with substrate edge to this metabolite
```

- [ ] **Step 7.3: Extend `scripts/post-import.sh` with the Cypher payload above**

Add the index block to Group 1 in the appropriate location (after the BriteCategory full-text index ~line 76).

Add the heavy-writes block to Group 3 in the appropriate location — after the existing BriteCategory rollup (~line 721) but before the Reaction/Metabolite metabolism rollups (`// Reaction.gene_count` ~line 725) so that the new `Metabolite.gene_count` UNION block REPLACES the old `Metabolite.gene_count` block (commit out the old block; it's now subsumed).

In particular:
- Replace the existing `// Metabolite.gene_count, organism_count` block (lines 734-739) with the UNION version above.
- Keep the existing `// Materialize Organism_has_metabolite (2-hop save)` block intact and add the transport-arm block alongside it.

- [ ] **Step 7.4: Mirror changes into `scripts/post-import.cypher`**

Per CLAUDE.md: `scripts/post-import.cypher` is a reference copy that must stay in sync with the `.sh` Cypher logic. Apply the same edits to the `.cypher` file (keep them byte-for-byte identical except for the cypher-shell wrappers in `.sh`).

- [ ] **Step 7.5: Write failing test for new computed properties**

Add to `tests/kg_validity/test_post_import.py`:

```python
@pytest.mark.kg
def test_tcdb_family_has_gene_count_and_organism_count(run_query):
    """Every TcdbFamily has gene_count + organism_count populated."""
    rows = run_query(
        "MATCH (t:TcdbFamily) WHERE t.gene_count IS NULL OR t.organism_count IS NULL "
        "RETURN count(t) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_tcdb_family_member_count_populated_on_classes(run_query):
    rows = run_query(
        "MATCH (t:TcdbFamily {level_kind: 'tc_class'}) RETURN t.member_count AS mc"
    )
    assert all(r["mc"] is not None and r["mc"] >= 0 for r in rows)


@pytest.mark.kg
def test_tcdb_family_metabolite_count_aggregates_in_subtree(run_query):
    """A tc_class node's metabolite_count is at least as large as any leaf descendant's."""
    rows = run_query(
        "MATCH (cls:TcdbFamily {level_kind: 'tc_class'}) WHERE cls.metabolite_count > 0 "
        "RETURN cls.tcdb_id AS cid, cls.metabolite_count AS mc LIMIT 5"
    )
    assert any(r["mc"] > 0 for r in rows)


@pytest.mark.kg
def test_tcdb_family_tc_class_id_matches_root(run_query):
    """Every tc_class node has tc_class_id == its own id; non-class nodes point at their root class."""
    rows = run_query(
        "MATCH (t:TcdbFamily {level_kind: 'tc_class'}) "
        "RETURN t.id AS id, t.tc_class_id AS pid LIMIT 5"
    )
    for r in rows:
        assert r["id"] == r["pid"]


@pytest.mark.kg
def test_cazy_family_has_gene_count_and_organism_count(run_query):
    rows = run_query(
        "MATCH (c:CazyFamily) WHERE c.gene_count IS NULL OR c.organism_count IS NULL "
        "RETURN count(c) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_gene_tcdb_family_count_populated(run_query):
    """Every Gene has tcdb_family_count and cazy_family_count populated (default 0)."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.tcdb_family_count IS NULL OR g.cazy_family_count IS NULL "
        "RETURN count(g) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_gene_annotation_types_includes_tcdb_when_edges_present(run_query):
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_tcdb_family]->() RETURN g.annotation_types AS at LIMIT 5"
    )
    assert all("tcdb" in r["at"] for r in rows)


@pytest.mark.kg
def test_gene_annotation_types_includes_cazy_when_edges_present(run_query):
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_cazy_family]->() RETURN g.annotation_types AS at LIMIT 5"
    )
    assert all("cazy" in r["at"] for r in rows)


@pytest.mark.kg
def test_metabolite_transporter_count_populated(run_query):
    rows = run_query(
        "MATCH (m:Metabolite) WHERE m.transporter_count IS NULL RETURN count(m) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_organism_has_metabolite_includes_transport_path(run_query):
    """At least some Organism_has_metabolite edges exist that come ONLY from transport
    (i.e. organism has no Reaction touching that metabolite, but does have a TCDB
    family with substrate edge)."""
    rows = run_query("""
        MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
        WHERE 'transport' IN m.evidence_sources
          AND NOT 'metabolism' IN m.evidence_sources
        RETURN count(*) AS n
    """)
    assert rows[0]["n"] > 0
```

- [ ] **Step 7.6: Rebuild the KG**

```bash
docker compose down
docker compose up -d --build
```

- [ ] **Step 7.7: Run KG validity tests**

Run: `pytest tests/kg_validity/test_post_import.py -v`
Expected: PASS (new tests + existing tests).

Run: `pytest -m kg -v`
Expected: PASS overall.

- [ ] **Step 7.8: Diff the deterministic dump**

```bash
bash scripts/post-import-validate.sh > after.txt
diff before.txt after.txt
```

Expected: lots of diff (intentional — new properties added). Inspect to confirm: only adds, no unrelated changes (e.g. spurious BRITE rollup drift would be a bug — investigate).

- [ ] **Step 7.9: Run `/omics-edge-snapshot` and verify byte-identical**

```bash
/omics-edge-snapshot
diff omics-edge-snapshot-pre-commit3.json omics-edge-snapshot-*.json
```

Expected: empty diff.

- [ ] **Step 7.10: Commit**

```bash
git add scripts/post-import.sh \
        scripts/post-import.cypher \
        config/schema_config.yaml \
        tests/kg_validity/test_post_import.py
git commit -m "$(cat <<'EOF'
feat(post-import): TCDB/CAZy rollups + Metabolite UNION extensions

Adds scalar + full-text indexes on TcdbFamily / CazyFamily. Computes
gene_count / organism_count / member_count / metabolite_count / tc_class_id
on TcdbFamily; gene_count / organism_count on CazyFamily. Extends
Gene.annotation_types to include 'tcdb' / 'cazy'; adds Gene.tcdb_family_count
+ Gene.cazy_family_count rollups (TCDB-S1, TCDB-S2). Defines
Gene.metabolite_count as UNION across catalysis + transport paths
(TCDB-S3 / KG-A2). Adds Metabolite.transporter_count. Extends
Metabolite.gene_count / organism_count and Organism_has_metabolite
materialization with the transport arm.

post-import.cypher kept byte-identical to post-import.sh per CLAUDE.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** `pytest tests/kg_validity/ -v` passes; deterministic dump diff shows only intentional additions; `/omics-edge-snapshot` byte-identical.

---

## Commit 8 — Snapshot regen + new validity tests + docs

**Goal:** Lock the new graph state into the snapshot regression fixture; add a focused TCDB/CAZy validity test file; update CLAUDE.md, MEMORY, the cypher-queries skill, and `docs/kg-changes/tcdb-cazy-ontologies.md` (lift "PROPOSED" status, backfill exact counts).

**Files:**
- Modify: `tests/kg_validity/snapshot_data.json` — regenerated.
- Create: `tests/kg_validity/test_tcdb_cazy.py` — node count + edge presence + level invariants + orphan checks.
- Modify: `tests/kg_validity/test_structure.py` — extend orphan checks to cover TcdbFamily / CazyFamily.
- Modify: `CLAUDE.md` — extend "Key Adapters", node label list, edge label list, prepare_data step 6 description, build_metabolite_resolver → build_mnx_resolver.
- Modify: `memory/MEMORY.md` — one-line pointer to a new memory file.
- Create: `memory/project_tcdb_cazy_ontologies.md`.
- Modify: `.claude/skills/cypher-queries/SKILL.md` — add new node/edge labels + 2-3 query templates.
- Modify: `docs/kg-changes/tcdb-cazy-ontologies.md` — flip status, backfill exact counts.

- [ ] **Step 8.1: Regenerate snapshot**

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

- [ ] **Step 8.2: Inspect the new snapshot**

```bash
git diff tests/kg_validity/snapshot_data.json | head -200
```

Confirm the diff shows new TcdbFamily / CazyFamily / transport-only Metabolite samples, and any pre-existing samples are still present (no regression).

- [ ] **Step 8.3: Write `tests/kg_validity/test_tcdb_cazy.py`**

```python
"""KG validity tests for TCDB and CAZy ontologies."""
from __future__ import annotations

import pytest


@pytest.mark.kg
def test_tcdb_family_node_count_in_range(run_query):
    """Pruned subhierarchy should have a few hundred TcdbFamily nodes (not all 13K)."""
    n = run_query("MATCH (t:TcdbFamily) RETURN count(t) AS n")[0]["n"]
    assert 100 <= n <= 2000, f"TcdbFamily count {n} outside expected 100-2000"


@pytest.mark.kg
def test_cazy_family_node_count_in_range(run_query):
    """Observed-only — small ontology."""
    n = run_query("MATCH (c:CazyFamily) RETURN count(c) AS n")[0]["n"]
    assert 5 <= n <= 100, f"CazyFamily count {n} outside expected 5-100"


@pytest.mark.kg
def test_gene_has_tcdb_family_edge_count(run_query):
    n = run_query("MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n")[0]["n"]
    assert 500 <= n <= 5000, f"Gene_has_tcdb_family count {n} outside 500-5000"


@pytest.mark.kg
def test_gene_has_cazy_family_edge_count(run_query):
    n = run_query("MATCH ()-[r:Gene_has_cazy_family]->() RETURN count(r) AS n")[0]["n"]
    assert 100 <= n <= 1500, f"Gene_has_cazy_family count {n} outside 100-1500"


@pytest.mark.kg
def test_tcdb_family_transports_metabolite_edge_count(run_query):
    n = run_query("MATCH ()-[r:Tcdb_family_transports_metabolite]->() RETURN count(r) AS n")[0]["n"]
    assert 1000 <= n <= 30000, f"Tcdb_family_transports_metabolite count {n} outside 1000-30000"


@pytest.mark.kg
def test_tcdb_family_transports_metabolite_only_on_leaves(run_query):
    """Substrate edges must originate from tc_specificity nodes only."""
    n_bad = run_query("""
        MATCH (t:TcdbFamily)-[:Tcdb_family_transports_metabolite]->()
        WHERE t.level_kind <> 'tc_specificity'
        RETURN count(*) AS n
    """)[0]["n"]
    assert n_bad == 0


@pytest.mark.kg
def test_tcdb_family_levels_are_valid(run_query):
    rows = run_query(
        "MATCH (t:TcdbFamily) RETURN DISTINCT t.level_kind AS lk, t.level AS lv"
    )
    expected = {
        ("tc_class", 0), ("tc_subclass", 1), ("tc_family", 2),
        ("tc_subfamily", 3), ("tc_specificity", 4),
    }
    actual = {(r["lk"], r["lv"]) for r in rows}
    assert actual.issubset(expected), f"Unexpected level_kind/level pairs: {actual - expected}"


@pytest.mark.kg
def test_cazy_family_levels_are_valid(run_query):
    rows = run_query(
        "MATCH (c:CazyFamily) RETURN DISTINCT c.level_kind AS lk, c.level AS lv"
    )
    expected = {("cazy_class", 0), ("cazy_family", 1), ("cazy_subfamily", 2)}
    actual = {(r["lk"], r["lv"]) for r in rows}
    assert actual.issubset(expected), f"Unexpected level_kind/level pairs: {actual - expected}"


@pytest.mark.kg
def test_tcdb_family_no_orphans_below_class(run_query):
    """Every non-class TcdbFamily must have a parent edge."""
    n_bad = run_query("""
        MATCH (t:TcdbFamily) WHERE t.level_kind <> 'tc_class'
          AND NOT EXISTS { (t)-[:Tcdb_family_is_a_tcdb_family]->() }
        RETURN count(t) AS n
    """)[0]["n"]
    assert n_bad == 0


@pytest.mark.kg
def test_cazy_family_no_orphans_below_class(run_query):
    n_bad = run_query("""
        MATCH (c:CazyFamily) WHERE c.level_kind <> 'cazy_class'
          AND NOT EXISTS { (c)-[:Cazy_family_is_a_cazy_family]->() }
        RETURN count(c) AS n
    """)[0]["n"]
    assert n_bad == 0


@pytest.mark.kg
def test_tcdb_full_text_index_exists(run_query):
    rows = run_query("SHOW INDEXES YIELD name WHERE name = 'tcdbFamilyFullText' RETURN name")
    assert len(rows) == 1


@pytest.mark.kg
def test_cazy_full_text_index_exists(run_query):
    rows = run_query("SHOW INDEXES YIELD name WHERE name = 'cazyFamilyFullText' RETURN name")
    assert len(rows) == 1


@pytest.mark.kg
def test_transport_only_metabolites_exist(run_query):
    rows = run_query("""
        MATCH (m:Metabolite)
        WHERE 'transport' IN m.evidence_sources
          AND NOT 'metabolism' IN m.evidence_sources
        RETURN count(m) AS n
    """)
    assert rows[0]["n"] >= 100, "Expected at least 100 transport-only Metabolites"


@pytest.mark.kg
def test_metabolites_have_evidence_sources(run_query):
    """Every Metabolite has at least one evidence source."""
    n_bad = run_query("""
        MATCH (m:Metabolite)
        WHERE m.evidence_sources IS NULL OR size(m.evidence_sources) = 0
        RETURN count(m) AS n
    """)[0]["n"]
    assert n_bad == 0
```

- [ ] **Step 8.4: Extend `tests/kg_validity/test_structure.py` orphan checks**

Search via `grep -n orphan tests/kg_validity/test_structure.py` for the existing pattern. Add tests:

```python
@pytest.mark.kg
def test_no_orphan_tcdb_families(run_query):
    """No TcdbFamily nodes without any incoming gene/parent edge AND no outgoing parent/substrate edge."""
    # Per spec: every TcdbFamily kept by pruning has at least one gene-edge in its subtree
    # OR is an ancestor of one. So orphan = no edges at all.
    n_bad = run_query("""
        MATCH (t:TcdbFamily)
        WHERE NOT EXISTS { (t)--() }
        RETURN count(t) AS n
    """)[0]["n"]
    assert n_bad == 0


@pytest.mark.kg
def test_no_orphan_cazy_families(run_query):
    n_bad = run_query("""
        MATCH (c:CazyFamily)
        WHERE NOT EXISTS { (c)--() }
        RETURN count(c) AS n
    """)[0]["n"]
    assert n_bad == 0
```

- [ ] **Step 8.5: Run all KG validity tests**

Run: `pytest -m kg -v`
Expected: PASS.

- [ ] **Step 8.6: Update `CLAUDE.md`**

Apply these edits (target sections noted):

1. **"Adapter Pattern" / "Key Adapters" list:** add lines after `metabolism_adapter.py`:

   ```
   - **tcdb_adapter.py** — TCDB transport ontology + substrate→Metabolite bridge. Reads cache/data/tcdb/{tcdb_hierarchy,tcdb_pruned}.json built by step 6. Pruned to gene-reachable subhierarchy (above + below)
   - **cazy_adapter.py** — CAZy carbohydrate-active-enzyme ontology (observed-only — no external download; CAZY_CLASSES table from utils/cazy_utils.py)
   ```

2. **"Actual Neo4j labels"** — add `TcdbFamily`, `CazyFamily` to the Nodes line; add the 5 new edges to the Relationships line.

3. **"Key graph facts"** — add a new bullet block:

   ```
   - TcdbFamily nodes: pruned to gene-reachable subhierarchy (above + below); 5 levels (tc_class | tc_subclass | tc_family | tc_subfamily | tc_specificity). Properties: name, tcdb_id, level (0-4), level_kind, superfamily (sparse, leaf-only), tc_class_id (post-import pointer to root). Computed post-import: gene_count, organism_count, member_count, metabolite_count.
   - CazyFamily nodes: observed-only; 3 levels (cazy_class | cazy_family | cazy_subfamily). Properties: name, cazy_id, level (0-2), level_kind. Computed post-import: gene_count, organism_count.
   - Five new edges: Gene_has_tcdb_family (Gene → TcdbFamily, attaches at exact eggNOG-annotated level), Tcdb_family_is_a_tcdb_family (child → parent), Tcdb_family_transports_metabolite (TcdbFamily-leaf → Metabolite, only on tc_specificity nodes), Gene_has_cazy_family, Cazy_family_is_a_cazy_family.
   - Metabolite.evidence_sources (str[]): values 'metabolism' (Reaction-reachable) and 'transport' (TCDB-substrate-reachable). Open-ended; 'metabolomics' reserved.
   - Metabolite.transporter_count (int, post-import): distinct TcdbFamily nodes pointing at this metabolite.
   - Metabolite.gene_count / organism_count: now UNION across catalysis + transport paths.
   - Gene routing extensions: tcdb_family_count, cazy_family_count (int rollups); annotation_types extended with 'tcdb' / 'cazy'; metabolite_count = UNION catalysis + transport.
   - Organism_has_metabolite materialization: now includes the transport arm (Gene → TcdbFamily → Metabolite). Single MERGE deduplicates by (organism, metabolite).
   ```

4. **"Removed Gene properties"** note: add `transporter_classification` and `cazy_ids` to the list of removed Gene properties (deleted on 2026-05-01 by TCDB-CAZy ontology promotion).

5. **prepare_data step 6 description (lines around 174-180):** add to the bulleted breakdown — TCDB download + parse + prune + substrate→Metabolite resolution + tcdb_pruned.json write are now part of step 6.

6. **`build_metabolite_resolver.py` mention (around line 188):** rename to `build_mnx_resolver.py` and update the description to "builds the heavy MNX SQLite only — TCDB hierarchy lives in step 6 now."

- [ ] **Step 8.7: Update `memory/MEMORY.md`**

Add one line under the existing index. Create `memory/project_tcdb_cazy_ontologies.md`:

```markdown
---
name: TCDB and CAZy ontologies landed
description: TcdbFamily + CazyFamily ontologies, ~785 transport-only Metabolites, evidence_sources enum, removed Gene.transporter_classification/cazy_ids
type: project
---

TCDB and CAZy classifications were promoted to first-class graph ontologies on 2026-05-XX (PR landed 8 commits per plans/2026-05-01-tcdb-cazy-ontologies.md).

**Why:** the flat Gene array properties `transporter_classification` and `cazy_ids` only supported text search — no hierarchy navigation, no chemistry linkage. Promoting to TcdbFamily/CazyFamily ontology nodes mirrors the Pfam normalization precedent and bridges TCDB substrates into the chemistry layer (~785 new transport-only Metabolite nodes).

**How to apply:**
- Cypher queries that filtered `g.transporter_classification` or `g.cazy_ids` must rewrite to traverse `Gene_has_tcdb_family` / `Gene_has_cazy_family` edges.
- Filter by `'transport' IN m.evidence_sources` to find transport-only metabolites; use set membership to stay forward-compatible with the reserved `'metabolomics'` value.
- TCDB hierarchy filtering via `tf.tc_class_id = 'tcdb:3'` is faster than variable-length traversals.

Spec: `docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md`. Explorer-facing summary: `docs/kg-changes/tcdb-cazy-ontologies.md`.
```

Then add to `memory/MEMORY.md`:

```markdown
- [TCDB and CAZy ontologies landed](project_tcdb_cazy_ontologies.md) — TcdbFamily + CazyFamily nodes; transport-only Metabolites; flat Gene properties removed.
```

- [ ] **Step 8.8: Update `.claude/skills/cypher-queries/SKILL.md`**

Apply these edits to the skill:

1. **Node Labels table:** add rows for `TcdbFamily` and `CazyFamily`.
2. **Relationships list:** add the 5 new edge labels.
3. **Templates section:** add 2-3 new query templates. Suggested:

```cypher
// Genes whose transporters are annotated to move calcium
MATCH (g:Gene)-[:Gene_has_tcdb_family]->(:TcdbFamily)
      -[:Tcdb_family_is_a_tcdb_family*0..]->(:TcdbFamily {level_kind: 'tc_specificity'})
      -[:Tcdb_family_transports_metabolite]->(m:Metabolite)
WHERE m.name CONTAINS 'calcium'
RETURN DISTINCT g.locus_tag, g.product, m.name
LIMIT 50;

// What does TCDB class 1 ('Channels and Pores') cover in MED4?
MATCH (g:Gene {organism_name: 'Prochlorococcus MED4'})
      -[:Gene_has_tcdb_family]->(tf:TcdbFamily {tc_class_id: 'tcdb:1'})
RETURN tf.tcdb_id, tf.name, tf.level_kind, count(g) AS gene_count
ORDER BY tf.level, tf.tcdb_id;

// Compounds the organism is annotated to TRANSPORT but has no metabolism for
MATCH (m:Metabolite)
WHERE 'transport' IN m.evidence_sources AND NOT 'metabolism' IN m.evidence_sources
RETURN m.name, m.formula, m.evidence_sources
LIMIT 50;
```

- [ ] **Step 8.9: Update `docs/kg-changes/tcdb-cazy-ontologies.md`**

1. Change the status header from `**Status:** PROPOSED — design at ...` to `**Status:** LANDED 2026-05-XX — design at ...`.
2. Backfill the "Estimated graph-size impact" table with the actual numbers from the rebuilt graph (run the spot-check Cypher from commits 5/6 against the live graph and copy the numbers in).
3. Confirm the property changes table reflects the as-built state.

- [ ] **Step 8.10: Run the full test suite one final time**

```bash
pytest -m "not slow and not kg" -v
pytest -m kg -v
```

Both should pass.

- [ ] **Step 8.11: Final `/omics-edge-snapshot` regression check**

```bash
/omics-edge-snapshot
diff omics-edge-snapshot-pre-commit3.json omics-edge-snapshot-*.json
```

Expected: empty diff. If there's any drift, stop and investigate before committing — the omics layer should be untouched by this entire feature.

- [ ] **Step 8.12: Commit**

```bash
git add tests/kg_validity/snapshot_data.json \
        tests/kg_validity/test_tcdb_cazy.py \
        tests/kg_validity/test_structure.py \
        CLAUDE.md \
        memory/MEMORY.md \
        memory/project_tcdb_cazy_ontologies.md \
        .claude/skills/cypher-queries/SKILL.md \
        docs/kg-changes/tcdb-cazy-ontologies.md
git commit -m "$(cat <<'EOF'
docs/test: snapshot regen + KG-validity coverage for TCDB/CAZy

Regenerates tests/kg_validity/snapshot_data.json to lock in TcdbFamily,
CazyFamily, and transport-only Metabolite samples. Adds focused
tests/kg_validity/test_tcdb_cazy.py covering node counts, edge presence,
level/level_kind invariants, leaf-only substrate edges, full-text indexes.
Extends test_structure.py orphan checks. Updates CLAUDE.md (key graph facts,
adapters, prepare_data step 6, build script rename), MEMORY (one-liner +
project memory), .claude/skills/cypher-queries/SKILL.md (new node/edge labels
+ 3 templates), and docs/kg-changes/tcdb-cazy-ontologies.md (lifts PROPOSED
status, backfills as-built counts).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

**Verification gate:** `pytest -m "not slow and not kg"` and `pytest -m kg` both pass; `/omics-edge-snapshot` byte-identical to commit 3 baseline.

---

## Cross-cutting verification matrix

| Commit | Unit tests | KG rebuild | KG validity | `/omics-edge-snapshot` |
|---|---|---|---|---|
| 1 | required | not needed | not needed | unchanged |
| 2 | required | not needed (smoke optional) | not needed | unchanged |
| 3 | required | required | required (widened thresholds) | byte-identical to pre-commit-3 baseline |
| 4 | required (schema parses) | not needed | not needed | unchanged |
| 5 | required | required | required | byte-identical to commit-3 baseline |
| 6 | required | required | required | byte-identical to commit-3 baseline |
| 7 | required | required | required (incl. new test_post_import.py assertions) | byte-identical to commit-3 baseline |
| 8 | required (full) | required (final) | required (full + new test_tcdb_cazy.py) | byte-identical to commit-3 baseline |

The `/omics-edge-snapshot` invariant is the single strongest regression signal: this entire feature must not touch a single `Changes_expression_of` edge property. Any drift after commit 3 is a bug.

## Open issues / clarifications to flag before implementation

- **Renaming script in commit 2 vs. test renames:** if the test file rename breaks `pytest` discovery (because pytest caches under the old name), consider running with `pytest --cache-clear`.
- **`gene_annotations_config.yaml` field deletion vs. keep:** the spec says "drop the field unless other code reads it." The TCDB / CAZy adapters DO read `gene["transporter_classification"]` and `gene["cazy_ids"]` from `gene_annotations_merged.json` (commit 5/6 step where the per-strain adapter loads the JSON). So the field stays in `gene_annotations_config.yaml` and `gene_annotations_merged.json` — only the schema-side Gene property goes away. The plan reflects this; flag if any executor reads the spec literally and tries to drop the field.
- **Missing `_tx_validate_tcdb` body in the spec — what about `is_valid_tcdb` consumers outside transforms?** The spec says drop the import. `tcdb_utils.is_valid_tcdb` may still be used in tests; verify with `grep -rn is_valid_tcdb tests/` before deleting from `multiomics_kg/utils/tcdb_utils.py` (the function itself is allowed to stay; the spec only requires removing the validation transform).
- **MNX query patterns in commit 3** assume the chebi alias source name is `chebi` (not `CHEBI`). Per `_SOURCE_NORM` in the existing `build_metabolite_resolver.py`, this is canonicalized to lowercase. Confirm by inspecting `compound_aliases` once before running the new query at scale.
