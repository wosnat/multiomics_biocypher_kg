# Unified Ontology `level` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a canonical `level: int` property to every ontology-term node (ten labels) and a sparse `level_is_best_effort: "true"` marker on GO DAG terms whose min-path ≠ max-path. Rename `KeggTerm.level` (string) to `KeggTerm.level_kind`.

**Architecture:** Pure-adapter implementation. No post-import Cypher changes. Each adapter already has the structural information needed — GO via in-memory BFS over `go_utils.load_go_data()`, EC via nested loops, CyanorakRole via parent pointers, others by kind. GO BFS uses the three canonical roots (`GO:0008150`, `GO:0003674`, `GO:0005575`) and tracks both min and max depth.

**Tech Stack:** Python (adapter generators), BioCypher schema YAML, pytest (unit + KG-validity), Neo4j (integration verification only).

**Spec:** [docs/superpowers/specs/2026-04-12-ontology-level-design.md](../specs/2026-04-12-ontology-level-design.md)

---

## File Structure

**New files:**
- `tests/test_go_levels.py` — unit tests for the GO BFS helper.
- `tests/kg_validity/test_ontology_level.py` — live-graph verification.

**Modified files:**
- `config/schema_config.yaml` — add `level` / `level_is_best_effort` to 10 node types; rename KEGG `level` → `level_kind`.
- `multiomics_kg/utils/go_utils.py` — add `compute_go_levels` helper + canonical root constants.
- `multiomics_kg/adapters/go_adapter.py` — call the helper and attach properties inside `get_go_nodes()`.
- `multiomics_kg/adapters/ec_adapter.py` — set `level` in each nested loop in `get_nodes()`.
- `multiomics_kg/adapters/functional_annotation_adapter.py` — KEGG rename + levels, CyanorakRole depth helper, TigrRole/COG/Pfam/PfamClan constants.
- `tests/test_kegg_annotation_adapter.py` — update KEGG assertions for `level_kind` + new int `level`.
- `CLAUDE.md` — one bullet under "Actual Neo4j labels" documenting the unified convention.

---

## Task 1: Schema additions in `schema_config.yaml`

**Files:**
- Modify: `config/schema_config.yaml:352-376` (GO labels)
- Modify: `config/schema_config.yaml:405-413` (KEGG)
- Modify: `config/schema_config.yaml:415-469` (EcNumber, COG, Cyanorak, Tigr, Pfam, PfamClan)

- [ ] **Step 1: Add `level` + `level_is_best_effort` to GO labels**

In `config/schema_config.yaml`, under each of `biological process`, `cellular component`, `molecular function`, add two lines to the `properties:` block:

```yaml
    level: int
    level_is_best_effort: str
```

Example (biological process, lines 352–359 become):

```yaml
biological process:
  is_a: gene ontology
  represented_as: node
  preferred_id: go
  label_in_input: biological process
  properties:
    name: str
    anc2vec_embedding: float[]
    level: int
    level_is_best_effort: str
```

Do the same for `cellular component` and `molecular function`.

- [ ] **Step 2: Rename KEGG `level` → `level_kind` and add int `level`**

In `config/schema_config.yaml:405-412`, replace:

```yaml
kegg term:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.orthology
  label_in_input: kegg_term
  properties:
    name: str
    level: str
```

with:

```yaml
kegg term:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.orthology
  label_in_input: kegg_term
  properties:
    name: str
    level_kind: str
    level: int
    level_is_best_effort: str
```

- [ ] **Step 3: Add `level` + `level_is_best_effort` to the remaining six labels**

For each of `ec number`, `cog functional category`, `cyanorak role`, `tigr role`, `pfam`, `pfam clan` (lines 415–469), append to the `properties:` block:

```yaml
    level: int
    level_is_best_effort: str
```

- [ ] **Step 4: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add unified level/level_is_best_effort properties to 10 ontology node types; rename KeggTerm.level → level_kind"
```

---

## Task 2: `compute_go_levels` helper in `go_utils.py` (unit-test-driven)

**Files:**
- Modify: `multiomics_kg/utils/go_utils.py` (add function at module level)
- Test: `tests/test_go_levels.py` (new file)

- [ ] **Step 1: Write failing unit tests**

Create `tests/test_go_levels.py`:

```python
"""
Unit tests for compute_go_levels (GO BFS level assignment).

Uses synthetic go_data dicts — does not hit the real OBO file or KG.
"""
from multiomics_kg.utils.go_utils import (
    CANONICAL_GO_ROOTS,
    compute_go_levels,
)


def _bp(parents):
    """Helper: build a biological_process term dict entry."""
    return {"name": "", "namespace": "biological_process", "parents": parents}


def test_root_only_gets_level_zero():
    go_data = {"GO:0008150": _bp([])}
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert orphans == []


def test_linear_chain_assigns_monotonic_levels():
    go_data = {
        "GO:0008150": _bp([]),                       # root
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),  # depth 1
        "GO:0000002": _bp([["GO:0000001", "is_a"]]),  # depth 2
    }
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert levels["GO:0000001"] == (1, False)
    assert levels["GO:0000002"] == (2, False)
    assert orphans == []


def test_diamond_equal_arms_not_best_effort():
    # Two parents both at depth 1 → confluence node at depth 2 with min==max.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),
        "GO:0000002": _bp([["GO:0008150", "is_a"]]),
        "GO:0000003": _bp([["GO:0000001", "is_a"], ["GO:0000002", "is_a"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000003"] == (2, False)


def test_diamond_unequal_arms_flags_best_effort():
    # One parent at depth 1, another at depth 2 → min=2, max=3, flagged.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),        # depth 1
        "GO:0000002": _bp([["GO:0000001", "is_a"]]),         # depth 2
        "GO:0000003": _bp([["GO:0008150", "is_a"], ["GO:0000002", "is_a"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000003"] == (1, True)  # min=1 via direct root, max=3 via chain


def test_part_of_is_traversed():
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "part_of"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000001"] == (1, False)


def test_regulates_is_not_traversed():
    # A term whose only parent is via 'regulates' must not be reachable
    # from a root → orphan.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "regulates"]]),
    }
    levels, orphans = compute_go_levels(go_data)
    assert "GO:0000001" not in levels
    assert "GO:0000001" in orphans


def test_orphan_not_canonical_root_is_surfaced():
    # Term with empty parents that isn't one of the three canonical roots.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0999999": _bp([]),
    }
    levels, orphans = compute_go_levels(go_data)
    assert "GO:0999999" not in levels
    assert orphans == ["GO:0999999"]


def test_all_three_namespaces_handled():
    go_data = {
        "GO:0008150": {"name": "", "namespace": "biological_process", "parents": []},
        "GO:0003674": {"name": "", "namespace": "molecular_function", "parents": []},
        "GO:0005575": {"name": "", "namespace": "cellular_component", "parents": []},
    }
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert levels["GO:0003674"] == (0, False)
    assert levels["GO:0005575"] == (0, False)
    assert orphans == []


def test_canonical_roots_constant():
    assert CANONICAL_GO_ROOTS == frozenset({
        "GO:0008150",  # biological_process
        "GO:0003674",  # molecular_function
        "GO:0005575",  # cellular_component
    })
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/test_go_levels.py -v
```

Expected: `ImportError: cannot import name 'CANONICAL_GO_ROOTS' from 'multiomics_kg.utils.go_utils'` (or similar).

- [ ] **Step 3: Add `CANONICAL_GO_ROOTS` constant and `compute_go_levels` helper**

In `multiomics_kg/utils/go_utils.py`, after the `NAMESPACE_TO_LABEL` constant (around line 28), add:

```python
# Canonical GO namespace root IDs. Any term with empty parents that is not
# one of these is treated as an orphan (data anomaly).
CANONICAL_GO_ROOTS: frozenset[str] = frozenset({
    "GO:0008150",  # biological_process
    "GO:0003674",  # molecular_function
    "GO:0005575",  # cellular_component
})

# Parent relations traversed for hierarchy-depth computation.
# `regulates` / `positively_regulates` / `negatively_regulates` are intentionally
# excluded — they don't carry hierarchy semantics for the explorer's level rollup.
_LEVEL_TRAVERSAL_RELATIONS: frozenset[str] = frozenset({"is_a", "part_of"})
```

Then append at the end of the module:

```python
def compute_go_levels(
    go_data: dict[str, dict],
) -> tuple[dict[str, tuple[int, bool]], list[str]]:
    """
    Compute hierarchy depth for every GO term reachable from its namespace root.

    Traverses ``is_a`` and ``part_of`` parents only. For each term, tracks the
    minimum and maximum depth observed along any root path; returns
    ``(min_depth, is_best_effort)`` where ``is_best_effort = min_depth != max_depth``.

    Args:
        go_data: ``{go_id: {name, namespace, parents: [[parent_id, relation], ...]}}``
            as produced by :func:`load_go_data`.

    Returns:
        A pair ``(levels, orphans)`` where:

        - ``levels`` maps every reachable go_id to ``(min_depth, is_best_effort)``.
        - ``orphans`` lists go_ids that are not one of the canonical roots and
          whose only path to a root is via a non-traversed relation (or none
          at all).
    """
    # Reverse adjacency: for each term, which terms list it as a parent
    # (via an allowed relation)?
    children: dict[str, list[str]] = {gid: [] for gid in go_data}
    for child_id, entry in go_data.items():
        for parent_id, relation in entry.get("parents", []):
            if relation not in _LEVEL_TRAVERSAL_RELATIONS:
                continue
            if parent_id in children:
                children[parent_id].append(child_id)

    min_depth: dict[str, int] = {}
    max_depth: dict[str, int] = {}

    # Seed: every canonical root present in go_data starts at depth 0.
    seeds = [gid for gid in CANONICAL_GO_ROOTS if gid in go_data]
    for root in seeds:
        min_depth[root] = 0
        max_depth[root] = 0

    # Iterative relaxation: each child's depth is at most/at least
    # parent_depth + 1. Loop until no further updates.
    # GO is a DAG with ~6k terms — converges in a handful of passes.
    changed = True
    while changed:
        changed = False
        # Snapshot IDs already reached to iterate stably.
        frontier = list(min_depth.keys())
        for parent_id in frontier:
            pd_min = min_depth[parent_id]
            pd_max = max_depth[parent_id]
            for child_id in children.get(parent_id, []):
                new_min = pd_min + 1
                new_max = pd_max + 1
                if child_id not in min_depth or new_min < min_depth[child_id]:
                    min_depth[child_id] = new_min
                    changed = True
                if child_id not in max_depth or new_max > max_depth[child_id]:
                    max_depth[child_id] = new_max
                    changed = True

    levels: dict[str, tuple[int, bool]] = {
        gid: (min_depth[gid], min_depth[gid] != max_depth[gid])
        for gid in min_depth
    }
    orphans: list[str] = sorted(
        gid for gid in go_data
        if gid not in levels and gid not in CANONICAL_GO_ROOTS
    )
    return levels, orphans
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
uv run pytest tests/test_go_levels.py -v
```

Expected: all 9 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/go_utils.py tests/test_go_levels.py
git commit -m "feat(go_utils): add compute_go_levels BFS helper + CANONICAL_GO_ROOTS"
```

---

## Task 3: GO adapter integration

**Files:**
- Modify: `multiomics_kg/adapters/go_adapter.py` (the `GO.get_go_nodes` method around lines 807–855)

- [ ] **Step 1: Import the helper**

Near the top of `go_adapter.py`, below the existing imports, add:

```python
from pathlib import Path

from multiomics_kg.utils.go_utils import compute_go_levels, load_go_data
```

(If imports from `multiomics_kg.utils.go_utils` already exist in the file, merge into the existing import statement.)

- [ ] **Step 2: Compute levels at node emission time**

In `go_adapter.py`, modify `get_go_nodes` (around line 807). After the existing `counter = 0` line (around line 823) but before the `for go_term in tqdm(...)` loop, insert:

```python
        # Compute hierarchy levels once, up-front. We load the compact cache
        # separately from pypath's GeneOntology — it's cheap (~5 MB JSON) and
        # keeps the level computation self-contained.
        go_data = load_go_data(Path("cache/data"))
        go_levels, go_orphans = compute_go_levels(go_data)
        if go_orphans:
            logger.warning(
                f"GO.get_go_nodes: {len(go_orphans)} orphan terms "
                f"(not reachable from canonical roots via is_a/part_of): "
                f"{go_orphans[:10]}{'...' if len(go_orphans) > 10 else ''}"
            )
```

Then inside the loop body, after the existing `node_props` assembly but before `node_list.append(...)` (around line 848), add:

```python
                # Unified hierarchy depth (spec: docs/superpowers/specs/2026-04-12-ontology-level-design.md)
                if go_term in go_levels:
                    depth, is_best_effort = go_levels[go_term]
                    node_props["level"] = depth
                    if is_best_effort:
                        node_props["level_is_best_effort"] = "true"
```

- [ ] **Step 3: Commit**

Since this is an adapter wiring change without a dedicated unit test (the BFS itself is unit-tested; integration is covered by the KG-validity tests in Task 8), commit now:

```bash
git add multiomics_kg/adapters/go_adapter.py
git commit -m "feat(go_adapter): emit unified level + sparse level_is_best_effort on GO nodes"
```

---

## Task 4: EcNumber adapter

**Files:**
- Modify: `multiomics_kg/adapters/ec_adapter.py:235-286`

- [ ] **Step 1: Set `level` inside the four nested loops**

In `ec_adapter.py`, at each of the four `node_list.append(...)` calls inside `get_nodes()`, add `"level"` to the `props` dict.

Specifically:

Line ~246 — after `props[ECNodeField.NAME.value] = level_1_dict["name"]` and before `node_list.append((level_1_id, label, props))`, insert:

```python
            props["level"] = 0
```

Line ~255 — after `props[ECNodeField.NAME.value] = level_2_dict["name"]` and before `node_list.append((level_2_id, label, props))`:

```python
                    props["level"] = 1
```

Line ~266 — after `props[ECNodeField.NAME.value] = level_3_dict["name"]`:

```python
                            props["level"] = 2
```

Line ~279 — inside the dict comprehension producing `props` for level 4, change the assignment to:

```python
                                    props = {
                                        field_name : self.clean_text(self.enzymes[level_4_entry][dict_key])
                                        for field_name, dict_key in field_name_to_dict_key.items()
                                        if (field_name in self.ec_node_fields) and (dict_key in self.enzymes[level_4_entry])
                                    }
                                    props["level"] = 3
```

(i.e. append one line after the dict-comp assignment.)

- [ ] **Step 2: Commit**

```bash
git add multiomics_kg/adapters/ec_adapter.py
git commit -m "feat(ec_adapter): emit level (0–3) on EcNumber nodes via nested-loop depth"
```

---

## Task 5: KEGG — rename `level` → `level_kind`, add int `level`

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:601-630`
- Modify: `tests/test_kegg_annotation_adapter.py:329-407`

- [ ] **Step 1: Update KEGG node emitters**

In `functional_annotation_adapter.py`, update the four KEGG `yield` statements:

Line 601 — change:

```python
            yield (_ko_node_id(ko_id), "kegg_term", {"name": name, "level": "ko"})
```

to:

```python
            yield (_ko_node_id(ko_id), "kegg_term", {"name": name, "level_kind": "ko", "level": 3})
```

Line 614 — change:

```python
            yield (_pathway_node_id(pw_id), "kegg_term", {"name": name, "level": "pathway"})
```

to:

```python
            yield (_pathway_node_id(pw_id), "kegg_term", {"name": name, "level_kind": "pathway", "level": 2})
```

Line 622 — change:

```python
            yield (_subcat_node_id(sc), "kegg_term", {"name": name, "level": "subcategory"})
```

to:

```python
            yield (_subcat_node_id(sc), "kegg_term", {"name": name, "level_kind": "subcategory", "level": 1})
```

Line 630 — change:

```python
            yield (_cat_node_id(cat), "kegg_term", {"name": name, "level": "category"})
```

to:

```python
            yield (_cat_node_id(cat), "kegg_term", {"name": name, "level_kind": "category", "level": 0})
```

- [ ] **Step 2: Update existing KEGG tests**

In `tests/test_kegg_annotation_adapter.py`, replace every occurrence of `"level"` (as a KEGG property key) with `"level_kind"`. The 5 relevant sites are on lines 329, 331, 334, 339, 363, 407 per the earlier grep.

Specifically, change:

Line 329-334:

```python
def test_multi_kegg_nodes_have_level(multi_kegg):
    nodes = list(multi_kegg.get_nodes())
    levels = {n[2]["level"] for n in nodes}
    assert levels == {"ko", "pathway", "subcategory", "category"}
    for _, _, props in nodes:
        assert "level" in props
```

to:

```python
def test_multi_kegg_nodes_have_level(multi_kegg):
    nodes = list(multi_kegg.get_nodes())
    level_kinds = {n[2]["level_kind"] for n in nodes}
    assert level_kinds == {"ko", "pathway", "subcategory", "category"}
    expected_ints = {"category": 0, "subcategory": 1, "pathway": 2, "ko": 3}
    for _, _, props in nodes:
        assert "level_kind" in props
        assert "level" in props
        assert props["level"] == expected_ints[props["level_kind"]]
```

And change every `n[2].get("level") == "ko"` / `"pathway"` / etc. (lines 339, 363, 407) to `n[2].get("level_kind") == "ko"` / `"pathway"` / etc. Leave the rest of each test body unchanged.

- [ ] **Step 3: Run existing KEGG tests**

```bash
uv run pytest tests/test_kegg_annotation_adapter.py -v
```

Expected: all tests PASS, including the enhanced `test_multi_kegg_nodes_have_level`.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py tests/test_kegg_annotation_adapter.py
git commit -m "feat(kegg): rename level→level_kind, add integer level (0=category..3=ko)"
```

---

## Task 6: CyanorakRole / TigrRole / COG levels

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:944-984` (in `MultiCogRoleAnnotationAdapter.get_nodes`)

- [ ] **Step 1: Add `_role_depth` module-level helper**

In `functional_annotation_adapter.py`, add this helper near the other `_*_node_id` helpers (around line 770):

```python
def _role_depth(code: str, role_tree: dict[str, dict]) -> int:
    """
    Compute hierarchy depth for a CyanorakRole code by walking ``parent`` pointers.

    Roots (no parent) → 0. Secondary roles → 1. Sub-roles → 2.
    """
    depth = 0
    cur = code
    while True:
        entry = role_tree.get(cur)
        if entry is None or entry.get("parent") is None:
            return depth
        cur = entry["parent"]
        depth += 1
```

- [ ] **Step 2: Set `level` in the three emitters**

In `MultiCogRoleAnnotationAdapter.get_nodes` (lines 944–984), update the three `yield` statements:

Line 955-959 — change:

```python
            yield (
                _cog_cat_node_id(letter),
                "cog functional category",
                {"code": letter, "name": _clean_str(name)},
            )
```

to:

```python
            yield (
                _cog_cat_node_id(letter),
                "cog functional category",
                {"code": letter, "name": _clean_str(name), "level": 0},
            )
```

Line 966-970 — change:

```python
            yield (
                _cyanorak_role_node_id(code),
                "cyanorak role",
                {"code": code, "name": _clean_str(full_role_description(code, self.role_tree))},
            )
```

to:

```python
            yield (
                _cyanorak_role_node_id(code),
                "cyanorak role",
                {
                    "code": code,
                    "name": _clean_str(full_role_description(code, self.role_tree)),
                    "level": _role_depth(code, self.role_tree),
                },
            )
```

Line 978-982 — change:

```python
            yield (
                _tigr_role_node_id(code),
                "tigr role",
                {"code": code, "name": _clean_str(desc)},
            )
```

to:

```python
            yield (
                _tigr_role_node_id(code),
                "tigr role",
                {"code": code, "name": _clean_str(desc), "level": 0},
            )
```

- [ ] **Step 3: Verify adapter-level behaviour**

Add a quick spot-check by running the existing functional-annotation adapter tests:

```bash
uv run pytest tests/test_functional_annotation_adapter.py -v
```

Expected: PASS (no existing test reads `level` on these nodes; the new key is additive).

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py
git commit -m "feat(cog/cyanorak/tigr): emit level (parent-chain depth for CyanorakRole, 0 for flat)"
```

---

## Task 7: Pfam / PfamClan levels

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:1177-1213` (in `MultiPfamAnnotationAdapter.get_nodes`)

- [ ] **Step 1: Set `level = 1` on Pfam nodes**

In `functional_annotation_adapter.py`, update the two Pfam `yield` statements inside `MultiPfamAnnotationAdapter.get_nodes`:

Lines 1177-1181 (missing-entry branch) — change:

```python
                yield (
                    _pfam_node_id(pf_id),
                    "pfam",
                    {"name": pf_id, "short_name": pf_id},
                )
```

to:

```python
                yield (
                    _pfam_node_id(pf_id),
                    "pfam",
                    {"name": pf_id, "short_name": pf_id, "level": 1},
                )
```

Lines 1185-1192 (normal branch) — change:

```python
            yield (
                _pfam_node_id(pf_id),
                "pfam",
                {
                    "name": _clean_str(entry.description),
                    "short_name": _clean_str(entry.shortname),
                },
            )
```

to:

```python
            yield (
                _pfam_node_id(pf_id),
                "pfam",
                {
                    "name": _clean_str(entry.description),
                    "short_name": _clean_str(entry.shortname),
                    "level": 1,
                },
            )
```

- [ ] **Step 2: Set `level = 0` on PfamClan nodes**

Lines 1208-1212 — change:

```python
            yield (
                _pfam_clan_node_id(clan_acc),
                "pfam_clan",
                {"name": _clean_str(clan_name)},
            )
```

to:

```python
            yield (
                _pfam_clan_node_id(clan_acc),
                "pfam_clan",
                {"name": _clean_str(clan_name), "level": 0},
            )
```

- [ ] **Step 3: Run existing Pfam adapter tests**

```bash
uv run pytest tests/test_functional_annotation_adapter.py -k pfam -v
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py
git commit -m "feat(pfam): emit level (1 for Pfam, 0 for PfamClan) including clan-less Pfams"
```

---

## Task 8: KG-validity tests

**Files:**
- Create: `tests/kg_validity/test_ontology_level.py`

These tests require a rebuilt KG and a running Neo4j. They auto-skip if Neo4j is unreachable, matching the pattern in the other `kg_validity/test_*.py` files.

- [ ] **Step 1: Read an existing KG-validity test to match its fixture pattern**

Read `tests/kg_validity/test_post_import.py` to see how the `driver` / `session` fixture is wired. Use the same `conftest.py`-provided fixture in the new file.

- [ ] **Step 2: Write `tests/kg_validity/test_ontology_level.py`**

Create the file with:

```python
"""
KG-validity tests for the unified `level` property on ontology-term nodes.

Spec: docs/superpowers/specs/2026-04-12-ontology-level-design.md
Change note: docs/kg-changes/ontology-level.md

Requires a running Neo4j with the rebuilt graph.
"""
import pytest

pytestmark = pytest.mark.kg


ONTOLOGY_LABELS = [
    "BiologicalProcess",
    "MolecularFunction",
    "CellularComponent",
    "EcNumber",
    "KeggTerm",
    "CyanorakRole",
    "TigrRole",
    "CogFunctionalCategory",
    "Pfam",
    "PfamClan",
]


# GO BP depth → expected term count, per the parent spec verification section
# (observed on live KG as of 2026-04-12). Tests warn on ≤5% drift, fail on more.
EXPECTED_BP_LEVEL_DIST = {
    1: 16, 2: 98, 3: 274, 4: 565, 5: 721,
    6: 650, 7: 462, 8: 211, 9: 42, 10: 10, 11: 2,
}


def test_every_ontology_term_has_level(session):
    """Coverage: no ontology-term node has a null `level`."""
    for label in ONTOLOGY_LABELS:
        result = session.run(
            f"MATCH (t:{label}) WHERE t.level IS NULL RETURN count(t) AS n"
        ).single()
        assert result["n"] == 0, f"{label}: {result['n']} nodes missing `level`"


def test_flat_ontologies_all_level_zero(session):
    """TigrRole and CogFunctionalCategory are flat — every term at level 0."""
    tigr = session.run(
        "MATCH (t:TigrRole) RETURN count(t) AS n, min(t.level) AS lo, max(t.level) AS hi"
    ).single()
    assert tigr["n"] == 114
    assert tigr["lo"] == 0 and tigr["hi"] == 0

    cog = session.run(
        "MATCH (t:CogFunctionalCategory) RETURN count(t) AS n, min(t.level) AS lo, max(t.level) AS hi"
    ).single()
    assert cog["n"] == 26
    assert cog["lo"] == 0 and cog["hi"] == 0


def test_pfam_and_pfam_clan_fixed_levels(session):
    """Pfam level=1 (including clan-less), PfamClan level=0."""
    pfam = session.run(
        "MATCH (t:Pfam) RETURN count(t) AS n, min(t.level) AS lo, max(t.level) AS hi"
    ).single()
    assert pfam["n"] == 5471, f"unexpected Pfam count {pfam['n']}"
    assert pfam["lo"] == 1 and pfam["hi"] == 1

    clan = session.run(
        "MATCH (t:PfamClan) RETURN count(t) AS n, min(t.level) AS lo, max(t.level) AS hi"
    ).single()
    assert clan["n"] == 509, f"unexpected PfamClan count {clan['n']}"
    assert clan["lo"] == 0 and clan["hi"] == 0


def test_tree_ontology_roots(session):
    """CyanorakRole / EcNumber / KeggTerm have the expected root counts."""
    cyr = session.run("MATCH (t:CyanorakRole {level: 0}) RETURN count(t) AS n").single()
    assert cyr["n"] == 19, f"expected 19 CyanorakRole roots, got {cyr['n']}"

    ec = session.run("MATCH (t:EcNumber {level: 0}) RETURN count(t) AS n").single()
    assert ec["n"] == 7, f"expected 7 EcNumber level-0 terms, got {ec['n']}"

    kegg = session.run("MATCH (t:KeggTerm {level: 0}) RETURN count(t) AS n").single()
    assert kegg["n"] == 6, f"expected 6 KeggTerm category nodes, got {kegg['n']}"


def test_go_exactly_one_level_zero_per_namespace(session):
    """GO has exactly one canonical root per namespace (BP/MF/CC)."""
    for label in ("BiologicalProcess", "MolecularFunction", "CellularComponent"):
        res = session.run(
            f"MATCH (t:{label} {{level: 0}}) RETURN count(t) AS n"
        ).single()
        assert res["n"] == 1, f"{label}: expected 1 level-0 root, got {res['n']}"


def test_go_min_depth_invariant(session):
    """
    For every non-root GO term, its level equals 1 + min(parent.level) over
    `is_a`/`part_of` parent edges.
    """
    cypher = """
    MATCH (t)-[:Biological_process_is_a_biological_process|Biological_process_part_of_biological_process]->(p)
    WHERE t.level > 0
    WITH t, min(p.level) AS parent_min
    WHERE t.level <> parent_min + 1
    RETURN count(t) AS violations
    """
    assert session.run(cypher).single()["violations"] == 0


def test_go_bp_level_distribution_matches_spec(session):
    """
    BP depth histogram should match the spec's expected table. Allow ±5% drift
    per bucket (warn band); >5% fails as a regression.
    """
    rows = session.run(
        "MATCH (t:BiologicalProcess) RETURN t.level AS level, count(*) AS n ORDER BY level"
    ).data()
    observed = {r["level"]: r["n"] for r in rows}

    for depth, expected in EXPECTED_BP_LEVEL_DIST.items():
        actual = observed.get(depth, 0)
        tolerance = max(1, int(round(expected * 0.05)))
        assert abs(actual - expected) <= tolerance, (
            f"BP depth {depth}: expected {expected}±{tolerance}, got {actual}"
        )


def test_go_best_effort_flag_has_both_states(session):
    """
    At least one GO term has level_is_best_effort set, at least one doesn't.
    Sanity check that the flag is emitted AND sparse.
    """
    flagged = session.run(
        "MATCH (t:BiologicalProcess) WHERE t.level_is_best_effort = 'true' RETURN count(t) AS n"
    ).single()["n"]
    unflagged = session.run(
        "MATCH (t:BiologicalProcess) WHERE t.level_is_best_effort IS NULL RETURN count(t) AS n"
    ).single()["n"]
    assert flagged > 0, "no GO BP terms flagged level_is_best_effort"
    assert unflagged > 0, "all GO BP terms flagged — BFS is wrong"


def test_kegg_level_kind_to_int_mapping(session):
    """Every KeggTerm has both level_kind and level, and the mapping is consistent."""
    expected = {"category": 0, "subcategory": 1, "pathway": 2, "ko": 3}
    rows = session.run(
        "MATCH (t:KeggTerm) RETURN t.level_kind AS kind, t.level AS level"
    ).data()
    assert rows, "no KeggTerm nodes in graph"
    for row in rows:
        kind = row["kind"]
        lvl = row["level"]
        assert kind in expected, f"unknown level_kind {kind!r}"
        assert lvl == expected[kind], f"{kind}: expected level {expected[kind]}, got {lvl}"


def test_non_go_labels_never_best_effort(session):
    """level_is_best_effort is only ever emitted on GO labels."""
    non_go = [
        "EcNumber", "KeggTerm", "CyanorakRole", "TigrRole",
        "CogFunctionalCategory", "Pfam", "PfamClan",
    ]
    for label in non_go:
        res = session.run(
            f"MATCH (t:{label}) WHERE t.level_is_best_effort IS NOT NULL RETURN count(t) AS n"
        ).single()
        assert res["n"] == 0, f"{label}: {res['n']} nodes with unexpected level_is_best_effort"
```

- [ ] **Step 3: Run the new KG-validity tests against the rebuilt KG**

Before running, rebuild the graph:

```bash
docker compose down
docker compose up -d --build
# wait for Neo4j to come up
uv run pytest tests/kg_validity/test_ontology_level.py -v
```

Expected: all 10 tests PASS.

If any fail, investigate whether it's an adapter bug or a spec-distribution drift; adjust the affected adapter or expectations, then re-run.

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/test_ontology_level.py
git commit -m "test(kg_validity): add unified ontology-level verification tests"
```

---

## Task 9: Document the change in CLAUDE.md

**Files:**
- Modify: `CLAUDE.md` (under "Actual Neo4j labels" section)

- [ ] **Step 1: Add a bullet to the key graph facts list**

In `CLAUDE.md`, locate the "Key graph facts" bullet list (the one with entries like "Gene↔Protein linkage", "Experiment nodes", "Expression edges"). Add this bullet near the end of the list:

```markdown
- Unified ontology `level` property: every term on `BiologicalProcess`, `MolecularFunction`, `CellularComponent`, `EcNumber`, `KeggTerm`, `CyanorakRole`, `TigrRole`, `CogFunctionalCategory`, `Pfam`, `PfamClan` carries `level: int` with `0` = broadest/root. GO DAG terms whose min-path ≠ max-path to the namespace root additionally carry `level_is_best_effort = "true"` (sparse — property absent otherwise; check with `t.level_is_best_effort IS NOT NULL`). KEGG's original string `level` has been renamed to `level_kind` (`category|subcategory|pathway|ko`); the new integer `level` maps `category=0, subcategory=1, pathway=2, ko=3`. See `docs/kg-changes/ontology-level.md`.
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document unified ontology level convention in CLAUDE.md"
```

---

## Integration verification (final gate)

- [ ] **Step 1: Full unit test suite**

```bash
uv run pytest -m "not slow and not kg" -v
```

Expected: PASS.

- [ ] **Step 2: KG-validity suite (requires running Neo4j)**

```bash
uv run pytest -m kg -v
```

Expected: PASS (including the new `test_ontology_level.py`).

- [ ] **Step 3: Edge snapshot comparison**

Run the `/omics-edge-snapshot` skill before and after the rebuild to confirm no expression-edge regressions caused by unrelated drift. Level changes should not affect edge counts.

- [ ] **Step 4: Final commit if anything else surfaces**

If any of the above surfaces a small fix, commit it with a focused message referencing the plan. Otherwise, the feature is complete.
