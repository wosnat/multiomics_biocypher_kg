# Chemistry slice-1 KG-side schema additions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add four schema properties (`Gene.reaction_count`, `Gene.metabolite_count`, `Metabolite.elements`, sparse pathway-level `KeggTerm.reaction_count` / `metabolite_count`) to support the chemistry MCP slice-1 explorer surface.

**Architecture:** Three of the four properties are post-import Cypher rollups (added to `scripts/post-import.sh` Group 3); one (`Metabolite.elements`) is a build-time list parsed from `Metabolite.formula` via the `chemparse` library. All schema declarations live in `config/schema_config.yaml`. `Gene.metabolite_count` is defined as **UNION across all gene-reaching paths** so the future TCDB-CAZy spec can extend the same property when transport edges land.

**Tech Stack:** Python 3.11, BioCypher, Neo4j 5, `chemparse` (new dep, pure Python), bash, Cypher.

**Spec:** [docs/superpowers/specs/2026-05-01-chemistry-slice1-kg-asks-design.md](../specs/2026-05-01-chemistry-slice1-kg-asks-design.md)

---

## File map

**Create:**
- `tests/test_metabolism_adapter_elements.py` — unit tests for `_parse_elements` helper.
- `tests/kg_validity/test_chemistry_slice1_rollups.py` — KG validity tests for the four new properties (separate file to keep `test_metabolism.py` focused on Reaction/Metabolite structural assertions).

**Modify:**
- `pyproject.toml` — add `chemparse>=0.3.1` to `dependencies`.
- `config/schema_config.yaml` — add 5 property declarations (gene: 2, metabolite: 1, kegg term: 2).
- `multiomics_kg/adapters/metabolism_adapter.py` — add `_parse_elements` helper + add `elements` to the metabolite props dict (around line 88-99).
- `scripts/post-import.sh` — add 4 Cypher statement blocks at the bottom of Group 3 (the existing "Metabolism rollups" subsection).
- `scripts/post-import.cypher` — mirror the same 4 statements (reference copy kept in sync per CLAUDE.md).
- `CLAUDE.md` — three appends in "Key graph facts" + one dependency note.
- `docs/kg-changes/metabolism-chemistry-layer.md` — append "Slice-1 explorer-coordinated additions" subsection.

**No changes:**
- `tests/kg_validity/snapshot_data.json` — these properties don't exist in the current snapshot, so no regen needed; the new tests are the regression net.
- `tests/test_kg_constants_drift.py` — verify in Task 8 whether it enumerates expected properties; only modify if it does.

---

## Task 1 — Add `chemparse` dependency

**Files:**
- Modify: `pyproject.toml` (the `dependencies = [...]` block, around line 16-50).

- [ ] **Step 1: Add chemparse to dependencies**

In `pyproject.toml`, find the `dependencies = [` block and add `"chemparse>=0.3.1",` alphabetically (between `"biopython>=1.86",` and `"bioregistry>=0.13.11",` is fine; use alphabetical or just append before the closing `]`).

Current snippet (around line 16-25):
```toml
dependencies = [
    "biocypher>=0.11.0",
    "biopython>=1.86",
    "bioregistry>=0.13.11",
    ...
]
```

Add `"chemparse>=0.3.1",` to the list.

- [ ] **Step 2: Run `uv sync` to install chemparse**

```bash
uv sync
```

Expected: chemparse pulled in successfully (~10 KB, pure Python).

- [ ] **Step 3: Verify chemparse imports and parses**

```bash
uv run python -c "import chemparse; print(chemparse.parse_formula('C10H12N5O13P3'))"
```

Expected: `{'C': 10.0, 'H': 12.0, 'N': 5.0, 'O': 13.0, 'P': 3.0}`

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "deps: add chemparse for Hill-notation formula parsing"
```

---

## Task 2 — Unit tests for `_parse_elements` helper

**Files:**
- Create: `tests/test_metabolism_adapter_elements.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_metabolism_adapter_elements.py`:

```python
"""Unit tests for the _parse_elements helper in metabolism_adapter."""
from __future__ import annotations

import pytest

from multiomics_kg.adapters.metabolism_adapter import _parse_elements


def test_parse_elements_water():
    assert _parse_elements("H2O") == ["H", "O"]


def test_parse_elements_complex_phosphate():
    assert _parse_elements("C10H12N5O13P3") == ["C", "H", "N", "O", "P"]


def test_parse_elements_two_letter_no_substring_clash():
    """Hill parsing must not split 'Na' into ['N','a'] or 'Cl' into ['C','l']."""
    assert _parse_elements("NaCl") == ["Cl", "Na"]


def test_parse_elements_kegg_charge_suffix():
    """KEGG uses '*N' suffix for charge/radical state; chemparse must tolerate it
    (parse the elements correctly, ignore the trailing *N)."""
    elts = _parse_elements("C42H44FeN8O8S2*4")
    assert "Fe" in elts
    assert "C" in elts
    assert "S" in elts
    assert "*" not in elts
    # Sorted, unique
    assert elts == sorted(set(elts))


def test_parse_elements_iron_sulfur_cluster():
    """Another KEGG-style charge suffix case."""
    elts = _parse_elements("Fe2S2*8")
    assert "Fe" in elts
    assert "S" in elts


def test_parse_elements_null_returns_empty():
    assert _parse_elements(None) == []


def test_parse_elements_empty_string_returns_empty():
    assert _parse_elements("") == []


def test_parse_elements_malformed_returns_empty():
    """Garbage input must not crash the adapter."""
    assert _parse_elements("???") == []


def test_parse_elements_returns_sorted_list():
    """Result is always sorted alphabetically by element symbol."""
    result = _parse_elements("ZnO")
    assert result == sorted(result)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/test_metabolism_adapter_elements.py -v
```

Expected: All tests FAIL with `ImportError: cannot import name '_parse_elements' from 'multiomics_kg.adapters.metabolism_adapter'`.

- [ ] **Step 3: Add `_parse_elements` helper to metabolism_adapter.py**

In `multiomics_kg/adapters/metabolism_adapter.py`, add the import and helper near the top, after the existing `_drop_nulls` helper (after line 35):

```python
import chemparse


def _parse_elements(formula: str | None) -> list[str]:
    """Sorted unique element symbols present in a Hill-notation formula.

    Empty list when formula is null/empty. Returns empty list (rather than
    raising) on any parse failure, so a malformed KEGG formula cannot break
    the build.
    """
    if not formula:
        return []
    try:
        return sorted(chemparse.parse_formula(formula).keys())
    except Exception:
        return []
```

The `import chemparse` line should go at the top of the file with the other imports (after `import logging`, line 12-13).

- [ ] **Step 4: Run unit tests to verify they pass**

```bash
uv run pytest tests/test_metabolism_adapter_elements.py -v
```

Expected: All 9 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/test_metabolism_adapter_elements.py multiomics_kg/adapters/metabolism_adapter.py
git commit -m "feat: add _parse_elements helper for Hill-notation formula parsing"
```

---

## Task 3 — Wire `elements` into the metabolite node yielder

**Files:**
- Modify: `multiomics_kg/adapters/metabolism_adapter.py` (around lines 84-100)

- [ ] **Step 1: Read the current metabolite yielder**

```bash
grep -n "for cpd_id, cpd in" /home/osnat/github/multiomics_biocypher_kg/multiomics_kg/adapters/metabolism_adapter.py
```

Confirm the loop at around line 84 emits a `props = _drop_nulls({...})` dict and yields `(node_id, "metabolite", props)`.

- [ ] **Step 2: Add `elements` to the props dict**

Find the `_drop_nulls({...})` call inside the metabolite yielder (currently lines ~88-99). Add one line:

```python
            props = _drop_nulls({
                "kegg_compound_id": cpd_id,
                "name": _clean_str(cpd.get("name", "")),
                "formula": _clean_str(cpd.get("formula")),
                "elements": _parse_elements(cpd.get("formula")),  # NEW
                "mass": cpd.get("mass"),
                "inchikey": _clean_str(cpd.get("inchikey")),
                "smiles": _clean_str(cpd.get("smiles")),
                "mnxm_id": cpd.get("mnxm_id"),
                "chebi_id": cpd.get("chebi_id"),
                "hmdb_id": cpd.get("hmdb_id"),
            })
```

Note: `_drop_nulls` will sparse-out the property for the ~31 metabolites without formula (it drops empty lists per the existing convention) — this is the correct semantic.

- [ ] **Step 3: Add an integration test for the metabolite yielder**

Append to `tests/test_metabolism_adapter_elements.py`:

```python
def test_metabolite_yielder_emits_elements_when_formula_present(tmp_path):
    """Smoke test: a metabolite with formula gets an 'elements' key in props."""
    import json
    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter

    # Minimal kegg_data.json fixture
    fake_data = {
        "reactions": {},
        "compounds": {
            "C00001": {"name": "H2O", "formula": "H2O"},
            "C99999": {"name": "Generic", "formula": None},
        },
    }
    fake_path = tmp_path / "kegg_data.json"
    fake_path.write_text(json.dumps(fake_data))

    adapter = MetabolismAdapter(kegg_data_path=fake_path, test_mode=False)
    nodes = list(adapter.get_nodes())

    # Find the H2O node — it should have elements=["H", "O"]
    h2o = next((p for nid, lbl, p in nodes if "C00001" in nid), None)
    assert h2o is not None
    assert h2o.get("elements") == ["H", "O"]

    # Find the no-formula node — elements should be absent (sparse, dropped by _drop_nulls)
    generic = next((p for nid, lbl, p in nodes if "C99999" in nid), None)
    assert generic is not None
    assert "elements" not in generic
```

If `MetabolismAdapter`'s constructor signature differs from `kegg_data_path=...`, adjust the fixture to match the actual constructor. Check by reading lines ~38-62 of `metabolism_adapter.py`.

- [ ] **Step 4: Run the new test**

```bash
uv run pytest tests/test_metabolism_adapter_elements.py -v
```

Expected: All tests pass, including the new yielder smoke test.

- [ ] **Step 5: Run the broader metabolism test suite to confirm no regression**

```bash
uv run pytest tests/test_metabolism_adapter.py tests/test_metabolism_smoke.py -v
```

Expected: All existing tests still pass.

- [ ] **Step 6: Commit**

```bash
git add tests/test_metabolism_adapter_elements.py multiomics_kg/adapters/metabolism_adapter.py
git commit -m "feat: emit Metabolite.elements from formula via Hill parsing"
```

---

## Task 4 — Schema declarations

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add `reaction_count` and `metabolite_count` to gene block**

Find the `gene:` block (around line 355). Locate the "Routing signals (pre-computed by post-import Cypher, not by adapter)" comment (around line 388). Add two lines after `expression_edge_count`:

Before (around lines 388-394):
```yaml
    # Routing signals (pre-computed by post-import Cypher, not by adapter)
    annotation_types: str[]           # ontology edge types with >0 edges, e.g. ["go_bp", "go_mf", "pfam"]
    expression_edge_count: int        # total Changes_expression_of edges
    significant_up_count: int         # expression edges where expression_status = 'significant_up'
    significant_down_count: int       # expression edges where expression_status = 'significant_down'
```

After:
```yaml
    # Routing signals (pre-computed by post-import Cypher, not by adapter)
    annotation_types: str[]           # ontology edge types with >0 edges, e.g. ["go_bp", "go_mf", "pfam"]
    expression_edge_count: int        # total Changes_expression_of edges
    significant_up_count: int         # expression edges where expression_status = 'significant_up'
    significant_down_count: int       # expression edges where expression_status = 'significant_down'
    reaction_count: int               # post-import: count of Gene_catalyzes_reaction edges (KG-A1)
    metabolite_count: int             # post-import: distinct Metabolite nodes reachable via any gene-reaching path; slice-1: catalysis only; UNION on TCDB-CAZy landing (KG-A2)
```

- [ ] **Step 2: Add `elements` to metabolite block**

Find the `metabolite:` block (around line 655). Add `elements: str[]` to the properties list:

Before:
```yaml
metabolite:
  ...
  properties:
    kegg_compound_id: str
    name: str
    formula: str
    mass: float
    inchikey: str
    ...
```

After (insert `elements: str[]` immediately after `formula: str`):
```yaml
metabolite:
  ...
  properties:
    kegg_compound_id: str
    name: str
    formula: str
    elements: str[]              # build-time: sorted unique element symbols, Hill-parsed via chemparse; sparse (absent on metabolites without formula) (KG-A3)
    mass: float
    inchikey: str
    ...
```

- [ ] **Step 3: Add `reaction_count` and `metabolite_count` to kegg term block**

Find the `kegg term:` block (around line 557). Add two properties:

Before:
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

After:
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
    reaction_count: int          # post-import (sparse, pathway-only): incoming Reaction_in_kegg_pathway edges (KG-A4)
    metabolite_count: int        # post-import (sparse, pathway-only): incoming Metabolite_in_pathway edges (KG-A4)
```

- [ ] **Step 4: Verify schema parses cleanly**

```bash
uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml').read()); print('OK')"
```

Expected: `OK`.

- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: declare Gene/Metabolite/KeggTerm chemistry-slice-1 properties"
```

---

## Task 5 — Add post-import Cypher rollups

**Files:**
- Modify: `scripts/post-import.sh` (Group 3, "Metabolism rollups" subsection at the bottom)
- Modify: `scripts/post-import.cypher` (mirror copy)

- [ ] **Step 1: Locate the metabolism rollups subsection**

```bash
grep -n "Metabolism rollups" /home/osnat/github/multiomics_biocypher_kg/scripts/post-import.sh
```

Expected: One match around line 723.

- [ ] **Step 2: Read the existing metabolism rollups block to confirm placement**

```bash
sed -n '720,765p' /home/osnat/github/multiomics_biocypher_kg/scripts/post-import.sh
```

Note the existing pattern: each rollup is wrapped in `CALL { ... } IN TRANSACTIONS OF N ROWS;`. The new statements go after the existing `Reaction.gene_count` / `Metabolite.gene_count` / `Organism_has_metabolite` / `OrganismTaxon.metabolite_count` / `OrganismTaxon.reaction_count` blocks, before the closing `CYPHER` heredoc terminator (around line 763).

- [ ] **Step 3: Add the four new Cypher statements**

In `scripts/post-import.sh`, after the last `} IN TRANSACTIONS OF 100 ROWS;` of the existing metabolism rollups (around line 762) and before the closing `CYPHER` line (around line 763), insert:

```cypher
// ── Chemistry slice-1 rollups (KG-A1, KG-A2, KG-A4) ───────────────────────

// KG-A1: Gene.reaction_count -- single-hop count of catalysis edges.
// count(r) returns 0 cleanly on no-match; no defaults pass needed.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_catalyzes_reaction]->()
  WITH g, count(r) AS rxn_count
  SET g.reaction_count = rxn_count
} IN TRANSACTIONS OF 1000 ROWS;

// KG-A2: Gene.metabolite_count -- 2-hop DISTINCT count, slice-1 form
// (catalysis only). The property is defined as UNION across all
// gene-reaching paths; on TCDB-CAZy landing the union expands to include
// the transport arm. See chemistry-slice-1 KG-A2 / TCDB-S3 coordination.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)
                 -[:Reaction_has_metabolite]->(m:Metabolite)
  WITH g, count(DISTINCT m) AS met_count
  SET g.metabolite_count = met_count
} IN TRANSACTIONS OF 1000 ROWS;

// KG-A4: KeggTerm pathway-level rollups (sparse on pathways only).
// level_kind = 'pathway' filter; KOs / categories left unset.
MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (r:Reaction)-[:Reaction_in_kegg_pathway]->(p)
  WITH p, count(r) AS rxn_count
  SET p.reaction_count = rxn_count
} IN TRANSACTIONS OF 100 ROWS;

MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p)
  WITH p, count(m) AS met_count
  SET p.metabolite_count = met_count
} IN TRANSACTIONS OF 100 ROWS;
```

- [ ] **Step 4: Mirror the same statements into `scripts/post-import.cypher`**

```bash
grep -n "Metabolism rollups\|Organism rollup props" /home/osnat/github/multiomics_biocypher_kg/scripts/post-import.cypher
```

Locate the equivalent metabolism-rollups section in `post-import.cypher` and append the same four Cypher statements at the end of that section. (The two files must stay in sync per CLAUDE.md.)

- [ ] **Step 5: Sanity-check the bash script syntax**

```bash
bash -n scripts/post-import.sh
```

Expected: No output (clean syntax).

- [ ] **Step 6: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "feat(post-import): add Gene/KeggTerm chemistry-slice-1 rollups"
```

---

## Task 6 — KG validity tests

**Files:**
- Create: `tests/kg_validity/test_chemistry_slice1_rollups.py`

- [ ] **Step 1: Write the test file**

Create `tests/kg_validity/test_chemistry_slice1_rollups.py`:

```python
"""KG validity tests for chemistry-slice-1 KG-side properties (KG-A1..A4).

Requires: a running Neo4j instance with the deployed graph (skipped otherwise).
Mark: @pytest.mark.kg
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.kg]


# ── KG-A1: Gene.reaction_count ────────────────────────────────────────────

def test_gene_reaction_count_present_on_all_genes(run_query):
    """Every Gene carries reaction_count (default 0)."""
    rows = run_query(
        "MATCH (g:Gene) RETURN count(g) AS total, count(g.reaction_count) AS with_prop"
    )
    assert rows[0]["with_prop"] == rows[0]["total"], (
        "reaction_count must be set on every Gene"
    )


def test_gene_reaction_count_total_positive(run_query):
    """Sum across all genes must be > 0 (KG has metabolism)."""
    n = run_query("MATCH (g:Gene) RETURN sum(g.reaction_count) AS total")[0]["total"]
    assert n > 0, f"Expected positive total reaction_count; got {n}"


def test_gene_reaction_count_consistent_with_edges(run_query):
    """Property matches actual outgoing edge count for sampled catalyzer genes."""
    rows = run_query("""
        MATCH (g:Gene) WHERE g.reaction_count > 0
        WITH g LIMIT 10
        OPTIONAL MATCH (g)-[r:Gene_catalyzes_reaction]->()
        WITH g, g.reaction_count AS prop, count(r) AS actual
        RETURN g.locus_tag AS lt, prop, actual
    """)
    assert len(rows) > 0, "No catalyzer genes found — KG may be empty"
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['lt']}: reaction_count={r['prop']} vs actual edges={r['actual']}"
        )


# ── KG-A2: Gene.metabolite_count ──────────────────────────────────────────

def test_gene_metabolite_count_present_on_all_genes(run_query):
    rows = run_query(
        "MATCH (g:Gene) RETURN count(g) AS total, count(g.metabolite_count) AS with_prop"
    )
    assert rows[0]["with_prop"] == rows[0]["total"]


def test_gene_metabolite_count_total_positive(run_query):
    n = run_query("MATCH (g:Gene) RETURN sum(g.metabolite_count) AS total")[0]["total"]
    assert n > 0


def test_gene_metabolite_count_consistent_with_2hop(run_query):
    """Property matches the actual 2-hop DISTINCT metabolite count for sampled genes."""
    rows = run_query("""
        MATCH (g:Gene) WHERE g.metabolite_count > 0
        WITH g LIMIT 10
        OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)
                       -[:Reaction_has_metabolite]->(m:Metabolite)
        WITH g, g.metabolite_count AS prop, count(DISTINCT m) AS actual
        RETURN g.locus_tag AS lt, prop, actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['lt']}: metabolite_count={r['prop']} vs actual 2-hop DISTINCT={r['actual']}"
        )


# ── KG-A3: Metabolite.elements ────────────────────────────────────────────

def test_metabolite_elements_populated_for_metabolites_with_formula(run_query):
    """Every metabolite that has a formula must have an elements list."""
    rows = run_query("""
        MATCH (m:Metabolite) WHERE m.formula IS NOT NULL
        RETURN count(m) AS total, count(m.elements) AS with_prop
    """)
    assert rows[0]["with_prop"] == rows[0]["total"], (
        f"Metabolites with formula but no elements: "
        f"{rows[0]['total'] - rows[0]['with_prop']}"
    )


def test_metabolite_elements_sparse_when_no_formula(run_query):
    """Metabolites without formula must NOT carry an elements property (sparse convention)."""
    n = run_query("""
        MATCH (m:Metabolite) WHERE m.formula IS NULL AND m.elements IS NOT NULL
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"Found {n} no-formula metabolites with stale elements property"


def test_metabolite_elements_h2o_correct(run_query):
    """Spot check: H2O must contain exactly H and O."""
    rows = run_query(
        "MATCH (m:Metabolite) WHERE m.formula = 'H2O' RETURN m.elements AS els LIMIT 1"
    )
    if rows:
        assert sorted(rows[0]["els"]) == ["H", "O"]


def test_metabolite_elements_no_substring_footgun_on_sodium(run_query):
    """A metabolite whose formula contains 'Na' but not standalone 'N' must NOT
    have 'N' in its elements list (the whole point of having a parsed list)."""
    n = run_query("""
        MATCH (m:Metabolite)
        WHERE m.formula CONTAINS 'Na'
          AND NOT m.formula =~ '.*N[^a].*'
          AND NOT m.formula ENDS WITH 'N'
          AND 'N' IN m.elements
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"{n} metabolites have false 'N' element from 'Na' substring"


def test_metabolite_elements_no_false_carbon_from_chlorine(run_query):
    """Metabolites whose formula contains 'Cl' but no real C must not have 'C' in elements."""
    n = run_query("""
        MATCH (m:Metabolite)
        WHERE m.formula CONTAINS 'Cl'
          AND NOT m.formula =~ '.*C[^l].*'
          AND NOT m.formula ENDS WITH 'C'
          AND 'C' IN m.elements
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0


# ── KG-A4: KeggTerm pathway-level rollups (sparse) ────────────────────────

def test_keggterm_pathway_rollups_present_on_pathways(run_query):
    """All pathway-level KeggTerms have both rollup props."""
    rows = run_query("""
        MATCH (k:KeggTerm) WHERE k.level_kind = 'pathway'
        RETURN count(k) AS total,
               count(k.reaction_count) AS with_rxn,
               count(k.metabolite_count) AS with_met
    """)
    assert rows[0]["with_rxn"] == rows[0]["total"], "all pathways need reaction_count"
    assert rows[0]["with_met"] == rows[0]["total"], "all pathways need metabolite_count"


def test_keggterm_non_pathway_rollups_sparse(run_query):
    """KOs and categories must NOT have the rollup props (sparse semantics)."""
    rows = run_query("""
        MATCH (k:KeggTerm) WHERE k.level_kind <> 'pathway'
        RETURN count(k.reaction_count) AS rxn_set,
               count(k.metabolite_count) AS met_set
    """)
    assert rows[0]["rxn_set"] == 0, (
        "Non-pathway KeggTerms must not carry reaction_count (sparse convention)"
    )
    assert rows[0]["met_set"] == 0


def test_keggterm_pathway_rollup_consistent_with_edges(run_query):
    """Spot-check: reaction_count matches actual incoming Reaction_in_kegg_pathway edges."""
    rows = run_query("""
        MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway' AND p.reaction_count > 0
        WITH p LIMIT 10
        OPTIONAL MATCH (r:Reaction)-[:Reaction_in_kegg_pathway]->(p)
        RETURN p.id AS id, p.reaction_count AS prop, count(r) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['id']}: reaction_count={r['prop']} vs actual={r['actual']}"
        )


def test_keggterm_pathway_metabolite_count_consistent(run_query):
    rows = run_query("""
        MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway' AND p.metabolite_count > 0
        WITH p LIMIT 10
        OPTIONAL MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p)
        RETURN p.id AS id, p.metabolite_count AS prop, count(m) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"]
```

- [ ] **Step 2: Run the new tests against the *current* (pre-rebuild) graph to confirm they fail**

```bash
uv run pytest tests/kg_validity/test_chemistry_slice1_rollups.py -v
```

Expected: All tests FAIL because the new properties do not exist yet in the deployed graph. (This proves the tests actually check the new properties — a green-on-day-one test is a worthless test.)

- [ ] **Step 3: Commit (test scaffold only; will go green after Task 7's rebuild)**

```bash
git add tests/kg_validity/test_chemistry_slice1_rollups.py
git commit -m "test(kg): KG-validity tests for chemistry-slice-1 rollups"
```

---

## Task 7 — Full Docker rebuild + verification

**Files:** none modified — this task runs the build pipeline and verifies green tests.

- [ ] **Step 1: Snapshot expression edges before rebuild**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py snapshot --name before_chemistry_slice1
```

(The `/omics-edge-snapshot` skill — see `.claude/skills/omics-edge-snapshot/SKILL.md` — wraps this script. Read the SKILL.md for current flags; the snapshot lives in `.claude/skills/omics-edge-snapshot/snapshots/`.) The point is to capture a pre-change baseline so we can prove no expression-edge regression.

- [ ] **Step 2: Build CSVs with the new adapter logic**

```bash
uv run python create_knowledge_graph.py
```

Expected: clean run; no errors related to `elements` or `chemparse`.

- [ ] **Step 3: Bring down deploy + app, run full Docker build**

```bash
docker compose down deploy app
docker compose up -d build
```

Wait for `build` to exit (check with `docker compose ps`). Then:

```bash
docker compose up -d import
docker compose logs -f import
```

Wait for the import container to finish (look for `import.status` artifact).

```bash
docker compose up -d post-process
docker compose logs -f post-process
```

Expected log output: `[timing]` lines for each Cypher group. The Group 3 timing should now include the four new statements.

```bash
docker compose up -d deploy app
```

- [ ] **Step 4: Wait for Neo4j to be ready**

```bash
for i in $(seq 1 60); do
  if docker exec deploy cypher-shell "RETURN 1" >/dev/null 2>&1; then
    echo "ready"; break
  fi
  sleep 1
done
```

- [ ] **Step 5: Run the new KG validity tests — they must now pass**

```bash
uv run pytest tests/kg_validity/test_chemistry_slice1_rollups.py -v
```

Expected: All 14 tests PASS.

- [ ] **Step 6: Run the full KG validity suite to confirm no regression**

```bash
uv run pytest -m kg -v
```

Expected: all green; no regressions in `test_metabolism.py`, `test_structure.py`, etc.

- [ ] **Step 7: Run unit tests**

```bash
uv run pytest -m "not slow and not kg" -v
```

Expected: all green.

- [ ] **Step 8: Snapshot expression edges after rebuild and diff**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py snapshot --name after_chemistry_slice1
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py compare \
  --before before_chemistry_slice1 --after after_chemistry_slice1
```

Expected: zero per-paper edge regressions. This PR doesn't touch expression edges; any change is a bug. (Adjust the exact CLI flags to match the script's current interface — see SKILL.md.)

- [ ] **Step 9: Spot-check the new properties in Cypher**

```bash
docker exec deploy cypher-shell "
  MATCH (g:Gene) WHERE g.reaction_count > 0
  RETURN g.locus_tag, g.reaction_count, g.metabolite_count
  ORDER BY g.reaction_count DESC LIMIT 5;
"
```

Expected: 5 rows with non-zero `reaction_count` and `metabolite_count`.

```bash
docker exec deploy cypher-shell "
  MATCH (m:Metabolite {formula: 'H2O'}) RETURN m.elements LIMIT 1;
"
```

Expected: `[\"H\", \"O\"]`.

```bash
docker exec deploy cypher-shell "
  MATCH (k:KeggTerm) WHERE k.level_kind = 'pathway'
  RETURN k.id, k.reaction_count, k.metabolite_count
  ORDER BY k.reaction_count DESC LIMIT 5;
"
```

Expected: 5 pathways with positive counts.

- [ ] **Step 10: No commit yet (this task only validates).**

If a verification step failed, treat it as a bug to fix in the relevant earlier task — do NOT add a "fix" commit on top of broken state.

---

## Task 8 — Documentation updates

**Files:**
- Modify: `CLAUDE.md` (in the "Key graph facts" section)
- Modify: `docs/kg-changes/metabolism-chemistry-layer.md`
- Verify: `tests/test_kg_constants_drift.py` (modify only if it enumerates expected properties)

- [ ] **Step 1: Update `CLAUDE.md` Gene routing-signals fact**

Find the line in `CLAUDE.md` that lists Gene routing signals (the long line under "Gene nodes carry ~27 properties." — search for `expression_edge_count`). Append `reaction_count`, `metabolite_count`, with brief notes:

```
... routing signals (annotation_types, expression_edge_count, significant_up_count, significant_down_count, reaction_count, metabolite_count, closest_ortholog_group_size, closest_ortholog_genera, cluster_membership_count, cluster_types — set by post-import Cypher, not adapter). reaction_count is the count of Gene_catalyzes_reaction edges; metabolite_count is the distinct count of Metabolite nodes reachable via any gene-reaching path (slice-1: catalysis only; UNION extends to transport on TCDB-CAZy landing).
```

- [ ] **Step 2: Update Metabolite block in `CLAUDE.md`**

Find the line about Metabolite nodes (search for `Metabolite nodes (~1.5K)`). Append `elements` to the property list:

```
... `chebi_id` (nullable), `hmdb_id` (nullable), `elements` (str[], sparse — sorted unique element symbols Hill-parsed from formula via chemparse; absent on the ~31 metabolites without formula). Computed post-import: gene_count, organism_count.
```

- [ ] **Step 3: Update KeggTerm pathway facts in `CLAUDE.md`**

There may already be a fact line about KeggTerm pathway nodes. If not, add to the relevant key-fact block:

```
- KeggTerm pathway-level nodes (level_kind = 'pathway', ~377 nodes) carry post-import `reaction_count: int` (incoming Reaction_in_kegg_pathway) and `metabolite_count: int` (incoming Metabolite_in_pathway). Sparse — only on pathways; KOs and categories have neither.
```

- [ ] **Step 4: Add chemparse to dependency note in `CLAUDE.md`**

If `CLAUDE.md` has a "Notes" or "Dependencies" section listing key non-obvious deps, add:

```
- `chemparse` — pure-Python Hill-notation formula parser (used in metabolism_adapter for Metabolite.elements)
```

If no such section exists, skip.

- [ ] **Step 5: Append slice-1 additions to chemistry-layer doc**

Append to `docs/kg-changes/metabolism-chemistry-layer.md`:

```markdown

## Slice-1 explorer-coordinated additions (2026-05-01)

Four properties added to support the chemistry MCP slice-1 explorer surface
(see [`docs/superpowers/specs/2026-05-01-chemistry-slice1-kg-asks-design.md`](../superpowers/specs/2026-05-01-chemistry-slice1-kg-asks-design.md)
and the explorer-side companion
`multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md`):

- **`Gene.reaction_count: int`** (post-import) — count of `Gene_catalyzes_reaction` edges.
  Routing signal for `gene_overview`.
- **`Gene.metabolite_count: int`** (post-import) — distinct Metabolite nodes
  reachable via *any* gene-reaching path. Slice-1 form is catalysis-only
  (2-hop `Gene → Reaction → Metabolite`); on TCDB-CAZy landing the same
  property expands to UNION across catalysis + transport. Coordinated with
  TCDB-CAZy spec via TCDB-S3.
- **`Metabolite.elements: list[str]`** (build-time) — sorted unique element
  symbols, Hill-parsed from `formula` via the `chemparse` library. Sparse —
  absent on the ~31 metabolites without `formula`. Eliminates the
  `formula CONTAINS 'N'` substring footgun (no false matches on `Na`/`Ne`/`Cl`).
- **`KeggTerm.reaction_count: int` / `KeggTerm.metabolite_count: int`** (post-import,
  sparse on `level_kind = 'pathway'`) — incoming `Reaction_in_kegg_pathway` /
  `Metabolite_in_pathway` edge counts per pathway. Filed for a Tier-2
  pathway-anchored explorer surface; not consumed by slice-1 directly.
```

- [ ] **Step 6: Check and update `tests/test_kg_constants_drift.py` if needed**

```bash
grep -n "expression_edge_count\|reaction_count\|metabolite_count\|elements" tests/test_kg_constants_drift.py 2>/dev/null
```

If the file exists and enumerates expected Gene/Metabolite/KeggTerm property names (typical pattern: a constant list compared against actual graph properties), add the four new properties to the relevant lists. If the file doesn't exist or doesn't enumerate properties, skip.

- [ ] **Step 7: Run unit + KG tests one more time**

```bash
uv run pytest -m "not slow and not kg" -v
uv run pytest -m kg -v
```

Expected: all green.

- [ ] **Step 8: Commit**

```bash
git add CLAUDE.md docs/kg-changes/metabolism-chemistry-layer.md
# only if changed:
git add tests/test_kg_constants_drift.py 2>/dev/null || true
git commit -m "docs: note chemistry-slice-1 KG-side property additions"
```

---

## Task 9 — Final verification

**Files:** none modified.

- [ ] **Step 1: Confirm all acceptance criteria**

Run through the 8 acceptance criteria from the spec doc:

1. `Gene.reaction_count` exists and matches edge count for sampled genes. ✓ (test_gene_reaction_count_consistent_with_edges)
2. `Gene.metabolite_count` exists and matches 2-hop DISTINCT for sampled genes. ✓ (test_gene_metabolite_count_consistent_with_2hop)
3. `Metabolite.elements` populated for ≥98% of metabolites + no substring footguns. ✓ (test_metabolite_elements_populated_for_metabolites_with_formula + the two no-substring-footgun tests)
4. KeggTerm pathway sparseness. ✓ (test_keggterm_pathway_rollups_present_on_pathways + test_keggterm_non_pathway_rollups_sparse)
5. `pytest -m kg -v` passes. ✓
6. `pytest -m "not slow and not kg" -v` passes. ✓
7. `/omics-edge-snapshot` before/after shows zero regression. ✓ (Task 7 step 8)
8. Schema declarations match property names in `schema_config.yaml`. ✓ (Task 4)

If any check fails, treat as a bug and revisit the relevant earlier task.

- [ ] **Step 2: Show the final git log for this branch**

```bash
git log --oneline main..HEAD
```

Expected: ~6-7 commits in topical order:
```
docs: note chemistry-slice-1 KG-side property additions
test(kg): KG-validity tests for chemistry-slice-1 rollups
feat(post-import): add Gene/KeggTerm chemistry-slice-1 rollups
schema: declare Gene/Metabolite/KeggTerm chemistry-slice-1 properties
feat: emit Metabolite.elements from formula via Hill parsing
feat: add _parse_elements helper for Hill-notation formula parsing
deps: add chemparse for Hill-notation formula parsing
spec: chemistry-slice-1 KG-side schema additions (KG-A1..A4)  # already on main
```

- [ ] **Step 3: Done. Ready for PR / merge.**

---

## Out-of-scope reminders

If during implementation any of these surface, **do not** add them:

- `Reaction_has_metabolite.role` (substrate vs product)
- `Metabolite.element_counts: map<str,int>`
- `Reaction.elements` / `OrganismTaxon.elements_observed`
- TCDB-CAZy properties — those belong in the TCDB-CAZy spec PR.
- The forward post-TCDB-landing UNION Cypher form for `Gene.metabolite_count` — written down in the spec for that future PR; not coded now.
