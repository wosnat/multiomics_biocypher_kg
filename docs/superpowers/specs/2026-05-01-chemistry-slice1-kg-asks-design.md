---
title: Chemistry slice-1 KG-side schema additions (KG-A1..A4)
date: 2026-05-01
status: design — ready for implementation
companion: multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md
coordinated_with: docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md (TCDB-S3)
---

# Chemistry slice-1 KG-side schema additions

## Goal

Add four small schema properties so the chemistry MCP slice-1 surface in
`multiomics_explorer` ships cleanly. All four are mechanical extensions of
patterns that already exist in `scripts/post-import.sh` and the metabolism
adapter — no new architectural surface.

The explorer side is forward-compatible via `coalesce(g.<prop>, 0)`, so this
PR can land before, in parallel with, or slightly after the explorer's slice-1
PRs without blocking either side.

## Scope

In scope:

- **KG-A1** — `Gene.reaction_count: int`, post-import rollup (single-hop count
  of `Gene_catalyzes_reaction` edges).
- **KG-A2** — `Gene.metabolite_count: int`, post-import rollup. Slice-1 form
  is 2-hop catalysis-only (`Gene → Reaction → Metabolite`, distinct).
  **Defined semantically as UNION across all gene-reaching paths**, so the
  TCDB-CAZy spec extends the same property when transport edges land. No
  second property like `Gene.transport_metabolite_count` is introduced.
- **KG-A3** — `Metabolite.elements: list[str]`, build-time Hill-parsed
  element-symbol list, populated in `metabolism_adapter.py` from `formula`.
  Eliminates the `formula CONTAINS 'N'` substring footgun.
- **KG-A4** — `KeggTerm.reaction_count: int` and `KeggTerm.metabolite_count: int`,
  pathway-level post-import rollups (sparse: filtered to
  `level_kind = 'pathway'`).

Out of scope (explicitly):

- TCDB-CAZy properties (`Gene.tcdb_family_count`, `Gene.cazy_family_count`,
  `TcdbFamily.metabolite_count`, etc.) — owned by the TCDB-CAZy spec.
- `Metabolite.element_counts: map<str,int>` — explorer explicitly defers;
  slice-1 needs presence only.
- `Reaction.elements` / `OrganismTaxon.elements_observed` rollups — derivable
  at query time, not worth materializing.
- The forward "Cypher with TCDB transport arm" for KG-A2 — written down for
  the day TCDB lands; not coded now.
- `Reaction_has_metabolite.role` (substrate vs product) — needs external
  authoritative direction (Rhea / curated thermodynamics); out of scope per
  the explorer ask doc.

## Coordination

- **TCDB-CAZy** — `Gene.metabolite_count` UNION semantics is already locked in
  writing on both sides:
  - `multiomics_biocypher_kg/docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md`
    §"Coordination with chemistry-slice-1"
  - `multiomics_biocypher_kg/docs/kg-changes/tcdb-cazy-ontologies.md` Gene
    property table
  - The explorer-side companion doc (cross-spec point TCDB-S3).
  When TCDB-CAZy lands, that PR replaces the slice-1 KG-A2 statement with a
  union form (catalysis arm + transport arm) and grows `metabolite_count`
  automatically. Consumers read one property regardless of source path.
- **Explorer slice-1** — independent. `coalesce(g.<prop>, 0)` keeps explorer
  code timing-resilient.

## Components & changes

### 1. `pyproject.toml`

Add `chemparse>=0.3.1` to project deps. Pure-Python (~10 KB), actively
maintained on PyPI. Handles nested parentheses and square brackets (added
2024) — useful if KEGG's formula notation grows.

Note: chemparse does **not** validate symbols against the real periodic
table — `parse_formula('R')` returns `{'R': 1.0}` and `parse_formula('(C5H10O5)n')`
captures lowercase `n` as an "element". Today this is non-issue because
the deployed KEGG cache has zero formulas with `R`/`X`/polymer placeholders;
if KEGG ever ships such formulas, a follow-up should wrap `_parse_elements`
with a real-element-symbol allow-list.

### 2. `config/schema_config.yaml`

Three edits:

```yaml
gene:
  ...
  properties:
    ...
    # Routing signals (pre-computed by post-import Cypher, not by adapter)
    ...
    reaction_count: int            # post-import: Gene_catalyzes_reaction edges (KG-A1)
    metabolite_count: int          # post-import: distinct Metabolite nodes
                                   # reachable via any gene-reaching path
                                   # (slice-1: catalysis only; UNION on TCDB landing) (KG-A2)
```

```yaml
metabolite:
  ...
  properties:
    ...
    elements: str[]                # build-time: sorted unique element symbols
                                   # Hill-parsed from formula via chemparse;
                                   # empty list when formula is null (KG-A3)
```

```yaml
kegg term:
  ...
  properties:
    ...
    reaction_count: int            # post-import (sparse, pathway-only):
                                   # incoming Reaction_in_kegg_pathway edges (KG-A4)
    metabolite_count: int          # post-import (sparse, pathway-only):
                                   # incoming Metabolite_in_pathway edges (KG-A4)
```

`reaction_count` and `metabolite_count` on `kegg term:` are declared as
properties but only populated on `level_kind = 'pathway'` nodes — same
sparseness pattern as the existing `level_is_best_effort`.

### 3. `multiomics_kg/adapters/metabolism_adapter.py`

Module-level helper near the top of the file:

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

In the `metabolite` node yielder (around lines 84-100 of the current file),
add `"elements": _parse_elements(cpd.get("formula"))` to the props dict.
Do NOT route through `_clean_str` — `elements` is a list of element symbols,
not a string.

**Sparse-output behavior.** The existing `_drop_nulls` helper (line 33-35
of the current adapter) drops keys whose value is `None` *or* `[]`, by
convention. This means the 31 metabolites without `formula` will not carry
an `elements` property at all (rather than carrying `[]`). This is the
correct semantics: those metabolites have no formula evidence, so explorer
queries like `WHERE 'N' IN m.elements` will correctly filter them out
(null comparison is false in Cypher predicates).

Schema-side, `elements: str[]` is declared but sparse — same pattern as
nullable Metabolite props (`mnxm_id`, `chebi_id`, `hmdb_id`).

### 4. `scripts/post-import.sh`

Three new statement blocks added to **Group 3** (the heavy `CALL { } IN
TRANSACTIONS` group), in the existing "Metabolism rollups" subsection (where
`Reaction.gene_count` and `Metabolite.gene_count` already live, lines
~723-762):

```cypher
// KG-A1: Gene.reaction_count -- single-hop count of catalysis edges.
// count(r) returns 0 cleanly on no-match, no defaults pass needed.
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

// KG-A4: pathway-level rollups, sparse on KeggTerm pathways only
// (level_kind = 'pathway' filter; KOs/categories left without these props).
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

### 5. `scripts/post-import.cypher`

Mirror the same Cypher statements (the reference copy kept in sync per
CLAUDE.md). Run `scripts/post-import-validate.sh` against the deployed graph
before and after to verify the dump is byte-identical for already-existing
properties (only the four new properties should appear in the diff).

## Build sequencing

Single PR, single full Docker rebuild:

1. `uv sync` — pulls in chemparse.
2. `uv run python create_knowledge_graph.py` — rebuilds CSVs with
   `Metabolite.elements`.
3. `docker compose down deploy app && docker compose up -d build import
   post-process deploy app` — full rebuild + import + post-process.
4. `pytest -m kg -v` — KG validity assertions (new tests below).
5. `pytest -m "not slow and not kg" -v` — unit tests including
   `_parse_elements` cases.
6. `/omics-edge-snapshot` before/after — must show 0 expression-edge
   regression (this PR doesn't touch expression edges; pure additions).

## Tests

### Unit tests (`tests/test_metabolism_adapter_elements.py` or extension of an
existing metabolism-adapter test file)

Cases:

- `_parse_elements("H2O")` → `["H", "O"]`
- `_parse_elements("C10H12N5O13P3")` → `["C", "H", "N", "O", "P"]`
- `_parse_elements("NaCl")` → `["Cl", "Na"]` (no false `["C", "H", "N"]` —
  the two-letter-element regression test).
- `_parse_elements("C42H44FeN8O8S2*4")` → contains `Fe`, no `*`, contains
  `C`, `H`, `N`, `O`, `S` (KEGG charge/radical suffix tolerance).
- `_parse_elements(None)` → `[]`
- `_parse_elements("")` → `[]`
- `_parse_elements("???")` → `[]` (malformed input must not crash).

### KG validity tests (`tests/kg_validity/test_metabolism_rollups.py` or
extension of `test_structure.py`)

Cases (all `@pytest.mark.kg`):

- `Gene.reaction_count` exists on all genes; `sum > 0`.
- Spot-check: `g.reaction_count == count(g-[:Gene_catalyzes_reaction]->())`
  for ≥5 sampled catalyzer genes.
- `Gene.metabolite_count` exists on all genes; `sum > 0`.
- Spot-check: `g.metabolite_count == count(DISTINCT m)` from explicit 2-hop
  traversal for ≥5 sampled genes.
- `Metabolite.elements` populated for ≥98% of metabolites (the ~31 without
  `formula` carry no `elements` property — sparse convention).
- No substring footgun: for every metabolite *with* an `elements` list,
  `'N' IN m.elements` iff `m.formula` actually contains the element N (not
  a substring of `Na`/`Ne`/`Cl`).
- KeggTerm pathway sparseness: `level_kind = 'pathway'` nodes carry
  `reaction_count`/`metabolite_count`; `level_kind = 'ko'` nodes do not.
- KeggTerm pathway consistency: spot-check `p.reaction_count == count(r)`
  via explicit edge query.

### Drift tests

If `tests/test_kg_constants_drift.py` enumerates expected Gene / Metabolite /
KeggTerm property names, add the four new ones. Verify during implementation.

### Snapshot

`tests/kg_validity/snapshot_data.json` doesn't sample any of these (they
don't exist yet); no regen required. The dedicated rollup tests above are
the regression net.

## Documentation updates

- **`CLAUDE.md`** — three appends in the "Key graph facts" section:
  - Gene routing-signals list: add `reaction_count` and `metabolite_count`
    (with UNION-semantics note for the latter).
  - Metabolite block: note `elements: list[str]` (Hill-parsed via chemparse).
  - KeggTerm/pathway: note sparse `reaction_count` / `metabolite_count` on
    `level_kind = 'pathway'` nodes.
  - Add `chemparse` to the dependency note.
- **`docs/kg-changes/metabolism-chemistry-layer.md`** — append a "Slice-1
  explorer-coordinated additions (2026-05-01)" subsection summarizing the
  four properties; cross-link the explorer-side asks doc.
- **Memory** — no entry. Properties are derivable from the schema file; the
  UNION semantics and chemparse choice are in this spec.
- **Plan file** — none needed; spec doc is the work record for a single
  small PR.

## Acceptance criteria

A reviewer can verify this PR is done by:

1. `Gene.reaction_count` exists and matches edge count for sampled genes.
2. `Gene.metabolite_count` exists and matches the 2-hop DISTINCT count for
   sampled catalyzer genes.
3. `Metabolite.elements` is populated for ≥98% of metabolites (sparse —
   the ~31 no-formula metabolites carry no `elements` property), and
   `'N' IN m.elements ↔ formula` actually contains element N for spot-checks
   on metabolites that do carry the property.
4. `KeggTerm {level_kind:'pathway'}` nodes have `reaction_count` and
   `metabolite_count`; KO/category nodes have neither.
5. `pytest -m kg -v` passes.
6. `pytest -m "not slow and not kg" -v` passes (chemparse unit tests
   included).
7. `/omics-edge-snapshot` before/after shows zero expression-edge regression.
8. Schema declarations match property names in `schema_config.yaml`.

## References

- Explorer-side companion: `multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-chemistry-slice1-asks.md`
- TCDB-CAZy KG-side spec (UNION coordination): `docs/superpowers/specs/2026-05-01-tcdb-cazy-ontologies-design.md`
- TCDB-CAZy explorer-facing summary: `docs/kg-changes/tcdb-cazy-ontologies.md`
- Chemistry layer reference: `docs/kg-changes/metabolism-chemistry-layer.md`
- Metabolism scaffold spec: `docs/superpowers/specs/2026-04-30-metabolite-reactions-scaffold-revised.md`
