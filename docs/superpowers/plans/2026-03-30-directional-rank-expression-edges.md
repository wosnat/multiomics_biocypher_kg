# Directional Rank on Expression Edges — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `rank_up` and `rank_down` precomputed properties to `Changes_expression_of` edges, ranking genes among significant-up or significant-down genes per experiment x timepoint.

**Architecture:** Two new Cypher blocks in the post-import script compute directional ranks after `expression_status` is set. `rank_up` ranks by |log2FC| descending among `significant_up` edges within each (experiment, time_point_order). `rank_down` does the same for `significant_down`. Non-significant edges and edges going the opposite direction get null. KG validity tests verify correctness.

**Tech Stack:** Neo4j Cypher (post-import script), pytest (KG validity tests)

**Spec reference:** Component 2 of `multiomics_research/docs/superpowers/specs/2026-03-30-gene-response-profile-and-analysis-utils-design.md`

**Note on edge type:** The spec mentions `AFFECTS_EXPRESSION` — the actual Neo4j label is `Changes_expression_of`.

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `config/schema_config.yaml` | Modify (~line 98) | Declare `rank_by_effect`, `rank_up`, `rank_down` on expression edge |
| `scripts/post-import.sh` | Modify (~line 228) | Add `rank_up` / `rank_down` Cypher computation block |
| `scripts/post-import.cypher` | Modify (~line 214) | Mirror the same Cypher (reference copy kept in sync) |
| `tests/kg_validity/test_expression.py` | Modify (append) | Add KG tests for `rank_up` / `rank_down` properties |
| `CLAUDE.md` | Modify | Document new properties in "Key graph facts" section |

---

### Task 1: Add `rank_up` / `rank_down` to schema and post-import scripts

**Files:**
- Modify: `config/schema_config.yaml:98` (expression edge properties)
- Modify: `scripts/post-import.sh:228` (after the existing `rank_by_effect` block)
- Modify: `scripts/post-import.cypher:214` (after the existing `rank_by_effect` block)

- [ ] **Step 1: Add rank properties to schema_config.yaml**

In `config/schema_config.yaml`, find the `experiment to gene expression association` edge properties block. After `expression_status: str`, add:

```yaml
    rank_by_effect: int   # post-import: rank by |log2FC| among all genes per experiment x timepoint (1 = strongest)
    rank_up: int          # post-import: rank by |log2FC| among significant_up genes per experiment x timepoint; null if not significant_up
    rank_down: int        # post-import: rank by |log2FC| among significant_down genes per experiment x timepoint; null if not significant_down
```

Note: `rank_by_effect` already exists in the post-import script but was never declared in the schema. Adding it here for completeness.

- [ ] **Step 2: Commit schema change**

```bash
git add config/schema_config.yaml
git commit -m "schema: declare rank_by_effect, rank_up, rank_down on expression edges"
```

- [ ] **Step 3: Add directional rank Cypher to `post-import.sh`**

Insert after the `rank_by_effect` block (after line 228), before `--- closest_ortholog_group_size ---`:

```bash
echo "--- rank_up (among significant_up per experiment x timepoint) ---"
cypher-shell <<'CYPHER'
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_up'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_up = i + 1
} IN TRANSACTIONS OF 10 ROWS;
CYPHER

echo "--- rank_down (among significant_down per experiment x timepoint) ---"
cypher-shell <<'CYPHER'
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_down'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_down = i + 1
} IN TRANSACTIONS OF 10 ROWS;
CYPHER
```

- [ ] **Step 4: Add the same Cypher to `post-import.cypher`**

Insert after the `rank_by_effect` block (after line 214), before `// closest_ortholog_group_size`:

```cypher
// rank_up: among significant_up edges per experiment + timepoint, rank by |log2FC| descending
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_up'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_up = i + 1
} IN TRANSACTIONS OF 10 ROWS;

// rank_down: among significant_down edges per experiment + timepoint, rank by |log2FC| descending
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_down'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_down = i + 1
} IN TRANSACTIONS OF 10 ROWS;
```

- [ ] **Step 5: Verify both files have identical Cypher logic**

Run:
```bash
diff <(grep -A5 'rank_up' scripts/post-import.sh | grep -v '^echo\|^cypher-shell\|^CYPHER') \
     <(grep -A5 'rank_up' scripts/post-import.cypher | grep -v '^//')
```

Expected: Only whitespace/comment differences, no logic divergence.

- [ ] **Step 6: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "feat: add rank_up and rank_down directional ranks on expression edges

Ranks genes by |log2FC| among significant_up (rank_up) or significant_down
(rank_down) genes within each experiment x timepoint. Non-significant edges
and edges in the opposite direction have null for the respective rank.
Computed after expression_status, same pattern as existing rank_by_effect."
```

---

### Task 2: Add KG validity tests for directional ranks

**Files:**
- Modify: `tests/kg_validity/test_expression.py` (append new test section)

- [ ] **Step 1: Add tests for `rank_up` / `rank_down` properties**

Append to `tests/kg_validity/test_expression.py`, after the `expression_status` test section:

```python
# ---------------------------------------------------------------------------
# Directional ranks (rank_up, rank_down — post-import computed)
# ---------------------------------------------------------------------------

def test_rank_up_only_on_significant_up(run_query):
    """rank_up must be null on edges that are not significant_up."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_up IS NOT NULL
          AND e.expression_status <> 'significant_up'
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have rank_up set but are not significant_up"
    )


def test_rank_down_only_on_significant_down(run_query):
    """rank_down must be null on edges that are not significant_down."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_down IS NOT NULL
          AND e.expression_status <> 'significant_down'
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have rank_down set but are not significant_down"
    )


def test_rank_up_covers_all_significant_up(run_query):
    """Every significant_up edge must have rank_up set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'significant_up'
          AND e.rank_up IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} significant_up edges are missing rank_up"
    )


def test_rank_down_covers_all_significant_down(run_query):
    """Every significant_down edge must have rank_down set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'significant_down'
          AND e.rank_down IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} significant_down edges are missing rank_down"
    )


def test_rank_up_starts_at_one(run_query):
    """Each experiment x timepoint must have a rank_up = 1 if any significant_up edges exist."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_up'
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_up) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_up not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_down_starts_at_one(run_query):
    """Each experiment x timepoint must have a rank_down = 1 if any significant_down edges exist."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_down'
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_down) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_down not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_up_contiguous(run_query):
    """rank_up values must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_up'
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_up) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_up: "
        f"{result[0]['examples']}"
    )


def test_rank_down_contiguous(run_query):
    """rank_down values must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_down'
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_down) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_down: "
        f"{result[0]['examples']}"
    )


def test_rank_up_rank_down_mutually_exclusive(run_query):
    """No edge should have both rank_up and rank_down set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_up IS NOT NULL AND e.rank_down IS NOT NULL
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have both rank_up and rank_down set"
    )


def test_not_significant_edges_have_no_directional_rank(run_query):
    """not_significant edges must have both rank_up and rank_down as null."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'not_significant'
          AND (e.rank_up IS NOT NULL OR e.rank_down IS NOT NULL)
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} not_significant edges have a directional rank set"
    )
```

- [ ] **Step 2: Run tests to verify they fail (no rank_up/rank_down in current graph)**

Run:
```bash
pytest tests/kg_validity/test_expression.py -v -k "rank_up or rank_down or directional_rank"
```

Expected: Tests that check for coverage (`test_rank_up_covers_all_significant_up`, `test_rank_down_covers_all_significant_down`) will FAIL because no edges have these properties yet. Tests for null-on-wrong-status will PASS vacuously. This confirms the tests are wired correctly.

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_expression.py
git commit -m "test: add KG validity tests for rank_up and rank_down directional ranks"
```

---

### Task 3: Add tests for existing `rank_by_effect` property

There are currently no KG tests for `rank_by_effect`. Add basic coverage alongside the directional rank tests.

**Files:**
- Modify: `tests/kg_validity/test_expression.py` (append)

- [ ] **Step 1: Add `rank_by_effect` tests**

Append to `tests/kg_validity/test_expression.py`, after the directional rank section:

```python
# ---------------------------------------------------------------------------
# rank_by_effect (post-import computed — all edges regardless of direction)
# ---------------------------------------------------------------------------

def test_rank_by_effect_populated(run_query):
    """Every expression edge must have rank_by_effect set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_by_effect IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} expression edges are missing rank_by_effect"
    )


def test_rank_by_effect_starts_at_one(run_query):
    """Each experiment x timepoint must have rank_by_effect = 1."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_by_effect) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_by_effect not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_by_effect_contiguous(run_query):
    """rank_by_effect must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_by_effect) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_by_effect: "
        f"{result[0]['examples']}"
    )
```

- [ ] **Step 2: Run tests against live graph**

Run:
```bash
pytest tests/kg_validity/test_expression.py -v -k "rank_by_effect"
```

Expected: PASS (these properties already exist in the deployed graph).

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_expression.py
git commit -m "test: add KG validity tests for existing rank_by_effect property"
```

---

### Task 4: Update CLAUDE.md documentation

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add `rank_up` / `rank_down` to the expression edge properties list**

In the "Key graph facts" section, find the line:
```
- Expression edges: ~188K `Changes_expression_of` edges ... Edge properties: `time_point` (str), `time_point_order` (int), `time_point_hours` (float), `log2_fold_change` (float), `adjusted_p_value` (float), `expression_direction` (str), `significant` (str), `expression_status` (str, post-import derived: `"significant_up"` | `"significant_down"` | `"not_significant"`).
```

Add after `expression_status` in the edge properties list:
```
`rank_by_effect` (int, post-import: rank by |log2FC| among all genes per experiment x timepoint), `rank_up` (int|null, post-import: rank by |log2FC| among significant_up genes per experiment x timepoint; null if not significant_up), `rank_down` (int|null, post-import: rank by |log2FC| among significant_down genes per experiment x timepoint; null if not significant_down)
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document rank_up and rank_down expression edge properties"
```

---

### Task 5: Rebuild KG and run tests

This task is manual — rebuild the Docker graph and validate.

- [ ] **Step 1: Rebuild KG with Docker**

```bash
docker compose down -v
docker compose up -d
```

Wait for all stages to complete. The post-import stage will compute `rank_up` and `rank_down`.

- [ ] **Step 2: Run all KG validity tests**

```bash
pytest tests/kg_validity/ -v
```

Expected: All tests PASS, including the new directional rank tests and existing tests.

- [ ] **Step 3: Spot-check directional ranks**

Run a quick sanity check query against Neo4j:

```bash
cypher-shell "
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.rank_up = 1
RETURN e.id AS experiment, r.time_point_order AS tp,
       g.locus_tag AS gene, r.log2_fold_change AS log2fc
ORDER BY e.id LIMIT 10
"
```

Expected: Each row shows the gene with the largest positive log2FC for its experiment x timepoint. The log2FC values should all be positive.

- [ ] **Step 4: Verify null partition**

```bash
cypher-shell "
MATCH ()-[r:Changes_expression_of]->()
RETURN
  count(r) AS total,
  count(r.rank_up) AS with_rank_up,
  count(r.rank_down) AS with_rank_down,
  count(CASE WHEN r.rank_up IS NOT NULL AND r.rank_down IS NOT NULL THEN 1 END) AS both,
  count(CASE WHEN r.expression_status = 'significant_up' THEN 1 END) AS sig_up,
  count(CASE WHEN r.expression_status = 'significant_down' THEN 1 END) AS sig_down
"
```

Expected: `with_rank_up` = `sig_up`, `with_rank_down` = `sig_down`, `both` = 0.

- [ ] **Step 5: Commit (if any adjustments were needed)**

Only if test failures required fixes — otherwise skip this step.
