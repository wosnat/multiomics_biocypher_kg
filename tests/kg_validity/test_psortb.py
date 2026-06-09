"""
PSORTb SubcellularLocalization KG-validity tests.

Validates the Phase-2 psortb integration in the live graph:
1. The 5 real PSORTb classes exist as SubcellularLocalization nodes (no "Unknown").
2. Nodes are flat (level 0) and carry post-import gene_count / organism_count.
3. Gene_has_subcellular_localization edges exist, 1:1, with a numeric `score`
   in the kept range, and a post-import rank_by_score.
4. No orphan edges; targets are always one of the 5 classes.
5. Gene.subcellular_localization routing string matches the linked node.
6. STRUCTURAL: psortb is deliberately NOT folded into Gene.annotation_types.
"""

import pytest

pytestmark = pytest.mark.kg

# The five real PSORTb (Gram-negative) localizations the KG keeps.
EXPECTED_CLASSES = {
    "Cytoplasmic",
    "CytoplasmicMembrane",
    "Periplasmic",
    "OuterMembrane",
    "Extracellular",
}


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------

def test_localization_nodes_are_the_real_classes(run_query):
    """Node psortb_ids must be a subset of the 5 real classes; ≥4 observed."""
    rows = run_query(
        "MATCH (n:SubcellularLocalization) RETURN n.psortb_id AS id"
    )
    ids = {r["id"] for r in rows}
    assert ids, "No SubcellularLocalization nodes found"
    assert ids <= EXPECTED_CLASSES, f"Unexpected localization ids: {ids - EXPECTED_CLASSES}"
    assert len(ids) >= 4, f"Only {len(ids)} localization classes present: {sorted(ids)}"


def test_no_unknown_localization_node(run_query):
    """The 'Unknown' sentinel must never be emitted as a node."""
    rows = run_query(
        "MATCH (n:SubcellularLocalization) WHERE n.psortb_id = 'Unknown' "
        "OR n.name = 'Unknown' RETURN count(n) AS cnt"
    )
    assert rows[0]["cnt"] == 0, "An 'Unknown' SubcellularLocalization node exists; it must be skipped"


def test_localization_nodes_are_flat_with_counts(run_query):
    """Flat ontology: level == 0 on every node; post-import counts present and ≥0."""
    rows = run_query(
        "MATCH (n:SubcellularLocalization) "
        "RETURN n.psortb_id AS id, n.level AS level, n.name AS name, "
        "n.gene_count AS gc, n.organism_count AS oc"
    )
    for r in rows:
        assert r["level"] == 0, f"{r['id']} has level {r['level']}, expected 0 (flat)"
        assert r["name"], f"{r['id']} missing name"
        assert r["gc"] is not None and r["gc"] >= 0, f"{r['id']} bad gene_count {r['gc']}"
        assert r["oc"] is not None and r["oc"] >= 0, f"{r['id']} bad organism_count {r['oc']}"
    assert sum(r["gc"] for r in rows) > 0, "No genes attached to any localization node"


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------

def test_localization_edges_exist(run_query):
    """Gene_has_subcellular_localization edges must exist in quantity."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_subcellular_localization]->(:SubcellularLocalization) "
        "RETURN count(r) AS cnt"
    )
    assert rows[0]["cnt"] > 1000, f"Only {rows[0]['cnt']} localization edges; expected > 1000"


def test_edge_score_is_numeric_and_in_kept_range(run_query):
    """Every edge carries a numeric score; kept calls sit in PSORTb's confident band."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_subcellular_localization]->() "
        "RETURN min(r.score) AS lo, max(r.score) AS hi, "
        "count(r) AS total, count(r.score) AS with_score"
    )
    row = rows[0]
    assert row["with_score"] == row["total"], (
        f"{row['total'] - row['with_score']} edges missing a score"
    )
    assert row["hi"] <= 10.0, f"Max score {row['hi']} exceeds PSORTb's 10.0 ceiling"
    assert row["lo"] >= 7.0, (
        f"Min kept score {row['lo']} below PSORTb's confident-call floor (~7.5)"
    )


def test_no_orphan_localization_edges(run_query):
    """Edge targets must always be one of the real classes (no dangling)."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_subcellular_localization]->(n:SubcellularLocalization) "
        "WHERE NOT n.psortb_id IN $classes RETURN count(r) AS bad",
        classes=sorted(EXPECTED_CLASSES),
    )
    assert rows[0]["bad"] == 0, f"{rows[0]['bad']} edges point at a non-class node"


def test_rank_by_score_present(run_query):
    """Post-import rank_by_score: present, ≥1, and rank 1 exists per populated node."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_subcellular_localization]->() "
        "RETURN count(r) AS total, count(r.rank_by_score) AS ranked, min(r.rank_by_score) AS lo"
    )
    row = rows[0]
    assert row["ranked"] == row["total"], "Some edges missing rank_by_score"
    assert row["lo"] == 1, f"Smallest rank is {row['lo']}, expected 1"


# ---------------------------------------------------------------------------
# Gene routing string + structural-not-functional invariant
# ---------------------------------------------------------------------------

def test_gene_routing_string_matches_edge(run_query):
    """Gene.subcellular_localization equals the linked node's psortb_id (1:1)."""
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_subcellular_localization]->(n:SubcellularLocalization) "
        "WHERE g.subcellular_localization <> n.psortb_id RETURN count(g) AS bad"
    )
    assert rows[0]["bad"] == 0, f"{rows[0]['bad']} genes have a mismatched subcellular_localization"


def test_psortb_not_folded_into_annotation_types(run_query):
    """Structural source: psortb must NOT appear in Gene.annotation_types."""
    rows = run_query(
        "MATCH (g:Gene) WHERE any(t IN g.annotation_types "
        "WHERE toLower(t) CONTAINS 'local' OR toLower(t) = 'psortb') "
        "RETURN count(g) AS bad"
    )
    assert rows[0]["bad"] == 0, (
        "psortb localization leaked into Gene.annotation_types (it is structural, not functional)"
    )
