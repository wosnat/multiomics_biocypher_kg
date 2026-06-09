"""
SignalP SignalPeptideType KG-validity tests.

Validates the Phase-2 signalp integration in the live graph:
1. The real SignalP 6.0 types exist as SignalPeptideType nodes (no "OTHER").
2. Nodes are flat (level 0) and carry post-import gene_count / organism_count.
3. Gene_has_signal_peptide_type edges exist, 1:1, with a numeric `probability`
   in [0,1], an int `cleavage_site` when present, and a post-import
   rank_by_probability.
4. No orphan edges; targets are always one of the real types.
5. Gene.signal_peptide_type routing string matches the linked node.
6. STRUCTURAL: signalp is deliberately NOT folded into Gene.annotation_types.
"""

import pytest

pytestmark = pytest.mark.kg

# The five real SignalP 6.0 signal-peptide types the KG keeps.
EXPECTED_TYPES = {"SP", "LIPO", "TAT", "TATLIPO", "PILIN"}


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------

def test_type_nodes_are_the_real_types(run_query):
    """Node signalp_ids must be a subset of the 5 real types; >=4 observed."""
    rows = run_query("MATCH (n:SignalPeptideType) RETURN n.signalp_id AS id")
    ids = {r["id"] for r in rows}
    assert ids, "No SignalPeptideType nodes found"
    assert ids <= EXPECTED_TYPES, f"Unexpected signalp ids: {ids - EXPECTED_TYPES}"
    assert len(ids) >= 4, f"Only {len(ids)} signal-peptide types present: {sorted(ids)}"


def test_no_other_type_node(run_query):
    """The 'OTHER' sentinel must never be emitted as a node."""
    rows = run_query(
        "MATCH (n:SignalPeptideType) WHERE n.signalp_id = 'OTHER' "
        "OR n.name = 'OTHER' RETURN count(n) AS cnt"
    )
    assert rows[0]["cnt"] == 0, "An 'OTHER' SignalPeptideType node exists; it must be skipped"


def test_type_nodes_are_flat_with_counts(run_query):
    """Flat ontology: level == 0 on every node; post-import counts present and >=0."""
    rows = run_query(
        "MATCH (n:SignalPeptideType) "
        "RETURN n.signalp_id AS id, n.level AS level, n.name AS name, "
        "n.gene_count AS gc, n.organism_count AS oc"
    )
    for r in rows:
        assert r["level"] == 0, f"{r['id']} has level {r['level']}, expected 0 (flat)"
        assert r["name"], f"{r['id']} missing name"
        assert r["gc"] is not None and r["gc"] >= 0, f"{r['id']} bad gene_count {r['gc']}"
        assert r["oc"] is not None and r["oc"] >= 0, f"{r['id']} bad organism_count {r['oc']}"
    assert sum(r["gc"] for r in rows) > 0, "No genes attached to any signal-peptide-type node"


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------

def test_type_edges_exist(run_query):
    """Gene_has_signal_peptide_type edges must exist in quantity."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_signal_peptide_type]->(:SignalPeptideType) "
        "RETURN count(r) AS cnt"
    )
    assert rows[0]["cnt"] > 1000, f"Only {rows[0]['cnt']} signal-peptide edges; expected > 1000"


def test_edge_probability_is_numeric_and_in_range(run_query):
    """Every edge carries a numeric probability in [0,1]."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_signal_peptide_type]->() "
        "RETURN min(r.probability) AS lo, max(r.probability) AS hi, "
        "count(r) AS total, count(r.probability) AS with_prob"
    )
    row = rows[0]
    assert row["with_prob"] == row["total"], (
        f"{row['total'] - row['with_prob']} edges missing a probability"
    )
    assert 0.0 <= row["lo"] <= row["hi"] <= 1.0, (
        f"probability out of [0,1]: lo={row['lo']} hi={row['hi']}"
    )


def test_cleavage_site_is_positive_int_when_present(run_query):
    """cleavage_site, where present, must be a positive integer."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_signal_peptide_type]->() "
        "WHERE r.cleavage_site IS NOT NULL "
        "RETURN count(r) AS cnt, min(r.cleavage_site) AS lo"
    )
    row = rows[0]
    assert row["cnt"] > 0, "No edges carry a cleavage_site (expected the common case)"
    assert row["lo"] >= 1, f"cleavage_site min {row['lo']} is not a positive residue position"


def test_no_orphan_type_edges(run_query):
    """Edge targets must always be one of the real types (no dangling)."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_signal_peptide_type]->(n:SignalPeptideType) "
        "WHERE NOT n.signalp_id IN $types RETURN count(r) AS bad",
        types=sorted(EXPECTED_TYPES),
    )
    assert rows[0]["bad"] == 0, f"{rows[0]['bad']} edges point at a non-type node"


def test_rank_by_probability_present(run_query):
    """Post-import rank_by_probability: present, >=1, rank 1 exists per populated node."""
    rows = run_query(
        "MATCH (:Gene)-[r:Gene_has_signal_peptide_type]->() "
        "RETURN count(r) AS total, count(r.rank_by_probability) AS ranked, "
        "min(r.rank_by_probability) AS lo"
    )
    row = rows[0]
    assert row["ranked"] == row["total"], "Some edges missing rank_by_probability"
    assert row["lo"] == 1, f"Smallest rank is {row['lo']}, expected 1"


# ---------------------------------------------------------------------------
# Gene routing string + structural-not-functional invariant
# ---------------------------------------------------------------------------

def test_gene_routing_string_matches_edge(run_query):
    """Gene.signal_peptide_type equals the linked node's signalp_id (1:1)."""
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_signal_peptide_type]->(n:SignalPeptideType) "
        "WHERE g.signal_peptide_type <> n.signalp_id RETURN count(g) AS bad"
    )
    assert rows[0]["bad"] == 0, f"{rows[0]['bad']} genes have a mismatched signal_peptide_type"


def test_signalp_not_folded_into_annotation_types(run_query):
    """Structural source: signalp must NOT appear in Gene.annotation_types."""
    rows = run_query(
        "MATCH (g:Gene) WHERE any(t IN g.annotation_types "
        "WHERE toLower(t) CONTAINS 'signal' OR toLower(t) = 'signalp') "
        "RETURN count(g) AS bad"
    )
    assert rows[0]["bad"] == 0, (
        "signalp leaked into Gene.annotation_types (it is structural, not functional)"
    )
