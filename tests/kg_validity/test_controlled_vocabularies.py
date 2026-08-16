"""Live-graph vocabulary contract checks.

The CSV gate (tests/test_create_knowledge_graph.py) is structurally blind to
post-import-computed properties (transport_substrate_resolution,
tcdb_evidence_score_max, source_agreement, pfam_support, go_support), so this
suite is not redundant with it.

Query construction note: each `VocabEntry` already declares its own
`value_type`, so branching is done in Python on `entry.value_type ==
"string_array"` (build an UNWIND query) vs anything else (build a scalar
query), rather than calling Neo4j's `valueType()` function inside the Cypher
-- this project runs Neo4j 5.15 and that function's availability there is
uncertain, and the VocabEntry field makes the version-dependent call
unnecessary.
"""

import pytest

from multiomics_kg.utils.controlled_vocab import load_vocabularies, vocabularies_hash

pytestmark = pytest.mark.kg


@pytest.fixture(scope="module")
def declared():
    return load_vocabularies()


def test_every_vocabulary_has_a_node(run_query, declared):
    rows = run_query("MATCH (v:ControlledVocabulary) RETURN v.id AS id")
    assert {r["id"] for r in rows} == set(declared)


def _seen_values(run_query, entry):
    """Return the distinct values observed for one VocabEntry's property."""
    if entry.applies_to_kind == "node":
        match_clause = f"MATCH (n:{entry.applies_to}) WHERE n.`{entry.property}` IS NOT NULL"
        var = "n"
    else:
        match_clause = (f"MATCH ()-[r:{entry.applies_to}]->() "
                         f"WHERE r.`{entry.property}` IS NOT NULL")
        var = "r"
    if entry.value_type == "string_array":
        q = f"{match_clause} UNWIND {var}.`{entry.property}` AS v RETURN DISTINCT v AS v"
    else:
        q = f"{match_clause} RETURN DISTINCT {var}.`{entry.property}` AS v"
    return {r["v"] for r in run_query(q)}


def test_observed_values_are_declared(run_query, declared):
    problems = []
    for e in declared.values():
        if not e.closed:
            continue
        seen = _seen_values(run_query, e)
        undeclared = seen - set(e.values)
        if undeclared:
            problems.append(f"{e.id}: undeclared {sorted(undeclared)}")
        if not e.expected_empty and seen:
            unseen = set(e.values) - seen
            if unseen:
                problems.append(f"{e.id}: declared but absent {sorted(unseen)}")
        if e.expected_empty:
            assert not seen, f"{e.id} is expected_empty but the graph has {seen}"
    assert not problems, "\n".join(problems)


def test_every_sources_value_joins_a_data_source(run_query):
    """R2: sources values are DataSource ids."""
    rows = run_query("""
        MATCH ()-[r]->() WHERE r.sources IS NOT NULL
        UNWIND r.sources AS s
        WITH DISTINCT s
        WHERE NOT EXISTS { MATCH (d:DataSource) WHERE d.id = s }
        RETURN collect(s) AS orphans
    """)
    assert rows[0]["orphans"] == [], (
        f"sources values with no DataSource node: {rows[0]['orphans']}")


def test_no_native_boolean_properties_remain(run_query):
    """R5. Checks the properties this change converted."""
    for label, prop in [("TcdbFamily", "is_promiscuous"),
                        ("InterproEntry", "is_promiscuous")]:
        rows = run_query(
            f"MATCH (n:{label}) WHERE n.`{prop}` IS NOT NULL RETURN count(n) AS n")
        assert rows[0]["n"] == 0, f"{label}.{prop} should have been deleted"


def test_layer_a_router_edges_carry_no_properties(run_query):
    for rel in ("Interpro_entry_related_to_ec_number",
                "Interpro_entry_related_to_cazy_family"):
        rows = run_query(
            f"MATCH ()-[r:{rel}]->() WHERE size(keys(r)) > 0 RETURN count(r) AS n")
        assert rows[0]["n"] == 0, f"{rel} edges still carry properties"


def test_schema_info_carries_the_vocabulary_hash(run_query, declared):
    rows = run_query(
        "MATCH (s:Schema_info {id:'schema_info'}) "
        "RETURN s.controlled_vocabularies_hash AS h")
    assert rows[0]["h"] == vocabularies_hash(list(declared.values())), (
        "Schema_info hash does not match the shipped config -- the build and the "
        "checkout disagree about the vocabulary set.")
