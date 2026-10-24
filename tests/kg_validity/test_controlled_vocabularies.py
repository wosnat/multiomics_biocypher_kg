"""Live-graph vocabulary contract checks.

The CSV gate (tests/test_create_knowledge_graph.py) is structurally blind to
post-import-computed properties (e.g. TcdbFamily.tcdb_substrate_count,
Gene_has_merops_family.pfam_support), so this suite is not redundant with it.

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
    """Return the distinct values observed for one VocabEntry's property.

    An empty string means "absent", not a vocabulary member -- mirrors the
    CSV-scanning gate in tests/test_create_knowledge_graph.py, which strips
    BioCypher's quoted-empty (`''`) encoding of an unset string field before
    comparing against declared values. Excluded here for both the scalar and
    the string_array (post-UNWIND) query paths.
    """
    if entry.applies_to_kind == "node":
        match_clause = f"MATCH (n:{entry.applies_to}) WHERE n.`{entry.property}` IS NOT NULL"
        var = "n"
    else:
        match_clause = (f"MATCH ()-[r:{entry.applies_to}]->() "
                         f"WHERE r.`{entry.property}` IS NOT NULL")
        var = "r"
    if entry.value_type == "string_array":
        q = (f"{match_clause} UNWIND {var}.`{entry.property}` AS v "
             f"WITH DISTINCT v WHERE v <> '' RETURN v AS v")
    else:
        q = f"{match_clause} AND {var}.`{entry.property}` <> '' RETURN DISTINCT {var}.`{entry.property}` AS v"
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
        if e.exhaustive and not e.expected_empty and seen:
            unseen = set(e.values) - seen
            if unseen:
                problems.append(f"{e.id}: declared but absent {sorted(unseen)}")
        if e.expected_empty:
            assert not seen, f"{e.id} is expected_empty but the graph has {seen}"
    assert not problems, "\n".join(problems)


def test_min_size_lists_are_dense_and_long_enough(run_query, declared):
    """A vocabulary with min_size N asserts the string_array property is
    present on EVERY carrier with at least N elements — no null, no []."""
    problems = []
    for e in declared.values():
        if e.min_size is None:
            continue
        if e.applies_to_kind == "node":
            q = (f"MATCH (n:{e.applies_to}) WHERE n.`{e.property}` IS NULL "
                 f"OR size(n.`{e.property}`) < {int(e.min_size)} RETURN count(n) AS bad")
        else:
            q = (f"MATCH ()-[r:{e.applies_to}]->() WHERE r.`{e.property}` IS NULL "
                 f"OR size(r.`{e.property}`) < {int(e.min_size)} RETURN count(r) AS bad")
        bad = run_query(q)[0]["bad"]
        if bad:
            problems.append(f"{e.id}: {bad} carriers null or shorter than min_size={e.min_size}")
    assert not problems, "\n".join(problems)


def test_every_sources_value_joins_a_data_source(run_query):
    """R2: every sources value corresponds to a DataSource node whose `id`
    property is `data_source:<value>` (the DataSource id is prefixed, the
    edge-level sources value is bare)."""
    rows = run_query("""
        MATCH ()-[r]->() WHERE r.sources IS NOT NULL
        UNWIND r.sources AS s
        WITH DISTINCT s
        WHERE NOT EXISTS { MATCH (d:DataSource) WHERE d.id = 'data_source:' + s }
        RETURN collect(s) AS orphans
    """)
    assert rows[0]["orphans"] == [], (
        f"sources values with no DataSource node: {rows[0]['orphans']}")


def test_no_native_boolean_properties_remain(run_query):
    """R5. Deleted properties stay gone, and the two-state / categorical
    properties are strings on the live graph, not native bool.

    TCDB cleanup (2026-10-23): source_agreement / pfam_support / go_support left
    Gene_has_tcdb_family (removed by (label, property) — MEROPS keeps
    pfam_support); tcdb_call / tcdb_gate are the new categorical strings."""
    for label, prop in [("TcdbFamily", "is_promiscuous"),
                        ("InterproEntry", "is_promiscuous"),
                        ("TcdbFamily", "substrate_breadth")]:
        rows = run_query(
            f"MATCH (n:{label}) WHERE n.`{prop}` IS NOT NULL RETURN count(n) AS n")
        assert rows[0]["n"] == 0, f"{label}.{prop} should have been deleted"

    for prop in ("source_agreement", "pfam_support", "go_support"):
        rows = run_query(
            f"MATCH ()-[r:Gene_has_tcdb_family]->() "
            f"WHERE r.`{prop}` IS NOT NULL RETURN count(r) AS n")
        assert rows[0]["n"] == 0, f"Gene_has_tcdb_family.{prop} should have been deleted"

    for rel, prop in (("Gene_has_tcdb_family", "tcdb_call"),
                      ("Gene_has_tcdb_family", "tcdb_gate"),
                      ("Gene_has_merops_family", "pfam_support")):
        rows = run_query(
            f"MATCH ()-[r:{rel}]->() "
            f"WHERE r.`{prop}` IS NOT NULL RETURN r.`{prop}` AS v LIMIT 5")
        assert rows, f"no {rel} edges carry {prop}"
        for row in rows:
            assert isinstance(row["v"], str), (
                f"{rel}.{prop} is {type(row['v'])!r}, expected str "
                f"(R5: no native bool)")

    rows = run_query("MATCH (t:TcdbFamily) WHERE t.direction IS NOT NULL "
                     "RETURN t.direction AS v LIMIT 5")
    assert rows and all(isinstance(r["v"], str) for r in rows)


def test_layer_a_router_edges_carry_no_properties(run_query):
    # BioCypher always writes an `id` column on relationships, so
    # size(keys(r)) can never be 0 here -- assert the two deleted
    # properties specifically stay gone instead.
    for rel in ("Interpro_entry_related_to_ec_number",
                "Interpro_entry_related_to_cazy_family"):
        rows = run_query(
            f"MATCH ()-[r:{rel}]->() "
            "WHERE r.ambiguous IS NOT NULL OR r.source_db IS NOT NULL "
            "RETURN count(r) AS n")
        assert rows[0]["n"] == 0, f"{rel} edges still carry ambiguous/source_db"


def test_schema_info_carries_the_vocabulary_hash(run_query, declared):
    rows = run_query(
        "MATCH (s:Schema_info {id:'schema_info'}) "
        "RETURN s.controlled_vocabularies_hash AS h")
    assert rows[0]["h"] == vocabularies_hash(list(declared.values())), (
        "Schema_info hash does not match the shipped config -- the build and the "
        "checkout disagree about the vocabulary set.")


def test_value_descriptions_split_back_to_values(run_query, declared):
    """`value_descriptions` is str[] of '<value>: <one line>' (explorer B1);
    the value prefixes must equal the node's own `values`, in order."""
    rows = run_query(
        "MATCH (v:ControlledVocabulary) WHERE v.value_descriptions IS NOT NULL "
        "RETURN v.id AS id, v.values AS vals, v.value_descriptions AS vd")
    if not any(e.value_descriptions for e in declared.values()):
        pytest.skip("no value_descriptions declared")
    if not rows:
        pytest.skip("value_descriptions declared but no live node carries them yet (rebuild pending)")
    for r in rows:
        assert [s.split(": ", 1)[0] for s in r["vd"]] == r["vals"], r["id"]
        assert all(s.split(": ", 1)[1].strip() for s in r["vd"]), r["id"]
