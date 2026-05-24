"""
KG-validity tests for the unified `level` property on ontology-term nodes.

Spec: docs/superpowers/specs/2026-04-12-ontology-level-design.md
Change note: docs/kg-changes/ontology-level.md

Requires a running Neo4j with the rebuilt graph. Uses the `run_query` fixture
from conftest.py, which returns query results as a list of dicts.
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


EXPECTED_BP_LEVEL_DIST = {
    1: 16, 2: 107, 3: 312, 4: 654, 5: 795,
    6: 700, 7: 501, 8: 205, 9: 42, 10: 11, 11: 2,
}


NON_GO_LABELS = [
    "EcNumber",
    "KeggTerm",
    "CyanorakRole",
    "TigrRole",
    "CogFunctionalCategory",
    "Pfam",
    "PfamClan",
]


GO_NAMESPACES = ["BiologicalProcess", "MolecularFunction", "CellularComponent"]


def test_every_ontology_term_has_level(run_query):
    """Every term node across all ten ontology labels must have a non-null level."""
    for label in ONTOLOGY_LABELS:
        rows = run_query(
            f"MATCH (t:{label}) WHERE t.level IS NULL RETURN count(t) AS n"
        )
        assert rows[0]["n"] == 0, (
            f"{rows[0]['n']} {label} nodes have null level"
        )


def test_flat_ontologies_all_level_zero(run_query):
    """TigrRole and CogFunctionalCategory are flat — all terms at level 0."""
    rows = run_query(
        "MATCH (t:TigrRole) RETURN count(t) AS n, min(t.level) AS mn, max(t.level) AS mx"
    )
    assert rows[0]["n"] == 114, f"TigrRole count {rows[0]['n']} != 114"
    assert rows[0]["mn"] == 0, f"TigrRole min level {rows[0]['mn']} != 0"
    assert rows[0]["mx"] == 0, f"TigrRole max level {rows[0]['mx']} != 0"

    rows = run_query(
        "MATCH (t:CogFunctionalCategory) RETURN count(t) AS n, min(t.level) AS mn, max(t.level) AS mx"
    )
    assert rows[0]["n"] == 26, f"CogFunctionalCategory count {rows[0]['n']} != 26"
    assert rows[0]["mn"] == 0, f"CogFunctionalCategory min level {rows[0]['mn']} != 0"
    assert rows[0]["mx"] == 0, f"CogFunctionalCategory max level {rows[0]['mx']} != 0"


def test_pfam_and_pfam_clan_fixed_levels(run_query):
    """Pfam domains are all level 1; PfamClans are all level 0."""
    rows = run_query(
        "MATCH (t:Pfam) RETURN count(t) AS n, min(t.level) AS mn, max(t.level) AS mx"
    )
    assert rows[0]["n"] == 5745, f"Pfam count {rows[0]['n']} != 5745"
    assert rows[0]["mn"] == 1, f"Pfam min level {rows[0]['mn']} != 1"
    assert rows[0]["mx"] == 1, f"Pfam max level {rows[0]['mx']} != 1"

    rows = run_query(
        "MATCH (t:PfamClan) RETURN count(t) AS n, min(t.level) AS mn, max(t.level) AS mx"
    )
    assert rows[0]["n"] == 516, f"PfamClan count {rows[0]['n']} != 516"
    assert rows[0]["mn"] == 0, f"PfamClan min level {rows[0]['mn']} != 0"
    assert rows[0]["mx"] == 0, f"PfamClan max level {rows[0]['mx']} != 0"


def test_tree_ontology_roots(run_query):
    """Tree ontologies have known root (level=0) counts."""
    rows = run_query("MATCH (t:CyanorakRole {level:0}) RETURN count(t) AS n")
    assert rows[0]["n"] == 19, f"CyanorakRole level-0 count {rows[0]['n']} != 19"

    rows = run_query("MATCH (t:EcNumber {level:0}) RETURN count(t) AS n")
    assert rows[0]["n"] == 7, f"EcNumber level-0 count {rows[0]['n']} != 7"

    rows = run_query("MATCH (t:KeggTerm {level:0}) RETURN count(t) AS n")
    assert rows[0]["n"] == 6, f"KeggTerm level-0 count {rows[0]['n']} != 6"


def test_go_exactly_one_level_zero_per_namespace(run_query):
    """Each GO namespace has exactly one level-0 (root) term."""
    for label in GO_NAMESPACES:
        rows = run_query(
            f"MATCH (t:{label} {{level:0}}) RETURN count(t) AS n"
        )
        assert rows[0]["n"] == 1, (
            f"{label} has {rows[0]['n']} level-0 terms, expected 1"
        )


def test_go_min_depth_invariant(run_query):
    """Every non-root GO term has level = 1 + min(parent.level) over is_a/part_of edges."""
    rows = run_query("""
        MATCH (t)-[:Biological_process_is_a_biological_process|Biological_process_part_of_biological_process]->(p)
        WHERE t.level > 0
        WITH t, min(p.level) AS parent_min
        WHERE t.level <> parent_min + 1
        RETURN count(t) AS violations
    """)
    assert rows[0]["violations"] == 0, (
        f"{rows[0]['violations']} BP terms violate the min-depth invariant"
    )


def test_go_bp_level_distribution_matches_spec(run_query):
    """BP level distribution matches spec within 5% (or +/- 1) per depth bucket."""
    rows = run_query(
        "MATCH (t:BiologicalProcess) RETURN t.level AS level, count(*) AS n ORDER BY level"
    )
    observed = {row["level"]: row["n"] for row in rows}
    for depth, expected in EXPECTED_BP_LEVEL_DIST.items():
        actual = observed.get(depth, 0)
        tolerance = max(1, int(round(expected * 0.05)))
        assert abs(actual - expected) <= tolerance, (
            f"BP level {depth}: observed {actual}, expected {expected} "
            f"(tolerance +/-{tolerance})"
        )


def test_go_best_effort_flag_has_both_states(run_query):
    """Some BP terms should be flagged best-effort, others not."""
    rows = run_query("""
        MATCH (t:BiologicalProcess)
        WHERE t.level_is_best_effort = 'true'
        RETURN count(t) AS n
    """)
    assert rows[0]["n"] > 0, (
        "Expected at least some BP terms with level_is_best_effort = 'true'"
    )

    rows = run_query("""
        MATCH (t:BiologicalProcess)
        WHERE t.level_is_best_effort IS NULL
        RETURN count(t) AS n
    """)
    assert rows[0]["n"] > 0, (
        "Expected at least some BP terms with level_is_best_effort IS NULL"
    )


def test_kegg_level_kind_to_int_mapping(run_query):
    """Every KeggTerm's level_kind maps deterministically to its integer level."""
    expected = {"category": 0, "subcategory": 1, "pathway": 2, "ko": 3}
    rows = run_query(
        "MATCH (t:KeggTerm) RETURN t.level_kind AS kind, t.level AS level"
    )
    assert len(rows) > 0, "No KeggTerm nodes found"
    for row in rows:
        kind = row["kind"]
        level = row["level"]
        assert kind in expected, (
            f"KeggTerm has unexpected level_kind {kind!r}"
        )
        assert level == expected[kind], (
            f"KeggTerm level_kind={kind!r} has level={level}, expected {expected[kind]}"
        )


def test_non_go_labels_never_best_effort(run_query):
    """Non-GO ontologies should never carry level_is_best_effort."""
    for label in NON_GO_LABELS:
        rows = run_query(
            f"MATCH (t:{label}) WHERE t.level_is_best_effort IS NOT NULL RETURN count(t) AS n"
        )
        assert rows[0]["n"] == 0, (
            f"{rows[0]['n']} {label} nodes have level_is_best_effort set; expected none"
        )
