"""KG validity: F1.2 + F1.3 annotation_quality + annotation_state."""

import pytest

pytestmark = pytest.mark.kg

VALID_STATES = {"no_evidence", "catch_all_only", "informative_single", "informative_multi"}
STATE_TO_QUALITY = {
    "no_evidence": 0,
    "catch_all_only": 1,
    "informative_single": 2,
    "informative_multi": 3,
}


def test_every_gene_has_valid_state(run_query):
    result = run_query(
        "MATCH (g:Gene) "
        "WHERE g.annotation_state IS NULL OR NOT g.annotation_state IN $valid "
        "RETURN count(*) AS n, collect(DISTINCT g.annotation_state)[..5] AS sample",
        valid=list(VALID_STATES),
    )
    assert result[0]["n"] == 0, f"Genes with invalid state: {result[0]}"


def test_quality_state_one_to_one(run_query):
    """annotation_quality must be the numeric encoding of annotation_state."""
    result = run_query(
        "MATCH (g:Gene) "
        "WITH g.annotation_state AS state, g.annotation_quality AS qual, count(*) AS n "
        "RETURN state, qual, n ORDER BY state, qual"
    )
    failures = []
    for row in result:
        expected = STATE_TO_QUALITY.get(row["state"])
        if row["qual"] != expected:
            failures.append(f"state={row['state']!r}: quality={row['qual']!r} (expected {expected!r}, n={row['n']})")
    assert not failures, "\n".join(failures)


def test_known_well_annotated_gene(run_query):
    """Spot-check: dnaA in MED4 (PMM0001) should be informative_multi."""
    result = run_query(
        "MATCH (g:Gene {locus_tag: 'PMM0001'}) "
        "RETURN g.annotation_state AS state, g.annotation_quality AS qual"
    )
    assert result, "PMM0001 not found"
    assert result[0]["state"] == "informative_multi"
    assert result[0]["qual"] == 3


def test_informative_subset_of_annotation_types(run_query):
    """For every gene, informative_annotation_types must be a subset of
    annotation_types — informativeness can only filter OUT, not add.
    Exception: 'reaction', 'transporter', 'cazy' are new informative-only
    sources not present in legacy annotation_types."""
    result = run_query(
        "MATCH (g:Gene) "
        "WITH g, [t IN g.informative_annotation_types "
        "         WHERE NOT t IN g.annotation_types AND NOT t IN $exempt | t] AS extra "
        "WHERE size(extra) > 0 "
        "RETURN count(*) AS n, collect(DISTINCT extra)[..3] AS sample",
        exempt=["reaction", "transporter", "cazy"],
    )
    assert result[0]["n"] == 0, f"Genes with extra informative-only types outside exempt list: {result[0]}"


def test_no_evidence_gene_has_empty_informative_types(run_query):
    result = run_query(
        "MATCH (g:Gene {annotation_state: 'no_evidence'}) "
        "WHERE size(g.informative_annotation_types) > 0 "
        "RETURN count(*) AS n"
    )
    assert result[0]["n"] == 0
