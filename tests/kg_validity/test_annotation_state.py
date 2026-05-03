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
        {"valid": list(VALID_STATES)},
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
