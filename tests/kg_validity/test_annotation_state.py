"""
Live invariants for the post-import Gene annotation_state / annotation_quality
computation (post-import.cypher F1.2) and its granular parallel
informative_annotation_types (F1.4).

These are invariants that hold on any healthy build, independent of what data
has landed — the complement of tests/kg_validity/capture_annotation_state.py,
which captures the (legitimately shifting) distribution as a baseline artifact,
and of tests/test_annotation_quality_buckets.py, which checks the bucket logic
statically in the script text. Here the *live result* is checked, so drift
between the two independently maintained post-import statements — or a stale
graph vs the current checkout — fails loudly.

The has_any_edge contract test is the one that would have caught the
2026-08-17 bug (genes whose only annotation was InterPro/NCBIfam miscounted
as no_evidence) before it shipped.
"""

import pytest


pytestmark = pytest.mark.kg


VALID_STATES = ["no_evidence", "catch_all_only", "informative_single", "informative_multi"]

# state <-> quality is a fixed 1:1 encoding (F1.2 sets both in one CASE)
STATE_TO_QUALITY = {
    "no_evidence": 0,
    "catch_all_only": 1,
    "informative_single": 2,
    "informative_multi": 3,
}

# The exact relationship list behind F1.2's has_any_edge predicate. Keep in
# sync with scripts/post-import.cypher / .sh (the SOURCE_BUCKETS block);
# tests/test_annotation_quality_buckets.py gates the two scripts against each
# other, this list gates the live graph against the checkout.
HAS_ANY_EDGE_RELS = (
    "Gene_involved_in_biological_process|Gene_enables_molecular_function"
    "|Gene_located_in_cellular_component"
    "|Gene_has_kegg_ko|Gene_has_pfam|Gene_catalyzes_ec_number"
    "|Gene_in_cog_category|Gene_has_cyanorak_role|Gene_has_tigr_role"
    "|Gene_catalyzes_reaction|Gene_has_tcdb_family|Gene_has_cazy_family"
    "|Gene_has_interpro_entry|Gene_has_ncbifam_family|Gene_has_merops_family"
)


def test_every_gene_has_valid_state_and_quality(run_query):
    """Every Gene carries annotation_state (4-value enum) and annotation_quality
    (0-3). A silently skipped post-import statement leaves these absent."""
    result = run_query(f"""
        MATCH (g:Gene)
        WHERE g.annotation_state IS NULL
           OR NOT g.annotation_state IN $states
           OR g.annotation_quality IS NULL
           OR NOT g.annotation_quality IN [0, 1, 2, 3]
        RETURN count(g) AS bad
    """, states=VALID_STATES)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Gene node(s) have a missing or out-of-vocabulary "
        f"annotation_state / annotation_quality"
    )


def test_state_quality_correspondence(run_query):
    """annotation_quality is the numeric encoding of annotation_state — the two
    are set by one CASE in F1.2 and must agree per gene, always."""
    result = run_query("""
        MATCH (g:Gene)
        WITH g, CASE g.annotation_state
                  WHEN 'no_evidence'        THEN 0
                  WHEN 'catch_all_only'     THEN 1
                  WHEN 'informative_single' THEN 2
                  WHEN 'informative_multi'  THEN 3
                END AS expected
        WHERE g.annotation_quality <> expected
        RETURN count(g) AS mismatches,
               collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatches"] == 0, (
        f"{result[0]['mismatches']} gene(s) where annotation_quality does not "
        f"encode annotation_state; examples: {result[0]['examples']}"
    )


def test_state_matches_recomputed_bucket_count(run_query):
    """Recompute F1.2's informative_count from informative_annotation_types
    (F1.4) and assert the clamp matches annotation_state.

    The 8 edge-based buckets use predicates identical to the granular list's
    entries (go = go_bp|go_mf|go_cc coarsened), so the list is a faithful
    witness for them. The 'role' bucket is the one deliberate divergence —
    F1.2 reads gene_category (a property), NOT the cog/cyanorak/tigr role
    edges the list reports — so role is recomputed from gene_category here.
    'brite' and 'merops' appear in the list but are NOT buckets and are
    ignored. Catches F1.2/F1.4 drift and stale-graph-vs-checkout skew.
    """
    result = run_query("""
        MATCH (g:Gene)
        WITH g, coalesce(g.informative_annotation_types, []) AS l
        WITH g,
             (CASE WHEN any(x IN l WHERE x IN ['go_bp', 'go_mf', 'go_cc']) THEN 1 ELSE 0 END
              + CASE WHEN 'kegg' IN l THEN 1 ELSE 0 END
              + CASE WHEN 'pfam' IN l THEN 1 ELSE 0 END
              + CASE WHEN 'ec' IN l THEN 1 ELSE 0 END
              + CASE WHEN g.gene_category IS NOT NULL AND g.gene_category <> 'Unknown'
                     THEN 1 ELSE 0 END
              + CASE WHEN 'reaction' IN l THEN 1 ELSE 0 END
              + CASE WHEN 'transporter' IN l THEN 1 ELSE 0 END
              + CASE WHEN 'cazy' IN l THEN 1 ELSE 0 END
              + CASE WHEN 'ncbifam' IN l THEN 1 ELSE 0 END) AS derived
        WHERE (derived >= 2 AND g.annotation_state <> 'informative_multi')
           OR (derived = 1 AND g.annotation_state <> 'informative_single')
           OR (derived = 0 AND g.annotation_state IN ['informative_single', 'informative_multi'])
        RETURN count(g) AS mismatches,
               collect({locus_tag: g.locus_tag, state: g.annotation_state,
                        derived: derived})[..5] AS examples
    """)
    assert result[0]["mismatches"] == 0, (
        f"{result[0]['mismatches']} gene(s) where annotation_state disagrees with "
        f"the bucket count recomputed from informative_annotation_types + "
        f"gene_category; examples: {result[0]['examples']}"
    )


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
    Exception: 'reaction' and 'transporter' are informative-only tokens
    (legacy uses 'tcdb' for the same TCDB edges and has no reaction token;
    verified against the 2026-08-18 baseline capture)."""
    result = run_query(
        "MATCH (g:Gene) "
        "WITH g, [t IN coalesce(g.informative_annotation_types, []) "
        "         WHERE NOT t IN coalesce(g.annotation_types, []) "
        "           AND NOT t IN $exempt | t] AS extra "
        "WHERE size(extra) > 0 "
        "RETURN count(*) AS n, collect(DISTINCT extra)[..3] AS sample",
        exempt=["reaction", "transporter"],
    )
    assert result[0]["n"] == 0, (
        f"Genes with informative-only types outside the exempt list: {result[0]}"
    )


def test_no_evidence_has_any_edge_contract(run_query):
    """The has_any_edge contract, checked live in both directions:
    no_evidence genes have ZERO edges of the counted relationship types (and
    an empty informative_annotation_types), catch_all_only genes have at
    least one. This is the test that would have caught the 2026-08-17
    InterPro/NCBIfam-only-gene miscount before it shipped."""
    leaked = run_query(f"""
        MATCH (g:Gene)
        WHERE g.annotation_state = 'no_evidence'
          AND (EXISTS {{ (g)-[:{HAS_ANY_EDGE_RELS}]->() }}
               OR size(coalesce(g.informative_annotation_types, [])) > 0)
        RETURN count(g) AS bad, collect(g.locus_tag)[..5] AS examples
    """)
    assert leaked[0]["bad"] == 0, (
        f"{leaked[0]['bad']} no_evidence gene(s) actually have counted "
        f"annotation edges (has_any_edge under-match — the 2026-08-17 bug "
        f"shape); examples: {leaked[0]['examples']}"
    )

    empty = run_query(f"""
        MATCH (g:Gene)
        WHERE g.annotation_state = 'catch_all_only'
          AND NOT EXISTS {{ (g)-[:{HAS_ANY_EDGE_RELS}]->() }}
        RETURN count(g) AS bad, collect(g.locus_tag)[..5] AS examples
    """)
    assert empty[0]["bad"] == 0, (
        f"{empty[0]['bad']} catch_all_only gene(s) have none of the counted "
        f"annotation edges (has_any_edge over-match); "
        f"examples: {empty[0]['examples']}"
    )
