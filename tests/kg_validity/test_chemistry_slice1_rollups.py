"""KG validity tests for chemistry-slice-1 KG-side properties (KG-A1..A4).

Requires: a running Neo4j instance with the deployed graph (skipped otherwise).
Mark: @pytest.mark.kg
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.kg]


# ── KG-A1: Gene.reaction_count ────────────────────────────────────────────

def test_gene_reaction_count_present_on_all_genes(run_query):
    """Every Gene carries reaction_count (default 0)."""
    rows = run_query(
        "MATCH (g:Gene) RETURN count(g) AS total, count(g.reaction_count) AS with_prop"
    )
    assert rows[0]["with_prop"] == rows[0]["total"], (
        "reaction_count must be set on every Gene"
    )


def test_gene_reaction_count_total_positive(run_query):
    """Sum across all genes must be > 0 (KG has metabolism)."""
    n = run_query("MATCH (g:Gene) RETURN sum(g.reaction_count) AS total")[0]["total"]
    assert n > 0, f"Expected positive total reaction_count; got {n}"


def test_gene_reaction_count_consistent_with_edges(run_query):
    """Property matches actual outgoing edge count for sampled catalyzer genes."""
    rows = run_query("""
        MATCH (g:Gene) WHERE g.reaction_count > 0
        WITH g LIMIT 10
        OPTIONAL MATCH (g)-[r:Gene_catalyzes_reaction]->()
        WITH g, g.reaction_count AS prop, count(r) AS actual
        RETURN g.locus_tag AS lt, prop, actual
    """)
    assert len(rows) > 0, "No catalyzer genes found — KG may be empty"
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['lt']}: reaction_count={r['prop']} vs actual edges={r['actual']}"
        )


# ── KG-A2: Gene.catalyzed_metabolite_count (renamed from metabolite_count, KG-SYNC-001) ──

def test_retired_catalysis_arm_names_are_absent(run_query):
    """KG-SYNC-001: the bare names were RETIRED with the arm split, not aliased.

    Label-scoped on purpose: metabolite_count survives on Experiment /
    Publication / KeggTerm / TcdbFamily etc., and gene_count survives on every
    ontology node — only the three renamed homes must be clean.
    """
    for label, prop in [
        ("Gene", "metabolite_count"),
        ("OrganismTaxon", "metabolite_count"),
        ("Metabolite", "gene_count"),
    ]:
        n = run_query(
            f"MATCH (x:{label}) WHERE x.{prop} IS NOT NULL RETURN count(x) AS n"
        )[0]["n"]
        assert n == 0, f"{n} {label} nodes still carry retired property {prop}"

def test_gene_catalyzed_metabolite_count_present_on_all_genes(run_query):
    rows = run_query(
        "MATCH (g:Gene) RETURN count(g) AS total, count(g.catalyzed_metabolite_count) AS with_prop"
    )
    assert rows[0]["with_prop"] == rows[0]["total"]


def test_gene_catalyzed_metabolite_count_total_positive(run_query):
    n = run_query("MATCH (g:Gene) RETURN sum(g.catalyzed_metabolite_count) AS total")[0]["total"]
    assert n > 0


def test_gene_catalyzed_metabolite_count_consistent_with_2hop(run_query):
    """`catalyzed_metabolite_count` matches the CATALYSIS 2-hop exactly.

    BREAKING (TCDB rollup fix): the old `metabolite_count` was the UNION of
    catalysis and transport. The arms have incomparable epistemics — catalysis
    p90 = 11 metabolites, transport p90 = 554, because step 6 rolls every
    descendant's substrates onto each ancestor — so they are now separate
    properties, and the catalysis arm was renamed to name its arm (KG-SYNC-001).
    The transport arm is asserted in `test_gene_transported_metabolite_count_*`.
    """
    rows = run_query("""
        MATCH (g:Gene) WHERE g.catalyzed_metabolite_count > 0
        WITH g LIMIT 10
        OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)
                       -[:Reaction_has_metabolite]->(m_cat:Metabolite)
        RETURN g.locus_tag AS lt,
               g.catalyzed_metabolite_count AS prop,
               count(DISTINCT m_cat) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['lt']}: catalyzed_metabolite_count={r['prop']} vs catalysis 2-hop DISTINCT={r['actual']}"
        )


def test_gene_transported_metabolite_count_consistent_with_deepest_2hop(run_query):
    """The transport arm: 2-hop from the gene's DEEPEST TC attachments only.

    Attaching at every level would let a gene annotated at both 3.A.1 and
    3.A.1.14 inherit the ABC superfamily's full 554-substrate rollup despite
    having a more specific call. 6,950 genes are in that position.
    """
    rows = run_query("""
        MATCH (g:Gene) WHERE g.transported_metabolite_count > 0
        WITH g LIMIT 10
        OPTIONAL MATCH (g)-[:Gene_has_tcdb_family]->(t:TcdbFamily)
        WHERE NOT EXISTS {
            MATCH (g)-[:Gene_has_tcdb_family]->(d:TcdbFamily)
            WHERE (d)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t)
        }
        OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m_tr:Metabolite)
        RETURN g.locus_tag AS lt,
               g.transported_metabolite_count AS prop,
               count(DISTINCT m_tr) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['lt']}: transported_metabolite_count={r['prop']} vs deepest 2-hop DISTINCT={r['actual']}"
        )


# ── KG-A3: Metabolite.elements ────────────────────────────────────────────

def test_metabolite_elements_populated_for_metabolites_with_formula(run_query):
    """Every metabolite whose formula contains at least one element symbol
    (any A-Z character) must have an elements list. Formulas like '*2'
    (charge-only generic placeholders, e.g. KEGG C01322 'RX') parse to no
    elements and are correctly dropped by the sparse-output convention."""
    rows = run_query("""
        MATCH (m:Metabolite)
        WHERE m.formula IS NOT NULL AND m.formula =~ '.*[A-Z].*'
        RETURN count(m) AS total, count(m.elements) AS with_prop
    """)
    assert rows[0]["with_prop"] == rows[0]["total"], (
        f"Metabolites with element-bearing formula but no elements: "
        f"{rows[0]['total'] - rows[0]['with_prop']}"
    )


def test_metabolite_elements_sparse_when_no_formula(run_query):
    """Metabolites without formula must NOT carry an elements property (sparse convention)."""
    n = run_query("""
        MATCH (m:Metabolite) WHERE m.formula IS NULL AND m.elements IS NOT NULL
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"Found {n} no-formula metabolites with stale elements property"


def test_metabolite_elements_h2o_correct(run_query):
    """Spot check: H2O must contain exactly H and O."""
    rows = run_query(
        "MATCH (m:Metabolite) WHERE m.formula = 'H2O' RETURN m.elements AS els LIMIT 1"
    )
    if rows:
        assert sorted(rows[0]["els"]) == ["H", "O"]


def test_metabolite_elements_no_substring_footgun_on_sodium(run_query):
    """A metabolite whose formula contains 'Na' but not standalone 'N' must NOT
    have 'N' in its elements list (the whole point of having a parsed list)."""
    n = run_query("""
        MATCH (m:Metabolite)
        WHERE m.formula CONTAINS 'Na'
          AND NOT m.formula =~ '.*N[^a].*'
          AND NOT m.formula ENDS WITH 'N'
          AND 'N' IN m.elements
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"{n} metabolites have false 'N' element from 'Na' substring"


def test_metabolite_elements_no_false_carbon_from_chlorine(run_query):
    """Metabolites whose formula contains 'Cl' but no real C must not have 'C' in elements."""
    n = run_query("""
        MATCH (m:Metabolite)
        WHERE m.formula CONTAINS 'Cl'
          AND NOT m.formula =~ '.*C[^l].*'
          AND NOT m.formula ENDS WITH 'C'
          AND 'C' IN m.elements
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0


# ── KG-A4: KeggTerm pathway-level rollups (sparse) ────────────────────────

def test_keggterm_pathway_rollups_present_on_pathways(run_query):
    """All pathway-level KeggTerms have both rollup props."""
    rows = run_query("""
        MATCH (k:KeggTerm) WHERE k.level_kind = 'pathway'
        RETURN count(k) AS total,
               count(k.reaction_count) AS with_rxn,
               count(k.metabolite_count) AS with_met
    """)
    assert rows[0]["with_rxn"] == rows[0]["total"], "all pathways need reaction_count"
    assert rows[0]["with_met"] == rows[0]["total"], "all pathways need metabolite_count"


def test_keggterm_non_pathway_rollups_sparse(run_query):
    """KOs and categories must NOT have the rollup props (sparse semantics)."""
    rows = run_query("""
        MATCH (k:KeggTerm) WHERE k.level_kind <> 'pathway'
        RETURN count(k.reaction_count) AS rxn_set,
               count(k.metabolite_count) AS met_set
    """)
    assert rows[0]["rxn_set"] == 0, (
        "Non-pathway KeggTerms must not carry reaction_count (sparse convention)"
    )
    assert rows[0]["met_set"] == 0


def test_keggterm_pathway_rollup_consistent_with_edges(run_query):
    """Spot-check: reaction_count matches actual incoming Reaction_in_kegg_pathway edges."""
    rows = run_query("""
        MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway' AND p.reaction_count > 0
        WITH p LIMIT 10
        OPTIONAL MATCH (r:Reaction)-[:Reaction_in_kegg_pathway]->(p)
        RETURN p.id AS id, p.reaction_count AS prop, count(r) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"], (
            f"{r['id']}: reaction_count={r['prop']} vs actual={r['actual']}"
        )


def test_keggterm_pathway_metabolite_count_consistent(run_query):
    rows = run_query("""
        MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway' AND p.metabolite_count > 0
        WITH p LIMIT 10
        OPTIONAL MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p)
        RETURN p.id AS id, p.metabolite_count AS prop, count(m) AS actual
    """)
    assert len(rows) > 0
    for r in rows:
        assert r["prop"] == r["actual"]
