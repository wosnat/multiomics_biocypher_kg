"""KG validity tests for the metabolism layer (Reaction + Metabolite + 4 edges).

Requires: a running Neo4j instance with the deployed graph (skipped otherwise).
Mark: @pytest.mark.kg
"""
from __future__ import annotations

import pytest


pytestmark = [pytest.mark.kg]


def test_reaction_count_in_expected_range(run_query):
    """Spec acceptance: >= 2,000 reactions (gene-reachable subset of KEGG ~14K)."""
    n = run_query("MATCH (r:Reaction) RETURN count(r) AS n")[0]["n"]
    assert 2000 <= n <= 6000, f"Reaction count {n} outside expected 2000-6000"


def test_metabolite_count_in_expected_range(run_query):
    """Spec acceptance: 1,000 <= count <= 6,000 (post-TCDB transport-only growth)."""
    n = run_query("MATCH (m:Metabolite) RETURN count(m) AS n")[0]["n"]
    assert 1000 <= n <= 6000, f"Metabolite count {n} outside expected 1000-6000"


def test_all_reactions_have_kegg_primary_id(run_query):
    """Primary IDs all start with kegg.reaction: per Spec 1.2."""
    n_bad = run_query(
        "MATCH (r:Reaction) WHERE NOT r.id STARTS WITH 'kegg.reaction:' RETURN count(r) AS n"
    )[0]["n"]
    assert n_bad == 0, f"{n_bad} Reaction nodes lack kegg.reaction: prefix"


def test_metabolite_primary_ids_predominantly_kegg(run_query):
    """Most metabolites should have kegg.compound: prefix; a few may fall back to chebi/mnx."""
    total = run_query("MATCH (m:Metabolite) RETURN count(m) AS n")[0]["n"]
    assert total > 0, "No Metabolite nodes found — graph not yet rebuilt with metabolism layer"
    kegg_count = run_query(
        "MATCH (m:Metabolite) WHERE m.id STARTS WITH 'kegg.compound:' RETURN count(m) AS n"
    )[0]["n"]
    assert kegg_count / total >= 0.95, (
        f"Only {kegg_count}/{total} metabolites have kegg.compound prefix"
    )


def test_kegg_reaction_id_property_matches_id_suffix(run_query):
    """Reaction.kegg_reaction_id should match the suffix of the primary id."""
    n_bad = run_query(
        "MATCH (r:Reaction) "
        "WHERE 'kegg.reaction:' + r.kegg_reaction_id <> r.id "
        "RETURN count(r) AS n"
    )[0]["n"]
    assert n_bad == 0


def test_every_reaction_has_at_least_one_metabolite(run_query):
    """Every Reaction should connect to >=1 Metabolite (some KEGG-only orphans tolerated)."""
    total = run_query("MATCH (r:Reaction) RETURN count(r) AS n")[0]["n"]
    assert total > 0, "No Reaction nodes found — graph not yet rebuilt with metabolism layer"
    n_bad = run_query("""
        MATCH (r:Reaction)
        WHERE NOT EXISTS { MATCH (r)-[:Reaction_has_metabolite]->() }
        RETURN count(r) AS n
    """)[0]["n"]
    assert n_bad / total < 0.05, f"{n_bad}/{total} reactions lack any metabolite edge"


def test_organism_has_metabolite_edges_materialized(run_query):
    """>=1K Organism_has_metabolite edges should exist after post-import."""
    n = run_query(
        "MATCH ()-[r:Organism_has_metabolite]->() RETURN count(r) AS n"
    )[0]["n"]
    assert n >= 1000, f"Only {n} Organism_has_metabolite edges (expected >= 1000)"


def test_pyruvate_kinase_spot_check(run_query):
    """Spec spot-check: kegg.reaction:R00200 (pyruvate kinase) exists with expected props."""
    results = run_query(
        "MATCH (r:Reaction {id: 'kegg.reaction:R00200'}) RETURN r"
    )
    assert len(results) > 0, "kegg.reaction:R00200 not found"
    rxn = results[0]["r"]
    assert rxn["kegg_reaction_id"] == "R00200"
    n_meta = run_query(
        "MATCH (:Reaction {id: 'kegg.reaction:R00200'})-[:Reaction_has_metabolite]->(m) "
        "RETURN count(m) AS n"
    )[0]["n"]
    assert n_meta == 4, f"R00200 has {n_meta} metabolites; expected 4"
    has_glycolysis = run_query(
        "MATCH (:Reaction {id: 'kegg.reaction:R00200'})-[:Reaction_in_kegg_pathway]->(p) "
        "WHERE p.id CONTAINS 'ko00010' RETURN count(p) AS n"
    )[0]["n"]
    assert has_glycolysis >= 1


def test_reaction_class_values_valid(run_query):
    """reaction_class enum: 'transport' or 'chemical'."""
    bad = run_query(
        "MATCH (r:Reaction) WHERE NOT r.reaction_class IN ['transport', 'chemical'] "
        "RETURN count(r) AS n"
    )[0]["n"]
    assert bad == 0


def test_mass_balance_values_valid(run_query):
    """mass_balance enum: 'balanced' or 'unbalanced'."""
    bad = run_query(
        "MATCH (r:Reaction) WHERE NOT r.mass_balance IN ['balanced', 'unbalanced'] "
        "RETURN count(r) AS n"
    )[0]["n"]
    assert bad == 0


def test_metabolite_in_pathway_edges_present(run_query):
    """Spec 1.2.1: ~5K-15K Metabolite_in_pathway edges (Option B-pruned)."""
    n = run_query(
        "MATCH ()-[r:Metabolite_in_pathway]->() RETURN count(r) AS n"
    )[0]["n"]
    assert 5000 <= n <= 15000, f"{n} Metabolite_in_pathway edges (expected 5000-15000)"


def test_metabolite_pathways_only_kg_evidenced(run_query):
    """Option B: Metabolite_in_pathway targets must also be reachable from genes."""
    # No metabolite-in-pathway edge should target a pathway with no Gene→KO edge to it
    # (since Option B prunes compound-only pathways at step 6).
    n_orphan = run_query("""
        MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p:KeggTerm)
        WHERE NOT EXISTS {
            MATCH (g:Gene)-[:Gene_has_kegg_ko]->(:KeggTerm)-[:Kegg_term_is_a_kegg_term]->(p)
        }
        AND NOT EXISTS {
            MATCH (g:Gene)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_in_kegg_pathway]->(p)
        }
        RETURN count(DISTINCT p) AS n
    """)[0]["n"]
    assert n_orphan == 0, f"{n_orphan} pathways have only Metabolite-in-pathway edges (Option B violation)"
