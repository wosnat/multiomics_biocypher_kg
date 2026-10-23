"""KG validity tests for GO transport edges (GO term -> Metabolite / ChemicalRole).

Requires: a running Neo4j instance with a graph built WITH the go-transport
feature (skipped if Neo4j is unreachable; fails on a pre-feature graph).
Mark: @pytest.mark.kg
"""
from __future__ import annotations

import pytest

pytestmark = [pytest.mark.kg]

# Baselines captured on the GO-refreshed graph BEFORE this feature
# (scratchpad baselines/counts.txt, Task 0 Step 10).
BASELINE_PATHWAY_COUNT = 452
BASELINE_SUM_REACTION_COUNT = 6408
BASELINE_SUM_METABOLITE_COUNT = 9342
# Task 4 measured exactly 6 GO-only compound links into kept pathways
# (controller ruling: GO-only compounds may link into kept pathways).
GO_ONLY_PATHWAY_LINKS = 6
# TCDB cleanup (2026-10-23, spec 2026-10-22 §5.4): losing the TCDB class-level
# substrate roll-up removed 3 gene-less pathways (ko04145 Phagosome, ko04510
# Focal adhesion, ko00532 chondroitin/dermatan sulfate GAG biosynthesis — see
# docs/kg-changes/tcdb-cleanup/removed_ids.tsv) and 107 Metabolite_in_pathway
# links (diff of cache/data/kegg/kegg_data.json main vs feat/tcdb-cleanup:
# 9,348 -> 9,241 compound->pathway pairs, 0 added). Reactions are untouched.
TCDB_CLEANUP_REMOVED_PATHWAYS = 3
TCDB_CLEANUP_LOST_PATHWAY_LINKS = 107

METABOLITE_EDGES = [
    "Molecular_function_transports_metabolite",
    "Biological_process_transports_metabolite",
]
ROLE_EDGES = [
    "Molecular_function_transports_chemical_role",
    "Biological_process_transports_chemical_role",
]
ALL_EDGES = METABOLITE_EDGES + ROLE_EDGES


@pytest.mark.parametrize("rel", ALL_EDGES)
def test_edge_type_present(run_query, rel):
    n = run_query(f"MATCH ()-[r:{rel}]->() RETURN count(r) AS n")[0]["n"]
    assert n > 0, f"no {rel} edges"


@pytest.mark.parametrize(
    "rel,target_label",
    [(r, "Metabolite") for r in METABOLITE_EDGES]
    + [(r, "ChemicalRole") for r in ROLE_EDGES],
)
def test_edge_endpoints_valid(run_query, rel, target_label):
    n = run_query(
        f"""
        MATCH (s)-[r:{rel}]->(t)
        WHERE NOT (s:MolecularFunction OR s:BiologicalProcess)
           OR NOT t:{target_label}
        RETURN count(r) AS n
        """
    )[0]["n"]
    assert n == 0, f"{n} {rel} edges with a non-GO source or non-{target_label} target"


@pytest.mark.parametrize("rel", ALL_EDGES)
def test_edges_carry_go_chebi_ids(run_query, rel):
    n = run_query(
        f"""
        MATCH ()-[r:{rel}]->()
        WHERE r.go_chebi_ids IS NULL OR size(r.go_chebi_ids) = 0
        RETURN count(r) AS n
        """
    )[0]["n"]
    assert n == 0, f"{n} {rel} edges lack go_chebi_ids"


def test_carbonate_like_transport_edge(run_query):
    rows = run_query(
        """
        MATCH (:MolecularFunction {id: 'go:0015110'})
              -[r:Molecular_function_transports_metabolite]->
              (m:Metabolite {id: 'kegg.compound:C01417'})
        RETURN r.go_chebi_ids AS ids
        """
    )
    assert rows, "go:0015110 -> kegg.compound:C01417 edge missing"
    assert "CHEBI:29195" in rows[0]["ids"]


def test_siderophore_transport_role_edge(run_query):
    n = run_query(
        """
        MATCH (:BiologicalProcess {id: 'go:0015891'})
              -[r:Biological_process_transports_chemical_role]->
              (:ChemicalRole {id: 'chebi.role:26672'})
        RETURN count(r) AS n
        """
    )[0]["n"]
    assert n == 1, "go:0015891 -> chebi.role:26672 edge missing"


def test_excluded_targets_recorded(run_query):
    rows = run_query(
        "MATCH (s:Schema_info {id: 'schema_info'}) RETURN s.go_transport_excluded_targets AS x"
    )
    assert rows and rows[0]["x"], "Schema_info.go_transport_excluded_targets missing"
    excluded = rows[0]["x"]
    assert len(excluded) >= 13, f"only {len(excluded)} excluded targets (13 measured)"
    assert "kegg.compound:C00080" in excluded, "H+ not in exclusion list"
    assert "kegg.compound:C01330" in excluded, "Na+ not in exclusion list"


def test_no_edge_targets_excluded_id(run_query):
    rows = run_query(
        "MATCH (s:Schema_info {id: 'schema_info'}) RETURN s.go_transport_excluded_targets AS x"
    )
    assert rows and rows[0]["x"] and len(rows[0]["x"]) >= 13, "exclusion list missing/short"
    n = run_query(
        """
        MATCH (s:Schema_info {id: 'schema_info'})
        WITH s.go_transport_excluded_targets AS ex
        MATCH ()-[r]->(t)
        WHERE type(r) IN $rels AND t.id IN ex
        RETURN count(r) AS n
        """,
        rels=ALL_EDGES,
    )[0]["n"]
    assert n == 0, f"{n} GO transport edges target an excluded id"


def test_fadD_not_a_transporter_via_go(run_query):
    has_go = run_query(
        """
        MATCH (:Gene {locus_tag: 'PMM0402'})
              -[:Gene_enables_molecular_function|Gene_involved_in_biological_process]->(t)
        RETURN count(t) AS n
        """
    )[0]["n"]
    assert has_go >= 1, "PMM0402 missing or has no GO edges (test would be vacuous)"
    n = run_query(
        """
        MATCH (g:Gene {locus_tag: 'PMM0402'})
              -[:Gene_enables_molecular_function|Gene_involved_in_biological_process]->
              (:MolecularFunction|BiologicalProcess)
              -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->
              (m:Metabolite)
        RETURN count(DISTINCT m) AS n
        """
    )[0]["n"]
    assert n == 0, f"fadD (PMM0402) reaches {n} metabolites via GO transport"


def test_med4_genes_reaching_metabolites_in_range(run_query):
    n = run_query(
        """
        MATCH (g:Gene {organism_name: 'Prochlorococcus MED4'})
              -[:Gene_enables_molecular_function|Gene_involved_in_biological_process]->
              (:MolecularFunction|BiologicalProcess)
              -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->
              (:Metabolite)
        RETURN count(DISTINCT g) AS n
        """
    )[0]["n"]
    assert 45 <= n <= 75, f"MED4 genes reaching a metabolite via GO route: {n} (measured 59)"


def test_go_transport_organism_props_consistent(run_query):
    bad = run_query(
        """
        MATCH (m:Metabolite)
        WHERE m.go_transport_organism_names IS NOT NULL
          AND m.go_transport_organism_count <> size(m.go_transport_organism_names)
        RETURN count(m) AS n
        """
    )[0]["n"]
    assert bad == 0, f"{bad} metabolites with count != size(names)"
    carriers = run_query(
        "MATCH (m:Metabolite) WHERE m.go_transport_organism_names IS NOT NULL RETURN count(m) AS n"
    )[0]["n"]
    assert carriers > 0, "no metabolite carries go_transport_organism_names"
    # sparse: absent elsewhere, and the two properties co-occur
    mismatch = run_query(
        """
        MATCH (m:Metabolite)
        WHERE (m.go_transport_organism_names IS NULL) <> (m.go_transport_organism_count IS NULL)
        RETURN count(m) AS n
        """
    )[0]["n"]
    assert mismatch == 0, f"{mismatch} metabolites carry only one of names/count"
    # carriers are exactly the metabolites reached by a GO transport edge
    off = run_query(
        """
        MATCH (m:Metabolite)
        WHERE m.go_transport_organism_names IS NOT NULL
          AND NOT EXISTS {
            MATCH (:MolecularFunction|BiologicalProcess)
                  -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m)
          }
        RETURN count(m) AS n
        """
    )[0]["n"]
    assert off == 0, f"{off} carriers with no incoming GO transport edge"
    stray = run_query(
        """
        MATCH (m:Metabolite)
        WHERE NOT EXISTS {
            MATCH (:MolecularFunction|BiologicalProcess)
                  -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m)
          }
          AND (m.go_transport_organism_count IS NOT NULL
               OR m.go_transport_organism_names IS NOT NULL)
        RETURN count(m) AS n
        """
    )[0]["n"]
    assert stray == 0, f"{stray} metabolites without GO transport edges carry go_transport props"


def test_go_transport_evidence_source_has_edge(run_query):
    bad = run_query(
        """
        MATCH (m:Metabolite)
        WHERE 'go_transport' IN m.evidence_sources
          AND NOT EXISTS {
            MATCH (:MolecularFunction|BiologicalProcess)
                  -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m)
          }
        RETURN count(m) AS n
        """
    )[0]["n"]
    assert bad == 0, f"{bad} go_transport metabolites with no incoming GO transport edge"


def test_pathway_and_reaction_counts_unchanged(run_query):
    n = run_query(
        "MATCH (k:KeggTerm {level_kind: 'pathway'}) RETURN count(k) AS n"
    )[0]["n"]
    expected_n = BASELINE_PATHWAY_COUNT - TCDB_CLEANUP_REMOVED_PATHWAYS   # 449
    assert n == expected_n, f"pathway count {n} != {expected_n}"
    s = run_query(
        "MATCH (k:KeggTerm {level_kind: 'pathway'}) "
        "RETURN sum(k.reaction_count) AS r, sum(k.metabolite_count) AS m"
    )[0]
    assert s["r"] == BASELINE_SUM_REACTION_COUNT, f"sum(reaction_count) {s['r']}"
    expected_m = (BASELINE_SUM_METABOLITE_COUNT + GO_ONLY_PATHWAY_LINKS
                  - TCDB_CLEANUP_LOST_PATHWAY_LINKS)   # 9,241
    assert s["m"] == expected_m, (
        f"sum(metabolite_count) {s['m']} != {BASELINE_SUM_METABOLITE_COUNT} + "
        f"{GO_ONLY_PATHWAY_LINKS} - {TCDB_CLEANUP_LOST_PATHWAY_LINKS}"
    )
