"""Live checks for the KEGG modules slice (spec docs/superpowers/specs/2026-11-18-kegg-modules-design.md).

Counts are pinned to the build-time dry run of 2026-11-18 (±5%, a KEGG refresh moves them a
little); invariants and biology hold on any healthy build.
"""
import re

import pytest

pytestmark = pytest.mark.kg

DRY_RUN = {"modules": 414, "classes": 57, "ko_edges": 1722, "pathway_edges": 1254,
           "organism_edges": 26004}
CYANO = ("Prochlorococcus", "Synechococcus", "Parasynechococcus", "Cyanobium", "Thermosynechococcus")
# One-step photosystem complexes these strains miss a required subunit KO for (annotation gaps,
# measured 2026-11-18): psbD K02706 in PCC 7942 and CB0205 (+ psbC K02707), psbB K02705 + psbD
# in the HNLC1 scaffold draft, PsaE K02693 in HNLC1.
PHOTOSYSTEM_GAPS = {"M00161": {"Synechococcus elongatus PCC 7942", "Synechococcus CB0205",
                               "Prochlorococcus HNLC1"},
                    "M00163": {"Prochlorococcus HNLC1"}}
# spec §2 table: steps_present / steps_total; None = no edge
BIOLOGY = {
    "Prochlorococcus MED4": {"M00009": (5, 8), "M00531": None, "M00002": (5, 5)},
    "Prochlorococcus MIT9313": {"M00009": (6, 8), "M00531": (1, 2), "M00002": (5, 5)},
    "Prochlorococcus NATL2A": {"M00009": (5, 8), "M00531": (1, 2), "M00002": (5, 5)},
    "Synechococcus WH8102": {"M00009": (6, 8), "M00531": (2, 2), "M00002": (5, 5)},
    "Alteromonas macleodii HOT1A3": {"M00009": (8, 8), "M00531": (1, 2), "M00002": (5, 5)},
    "Pseudomonas putida KT2440": {"M00009": (8, 8), "M00531": (1, 2), "M00002": (5, 5)},
}


def _near(actual, expected):
    assert abs(actual - expected) <= 0.05 * expected, (actual, expected)


def test_counts_near_the_dry_run(run_query):
    r = run_query("""
        RETURN COUNT { (m:KeggModule WHERE m.level = 3) } AS modules,
               COUNT { (m:KeggModule WHERE m.level < 3) } AS classes,
               COUNT { ()-[:Kegg_term_in_kegg_module]->() } AS ko_edges,
               COUNT { ()-[:Kegg_module_in_kegg_pathway]->() } AS pathway_edges,
               COUNT { ()-[:Organism_has_kegg_module]->() } AS organism_edges""")[0]
    for k, v in DRY_RUN.items():
        _near(r[k], v)


def test_organism_edge_invariants(run_query):
    bad = run_query("""
        MATCH (:OrganismTaxon)-[r:Organism_has_kegg_module]->(m:KeggModule)
        WITH r, m, r.steps_total - r.steps_present AS missing
        WHERE size(r.missing_step_indices) <> missing
           OR size(r.missing_step_definitions) <> missing
           OR coalesce(size(r.ec_genes), 0) <> r.steps_present_via_ec
           OR coalesce(size(r.ec_filled_step_indices), 0) <> r.steps_present_via_ec
           OR r.missing_step_indices IS NULL OR r.missing_step_definitions IS NULL
           OR r.missing_kos IS NULL
           OR r.steps_total <> m.step_count OR m.level <> 3 OR r.steps_present < 1
           OR r.status <> CASE missing WHEN 0 THEN 'complete' WHEN 1 THEN 'one_step_missing'
                                       WHEN 2 THEN 'two_steps_missing' ELSE 'partial' END
        RETURN m.module_id AS m LIMIT 5""")
    assert bad == []


def test_every_node_has_module_type_and_modules_one_parent(run_query):
    assert run_query("MATCH (m:KeggModule) WHERE m.module_type IS NULL RETURN m.id LIMIT 5") == []
    bad = run_query("""
        MATCH (m:KeggModule {level: 3})
        WITH m, COUNT { (m)-[:Kegg_module_is_a_kegg_module]->() } AS parents
        WHERE parents <> 1 RETURN m.module_id AS m, parents""")
    assert bad == []


def test_rollups(run_query):
    bad = run_query("""
        MATCH (m:KeggModule)
        WHERE m.gene_count IS NULL OR m.organism_count IS NULL OR m.complete_organism_count IS NULL
           OR (m.level = 3) <> (m.direct_gene_count IS NOT NULL)
        RETURN m.id LIMIT 5""")
    assert bad == []
    r = run_query("""MATCH (m:KeggModule {id: 'kegg.module:pathway_modules'})
                     RETURN m.complete_organism_count AS c, m.gene_count AS g""")[0]
    assert r["c"] >= 100 and r["g"] > 10000


def test_photosystems_complete_in_cyanobacteria_only(run_query):
    rows = run_query("""
        MATCH (o:OrganismTaxon {organism_type: 'genome_strain'})
        UNWIND ['M00161', 'M00163'] AS mid
        OPTIONAL MATCH (o)-[r:Organism_has_kegg_module]->(:KeggModule {module_id: mid})
        RETURN o.preferred_name AS org, mid, r.status AS status""")
    for row in rows:
        cyano = row["org"].startswith(CYANO)
        if cyano and row["org"] not in PHOTOSYSTEM_GAPS[row["mid"]]:
            assert row["status"] == "complete", row
        else:
            assert row["status"] != "complete", row


def test_spec_biology_table(run_query):
    rows = run_query("""
        MATCH (o:OrganismTaxon) WHERE o.preferred_name IN $orgs
        UNWIND ['M00009', 'M00531', 'M00002'] AS mid
        OPTIONAL MATCH (o)-[r:Organism_has_kegg_module]->(:KeggModule {module_id: mid})
        RETURN o.preferred_name AS org, mid, r.steps_present AS p, r.steps_total AS t""",
                     orgs=list(BIOLOGY))
    assert len(rows) == 3 * len(BIOLOGY)
    for row in rows:
        got = None if row["p"] is None else (row["p"], row["t"])
        assert got == BIOLOGY[row["org"]][row["mid"]], row


def test_natl2a_pyruvate_kinase_through_ec(run_query):
    r = run_query("""
        MATCH (:OrganismTaxon {preferred_name: 'Prochlorococcus NATL2A'})
              -[r:Organism_has_kegg_module]->(:KeggModule {module_id: 'M00002'})
        RETURN r.ec_genes AS g, r.ec_filled_step_indices AS i, r.status AS s""")[0]
    assert (r["g"], r["i"], r["s"]) == (["PMN2A_0295"], [5], "complete")


def test_one_source_of_kos(run_query):
    """The evaluator's KO-only steps equal the steps Gene_has_kegg_ko reaches for modules whose
    steps are plain KO alternatives, and never exceed them for KO-only modules otherwise."""
    defs = run_query("MATCH (m:KeggModule {level: 3}) RETURN m.module_id AS m, m.definition AS d")
    plain_step = re.compile(r"^\(?K\d{5}(,K\d{5})*\)?$")
    plain = [r["m"] for r in defs
             if all(plain_step.match(s) for s in r["d"].replace("\n", " ").split() if s != "--")]
    no_refs = [r["m"] for r in defs if not re.search(r"M\d{5}", r["d"])]
    rows = run_query("""
        MATCH (o:OrganismTaxon)-[r:Organism_has_kegg_module]->(m:KeggModule)
        WHERE m.module_id IN $mods
        CALL {
          WITH o, m
          OPTIONAL MATCH (o)<-[:Gene_belongs_to_organism]-(:Gene)-[:Gene_has_kegg_ko]->(:KeggTerm)
                         -[k:Kegg_term_in_kegg_module]->(m)
          UNWIND CASE WHEN k IS NULL THEN [] ELSE k.steps END AS st
          RETURN count(DISTINCT st) AS reached
        }
        RETURN m.module_id AS m, o.preferred_name AS o,
               r.steps_present - r.steps_present_via_ec AS ko_only, reached""", mods=no_refs)
    plain_set = set(plain)
    assert len(plain_set) > 100 and rows
    for row in rows:
        if row["m"] in plain_set:
            assert row["ko_only"] == row["reached"], row
        else:
            assert row["ko_only"] <= row["reached"], row
