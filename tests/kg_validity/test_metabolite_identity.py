"""Live-graph invariants for the metabolite identity + classes change.

Spec: docs/superpowers/specs/2026-10-13-metabolite-identity-and-classes-design.md §11
(KG validity). Every test is read-only. They describe the POST-rebuild graph; run
against a pre-change graph most of them fail (that is the RED evidence).
"""
import csv
import json
import pathlib

import pytest

pytestmark = pytest.mark.kg

ROOT = pathlib.Path(__file__).resolve().parents[2]
REMAP = ROOT / "docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv"
PRE_ELEMENTS = ROOT / "tests/kg_validity/fixtures/metabolite_elements_pre_identity.json"
PRE_COUNTS = ROOT / "tests/kg_validity/fixtures/metabolite_identity_pre_counts.json"
TCDB_PRUNED = ROOT / "cache/data/tcdb/tcdb_pruned.json"

# all_identifiers CURIE forms (spec §3.2). The mnx local part may be non-MNXM (e.g. mnx:WATER).
CURIE_RE = (r"kegg\.compound:C\d{5}|chebi:\d+|mnx:\S+|hmdb:HMDB\d{7}"
            r"|inchikey:[A-Z]{14}-[A-Z]{10}-[A-Z]")
# Placeholder / generic formula: '*', an R group (R not starting Ra/Rb/Re/Rf/Rg/Rh/Ru), or a ")n" repeat.
PLACEHOLDER_FORMULA_RE = r".*(\*|R(?![abefghu])|\)n).*"

IDENTITY_BASES = {"mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual"}
# remap.tsv bases that name a merge / re-key target (the new_id carries identity_basis)
MERGE_BASES = {"mnx_same_entity", "chebi_conjugate", "chebi_tautomer", "manual"}
REKEY_BASIS = "rekey_chebi_kegg_xref"


def _need(path: pathlib.Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(ROOT)} missing; rebuild prepare_data step 6 / restore the fixture")
    return path


def _one(run_query, q, **kw):
    return run_query(q, **kw)[0]


def _remap_rows():
    with _need(REMAP).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- identity: one node per compound ------------------------------------------

def test_no_duplicate_identifiers(run_query):
    for prop in ("kegg_compound_id", "chebi_id", "mnxm_id"):
        dup = run_query(f"MATCH (m:Metabolite) WHERE m.{prop} IS NOT NULL WITH m.{prop} AS v, count(*) AS n "
                        f"WHERE n > 1 RETURN v, n LIMIT 5")
        assert not dup, (prop, dup)
    dup = run_query("""
        CALL { MATCH (m:Metabolite) UNWIND m.all_identifiers AS x RETURN x
               UNION ALL MATCH (r:ChemicalRole) UNWIND r.all_identifiers AS x RETURN x }
        WITH x, count(*) AS n WHERE n > 1 RETURN x, n LIMIT 5""")
    assert not dup, dup


def test_no_inchikey_on_two_nodes(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE m.inchikey IS NOT NULL WITH m.inchikey AS k, count(*) AS n "
                         "WHERE n > 1 RETURN k LIMIT 5")


def test_chebi_keyed_nodes_chebi_id_matches(run_query):
    bad = run_query("MATCH (m:Metabolite) WHERE m.id STARTS WITH 'chebi:' "
                    "AND coalesce(m.chebi_id, '') <> split(m.id, ':')[1] RETURN m.id, m.chebi_id LIMIT 5")
    assert not bad, bad


def test_all_identifiers_are_curies(run_query):
    missing = run_query("MATCH (m:Metabolite) WHERE m.all_identifiers IS NULL OR size(m.all_identifiers) = 0 "
                        "RETURN m.id LIMIT 5")
    assert not missing, missing
    bad = run_query("MATCH (m:Metabolite) UNWIND m.all_identifiers AS x WITH m, x "
                    "WHERE NOT x =~ $re RETURN m.id, x LIMIT 5", re=f"(?:{CURIE_RE})")
    assert not bad, bad
    # every node's own id is one of its identifiers
    assert not run_query("MATCH (m:Metabolite) WHERE NOT m.id IN m.all_identifiers RETURN m.id LIMIT 5")


# --- ChemicalRole --------------------------------------------------------------

def test_no_role_metabolites_and_role_ids(run_query):
    role_ids = [r["id"] for r in run_query("MATCH (r:ChemicalRole) RETURN r.id AS id")]
    assert role_ids, "no ChemicalRole nodes"
    assert all(i.startswith("chebi.role:") and i.split(":")[1].isdigit() for i in role_ids)
    # the role's own ChEBI term is one of its identifiers
    assert not run_query("MATCH (r:ChemicalRole) WHERE r.all_identifiers IS NULL "
                         "OR NOT 'chebi:' + split(r.id, ':')[1] IN r.all_identifiers RETURN r.id LIMIT 5")
    secondary_hits = run_query("""MATCH (r:ChemicalRole) UNWIND r.all_identifiers AS x
        WITH x WHERE x STARTS WITH 'chebi:' MATCH (m:Metabolite {id: x}) RETURN m.id LIMIT 5""")
    assert not secondary_hits, secondary_hits
    # every remap row with basis removed_role is gone as a Metabolite
    removed = [r["old_id"] for r in _remap_rows() if r["basis"] == "removed_role"]
    assert removed
    assert not run_query("MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id LIMIT 5", ids=removed)
    assert not run_query("MATCH (r:ChemicalRole) WHERE r.level_is_best_effort IS NOT NULL RETURN r.id LIMIT 1")


def test_old_secondary_role_id_resolves(run_query):
    rows = run_query("MATCH (r:ChemicalRole) WHERE 'chebi:10004' IN r.all_identifiers RETURN r.id AS id")
    assert [r["id"] for r in rows] == ["chebi.role:33229"]


def test_role_subtree_counts_distinct_and_ge_direct(run_query):
    bad = run_query("MATCH (r:ChemicalRole) WHERE r.tcdb_family_count IS NULL OR r.transporter_count IS NULL "
                    "OR r.tcdb_family_count < r.direct_tcdb_family_count "
                    "OR r.transporter_count < r.direct_transporter_count RETURN r.id LIMIT 5")
    assert not bad, bad
    row = _one(run_query, """MATCH (r:ChemicalRole) WHERE r.level = 0
        OPTIONAL MATCH (r)<-[:Chemical_role_is_a_chemical_role*0..]-(d:ChemicalRole)
                       <-[:Tcdb_family_transports_chemical_role]-(t:TcdbFamily)
        WITH r, count(DISTINCT t) AS n RETURN collect(coalesce(r.tcdb_family_count = n, false)) AS ok""")
    assert row["ok"], "no level-0 ChemicalRole"
    assert all(row["ok"])
    # transporter_count = DISTINCT families via most_specific role edges over the subtree
    row = _one(run_query, """MATCH (r:ChemicalRole) WHERE r.level = 0
        OPTIONAL MATCH (r)<-[:Chemical_role_is_a_chemical_role*0..]-(d:ChemicalRole)
                       <-[e:Tcdb_family_transports_chemical_role]-(t:TcdbFamily)
        WHERE e.substrate_depth = 'most_specific'
        WITH r, count(DISTINCT t) AS n RETURN collect(coalesce(r.transporter_count = n, false)) AS ok""")
    assert row["ok"] and all(row["ok"])


def test_roles_not_propagated(run_query):
    """Edge-level proxy for "roles are never propagated": no ChemicalRole node is adjacent to
    a Gene, OrganismTaxon or Metabolite. Role facts reach those only by traversal."""
    assert not run_query("MATCH (r:ChemicalRole)--(x) WHERE x:Gene OR x:OrganismTaxon OR x:Metabolite "
                         "RETURN x LIMIT 1")


# --- Metabolite set ---------------------------------------------------------------

def test_no_isolated_metabolite_and_metabolomics_tag_matches_assays(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE NOT EXISTS { (m)--() } RETURN m.id LIMIT 5")
    assert not run_query("""MATCH (m:Metabolite) WHERE 'metabolomics' IN m.evidence_sources
        AND NOT EXISTS { (:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m) }
        RETURN m.id LIMIT 5""")
    assert not run_query("""MATCH (:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
        WHERE NOT 'metabolomics' IN coalesce(m.evidence_sources, []) RETURN DISTINCT m.id LIMIT 5""")


def test_recurated_ids_are_not_identifiers(run_query):
    assert run_query("MATCH (m:Metabolite) WHERE m.recurated_from_ids IS NOT NULL RETURN m.id LIMIT 1"), \
        "no recurated_from_ids anywhere"
    assert not run_query("MATCH (m:Metabolite) UNWIND coalesce(m.recurated_from_ids, []) AS x "
                         "MATCH (n:Metabolite) WHERE x IN n.all_identifiers OR n.id = x RETURN x LIMIT 5")
    dup = run_query("MATCH (m:Metabolite) UNWIND coalesce(m.recurated_from_ids, []) AS x WITH x, count(*) AS n "
                    "WHERE n > 1 RETURN x")
    assert not dup, dup


def test_structure_and_elements(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind IS NULL RETURN m.id LIMIT 5")
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind <> 'defined' AND m.inchikey IS NOT NULL "
                         "RETURN m.id LIMIT 5")
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind <> 'defined' AND m.mass IS NOT NULL "
                         "RETURN m.id LIMIT 5")
    # direct check, independent of structure_kind: a placeholder formula never carries a mass
    bad = run_query("MATCH (m:Metabolite) WHERE m.formula =~ $re AND m.mass IS NOT NULL "
                    "RETURN m.id, m.formula LIMIT 5", re=PLACEHOLDER_FORMULA_RE)
    assert not bad, bad


def test_no_node_loses_elements(run_query):
    """Every pre-change metabolite with `elements` still has them on its (possibly remapped) node.

    Skipped: remap rows with an EMPTY new_id — removed_role (now a ChemicalRole, not a
    Metabolite) and no_longer_reachable (no node at all). Every other basis is checked,
    including class_term_rejected (its replacement chebi node is generic/defined but
    still carries elements — verified against cache/data/kegg/kegg_data.json on
    2026-10-14) and removed_placeholder (remapped onto the real KEGG compound). The
    remap TSV has no isotopologue_rejected rows at the time of writing.
    """
    pre = json.loads(_need(PRE_ELEMENTS).read_text())
    remap = {r["old_id"]: r["new_id"] for r in _remap_rows()}
    with_el = {r["id"] for r in run_query(
        "MATCH (m:Metabolite) WHERE m.elements IS NOT NULL AND size(m.elements) > 0 RETURN m.id AS id")}
    lost, skipped = [], []
    for old in pre:
        new = remap.get(old, old)
        if not new:          # removed role / no longer reachable: not a metabolite any more
            skipped.append(old)
            continue
        if new not in with_el:
            lost.append((old, new))
    assert not lost, lost[:10]
    # ...and the skipped ids really are gone from the Metabolite id set
    still = run_query("MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id AS id", ids=skipped)
    assert not still, still[:5]


# Known, explained divergences between the remap TSV and identity_basis (as built 2026-10-14,
# checked against cache/data/kegg/kegg_data.json). The remap TSV covers only ids that were on
# the PRE-change graph, and its fallback rows (outputs._remap_fallback) label a relocation
# "rekey_chebi_kegg_xref" / "mnx_same_entity" approximately — the old term is now an identifier
# of that node, not necessarily a re-key or merge. Pinned exactly, so any drift fails.
#  - remap target, no identity_basis: chebi:176600 -> chebi:49470 is a fallback "mnx_same_entity"
#    row (the old substrate string now resolves there by MNX id); nothing was merged.
REMAP_TARGETS_WITHOUT_BASIS = {"chebi:49470"}
#  - rekey targets whose identity_basis lacks chebi_kegg_xref: the old chebi id is a member of the
#    node's MNX / conjugate group (or a fallback relocation), so the basis names that instead.
REKEY_TARGETS_WITHOUT_XREF_BASIS = {
    "kegg.compound:C00261", "kegg.compound:C01197", "kegg.compound:C05345", "kegg.compound:C05984"}
#  - identity_basis nodes that are no remap merge/rekey target: everything merged or re-keyed into
#    them was never a node of the pre-change graph (or, C00197, the pre-change id is remapped as
#    remeasured_elsewhere).
BASIS_NODES_NOT_IN_REMAP = {
    "chebi:104011", "kegg.compound:C00076", "kegg.compound:C00117", "kegg.compound:C00197",
    "kegg.compound:C00760", "kegg.compound:C01413", "kegg.compound:C01935", "kegg.compound:C06269",
    "kegg.compound:C10309"}


def test_identity_basis_only_on_merged_or_rekeyed(run_query):
    rows = _remap_rows()
    basis = {r["id"]: r["b"] for r in run_query(
        "MATCH (m:Metabolite) WHERE m.identity_basis IS NOT NULL RETURN m.id AS id, m.identity_basis AS b")}
    assert basis, "no identity_basis anywhere"
    values = {v for b in basis.values() for v in b}
    assert values <= IDENTITY_BASES, values - IDENTITY_BASES
    assert not [i for i, b in basis.items() if not b], "empty identity_basis list"

    rekey_targets = {r["new_id"] for r in rows if r["basis"] == REKEY_BASIS and r["new_id"]}
    assert rekey_targets
    no_xref = {t for t in rekey_targets if "chebi_kegg_xref" not in basis.get(t, [])}
    assert no_xref == REKEY_TARGETS_WITHOUT_XREF_BASIS, sorted(no_xref ^ REKEY_TARGETS_WITHOUT_XREF_BASIS)

    expected = rekey_targets | {r["new_id"] for r in rows if r["basis"] in MERGE_BASES and r["new_id"]}
    assert expected - set(basis) == REMAP_TARGETS_WITHOUT_BASIS, sorted(expected - set(basis))
    assert set(basis) - expected == BASIS_NODES_NOT_IN_REMAP, sorted(set(basis) - expected)


def test_every_rekeyed_node_has_kegg_enrichment(run_query):
    targets = sorted({r["new_id"] for r in _remap_rows() if r["basis"] == REKEY_BASIS and r["new_id"]})
    assert targets
    live = {r["id"]: r for r in run_query(
        "MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id AS id, m.name AS name, "
        "m.evidence_sources AS ev, m.structure_kind AS sk, m.formula AS formula", ids=targets)}
    missing = [t for t in targets if t not in live]
    assert not missing, missing[:10]
    bad = [r for r in live.values()
           if not r["name"] or not r["ev"] or (r["sk"] == "defined" and not r["formula"])]
    assert not bad, bad[:10]


# --- stereo / enantiomer edges -------------------------------------------------------

def test_stereo_edges(run_query):
    n = _one(run_query, "MATCH ()-[e:Metabolite_is_stereo_form_of|Metabolite_is_enantiomer_of]->() "
                        "RETURN count(e) AS n")["n"]
    assert n > 0, "no stereo/enantiomer edges"
    assert not run_query("MATCH (a)-[e:Metabolite_is_stereo_form_of|Metabolite_is_enantiomer_of]->(b) "
                         "WHERE a = b OR e.sources IS NULL OR e.evidence IS NULL "
                         "OR e.sources <> ['chebi'] OR e.evidence <> 'curated' RETURN a.id LIMIT 5")
    # symmetric relation stored once, lower -> higher id
    assert not run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b) WHERE a.id > b.id RETURN a.id LIMIT 5")
    assert not run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b)-[:Metabolite_is_enantiomer_of]->(a) "
                         "RETURN a.id LIMIT 1")


def test_no_anomer_enantiomer_edges(run_query):
    from multiomics_kg.utils.metabolite_identity.names import strip_anomer_tokens
    names = run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b) RETURN a.name AS a, b.name AS b")
    assert names, "no enantiomer edges"
    bad = [r for r in names if strip_anomer_tokens(r["a"]) == strip_anomer_tokens(r["b"])]
    assert not bad, bad[:5]


# --- assay edges ------------------------------------------------------------------------

def test_assay_edges(run_query):
    assert not run_query("MATCH ()-[e:Assay_quantifies_metabolite|Assay_flags_metabolite]->() "
                         "WHERE e.stereo_assignment IS NULL RETURN e LIMIT 1")
    dup = run_query("""MATCH (a:MetaboliteAssay)-[e:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m)
        WITH a, m, e.condition_label AS c, e.time_point AS t, count(e) AS n WHERE n > 1
        RETURN a.id AS assay, m.id AS metabolite, n LIMIT 5""")
    assert not dup, dup


# --- MetaboliteClass ------------------------------------------------------------------

def test_metabolite_class(run_query):
    assert not run_query("MATCH (c:MetaboliteClass) WHERE c.metabolite_count IS NULL "
                         "OR c.metabolite_count < c.direct_metabolite_count RETURN c.id LIMIT 5")
    # is_a is child -> parent: a child is never shallower than its parent
    assert not run_query("MATCH (c:MetaboliteClass)-[:Metabolite_class_is_a_metabolite_class]->(p:MetaboliteClass) "
                         "WHERE c.level <= p.level RETURN c.id, p.id LIMIT 5")
    from multiomics_kg.utils.metabolite_identity.classes import BR08001_LEVEL0
    cats = {x["c"] for x in run_query("MATCH (m:Metabolite) UNWIND m.kegg_compound_categories AS c "
                                      "RETURN DISTINCT c")}
    assert cats, "no kegg_compound_categories"
    assert cats <= set(BR08001_LEVEL0), cats - set(BR08001_LEVEL0)


def test_new_edge_sources_join_data_sources(run_query):
    n = _one(run_query, "MATCH ()-[e:Metabolite_in_metabolite_class|Metabolite_is_stereo_form_of"
                        "|Metabolite_is_enantiomer_of]->() RETURN count(e) AS n")["n"]
    assert n > 0
    # UNWIND of a null/empty list emits no rows, so check those explicitly
    assert not run_query("MATCH ()-[e:Metabolite_in_metabolite_class|Metabolite_is_stereo_form_of"
                         "|Metabolite_is_enantiomer_of]->() WHERE e.sources IS NULL OR size(e.sources) = 0 "
                         "RETURN type(e) AS t LIMIT 5")
    bad = run_query("""MATCH ()-[e:Metabolite_in_metabolite_class|Metabolite_is_stereo_form_of|Metabolite_is_enantiomer_of]->()
        UNWIND e.sources AS s WITH DISTINCT s WHERE NOT EXISTS { MATCH (d:DataSource {id: 'data_source:' + s}) }
        RETURN s""")
    assert not bad, bad


# --- full-text --------------------------------------------------------------------------

# Bare ids checked present in cache/data/kegg/kegg_data.json all_identifiers (2026-10-14):
#   C01353      -> kegg.compound:C00288 (merged in, mnx_same_entity)
#   MNXM732398  -> kegg.compound:C00244 (nitrate)
#   HMDB0000148 -> kegg.compound:C00025 (L-glutamate). The brief's HMDB0002878 is on no
#                  node after the identity change, so an HMDB id that is was substituted.
FULLTEXT_PROBES = {
    "C01353": "kegg.compound:C00288",
    "MNXM732398": "kegg.compound:C00244",
    "HMDB0000148": "kegg.compound:C00025",
}


def test_fulltext_hits(run_query):
    for q, expected in FULLTEXT_PROBES.items():
        hits = [r["id"] for r in run_query(
            "CALL db.index.fulltext.queryNodes('metaboliteFullText', $q) YIELD node "
            "RETURN node.id AS id LIMIT 5", q=q)]
        assert expected in hits, (q, hits)
    assert run_query("CALL db.index.fulltext.queryNodes('chemicalRoleFullText', 'siderophore') YIELD node "
                     "RETURN node.id")


# --- TCDB substrate_breadth pin -------------------------------------------------------

def test_lumping_pin(run_query):
    pin = json.loads(_need(TCDB_PRUNED).read_text())
    expected = set(pin["lumping_pin"]["lumping_tcdb_ids"]) & set(pin["kept_tcdb_ids"])
    assert expected
    live = {r["id"] for r in run_query("MATCH (t:TcdbFamily {substrate_breadth: 'lumping'}) RETURN t.tcdb_id AS id")}
    assert live == expected, (sorted(expected - live)[:5], sorted(live - expected)[:5])
    # dense
    assert not run_query("MATCH (t:TcdbFamily) WHERE t.substrate_breadth IS NULL RETURN t.tcdb_id LIMIT 3")
    assert not run_query("MATCH (t:TcdbFamily) WHERE t.level < 2 AND t.substrate_breadth <> 'not_assessed' "
                         "RETURN t.tcdb_id LIMIT 3")


def test_transport_substrate_resolution_unchanged(run_query):
    baseline = json.loads(_need(PRE_COUNTS).read_text())["transport_substrate_resolution"]
    expected = {k: v for k, v in baseline.items() if k != "None"}   # null stored as "None"
    dist = {r["v"]: r["n"] for r in run_query("MATCH (g:Gene) WHERE g.transport_substrate_resolution IS NOT NULL "
                                              "RETURN g.transport_substrate_resolution AS v, count(*) AS n")}
    assert dist == expected
