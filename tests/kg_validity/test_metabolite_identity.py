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
IDENTITY_REPORT = ROOT / "cache/data/metabolomics/metabolite_identity_report.json"

# all_identifiers CURIE forms (spec §3.2). The mnx local part may be non-MNXM (e.g. mnx:WATER).
CURIE_RE = (r"kegg\.compound:C\d{5}|chebi:\d+|mnx:\S+|hmdb:HMDB\d{7}"
            r"|inchikey:[A-Z]{14}-[A-Z]{10}-[A-Z]")
# Placeholder / generic formula: '*', an R group (R not starting Ra/Rb/Re/Rf/Rg/Rh/Ru), or a ")n" repeat.
PLACEHOLDER_FORMULA_RE = r".*(\*|R(?![abefghu])|\)n).*"

IDENTITY_BASES = {"mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual"}


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


def _identity_ground_truth(run_query):
    """Expected identity_basis carriers from the step-6 identity report (ground truth; the remap
    TSV's fallback labels are approximate release-file labels, not identity facts).

    Returns (merge_groups, rekey_targets, rekeys_by_target), restricted to WRITTEN nodes: a
    group / re-key whose canonical node is pruned at stage K (report `dropped_without_rows`,
    e.g. chebi:33283, kegg.compound:C05776) is not in the graph. A re-key target that was
    itself merged into another node is followed to that group's canonical id."""
    report = json.loads(_need(IDENTITY_REPORT).read_text())
    groups = report["merge_groups"]                       # {canonical: {basis: [...], members: [...]}}
    canon = {m: c for c, g in groups.items() for m in g["members"]}
    written = {r["id"] for r in run_query("MATCH (m:Metabolite) RETURN m.id AS id")}
    groups = {c: g for c, g in groups.items() if c in written}
    by_target: dict[str, list[str]] = {}
    for old, new in report["rekeys"].items():             # {chebi:<n>: kegg.compound:C...}
        t = canon.get(new, new)
        if t in written:
            by_target.setdefault(t, []).append(old)
    return groups, set(by_target), by_target


def test_identity_basis_only_on_merged_or_rekeyed(run_query):
    groups, rekey_targets, rekeys_by_target = _identity_ground_truth(run_query)
    assert groups and rekey_targets
    basis = {r["id"]: r["b"] for r in run_query(
        "MATCH (m:Metabolite) WHERE m.identity_basis IS NOT NULL RETURN m.id AS id, m.identity_basis AS b")}
    values = {v for b in basis.values() for v in b}
    assert values <= IDENTITY_BASES, values - IDENTITY_BASES
    assert not [i for i, b in basis.items() if not b], "empty identity_basis list"
    # exactly the merged and re-keyed nodes carry identity_basis
    expected = set(groups) | rekey_targets
    assert set(basis) == expected, (sorted(expected - set(basis))[:10], sorted(set(basis) - expected)[:10])
    # each merge group's canonical carries that group's bases
    bad_groups = [(c, g["basis"], basis.get(c)) for c, g in groups.items()
                  if not set(g["basis"]) <= set(basis.get(c, []))]
    assert not bad_groups, bad_groups[:10]
    # each re-key target carries chebi_kegg_xref (spec §3.2)
    no_xref = sorted(t for t in rekey_targets if "chebi_kegg_xref" not in basis.get(t, []))
    assert not no_xref, no_xref[:10]
    # every re-keyed old id is an identifier of its target (spec §3.2 all_identifiers)
    idents = {r["id"]: set(r["x"]) for r in run_query(
        "MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id AS id, m.all_identifiers AS x",
        ids=sorted(rekey_targets))}
    missing = [(o, t) for t, olds in rekeys_by_target.items() for o in olds if o not in idents.get(t, set())]
    assert not missing, missing[:10]


def test_every_rekeyed_node_has_kegg_enrichment(run_query):
    _, rekey_targets, _ = _identity_ground_truth(run_query)
    targets = sorted(rekey_targets)
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
    cats = {x["c"] for x in run_query("MATCH (m:Metabolite) UNWIND m.kegg_functional_categories AS c "
                                      "RETURN DISTINCT c")}
    assert cats, "no kegg_functional_categories"
    assert cats <= set(BR08001_LEVEL0), cats - set(BR08001_LEVEL0)
    assert not run_query("MATCH (m:Metabolite) WHERE m.kegg_compound_categories IS NOT NULL "
                         "RETURN m.id LIMIT 1"), "retired kegg_compound_categories still set"


def test_kegg_brite_classes_match_derivation(run_query):
    """Explorer ask 6 (2026-10-15): the property equals the derivation from MetaboliteClass,
    for every metabolite (absent when the metabolite is in no non-br08001 tree)."""
    import re
    rows = run_query("""
        MATCH (m:Metabolite)
        OPTIONAL MATCH (m)-[:Metabolite_in_metabolite_class]->(:MetaboliteClass)
              -[:Metabolite_class_is_a_metabolite_class*0..]->(top:MetaboliteClass {level: 0})
        WHERE top.tree_code <> 'br08001'
        RETURN m.id AS id, m.kegg_brite_classes AS got, collect(DISTINCT [top.tree, top.name]) AS tops""")
    bad = []
    for r in rows:
        want = set()
        for tree, name in r["tops"]:
            if tree is None:
                continue
            title = tree.replace("_", " ")
            title = title[0].upper() + title[1:]
            name = name[2:].strip() if re.fullmatch(r"[A-Z]{2}  .+", name) else name
            want.add(f"{title}: {name}")
        got = r["got"]
        if (sorted(want) or None) != got:
            bad.append((r["id"], got, sorted(want)))
    assert not bad, bad[:5]
    assert sum(1 for r in rows if r["got"]) > 500


# --- MNX stereo guard (explorer ask 2026-10-15) ------------------------------------------

# Pairs MNX called one entity that must stay two nodes: stereo-unspecified vs specific,
# generic sugar vs anomer, class vs member, two protein-bound species.
SPLIT_PAIRS = [("C00085", "C05345"), ("C00597", "C00197"), ("C00900", "C06010"),
               ("C05702", "C01102"), ("C00354", "C05378"), ("C03688", "C16240"),
               ("C04039", "C04272"), ("C04644", "C17345"), ("C00725", "C16241"),
               ("C02794", "C03227"), ("C00193", "C00261")]


@pytest.mark.parametrize("a,b", SPLIT_PAIRS)
def test_mnx_stereo_pairs_are_separate_nodes(run_query, a, b):
    rows = run_query("MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id AS id, m.kegg_compound_ids AS k",
                     ids=[f"kegg.compound:{a}", f"kegg.compound:{b}"])
    assert sorted(r["id"] for r in rows) == sorted([f"kegg.compound:{a}", f"kegg.compound:{b}"])
    for r in rows:
        assert r["k"] == [r["id"].split(":")[1]], r


@pytest.mark.parametrize("specific,generic", [("C05345", "C00085"), ("C05378", "C00354"),
                                              ("C00197", "C00597"), ("C06010", "C00900")])
def test_split_pairs_are_stereo_linked(run_query, specific, generic):
    assert run_query("MATCH (:Metabolite {id: $s})-[:Metabolite_is_stereo_form_of]->(:Metabolite {id: $g}) "
                     "RETURN 1 AS x", s=f"kegg.compound:{specific}", g=f"kegg.compound:{generic}")


def test_halide_pins_name_the_ion(run_query):
    rows = {r["id"]: r["k"] for r in run_query(
        "MATCH (m:Metabolite) WHERE m.id IN $ids RETURN m.id AS id, m.kegg_compound_ids AS k",
        ids=["kegg.compound:C00698", "kegg.compound:C01324", "kegg.compound:C00742"])}
    assert rows == {"kegg.compound:C00698": ["C00698", "C01327"],
                    "kegg.compound:C01324": ["C01324", "C13645"],
                    "kegg.compound:C00742": ["C00742", "C16487"]}


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


# Curated xref: overrides (explorer ask 5): generic HMDB / ChEBI records sit on the generic node.
@pytest.mark.parametrize("ident,node", [
    ("hmdb:HMDB0000122", "kegg.compound:C00031"), ("chebi:17634", "kegg.compound:C00031"),
    ("hmdb:HMDB0000143", "kegg.compound:C00124"), ("hmdb:HMDB0000660", "kegg.compound:C00095"),
    ("hmdb:HMDB0001058", "kegg.compound:C00354"), ("hmdb:HMDB0000124", "kegg.compound:C00085")])
def test_generic_xrefs_on_generic_node(run_query, ident, node):
    rows = run_query("MATCH (m:Metabolite) WHERE $i IN m.all_identifiers RETURN m.id AS id", i=ident)
    assert [r["id"] for r in rows] == [node]
