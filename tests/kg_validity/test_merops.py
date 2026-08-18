"""KG validity: MEROPS peptidase ontology (merops-diamond Phase-2).

Expected figures computed from the committed 42-strain calls.json corpus +
merops_reference.json before the rebuild (see
docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md): 155
MeropsFamily nodes (41 clans / 97 families / 17 subfamilies), 108 is-a edges,
<=4,254 Gene_has_merops_family edges (candidates whose WP_ maps to a gene),
call_class split ~3,425 peptidase / 777 nonpeptidase_homolog / 52 inhibitor.
"""

import pytest

pytestmark = pytest.mark.kg


# ── nodes ────────────────────────────────────────────────────────────────────

def test_merops_family_count_in_range(run_query):
    n = run_query("MATCH (m:MeropsFamily) RETURN count(m) AS n")[0]["n"]
    assert 140 <= n <= 175, f"expected MeropsFamily count in [140, 175], got {n}"


def test_merops_level_distribution(run_query):
    rows = run_query("""
        MATCH (m:MeropsFamily)
        RETURN m.level AS level, m.level_kind AS kind, count(*) AS n
        ORDER BY level
    """)
    by_level = {r["level"]: (r["kind"], r["n"]) for r in rows}
    assert set(by_level) == {0, 1, 2}
    assert by_level[0][0] == "merops_clan"
    assert by_level[1][0] == "merops_family"
    assert by_level[2][0] == "merops_subfamily"
    assert by_level[1][1] > by_level[2][1], "families must outnumber subfamilies"


def test_merops_node_ids_use_registered_prefixes(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily)
        WHERE NOT (m.id STARTS WITH 'merops.clan:' OR m.id STARTS WITH 'merops.family:')
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, "all MeropsFamily ids must use merops.clan:/merops.family: CURIEs"


def test_every_merops_node_named(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily)
        WHERE m.name IS NULL OR m.name = ''
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0


def test_s14_node_is_serine_clp(run_query):
    """S14 = ClpP endopeptidase family: serine peptidase, clan SK."""
    row = run_query("""
        MATCH (m:MeropsFamily {merops_id: 'S14'})
        OPTIONAL MATCH (m)-[:Merops_family_is_a_merops_family]->(p:MeropsFamily)
        RETURN m.catalytic_type AS ct, m.family_type AS ft, m.level AS level,
               m.name AS name, p.merops_id AS clan
    """)
    assert row, "S14 node not found"
    assert row[0]["ct"] == "serine"
    assert row[0]["ft"] == "peptidase"
    assert row[0]["level"] == 1
    assert row[0]["clan"] == "SK"
    assert "clp" in (row[0]["name"] or "").lower()


def test_inhibitor_families_typed_with_null_catalytic(run_query):
    rows = run_query("""
        MATCH (m:MeropsFamily)
        WHERE m.merops_id STARTS WITH 'I'
        RETURN count(*) AS n,
               sum(CASE WHEN m.family_type = 'inhibitor' THEN 1 ELSE 0 END) AS inhib,
               sum(CASE WHEN m.catalytic_type IS NULL THEN 1 ELSE 0 END) AS null_ct
    """)
    assert rows[0]["n"] > 0, "expected at least one inhibitor family/clan"
    assert rows[0]["inhib"] == rows[0]["n"]
    assert rows[0]["null_ct"] == rows[0]["n"]


def test_catalytic_type_vocabulary(run_query):
    rows = run_query("""
        MATCH (m:MeropsFamily)
        WHERE m.catalytic_type IS NOT NULL
        RETURN DISTINCT m.catalytic_type AS ct
    """)
    allowed = {"serine", "cysteine", "metallo", "aspartic", "threonine",
               "glutamic", "asparagine_lyase", "mixed", "unknown"}
    got = {r["ct"] for r in rows}
    assert got <= allowed, f"unexpected catalytic_type values: {got - allowed}"


def test_clan_descriptions_present(run_query):
    """Most clans carry the clan.txt fold/mechanism description (sparse prop)."""
    rows = run_query("""
        MATCH (m:MeropsFamily {level: 0})
        RETURN count(*) AS n,
               sum(CASE WHEN m.description IS NOT NULL AND m.description <> ''
                   THEN 1 ELSE 0 END) AS described
    """)
    assert rows[0]["described"] >= rows[0]["n"] * 0.8


# ── hierarchy ────────────────────────────────────────────────────────────────

def test_hierarchy_child_level_greater_than_parent(run_query):
    n = run_query("""
        MATCH (c:MeropsFamily)-[:Merops_family_is_a_merops_family]->(p:MeropsFamily)
        WHERE c.level <= p.level
        RETURN count(*) AS n
    """)[0]["n"]
    assert n == 0


def test_no_orphan_hierarchy_edges(run_query):
    n = run_query("""
        MATCH (a)-[r:Merops_family_is_a_merops_family]->(b)
        WHERE NOT a:MeropsFamily OR NOT b:MeropsFamily
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_subfamilies_all_have_parents(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily {level: 2})
        WHERE NOT EXISTS { (m)-[:Merops_family_is_a_merops_family]->() }
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, "every subfamily must link to its family"


# ── Gene_has_merops_family edges ─────────────────────────────────────────────

def test_gene_merops_edges_exist(run_query):
    n = run_query("MATCH ()-[r:Gene_has_merops_family]->() RETURN count(r) AS n")[0]["n"]
    assert 3800 <= n <= 4300, f"expected ~4.2K Gene_has_merops_family edges, got {n}"


def test_no_orphan_gene_merops_edges(run_query):
    n = run_query("""
        MATCH (g)-[r:Gene_has_merops_family]->(m)
        WHERE NOT g:Gene OR NOT m:MeropsFamily
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_edge_call_class_vocabulary_and_split(run_query):
    rows = run_query("""
        MATCH ()-[r:Gene_has_merops_family]->()
        RETURN r.call_class AS c, count(*) AS n
    """)
    by_class = {r["c"]: r["n"] for r in rows}
    assert set(by_class) <= {"peptidase", "inhibitor", "nonpeptidase_homolog"}
    assert by_class.get("peptidase", 0) > by_class.get("nonpeptidase_homolog", 0)
    assert by_class.get("inhibitor", 0) > 0


def test_edge_numeric_ranges(run_query):
    n = run_query("""
        MATCH ()-[r:Gene_has_merops_family]->()
        WHERE NOT (r.tier IN [1, 2, 3])
           OR r.confidence_score <= 0 OR r.confidence_score > 1
           OR r.identity <= 0 OR r.identity > 100
           OR r.qcov <= 0 OR r.qcov > 100
           OR r.evalue > 0.001
           OR r.consensus_n < 1
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_edge_best_hit_kind_vocabulary(run_query):
    rows = run_query("""
        MATCH ()-[r:Gene_has_merops_family]->()
        RETURN DISTINCT r.best_hit_kind AS k
    """)
    assert {r["k"] for r in rows} <= {"holotype", "putative", "nonpeptidase_homolog"}


def test_tier3_dominates(run_query):
    """~92% of candidates are tier-3 remote homology — expected shape, not a bug."""
    rows = run_query("""
        MATCH ()-[r:Gene_has_merops_family]->()
        RETURN r.tier AS tier, count(*) AS n
    """)
    by_tier = {r["tier"]: r["n"] for r in rows}
    total = sum(by_tier.values())
    assert by_tier.get(3, 0) / total > 0.8


def test_med4_clpp_spot_check(run_query):
    """PMM1313 (ClpP) → S14, tier 3, holotype-anchored peptidase call."""
    row = run_query("""
        MATCH (g:Gene {locus_tag: 'PMM1313'})-[r:Gene_has_merops_family]->(m:MeropsFamily)
        RETURN m.merops_id AS code, r.call_class AS cc, r.best_hit_id AS bh
    """)
    assert row, "PMM1313 has no Gene_has_merops_family edge"
    assert row[0]["code"] == "S14"
    assert row[0]["cc"] == "peptidase"
    assert row[0]["bh"].startswith("S14.")


# ── rollups ──────────────────────────────────────────────────────────────────

def test_rollup_counts_present_and_sane(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily)
        WHERE m.gene_count IS NULL OR m.organism_count IS NULL
           OR m.peptidase_gene_count IS NULL OR m.member_count IS NULL
           OR m.peptidase_gene_count > m.gene_count
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0


def test_every_node_reaches_a_gene(run_query):
    """Observed-only pruning + subtree rollup ⇒ every node has gene_count > 0."""
    n = run_query("""
        MATCH (m:MeropsFamily)
        WHERE coalesce(m.gene_count, 0) = 0
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0


def test_clan_subtree_geq_family(run_query):
    """A clan's subtree gene_count must be >= any child family's."""
    n = run_query("""
        MATCH (f:MeropsFamily)-[:Merops_family_is_a_merops_family]->(c:MeropsFamily)
        WHERE f.gene_count > c.gene_count
        RETURN count(*) AS n
    """)[0]["n"]
    assert n == 0


def test_homolog_prone_family_shows_count_gap(run_query):
    """C26 (gamma-glutamyl-hydrolase-like) is dominated by dead homologs —
    the peptidase_gene_count / gene_count gap must be visible."""
    row = run_query("""
        MATCH (m:MeropsFamily {merops_id: 'C26'})
        RETURN m.gene_count AS gc, m.peptidase_gene_count AS pgc
    """)
    assert row, "C26 node not found"
    assert row[0]["pgc"] < row[0]["gc"]


# ── Gene routing signals ─────────────────────────────────────────────────────

def test_gene_merops_count_and_classes_consistent(run_query):
    n = run_query("""
        MATCH (g:Gene)
        WHERE (g.merops_family_count > 0) <> (size(coalesce(g.merops_classes, [])) > 0)
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0


def test_gene_merops_classes_vocabulary(run_query):
    rows = run_query("""
        MATCH (g:Gene)
        WHERE size(coalesce(g.merops_classes, [])) > 0
        UNWIND g.merops_classes AS c
        RETURN DISTINCT c
    """)
    assert {r["c"] for r in rows} <= {"peptidase", "inhibitor", "nonpeptidase_homolog"}


def test_merops_in_annotation_types_is_tier_gated(run_query):
    """'merops' appears in annotation_types iff the gene has a tier<=2 edge."""
    n = run_query("""
        MATCH (g:Gene)
        WHERE ('merops' IN coalesce(g.annotation_types, []))
              <> EXISTS { MATCH (g)-[r:Gene_has_merops_family]->() WHERE r.tier <= 2 }
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0


def test_heterotrophs_carry_more_merops_genes_than_pro(run_query):
    """Lifestyle gradient (Phase-1 QC): heterotroph MIT1002 >> MED4."""
    rows = run_query("""
        MATCH (g:Gene)
        WHERE g.merops_family_count > 0
          AND g.organism_name IN ['Prochlorococcus MED4', 'Alteromonas macleodii MIT1002']
        RETURN g.organism_name AS org, count(g) AS n
    """)
    by_org = {r["org"]: r["n"] for r in rows}
    med4 = by_org.get("Prochlorococcus MED4", 0)
    mit1002 = by_org.get("Alteromonas macleodii MIT1002", 0)
    assert med4 > 0 and mit1002 > 0
    assert mit1002 > med4 * 2
