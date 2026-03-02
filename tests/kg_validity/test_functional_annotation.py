"""
Functional annotation tests for the multi-omics knowledge graph.

Validates GO, EC, and KEGG nodes and edges produced by
MultiGoAnnotationAdapter, MultiEcAnnotationAdapter (via ec_adapter), and
MultiKeggAnnotationAdapter (functional_annotation_adapter.py).

Checks:
- Node type presence and minimum counts
- Name property populated on all annotation nodes
- Hierarchy completeness (GO is-a, EC is-a, KEGG BRITE hierarchy)
- Gene annotation edge counts and per-gene coverage thresholds
- Correct node types connected by each edge type
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# GO: node counts and properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,min_count", [
    ("BiologicalProcess", 1500),
    ("CellularComponent", 200),
    ("MolecularFunction", 5000),
])
def test_go_node_count_minimum(run_query, label, min_count):
    """GO term nodes must meet minimum count (closure includes ancestors)."""
    result = run_query(f"MATCH (n:{label}) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt >= min_count, (
        f"Only {cnt} {label} nodes found; expected >= {min_count}"
    )


def test_go_root_terms_present(run_query):
    """The three GO root terms must always be present as ancestor closure nodes."""
    roots = [
        ("BiologicalProcess", "go:0008150", "biological_process"),
        ("CellularComponent", "go:0005575", "cellular_component"),
        ("MolecularFunction", "go:0003674", "molecular_function"),
    ]
    for label, go_id, expected_name in roots:
        result = run_query(
            f"MATCH (n:{label}) WHERE n.id = $nid RETURN n.name AS name",
            nid=go_id,
        )
        assert len(result) == 1, (
            f"GO root term {go_id} ({label}) not found in graph"
        )
        assert result[0]["name"] == expected_name, (
            f"GO root {go_id}: expected name '{expected_name}', got '{result[0]['name']}'"
        )


@pytest.mark.parametrize("label", [
    "BiologicalProcess",
    "CellularComponent",
    "MolecularFunction",
])
def test_go_nodes_have_name(run_query, label):
    """All GO term nodes must have a name property."""
    result = run_query(
        f"MATCH (n:{label}) WHERE n.name IS NULL RETURN count(n) AS missing"
    )
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} {label} nodes are missing the name property"


# ---------------------------------------------------------------------------
# GO: hierarchy edges
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel_type,source_label,target_label,min_count", [
    ("Biological_process_is_a_biological_process", "BiologicalProcess", "BiologicalProcess", 2000),
    ("Cellular_component_is_a_cellular_component", "CellularComponent",  "CellularComponent",  200),
    ("Molecular_function_is_a_molecular_function", "MolecularFunction",  "MolecularFunction",  1000),
])
def test_go_hierarchy_edges(run_query, rel_type, source_label, target_label, min_count):
    """GO is-a hierarchy edges must exist with sufficient count."""
    result = run_query(
        f"MATCH ()-[r:`{rel_type}`]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= min_count, (
        f"Only {cnt} {rel_type} edges found; expected >= {min_count}"
    )


@pytest.mark.parametrize("rel_type,source_label,target_label", [
    ("Biological_process_is_a_biological_process", "BiologicalProcess", "BiologicalProcess"),
    ("Cellular_component_is_a_cellular_component", "CellularComponent",  "CellularComponent"),
    ("Molecular_function_is_a_molecular_function", "MolecularFunction",  "MolecularFunction"),
])
def test_go_hierarchy_connects_correct_types(run_query, rel_type, source_label, target_label):
    """GO hierarchy edges must connect nodes of the same GO namespace."""
    result = run_query(
        f"MATCH (src)-[r:`{rel_type}`]->(tgt) "
        f"WHERE NOT src:{source_label} OR NOT tgt:{target_label} "
        f"RETURN count(r) AS bad"
    )
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} {rel_type} edges connect wrong node types "
        f"(expected {source_label} -> {target_label})"
    )


# ---------------------------------------------------------------------------
# GO: gene annotation edges
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel_type,target_label,min_count", [
    ("Gene_involved_in_biological_process", "BiologicalProcess", 50_000),
    ("Gene_located_in_cellular_component",  "CellularComponent",  20_000),
    ("Gene_enables_molecular_function",      "MolecularFunction",  50_000),
])
def test_go_gene_annotation_edge_counts(run_query, rel_type, target_label, min_count):
    """Gene→GO annotation edges must meet minimum count across all strains."""
    result = run_query(
        f"MATCH ()-[r:`{rel_type}`]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= min_count, (
        f"Only {cnt} {rel_type} edges found; expected >= {min_count}"
    )


@pytest.mark.parametrize("rel_type,target_label", [
    ("Gene_involved_in_biological_process", "BiologicalProcess"),
    ("Gene_located_in_cellular_component",  "CellularComponent"),
    ("Gene_enables_molecular_function",      "MolecularFunction"),
])
def test_go_gene_edges_connect_correct_types(run_query, rel_type, target_label):
    """Gene→GO edges must have Gene as source and the correct GO namespace as target."""
    result = run_query(
        f"MATCH (src)-[r:`{rel_type}`]->(tgt) "
        f"WHERE NOT src:Gene OR NOT tgt:{target_label} "
        f"RETURN count(r) AS bad"
    )
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} {rel_type} edges connect wrong node types "
        f"(expected Gene -> {target_label})"
    )


def test_go_gene_coverage(run_query):
    """
    At least 50% of all Gene nodes should carry at least one GO annotation
    (across BP, CC, or MF). A lower fraction indicates a data pipeline failure.
    Actual coverage is ~65%.
    """
    result = run_query("""
        MATCH (g:Gene)
        WITH count(g) AS total,
             count(CASE WHEN exists((g)-[:Gene_involved_in_biological_process|
                                        Gene_located_in_cellular_component|
                                        Gene_enables_molecular_function]->())
                        THEN 1 END) AS annotated
        RETURN total, annotated
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Gene nodes found")
    coverage = row["annotated"] / row["total"]
    assert coverage >= 0.50, (
        f"GO annotation coverage is {coverage:.1%} "
        f"({row['annotated']}/{row['total']}); expected >= 50%"
    )


# ---------------------------------------------------------------------------
# EC: node counts and properties
# ---------------------------------------------------------------------------

def test_ec_node_count_minimum(run_query):
    """EcNumber hierarchy nodes must exceed 5000 (full Expasy hierarchy)."""
    result = run_query("MATCH (n:EcNumber) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt >= 5000, (
        f"Only {cnt} EcNumber nodes found; expected >= 5000"
    )


def test_ec_nodes_have_name(run_query):
    """All EcNumber nodes must have a name property."""
    result = run_query(
        "MATCH (n:EcNumber) WHERE n.name IS NULL RETURN count(n) AS missing"
    )
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} EcNumber nodes are missing the name property"


@pytest.mark.parametrize("ec_id", [
    "ec:1.-.-.-",  # Oxidoreductases
    "ec:2.-.-.-",  # Transferases
    "ec:3.-.-.-",  # Hydrolases
    "ec:4.-.-.-",  # Lyases
    "ec:5.-.-.-",  # Isomerases
    "ec:6.-.-.-",  # Ligases
    "ec:7.-.-.-",  # Translocases
])
def test_ec_root_classes_present(run_query, ec_id):
    """All 7 top-level EC classes must be present in the hierarchy."""
    result = run_query(
        "MATCH (n:EcNumber) WHERE n.id = $nid RETURN count(n) AS cnt",
        nid=ec_id,
    )
    assert result[0]["cnt"] == 1, (
        f"EC root class {ec_id} not found in graph"
    )


# ---------------------------------------------------------------------------
# EC: hierarchy and gene annotation edges
# ---------------------------------------------------------------------------

def test_ec_hierarchy_edges_exist(run_query):
    """Ec_number_is_a_ec_number hierarchy edges must exceed 5000."""
    result = run_query(
        "MATCH ()-[r:Ec_number_is_a_ec_number]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 5000, (
        f"Only {cnt} Ec_number_is_a_ec_number edges found; expected >= 5000"
    )


def test_ec_hierarchy_connects_correct_types(run_query):
    """Ec_number_is_a_ec_number must connect EcNumber -> EcNumber."""
    result = run_query("""
        MATCH (src)-[r:Ec_number_is_a_ec_number]->(tgt)
        WHERE NOT src:EcNumber OR NOT tgt:EcNumber
        RETURN count(r) AS bad
    """)
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} Ec_number_is_a_ec_number edges connect wrong node types"
    )


def test_gene_catalyzes_ec_edge_count(run_query):
    """Gene_catalyzes_ec_number edges must exceed 10,000."""
    result = run_query(
        "MATCH ()-[r:Gene_catalyzes_ec_number]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 10_000, (
        f"Only {cnt} Gene_catalyzes_ec_number edges found; expected >= 10,000"
    )


def test_gene_catalyzes_ec_connects_correct_types(run_query):
    """Gene_catalyzes_ec_number must connect Gene -> EcNumber."""
    result = run_query("""
        MATCH (src)-[r:Gene_catalyzes_ec_number]->(tgt)
        WHERE NOT src:Gene OR NOT tgt:EcNumber
        RETURN count(r) AS bad
    """)
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} Gene_catalyzes_ec_number edges connect wrong node types "
        f"(expected Gene -> EcNumber)"
    )


def test_ec_gene_coverage(run_query):
    """
    At least 30% of Gene nodes should have an EC annotation.
    Actual coverage is ~36%. Genes lacking eggNOG annotations have no EC.
    """
    result = run_query("""
        MATCH (g:Gene)
        WITH count(g) AS total,
             count(CASE WHEN exists((g)-[:Gene_catalyzes_ec_number]->()) THEN 1 END) AS annotated
        RETURN total, annotated
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Gene nodes found")
    coverage = row["annotated"] / row["total"]
    assert coverage >= 0.30, (
        f"EC annotation coverage is {coverage:.1%} "
        f"({row['annotated']}/{row['total']}); expected >= 30%"
    )


# ---------------------------------------------------------------------------
# KEGG: node counts and properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,min_count", [
    ("KeggOrthologousGroup", 1500),
    ("KeggPathway",          100),
    ("KeggSubcategory",      30),
    ("KeggCategory",         5),
])
def test_kegg_node_count_minimum(run_query, label, min_count):
    """KEGG hierarchy nodes must meet minimum count."""
    result = run_query(f"MATCH (n:{label}) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt >= min_count, (
        f"Only {cnt} {label} nodes found; expected >= {min_count}"
    )


@pytest.mark.parametrize("label", [
    "KeggOrthologousGroup",
    "KeggPathway",
    "KeggSubcategory",
    "KeggCategory",
])
def test_kegg_nodes_have_name(run_query, label):
    """All KEGG hierarchy nodes must have a name property."""
    result = run_query(
        f"MATCH (n:{label}) WHERE n.name IS NULL RETURN count(n) AS missing"
    )
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} {label} nodes are missing the name property"


# ---------------------------------------------------------------------------
# KEGG: hierarchy edges
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel_type,source_label,target_label,min_count", [
    ("Ko_in_kegg_pathway",                 "KeggOrthologousGroup", "KeggPathway",          3000),
    ("Kegg_pathway_in_kegg_subcategory",   "KeggPathway",          "KeggSubcategory",       100),
    ("Kegg_subcategory_in_kegg_category",  "KeggSubcategory",      "KeggCategory",           30),
])
def test_kegg_hierarchy_edges(run_query, rel_type, source_label, target_label, min_count):
    """KEGG BRITE hierarchy edges must exist with sufficient count."""
    result = run_query(
        f"MATCH ()-[r:`{rel_type}`]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= min_count, (
        f"Only {cnt} {rel_type} edges found; expected >= {min_count}"
    )


@pytest.mark.parametrize("rel_type,source_label,target_label", [
    ("Ko_in_kegg_pathway",                 "KeggOrthologousGroup", "KeggPathway"),
    ("Kegg_pathway_in_kegg_subcategory",   "KeggPathway",          "KeggSubcategory"),
    ("Kegg_subcategory_in_kegg_category",  "KeggSubcategory",      "KeggCategory"),
])
def test_kegg_hierarchy_connects_correct_types(run_query, rel_type, source_label, target_label):
    """KEGG hierarchy edges must connect the correct node types."""
    result = run_query(
        f"MATCH (src)-[r:`{rel_type}`]->(tgt) "
        f"WHERE NOT src:{source_label} OR NOT tgt:{target_label} "
        f"RETURN count(r) AS bad"
    )
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} {rel_type} edges connect wrong node types "
        f"(expected {source_label} -> {target_label})"
    )


# ---------------------------------------------------------------------------
# KEGG: gene annotation edges
# ---------------------------------------------------------------------------

def test_gene_kegg_ko_edge_count(run_query):
    """Gene_has_kegg_ko edges must exceed 10,000."""
    result = run_query(
        "MATCH ()-[r:Gene_has_kegg_ko]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 10_000, (
        f"Only {cnt} Gene_has_kegg_ko edges found; expected >= 10,000"
    )


def test_gene_kegg_ko_connects_correct_types(run_query):
    """Gene_has_kegg_ko must connect Gene -> KeggOrthologousGroup."""
    result = run_query("""
        MATCH (src)-[r:Gene_has_kegg_ko]->(tgt)
        WHERE NOT src:Gene OR NOT tgt:KeggOrthologousGroup
        RETURN count(r) AS bad
    """)
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} Gene_has_kegg_ko edges connect wrong node types "
        f"(expected Gene -> KeggOrthologousGroup)"
    )


def test_kegg_gene_coverage(run_query):
    """
    At least 40% of Gene nodes should have a KEGG KO annotation.
    Actual coverage is ~51%. Genes lacking eggNOG annotations have no KO.
    """
    result = run_query("""
        MATCH (g:Gene)
        WITH count(g) AS total,
             count(CASE WHEN exists((g)-[:Gene_has_kegg_ko]->()) THEN 1 END) AS annotated
        RETURN total, annotated
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Gene nodes found")
    coverage = row["annotated"] / row["total"]
    assert coverage >= 0.40, (
        f"KEGG KO annotation coverage is {coverage:.1%} "
        f"({row['annotated']}/{row['total']}); expected >= 40%"
    )


def test_kegg_ko_pathway_coverage(run_query):
    """
    Most KEGG KO nodes observed in gene annotations should be linked to
    at least one KEGG pathway. Require >= 50% pathway linkage.
    Unmapped KOs exist in KEGG (e.g., hypothetical proteins) and are expected.
    """
    result = run_query("""
        MATCH (ko:KeggOrthologousGroup)
        WITH count(ko) AS total,
             count(CASE WHEN exists((ko)-[:Ko_in_kegg_pathway]->()) THEN 1 END) AS linked
        RETURN total, linked
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No KeggOrthologousGroup nodes found")
    linked_fraction = row["linked"] / row["total"]
    assert linked_fraction >= 0.50, (
        f"Only {row['linked']}/{row['total']} KO nodes ({linked_fraction:.1%}) "
        f"are linked to a pathway; expected >= 50%"
    )


# ---------------------------------------------------------------------------
# COG functional categories: node counts, properties, and gene edges
# ---------------------------------------------------------------------------

def test_cog_node_count(run_query):
    """All 26 standard COG functional category letters (A-Z) must be present."""
    result = run_query("MATCH (n:CogFunctionalCategory) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt == 26, (
        f"Expected exactly 26 CogFunctionalCategory nodes (A-Z), found {cnt}"
    )


@pytest.mark.parametrize("letter,expected_name", [
    ("J", "Translation, ribosomal structure and biogenesis"),
    ("C", "Energy production and conversion"),
    ("S", "Function unknown"),
    ("R", "General function prediction only"),
])
def test_cog_known_categories_present(run_query, letter, expected_name):
    """Spot-check well-known COG category nodes by letter code."""
    result = run_query(
        "MATCH (n:CogFunctionalCategory) WHERE n.code = $code RETURN n.id AS nid, n.name AS name",
        code=letter,
    )
    assert len(result) == 1, f"COG category '{letter}' not found"
    assert result[0]["name"] == expected_name, (
        f"COG {letter}: expected name '{expected_name}', got '{result[0]['name']}'"
    )


def test_cog_nodes_have_code_and_name(run_query):
    """All CogFunctionalCategory nodes must have both code and name properties."""
    result = run_query("""
        MATCH (n:CogFunctionalCategory)
        WHERE n.code IS NULL OR n.name IS NULL
        RETURN count(n) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} CogFunctionalCategory nodes are missing code or name"
    )


def test_gene_cog_category_edge_count(run_query):
    """Gene_in_cog_category edges must exceed 20,000 (all 13 strains annotated)."""
    result = run_query(
        "MATCH ()-[r:Gene_in_cog_category]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 20_000, (
        f"Only {cnt} Gene_in_cog_category edges found; expected >= 20,000"
    )


def test_gene_cog_edges_connect_correct_types(run_query):
    """Gene_in_cog_category must connect Gene -> CogFunctionalCategory."""
    result = run_query("""
        MATCH (src)-[r:Gene_in_cog_category]->(tgt)
        WHERE NOT src:Gene OR NOT tgt:CogFunctionalCategory
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Gene_in_cog_category edges connect wrong node types "
        f"(expected Gene -> CogFunctionalCategory)"
    )


def test_cog_gene_coverage(run_query):
    """
    At least 60% of Gene nodes should have a COG category annotation.
    All 13 strains are annotated via eggNOG; actual coverage is ~82%.
    """
    result = run_query("""
        MATCH (g:Gene)
        WITH count(g) AS total,
             count(CASE WHEN exists((g)-[:Gene_in_cog_category]->()) THEN 1 END) AS annotated
        RETURN total, annotated
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Gene nodes found")
    coverage = row["annotated"] / row["total"]
    assert coverage >= 0.60, (
        f"COG category coverage is {coverage:.1%} "
        f"({row['annotated']}/{row['total']}); expected >= 60%"
    )


# ---------------------------------------------------------------------------
# CyanorakRole: node counts, hierarchy, properties, and gene edges
# ---------------------------------------------------------------------------

def test_cyanorak_role_node_count(run_query):
    """The full Cyanorak role tree has ~172 nodes; expect at least 150."""
    result = run_query("MATCH (n:CyanorakRole) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt >= 150, (
        f"Only {cnt} CyanorakRole nodes found; expected >= 150 (full tree is 172)"
    )


def test_cyanorak_role_nodes_have_properties(run_query):
    """All CyanorakRole nodes must have both code and description properties."""
    result = run_query("""
        MATCH (n:CyanorakRole)
        WHERE n.code IS NULL OR n.description IS NULL
        RETURN count(n) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} CyanorakRole nodes are missing code or description"
    )


def test_cyanorak_role_hierarchy_edges(run_query):
    """Cyanorak_role_is_a_cyanorak_role hierarchy edges must cover most of the tree."""
    result = run_query(
        "MATCH ()-[r:Cyanorak_role_is_a_cyanorak_role]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 100, (
        f"Only {cnt} Cyanorak_role_is_a_cyanorak_role edges; expected >= 100"
    )


def test_cyanorak_role_hierarchy_connects_correct_types(run_query):
    """Cyanorak_role_is_a_cyanorak_role must connect CyanorakRole -> CyanorakRole."""
    result = run_query("""
        MATCH (src)-[r:Cyanorak_role_is_a_cyanorak_role]->(tgt)
        WHERE NOT src:CyanorakRole OR NOT tgt:CyanorakRole
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Cyanorak_role_is_a_cyanorak_role edges connect wrong node types"
    )


def test_gene_cyanorak_role_edge_count(run_query):
    """Gene_has_cyanorak_role edges must exceed 10,000 (Pro/Syn strains only)."""
    result = run_query(
        "MATCH ()-[r:Gene_has_cyanorak_role]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 10_000, (
        f"Only {cnt} Gene_has_cyanorak_role edges found; expected >= 10,000"
    )


def test_gene_cyanorak_role_edges_connect_correct_types(run_query):
    """Gene_has_cyanorak_role must connect Gene -> CyanorakRole."""
    result = run_query("""
        MATCH (src)-[r:Gene_has_cyanorak_role]->(tgt)
        WHERE NOT src:Gene OR NOT tgt:CyanorakRole
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Gene_has_cyanorak_role edges connect wrong node types "
        f"(expected Gene -> CyanorakRole)"
    )


# ---------------------------------------------------------------------------
# TigrRole: node counts, properties, and gene edges
# ---------------------------------------------------------------------------

def test_tigr_role_node_count(run_query):
    """TigrRole nodes (only codes present in data) must exceed 50."""
    result = run_query("MATCH (n:TigrRole) RETURN count(n) AS cnt")
    cnt = result[0]["cnt"]
    assert cnt >= 50, (
        f"Only {cnt} TigrRole nodes found; expected >= 50"
    )


def test_tigr_role_nodes_have_properties(run_query):
    """All TigrRole nodes must have both code and description properties."""
    result = run_query("""
        MATCH (n:TigrRole)
        WHERE n.code IS NULL OR n.description IS NULL
        RETURN count(n) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} TigrRole nodes are missing code or description"
    )


def test_gene_tigr_role_edge_count(run_query):
    """Gene_has_tigr_role edges must exceed 10,000 (Pro/Syn strains only)."""
    result = run_query(
        "MATCH ()-[r:Gene_has_tigr_role]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt >= 10_000, (
        f"Only {cnt} Gene_has_tigr_role edges found; expected >= 10,000"
    )


def test_gene_tigr_role_edges_connect_correct_types(run_query):
    """Gene_has_tigr_role must connect Gene -> TigrRole."""
    result = run_query("""
        MATCH (src)-[r:Gene_has_tigr_role]->(tgt)
        WHERE NOT src:Gene OR NOT tgt:TigrRole
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Gene_has_tigr_role edges connect wrong node types "
        f"(expected Gene -> TigrRole)"
    )
