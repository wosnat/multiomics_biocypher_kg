"""Fast schema validation tests for config/schema_config.yaml.

Verifies that the BioCypher schema loads without errors and contains
all expected node types and edge types.  Runs in < 2 seconds with no
network access and no file I/O beyond reading the two config YAML files.
"""

import pytest
from biocypher import BioCypher


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def schema():
    """Load the schema once for all tests in this module."""
    bc = BioCypher(
        biocypher_config_path="config/biocypher_config.yaml",
        schema_config_path="config/schema_config.yaml",
    )
    mapping = bc._get_ontology_mapping()
    return mapping.extended_schema


@pytest.fixture(scope="module")
def node_labels(schema):
    """Set of label_in_input values for all node entries."""
    return {
        v["label_in_input"]
        for v in schema.values()
        if v.get("represented_as") == "node" and "label_in_input" in v
    }


@pytest.fixture(scope="module")
def edge_labels(schema):
    """Set of label_in_input values for all edge entries."""
    return {
        v["label_in_input"]
        for v in schema.values()
        if v.get("represented_as") == "edge" and "label_in_input" in v
    }


# ---------------------------------------------------------------------------
# 1. Schema loads without errors
# ---------------------------------------------------------------------------

class TestSchemaLoads:
    def test_schema_is_not_empty(self, schema):
        assert len(schema) > 0, "Schema should contain at least one entry"

    def test_schema_has_nodes_and_edges(self, schema):
        reps = {v.get("represented_as") for v in schema.values()}
        assert "node" in reps, "Schema should define at least one node type"
        assert "edge" in reps, "Schema should define at least one edge type"


# ---------------------------------------------------------------------------
# 2. Key node types exist
# ---------------------------------------------------------------------------

class TestNodeTypes:
    @pytest.mark.parametrize("label", [
        "gene",
        "protein",
        "organism",
        "publication",
        "experiment",
        "ortholog_group",
        "pfam",
        "pfam_clan",
    ])
    def test_node_type_present(self, node_labels, label):
        assert label in node_labels, (
            f"Expected node label_in_input '{label}' not found in schema. "
            f"Available: {sorted(node_labels)}"
        )


# ---------------------------------------------------------------------------
# 3. Key edge types exist
# ---------------------------------------------------------------------------

class TestEdgeTypes:
    @pytest.mark.parametrize("label", [
        "changes_expression_of",
        "has_experiment",
        "tests_coculture_with",
        "gene_belongs_to_organism",
        "Gene_encodes_protein",
        "gene_in_ortholog_group",
        "gene_has_pfam",
        "pfam_in_pfam_clan",
    ])
    def test_edge_type_present(self, edge_labels, label):
        assert label in edge_labels, (
            f"Expected edge label_in_input '{label}' not found in schema. "
            f"Available: {sorted(edge_labels)}"
        )


# ---------------------------------------------------------------------------
# 4. Spot-check critical properties
# ---------------------------------------------------------------------------

class TestSchemaProperties:
    def test_gene_has_locus_tag(self, schema):
        gene = schema.get("gene", {})
        props = gene.get("properties", {})
        assert "locus_tag" in props, "Gene node must have locus_tag property"

    def test_expression_edge_has_log2fc(self, schema):
        for key, val in schema.items():
            if val.get("label_in_input") == "changes_expression_of":
                props = val.get("properties", {})
                assert "log2_fold_change" in props
                assert "adjusted_p_value" in props
                assert "expression_direction" in props
                assert "time_point" in props
                assert "time_point_order" in props
                assert "time_point_hours" in props
                break
        else:
            pytest.fail("changes_expression_of not found")

    def test_experiment_node_properties(self, schema):
        for key, val in schema.items():
            if val.get("label_in_input") == "experiment":
                props = val.get("properties", {})
                assert "name" in props
                assert "organism_name" in props
                assert "treatment_type" in props
                assert "treatment" in props
                assert "control" in props
                assert "omics_type" in props
                assert "is_time_course" in props
                assert "table_scope" in props
                assert "table_scope_detail" in props
                break
        else:
            pytest.fail("experiment node not found")

    def test_publication_has_doi(self, schema):
        for key, val in schema.items():
            if val.get("label_in_input") == "publication":
                assert val.get("preferred_id") == "doi"
                break
        else:
            pytest.fail("publication node not found")

    def test_organism_preferred_id(self, schema):
        for key, val in schema.items():
            if val.get("label_in_input") == "organism":
                assert val.get("preferred_id") == "insdc.gcf"
                break
        else:
            pytest.fail("organism node not found")


# ---------------------------------------------------------------------------
# 5. Removed types must NOT be in schema
# ---------------------------------------------------------------------------

class TestRemovedTypes:
    """Verify that old schema types have been removed."""

    def test_no_environmental_condition_node(self, node_labels):
        assert "environmental_condition" not in node_labels, (
            "environmental_condition node type should have been removed from schema"
        )

    def test_no_condition_changes_expression_of_edge(self, edge_labels):
        assert "condition_changes_expression_of" not in edge_labels, (
            "condition_changes_expression_of edge type should have been removed from schema"
        )

    def test_no_coculture_changes_expression_of_edge(self, edge_labels):
        assert "coculture_changes_expression_of" not in edge_labels, (
            "coculture_changes_expression_of edge type should have been removed from schema"
        )

    def test_no_published_expression_data_about_edge(self, edge_labels):
        assert "published_expression_data_about" not in edge_labels, (
            "published_expression_data_about edge type should have been removed from schema"
        )
