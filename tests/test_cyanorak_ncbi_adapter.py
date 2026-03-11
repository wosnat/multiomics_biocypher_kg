"""Unit tests for the CyanorakNcbi adapter."""

import json
import os
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest
from pydantic import ValidationError

from multiomics_kg.adapters.cyanorak_ncbi_adapter import (
    ClusterNodeField,
    CyanorakNcbi,
    GeneEdgeType,
    GeneModel,
    GeneNodeField,
    MultiCyanorakNcbi,
    _fetch_ncbi_taxonomy,
)


# ---------------------------------------------------------------------------
# Minimal merged JSON test data
# ---------------------------------------------------------------------------

MERGED_JSON_DATA = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "locus_tag_ncbi": "TX50_RS00020",
        "locus_tag_cyanorak": "CK_Pro_MED4_00001",
        "protein_id": "WP_011131639.1",
        "old_locus_tags": ["PMM0001"],
        "start": 174,
        "end": 1331,
        "strand": "+",
        "start_cyanorak": 1331,
        "end_cyanorak": 1331,
        "strand_cyanorak": "+",
        "gene_name": "dnaN",
        "gene_synonyms": ["PMM0001"],
        "gene_name_source": "cyanorak",
        "product": "DNA polymerase III, beta subunit",
        "product_cyanorak": "DNA polymerase III, beta subunit",
        "product_source": "cyanorak",
        "function_description": "Replicative DNA polymerase beta subunit",
        "function_description_source": "uniprot",
        "cluster_number": "CK_00000364",
        "cyanorak_Role": ["F.1"],
        "cyanorak_Role_description": ["DNA replication"],
        "tIGR_Role": ["132"],
        "tIGR_Role_description": ["DNA metabolism"],
        "cog_category": "L",
        "protein_family": "Beta sliding clamp family",
        "go_terms": ["GO:0003887", "GO:0005737"],
        "go_term_descriptions": "DNA-directed DNA polymerase activity|cytoplasm",
        "ec_numbers": ["2.7.7.7"],
        "kegg_ko": ["K02338"],
        "kegg_ko_descriptions": ["DNA polymerase III, beta subunit"],
        "kegg_pathway": ["map00230", "map03030"],
        "kegg_module": ["M00260"],
        "kegg_reaction": ["R00375"],
        "kegg_brite": ["ko00001"],
        "eggnog_ogs": ["COG0592", "bactNOG00989"],
        "eggnog_og_descriptions": ["Replication, recombination"],
        "seed_ortholog": "59919.PMM0001",
        "max_annot_lvl": "1117|Cyanobacteria",
        "seed_ortholog_evalue": 1.12e-267,
        "pfam_ids": ["PF02768"],
        "pfam_descriptions": ["Beta sliding clamp, N-terminal domain"],
        "annotation_quality": 3,
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "locus_tag_ncbi": "TX50_RS00025",
        "protein_id": "WP_011131640.1",
        "start": 1333,
        "end": 2040,
        "strand": "+",
        "gene_name": "PMM0002",
        "gene_synonyms": ["PMM0002"],
        "product": "hypothetical protein",
        "cluster_number": "CK_00000363",
        "cog_category": "S",
        "go_terms": [],
        "ec_numbers": [],
        "kegg_ko": [],
        "eggnog_ogs": ["COG0243"],
        "annotation_quality": 1,
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "locus_tag_ncbi": "TX50_RS00030",
        "protein_id": "WP_011131641.1",
        "start": 2044,
        "end": 4383,
        "strand": "+",
        "gene_name": "purL",
        "gene_synonyms": ["PMM0003"],
        "product": "phosphoribosylformylglycinamidine synthase",
        "cluster_number": "CK_00000362",
        "cog_category": "F",
        "go_terms": ["GO:0009152"],
        "ec_numbers": ["6.3.5.3"],
        "kegg_ko": ["K01952"],
        "eggnog_ogs": ["COG0046"],
        "annotation_quality": 2,
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_data_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def merged_json_path(temp_data_dir):
    """Write MERGED_JSON_DATA to a temp gene_annotations_merged.json."""
    path = os.path.join(temp_data_dir, "gene_annotations_merged.json")
    with open(path, "w") as f:
        json.dump(MERGED_JSON_DATA, f)
    return path


@pytest.fixture
def adapter(temp_data_dir, merged_json_path):
    """CyanorakNcbi loaded from merged JSON, no taxonomy fetch."""
    a = CyanorakNcbi(
        data_dir=temp_data_dir,
        strain_name="MED4",
        ncbi_accession="GCF_000011465.1",
    )
    a.download_data()
    return a


@pytest.fixture
def adapter_no_prefix(temp_data_dir, merged_json_path):
    a = CyanorakNcbi(
        data_dir=temp_data_dir,
        strain_name="MED4",
        ncbi_accession="GCF_000011465.1",
        add_prefix=False,
    )
    a.download_data()
    return a


@pytest.fixture
def adapter_test_mode(temp_data_dir, merged_json_path):
    a = CyanorakNcbi(
        data_dir=temp_data_dir,
        strain_name="MED4",
        ncbi_accession="GCF_000011465.1",
        test_mode=True,
    )
    a.download_data()
    return a


@pytest.fixture
def adapter_subset_fields(temp_data_dir, merged_json_path):
    a = CyanorakNcbi(
        data_dir=temp_data_dir,
        strain_name="MED4",
        ncbi_accession="GCF_000011465.1",
        gene_node_fields=[
            GeneNodeField.LOCUS_TAG,
            GeneNodeField.PRODUCT,
            GeneNodeField.GENE_NAME,
        ],
    )
    a.download_data()
    return a


# ---------------------------------------------------------------------------
# GeneNodeField enum
# ---------------------------------------------------------------------------


class TestGeneNodeFieldEnum:
    def test_core_fields_exist(self):
        assert GeneNodeField.LOCUS_TAG.value == "locus_tag"
        assert GeneNodeField.PRODUCT.value == "product"
        assert GeneNodeField.EC_NUMBERS.value == "ec_numbers"

    def test_new_fields_exist(self):
        assert GeneNodeField.GENE_NAME.value == "gene_name"
        assert GeneNodeField.GENE_SYNONYMS.value == "gene_synonyms"
        assert GeneNodeField.GO_TERMS.value == "go_terms"
        assert GeneNodeField.KEGG_KO.value == "kegg_ko"
        assert GeneNodeField.EGGNOG_OGS.value == "eggnog_ogs"
        assert GeneNodeField.PFAM_IDS.value == "pfam_ids"
        assert GeneNodeField.COG_CATEGORY.value == "cog_category"
        assert GeneNodeField.ANNOTATION_QUALITY.value == "annotation_quality"
        assert GeneNodeField.SEED_ORTHOLOG_EVALUE.value == "seed_ortholog_evalue"

    def test_old_fields_removed(self):
        field_values = {f.value for f in GeneNodeField}
        assert "gene_names" not in field_values
        assert "eggNOG" not in field_values
        assert "kegg" not in field_values
        assert "protein_domains" not in field_values
        assert "go_component" not in field_values
        assert "Ontology_term" not in field_values

    def test_enum_contains_check(self):
        # GeneEnumMeta.__contains__ checks member names (uppercase), not values
        assert "LOCUS_TAG" in GeneNodeField
        assert "GENE_NAME" in GeneNodeField

    def test_enum_contains_nonexistent(self):
        assert "nonexistent_field" not in GeneNodeField

    def test_missing_case_insensitive_lookup(self):
        result = GeneNodeField._missing_("LOCUS_TAG")
        assert result == GeneNodeField.LOCUS_TAG

    def test_missing_returns_none_for_invalid(self):
        result = GeneNodeField._missing_("definitely_not_a_field")
        assert result is None


# ---------------------------------------------------------------------------
# GeneModel validation
# ---------------------------------------------------------------------------


class TestGeneModelValidation:
    def test_default_construction(self):
        model = GeneModel()
        assert model.gene_node_fields is None
        assert model.test_mode is False

    def test_test_mode_flag(self):
        model = GeneModel(test_mode=True)
        assert model.test_mode is True

    def test_organism_int(self):
        model = GeneModel(organism=59919)
        assert model.organism == 59919

    def test_organism_star(self):
        model = GeneModel(organism="*")
        assert model.organism == "*"

    def test_organism_none(self):
        model = GeneModel(organism=None)
        assert model.organism is None

    def test_invalid_output_dir(self):
        with pytest.raises(ValidationError):
            GeneModel(output_dir="/nonexistent/path/that/does/not/exist")

    def test_valid_output_dir(self, temp_data_dir):
        model = GeneModel(output_dir=temp_data_dir)
        assert model.output_dir is not None

    def test_gene_node_fields_list(self):
        model = GeneModel(gene_node_fields=[GeneNodeField.LOCUS_TAG, GeneNodeField.GENE_NAME])
        assert len(model.gene_node_fields) == 2


# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------


class TestConstructor:
    def test_default_node_fields(self, adapter):
        expected = {field.value for field in GeneNodeField}
        assert expected == set(adapter.gene_node_fields)

    def test_custom_node_fields(self, adapter_subset_fields):
        assert adapter_subset_fields.gene_node_fields == [
            GeneNodeField.LOCUS_TAG.value,
            GeneNodeField.PRODUCT.value,
            GeneNodeField.GENE_NAME.value,
        ]

    def test_default_edge_types(self, adapter):
        assert GeneEdgeType.GENE_IN_CLUSTER in adapter.edge_types

    def test_add_prefix_default_true(self, adapter):
        assert adapter.add_prefix is True

    def test_add_prefix_false(self, adapter_no_prefix):
        assert adapter_no_prefix.add_prefix is False

    def test_test_mode_sets_early_stopping(self, adapter_test_mode):
        assert adapter_test_mode.early_stopping == 100

    def test_normal_mode_no_early_stopping(self, adapter):
        assert adapter.early_stopping is None

    def test_strain_name_stored(self, adapter):
        assert adapter.strain_name == "MED4"

    def test_ncbi_accession_stored(self, adapter):
        assert adapter.ncbi_accession == "GCF_000011465.1"


# ---------------------------------------------------------------------------
# download_data()
# ---------------------------------------------------------------------------


class TestDownloadData:
    def test_loads_merged_json(self, adapter):
        assert adapter.data_df is not None
        assert isinstance(adapter.data_df, pd.DataFrame)

    def test_correct_row_count(self, adapter):
        assert len(adapter.data_df) == 3

    def test_locus_tag_column_present(self, adapter):
        assert "locus_tag" in adapter.data_df.columns

    def test_list_columns_preserved(self, adapter):
        row = adapter.data_df[adapter.data_df["locus_tag"] == "PMM0001"].iloc[0]
        assert isinstance(row["go_terms"], list)
        assert isinstance(row["kegg_ko"], list)
        assert isinstance(row["gene_synonyms"], list)

    def test_numeric_columns_preserved(self, adapter):
        row = adapter.data_df[adapter.data_df["locus_tag"] == "PMM0001"].iloc[0]
        assert row["annotation_quality"] == 3
        assert isinstance(row["seed_ortholog_evalue"], float)

    def test_missing_json_raises_error(self, temp_data_dir):
        a = CyanorakNcbi(data_dir=temp_data_dir, strain_name="TEST")
        with pytest.raises(FileNotFoundError, match="gene_annotations_merged.json"):
            a.download_data()

    def test_missing_data_dir_raises_error(self):
        a = CyanorakNcbi(strain_name="TEST")
        with pytest.raises(ValueError, match="data_dir is required"):
            a.download_data()


# ---------------------------------------------------------------------------
# clean_text
# ---------------------------------------------------------------------------


class TestCleanText:
    def test_pipe_replaced_with_comma(self, adapter):
        assert adapter.clean_text("a|b") == "a,b"

    def test_single_quote_replaced_with_caret(self, adapter):
        assert adapter.clean_text("it's") == "it^s"

    def test_both_special_chars(self, adapter):
        assert adapter.clean_text("a|b'c") == "a,b^c"

    def test_no_special_chars(self, adapter):
        assert adapter.clean_text("hello") == "hello"

    def test_list_input(self, adapter):
        result = adapter.clean_text(["a|b", "c'd"])
        assert result == ["a,b", "c^d"]

    def test_non_string_passthrough(self, adapter):
        assert adapter.clean_text(42) == 42

    def test_none_input(self, adapter):
        assert adapter.clean_text(None) is None

    def test_empty_string(self, adapter):
        assert adapter.clean_text("") == ""


# ---------------------------------------------------------------------------
# _get_split_character / _split_field
# ---------------------------------------------------------------------------


class TestSplitField:
    def test_cyanorak_role_comma_split(self, adapter):
        assert adapter._get_split_character("cyanorak_Role") == ","

    def test_tigr_role_comma_split(self, adapter):
        assert adapter._get_split_character("tIGR_Role") == ","

    def test_ec_numbers_comma_split(self, adapter):
        assert adapter._get_split_character("ec_numbers") == ","

    def test_product_no_split(self, adapter):
        assert adapter._get_split_character("product") is None

    def test_gene_name_no_split(self, adapter):
        assert adapter._get_split_character("gene_name") is None

    def test_go_terms_no_split(self, adapter):
        # go_terms come as lists from merged JSON, no CSV split needed
        assert adapter._get_split_character("go_terms") is None

    def test_split_comma_separated_string(self, adapter):
        result = adapter._split_field("cyanorak_Role", "F.1,F.2")
        assert result == ["F.1", "F.2"]

    def test_no_split_returns_string(self, adapter):
        result = adapter._split_field("product", "DNA polymerase")
        assert result == "DNA polymerase"

    def test_comma_followed_by_space_not_split(self, adapter):
        # Regex r',(?! )' does not split on comma+space, so "a, b, c" stays as one element
        result = adapter._split_field("cyanorak_Role_description", "a, b, c")
        assert result == ["a, b, c"]


# ---------------------------------------------------------------------------
# add_prefix_to_id
# ---------------------------------------------------------------------------


class TestAddPrefixToId:
    def test_add_prefix_enabled(self, adapter):
        result = adapter.add_prefix_to_id(prefix="ncbigene", identifier="PMM0001")
        assert "PMM0001" in result
        assert "ncbigene" in result

    def test_add_prefix_disabled(self, adapter_no_prefix):
        result = adapter_no_prefix.add_prefix_to_id(prefix="ncbigene", identifier="PMM0001")
        assert result == "PMM0001"

    def test_none_identifier(self, adapter):
        result = adapter.add_prefix_to_id(prefix="ncbigene", identifier=None)
        assert result is None

    def test_empty_identifier(self, adapter):
        result = adapter.add_prefix_to_id(prefix="ncbigene", identifier="")
        assert result == ""


# ---------------------------------------------------------------------------
# set_node_fields / set_edge_types
# ---------------------------------------------------------------------------


class TestSetNodeFields:
    def test_set_node_fields_with_list(self, adapter):
        adapter.set_node_fields([GeneNodeField.LOCUS_TAG, GeneNodeField.GENE_NAME])
        assert adapter.gene_node_fields == ["locus_tag", "gene_name"]

    def test_set_node_fields_none(self, adapter):
        adapter.set_node_fields(None)
        assert set(adapter.gene_node_fields) == {f.value for f in GeneNodeField}

    def test_set_node_fields_empty_list(self, adapter):
        adapter.set_node_fields([])
        assert set(adapter.gene_node_fields) == {f.value for f in GeneNodeField}

    def test_set_node_fields_single(self, adapter):
        adapter.set_node_fields([GeneNodeField.LOCUS_TAG])
        assert adapter.gene_node_fields == ["locus_tag"]


class TestSetEdgeTypes:
    def test_set_edge_types_none(self, adapter):
        adapter.set_edge_types(None)
        assert adapter.edge_types == list(GeneEdgeType)

    def test_set_edge_types_specific(self, adapter):
        adapter.set_edge_types([GeneEdgeType.GENE_IN_CLUSTER])
        assert GeneEdgeType.GENE_IN_CLUSTER in adapter.edge_types
        assert len(adapter.edge_types) == 1


# ---------------------------------------------------------------------------
# get_nodes() — gene nodes
# ---------------------------------------------------------------------------


class TestGetGeneNodes:
    def test_returns_list(self, adapter):
        nodes = adapter.get_nodes()
        assert isinstance(nodes, list)

    def test_node_counts(self, adapter):
        nodes = adapter.get_nodes()
        gene_nodes = [n for n in nodes if n[1] == "gene"]
        cluster_nodes = [n for n in nodes if n[1] == "cyanorak_cluster"]
        organism_nodes = [n for n in nodes if n[1] == "organism"]
        assert len(gene_nodes) == 3
        assert len(cluster_nodes) == 3
        assert len(organism_nodes) == 1

    def test_gene_node_structure(self, adapter):
        nodes = adapter.get_nodes()
        for node_id, label, props in [n for n in nodes if n[1] == "gene"]:
            assert isinstance(node_id, str)
            assert label == "gene"
            assert isinstance(props, dict)
            assert "locus_tag" in props

    def test_gene_id_uses_ncbigene_prefix(self, adapter):
        nodes = adapter.get_nodes()
        ids = [n[0] for n in nodes if n[1] == "gene"]
        assert all("PMM" in nid for nid in ids)

    def test_list_fields_preserved_as_lists(self, adapter):
        nodes = adapter.get_nodes()
        pmm0001 = next(n for n in nodes if n[1] == "gene" and "PMM0001" in n[0])
        props = pmm0001[2]
        assert isinstance(props.get("go_terms"), list)
        assert isinstance(props.get("kegg_ko"), list)
        assert isinstance(props.get("eggnog_ogs"), list)

    def test_numeric_fields_correct_types(self, adapter):
        nodes = adapter.get_nodes()
        pmm0001 = next(n for n in nodes if n[1] == "gene" and "PMM0001" in n[0])
        props = pmm0001[2]
        assert isinstance(props["start"], int)
        assert isinstance(props["end"], int)
        assert isinstance(props["annotation_quality"], int)
        assert isinstance(props["seed_ortholog_evalue"], float)

    def test_string_fields_correct_types(self, adapter):
        nodes = adapter.get_nodes()
        pmm0001 = next(n for n in nodes if n[1] == "gene" and "PMM0001" in n[0])
        props = pmm0001[2]
        assert isinstance(props["gene_name"], str)
        assert isinstance(props["cog_category"], str)

    def test_new_fields_populated(self, adapter):
        nodes = adapter.get_nodes()
        pmm0001 = next(n for n in nodes if n[1] == "gene" and "PMM0001" in n[0])
        props = pmm0001[2]
        assert props["gene_name"] == "dnaN"
        assert "K02338" in props["kegg_ko"]
        assert "GO:0003887" in props["go_terms"]
        assert props["cog_category"] == "L"
        assert props["seed_ortholog"] == "59919.PMM0001"

    def test_pipe_cleaned_from_string_values(self, adapter):
        nodes = adapter.get_nodes()
        pmm0001 = next(n for n in nodes if n[1] == "gene" and "PMM0001" in n[0])
        props = pmm0001[2]
        desc = props.get("go_term_descriptions", "")
        assert "|" not in desc

    def test_subset_fields_only_returns_requested_fields(self, adapter_subset_fields):
        nodes = adapter_subset_fields.get_nodes()
        for _, _, props in [n for n in nodes if n[1] == "gene"]:
            assert "locus_tag" in props
            assert "product" in props
            assert "gene_name" in props
            assert "cog_category" not in props


# ---------------------------------------------------------------------------
# get_edges()
# ---------------------------------------------------------------------------


class TestGetEdges:
    def test_gene_cluster_edges_present(self, adapter):
        edges = adapter.get_edges()
        gene_cluster = [e for e in edges if e[3] == "gene_in_cyanorak_cluster"]
        assert len(gene_cluster) == 3

    def test_gene_organism_edges_present(self, adapter):
        edges = adapter.get_edges()
        gene_org = [e for e in edges if e[3] == "gene_belongs_to_organism"]
        assert len(gene_org) == 3

    def test_edge_structure(self, adapter):
        for edge_id, source, target, edge_type, props in adapter.get_edges():
            assert isinstance(edge_id, str)
            assert isinstance(source, str)
            assert isinstance(target, str)
            assert isinstance(edge_type, str)
            assert isinstance(props, dict)


# ---------------------------------------------------------------------------
# Organism node
# ---------------------------------------------------------------------------


class TestOrganismNode:
    def test_organism_node_created(self, adapter):
        nodes = adapter.get_nodes()
        org_nodes = [n for n in nodes if n[1] == "organism"]
        assert len(org_nodes) == 1

    def test_organism_node_id_contains_accession(self, adapter):
        nodes = adapter.get_nodes()
        org_node = next(n for n in nodes if n[1] == "organism")
        assert "GCF_000011465.1" in org_node[0]

    def test_organism_node_strain_name(self, adapter):
        nodes = adapter.get_nodes()
        org_node = next(n for n in nodes if n[1] == "organism")
        assert org_node[2]["strain_name"] == "MED4"

    def test_organism_node_preferred_name(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
            preferred_name="Prochlorococcus MED4",
        )
        a.taxonomy = {}
        nodes = a._get_organism_node()
        assert nodes[0][2]["preferred_name"] == "Prochlorococcus MED4"

    def test_organism_node_no_preferred_name_when_not_set(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.taxonomy = {}
        nodes = a._get_organism_node()
        assert "preferred_name" not in nodes[0][2]

    def test_organism_node_taxonomy_from_download_data(self, temp_data_dir, merged_json_path):
        mock_tax = {"lineage": "Bacteria;Cyanobacteria", "genus": "Prochlorococcus"}
        with patch(
            "multiomics_kg.adapters.cyanorak_ncbi_adapter._fetch_ncbi_taxonomy",
            return_value=mock_tax,
        ):
            a = CyanorakNcbi(
                data_dir=temp_data_dir,
                strain_name="MED4",
                ncbi_accession="GCF_000011465.1",
                ncbi_taxon_id=59919,
            )
            a.download_data()
        nodes = a.get_nodes()
        org_node = next(n for n in nodes if n[1] == "organism")
        assert org_node[2].get("genus") == "Prochlorococcus"


# ---------------------------------------------------------------------------
# Organism node taxonomy
# ---------------------------------------------------------------------------


class TestOrganismNodeTaxonomy:
    def test_clade_in_organism_node(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
            clade="HLI",
        )
        a.taxonomy = {}
        nodes = a._get_organism_node()
        assert nodes[0][2]["clade"] == "HLI"

    def test_no_clade_not_in_props(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.taxonomy = {}
        nodes = a._get_organism_node()
        assert "clade" not in nodes[0][2]

    def test_taxonomy_ranks_propagated(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.taxonomy = {
            "lineage": "Bacteria;Cyanobacteria",
            "phylum": "Cyanobacteria",
            "genus": "Prochlorococcus",
            "species": "Prochlorococcus marinus",
        }
        nodes = a._get_organism_node()
        props = nodes[0][2]
        assert props["phylum"] == "Cyanobacteria"
        assert props["genus"] == "Prochlorococcus"
        assert props["species"] == "Prochlorococcus marinus"

    def test_empty_taxonomy_no_extra_props(self):
        a = CyanorakNcbi(
            data_dir="/fake/dir",
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.taxonomy = {}
        nodes = a._get_organism_node()
        props = nodes[0][2]
        assert "lineage" not in props
        assert "genus" not in props

    def test_no_ncbi_accession_returns_empty(self):
        a = CyanorakNcbi(data_dir="/fake/dir")
        a.taxonomy = {}
        nodes = a._get_organism_node()
        assert nodes == []


# ---------------------------------------------------------------------------
# _fetch_ncbi_taxonomy
# ---------------------------------------------------------------------------


MOCK_TAX_XML = """<?xml version="1.0" ?>
<TaxaSet>
  <Taxon>
    <Lineage>cellular organisms; Bacteria; Cyanobacteria/Melainabacteria group</Lineage>
    <LineageEx>
      <Taxon>
        <ScientificName>Bacteria</ScientificName>
        <Rank>superkingdom</Rank>
      </Taxon>
      <Taxon>
        <ScientificName>Cyanobacteria</ScientificName>
        <Rank>phylum</Rank>
      </Taxon>
      <Taxon>
        <ScientificName>Prochlorococcus</ScientificName>
        <Rank>genus</Rank>
      </Taxon>
      <Taxon>
        <ScientificName>Prochlorococcus marinus</ScientificName>
        <Rank>species</Rank>
      </Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>
"""


@pytest.fixture
def mock_curl_ok():
    mock = MagicMock()
    mock.result = MOCK_TAX_XML
    with patch("multiomics_kg.adapters.cyanorak_ncbi_adapter.curl.Curl", return_value=mock):
        yield mock


class TestFetchNcbiTaxonomy:
    def test_parses_lineage_string(self, tmp_path, mock_curl_ok):
        result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert "Bacteria" in result["lineage"]

    def test_superkingdom_parsed(self, tmp_path, mock_curl_ok):
        result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert result["superkingdom"] == "Bacteria"

    def test_phylum_parsed(self, tmp_path, mock_curl_ok):
        result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert result["phylum"] == "Cyanobacteria"

    def test_genus_and_species_parsed(self, tmp_path, mock_curl_ok):
        result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert result["genus"] == "Prochlorococcus"
        assert result["species"] == "Prochlorococcus marinus"

    def test_writes_cache_file(self, tmp_path, mock_curl_ok):
        _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert (tmp_path / "taxonomy_59919.json").exists()

    def test_cache_hit_skips_curl(self, tmp_path):
        cache_path = tmp_path / "taxonomy_59919.json"
        with open(cache_path, "w") as f:
            json.dump({"lineage": "cached", "genus": "CachedGenus"}, f)
        with patch("multiomics_kg.adapters.cyanorak_ncbi_adapter.curl.Curl") as mock_curl:
            result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
            mock_curl.assert_not_called()
        assert result["genus"] == "CachedGenus"

    def test_curl_failure_returns_empty_dict(self, tmp_path):
        mock = MagicMock()
        mock.result = None
        with patch("multiomics_kg.adapters.cyanorak_ncbi_adapter.curl.Curl", return_value=mock):
            result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert result == {}

    def test_malformed_xml_returns_empty_dict(self, tmp_path):
        mock = MagicMock()
        mock.result = "not valid xml <<<"
        with patch("multiomics_kg.adapters.cyanorak_ncbi_adapter.curl.Curl", return_value=mock):
            result = _fetch_ncbi_taxonomy(59919, str(tmp_path))
        assert result == {}


# ---------------------------------------------------------------------------
# MultiCyanorakNcbi
# ---------------------------------------------------------------------------


@pytest.fixture
def multi_data_dir(tmp_path):
    """Create two per-strain dirs with merged JSON files."""
    for strain in ["MED4", "MIT9312"]:
        strain_dir = tmp_path / strain
        strain_dir.mkdir()
        data = {
            f"{strain}_0001": {
                "locus_tag": f"{strain}_0001",
                "gene_name": f"gene_{strain}",
                "start": 100,
                "end": 500,
                "strand": "+",
                "product": "hypothetical protein",
                "cluster_number": f"CK_{strain}_001",
                "go_terms": [],
                "ec_numbers": [],
                "kegg_ko": [],
                "annotation_quality": 1,
            }
        }
        with open(strain_dir / "gene_annotations_merged.json", "w") as f:
            json.dump(data, f)
    return tmp_path


@pytest.fixture
def multi_config_csv(multi_data_dir, tmp_path):
    csv_path = tmp_path / "genomes.csv"
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade,preferred_name\n")
        f.write(f"GCF_000011465.1,{multi_data_dir}/MED4/,MED4,59919,HLI,Prochlorococcus MED4\n")
        f.write(f"GCF_000015645.1,{multi_data_dir}/MIT9312/,MIT9312,74546,HLII,Prochlorococcus MIT9312\n")
    return str(csv_path)


class TestMultiCyanorakNcbiConstruction:
    def test_loads_correct_number_of_adapters(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        assert len(multi.adapters) == 2

    def test_adapters_are_cyanorak_instances(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        for a in multi.adapters:
            assert isinstance(a, CyanorakNcbi)

    def test_kwargs_passed_to_adapters(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv, add_prefix=False)
        for a in multi.adapters:
            assert a.add_prefix is False

    def test_clade_set_on_adapters(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        clades = {a.clade for a in multi.adapters}
        assert "HLI" in clades
        assert "HLII" in clades

    def test_comment_lines_skipped(self, tmp_path, multi_data_dir):
        csv_path = tmp_path / "genomes_with_comments.csv"
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
            f.write("# This is a comment\n")
            f.write(f"GCF_000011465.1,{multi_data_dir}/MED4/,MED4,59919,HLI\n")
        multi = MultiCyanorakNcbi(config_list_file=str(csv_path))
        assert len(multi.adapters) == 1

    def test_empty_clade_sets_none(self, tmp_path, multi_data_dir):
        csv_path = tmp_path / "genomes_no_clade.csv"
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
            f.write(f"GCF_000011465.1,{multi_data_dir}/MED4/,MED4,59919,\n")
        multi = MultiCyanorakNcbi(config_list_file=str(csv_path))
        assert multi.adapters[0].clade is None

    def test_preferred_name_passed_from_csv(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        preferred_names = {a.strain_name: a.preferred_name for a in multi.adapters}
        assert preferred_names["MED4"] == "Prochlorococcus MED4"
        assert preferred_names["MIT9312"] == "Prochlorococcus MIT9312"


class TestMultiCyanorakNcbiGetNodes:
    def test_aggregates_nodes_from_all_genomes(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        multi.download_data()
        nodes = multi.get_nodes()
        gene_nodes = [n for n in nodes if n[1] == "gene"]
        assert len(gene_nodes) == 2

    def test_organism_nodes_present(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        multi.download_data()
        nodes = multi.get_nodes()
        org_nodes = [n for n in nodes if n[1] == "organism"]
        assert len(org_nodes) == 2


class TestTreatmentOrganismNodes:
    @pytest.fixture
    def treatment_csv(self, tmp_path):
        path = tmp_path / "treatment_organisms.csv"
        with open(path, "w") as f:
            f.write("ncbi_taxon_id,organism_name\n")
            f.write("10754,Phage\n")
            f.write("# comment\n")
            f.write("413470,Marinobacter\n")
        return str(path)

    def test_no_treatment_nodes_without_file(self, multi_config_csv):
        multi = MultiCyanorakNcbi(config_list_file=multi_config_csv)
        multi.download_data()
        nodes = multi.get_nodes()
        org_nodes = [n for n in nodes if n[1] == "organism"]
        assert len(org_nodes) == 2

    def test_treatment_nodes_created(self, multi_config_csv, treatment_csv):
        with patch(
            "multiomics_kg.adapters.cyanorak_ncbi_adapter._fetch_ncbi_taxonomy",
            return_value={},
        ):
            multi = MultiCyanorakNcbi(
                config_list_file=multi_config_csv,
                treatment_organisms_file=treatment_csv,
            )
            multi.download_data()
            nodes = multi.get_nodes()
        treatment_nodes = [n for n in nodes if "ncbitaxon" in n[0]]
        assert len(treatment_nodes) == 2

    def test_treatment_node_ids(self, multi_config_csv, treatment_csv):
        with patch(
            "multiomics_kg.adapters.cyanorak_ncbi_adapter._fetch_ncbi_taxonomy",
            return_value={},
        ):
            multi = MultiCyanorakNcbi(
                config_list_file=multi_config_csv,
                treatment_organisms_file=treatment_csv,
            )
            multi.download_data()
            nodes = multi.get_nodes()
        ids = [n[0] for n in nodes if "ncbitaxon" in n[0]]
        assert any("10754" in nid for nid in ids)
        assert any("413470" in nid for nid in ids)

    def test_treatment_node_preferred_name(self, multi_config_csv, treatment_csv):
        with patch(
            "multiomics_kg.adapters.cyanorak_ncbi_adapter._fetch_ncbi_taxonomy",
            return_value={},
        ):
            multi = MultiCyanorakNcbi(
                config_list_file=multi_config_csv,
                treatment_organisms_file=treatment_csv,
            )
            multi.download_data()
            nodes = multi.get_nodes()
        treatment_nodes = {n[0]: n[2] for n in nodes if "ncbitaxon" in n[0]}
        assert treatment_nodes["ncbitaxon:10754"]["preferred_name"] == "Phage"
        assert treatment_nodes["ncbitaxon:413470"]["preferred_name"] == "Marinobacter"

    def test_comments_skipped_in_treatment_csv(self, multi_config_csv, treatment_csv):
        with patch(
            "multiomics_kg.adapters.cyanorak_ncbi_adapter._fetch_ncbi_taxonomy",
            return_value={},
        ):
            multi = MultiCyanorakNcbi(
                config_list_file=multi_config_csv,
                treatment_organisms_file=treatment_csv,
            )
            multi.download_data()
            nodes = multi.get_nodes()
        treatment_nodes = [n for n in nodes if "ncbitaxon" in n[0]]
        assert len(treatment_nodes) == 2


# ---------------------------------------------------------------------------
# Integration: real cache data (skipped if not available)
# ---------------------------------------------------------------------------


@pytest.fixture
def real_data_dir():
    path = "cache/data/Prochlorococcus/genomes/MED4"
    if not os.path.exists(os.path.join(path, "gene_annotations_merged.json")):
        pytest.skip("Real MED4 cache data not available")
    return path


class TestIntegrationWithRealData:
    def test_real_data_loads(self, real_data_dir):
        a = CyanorakNcbi(
            data_dir=real_data_dir,
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.download_data()
        assert len(a.data_df) > 1000

    def test_real_data_gene_nodes(self, real_data_dir):
        a = CyanorakNcbi(
            data_dir=real_data_dir,
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.download_data()
        nodes = a.get_nodes()
        gene_nodes = [n for n in nodes if n[1] == "gene"]
        assert len(gene_nodes) > 1000

    def test_real_data_known_locus_tag(self, real_data_dir):
        a = CyanorakNcbi(
            data_dir=real_data_dir,
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.download_data()
        nodes = a.get_nodes()
        gene_nodes = {n[0]: n[2] for n in nodes if n[1] == "gene"}
        pmm0001 = next(
            (props for _, props in gene_nodes.items() if props.get("locus_tag") == "PMM0001"),
            None,
        )
        assert pmm0001 is not None
        assert pmm0001.get("gene_name") == "dnaN"

    def test_real_data_list_fields_are_lists(self, real_data_dir):
        a = CyanorakNcbi(
            data_dir=real_data_dir,
            strain_name="MED4",
            ncbi_accession="GCF_000011465.1",
        )
        a.download_data()
        nodes = a.get_nodes()
        for _, _, props in [n for n in nodes if n[1] == "gene"][:10]:
            if "go_terms" in props:
                assert isinstance(props["go_terms"], list)
            if "kegg_ko" in props:
                assert isinstance(props["kegg_ko"], list)
