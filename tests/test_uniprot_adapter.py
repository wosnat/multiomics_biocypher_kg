"""Unit tests for the MultiUniprot adapter wrapper."""

import os
import tempfile
from unittest.mock import patch, MagicMock

import pytest

from multiomics_kg.adapters.uniprot_adapter import (
    MultiUniprot,
    Uniprot,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def config_csv(temp_dir):
    """Create a CSV config file with multiple organisms."""
    csv_path = os.path.join(temp_dir, "cyanobacteria_genomes.csv")
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
        f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/genomes/MED4/\n")
        f.write("GCF_000015645.1,Pro_AS9601,146891,AS9601,cache/genomes/AS9601/\n")
    return csv_path


@pytest.fixture
def single_organism_csv(temp_dir):
    """Create a CSV config file with a single organism."""
    csv_path = os.path.join(temp_dir, "single_genome.csv")
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
        f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/genomes/MED4/\n")
    return csv_path


@pytest.fixture
def empty_csv(temp_dir):
    """Create an empty CSV config file (header only)."""
    csv_path = os.path.join(temp_dir, "empty.csv")
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
    return csv_path


@pytest.fixture
def csv_with_comments(temp_dir):
    """Create a CSV config file with commented lines."""
    csv_path = os.path.join(temp_dir, "commented_genomes.csv")
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
        f.write("# This is a comment\n")
        f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/genomes/MED4/\n")
        f.write("# GCF_000015645.1,Pro_AS9601,146891,AS9601,cache/genomes/AS9601/\n")
    return csv_path


@pytest.fixture
def mock_uniprot_data():
    """Mock data that would be returned by Uniprot adapter (assembly-based IDs)."""
    return {
        "protein_nodes": [
            ("uniprot:P12345", "protein", {"length": 100, "organism_id": 59919}),
            ("uniprot:P12346", "protein", {"length": 200, "organism_id": 59919}),
        ],
        "organism_nodes": [
            ("insdc.gcf:GCF_000011465.1", "organism", {
                "organism_name": "Prochlorococcus marinus",
                "strain_name": "MED4",
                "ncbi_taxon_id": 59919,
            }),
        ],
        "edges": [
            (None, "uniprot:P12345", "insdc.gcf:GCF_000011465.1", "Protein_belongs_to_organism", {}),
            (None, "uniprot:P12346", "insdc.gcf:GCF_000011465.1", "Protein_belongs_to_organism", {}),
        ],
    }


@pytest.fixture
def mock_uniprot_data_2():
    """Mock data for a second organism (assembly-based IDs)."""
    return {
        "protein_nodes": [
            ("uniprot:Q99999", "protein", {"length": 150, "organism_id": 146891}),
        ],
        "organism_nodes": [
            ("insdc.gcf:GCF_000015645.1", "organism", {
                "organism_name": "Prochlorococcus marinus AS9601",
                "strain_name": "AS9601",
                "ncbi_taxon_id": 146891,
            }),
        ],
        "edges": [
            (None, "uniprot:Q99999", "insdc.gcf:GCF_000015645.1", "Protein_belongs_to_organism", {}),
        ],
    }


# ---------------------------------------------------------------------------
# Tests: MultiUniprot Construction
# ---------------------------------------------------------------------------


class TestMultiUniprotConstruction:
    def test_loads_correct_number_of_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert len(wrapper.adapters) == 2

    def test_single_organism_config(self, single_organism_csv):
        wrapper = MultiUniprot(config_list_file=single_organism_csv)
        assert len(wrapper.adapters) == 1

    def test_empty_csv_creates_no_adapters(self, empty_csv):
        wrapper = MultiUniprot(config_list_file=empty_csv)
        assert len(wrapper.adapters) == 0

    def test_adapters_are_uniprot_instances(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        for adapter in wrapper.adapters:
            assert isinstance(adapter, Uniprot)

    def test_organism_ids_extracted_correctly(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert wrapper.organism_ids == [59919, 146891]

    def test_adapters_have_correct_organism(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert wrapper.adapters[0].organism == 59919
        assert wrapper.adapters[1].organism == 146891

    def test_kwargs_passed_to_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv, test_mode=True)
        for adapter in wrapper.adapters:
            assert adapter.test_mode is True

    def test_rev_kwarg_passed_to_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv, rev=False)
        for adapter in wrapper.adapters:
            assert adapter.rev is False

    def test_add_prefix_kwarg_passed_to_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv, add_prefix=False)
        for adapter in wrapper.adapters:
            assert adapter.add_prefix is False

    def test_assembly_info_passed_to_adapters(self, config_csv):
        """Each adapter should receive assembly_info with accession and strain_name."""
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert len(wrapper.adapters[0].assembly_info) == 1
        assert wrapper.adapters[0].assembly_info[0]['accession'] == 'GCF_000011465.1'
        assert wrapper.adapters[0].assembly_info[0]['strain_name'] == 'MED4'

    def test_same_taxid_dedup(self, temp_dir):
        """Two rows with same taxid create one adapter with both assemblies."""
        csv_path = os.path.join(temp_dir, "dedup.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
            f.write("GCF_001077695.1,,28108,MIT1002,cache/MIT1002/\n")
            f.write("GCF_901457815.2,,28108,EZ55,cache/EZ55/\n")

        wrapper = MultiUniprot(config_list_file=csv_path)
        # One adapter (deduped by taxid)
        assert len(wrapper.adapters) == 1
        assert wrapper.organism_ids == [28108]
        assert wrapper.adapters[0].organism == 28108
        # But assembly_info has both assemblies
        assert len(wrapper.adapters[0].assembly_info) == 2
        accessions = [info['accession'] for info in wrapper.adapters[0].assembly_info]
        assert 'GCF_001077695.1' in accessions
        assert 'GCF_901457815.2' in accessions


# ---------------------------------------------------------------------------
# Tests: Comment line skipping
# ---------------------------------------------------------------------------


class TestMultiUniprotCommentSkipping:
    def test_skips_comment_lines(self, csv_with_comments):
        wrapper = MultiUniprot(config_list_file=csv_with_comments)
        # Only one organism should be loaded (the second is commented out)
        assert len(wrapper.adapters) == 1
        assert wrapper.organism_ids == [59919]

    def test_whitespace_before_comment_skipped(self, temp_dir):
        """Lines with whitespace before # are still treated as comments."""
        csv_path = os.path.join(temp_dir, "whitespace_comment.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,data_dir\n")
            f.write("  # This is a comment with leading whitespace\n")
            f.write("GCF_000011465.1,Pro_MED4,59919,cache/genomes/MED4/\n")
        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 1

    def test_all_commented_creates_no_adapters(self, temp_dir):
        """When all data rows are commented, no adapters are created."""
        csv_path = os.path.join(temp_dir, "all_commented.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,data_dir\n")
            f.write("# GCF_000011465.1,Pro_MED4,59919,cache/genomes/MED4/\n")
            f.write("# GCF_000015645.1,Pro_AS9601,146891,cache/genomes/AS9601/\n")
        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 0


# ---------------------------------------------------------------------------
# Tests: download_uniprot_data
# ---------------------------------------------------------------------------


class TestMultiUniprotDownloadData:
    def test_calls_download_on_all_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        # Mock the download_uniprot_data method on each adapter
        for adapter in wrapper.adapters:
            adapter.download_uniprot_data = MagicMock()

        wrapper.download_uniprot_data(cache=True)

        for adapter in wrapper.adapters:
            adapter.download_uniprot_data.assert_called_once_with(cache=True)

    def test_passes_kwargs_to_download(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        for adapter in wrapper.adapters:
            adapter.download_uniprot_data = MagicMock()

        wrapper.download_uniprot_data(cache=False, debug=True)

        for adapter in wrapper.adapters:
            adapter.download_uniprot_data.assert_called_once_with(cache=False, debug=True)


# ---------------------------------------------------------------------------
# Tests: get_nodes
# ---------------------------------------------------------------------------


class TestMultiUniprotGetNodes:
    def test_returns_list(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        # Mock get_nodes on each adapter
        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=mock_uniprot_data["protein_nodes"] + mock_uniprot_data["organism_nodes"]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=mock_uniprot_data_2["protein_nodes"] + mock_uniprot_data_2["organism_nodes"]
        )

        nodes = wrapper.get_nodes()
        assert isinstance(nodes, list)

    def test_aggregates_nodes_from_all_adapters(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=mock_uniprot_data["protein_nodes"] + mock_uniprot_data["organism_nodes"]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=mock_uniprot_data_2["protein_nodes"] + mock_uniprot_data_2["organism_nodes"]
        )

        nodes = wrapper.get_nodes()
        # 2 proteins + 1 organism from adapter 1, 1 protein + 1 organism from adapter 2
        assert len(nodes) == 5

    def test_deduplicates_organism_nodes(self, config_csv):
        """Test that duplicate organism nodes are removed."""
        wrapper = MultiUniprot(config_list_file=config_csv)

        # Both adapters return the same organism node
        shared_organism = ("insdc.gcf:GCF_000011465.1", "organism", {"organism_name": "Prochlorococcus"})
        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=[
                ("uniprot:P12345", "protein", {}),
                shared_organism,
            ]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=[
                ("uniprot:Q99999", "protein", {}),
                shared_organism,  # Same organism node
            ]
        )

        nodes = wrapper.get_nodes()
        # 2 proteins + 1 organism (deduplicated)
        assert len(nodes) == 3

        # Count organism nodes
        organism_count = sum(1 for _, label, _ in nodes if label == "organism")
        assert organism_count == 1

    def test_protein_nodes_not_deduplicated(self, config_csv):
        """Test that protein nodes are NOT deduplicated (they should be unique anyway)."""
        wrapper = MultiUniprot(config_list_file=config_csv)

        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=[
                ("uniprot:P12345", "protein", {}),
                ("uniprot:P12346", "protein", {}),
            ]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=[
                ("uniprot:Q99999", "protein", {}),
            ]
        )

        nodes = wrapper.get_nodes()
        assert len(nodes) == 3

    def test_passes_kwargs_to_get_nodes(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_nodes = MagicMock(return_value=[])
        wrapper.adapters[1].get_nodes = MagicMock(return_value=[])

        wrapper.get_nodes(protein_label="custom_protein")

        for adapter in wrapper.adapters:
            adapter.get_nodes.assert_called_once_with(protein_label="custom_protein")

    def test_empty_adapters_returns_empty_list(self, empty_csv):
        wrapper = MultiUniprot(config_list_file=empty_csv)
        nodes = wrapper.get_nodes()
        assert nodes == []


# ---------------------------------------------------------------------------
# Tests: get_edges
# ---------------------------------------------------------------------------


class TestMultiUniprotGetEdges:
    def test_returns_list(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_edges = MagicMock(return_value=mock_uniprot_data["edges"])
        wrapper.adapters[1].get_edges = MagicMock(return_value=mock_uniprot_data_2["edges"])

        edges = wrapper.get_edges()
        assert isinstance(edges, list)

    def test_aggregates_edges_from_all_adapters(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_edges = MagicMock(return_value=mock_uniprot_data["edges"])
        wrapper.adapters[1].get_edges = MagicMock(return_value=mock_uniprot_data_2["edges"])

        edges = wrapper.get_edges()
        # 2 edges from adapter 1, 1 edge from adapter 2
        assert len(edges) == 3

    def test_passes_kwargs_to_get_edges(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_edges = MagicMock(return_value=[])
        wrapper.adapters[1].get_edges = MagicMock(return_value=[])

        wrapper.get_edges(some_param="value")

        for adapter in wrapper.adapters:
            adapter.get_edges.assert_called_once_with(some_param="value")

    def test_empty_adapters_returns_empty_list(self, empty_csv):
        wrapper = MultiUniprot(config_list_file=empty_csv)
        edges = wrapper.get_edges()
        assert edges == []


# ---------------------------------------------------------------------------
# Tests: Integration with real config file format
# ---------------------------------------------------------------------------


class TestMultiUniprotRealConfigFormat:
    def test_parses_real_config_format(self, temp_dir):
        """Test parsing the actual config file format used in the project."""
        csv_path = os.path.join(temp_dir, "cyanobacteria_genomes.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
            f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/data/Prochlorococcus/genomes/MED4/\n")
            f.write("GCF_000015645.1,Pro_AS9601,146891,AS9601,cache/data/Prochlorococcus/genomes/AS9601/\n")
            f.write("GCF_000015965.1,Pro_MIT9301,167546,MIT9301,cache/data/Prochlorococcus/genomes/MIT9301/\n")
            f.write("# GCF_001989415.1,Pro_RSP50,1924285,RSP50,cache/data/Prochlorococcus/genomes/RSP50/\n")
            f.write("GCF_000014585.1,Syn_CC9311,64471,CC9311,cache/data/Synechococcus/genomes/CC9311/\n")

        wrapper = MultiUniprot(config_list_file=csv_path)

        # Should have 4 adapters (one commented out)
        assert len(wrapper.adapters) == 4
        assert wrapper.organism_ids == [59919, 146891, 167546, 64471]

    def test_handles_mixed_prochlorococcus_synechococcus(self, temp_dir):
        """Test that both Prochlorococcus and Synechococcus organisms are handled."""
        csv_path = os.path.join(temp_dir, "mixed.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
            f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/Pro/MED4/\n")
            f.write("GCF_000014585.1,Syn_CC9311,64471,CC9311,cache/Syn/CC9311/\n")

        wrapper = MultiUniprot(config_list_file=csv_path)

        assert len(wrapper.adapters) == 2
        assert wrapper.adapters[0].organism == 59919  # Prochlorococcus
        assert wrapper.adapters[1].organism == 64471  # Synechococcus
