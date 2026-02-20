"""Unit tests for gene-to-protein edge creation via RefSeq protein_id matching."""

import os
import tempfile

import pytest

from multiomics_kg.adapters.uniprot_adapter import (
    MultiUniprot,
    Uniprot,
    UniprotEdgeType,
    UniprotNodeField,
    UniprotNodeType,
    UniprotIDField,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def gene_mapping_med4(temp_dir):
    """Create a gene_mapping.csv for MED4 with a few genes."""
    data_dir = os.path.join(temp_dir, "MED4")
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, "gene_mapping.csv")
    with open(path, "w") as f:
        f.write("gene_names,locus_tag,protein_id,locus_tag_ncbi\n")
        f.write("dnaN,PMM0001,WP_011131639.1,TX50_RS00020\n")
        f.write(",PMM0002,WP_011131640.1,TX50_RS00025\n")
        f.write("rnpA,PMM0003,WP_011131641.1,TX50_RS00030\n")
    return data_dir


@pytest.fixture
def gene_mapping_mit1002(temp_dir):
    """Create a gene_mapping.csv for MIT1002 (Alteromonas)."""
    data_dir = os.path.join(temp_dir, "MIT1002")
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, "gene_mapping.csv")
    with open(path, "w") as f:
        f.write("gene_names,locus_tag,protein_id,locus_tag_ncbi\n")
        f.write("dnaN,MIT1002_00002,WP_049586368.1,ALT831_RS00010\n")
        f.write("recF,MIT1002_00003,WP_049586367.1,ALT831_RS00015\n")
        # Shared WP_ ID with EZ55 (same protein across strains)
        f.write("gyrB,MIT1002_00004,WP_SHARED_001.1,ALT831_RS00020\n")
    return data_dir


@pytest.fixture
def gene_mapping_ez55(temp_dir):
    """Create a gene_mapping.csv for EZ55 (Alteromonas, same taxid as MIT1002)."""
    data_dir = os.path.join(temp_dir, "EZ55")
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, "gene_mapping.csv")
    with open(path, "w") as f:
        f.write("gene_names,locus_tag,protein_id,locus_tag_ncbi\n")
        f.write("gyrB,EZ55_00010,WP_SHARED_001.1,BSR22_RS00050\n")
        f.write("dnaA,EZ55_00011,WP_061093873.1,BSR22_RS00055\n")
    return data_dir


def _make_adapter_with_mock_data(
    gene_mapping_paths, refseq_data, organism=59919
):
    """Create a Uniprot adapter with mocked data for testing edge creation.

    Args:
        gene_mapping_paths: List of gene_mapping.csv file paths.
        refseq_data: Dict mapping uniprot accession -> RefSeq protein ID.
        organism: NCBI taxid.

    Returns:
        Configured Uniprot adapter with mock data loaded.
    """
    adapter = Uniprot(
        organism=organism,
        edge_types=[UniprotEdgeType.GENE_TO_PROTEIN],
        node_fields=[UniprotNodeField.xref_refseq],
        gene_mapping_paths=gene_mapping_paths,
    )
    # Inject mock data directly (skip API download)
    adapter.data = {
        UniprotNodeField.xref_refseq.value: refseq_data,
    }
    return adapter


# ---------------------------------------------------------------------------
# Tests: _build_refseq_to_gene_map
# ---------------------------------------------------------------------------


class TestBuildRefseqToGeneMap:
    def test_builds_map_from_single_file(self, gene_mapping_med4):
        adapter = Uniprot(
            organism=59919,
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
        )
        refseq_map = adapter._build_refseq_to_gene_map()
        assert "WP_011131639.1" in refseq_map
        assert refseq_map["WP_011131639.1"] == ["PMM0001"]
        assert len(refseq_map) == 3

    def test_builds_map_from_multiple_files(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        adapter = Uniprot(
            organism=28108,
            gene_mapping_paths=[
                os.path.join(gene_mapping_mit1002, "gene_mapping.csv"),
                os.path.join(gene_mapping_ez55, "gene_mapping.csv"),
            ],
        )
        refseq_map = adapter._build_refseq_to_gene_map()
        # MIT1002 has 3 genes, EZ55 has 2, but WP_SHARED_001.1 is shared
        assert len(refseq_map) == 4  # 4 unique protein_ids

    def test_shared_refseq_maps_to_multiple_locus_tags(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        adapter = Uniprot(
            organism=28108,
            gene_mapping_paths=[
                os.path.join(gene_mapping_mit1002, "gene_mapping.csv"),
                os.path.join(gene_mapping_ez55, "gene_mapping.csv"),
            ],
        )
        refseq_map = adapter._build_refseq_to_gene_map()
        shared = refseq_map["WP_SHARED_001.1"]
        assert len(shared) == 2
        assert "MIT1002_00004" in shared
        assert "EZ55_00010" in shared

    def test_empty_paths_returns_empty_map(self):
        adapter = Uniprot(organism=59919, gene_mapping_paths=[])
        refseq_map = adapter._build_refseq_to_gene_map()
        assert refseq_map == {}

    def test_missing_file_logs_warning(self, temp_dir):
        adapter = Uniprot(
            organism=59919,
            gene_mapping_paths=[os.path.join(temp_dir, "nonexistent.csv")],
        )
        refseq_map = adapter._build_refseq_to_gene_map()
        assert refseq_map == {}

    def test_map_is_cached(self, gene_mapping_med4):
        adapter = Uniprot(
            organism=59919,
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
        )
        map1 = adapter._build_refseq_to_gene_map()
        map2 = adapter._build_refseq_to_gene_map()
        assert map1 is map2  # Same object, not rebuilt


# ---------------------------------------------------------------------------
# Tests: Gene-to-protein edge creation
# ---------------------------------------------------------------------------


class TestGeneToProteinEdges:
    def test_basic_matching(self, gene_mapping_med4):
        """Protein with matching RefSeq ID produces an edge."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
            refseq_data={
                "P12345": "WP_011131639.1",  # matches PMM0001
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 1
        _, source, target, label, props = edges[0]
        assert source == "uniprot:P12345"
        assert target == "ncbigene:PMM0001"
        assert label == "Gene_encodes_protein"

    def test_no_match(self, gene_mapping_med4):
        """Protein with non-matching RefSeq ID produces no edge."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
            refseq_data={
                "P99999": "WP_000000000.1",  # not in gene_mapping
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 0

    def test_missing_refseq(self, gene_mapping_med4):
        """Protein with empty RefSeq produces no edge."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
            refseq_data={
                "P12345": "",
                "P12346": None,
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 0

    def test_multiple_proteins_matched(self, gene_mapping_med4):
        """Multiple proteins each matching different genes."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
            refseq_data={
                "P12345": "WP_011131639.1",  # PMM0001
                "P12346": "WP_011131640.1",  # PMM0002
                "P99999": "WP_000000000.1",  # no match
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 2
        targets = {e[2] for e in edges}
        assert targets == {"ncbigene:PMM0001", "ncbigene:PMM0002"}

    def test_multi_strain_produces_multiple_edges(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        """Same WP_ ID in multiple strains creates one edge per gene."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_mit1002, "gene_mapping.csv"),
                os.path.join(gene_mapping_ez55, "gene_mapping.csv"),
            ],
            refseq_data={
                "A0A_SHARED": "WP_SHARED_001.1",
            },
            organism=28108,
        )
        edges = adapter.get_edges()
        assert len(edges) == 2
        targets = {e[2] for e in edges}
        assert targets == {"ncbigene:MIT1002_00004", "ncbigene:EZ55_00010"}

    def test_no_gene_mapping_paths_skips_gracefully(self):
        """GENE_TO_PROTEIN enabled but no gene_mapping_paths produces no edges."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[],
            refseq_data={
                "P12345": "WP_011131639.1",
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 0

    def test_edge_properties(self, gene_mapping_med4):
        """Edges have correct provenance properties."""
        adapter = _make_adapter_with_mock_data(
            gene_mapping_paths=[
                os.path.join(gene_mapping_med4, "gene_mapping.csv")
            ],
            refseq_data={
                "P12345": "WP_011131639.1",
            },
        )
        edges = adapter.get_edges()
        assert len(edges) == 1
        props = edges[0][4]
        assert props["source"] == "uniprot"
        assert "licence" in props


# ---------------------------------------------------------------------------
# Tests: MultiUniprot gene_mapping_paths passthrough
# ---------------------------------------------------------------------------


class TestMultiUniprotGeneMappingPaths:
    def test_gene_mapping_paths_passed_to_adapters(self, temp_dir):
        """MultiUniprot constructs gene_mapping paths from data_dir column."""
        csv_path = os.path.join(temp_dir, "genomes.csv")
        with open(csv_path, "w") as f:
            f.write(
                "ncbi_accession,cyanorak_organism,ncbi_taxon_id,"
                "strain_name,data_dir\n"
            )
            f.write(
                "GCF_000011465.1,Pro_MED4,59919,MED4,"
                "cache/data/Pro/MED4/\n"
            )

        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 1
        expected = ["cache/data/Pro/MED4/gene_mapping.csv"]
        assert wrapper.adapters[0].gene_mapping_paths == expected

    def test_same_taxid_collects_all_data_dirs(self, temp_dir):
        """Multiple strains with same taxid get all gene_mapping paths."""
        csv_path = os.path.join(temp_dir, "alt.csv")
        with open(csv_path, "w") as f:
            f.write(
                "ncbi_accession,cyanorak_organism,ncbi_taxon_id,"
                "strain_name,data_dir\n"
            )
            f.write("GCF_901457835.2,,28108,MIT1002,cache/MIT1002/\n")
            f.write("GCF_901457815.2,,28108,EZ55,cache/EZ55/\n")

        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 1
        paths = wrapper.adapters[0].gene_mapping_paths
        assert len(paths) == 2
        assert "cache/MIT1002/gene_mapping.csv" in paths
        assert "cache/EZ55/gene_mapping.csv" in paths

    def test_different_taxids_get_separate_paths(self, temp_dir):
        """Each unique taxid adapter gets only its own data_dir paths."""
        csv_path = os.path.join(temp_dir, "mixed.csv")
        with open(csv_path, "w") as f:
            f.write(
                "ncbi_accession,cyanorak_organism,ncbi_taxon_id,"
                "strain_name,data_dir\n"
            )
            f.write(
                "GCF_000011465.1,Pro_MED4,59919,MED4,"
                "cache/MED4/\n"
            )
            f.write("GCF_901457835.2,,28108,MIT1002,cache/MIT1002/\n")

        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 2
        # MED4 adapter
        assert wrapper.adapters[0].gene_mapping_paths == [
            "cache/MED4/gene_mapping.csv"
        ]
        # Alteromonas adapter
        assert wrapper.adapters[1].gene_mapping_paths == [
            "cache/MIT1002/gene_mapping.csv"
        ]
