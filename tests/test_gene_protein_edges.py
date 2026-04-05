"""Unit tests for gene-to-protein edge creation via RefSeq protein_id matching."""

import os
import tempfile

import pytest

from multiomics_kg.adapters.uniprot_adapter import (
    MultiUniprot,
    UniprotAdapter,
)

# Dummy assembly info used when we only care about gene_mapping loading
_MED4_INFO = [{"accession": "GCF_000011465.1", "strain_name": "MED4", "ncbi_taxon_id": 59919}]
_MIT1002_INFO = [{"accession": "GCF_901457835.2", "strain_name": "MIT1002", "ncbi_taxon_id": 28108}]
_EZ55_INFO = [{"accession": "GCF_901457815.2", "strain_name": "EZ55", "ncbi_taxon_id": 28108}]
_ALT_INFO = _MIT1002_INFO + _EZ55_INFO


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


def _make_adapter(data_dirs, assembly_info, refseq_data, organism=59919):
    """Create a UniprotAdapter with mocked protein data for testing edge creation.

    Args:
        data_dirs: List of data_dir paths (for gene_mapping.csv lookup).
        assembly_info: Parallel list of assembly info dicts.
        refseq_data: Dict mapping uniprot accession -> RefSeq WP_ id (str or None/empty).
        organism: NCBI taxid.

    Returns:
        UniprotAdapter with _data and _refseq_to_strains populated (no API download).
    """
    adapter = UniprotAdapter(
        organism_group="Prochlorococcus",
        ncbi_taxon_id=organism,
        assembly_info=assembly_info,
        data_dirs=data_dirs,
    )
    # Build _data: {uniprot_acc: {refseq_ids: [WP_...]}}
    adapter._data = {}
    for acc, rs in refseq_data.items():
        if rs:
            adapter._data[acc] = {"refseq_ids": [rs]}
        else:
            adapter._data[acc] = {"refseq_ids": []}
    # Build _refseq_to_strains from gene_mapping.csv files
    adapter._refseq_to_strains = adapter._load_gene_mapping()
    return adapter


def _gene_edges(adapter):
    """Return only Gene_encodes_protein edges as a list."""
    return [e for e in adapter.get_edges() if e[3] == "Gene_encodes_protein"]


# ---------------------------------------------------------------------------
# Tests: _load_gene_mapping
# ---------------------------------------------------------------------------


class TestLoadGeneMapping:
    def test_builds_map_from_single_file(self, gene_mapping_med4):
        adapter = UniprotAdapter(
            organism_group="Prochlorococcus",
            ncbi_taxon_id=59919,
            assembly_info=_MED4_INFO,
            data_dirs=[gene_mapping_med4],
        )
        refseq_map = adapter._load_gene_mapping()
        assert "WP_011131639.1" in refseq_map
        locus_tags = [lt for lt, _ in refseq_map["WP_011131639.1"]]
        assert locus_tags == ["PMM0001"]
        assert len(refseq_map) == 3

    def test_builds_map_from_multiple_files(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        adapter = UniprotAdapter(
            organism_group="Alteromonas",
            ncbi_taxon_id=28108,
            assembly_info=_ALT_INFO,
            data_dirs=[gene_mapping_mit1002, gene_mapping_ez55],
        )
        refseq_map = adapter._load_gene_mapping()
        # MIT1002 has 3 genes, EZ55 has 2, but WP_SHARED_001.1 is shared
        assert len(refseq_map) == 4  # 4 unique protein_ids

    def test_shared_refseq_maps_to_multiple_locus_tags(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        adapter = UniprotAdapter(
            organism_group="Alteromonas",
            ncbi_taxon_id=28108,
            assembly_info=_ALT_INFO,
            data_dirs=[gene_mapping_mit1002, gene_mapping_ez55],
        )
        refseq_map = adapter._load_gene_mapping()
        shared_locus_tags = [lt for lt, _ in refseq_map["WP_SHARED_001.1"]]
        assert len(shared_locus_tags) == 2
        assert "MIT1002_00004" in shared_locus_tags
        assert "EZ55_00010" in shared_locus_tags

    def test_empty_data_dirs_returns_empty_map(self):
        adapter = UniprotAdapter(
            organism_group="Prochlorococcus",
            ncbi_taxon_id=59919,
            assembly_info=[],
            data_dirs=[],
        )
        refseq_map = adapter._load_gene_mapping()
        assert refseq_map == {}

    def test_missing_file_logs_warning(self, temp_dir):
        adapter = UniprotAdapter(
            organism_group="Prochlorococcus",
            ncbi_taxon_id=59919,
            assembly_info=_MED4_INFO,
            data_dirs=[os.path.join(temp_dir, "nonexistent")],
        )
        refseq_map = adapter._load_gene_mapping()
        assert refseq_map == {}


# ---------------------------------------------------------------------------
# Tests: Gene-to-protein edge creation
# ---------------------------------------------------------------------------


class TestGeneToProteinEdges:
    def test_basic_matching(self, gene_mapping_med4):
        """Protein with matching RefSeq ID produces an edge."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_med4],
            assembly_info=_MED4_INFO,
            refseq_data={"P12345": "WP_011131639.1"},  # matches PMM0001
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 1
        _, source, target, label, props = edges[0]
        assert source == "ncbigene:PMM0001"
        assert target == "uniprot:P12345"
        assert label == "Gene_encodes_protein"

    def test_no_match(self, gene_mapping_med4):
        """Protein with non-matching RefSeq ID produces no edge."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_med4],
            assembly_info=_MED4_INFO,
            refseq_data={"P99999": "WP_000000000.1"},  # not in gene_mapping
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 0

    def test_missing_refseq(self, gene_mapping_med4):
        """Protein with empty RefSeq produces no edge."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_med4],
            assembly_info=_MED4_INFO,
            refseq_data={"P12345": "", "P12346": None},
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 0

    def test_multiple_proteins_matched(self, gene_mapping_med4):
        """Multiple proteins each matching different genes."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_med4],
            assembly_info=_MED4_INFO,
            refseq_data={
                "P12345": "WP_011131639.1",  # PMM0001
                "P12346": "WP_011131640.1",  # PMM0002
                "P99999": "WP_000000000.1",  # no match
            },
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 2
        sources = {e[1] for e in edges}
        assert sources == {"ncbigene:PMM0001", "ncbigene:PMM0002"}

    def test_multi_strain_produces_multiple_edges(
        self, gene_mapping_mit1002, gene_mapping_ez55
    ):
        """Same WP_ ID in multiple strains creates one edge per gene."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_mit1002, gene_mapping_ez55],
            assembly_info=_ALT_INFO,
            refseq_data={"A0A_SHARED": "WP_SHARED_001.1"},
            organism=28108,
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 2
        sources = {e[1] for e in edges}
        assert sources == {"ncbigene:MIT1002_00004", "ncbigene:EZ55_00010"}

    def test_no_data_dirs_skips_gracefully(self):
        """No data_dirs means no gene mapping, so no Gene_encodes_protein edges."""
        adapter = _make_adapter(
            data_dirs=[],
            assembly_info=[],
            refseq_data={"P12345": "WP_011131639.1"},
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 0

    def test_edge_properties(self, gene_mapping_med4):
        """Edges have correct provenance properties."""
        adapter = _make_adapter(
            data_dirs=[gene_mapping_med4],
            assembly_info=_MED4_INFO,
            refseq_data={"P12345": "WP_011131639.1"},
        )
        edges = _gene_edges(adapter)
        assert len(edges) == 1
        props = edges[0][4]
        assert props["source"] == "uniprot"
        assert "licence" in props


# ---------------------------------------------------------------------------
# Tests: MultiUniprot data_dirs passthrough
# ---------------------------------------------------------------------------


class TestMultiUniprotDataDirs:
    def test_data_dirs_passed_to_adapters(self, temp_dir):
        """MultiUniprot collects data_dir values from the genomes CSV."""
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
        assert wrapper.adapters[0].data_dirs == ["cache/data/Pro/MED4/"]

    def test_same_taxid_collects_all_data_dirs(self, temp_dir):
        """Multiple strains with same taxid get all data_dirs in one adapter."""
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
        dirs = wrapper.adapters[0].data_dirs
        assert len(dirs) == 2
        assert "cache/MIT1002/" in dirs
        assert "cache/EZ55/" in dirs

    def test_different_taxids_get_separate_adapters(self, temp_dir):
        """Each unique taxid gets its own adapter with only its data_dirs."""
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
        assert wrapper.adapters[0].data_dirs == ["cache/MED4/"]
        assert wrapper.adapters[1].data_dirs == ["cache/MIT1002/"]
