"""Unit tests for the MultiUniprot adapter wrapper."""

import os
import tempfile
from unittest.mock import patch, MagicMock

import pytest

from multiomics_kg.adapters.uniprot_adapter import (
    MultiUniprot,
    UniprotAdapter,
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
    """Mock data that would be returned by UniprotAdapter (assembly-based IDs)."""
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
            assert isinstance(adapter, UniprotAdapter)

    def test_organism_ids_extracted_correctly(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert wrapper.organism_ids == [59919, 146891]

    def test_adapters_have_correct_organism(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        assert wrapper.adapters[0].ncbi_taxon_id == 59919
        assert wrapper.adapters[1].ncbi_taxon_id == 146891

    def test_kwargs_passed_to_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv, test_mode=True)
        for adapter in wrapper.adapters:
            assert adapter.test_mode is True

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
        assert wrapper.adapters[0].ncbi_taxon_id == 28108
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
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
            f.write("  # This is a comment with leading whitespace\n")
            f.write("GCF_000011465.1,Pro_MED4,59919,MED4,cache/genomes/MED4/\n")
        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 1

    def test_all_commented_creates_no_adapters(self, temp_dir):
        """When all data rows are commented, no adapters are created."""
        csv_path = os.path.join(temp_dir, "all_commented.csv")
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir\n")
            f.write("# GCF_000011465.1,Pro_MED4,59919,MED4,cache/genomes/MED4/\n")
            f.write("# GCF_000015645.1,Pro_AS9601,146891,AS9601,cache/genomes/AS9601/\n")
        wrapper = MultiUniprot(config_list_file=csv_path)
        assert len(wrapper.adapters) == 0


# ---------------------------------------------------------------------------
# Tests: download_data
# ---------------------------------------------------------------------------


class TestMultiUniprotDownloadData:
    def test_calls_download_on_all_adapters(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        for adapter in wrapper.adapters:
            adapter.download_data = MagicMock()

        wrapper.download_data(cache=True)

        for adapter in wrapper.adapters:
            adapter.download_data.assert_called_once_with(cache=True)

    def test_passes_kwargs_to_download(self, config_csv):
        wrapper = MultiUniprot(config_list_file=config_csv)
        for adapter in wrapper.adapters:
            adapter.download_data = MagicMock()

        wrapper.download_data(cache=False)

        for adapter in wrapper.adapters:
            adapter.download_data.assert_called_once_with(cache=False)


# ---------------------------------------------------------------------------
# Tests: get_nodes
# ---------------------------------------------------------------------------


class TestMultiUniprotGetNodes:
    def test_returns_iterable(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=mock_uniprot_data["protein_nodes"] + mock_uniprot_data["organism_nodes"]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=mock_uniprot_data_2["protein_nodes"] + mock_uniprot_data_2["organism_nodes"]
        )

        nodes = wrapper.get_nodes()
        assert hasattr(nodes, '__iter__')

    def test_aggregates_nodes_from_all_adapters(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_nodes = MagicMock(
            return_value=mock_uniprot_data["protein_nodes"] + mock_uniprot_data["organism_nodes"]
        )
        wrapper.adapters[1].get_nodes = MagicMock(
            return_value=mock_uniprot_data_2["protein_nodes"] + mock_uniprot_data_2["organism_nodes"]
        )

        nodes = list(wrapper.get_nodes())
        # 2 proteins + 1 organism from adapter 1, 1 protein + 1 organism from adapter 2
        assert len(nodes) == 5

    def test_all_nodes_returned_without_deduplication(self, config_csv):
        """get_nodes yields all nodes; no deduplication is performed."""
        wrapper = MultiUniprot(config_list_file=config_csv)

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
                shared_organism,
            ]
        )

        nodes = list(wrapper.get_nodes())
        # 2 proteins + 2 organism entries (no dedup)
        assert len(nodes) == 4

    def test_protein_nodes_aggregated(self, config_csv):
        """Protein nodes from all adapters are all returned."""
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

        nodes = list(wrapper.get_nodes())
        assert len(nodes) == 3

    def test_empty_adapters_returns_empty(self, empty_csv):
        wrapper = MultiUniprot(config_list_file=empty_csv)
        nodes = list(wrapper.get_nodes())
        assert nodes == []


# ---------------------------------------------------------------------------
# Tests: get_edges
# ---------------------------------------------------------------------------


class TestMultiUniprotGetEdges:
    def test_returns_iterable(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_edges = MagicMock(return_value=mock_uniprot_data["edges"])
        wrapper.adapters[1].get_edges = MagicMock(return_value=mock_uniprot_data_2["edges"])

        edges = wrapper.get_edges()
        assert hasattr(edges, '__iter__')

    def test_aggregates_edges_from_all_adapters(self, config_csv, mock_uniprot_data, mock_uniprot_data_2):
        wrapper = MultiUniprot(config_list_file=config_csv)
        wrapper.adapters[0].get_edges = MagicMock(return_value=mock_uniprot_data["edges"])
        wrapper.adapters[1].get_edges = MagicMock(return_value=mock_uniprot_data_2["edges"])

        edges = list(wrapper.get_edges())
        # 2 edges from adapter 1, 1 edge from adapter 2
        assert len(edges) == 3

    def test_empty_adapters_returns_empty(self, empty_csv):
        wrapper = MultiUniprot(config_list_file=empty_csv)
        edges = list(wrapper.get_edges())
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
        assert wrapper.adapters[0].ncbi_taxon_id == 59919  # Prochlorococcus
        assert wrapper.adapters[1].ncbi_taxon_id == 64471  # Synechococcus


# ---------------------------------------------------------------------------
# Regression: edge cardinality, provenance, and property sanitisation
# ---------------------------------------------------------------------------


def _adapter_with(data, refseq_to_strains, locus_to_strains=None, own_proteomes=None):
    """Build a UniprotAdapter with its loaded state injected directly."""
    adapter = UniprotAdapter(
        organism_group="Prochlorococcus",
        ncbi_taxon_id=59919,
        assembly_info=[{"accession": "GCF_000011465.1", "strain_name": "MED4",
                        "ncbi_taxon_id": 59919}],
        data_dirs=["cache/genomes/MED4/"],
    )
    adapter._data = data
    adapter._refseq_to_strains = refseq_to_strains
    adapter._locus_to_strains = locus_to_strains or {}
    adapter._own_proteomes = own_proteomes if own_proteomes is not None \
        else adapter._derive_own_proteomes()
    return adapter


class TestEdgeCardinality:
    """A WP_ accession mapping to several locus tags in ONE assembly yields several
    Gene_encodes_protein edges but only ONE Protein_belongs_to_organism edge.

    80 real proteins hit this (e.g. uniprot:A8WIB5 -> PMM1896 + PMM2004). Keying
    both edges on (locus_tag, ncbi_acc) emitted ~131 duplicate organism edges per
    build, silently absorbed by BioCypher's Deduplicator.
    """

    def test_paralogs_yield_one_organism_edge_and_two_gene_edges(self):
        adapter = _adapter_with(
            {"A8WIB5": {"refseq_ids": ["WP_011132000.1"]}},
            {"WP_011132000.1": [("PMM1896", "GCF_000011465.1"),
                                ("PMM2004", "GCF_000011465.1")]},
        )
        edges = list(adapter.get_edges())
        gene_edges = [e for e in edges if e[3] == "Gene_encodes_protein"]
        org_edges = [e for e in edges if e[3] == "Protein_belongs_to_organism"]

        assert len(gene_edges) == 2
        assert {e[1] for e in gene_edges} == {"ncbigene:PMM1896", "ncbigene:PMM2004"}
        assert len(org_edges) == 1
        assert org_edges[0][2] == "insdc.gcf:GCF_000011465.1"

    def test_distinct_assemblies_still_get_one_organism_edge_each(self):
        adapter = _adapter_with(
            {"P00001": {"refseq_ids": ["WP_1"]}},
            {"WP_1": [("PMM0001", "GCF_000011465.1"), ("A9601_1", "GCF_000015645.1")]},
        )
        org_edges = [e for e in adapter.get_edges()
                     if e[3] == "Protein_belongs_to_organism"]
        assert len(org_edges) == 2
        assert {e[2] for e in org_edges} == {
            "insdc.gcf:GCF_000011465.1", "insdc.gcf:GCF_000015645.1"}

    def test_multiple_refseq_ids_mapping_to_same_gene_are_deduped(self):
        adapter = _adapter_with(
            {"P00002": {"refseq_ids": ["WP_1", "WP_2"]}},
            {"WP_1": [("PMM0001", "GCF_000011465.1")],
             "WP_2": [("PMM0001", "GCF_000011465.1")]},
        )
        edges = list(adapter.get_edges())
        assert len([e for e in edges if e[3] == "Gene_encodes_protein"]) == 1
        assert len([e for e in edges if e[3] == "Protein_belongs_to_organism"]) == 1


class TestProvenance:
    def test_no_version_is_claimed(self):
        """A UniProt release is never captured at download time, so asserting one
        on ~89K edges was a fabrication (the hardcoded '2024_03' was stamped on
        data actually fetched in 2026)."""
        adapter = _adapter_with({}, {})
        prov = adapter._provenance()
        assert "version" not in prov
        assert not hasattr(adapter, "data_version")

    def test_source_and_licence_are_still_reported(self):
        prov = _adapter_with({}, {})._provenance()
        assert prov["source"] == "uniprot"
        assert prov["licence"] == "CC BY 4.0"


class TestNodePropertySanitisation:
    def test_pipes_and_quotes_are_stripped_from_list_and_scalar_props(self):
        adapter = _adapter_with(
            {"P00003": {
                "catalytic_activities": ["Reaction=Xaa-|-Yaa", "Reaction=2'-deoxy"],
                "function_description": "cleaves X-|-Y using Mg'2+",
                "sequence_length": 367,
                "refseq_ids": ["WP_1"],
            }},
            {"WP_1": [("PMM0001", "GCF_000011465.1")]},
        )
        _, _, props = next(iter(adapter.get_nodes()))
        assert props["catalytic_activities"] == ["Reaction=Xaa--Yaa", "Reaction=2^-deoxy"]
        assert props["function_description"] == "cleaves X--Y using Mg^2+"
        # non-string properties pass through untouched
        assert props["sequence_length"] == 367


# ---------------------------------------------------------------------------
# Orphan-protein fix (plans/orphan_proteins.md, 2026-08-27): proteins are kept
# only when linkable to one of OUR assemblies; the organism edge no longer
# depends on the RefSeq join; gene_oln is a fallback gene join.
# ---------------------------------------------------------------------------

ACC = "GCF_000011465.1"


def _by_type(edges, label):
    return [e for e in edges if e[3] == label]


class TestOwnProteomeDerivation:
    def test_own_proteome_is_one_our_assembly_accounts_for(self):
        """An assembly owns a UniProt proteome when its WP_-matched proteins cover
        a MAJORITY of that proteome's entries (and >= MIN_OWN_PROTEOME_SUPPORT).
        The inverse metric (share of OUR matched proteins listing the proteome)
        is fooled by conserved proteins: MIT1002's matched proteins list
        UP000063991 59% of the time, yet 80% of that proteome's WP_ are absent
        from MIT1002's FASTA — it is another isolate's build."""
        data = {}
        for i in range(30):   # matched; list both proteomes
            data[f"P{i:05d}"] = {"refseq_ids": [f"WP_{i}"],
                                 "proteome_ids": ["UP000001026: Chromosome",
                                                  "UP000063991: Chromosome"]}
        for i in range(30, 130):   # unmatched; only the foreign proteome
            data[f"P{i:05d}"] = {"refseq_ids": [f"WP_{i}"],
                                 "proteome_ids": ["UP000063991: Chromosome"]}
        refseq = {f"WP_{i}": [(f"PMM{i:04d}", ACC)] for i in range(30)}
        adapter = _adapter_with(data, refseq)
        assert adapter._own_proteomes == {ACC: {"UP000001026"}}

    def test_no_matched_proteins_means_no_own_proteome(self):
        adapter = _adapter_with({"P1": {"proteome_ids": ["UP1: Chromosome"]}}, {})
        assert adapter._own_proteomes == {}

    def test_proteome_that_is_majority_for_several_assemblies_is_adopted_by_none(self):
        """Species-level taxid 28108: UniProt's UP000095392 is the majority
        proteome of all four A. macleodii assemblies (conserved proteins share one
        WP_), so it identifies none of them — proteome-only proteins there are
        unattributable and must not be pinned on every strain."""
        acc_a, acc_b = "GCF_A", "GCF_B"
        adapter = UniprotAdapter(
            organism_group="Alteromonas", ncbi_taxon_id=28108,
            assembly_info=[{"accession": acc_a, "strain_name": "A", "ncbi_taxon_id": 28108},
                           {"accession": acc_b, "strain_name": "B", "ncbi_taxon_id": 28108}],
            data_dirs=["x/A", "x/B"],
        )
        adapter._data = {f"P{i}": {"refseq_ids": [f"WP_{i}"],
                                   "proteome_ids": ["UP000095392: Unassembled WGS sequence"]}
                         for i in range(40)}
        adapter._refseq_to_strains = {f"WP_{i}": [(f"A_{i}", acc_a), (f"B_{i}", acc_b)]
                                      for i in range(40)}
        adapter._locus_to_strains = {}
        assert adapter._derive_own_proteomes() == {}


class TestProteinKeepRule:
    def test_foreign_proteome_protein_is_dropped_entirely(self):
        """A WP_ that matches nothing AND a proteome that is not ours = another
        isolate under a shared taxid (e.g. 28108). No node, no edges."""
        adapter = _adapter_with(
            {"A0A126Q220": {"refseq_ids": ["WP_061095617.1"],
                            "proteome_ids": ["UP000509458: Chromosome"]}},
            {}, own_proteomes={ACC: {"UP000001026"}},
        )
        assert list(adapter.get_nodes()) == []
        assert list(adapter.get_edges()) == []

    def test_unlinkable_protein_without_proteome_is_dropped(self):
        adapter = _adapter_with({"P1": {"gene_names": ["foo"]}}, {},
                                own_proteomes={ACC: {"UP000001026"}})
        assert list(adapter.get_nodes()) == []

    def test_own_proteome_protein_keeps_node_and_gets_organism_edge_only(self):
        """Same proteome as our matched proteins but a WP_ absent from our
        annotation build: it IS our organism's protein, so it gets the organism
        edge (pre-fe5c2bb behaviour) but no gene edge."""
        adapter = _adapter_with(
            {"A0A399DPR3": {"refseq_ids": ["WP_119361850.1"],
                            "proteome_ids": ["UP000001026: Chromosome"]}},
            {}, own_proteomes={ACC: {"UP000001026"}},
        )
        nodes = list(adapter.get_nodes())
        assert [n[0] for n in nodes] == ["uniprot:A0A399DPR3"]
        edges = list(adapter.get_edges())
        assert _by_type(edges, "Gene_encodes_protein") == []
        org = _by_type(edges, "Protein_belongs_to_organism")
        assert [(e[1], e[2]) for e in org] == [("uniprot:A0A399DPR3", f"insdc.gcf:{ACC}")]

    def test_wp_matched_protein_is_kept_even_without_proteome(self):
        adapter = _adapter_with(
            {"P1": {"refseq_ids": ["WP_1"]}},
            {"WP_1": [("PMM0001", ACC)]}, own_proteomes={},
        )
        assert len(list(adapter.get_nodes())) == 1


class TestLocusTagFallbackJoin:
    def test_gene_oln_links_protein_with_no_refseq_xref(self):
        """7,596 real proteins (e.g. A8WI08 ↔ PMM1805) carry no RefSeq xref but
        name a locus tag we already have in gene_mapping.csv."""
        adapter = _adapter_with(
            {"A8WI08": {"locus_tag": "PMM1805", "gene_names": ["PMM1805"],
                        "proteome_ids": ["UP000001026: Chromosome"]}},
            {}, locus_to_strains={"PMM1805": [("PMM1805", ACC)]},
            own_proteomes={ACC: {"UP000001026"}},
        )
        edges = list(adapter.get_edges())
        gene = _by_type(edges, "Gene_encodes_protein")
        assert [(e[1], e[2]) for e in gene] == [("ncbigene:PMM1805", "uniprot:A8WI08")]
        assert len(_by_type(edges, "Protein_belongs_to_organism")) == 1

    def test_gene_names_tokens_are_tried_when_locus_tag_field_is_absent(self):
        adapter = _adapter_with(
            {"P1": {"gene_names": ["nnrD", "PMM1290"]}},
            {}, locus_to_strains={"PMM1290": [("PMM1290", ACC)]}, own_proteomes={},
        )
        gene = _by_type(list(adapter.get_edges()), "Gene_encodes_protein")
        assert [e[1] for e in gene] == ["ncbigene:PMM1290"]

    def test_locus_join_is_fallback_only_when_refseq_hits(self):
        """When the WP_ join succeeds the locus tag is not consulted, so a
        stale gene_oln cannot add a second, contradictory gene edge."""
        adapter = _adapter_with(
            {"P1": {"refseq_ids": ["WP_1"], "locus_tag": "PMM0002"}},
            {"WP_1": [("PMM0001", ACC)]},
            locus_to_strains={"PMM0002": [("PMM0002", ACC)]}, own_proteomes={},
        )
        gene = _by_type(list(adapter.get_edges()), "Gene_encodes_protein")
        assert [e[1] for e in gene] == ["ncbigene:PMM0001"]


class TestGeneMappingLocusIndex:
    def test_load_gene_mapping_indexes_current_and_old_locus_tags(self, temp_dir):
        d = os.path.join(temp_dir, "MED4"); os.makedirs(d)
        with open(os.path.join(d, "gene_mapping.csv"), "w") as f:
            f.write("locus_tag,locus_tag_ncbi,old_locus_tags,locus_tag_cyanorak,protein_id\n")
            f.write("PMM0001,TX50_RS00020,PMM0001_old,CK_Pro_MED4_00001,WP_1\n")
            f.write("PMM0002,TX50_RS00025,,,\n")
        adapter = UniprotAdapter("Prochlorococcus", 59919,
                                 [{"accession": ACC, "strain_name": "MED4", "ncbi_taxon_id": 59919}],
                                 [d])
        refseq, locus = adapter._load_gene_mapping()
        assert refseq == {"WP_1": [("PMM0001", ACC)]}
        assert locus["PMM0001"] == [("PMM0001", ACC)]
        assert locus["TX50_RS00020"] == [("PMM0001", ACC)]
        assert locus["PMM0001_old"] == [("PMM0001", ACC)]
        assert locus["CK_Pro_MED4_00001"] == [("PMM0001", ACC)]
        assert locus["TX50_RS00025"] == [("PMM0002", ACC)]
