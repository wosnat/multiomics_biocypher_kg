"""
Unit tests for the Pfam section of
multiomics_kg/adapters/functional_annotation_adapter.py.

All tests are offline -- no downloads, no live Neo4j required.
"""

import json
import pytest
from pathlib import Path

from multiomics_kg.utils.pfam_utils import PfamData, PfamEntry
from multiomics_kg.adapters.functional_annotation_adapter import (
    PfamAnnotationAdapter,
    MultiPfamAnnotationAdapter,
    _pfam_node_id,
    _pfam_clan_node_id,
)


# ---------------------------------------------------------------------------
# Minimal test data
# ---------------------------------------------------------------------------

# PfamData with 4 entries: 3 with clan (CL0060), 1 without clan
MINI_PFAM_DATA = PfamData(
    by_accession={
        "PF00712": PfamEntry(
            accession="PF00712",
            shortname="DNA_pol3_beta",
            description="DNA polymerase III beta subunit, N-terminal domain",
            clan_accession="CL0060",
            clan_name="DNA_clamp",
        ),
        "PF02768": PfamEntry(
            accession="PF02768",
            shortname="DNA_pol3_beta_2",
            description="DNA polymerase III beta subunit, central domain",
            clan_accession="CL0060",
            clan_name="DNA_clamp",
        ),
        "PF00069": PfamEntry(
            accession="PF00069",
            shortname="Pkinase",
            description="Protein kinase domain",
            clan_accession="",
            clan_name="",
        ),
        "PF99999": PfamEntry(
            accession="PF99999",
            shortname="Test_entry",
            description="Entry with 'quotes' and | pipes",
            clan_accession="CL9999",
            clan_name="Test_clan",
        ),
    },
    by_shortname={
        "DNA_pol3_beta": "PF00712",
        "DNA_pol3_beta_2": "PF02768",
        "Pkinase": "PF00069",
        "Test_entry": "PF99999",
    },
    clans={
        "CL0060": "DNA_clamp",
        "CL9999": "Test_clan",
    },
)


# Gene annotations referencing these Pfam IDs
MINI_GENE_DATA = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "pfam_ids": ["PF00712", "PF02768"],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "pfam_ids": ["PF00069"],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "pfam_ids": [],
    },
    "PMM0004": {
        "locus_tag": "PMM0004",
        "pfam_ids": None,
    },
    "PMM0005": {
        "locus_tag": "PMM0005",
        # No pfam_ids key at all
    },
    "PMM0006": {
        "locus_tag": "PMM0006",
        "pfam_ids": ["PF99999"],
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def genome_dir(tmp_path):
    """Write MINI_GENE_DATA as gene_annotations_merged.json; return dir."""
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_GENE_DATA), encoding="utf-8"
    )
    return d


@pytest.fixture
def pfam_adapter(genome_dir):
    return PfamAnnotationAdapter(genome_dir=genome_dir)


@pytest.fixture
def genome_config_csv(tmp_path, genome_dir):
    """Write a minimal cyanobacteria_genomes.csv with one strain row."""
    csv_path = tmp_path / "genomes.csv"
    csv_path.write_text(
        "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
        f"GCF_000011465.1,Pro_MED4,59919,MED4,{genome_dir},HLI\n"
        "# this comment should be skipped\n",
        encoding="utf-8",
    )
    return str(csv_path)


@pytest.fixture
def pfam_cache_dir(tmp_path):
    """Create a cache dir with a pre-built pfam_reference.json (no downloads)."""
    pfam_dir = tmp_path / "cache" / "pfam"
    pfam_dir.mkdir(parents=True)
    from multiomics_kg.utils.pfam_utils import _pfam_data_to_json
    cache_file = pfam_dir / "pfam_reference.json"
    cache_file.write_text(
        json.dumps(_pfam_data_to_json(MINI_PFAM_DATA)),
        encoding="utf-8",
    )
    return tmp_path / "cache"


@pytest.fixture
def multi_adapter(genome_config_csv, pfam_cache_dir):
    return MultiPfamAnnotationAdapter(
        genome_config_file=genome_config_csv,
        cache_root=pfam_cache_dir,
    )


# ===========================================================================
# Tests for helper ID functions
# ===========================================================================


class TestPfamNodeIdHelpers:
    def test_pfam_node_id_format(self):
        nid = _pfam_node_id("PF00712")
        assert "pfam" in nid.lower()
        assert "PF00712" in nid or "00712" in nid

    def test_pfam_clan_node_id_format(self):
        nid = _pfam_clan_node_id("CL0060")
        assert "pfam.clan" in nid.lower() or "pfam_clan" in nid.lower() or "CL0060" in nid


# ===========================================================================
# Tests for PfamAnnotationAdapter (per-strain)
# ===========================================================================


class TestPfamAdapterGetAllPfamIds:
    def test_returns_correct_ids(self, pfam_adapter):
        ids = pfam_adapter.get_all_pfam_ids()
        assert ids == {"PF00712", "PF02768", "PF00069", "PF99999"}

    def test_empty_and_null_excluded(self, pfam_adapter):
        ids = pfam_adapter.get_all_pfam_ids()
        assert "" not in ids
        assert None not in ids

    def test_missing_key_excluded(self, pfam_adapter):
        """PMM0005 has no pfam_ids key at all."""
        ids = pfam_adapter.get_all_pfam_ids()
        # Should not crash; just skip PMM0005
        assert len(ids) == 4


class TestPfamAdapterEdgeStructure:
    def test_all_edges_are_5_tuples(self, pfam_adapter):
        for edge in pfam_adapter.get_edges():
            assert len(edge) == 5, f"Edge is not a 5-tuple: {edge}"
            edge_id, src, tgt, label, props = edge
            assert isinstance(edge_id, str)
            assert isinstance(src, str)
            assert isinstance(tgt, str)
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_gene_node_id_format(self, pfam_adapter):
        for _, src, _, _, _ in pfam_adapter.get_edges():
            assert src.startswith("ncbigene:"), f"Unexpected source ID: {src}"

    def test_pfam_target_id_format(self, pfam_adapter):
        for _, _, tgt, _, _ in pfam_adapter.get_edges():
            assert "pfam" in tgt.lower(), f"Unexpected target ID: {tgt}"

    def test_edge_label_is_gene_has_pfam(self, pfam_adapter):
        for _, _, _, label, _ in pfam_adapter.get_edges():
            assert label == "gene_has_pfam"


class TestPfamAdapterEdgeCounts:
    def test_pmm0001_two_pfam_edges(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0001" in e[1]]
        assert len(edges) == 2  # PF00712 and PF02768

    def test_pmm0002_one_pfam_edge(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0002" in e[1]]
        assert len(edges) == 1  # PF00069

    def test_pmm0003_empty_list_no_edges(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0003" in e[1]]
        assert len(edges) == 0

    def test_pmm0004_null_no_edges(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0004" in e[1]]
        assert len(edges) == 0

    def test_pmm0005_missing_key_no_edges(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0005" in e[1]]
        assert len(edges) == 0

    def test_pmm0006_one_pfam_edge(self, pfam_adapter):
        edges = [e for e in pfam_adapter.get_edges() if "PMM0006" in e[1]]
        assert len(edges) == 1  # PF99999

    def test_total_edge_count(self, pfam_adapter):
        edges = list(pfam_adapter.get_edges())
        # PMM0001: 2, PMM0002: 1, PMM0006: 1 = 4
        assert len(edges) == 4


class TestPfamAdapterEdgeIds:
    def test_edge_id_contains_pfam(self, pfam_adapter):
        for edge_id, _, _, _, _ in pfam_adapter.get_edges():
            assert "pfam" in edge_id.lower() or "has_pfam" in edge_id

    def test_no_duplicate_edge_ids(self, pfam_adapter):
        edges = list(pfam_adapter.get_edges())
        ids = [e[0] for e in edges]
        assert len(ids) == len(set(ids)), "Duplicate edge IDs in PfamAnnotationAdapter"


class TestPfamAdapterMissingFile:
    def test_missing_json_yields_no_edges(self, tmp_path):
        empty_dir = tmp_path / "empty_strain"
        empty_dir.mkdir()
        adapter = PfamAnnotationAdapter(genome_dir=empty_dir)
        assert list(adapter.get_edges()) == []

    def test_missing_json_get_all_pfam_ids_returns_empty(self, tmp_path):
        empty_dir = tmp_path / "empty_strain"
        empty_dir.mkdir()
        adapter = PfamAnnotationAdapter(genome_dir=empty_dir)
        assert adapter.get_all_pfam_ids() == set()


class TestPfamAdapterTestMode:
    def test_test_mode_limits_edges(self, tmp_path):
        big_genes = {
            f"PMM{i:04d}": {
                "locus_tag": f"PMM{i:04d}",
                "pfam_ids": ["PF00712"],
            }
            for i in range(200)
        }
        gene_file = tmp_path / "gene_annotations_merged.json"
        gene_file.write_text(json.dumps(big_genes), encoding="utf-8")
        adapter = PfamAnnotationAdapter(genome_dir=tmp_path, test_mode=True)
        edges = list(adapter.get_edges())
        assert len(edges) <= 100


# ===========================================================================
# Tests for MultiPfamAnnotationAdapter.get_nodes()
# ===========================================================================


class TestMultiPfamAdapterPfamNodes:
    def test_pfam_node_count(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        pfam_nodes = [n for n in nodes if n[1] == "pfam"]
        # 4 unique PF IDs: PF00712, PF02768, PF00069, PF99999
        assert len(pfam_nodes) == 4

    def test_pfam_node_ids_use_prefix(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        pfam_nodes = [n for n in nodes if n[1] == "pfam"]
        for node_id, _, _ in pfam_nodes:
            assert "pfam" in node_id.lower(), f"Pfam node ID missing prefix: {node_id}"

    def test_pfam_node_has_name_property(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "pfam":
                assert "name" in props, f"Pfam node {node_id} missing 'name'"
                assert isinstance(props["name"], str)
                assert len(props["name"]) > 0

    def test_pfam_node_has_short_name_property(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "pfam":
                assert "short_name" in props, f"Pfam node {node_id} missing 'short_name'"
                assert isinstance(props["short_name"], str)
                assert len(props["short_name"]) > 0

    def test_only_referenced_pfam_ids_emitted(self, multi_adapter):
        """Only PF IDs that genes reference should be emitted as nodes."""
        nodes = list(multi_adapter.get_nodes())
        pfam_nodes = [n for n in nodes if n[1] == "pfam"]
        # Gene data has 4 unique PF IDs; MINI_PFAM_DATA has 4 entries
        # All 4 should be emitted because genes reference all 4
        assert len(pfam_nodes) == 4


class TestMultiPfamAdapterPfamClanNodes:
    def test_pfam_clan_node_count(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        clan_nodes = [n for n in nodes if n[1] == "pfam clan"]
        # CL0060 (DNA_clamp, from PF00712 + PF02768) and CL9999 (Test_clan, from PF99999)
        # PF00069 has no clan
        assert len(clan_nodes) == 2

    def test_pfam_clan_node_ids_use_prefix(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        clan_nodes = [n for n in nodes if n[1] == "pfam clan"]
        for node_id, _, _ in clan_nodes:
            assert "pfam.clan" in node_id.lower() or "CL" in node_id, (
                f"PfamClan node ID format unexpected: {node_id}"
            )

    def test_pfam_clan_node_has_name(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "pfam clan":
                assert "name" in props, f"PfamClan node {node_id} missing 'name'"
                assert isinstance(props["name"], str)
                assert len(props["name"]) > 0

    def test_only_referenced_clans_emitted(self, multi_adapter):
        """Only clans whose domains are referenced by genes should appear."""
        nodes = list(multi_adapter.get_nodes())
        clan_nodes = [n for n in nodes if n[1] == "pfam clan"]
        clan_names = {n[2]["name"] for n in clan_nodes}
        assert "DNA_clamp" in clan_names
        assert "Test_clan" in clan_names


class TestMultiPfamAdapterNoDuplicateNodes:
    def test_no_duplicate_node_ids(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        node_ids = [n[0] for n in nodes]
        assert len(node_ids) == len(set(node_ids)), "Duplicate node IDs in get_nodes()"


class TestMultiPfamAdapterDeduplicationAcrossStrains:
    def test_pfam_nodes_deduplicated(self, tmp_path, pfam_cache_dir):
        """Two strains sharing the same Pfam ID should yield only one Pfam node."""
        shared_genes = {
            "PMM0001": {
                "locus_tag": "PMM0001",
                "pfam_ids": ["PF00712"],
            }
        }
        strain1 = tmp_path / "MED4"
        strain1.mkdir()
        (strain1 / "gene_annotations_merged.json").write_text(
            json.dumps(shared_genes), encoding="utf-8"
        )
        strain2 = tmp_path / "MIT9312"
        strain2.mkdir()
        (strain2 / "gene_annotations_merged.json").write_text(
            json.dumps(shared_genes), encoding="utf-8"
        )

        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            f"GCF_001,59919,MED4,{strain1},HLI\n"
            f"GCF_002,59919,MIT9312,{strain2},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiPfamAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_root=pfam_cache_dir,
        )
        nodes = list(adapter.get_nodes())
        pfam_nodes = [n for n in nodes if n[1] == "pfam"]
        pfam_ids = [n[0] for n in pfam_nodes]
        assert len(pfam_ids) == len(set(pfam_ids)), "Duplicate Pfam nodes across strains"


# ===========================================================================
# Tests for MultiPfamAnnotationAdapter.get_edges()
# ===========================================================================


class TestMultiPfamAdapterEdges:
    def test_gene_has_pfam_edges_present(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        gene_pfam = [e for e in edges if e[3] == "gene_has_pfam"]
        assert len(gene_pfam) == 4  # PMM0001: 2, PMM0002: 1, PMM0006: 1

    def test_pfam_in_pfam_clan_edges_present(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        clan_edges = [e for e in edges if e[3] == "pfam_in_pfam_clan"]
        # PF00712 → CL0060, PF02768 → CL0060, PF99999 → CL9999
        # PF00069 has no clan
        assert len(clan_edges) == 3

    def test_no_clan_edge_for_domain_without_clan(self, multi_adapter):
        """PF00069 (Pkinase) has no clan; should not generate a clan edge."""
        edges = list(multi_adapter.get_edges())
        clan_edges = [e for e in edges if e[3] == "pfam_in_pfam_clan"]
        pfam_sources = {e[1] for e in clan_edges}
        pf00069_node = _pfam_node_id("PF00069")
        assert pf00069_node not in pfam_sources

    def test_clan_edge_source_is_pfam(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if e[3] == "pfam_in_pfam_clan":
                src = e[1]
                assert "pfam" in src.lower(), f"Clan edge source is not a Pfam node: {src}"

    def test_clan_edge_target_is_pfam_clan(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if e[3] == "pfam_in_pfam_clan":
                tgt = e[2]
                assert "pfam.clan" in tgt.lower() or "CL" in tgt, (
                    f"Clan edge target is not a PfamClan node: {tgt}"
                )

    def test_gene_pfam_edge_targets_are_known_nodes(self, multi_adapter):
        """All gene→Pfam edges must point to nodes returned by get_nodes()."""
        node_ids = {n[0] for n in multi_adapter.get_nodes()}
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if e[3] == "gene_has_pfam":
                assert e[2] in node_ids, f"Pfam target {e[2]} not in node set"

    def test_clan_edge_targets_are_known_nodes(self, multi_adapter):
        """All Pfam→PfamClan edges must point to nodes returned by get_nodes()."""
        node_ids = {n[0] for n in multi_adapter.get_nodes()}
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if e[3] == "pfam_in_pfam_clan":
                assert e[2] in node_ids, f"PfamClan target {e[2]} not in node set"


# ===========================================================================
# Tests for string sanitization
# ===========================================================================


class TestPfamStringSanitization:
    def test_single_quotes_sanitized_in_pfam_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for _, label, props in nodes:
            if label == "pfam":
                for val in props.values():
                    assert "'" not in str(val), (
                        f"Raw single-quote in pfam property: {val!r}"
                    )

    def test_pipe_sanitized_in_pfam_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for _, label, props in nodes:
            if label == "pfam":
                for val in props.values():
                    assert "|" not in str(val), (
                        f"Pipe character in pfam property: {val!r}"
                    )

    def test_pf99999_description_cleaned(self, multi_adapter):
        """PF99999 has quotes and pipes in description; must be sanitized."""
        nodes = list(multi_adapter.get_nodes())
        pf99999 = next(
            (n for n in nodes if n[1] == "pfam" and "PF99999" in n[0]),
            None,
        )
        assert pf99999 is not None, "PF99999 node not found"
        desc = pf99999[2]["name"]
        assert "'" not in desc
        assert "|" not in desc

    def test_single_quotes_sanitized_in_clan_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for _, label, props in nodes:
            if label == "pfam clan":
                for val in props.values():
                    assert "'" not in str(val), (
                        f"Raw single-quote in pfam clan property: {val!r}"
                    )


# ===========================================================================
# Tests for CSV parsing edge cases
# ===========================================================================


class TestMultiPfamCsvParsing:
    def test_comment_lines_skipped(self, tmp_path, genome_dir, pfam_cache_dir):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "# this is a comment\n"
            f"GCF_000011465.1,59919,MED4,{genome_dir},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiPfamAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_root=pfam_cache_dir,
        )
        assert len(adapter._strain_adapters) == 1

    def test_empty_data_dir_skipped(self, tmp_path, pfam_cache_dir):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "GCF_000011465.1,59919,MED4,,HLI\n",
            encoding="utf-8",
        )
        adapter = MultiPfamAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_root=pfam_cache_dir,
        )
        assert len(adapter._strain_adapters) == 0
