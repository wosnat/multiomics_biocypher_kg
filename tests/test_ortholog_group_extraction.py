"""Tests for ortholog group extraction utilities."""

import pytest

from multiomics_kg.download.utils.ortholog_group_utils import (
    ORGANISM_GROUP_LEVELS,
    _parse_eggnog_ogs,
    extract_ortholog_groups,
    organism_group_from_path,
)


# ── organism_group_from_path ─────────────────────────────────────────────────


class TestOrganismGroupFromPath:
    def test_prochlorococcus(self):
        assert organism_group_from_path("cache/data/Prochlorococcus/genomes/MED4/") == "Prochlorococcus"

    def test_synechococcus(self):
        assert organism_group_from_path("cache/data/Synechococcus/genomes/CC9311/") == "Synechococcus"

    def test_alteromonas(self):
        assert organism_group_from_path("cache/data/Alteromonas/genomes/MIT1002/") == "Alteromonas"

    def test_unknown(self):
        assert organism_group_from_path("/some/other/path/") == "unknown"


# ── _parse_eggnog_ogs ────────────────────────────────────────────────────────


class TestParseEggnogOgs:
    def test_standard_entries(self):
        ogs = ["COG0592@2|Bacteria", "1MKTR@1212|Prochloraceae"]
        parsed = _parse_eggnog_ogs(ogs)
        assert parsed[2] == ("COG0592", "Bacteria")
        assert parsed[1212] == ("1MKTR", "Prochloraceae")

    def test_legacy_entries_skipped(self):
        """Entries without '@' are legacy short names and should be skipped."""
        ogs = ["COG0592", "1MKTR"]
        parsed = _parse_eggnog_ogs(ogs)
        assert parsed == {}

    def test_missing_pipe_skipped(self):
        """Entries with '@' but no '|' are malformed and skipped."""
        parsed = _parse_eggnog_ogs(["COG0592@2"])
        assert parsed == {}

    def test_non_numeric_taxon_skipped(self):
        parsed = _parse_eggnog_ogs(["COG0592@abc|Bacteria"])
        assert parsed == {}

    def test_empty_list(self):
        assert _parse_eggnog_ogs([]) == {}


# ── extract_ortholog_groups ──────────────────────────────────────────────────


class TestExtractOrthologGroups:
    def test_pro_gene_three_groups(self):
        """Prochlorococcus gene with cluster_number + eggnog_ogs → 3 groups."""
        gene = {
            "cluster_number": "CK_00000001",
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "1MKTR@1212|Prochloraceae",
                "3Q4X@1117|Cyanobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 3

        # Cyanorak
        assert groups[0] == {
            "og_id": "cyanorak:CK_00000001",
            "source": "cyanorak",
            "taxonomic_level": "curated",
            "taxon_id": 0,
        }
        # Bacteria COG
        assert groups[1] == {
            "og_id": "eggnog:COG0592@2",
            "source": "eggnog",
            "taxonomic_level": "Bacteria",
            "taxon_id": 2,
        }
        # Prochloraceae lowest-level
        assert groups[2] == {
            "og_id": "eggnog:1MKTR@1212",
            "source": "eggnog",
            "taxonomic_level": "Prochloraceae",
            "taxon_id": 1212,
        }

    def test_alt_gene_two_groups(self):
        """Alteromonas gene without cluster → 2 groups (bacteria + lowest)."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "4648R@72275|Alteromonadaceae",
            ],
        }
        groups = extract_ortholog_groups(gene, "Alteromonas")
        assert len(groups) == 2
        assert groups[0]["og_id"] == "eggnog:COG0592@2"
        assert groups[1]["og_id"] == "eggnog:4648R@72275"
        assert groups[1]["taxonomic_level"] == "Alteromonadaceae"

    def test_no_eggnog(self):
        """Gene with only cluster_number → 1 group."""
        gene = {"cluster_number": "CK_00000001"}
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 1
        assert groups[0]["source"] == "cyanorak"

    def test_empty_gene(self):
        """Gene with no annotations → empty list."""
        groups = extract_ortholog_groups({}, "Prochlorococcus")
        assert groups == []

    def test_legacy_format_skipped(self):
        """Entries without '@' in eggnog_ogs are ignored."""
        gene = {"eggnog_ogs": ["COG0592", "1MKTR"]}
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert groups == []

    def test_lowest_level_uses_whitelist_not_max_taxon_id(self):
        """Gene with cross-lineage OG (Pleurocapsales@52604) → picks Prochloraceae@1212."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "AAAA@52604|Pleurocapsales",     # cross-lineage, should NOT be picked
                "1MKTR@1212|Prochloraceae",       # target level, should be picked
                "3Q4X@1117|Cyanobacteria",        # fallback, not needed
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        # Should have bacteria + Prochloraceae (NOT Pleurocapsales)
        og_ids = [g["og_id"] for g in groups]
        assert "eggnog:1MKTR@1212" in og_ids
        assert "eggnog:AAAA@52604" not in og_ids

    def test_lowest_level_falls_back_to_cyanobacteria(self):
        """Gene missing target level but has Cyanobacteria@1117 → uses fallback."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "3Q4X@1117|Cyanobacteria",       # fallback level
                # No Prochloraceae@1212 entry
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:3Q4X@1117"
        assert groups[1]["taxonomic_level"] == "Cyanobacteria"
        assert groups[1]["taxon_id"] == 1117

    def test_dedup_og_ids(self):
        """Gene with duplicate OG entries → no duplicate in output."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "COG0592@2|Bacteria",   # duplicate
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        og_ids = [g["og_id"] for g in groups]
        assert og_ids.count("eggnog:COG0592@2") == 1

    def test_synechococcus_target_level(self):
        """Synechococcus gene uses Synechococcus@1129 as target level."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "XXXX@1129|Synechococcus",
            ],
        }
        groups = extract_ortholog_groups(gene, "Synechococcus")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:XXXX@1129"
        assert groups[1]["taxon_id"] == 1129

    def test_unknown_organism_group(self):
        """Unknown organism group → only bacteria-level, no lowest-level."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "1MKTR@1212|Prochloraceae",
            ],
        }
        groups = extract_ortholog_groups(gene, "unknown")
        assert len(groups) == 1
        assert groups[0]["og_id"] == "eggnog:COG0592@2"

    def test_bacteria_only_gene(self):
        """Gene with only bacteria-level OG → 1 eggnog group."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 1
        assert groups[0]["og_id"] == "eggnog:COG0592@2"

    def test_none_eggnog_ogs(self):
        """Gene with eggnog_ogs=None → handled gracefully."""
        gene = {"eggnog_ogs": None}
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert groups == []


# ── ORGANISM_GROUP_LEVELS coverage ───────────────────────────────────────────


class TestOrganismGroupLevels:
    def test_all_three_groups_defined(self):
        assert "Prochlorococcus" in ORGANISM_GROUP_LEVELS
        assert "Synechococcus" in ORGANISM_GROUP_LEVELS
        assert "Alteromonas" in ORGANISM_GROUP_LEVELS

    def test_taxon_ids_are_positive_integers(self):
        for group, (target, fallback) in ORGANISM_GROUP_LEVELS.items():
            assert isinstance(target, int) and target > 0, f"{group} target"
            assert isinstance(fallback, int) and fallback > 0, f"{group} fallback"
