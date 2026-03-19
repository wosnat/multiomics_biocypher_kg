"""
Tests for multiomics_kg/utils/paperconfig_utils.py.

Covers:
- parse_timepoint_hours(): all normalization patterns
- iter_csv_tables() / iter_analyses(): traversal of multi-table configs
- get_organism_for_entry(): organism extraction with fallbacks
- get_publication() / get_paper_name() / get_supplementary_materials()
- load_all_paperconfigs(): skips comments, blank lines, missing files
"""

import pytest
import yaml
from pathlib import Path

from multiomics_kg.utils.paperconfig_utils import (
    parse_timepoint_hours,
    iter_csv_tables,
    iter_analyses,
    get_organism_for_entry,
    get_publication,
    get_paper_name,
    get_supplementary_materials,
    get_experiments,
    get_experiment_for_analysis,
    get_organism_for_analysis,
    load_all_paperconfigs,
    load_paperconfig,
)


# ─── parse_timepoint_hours ─────────────────────────────────────────────


class TestParseTimepointHours:
    """Tests for parse_timepoint_hours() covering all documented patterns."""

    @pytest.mark.parametrize(
        "input_val, expected",
        [
            # Simple hours
            ("4h", 4.0),
            ("0.5h", 0.5),
            ("12h", 12.0),
            ("100h", 100.0),
            # Negative hours
            ("-12h", -12.0),
            ("-0.5h", -0.5),
            # Space between number and h
            ("4 h", 4.0),
            ("0.5 h", 0.5),
            # Day patterns
            ("day 18", 432.0),
            ("Day 2", 48.0),
            ("day 1", 24.0),
            ("Day 0", 0.0),
            # Hours with parenthetical annotation
            ("50h (P added)", 50.0),
            ("4h (some note)", 4.0),
            # Post-inoculation
            ("0.5h post-inoculation", 0.5),
            ("4h post-inoculation", 4.0),
            # Extended darkness — uses absolute time from parenthetical
            ("1h extended darkness (36h)", 36.0),
            ("2h extended darkness (48h)", 48.0),
            # Rescue — returns None
            ("R (rescue: re-fed)", None),
            ("R", None),
            # Pooled days — returns None
            ("days 60+89", None),
            # None and empty
            (None, None),
            ("", None),
            ("  ", None),
        ],
    )
    def test_parse_timepoint_hours(self, input_val, expected):
        result = parse_timepoint_hours(input_val)
        if expected is None:
            assert result is None, f"Expected None for {input_val!r}, got {result}"
        else:
            assert result == pytest.approx(expected), (
                f"Expected {expected} for {input_val!r}, got {result}"
            )


# ─── Sample config fixtures ────────────────────────────────────────────


@pytest.fixture
def multi_table_config():
    """Config with multiple supplementary table types."""
    return {
        "publication": {
            "papername": "Test Paper 2025",
            "papermainpdf": "data/test/paper.pdf",
            "supplementary_materials": {
                "supp_table_1": {
                    "type": "csv",
                    "filename": "data/test/table1.csv",
                    "statistical_analyses": [
                        {
                            "id": "analysis_1",
                            "type": "RNASEQ",
                            "organism": "Prochlorococcus MED4",
                            "name_col": "Gene",
                            "logfc_col": "log2FC",
                        },
                        {
                            "id": "analysis_2",
                            "type": "PROTEOMICS",
                            "organism": "Prochlorococcus MED4",
                            "name_col": "Protein",
                            "logfc_col": "ratio",
                        },
                    ],
                },
                "supp_table_2": {
                    "type": "csv",
                    "filename": "data/test/table2.csv",
                    "statistical_analyses": [
                        {
                            "id": "analysis_3",
                            "type": "RNASEQ",
                            "organism": "Alteromonas EZ55",
                            "name_col": "Gene",
                            "logfc_col": "log2FC",
                        },
                    ],
                },
                "id_map": {
                    "type": "id_translation",
                    "filename": "data/test/id_map.csv",
                    "organism": "Prochlorococcus MED4",
                },
                "gff_entry": {
                    "type": "annotation_gff",
                    "filename": "data/test/genes.gff",
                    "organism": "Prochlorococcus MED4",
                },
            },
        }
    }


@pytest.fixture
def no_type_config():
    """Config where csv tables lack an explicit type field (defaults to csv)."""
    return {
        "publication": {
            "papername": "NoType Paper",
            "supplementary_materials": {
                "table_a": {
                    "filename": "data/test/a.csv",
                    "statistical_analyses": [
                        {"id": "a1", "type": "RNASEQ", "name_col": "Gene"},
                    ],
                },
            },
        }
    }


@pytest.fixture
def empty_config():
    """Config with no supplementary_materials."""
    return {"publication": {"papername": "Empty Paper"}}


@pytest.fixture
def minimal_config():
    """Bare minimum config — no publication block at all."""
    return {}


# ─── iter_csv_tables ────────────────────────────────────────────────────


class TestIterCsvTables:
    def test_yields_only_csv_tables(self, multi_table_config):
        tables = list(iter_csv_tables(multi_table_config))
        keys = [k for k, _ in tables]
        assert "supp_table_1" in keys
        assert "supp_table_2" in keys
        assert "id_map" not in keys
        assert "gff_entry" not in keys
        assert len(tables) == 2

    def test_no_supplementary_materials(self, empty_config):
        assert list(iter_csv_tables(empty_config)) == []

    def test_no_publication(self, minimal_config):
        assert list(iter_csv_tables(minimal_config)) == []

    def test_default_type_is_csv(self, no_type_config):
        tables = list(iter_csv_tables(no_type_config))
        assert len(tables) == 1
        assert tables[0][0] == "table_a"


# ─── iter_analyses ──────────────────────────────────────────────────────


class TestIterAnalyses:
    def test_yields_all_analyses(self, multi_table_config):
        analyses = list(iter_analyses(multi_table_config))
        ids = [a["id"] for _, _, a in analyses]
        assert ids == ["analysis_1", "analysis_2", "analysis_3"]

    def test_table_key_matches(self, multi_table_config):
        analyses = list(iter_analyses(multi_table_config))
        # First two from supp_table_1, third from supp_table_2
        assert analyses[0][0] == "supp_table_1"
        assert analyses[1][0] == "supp_table_1"
        assert analyses[2][0] == "supp_table_2"

    def test_empty_analyses_list(self):
        config = {
            "publication": {
                "supplementary_materials": {
                    "table_x": {
                        "type": "csv",
                        "statistical_analyses": [],
                    },
                }
            }
        }
        assert list(iter_analyses(config)) == []

    def test_no_analyses_key(self):
        config = {
            "publication": {
                "supplementary_materials": {
                    "table_x": {
                        "type": "csv",
                        "filename": "data/test/x.csv",
                    },
                }
            }
        }
        assert list(iter_analyses(config)) == []

    def test_no_tables(self, empty_config):
        assert list(iter_analyses(empty_config)) == []


# ─── get_organism_for_entry ─────────────────────────────────────────────


class TestGetOrganismForEntry:
    def test_direct_organism_field(self, multi_table_config):
        entry = {"organism": "Prochlorococcus MED4"}
        result = get_organism_for_entry(multi_table_config, entry)
        assert result == "Prochlorococcus MED4"

    def test_organism_from_first_analysis(self, multi_table_config):
        entry = {
            "statistical_analyses": [
                {"organism": "Alteromonas MIT1002"},
                {"organism": "Alteromonas EZ55"},
            ]
        }
        result = get_organism_for_entry(multi_table_config, entry)
        assert result == "Alteromonas MIT1002"

    def test_no_organism_anywhere(self, multi_table_config):
        entry = {"statistical_analyses": [{"id": "x"}]}
        result = get_organism_for_entry(multi_table_config, entry)
        assert result is None

    def test_empty_entry(self, multi_table_config):
        result = get_organism_for_entry(multi_table_config, {})
        assert result is None

    def test_strips_quotes_and_whitespace(self, multi_table_config):
        entry = {"organism": '  "Prochlorococcus MED4"  '}
        result = get_organism_for_entry(multi_table_config, entry)
        assert result == "Prochlorococcus MED4"

    def test_id_translation_entry(self, multi_table_config):
        """id_translation entries have organism directly on the entry."""
        entry = {
            "type": "id_translation",
            "organism": "Prochlorococcus MIT9301",
        }
        result = get_organism_for_entry(multi_table_config, entry)
        assert result == "Prochlorococcus MIT9301"


# ─── get_publication / get_paper_name / get_supplementary_materials ─────


class TestPublicationHelpers:
    def test_get_publication(self, multi_table_config):
        pub = get_publication(multi_table_config)
        assert pub["papername"] == "Test Paper 2025"

    def test_get_publication_missing(self, minimal_config):
        assert get_publication(minimal_config) == {}

    def test_get_paper_name(self, multi_table_config):
        assert get_paper_name(multi_table_config) == "Test Paper 2025"

    def test_get_paper_name_fallback_to_dir(self, minimal_config):
        fallback = Path("/data/papers/AuthorYear2025/paperconfig.yaml")
        assert get_paper_name(minimal_config, fallback_path=fallback) == "AuthorYear2025"

    def test_get_paper_name_no_fallback(self, minimal_config):
        assert get_paper_name(minimal_config) == "unknown"

    def test_get_paper_name_prefers_config_over_fallback(self, multi_table_config):
        fallback = Path("/data/papers/Other/paperconfig.yaml")
        assert get_paper_name(multi_table_config, fallback_path=fallback) == "Test Paper 2025"

    def test_get_supplementary_materials(self, multi_table_config):
        sm = get_supplementary_materials(multi_table_config)
        assert "supp_table_1" in sm
        assert "id_map" in sm

    def test_get_supplementary_materials_missing(self, minimal_config):
        assert get_supplementary_materials(minimal_config) == {}

    def test_get_supplementary_materials_top_level(self):
        """Strain-level resource configs have supplementary_materials at top level (no publication block)."""
        config = {
            "supplementary_materials": {
                "id_trans": {
                    "type": "id_translation",
                    "organism": "Prochlorococcus MIT9313",
                }
            }
        }
        sm = get_supplementary_materials(config)
        assert "id_trans" in sm
        assert sm["id_trans"]["type"] == "id_translation"


# ─── load_paperconfig / load_all_paperconfigs ───────────────────────────


class TestLoadPaperconfigs:
    def test_load_single_paperconfig(self, tmp_path):
        pc = tmp_path / "paperconfig.yaml"
        pc.write_text(
            yaml.dump({"publication": {"papername": "TestLoad"}})
        )
        config = load_paperconfig(pc)
        assert config["publication"]["papername"] == "TestLoad"

    def test_load_empty_yaml(self, tmp_path):
        pc = tmp_path / "empty.yaml"
        pc.write_text("")
        config = load_paperconfig(pc)
        assert config == {}

    def test_load_all_skips_comments_and_blanks(self, tmp_path):
        # Create two valid paperconfigs
        dir1 = tmp_path / "Paper1"
        dir1.mkdir()
        pc1 = dir1 / "paperconfig.yaml"
        pc1.write_text(yaml.dump({"publication": {"papername": "Paper1"}}))

        dir2 = tmp_path / "Paper2"
        dir2.mkdir()
        pc2 = dir2 / "paperconfig.yaml"
        pc2.write_text(yaml.dump({"publication": {"papername": "Paper2"}}))

        # Write list file with comments, blank lines, and a missing path
        # Paths must be relative to PROJECT_ROOT, so we use absolute trick
        # by monkeypatching. Instead, write absolute paths and patch.
        list_file = tmp_path / "paperconfig_files.txt"
        list_file.write_text(
            f"# This is a comment\n"
            f"\n"
            f"  \n"
            f"# Another comment\n"
            f"nonexistent/path/paperconfig.yaml\n"
        )

        # We need to handle the PROJECT_ROOT prefix.
        # load_all_paperconfigs prepends PROJECT_ROOT to each path.
        # Use monkeypatch to temporarily set PROJECT_ROOT.
        import multiomics_kg.utils.paperconfig_utils as pcu
        original_root = pcu.PROJECT_ROOT
        try:
            pcu.PROJECT_ROOT = tmp_path
            # Rewrite list file with paths relative to tmp_path
            list_file.write_text(
                f"# comment line\n"
                f"\n"
                f"Paper1/paperconfig.yaml\n"
                f"  \n"
                f"# another comment\n"
                f"Paper2/paperconfig.yaml\n"
                f"missing/paperconfig.yaml\n"
            )
            results = load_all_paperconfigs(list_file)
            assert len(results) == 2
            names = [get_paper_name(cfg) for _, cfg in results]
            assert "Paper1" in names
            assert "Paper2" in names
        finally:
            pcu.PROJECT_ROOT = original_root

    def test_load_all_warns_on_missing(self, tmp_path, capsys):
        """Missing files produce a warning but do not raise."""
        import multiomics_kg.utils.paperconfig_utils as pcu
        original_root = pcu.PROJECT_ROOT
        try:
            pcu.PROJECT_ROOT = tmp_path
            list_file = tmp_path / "list.txt"
            list_file.write_text("nonexistent/paperconfig.yaml\n")
            results = load_all_paperconfigs(list_file)
            assert results == []
            captured = capsys.readouterr()
            assert "warn" in captured.out.lower()
        finally:
            pcu.PROJECT_ROOT = original_root


# ─── Experiment lookup (new format) ──────────────────────────────────


@pytest.fixture
def new_format_config():
    """A paperconfig in new format with experiments block."""
    return {
        "publication": {
            "papername": "Test 2025",
            "experiments": {
                "starvation_med4_rnaseq": {
                    "name": "MED4 N starvation vs replete (RNASEQ)",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": "nitrogen_stress",
                    "treatment_condition": "N starvation",
                    "control_condition": "N replete",
                },
                "coculture_hot1a3_med4_rnaseq": {
                    "name": "MED4 coculture vs axenic (RNASEQ)",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": "coculture",
                    "treatment_condition": "Coculture with HOT1A3",
                    "control_condition": "Axenic",
                    "treatment_organism": "Alteromonas macleodii HOT1A3",
                },
            },
            "supplementary_materials": {
                "table1": {
                    "type": "csv",
                    "filename": "data.csv",
                    "statistical_analyses": [
                        {
                            "id": "DE_starv_4h",
                            "experiment": "starvation_med4_rnaseq",
                            "timepoint": "4h",
                            "timepoint_hours": 4.0,
                            "name_col": "gene",
                            "logfc_col": "logFC",
                        },
                        {
                            "id": "DE_starv_24h",
                            "experiment": "starvation_med4_rnaseq",
                            "timepoint": "24h",
                            "timepoint_hours": 24.0,
                            "name_col": "gene",
                            "logfc_col": "logFC",
                        },
                    ],
                },
                "table2": {
                    "type": "csv",
                    "filename": "coculture.csv",
                    "statistical_analyses": [
                        {
                            "id": "DE_cocult",
                            "experiment": "coculture_hot1a3_med4_rnaseq",
                            "name_col": "gene",
                            "logfc_col": "logFC",
                        },
                    ],
                },
            },
        }
    }


class TestExperimentLookup:
    def test_get_experiments(self, new_format_config):
        exps = get_experiments(new_format_config)
        assert len(exps) == 2
        assert "starvation_med4_rnaseq" in exps
        assert "coculture_hot1a3_med4_rnaseq" in exps

    def test_get_experiments_missing(self, minimal_config):
        assert get_experiments(minimal_config) == {}

    def test_get_experiment_for_analysis(self, new_format_config):
        analysis = {"id": "test", "experiment": "starvation_med4_rnaseq"}
        exp = get_experiment_for_analysis(new_format_config, analysis)
        assert exp["organism"] == "Prochlorococcus MED4"
        assert exp["omics_type"] == "RNASEQ"

    def test_get_experiment_for_analysis_missing_ref(self, new_format_config):
        analysis = {"id": "test"}
        with pytest.raises(ValueError, match="missing 'experiment'"):
            get_experiment_for_analysis(new_format_config, analysis)

    def test_get_experiment_for_analysis_bad_ref(self, new_format_config):
        analysis = {"id": "test", "experiment": "nonexistent"}
        with pytest.raises(ValueError, match="unknown experiment"):
            get_experiment_for_analysis(new_format_config, analysis)

    def test_get_organism_for_analysis(self, new_format_config):
        analysis = {"id": "test", "experiment": "starvation_med4_rnaseq"}
        assert get_organism_for_analysis(new_format_config, analysis) == "Prochlorococcus MED4"

    def test_get_organism_for_entry_new_format(self, new_format_config):
        """get_organism_for_entry works with new-format analyses (via experiment block)."""
        table = new_format_config["publication"]["supplementary_materials"]["table1"]
        org = get_organism_for_entry(new_format_config, table)
        assert org == "Prochlorococcus MED4"

    def test_get_organism_for_entry_id_translation_unchanged(self):
        """id_translation entries still use direct organism field."""
        config = {"publication": {}}
        entry = {"type": "id_translation", "organism": "Prochlorococcus MIT9313"}
        assert get_organism_for_entry(config, entry) == "Prochlorococcus MIT9313"
