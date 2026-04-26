"""Unit tests for multiomics_kg/utils/gene_id_utils.py.

Coverage
--------
- ORGANISM_TO_GENOME_DIR / SKIP_PATTERNS constants
- get_genome_dir: known organisms, case/whitespace normalization, missing dirs
- is_organism_loaded: known, partial-match, unknown, empty
- load_gene_annotations: merged.json, fallback to wide.json, both missing, malformed JSON
- build_id_lookup: identity mapping, single-value fields, list fields, string fallback,
  supp CSV merge, missing annotations
- build_field_lookups: per-field dicts, list fields, missing annotations
- load_supp_lookup: well-formed CSV, missing file, malformed CSV, duplicate alt_ids
- map_gene_id: direct, lookup, supp, repadded (zero-padding), composite (comma-split),
  no match, empty/whitespace input
- is_id_like_column: keyword match, excluded names, numeric columns,
  low-cardinality, long-text, float-like, space-heavy, no keyword short values
- load_import_report: explicit path, missing file, Docker fallback skipped,
  return_by_source mode
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from multiomics_kg.utils.gene_id_utils import (
    ANNOTATION_ID_FIELDS,
    DESCRIPTION_COL_KEYWORDS,
    ID_COL_KEYWORDS,
    ORGANISM_TO_GENOME_DIR,
    SKIP_PATTERNS,
    build_field_lookups,
    build_id_lookup,
    get_genome_dir,
    is_id_like_column,
    is_organism_loaded,
    load_gene_annotations,
    load_import_report,
    load_supp_lookup,
    map_gene_id,
)


# ─── Fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_genome_dir(tmp_path):
    """Return a temp directory acting as a genome cache dir."""
    return tmp_path


@pytest.fixture
def minimal_annotations():
    """Minimal gene_annotations_merged.json content (two entries)."""
    return {
        "PMM0001": {
            "locus_tag_ncbi": "A9601_RS00100",
            "locus_tag_cyanorak": "CK_Pro_MED4_00001",
            "protein_id": "WP_011128533.1",
            "gene_name": "rpsB",
            "gene_synonyms": ["rps2", "PMM_0001"],
            "old_locus_tags": ["PMM1344"],
        },
        "PMM0002": {
            "locus_tag_ncbi": "A9601_RS00200",
            "locus_tag_cyanorak": "CK_Pro_MED4_00002",
            "protein_id": "WP_011128534.1",
            "gene_name": "rpsC",
            "gene_synonyms": [],
            "old_locus_tags": [],
        },
    }


@pytest.fixture
def annotations_file(tmp_genome_dir, minimal_annotations):
    """Write gene_annotations_merged.json and return (genome_dir, annotations)."""
    path = tmp_genome_dir / "gene_annotations_merged.json"
    path.write_text(json.dumps(minimal_annotations))
    return tmp_genome_dir, minimal_annotations


@pytest.fixture
def supp_csv(tmp_genome_dir):
    """Write a gene_mapping_supp.csv with two alt-ID rows."""
    df = pd.DataFrame({
        "locus_tag": ["PMM0001", "PMM0002"],
        "alt_id": ["old_gene_1", "old_gene_2"],
    })
    path = tmp_genome_dir / "gene_mapping_supp.csv"
    df.to_csv(path, index=False)
    return tmp_genome_dir


# ─── Constants ────────────────────────────────────────────────────────────────

class TestConstants:
    def test_organism_map_has_prochlorococcus_strains(self):
        keys = list(ORGANISM_TO_GENOME_DIR.keys())
        assert any("prochlorococcus med4" in k for k in keys)
        assert any("mit9313" in k for k in keys)

    def test_organism_map_has_alteromonas_aliases(self):
        # Both "alteromonas mit1002" and "alteromonas macleodii mit1002" should map
        # to the same directory.
        path1 = ORGANISM_TO_GENOME_DIR["alteromonas mit1002"]
        path2 = ORGANISM_TO_GENOME_DIR["alteromonas macleodii mit1002"]
        assert path1 == path2

    def test_skip_patterns_match_rna_features(self):
        for feature in ["tRNA_1", "ncRNA_2", "rRNA_3", "Yfr12", "tmRNA"]:
            assert SKIP_PATTERNS.match(feature), f"Expected {feature!r} to match SKIP_PATTERNS"

    def test_skip_patterns_do_not_match_protein_coding(self):
        for gene in ["PMM0001", "rpsB", "MIT9313_0001", "A9601_RS00100"]:
            assert not SKIP_PATTERNS.match(gene), f"Expected {gene!r} NOT to match SKIP_PATTERNS"

    def test_annotation_id_fields_contains_required_keys(self):
        required = {"locus_tag", "gene_name", "protein_id", "gene_synonyms", "old_locus_tags"}
        assert required.issubset(set(ANNOTATION_ID_FIELDS))


# ─── get_genome_dir ───────────────────────────────────────────────────────────

class TestGetGenomeDir:
    def test_exact_known_organism(self, tmp_path):
        # Create the expected sub-directory so it exists.
        rel = ORGANISM_TO_GENOME_DIR["prochlorococcus med4"]
        full = tmp_path / rel
        full.mkdir(parents=True)
        result = get_genome_dir("Prochlorococcus MED4", str(tmp_path))
        assert result == str(full)

    def test_case_insensitive(self, tmp_path):
        rel = ORGANISM_TO_GENOME_DIR["prochlorococcus med4"]
        full = tmp_path / rel
        full.mkdir(parents=True)
        result = get_genome_dir("PROCHLOROCOCCUS MED4", str(tmp_path))
        assert result == str(full)

    def test_strips_quotes_and_whitespace(self, tmp_path):
        rel = ORGANISM_TO_GENOME_DIR["prochlorococcus med4"]
        full = tmp_path / rel
        full.mkdir(parents=True)
        result = get_genome_dir('  "Prochlorococcus MED4"  ', str(tmp_path))
        assert result == str(full)

    def test_returns_none_for_unknown_organism(self, tmp_path):
        result = get_genome_dir("Escherichia coli K12", str(tmp_path))
        assert result is None

    def test_returns_none_when_dir_does_not_exist(self, tmp_path):
        # Known organism but directory not created.
        result = get_genome_dir("Prochlorococcus MED4", str(tmp_path))
        assert result is None

    def test_returns_none_for_empty_string(self, tmp_path):
        assert get_genome_dir("", str(tmp_path)) is None

    def test_returns_none_for_none(self, tmp_path):
        assert get_genome_dir(None, str(tmp_path)) is None

    def test_partial_match_organism(self, tmp_path):
        # "Prochlorococcus MIT9313" is a known key; passing just "MIT9313" should still
        # resolve because of the substring check.
        rel = ORGANISM_TO_GENOME_DIR["prochlorococcus mit9313"]
        full = tmp_path / rel
        full.mkdir(parents=True)
        result = get_genome_dir("prochlorococcus mit9313 strain X", str(tmp_path))
        assert result == str(full)


# ─── is_organism_loaded ───────────────────────────────────────────────────────

class TestIsOrganismLoaded:
    def test_exact_match(self):
        assert is_organism_loaded("prochlorococcus med4") is True

    def test_case_insensitive(self):
        assert is_organism_loaded("Prochlorococcus MED4") is True

    def test_partial_match_superset(self):
        # Longer string that contains a known key
        assert is_organism_loaded("Prochlorococcus MED4 strain") is True

    def test_partial_match_subset(self):
        # Shorter string contained in a known key — e.g. "med4"
        # Known key "prochlorococcus med4" contains "med4"
        assert is_organism_loaded("med4") is True

    def test_unknown_organism(self):
        assert is_organism_loaded("Bacillus subtilis") is False

    def test_empty_string(self):
        assert is_organism_loaded("") is False

    def test_none(self):
        assert is_organism_loaded(None) is False

    def test_alteromonas_alias(self):
        assert is_organism_loaded("Alteromonas macleodii EZ55") is True
        assert is_organism_loaded("Alteromonas EZ55") is True


# ─── load_gene_annotations ────────────────────────────────────────────────────

class TestLoadGeneAnnotations:
    def test_loads_merged_json(self, annotations_file):
        genome_dir, expected = annotations_file
        result = load_gene_annotations(str(genome_dir))
        assert result == expected

    def test_falls_back_to_wide_json(self, tmp_genome_dir, minimal_annotations):
        wide_path = tmp_genome_dir / "gene_annotations_wide.json"
        wide_path.write_text(json.dumps(minimal_annotations))
        result = load_gene_annotations(str(tmp_genome_dir))
        assert result == minimal_annotations

    def test_merged_takes_priority_over_wide(self, tmp_genome_dir, minimal_annotations):
        merged_data = {"PMM9999": {"locus_tag_ncbi": "X"}}
        (tmp_genome_dir / "gene_annotations_merged.json").write_text(json.dumps(merged_data))
        (tmp_genome_dir / "gene_annotations_wide.json").write_text(json.dumps(minimal_annotations))
        result = load_gene_annotations(str(tmp_genome_dir))
        assert result == merged_data

    def test_returns_none_when_no_file(self, tmp_genome_dir):
        result = load_gene_annotations(str(tmp_genome_dir))
        assert result is None

    def test_returns_none_on_malformed_json(self, tmp_genome_dir):
        (tmp_genome_dir / "gene_annotations_merged.json").write_text("{broken json")
        result = load_gene_annotations(str(tmp_genome_dir))
        assert result is None


# ─── build_id_lookup ─────────────────────────────────────────────────────────

class TestBuildIdLookup:
    def test_identity_mapping(self, annotations_file):
        genome_dir, _ = annotations_file
        lookup, locus_tags, _ = build_id_lookup(str(genome_dir))
        assert "PMM0001" in locus_tags
        assert lookup["PMM0001"] == "PMM0001"

    def test_single_value_fields_mapped(self, annotations_file):
        genome_dir, _ = annotations_file
        lookup, _, _ = build_id_lookup(str(genome_dir))
        assert lookup["A9601_RS00100"] == "PMM0001"
        assert lookup["CK_Pro_MED4_00001"] == "PMM0001"
        assert lookup["WP_011128533.1"] == "PMM0001"
        assert lookup["rpsB"] == "PMM0001"

    def test_list_field_gene_synonyms(self, annotations_file):
        genome_dir, _ = annotations_file
        lookup, _, _ = build_id_lookup(str(genome_dir))
        assert lookup["rps2"] == "PMM0001"
        assert lookup["PMM_0001"] == "PMM0001"

    def test_list_field_old_locus_tags(self, annotations_file):
        genome_dir, _ = annotations_file
        lookup, _, _ = build_id_lookup(str(genome_dir))
        assert lookup["PMM1344"] == "PMM0001"

    def test_empty_list_fields_ignored(self, annotations_file):
        # PMM0002 has empty gene_synonyms and old_locus_tags — should not cause errors
        genome_dir, _ = annotations_file
        lookup, locus_tags, _ = build_id_lookup(str(genome_dir))
        assert "PMM0002" in locus_tags

    def test_returns_none_when_no_annotations(self, tmp_genome_dir):
        result = build_id_lookup(str(tmp_genome_dir))
        assert result == (None, None, set())

    def test_supp_csv_merged_as_fallback(self, tmp_genome_dir, minimal_annotations):
        (tmp_genome_dir / "gene_annotations_merged.json").write_text(
            json.dumps(minimal_annotations)
        )
        df = pd.DataFrame({"locus_tag": ["PMM0001"], "alt_id": ["extra_alias"]})
        df.to_csv(tmp_genome_dir / "gene_mapping_supp.csv", index=False)

        lookup, _, supp_keys = build_id_lookup(str(tmp_genome_dir))
        assert lookup["extra_alias"] == "PMM0001"
        assert "extra_alias" in supp_keys

    def test_supp_csv_does_not_overwrite_primary(self, tmp_genome_dir, minimal_annotations):
        # "rpsB" is already mapped via gene_name; supp should not overwrite.
        (tmp_genome_dir / "gene_annotations_merged.json").write_text(
            json.dumps(minimal_annotations)
        )
        df = pd.DataFrame({"locus_tag": ["PMM0002"], "alt_id": ["rpsB"]})
        df.to_csv(tmp_genome_dir / "gene_mapping_supp.csv", index=False)

        lookup, _, supp_keys = build_id_lookup(str(tmp_genome_dir))
        assert lookup["rpsB"] == "PMM0001"   # original wins
        assert "rpsB" not in supp_keys

    def test_string_gene_synonyms_comma_split(self, tmp_genome_dir):
        """Fallback path: gene_synonyms stored as comma-delimited string."""
        data = {
            "PMM0005": {
                "gene_synonyms": "synA, synB",
                "old_locus_tags": [],
            }
        }
        (tmp_genome_dir / "gene_annotations_merged.json").write_text(json.dumps(data))
        lookup, _, _ = build_id_lookup(str(tmp_genome_dir))
        assert lookup.get("synA") == "PMM0005"
        assert lookup.get("synB") == "PMM0005"


# ─── build_field_lookups ─────────────────────────────────────────────────────

class TestBuildFieldLookups:
    def test_returns_none_when_no_annotations(self, tmp_genome_dir):
        result = build_field_lookups(str(tmp_genome_dir))
        assert result is None

    def test_single_value_fields(self, annotations_file):
        genome_dir, _ = annotations_file
        lookups = build_field_lookups(str(genome_dir))
        assert lookups is not None
        assert lookups["locus_tag_ncbi"]["A9601_RS00100"] == "PMM0001"
        assert lookups["gene_name"]["rpsB"] == "PMM0001"
        assert lookups["protein_id"]["WP_011128533.1"] == "PMM0001"

    def test_list_field_gene_synonyms(self, annotations_file):
        genome_dir, _ = annotations_file
        lookups = build_field_lookups(str(genome_dir))
        assert lookups["gene_synonyms"]["rps2"] == "PMM0001"
        assert lookups["gene_synonyms"]["PMM_0001"] == "PMM0001"

    def test_list_field_old_locus_tags(self, annotations_file):
        genome_dir, _ = annotations_file
        lookups = build_field_lookups(str(genome_dir))
        assert lookups["old_locus_tags"]["PMM1344"] == "PMM0001"

    def test_empty_list_does_not_add_entry(self, annotations_file):
        genome_dir, _ = annotations_file
        lookups = build_field_lookups(str(genome_dir))
        # PMM0002 has empty gene_synonyms — should not add empty key
        if "gene_synonyms" in lookups:
            assert "" not in lookups["gene_synonyms"]


# ─── load_supp_lookup ─────────────────────────────────────────────────────────

class TestLoadSuppLookup:
    def test_loads_valid_csv(self, supp_csv):
        lookup = load_supp_lookup(str(supp_csv))
        assert lookup["old_gene_1"] == "PMM0001"
        assert lookup["old_gene_2"] == "PMM0002"

    def test_returns_empty_when_file_missing(self, tmp_genome_dir):
        result = load_supp_lookup(str(tmp_genome_dir))
        assert result == {}

    def test_skips_nan_rows(self, tmp_genome_dir):
        df = pd.DataFrame({"locus_tag": ["PMM0001", None], "alt_id": ["ok_id", "bad_id"]})
        df.to_csv(tmp_genome_dir / "gene_mapping_supp.csv", index=False)
        lookup = load_supp_lookup(str(tmp_genome_dir))
        assert "ok_id" in lookup
        assert "bad_id" not in lookup

    def test_first_occurrence_wins_on_duplicates(self, tmp_genome_dir):
        df = pd.DataFrame({
            "locus_tag": ["PMM0001", "PMM0002"],
            "alt_id": ["dup_id", "dup_id"],
        })
        df.to_csv(tmp_genome_dir / "gene_mapping_supp.csv", index=False)
        lookup = load_supp_lookup(str(tmp_genome_dir))
        assert lookup["dup_id"] == "PMM0001"

    def test_strips_whitespace(self, tmp_genome_dir):
        df = pd.DataFrame({"locus_tag": ["  PMM0001  "], "alt_id": ["  alt  "]})
        df.to_csv(tmp_genome_dir / "gene_mapping_supp.csv", index=False)
        lookup = load_supp_lookup(str(tmp_genome_dir))
        assert "alt" in lookup
        assert lookup["alt"] == "PMM0001"


# ─── map_gene_id ─────────────────────────────────────────────────────────────

@pytest.fixture
def lookup_and_tags():
    """A small lookup dict and locus_tag set for map_gene_id tests."""
    locus_tags = {"PMM0001", "PMM0002", "PMM00010"}
    lookup = {
        "PMM0001": "PMM0001",   # identity
        "PMM0002": "PMM0002",
        "PMM00010": "PMM00010",
        "rpsB": "PMM0001",      # gene_name alias
        "WP_011128": "PMM0002", # protein_id alias
        "supp_alias": "PMM0001",
    }
    supp_keys = {"supp_alias"}
    return lookup, locus_tags, supp_keys


class TestMapGeneId:
    def test_direct_match(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("PMM0001", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "direct"

    def test_lookup_match(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("rpsB", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "lookup"

    def test_supp_method(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("supp_alias", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "supp"

    def test_no_match_returns_none(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("UNKNOWN_GENE", lookup, locus_tags, supp_keys)
        assert lt is None
        assert method is None

    def test_whitespace_is_stripped(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("  PMM0001  ", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "direct"

    def test_zero_padding_single_pad(self, lookup_and_tags):
        # The regex requires PREFIX_DIGITS format (underscore separator).
        # "MIT1002_001" (3-digit) → padded to "MIT1002_0001" (4-digit)
        locus_tags_ext = {"MIT1002_0001"}
        lookup_ext = {"MIT1002_0001": "MIT1002_0001"}
        lt, method = map_gene_id("MIT1002_001", lookup_ext, locus_tags_ext, set())
        assert lt == "MIT1002_0001"
        assert method == "repadded"

    def test_zero_padding_via_lookup(self, lookup_and_tags):
        # If the padded version is in lookup (not locus_tags), it should still resolve.
        # "MIT1002_0010" (4-digit) → padded to "MIT1002_00010" (5-digit)
        locus_tags = set()
        lookup = {"MIT1002_00010": "PMM_TARGET"}
        lt, method = map_gene_id("MIT1002_0010", lookup, locus_tags, set())
        assert lt == "PMM_TARGET"
        assert method == "repadded"

    def test_composite_comma_split_direct(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("PMM0001,PMM0002", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "composite_direct"

    def test_composite_comma_split_lookup(self, lookup_and_tags):
        lookup, locus_tags, supp_keys = lookup_and_tags
        lt, method = map_gene_id("unknown_part,rpsB", lookup, locus_tags, supp_keys)
        assert lt == "PMM0001"
        assert method == "composite_lookup"

    def test_supp_keys_default_to_empty_set(self, lookup_and_tags):
        lookup, locus_tags, _ = lookup_and_tags
        # Calling without supp_keys argument should not raise
        lt, method = map_gene_id("PMM0001", lookup, locus_tags)
        assert lt == "PMM0001"


# ─── is_id_like_column ───────────────────────────────────────────────────────

def _series(*values):
    return pd.Series(list(values))


class TestIsIdLikeColumn:
    """Tests for is_id_like_column heuristics."""

    def test_gene_keyword_col_accepted(self):
        s = _series("PMM0001", "PMM0002", "PMM0003")
        assert is_id_like_column(s, "Gene_ID", exclude_cols=set()) is True

    def test_locus_keyword_col_accepted(self):
        s = _series("A9601_RS00100", "A9601_RS00200", "A9601_RS00300")
        assert is_id_like_column(s, "locus_tag", exclude_cols=set()) is True

    def test_description_keyword_excluded(self):
        s = _series("hypothetical protein", "ribosomal subunit", "chaperone")
        assert is_id_like_column(s, "product", exclude_cols=set()) is False

    def test_excluded_col_name_rejected(self):
        s = _series("PMM0001", "PMM0002", "PMM0003")
        assert is_id_like_column(s, "Gene_ID", exclude_cols={"Gene_ID"}) is False

    def test_purely_numeric_rejected(self):
        s = _series("1", "2", "3", "4", "5")
        assert is_id_like_column(s, "count", exclude_cols=set()) is False

    def test_low_cardinality_categorical_rejected(self):
        # 2 unique values, 10 total — ratio ≥ 5 → categorical
        s = _series(*["Up", "Down"] * 5)
        assert is_id_like_column(s, "direction", exclude_cols=set()) is False

    def test_long_text_rejected(self):
        long_desc = "A" * 60
        s = _series(long_desc, long_desc, long_desc)
        assert is_id_like_column(s, "description", exclude_cols=set()) is False

    def test_float_like_fold_change_rejected(self):
        s = _series("-0.54", "1.8", "0.25", "2.1", "-3.0")
        assert is_id_like_column(s, "log2FC", exclude_cols=set()) is False

    def test_float_with_asterisk_rejected(self):
        s = _series("-0.54*", "1.8*", "0.25", "2.1", "-3.0")
        assert is_id_like_column(s, "fold_change", exclude_cols=set()) is False

    def test_empty_series_rejected(self):
        s = pd.Series([], dtype=str)
        assert is_id_like_column(s, "Gene", exclude_cols=set()) is False

    def test_all_nan_rejected(self):
        s = pd.Series([None, None, None])
        assert is_id_like_column(s, "Gene", exclude_cols=set()) is False

    def test_positional_keyword_excluded(self):
        s = _series("100", "200", "300")
        assert is_id_like_column(s, "start_position", exclude_cols=set()) is False

    def test_space_heavy_values_without_keyword_rejected(self):
        # Values with spaces >30% — no ID keyword
        s = _series("gene one", "gene two", "gene three", "gene four")
        assert is_id_like_column(s, "some_col", exclude_cols=set()) is False

    def test_short_values_without_keyword_accepted(self):
        # Values without spaces, no keyword in column name
        s = _series("PMM0001", "PMM0002", "PMM0003", "PMM0004")
        assert is_id_like_column(s, "entry", exclude_cols=set()) is True


# ─── load_import_report ───────────────────────────────────────────────────────

SAMPLE_REPORT = (
    "doi:10.1234/paper1 (global id space)-[Changes_expression_of]->ncbigene:PMM0001 (global id space)\n"
    "doi:10.1234/paper1 (global id space)-[Changes_expression_of]->ncbigene:PMM0002 (global id space)\n"
    "doi:10.1234/paper2 (global id space)-[Changes_expression_of]->ncbigene:PMM0003 (global id space)\n"
    "unrelated line that should be ignored\n"
)


class TestLoadImportReport:
    def test_explicit_path_returns_set(self, tmp_path):
        report = tmp_path / "import.report"
        report.write_text(SAMPLE_REPORT)
        result = load_import_report(import_report_path=str(report))
        assert result == {"PMM0001", "PMM0002", "PMM0003"}

    def test_explicit_path_return_by_source(self, tmp_path):
        report = tmp_path / "import.report"
        report.write_text(SAMPLE_REPORT)
        result = load_import_report(import_report_path=str(report), return_by_source=True)
        assert result is not None
        assert "all_missing" in result
        assert "by_source" in result
        assert result["all_missing"] == {"PMM0001", "PMM0002", "PMM0003"}
        by_src = result["by_source"]
        assert "doi:10.1234/paper1" in by_src
        assert {"PMM0001", "PMM0002"} == by_src["doi:10.1234/paper1"]
        assert {"PMM0003"} == by_src["doi:10.1234/paper2"]

    def test_returns_none_when_file_missing_and_no_docker(self, tmp_path):
        # No local report, no docker — patch subprocess to fail
        with patch("multiomics_kg.utils.gene_id_utils.subprocess.run") as mock_run:
            mock_run.side_effect = Exception("no docker")
            result = load_import_report(
                import_report_path=str(tmp_path / "nonexistent.report")
            )
        assert result is None

    def test_returns_none_on_empty_report(self, tmp_path):
        report = tmp_path / "import.report"
        report.write_text("nothing relevant here\n")
        result = load_import_report(import_report_path=str(report))
        assert result is None

    def test_local_fallback_path(self, tmp_path, monkeypatch):
        # Patch Path("output/import.report").exists() by writing a real file
        # and monkeypatching the working dir.
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        report = output_dir / "import.report"
        report.write_text(SAMPLE_REPORT)
        monkeypatch.chdir(tmp_path)
        # No explicit path → should fall back to output/import.report
        result = load_import_report()
        assert result == {"PMM0001", "PMM0002", "PMM0003"}

    def test_unrelated_lines_ignored(self, tmp_path):
        report = tmp_path / "import.report"
        report.write_text(
            "WARNING: something\n"
            "doi:10.1234/paper1 (global id space)-[Changes_expression_of]->ncbigene:PMM0001 (global id space)\n"
            "ERROR: something else\n"
        )
        result = load_import_report(import_report_path=str(report))
        assert result == {"PMM0001"}


class TestExtractUniprotAnnotationTokens:
    """Unit tests for extract_uniprot_annotation_tokens (Biller 2022 driver case)."""

    def test_full_input_entry_and_gn_locus_tag(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = ("Q31DF2_PROM9 Type II secretion system protein-like protein "
             "OS=Prochlorococcus marinus (strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")
        assert extract_uniprot_annotation_tokens(s) == [
            ("Q31DF2_PROM9", "uniprot_entry_name"),
            ("PMT9312_0032", "gene_name"),
        ]

    def test_full_input_entry_and_gn_gene_symbol(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = ("RL33_PROM9 50S ribosomal protein L33 "
             "OS=Prochlorococcus marinus (strain MIT 9312) GN=rpmG PE=3 SV=2")
        assert extract_uniprot_annotation_tokens(s) == [
            ("RL33_PROM9", "uniprot_entry_name"),
            ("rpmG", "gene_name"),
        ]

    def test_entry_only_no_gn(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = "Q7V5C8_PROMM Possible serine protease OS=Prochlorococcus marinus"
        assert extract_uniprot_annotation_tokens(s) == [
            ("Q7V5C8_PROMM", "uniprot_entry_name"),
        ]

    def test_gn_only_no_entry_prefix(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        # No leading <entry>_<ORG> shape; only GN= token
        s = "putative serine protease GN=PMT_1636 PE=4 SV=1"
        assert extract_uniprot_annotation_tokens(s) == [
            ("PMT_1636", "gene_name"),
        ]

    def test_neither_pattern_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        assert extract_uniprot_annotation_tokens("hypothetical protein") == []

    def test_empty_and_non_string_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        assert extract_uniprot_annotation_tokens("") == []
        assert extract_uniprot_annotation_tokens("   ") == []
        assert extract_uniprot_annotation_tokens(None) == []  # type: ignore[arg-type]
        assert extract_uniprot_annotation_tokens(123) == []   # type: ignore[arg-type]

    def test_first_gn_wins_when_multiple(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_annotation_tokens
        s = "X_Y abcd GN=first PE=1 GN=second"
        result = extract_uniprot_annotation_tokens(s)
        # First entry-name match + first GN= match
        assert ("first", "gene_name") in result
        assert ("second", "gene_name") not in result


class TestUniprotAnnotationStringValidatorAccepts:
    """Validator must accept the new id_type; the free-text guard must NOT fire."""

    def test_id_type_in_valid_set(self):
        from scripts.validate_paperconfig import VALID_ID_TYPES
        assert "uniprot_annotation_string" in VALID_ID_TYPES


class TestExtractNcbiDeflineTokens:
    """Unit tests for extract_ncbi_defline_tokens (Biller 2014/2022 + many older
    papers' supplementary tables embed BLAST-style deflines)."""

    def test_refseq_np_accession(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        s = "gi|33860650|ref|NP_892211.1| serine protease"
        assert extract_ncbi_defline_tokens(s) == [("NP_892211.1", "protein_id_refseq")]

    def test_refseq_wp_accession(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        s = "gi|494836200|ref|WP_011131639.1| DNA polymerase III subunit beta"
        assert extract_ncbi_defline_tokens(s) == [("WP_011131639.1", "protein_id_refseq")]

    def test_insdc_gb_accession(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        s = "gi|123|gb|CAE18549.1| possible serine protease"
        assert extract_ncbi_defline_tokens(s) == [("CAE18549.1", "protein_id_refseq")]

    def test_emb_and_dbj_sources(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        assert extract_ncbi_defline_tokens("gi|9|emb|CAA12345.1| desc") == [
            ("CAA12345.1", "protein_id_refseq"),
        ]
        assert extract_ncbi_defline_tokens("gi|9|dbj|BAA12345.1| desc") == [
            ("BAA12345.1", "protein_id_refseq"),
        ]

    def test_multiple_deflines_in_one_cell(self):
        """Some cells have multiple BLAST hits separated by commas/semicolons."""
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        s = ("gi|1|ref|NP_111.1| first hit; "
             "gi|2|ref|NP_222.1| second hit")
        assert extract_ncbi_defline_tokens(s) == [
            ("NP_111.1", "protein_id_refseq"),
            ("NP_222.1", "protein_id_refseq"),
        ]

    def test_no_defline_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        assert extract_ncbi_defline_tokens("plain product description") == []
        assert extract_ncbi_defline_tokens("NP_892211.1 alone") == []  # no gi|...| prefix

    def test_empty_and_non_string_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_ncbi_defline_tokens
        assert extract_ncbi_defline_tokens("") == []
        assert extract_ncbi_defline_tokens(None) == []  # type: ignore[arg-type]
        assert extract_ncbi_defline_tokens(42) == []    # type: ignore[arg-type]


class TestNcbiProteinDeflineValidatorAccepts:
    def test_id_type_in_valid_set(self):
        from scripts.validate_paperconfig import VALID_ID_TYPES
        assert "ncbi_protein_defline" in VALID_ID_TYPES


class TestExtractUniprotDeflineTokens:
    """Unit tests for extract_uniprot_defline_tokens (Kratzl 2024 style:
    UniProt FASTA deflines pasted into a CSV column)."""

    def test_swissprot_defline(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_defline_tokens
        s = "sp|Q31L36|RF1_SYNE7"
        assert extract_uniprot_defline_tokens(s) == [
            ("Q31L36", "uniprot_accession"),
            ("RF1_SYNE7", "uniprot_entry_name"),
        ]

    def test_trembl_defline_with_description(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_defline_tokens
        s = "tr|E0IXR1|E0IXR1_ECOLW Sucrose permease OS=Escherichia coli"
        assert extract_uniprot_defline_tokens(s) == [
            ("E0IXR1", "uniprot_accession"),
            ("E0IXR1_ECOLW", "uniprot_entry_name"),
        ]

    def test_multiple_deflines_in_one_cell(self):
        """Real Kratzl 2024 cell: free-text + ;sp|...|... + free-text + ;sp|...|..."""
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_defline_tokens
        s = ("C-phycocyanin beta chain OS=Synechococcus elongatus GN=cpcB1 "
             "PE=1 SV=2;sp|P00312|PHCB1_SYNP6 C-phycocyanin beta chain "
             "OS=Synechococcus PE=1 SV=2;sp|P00308|PHCA1_SYNP6 alpha chain")
        assert extract_uniprot_defline_tokens(s) == [
            ("P00312", "uniprot_accession"),
            ("PHCB1_SYNP6", "uniprot_entry_name"),
            ("P00308", "uniprot_accession"),
            ("PHCA1_SYNP6", "uniprot_entry_name"),
        ]

    def test_no_defline_returns_empty(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_defline_tokens
        assert extract_uniprot_defline_tokens("amt1") == []
        assert extract_uniprot_defline_tokens("plain product description") == []
        # Plain accession without sp|/tr| prefix is NOT extracted
        assert extract_uniprot_defline_tokens("Q31L36") == []

    def test_empty_and_non_string(self):
        from multiomics_kg.utils.gene_id_utils import extract_uniprot_defline_tokens
        assert extract_uniprot_defline_tokens("") == []
        assert extract_uniprot_defline_tokens(None) == []  # type: ignore[arg-type]
        assert extract_uniprot_defline_tokens(42) == []    # type: ignore[arg-type]


class TestUniprotDeflineValidatorAccepts:
    def test_id_type_in_valid_set(self):
        from scripts.validate_paperconfig import VALID_ID_TYPES
        assert "uniprot_defline" in VALID_ID_TYPES


class TestResolveRowUniprotAnnotation:
    """resolve_row must extract uniprot_annotation_string columns at lookup time."""

    @pytest.fixture
    def mapping_data(self):
        from multiomics_kg.utils.gene_id_utils import MappingData
        # PMT9312_0032 is a real locus tag; specific_lookup also has the
        # uniprot_entry_name token Q31DF2_PROM9.
        return MappingData(
            locus_tags={"PMT9312_0032"},
            specific_lookup={
                "Q31DF2_PROM9": "PMT9312_0032",
                "PMT9312_0032": "PMT9312_0032",
            },
            multi_lookup={"rpmG": ["PMT9312_0107"]},
            conflicts={},
        )

    def test_resolves_via_entry_name(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        row = {"Annot": ("Q31DF2_PROM9 protein OS=... "
                        "GN=NOT_A_KNOWN_TAG PE=4 SV=1")}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"
        assert "Annot" in method  # method string mentions the column

    def test_resolves_via_gn_locus_tag_when_entry_missing(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        # No leading entry name → only GN= token; that token is itself a real
        # locus_tag (Tier 1 specific_lookup hit).
        row = {"Annot": "putative protein GN=PMT9312_0032 PE=3 SV=1"}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, _ = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"

    def test_unresolvable_when_no_tokens_match(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        row = {"Annot": "ZZZZZ_NOTREAL protein OS=Foo GN=NOT_A_GENE PE=4 SV=1"}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt is None
        assert method == "unresolved"

    def test_skips_cell_with_no_extractable_tokens(self, mapping_data):
        from multiomics_kg.utils.gene_id_utils import resolve_row
        # When the cell has no entry/GN= tokens, that column contributes
        # nothing — but other columns can still resolve.
        row = {
            "Annot": "hypothetical protein",
            "Other": "PMT9312_0032",
        }
        id_columns = [
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
            {"column": "Other", "id_type": "locus_tag"},
        ]
        lt, _ = resolve_row(row, "Annot", id_columns, mapping_data)
        assert lt == "PMT9312_0032"

    def test_does_not_match_long_annotation_as_literal(self, mapping_data):
        """Regression guard: the raw long annotation string must NOT be looked up
        verbatim — that would either miss (best case) or pollute the fallback
        path (if the literal happens to coincide with a multi_lookup key)."""
        from multiomics_kg.utils.gene_id_utils import resolve_row
        long_annot = ("Q31DF2_PROM9 protein OS=Prochlorococcus marinus "
                      "(strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")
        row = {"Annot": long_annot}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, method = resolve_row(row, "Annot", id_columns, mapping_data)
        # Resolves via the extracted entry_name, not the literal long string
        assert lt == "PMT9312_0032"

    def test_excludes_noise_tokens_from_multi_lookup(self):
        """Regression guard: noise tokens like `OS`, `PE`, `protein` from the
        long annotation string must NOT feed Tier 2/3 multi_lookup fallback.
        Without _candidate_values narrowing, expand_list's word-tokenizer
        would surface `OS` and produce a spurious singleton resolution."""
        from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
        md = MappingData(
            locus_tags={"PMT9312_0032"},
            specific_lookup={},  # entry name not registered, so Tier 1 cannot help
            multi_lookup={"OS": ["WRONG_LT"]},  # noise token would catch here under old logic
            conflicts={},
        )
        row = {"Annot": ("Q31DF2_PROM9 protein OS=Prochlorococcus marinus "
                        "(strain MIT 9312) GN=somekey PE=4 SV=1")}
        id_columns = [{"column": "Annot", "id_type": "uniprot_annotation_string"}]
        lt, _ = resolve_row(row, "Annot", id_columns, md)
        # Must NOT resolve to WRONG_LT via the noise token. Acceptable outcomes:
        #   - lt is None (no other token resolves)
        #   - lt is something else if some legitimate token happens to match
        # The critical assertion is the absence of WRONG_LT.
        assert lt != "WRONG_LT"
