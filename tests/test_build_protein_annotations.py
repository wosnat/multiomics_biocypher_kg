"""Unit tests for multiomics_kg/download/build_protein_annotations.py.

Coverage
--------
- Helper functions: _nonempty, _split, _coerce_to_tokens
- Named transforms: _tx_add_go_prefix, _tx_strip_function_prefix,
                    _tx_extract_go_from_brackets
- load_uniprot_columnar
- ProteinAnnotationBuilder.build_merged  (all field types incl. bool edge cases)
- process_taxid: skip/force/output-file behaviour
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.build_protein_annotations import (
    ProteinAnnotationBuilder,
    load_uniprot_columnar,
    process_taxid,
)
from multiomics_kg.download.utils.annotation_helpers import (
    _coerce_to_tokens,
    _nonempty,
    _split,
)
from multiomics_kg.download.utils.annotation_transforms import (
    _tx_add_go_prefix,
    _tx_clean_catalytic_activity,
    _tx_clean_function_description,
    _tx_extract_cofactor_name,
    _tx_extract_go_from_brackets,
    _tx_extract_pathway_name,
    _tx_extract_signal_range,
    _tx_extract_tm_range,
    _tx_strip_function_prefix,
)

_PROCESS_TAXID_PROJECT_ROOT = "multiomics_kg.download.build_protein_annotations.PROJECT_ROOT"


# ─── _nonempty ────────────────────────────────────────────────────────────────

class TestNonempty:
    def test_none_is_empty(self):
        assert _nonempty(None) is False

    def test_empty_string_is_empty(self):
        assert _nonempty("") is False

    def test_whitespace_only_is_empty(self):
        assert _nonempty("   ") is False

    def test_dash_sentinel_is_empty(self):
        assert _nonempty("-") is False

    def test_dash_with_spaces_is_empty(self):
        assert _nonempty("  -  ") is False

    def test_non_empty_string_is_true(self):
        assert _nonempty("dnaN") is True

    def test_integer_zero_is_true(self):
        assert _nonempty(0) is True

    def test_false_bool_is_true(self):
        # False is a real value; must NOT be treated as empty
        assert _nonempty(False) is True

    def test_empty_list_is_empty(self):
        assert _nonempty([]) is False

    def test_list_with_only_empty_items_is_empty(self):
        assert _nonempty(["", "-", None]) is False

    def test_list_with_one_non_empty_item_is_true(self):
        assert _nonempty(["", "Q7V0G8"]) is True


# ─── _split ───────────────────────────────────────────────────────────────────

class TestSplit:
    def test_splits_on_delimiter(self):
        assert _split("GO:0009360;GO:0003677", ";") == ["GO:0009360", "GO:0003677"]

    def test_strips_whitespace_from_tokens(self):
        assert _split(" asd ; nnrD ", ";") == ["asd", "nnrD"]

    def test_skips_empty_tokens(self):
        assert _split("a;;b", ";") == ["a", "b"]

    def test_skips_dash_sentinel_tokens(self):
        assert _split("a;-;b", ";") == ["a", "b"]

    def test_empty_input_returns_empty_list(self):
        assert _split("", ";") == []

    def test_non_string_returns_empty_list(self):
        assert _split(None, ";") == []

    def test_single_token(self):
        assert _split("Q7TU21", ";") == ["Q7TU21"]


# ─── _coerce_to_tokens ────────────────────────────────────────────────────────

class TestCoerceToTokens:
    def test_string_splits_on_delimiter(self):
        assert _coerce_to_tokens("Q7TU21;Q7V0G8", ";") == ["Q7TU21", "Q7V0G8"]

    def test_list_returns_items_as_tokens(self):
        assert _coerce_to_tokens(["Q7TU21", "Q7V0G8"], ";") == ["Q7TU21", "Q7V0G8"]

    def test_list_skips_empty_and_none(self):
        assert _coerce_to_tokens(["Q7TU21", "", None, "-"], ";") == ["Q7TU21"]

    def test_non_string_non_list_returns_empty(self):
        assert _coerce_to_tokens(42, ";") == []

    def test_empty_list_returns_empty(self):
        assert _coerce_to_tokens([], ";") == []


# ─── _tx_add_go_prefix ────────────────────────────────────────────────────────

class TestTxAddGoPrefix:
    def test_prepends_go_to_bare_7_digit(self):
        # Real GO id stored in go_c_id: "0009360" → "GO:0009360"
        assert _tx_add_go_prefix("0009360") == "GO:0009360"

    def test_already_prefixed_left_unchanged(self):
        assert _tx_add_go_prefix("GO:0009360") == "GO:0009360"

    def test_non_go_token_left_unchanged(self):
        assert _tx_add_go_prefix("PF00712") == "PF00712"

    def test_dash_sentinel_returns_empty(self):
        assert _tx_add_go_prefix("-") == ""

    def test_empty_string_returns_empty(self):
        assert _tx_add_go_prefix("") == ""

    def test_strips_whitespace_before_check(self):
        assert _tx_add_go_prefix("  0009360  ") == "GO:0009360"


# ─── _tx_strip_function_prefix ───────────────────────────────────────────────

class TestTxStripFunctionPrefix:
    def test_strips_function_prefix(self):
        assert _tx_strip_function_prefix(
            "FUNCTION: Catalyzes the replication of chromosomal DNA."
        ) == "Catalyzes the replication of chromosomal DNA."

    def test_case_insensitive(self):
        assert _tx_strip_function_prefix("function: catalyzes") == "catalyzes"

    def test_no_prefix_passes_through(self):
        assert _tx_strip_function_prefix("Catalyzes something.") == "Catalyzes something."

    def test_non_string_returns_empty(self):
        assert _tx_strip_function_prefix(None) == ""

    def test_empty_string_returns_empty(self):
        assert _tx_strip_function_prefix("") == ""


# ─── _tx_extract_go_from_brackets ────────────────────────────────────────────

class TestTxExtractGoFromBrackets:
    def test_extracts_go_id_from_bracket_notation(self):
        term = "DNA polymerase III complex [GO:0009360]"
        assert _tx_extract_go_from_brackets(term) == "GO:0009360"

    def test_extracts_go_id_with_longer_term_name(self):
        term = "aspartate-semialdehyde dehydrogenase activity [GO:0004073]"
        assert _tx_extract_go_from_brackets(term) == "GO:0004073"

    def test_no_go_term_returns_empty(self):
        assert _tx_extract_go_from_brackets("no GO term here") == ""

    def test_empty_string_returns_empty(self):
        assert _tx_extract_go_from_brackets("") == ""

    def test_dash_sentinel_returns_empty(self):
        assert _tx_extract_go_from_brackets("-") == ""


# ─── load_uniprot_columnar ────────────────────────────────────────────────────

# Minimal column-oriented JSON matching real uniprot_raw_data.json format.
# Real UniProt accessions from taxid 59919 (Prochlorococcus MED4):
UNIPROT_COL_DATA = {
    "gene_primary": {"Q7TU21": "asd", "Q7V0G8": "nnrD", "Q7V0H7": "trmD"},
    "organism_id": {"Q7TU21": 59919, "Q7V0G8": 59919, "Q7V0H7": 59919},
    "reviewed": {"Q7TU21": "Reviewed", "Q7V0G8": "Unreviewed"},
    "length": {"Q7TU21": 367, "Q7V0H7": 238},
    # Q7V0G8 deliberately has no 'length' entry (field missing for that accession)
}


class TestLoadUniprotColumnar:
    @pytest.fixture
    def json_path(self, tmp_path):
        p = tmp_path / "uniprot_raw_data.json"
        p.write_text(json.dumps(UNIPROT_COL_DATA))
        return str(p)

    def test_all_accessions_present_as_keys(self, json_path):
        result = load_uniprot_columnar(json_path)
        assert set(result.keys()) == {"Q7TU21", "Q7V0G8", "Q7V0H7"}

    def test_row_contains_correct_field_values(self, json_path):
        result = load_uniprot_columnar(json_path)
        assert result["Q7TU21"]["gene_primary"] == "asd"
        assert result["Q7TU21"]["organism_id"] == 59919
        assert result["Q7TU21"]["reviewed"] == "Reviewed"

    def test_field_missing_for_accession_is_absent(self, json_path):
        # Q7V0G8 has no 'length' in the JSON
        result = load_uniprot_columnar(json_path)
        assert "length" not in result["Q7V0G8"]

    def test_accession_only_in_one_field_still_appears(self, json_path):
        # Q7V0H7 is in gene_primary but not in reviewed
        result = load_uniprot_columnar(json_path)
        assert "Q7V0H7" in result
        assert "reviewed" not in result["Q7V0H7"]


# ─── ProteinAnnotationBuilder.build_merged ───────────────────────────────────

MINIMAL_PROTEIN_CONFIG = {
    "fields": {
        "gene_symbol": {
            "type": "passthrough",
            "field": "gene_primary",
        },
        "is_reviewed": {
            "type": "reviewed_status",
            "field": "reviewed",
        },
        "sequence_length": {
            "type": "integer",
            "field": "length",
        },
        "annotation_score": {
            "type": "float",
            "field": "annotation_score",
        },
        "go_cellular_components": {
            "type": "passthrough_list",
            "field": "go_c",
            "delimiter": ";",
            "transform": "extract_go_from_brackets",
        },
        "function_description": {
            "type": "passthrough",
            "field": "cc_function",
            "transform": "strip_function_prefix",
        },
    }
}

# Representative UniProt row (raw data format, row-oriented)
UP_ROW = {
    "gene_primary": "asd",
    "reviewed": "Reviewed",
    "length": 367,
    "annotation_score": "4.0",
    # Raw GO terms: semicolon-separated "term description [GO:XXXXXXX]" strings
    "go_c": "DNA polymerase III complex [GO:0009360];cytoplasm [GO:0005737]",
    "cc_function": "FUNCTION: Catalyzes the reductive amination of aspartate semialdehyde.",
}


class TestProteinAnnotationBuilderBuildMerged:
    def setup_method(self):
        self.builder = ProteinAnnotationBuilder(MINIMAL_PROTEIN_CONFIG)

    # ── passthrough ──────────────────────────────────────────────────────────

    def test_passthrough_copies_value(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        assert result["gene_symbol"] == "asd"

    def test_passthrough_missing_field_omitted(self):
        result = self.builder.build_merged("Q7TU21", {})
        assert "gene_symbol" not in result

    def test_passthrough_with_strip_function_transform(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        assert result["function_description"].startswith("Catalyzes")
        assert "FUNCTION:" not in result["function_description"]

    # ── reviewed_status ────────────────────────────────────────────────────────

    def test_reviewed_string_gives_reviewed(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        assert result["is_reviewed"] == "reviewed"

    def test_unreviewed_string_gives_not_reviewed(self):
        row = {**UP_ROW, "reviewed": "Unreviewed"}
        result = self.builder.build_merged("Q7V0G8", row)
        assert "is_reviewed" in result
        assert result["is_reviewed"] == "not reviewed"

    def test_none_reviewed_omits_field(self):
        row = {**UP_ROW, "reviewed": None}
        result = self.builder.build_merged("Q7TU21", row)
        assert "is_reviewed" not in result

    def test_unknown_reviewed_value_omits_field(self):
        row = {**UP_ROW, "reviewed": "unknown"}
        result = self.builder.build_merged("Q7TU21", row)
        assert "is_reviewed" not in result

    # ── integer ──────────────────────────────────────────────────────────────

    def test_integer_converts_from_int(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        assert result["sequence_length"] == 367
        assert isinstance(result["sequence_length"], int)

    def test_integer_converts_from_float_string(self):
        row = {**UP_ROW, "length": "1234.0"}
        result = self.builder.build_merged("Q7TU21", row)
        assert result["sequence_length"] == 1234

    def test_integer_non_parseable_omits_field(self):
        row = {**UP_ROW, "length": "notanumber"}
        result = self.builder.build_merged("Q7TU21", row)
        assert "sequence_length" not in result

    # ── float ────────────────────────────────────────────────────────────────

    def test_float_converts_from_string(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        assert abs(result["annotation_score"] - 4.0) < 1e-9
        assert isinstance(result["annotation_score"], float)

    def test_float_non_parseable_omits_field(self):
        row = {**UP_ROW, "annotation_score": "bad"}
        result = self.builder.build_merged("Q7TU21", row)
        assert "annotation_score" not in result

    # ── passthrough_list ──────────────────────────────────────────────────────

    def test_passthrough_list_applies_extract_go_from_brackets_transform(self):
        result = self.builder.build_merged("Q7TU21", UP_ROW)
        go_terms = result["go_cellular_components"]
        assert "GO:0009360" in go_terms
        assert "GO:0005737" in go_terms

    def test_passthrough_list_from_semicolon_delimited_string(self):
        config = {
            "fields": {
                "gene_names": {
                    "type": "passthrough_list",
                    "field": "gene_names",
                    "delimiter": ";",
                }
            }
        }
        builder = ProteinAnnotationBuilder(config)
        row = {"gene_names": "asd;nnrD;trmD"}
        result = builder.build_merged("Q7TU21", row)
        assert result["gene_names"] == ["asd", "nnrD", "trmD"]

    def test_passthrough_list_pre_existing_list_passes_through(self):
        config = {
            "fields": {
                "proteome_ids": {
                    "type": "passthrough_list",
                    "field": "xref_proteomes",
                    "delimiter": ";",
                }
            }
        }
        builder = ProteinAnnotationBuilder(config)
        # Real example: proteomes are a comma-separated string in preprocess data
        row = {"xref_proteomes": ["UP000000625", "UP000002311"]}
        result = builder.build_merged("Q7TU21", row)
        assert result["proteome_ids"] == ["UP000000625", "UP000002311"]

    def test_passthrough_list_missing_field_omitted(self):
        result = self.builder.build_merged("Q7TU21", {})
        assert "go_cellular_components" not in result

    def test_passthrough_list_empty_after_filtering_omitted(self):
        # Terms with no [GO:...] notation → extract_go_from_brackets returns "" → dropped
        row = {**UP_ROW, "go_c": "term one;term two"}
        result = self.builder.build_merged("Q7TU21", row)
        assert "go_cellular_components" not in result

    def test_passthrough_list_dash_sentinel_tokens_excluded(self):
        # Second term has no GO id → filtered out; only GO:0009360 remains
        row = {**UP_ROW, "go_c": "DNA polymerase III complex [GO:0009360];no GO term here"}
        result = self.builder.build_merged("Q7TU21", row)
        assert all(t.startswith("GO:") for t in result["go_cellular_components"])


# ─── process_taxid ────────────────────────────────────────────────────────────

class TestProcessTaxid:
    def _write_raw_data(self, path: Path, data: dict | None = None):
        if data is None:
            data = {
                "gene_primary": {"Q7TU21": "asd"},
                "reviewed": {"Q7TU21": "Reviewed"},
                "length": {"Q7TU21": 367},
            }
        path.write_text(json.dumps(data))

    def test_skips_when_input_missing(self, tmp_path, capsys):
        # No uniprot_raw_data.json → should print skip and not create output
        with patch(_PROCESS_TAXID_PROJECT_ROOT, tmp_path):
            process_taxid("Prochlorococcus", 59919, MINIMAL_PROTEIN_CONFIG, force=False)
        out = capsys.readouterr().out
        assert "skipping" in out.lower()
        output_path = (
            tmp_path / "cache" / "data" / "Prochlorococcus"
            / "uniprot" / "59919" / "protein_annotations.json"
        )
        assert not output_path.exists()

    def test_skips_when_output_exists_and_no_force(self, tmp_path, capsys):
        uniprot_dir = (
            tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        )
        uniprot_dir.mkdir(parents=True)
        self._write_raw_data(uniprot_dir / "uniprot_raw_data.json")
        output_path = uniprot_dir / "protein_annotations.json"
        output_path.write_text('{"old": true}')
        mtime_before = output_path.stat().st_mtime

        with patch(_PROCESS_TAXID_PROJECT_ROOT, tmp_path):
            process_taxid("Prochlorococcus", 59919, MINIMAL_PROTEIN_CONFIG, force=False)

        assert output_path.stat().st_mtime == mtime_before
        assert "Skipping" in capsys.readouterr().out

    def test_force_overwrites_existing_output(self, tmp_path):
        uniprot_dir = (
            tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        )
        uniprot_dir.mkdir(parents=True)
        self._write_raw_data(uniprot_dir / "uniprot_raw_data.json")
        output_path = uniprot_dir / "protein_annotations.json"
        output_path.write_text('{"old": true}')

        with patch(_PROCESS_TAXID_PROJECT_ROOT, tmp_path):
            process_taxid("Prochlorococcus", 59919, MINIMAL_PROTEIN_CONFIG, force=True)

        saved = json.loads(output_path.read_text())
        assert "old" not in saved
        assert "Q7TU21" in saved

    def test_normal_run_creates_output_file(self, tmp_path):
        uniprot_dir = (
            tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        )
        uniprot_dir.mkdir(parents=True)
        self._write_raw_data(uniprot_dir / "uniprot_raw_data.json")
        output_path = uniprot_dir / "protein_annotations.json"

        with patch(_PROCESS_TAXID_PROJECT_ROOT, tmp_path):
            process_taxid("Prochlorococcus", 59919, MINIMAL_PROTEIN_CONFIG, force=False)

        assert output_path.exists()

    def test_output_keyed_by_uniprot_accession(self, tmp_path):
        uniprot_dir = (
            tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        )
        uniprot_dir.mkdir(parents=True)
        self._write_raw_data(uniprot_dir / "uniprot_raw_data.json")

        with patch(_PROCESS_TAXID_PROJECT_ROOT, tmp_path):
            process_taxid("Prochlorococcus", 59919, MINIMAL_PROTEIN_CONFIG, force=False)

        output = json.loads(
            (uniprot_dir / "protein_annotations.json").read_text()
        )
        assert "Q7TU21" in output
        assert output["Q7TU21"]["gene_symbol"] == "asd"


# ─── New transform unit tests ─────────────────────────────────────────────────

class TestTxCleanFunctionDescription:
    def test_strips_function_prefix(self):
        s = "FUNCTION: Catalyzes the reduction. {ECO:0000256,PIRNR:PIRNR000001}."
        assert _tx_clean_function_description(s) == "Catalyzes the reduction"

    def test_strips_eco_tags_only_no_prefix(self):
        s = "Catalyzes the reduction. {ECO:0000256,PIRNR:PIRNR000001}."
        assert _tx_clean_function_description(s) == "Catalyzes the reduction"

    def test_no_eco_passthrough(self):
        assert _tx_clean_function_description("Catalyzes the reduction.") \
               == "Catalyzes the reduction"

    def test_empty_returns_empty(self):
        assert _tx_clean_function_description("") == ""


class TestTxCleanCatalyticActivity:
    def test_strips_catalytic_activity_prefix(self):
        s = "CATALYTIC ACTIVITY: Reaction=ATP + H2O = ADP + phosphate; EC=3.6.1.3; Evidence={ECO:0000256,ARBA:ARB00001}."
        result = _tx_clean_catalytic_activity(s)
        assert "CATALYTIC ACTIVITY:" not in result
        assert "Reaction=ATP" in result

    def test_strips_eco_evidence(self):
        s = "CATALYTIC ACTIVITY: Reaction=A = B; Evidence={ECO:0000256,ARBA:ARBA00001}."
        assert "{ECO:" not in _tx_clean_catalytic_activity(s)

    def test_empty_chunk_returns_empty(self):
        assert _tx_clean_catalytic_activity("") == ""


class TestTxExtractCofactorName:
    def test_extracts_simple_name(self):
        s = "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; Evidence={ECO:0000256|HAMAP-Rule:MF_02225}"
        assert _tx_extract_cofactor_name(s) == "FMN"

    def test_extracts_metal_ion_name(self):
        s = "COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence={ECO:0000256}"
        assert _tx_extract_cofactor_name(s) == "Mg(2+)"

    def test_no_match_returns_empty(self):
        assert _tx_extract_cofactor_name("unrelated text") == ""

    def test_empty_returns_empty(self):
        assert _tx_extract_cofactor_name("") == ""


class TestTxExtractPathwayName:
    def test_strips_prefix_and_eco(self):
        s = "PATHWAY: Energy metabolism; oxidative phosphorylation. {ECO:0000256,ARBA:ARBA00004673}."
        result = _tx_extract_pathway_name(s)
        assert result == "Energy metabolism; oxidative phosphorylation"

    def test_no_eco_still_strips_prefix_and_period(self):
        s = "PATHWAY: Carbohydrate biosynthesis."
        assert _tx_extract_pathway_name(s) == "Carbohydrate biosynthesis"

    def test_empty_returns_empty(self):
        assert _tx_extract_pathway_name("") == ""


class TestTxExtractTmRange:
    def test_extracts_single_range(self):
        s = 'TRANSMEM 32..50; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"'
        assert _tx_extract_tm_range(s) == "32..50"

    def test_extracts_range_ignoring_trailing_annotation(self):
        s = 'TRANSMEM 7..24; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"; TRANSMEM 60..82'
        # Only the first TRANSMEM match is extracted (single-token use)
        assert _tx_extract_tm_range(s) == "7..24"

    def test_no_transmem_returns_empty(self):
        assert _tx_extract_tm_range("no TM here") == ""


class TestTxExtractSignalRange:
    def test_extracts_range(self):
        s = 'SIGNAL 1..26; /evidence="ECO:0000256|HAMAP-Rule:MF_01961"'
        assert _tx_extract_signal_range(s) == "1..26"

    def test_no_signal_returns_empty(self):
        assert _tx_extract_signal_range("TRANSMEM 1..20") == ""

    def test_empty_returns_empty(self):
        assert _tx_extract_signal_range("") == ""


# ─── split_pattern support in _resolve_passthrough_list ───────────────────────

class TestPassthroughListSplitPattern:
    """Tests for the split_pattern parameter in passthrough_list fields."""

    def _make_builder(self, field_cfg: dict) -> ProteinAnnotationBuilder:
        return ProteinAnnotationBuilder({"fields": {"out": field_cfg}})

    def test_single_entry_no_split_needed(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
            "transform": "clean_catalytic_activity",
        })
        row = {"data": "CATALYTIC ACTIVITY: Reaction=A = B; EC=1.1.1.1; Evidence={ECO:0000256}."}
        result = builder.build_merged("X", row)
        assert isinstance(result["out"], list)
        assert len(result["out"]) == 1
        assert "Reaction=A = B" in result["out"][0]

    def test_multi_entry_split_into_list(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
        })
        row = {"data": "CATALYTIC ACTIVITY: Reaction=A = B; CATALYTIC ACTIVITY: Reaction=C = D;"}
        result = builder.build_merged("X", row)
        assert len(result["out"]) == 2

    def test_empty_chunks_after_split_are_dropped(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=COFACTOR:)",
            "transform": "extract_cofactor_name",
        })
        row = {"data": "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; COFACTOR: Name=Mg(2+);"}
        result = builder.build_merged("X", row)
        assert result["out"] == ["FMN", "Mg(2+)"]

    def test_split_pattern_takes_priority_over_delimiter(self):
        """When split_pattern is set, delimiter is ignored."""
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=TRANSMEM)",
            "transform": "extract_tm_range",
            "delimiter": ";",   # should be ignored
        })
        row = {"data": 'TRANSMEM 12..31; /note="Helical"; TRANSMEM 37..57; /note="Helical"'}
        result = builder.build_merged("X", row)
        assert result["out"] == ["12..31", "37..57"]


# ─── Integration tests for new field configs ──────────────────────────────────

RAW_MULTI_CATALYTIC = (
    "CATALYTIC ACTIVITY: Reaction=ATP + H2O = ADP + phosphate; EC=3.6.1.3;"
    " Evidence={ECO:0000256,ARBA:ARBA00048988};"
    " CATALYTIC ACTIVITY: Reaction=C + D = E; EC=1.2.3.4;"
    " Evidence={ECO:0000256,HAMAP-Rule:MF_01920}"
)
RAW_MULTI_COFACTOR = (
    "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; Evidence={ECO:0000256|HAMAP};"
    " Note=Binds 1 FMN.; COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420;"
    " Evidence={ECO:0000256|HAMAP}"
)
RAW_MULTI_PATHWAY = (
    "PATHWAY: Cofactor biosynthesis; CoA: step 2/5. {ECO:0000256|HAMAP}.;"
    " PATHWAY: Cofactor biosynthesis; CoA: step 3/5. {ECO:0000256|HAMAP}."
)
RAW_MULTI_TM = (
    'TRANSMEM 12..31; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius";'
    ' TRANSMEM 37..57; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"'
)
RAW_SIGNAL = 'SIGNAL 1..26; /evidence="ECO:0000256|HAMAP-Rule:MF_01961"'
RAW_KEYWORDS = "ATP-binding;DNA repair;Helicase;Hydrolase"

NEW_FIELDS_CONFIG = {
    "fields": {
        "catalytic_activities": {
            "type": "passthrough_list",
            "field": "cc_catalytic_activity",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
            "transform": "clean_catalytic_activity",
        },
        "cofactor_names": {
            "type": "passthrough_list",
            "field": "cc_cofactor",
            "split_pattern": r"(?=COFACTOR:)",
            "transform": "extract_cofactor_name",
        },
        "pathways": {
            "type": "passthrough_list",
            "field": "cc_pathway",
            "split_pattern": r"(?=PATHWAY:)",
            "transform": "extract_pathway_name",
        },
        "transmembrane_regions": {
            "type": "passthrough_list",
            "field": "ft_transmem",
            "split_pattern": r"(?=TRANSMEM)",
            "transform": "extract_tm_range",
        },
        "signal_peptide": {
            "type": "passthrough",
            "field": "ft_signal",
            "transform": "extract_signal_range",
        },
        "keywords": {
            "type": "passthrough_list",
            "field": "keyword",
            "delimiter": ";",
        },
    }
}


class TestNewFieldConfigs:
    def setup_method(self):
        self.builder = ProteinAnnotationBuilder(NEW_FIELDS_CONFIG)

    def test_catalytic_activities_is_list(self):
        result = self.builder.build_merged("X", {"cc_catalytic_activity": RAW_MULTI_CATALYTIC})
        assert isinstance(result["catalytic_activities"], list)
        assert len(result["catalytic_activities"]) == 2

    def test_catalytic_activities_clean_no_prefix_or_eco(self):
        result = self.builder.build_merged("X", {"cc_catalytic_activity": RAW_MULTI_CATALYTIC})
        for item in result["catalytic_activities"]:
            assert "CATALYTIC ACTIVITY:" not in item
            assert "{ECO:" not in item

    def test_cofactor_names_extracts_names(self):
        result = self.builder.build_merged("X", {"cc_cofactor": RAW_MULTI_COFACTOR})
        assert result["cofactor_names"] == ["FMN", "Mg(2+)"]

    def test_pathways_is_clean_list(self):
        result = self.builder.build_merged("X", {"cc_pathway": RAW_MULTI_PATHWAY})
        pathways = result["pathways"]
        assert isinstance(pathways, list)
        assert len(pathways) == 2
        assert all("{ECO:" not in p for p in pathways)
        assert all("PATHWAY:" not in p for p in pathways)

    def test_transmembrane_regions_is_list_of_ranges(self):
        result = self.builder.build_merged("X", {"ft_transmem": RAW_MULTI_TM})
        assert result["transmembrane_regions"] == ["12..31", "37..57"]

    def test_signal_peptide_is_range_string(self):
        result = self.builder.build_merged("X", {"ft_signal": RAW_SIGNAL})
        assert result["signal_peptide"] == "1..26"

    def test_keywords_is_list(self):
        result = self.builder.build_merged("X", {"keyword": RAW_KEYWORDS})
        assert result["keywords"] == ["ATP-binding", "DNA repair", "Helicase", "Hydrolase"]

    def test_missing_fields_omitted_not_none(self):
        result = self.builder.build_merged("X", {})
        assert "catalytic_activities" not in result
        assert "cofactor_names" not in result
        assert "transmembrane_regions" not in result
