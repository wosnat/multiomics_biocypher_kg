"""Unit tests for multiomics_kg/download/download_uniprot.py.

Coverage
--------
- _extract_go_id
- _split_field
- preprocess_uniprot_data  (integer conversion, string sanitisation,
                            GO ID extraction, KEGG split, ENTREZ take-first,
                            single-element list unwrapping)
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.download_uniprot import (
    UniprotNodeField,
    _extract_go_id,
    _split_field,
    preprocess_uniprot_data,
)


# ─── _extract_go_id ──────────────────────────────────────────────────────────

class TestExtractGoId:
    def test_extracts_go_id_from_bracket_notation(self):
        term = "aspartate-semialdehyde dehydrogenase activity [GO:0004073]"
        assert _extract_go_id(term) == "0004073"

    def test_extracts_go_id_from_real_cellular_component(self):
        term = "DNA polymerase III complex [GO:0009360]"
        assert _extract_go_id(term) == "0009360"

    def test_returns_none_when_no_go_term(self):
        assert _extract_go_id("no GO term here") is None

    def test_returns_none_for_empty_string(self):
        assert _extract_go_id("") is None

    def test_handles_multiple_go_colons_takes_first(self):
        # Unlikely in practice but tests the split logic
        term = "some term [GO:0001234]"
        assert _extract_go_id(term) == "0001234"


# ─── _split_field ─────────────────────────────────────────────────────────────

class TestSplitField:
    def test_proteomes_split_on_comma(self):
        result = _split_field(
            UniprotNodeField.PROTEOME.value,
            "UP000000625,UP000002311",
        )
        assert result == ["UP000000625", "UP000002311"]

    def test_gene_names_split_on_space(self):
        # Real gene_names value: "dnaN repA" → ["dnaN", "repA"]
        result = _split_field(
            UniprotNodeField.PROTEIN_GENE_NAMES.value,
            "dnaN repA",
        )
        assert result == ["dnaN", "repA"]

    def test_kegg_splits_colon_and_takes_id_part(self):
        # Real KEGG value: "pro:PMM0001" → "PMM0001"
        result = _split_field(
            UniprotNodeField.KEGG_IDS.value,
            "pro:PMM0001",
        )
        assert result == "PMM0001"

    def test_multiple_kegg_entries_colon_split_each(self):
        result = _split_field(
            UniprotNodeField.KEGG_IDS.value,
            "pro:PMM0001;syn:slr0001",
        )
        assert result == ["PMM0001", "slr0001"]

    def test_entrez_gene_ids_takes_first_only(self):
        # Multiple GeneIDs → first only (scalar string)
        result = _split_field(
            UniprotNodeField.ENTREZ_GENE_IDS.value,
            "1234567;7654321",
        )
        assert result == "1234567"

    def test_single_element_result_unwrapped_to_scalar(self):
        # proteomes with one entry → scalar, not list
        result = _split_field(
            UniprotNodeField.PROTEOME.value,
            "UP000000625",
        )
        assert result == "UP000000625"

    def test_none_returns_none(self):
        assert _split_field(UniprotNodeField.PROTEOME.value, None) is None

    def test_empty_string_returns_empty_string(self):
        result = _split_field(UniprotNodeField.PROTEOME.value, "")
        assert result == ""

    def test_pipe_replaced_with_comma_before_split(self):
        # Field sanitisation: | → , (for admin-import compatibility)
        result = _split_field(
            UniprotNodeField.PROTEOME.value,
            "UP000000625|UP000002311",
        )
        # After | → , substitution and split on comma
        assert "UP000000625" in result
        assert "UP000002311" in result

    def test_apostrophe_replaced_with_caret(self):
        result = _split_field(
            UniprotNodeField.PROTEIN_GENE_NAMES.value,
            "O'Brien",
        )
        assert "'" not in result
        assert "^" in result


# ─── preprocess_uniprot_data ──────────────────────────────────────────────────

class TestPreprocessUniprotData:
    """Tests for preprocess_uniprot_data.

    The function modifies `data` in-place and returns (data, locations_set).
    We build minimal `data` dicts targeting one code path per test.
    """

    def _run(self, data: dict, fields: list[str] | None = None):
        from tqdm import tqdm as _tqdm  # imported inside test to avoid top-level dep issues
        if fields is None:
            fields = list(data.keys())
        result_data, locations = preprocess_uniprot_data(data, fields)
        return result_data, locations

    # ── integer conversion ────────────────────────────────────────────────────

    def test_integer_field_length_converted(self):
        data = {UniprotNodeField.LENGTH.value: {"Q7TU21": "367"}}
        result, _ = self._run(data)
        assert result[UniprotNodeField.LENGTH.value]["Q7TU21"] == 367
        assert isinstance(result[UniprotNodeField.LENGTH.value]["Q7TU21"], int)

    def test_integer_field_mass_removes_comma(self):
        data = {UniprotNodeField.MASS.value: {"Q7TU21": "38,274"}}
        result, _ = self._run(data)
        assert result[UniprotNodeField.MASS.value]["Q7TU21"] == 38274

    def test_integer_field_organism_id_converted(self):
        data = {UniprotNodeField.ORGANISM_ID.value: {"Q7TU21": "59919"}}
        result, _ = self._run(data)
        assert result[UniprotNodeField.ORGANISM_ID.value]["Q7TU21"] == 59919

    # ── string sanitisation ───────────────────────────────────────────────────

    def test_pipe_replaced_with_comma_in_string_fields(self):
        # cc_function is a string field (not in split_fields)
        data = {UniprotNodeField.cc_function.value: {"Q7TU21": "domain A|domain B"}}
        result, _ = self._run(data)
        assert "|" not in result[UniprotNodeField.cc_function.value]["Q7TU21"]
        assert "," in result[UniprotNodeField.cc_function.value]["Q7TU21"]

    def test_apostrophe_replaced_with_caret_in_string_fields(self):
        data = {UniprotNodeField.cc_function.value: {"Q7TU21": "O'something"}}
        result, _ = self._run(data)
        val = result[UniprotNodeField.cc_function.value]["Q7TU21"]
        assert "'" not in val
        assert "^" in val

    # ── GO ID extraction ──────────────────────────────────────────────────────

    def test_go_c_id_field_created_from_go_c(self):
        # go_c raw value is a semicolon-delimited string (from UniProt API).
        # preprocess_uniprot_data first splits it via _split_field, then
        # extracts GO IDs into go_c_id.
        data = {
            UniprotNodeField.CELLULAR_COMPONENT.value: {
                "Q7TU21": (
                    "DNA polymerase III complex [GO:0009360];"
                    "cytoplasm [GO:0005737]"
                ),
            }
        }
        result, _ = self._run(data)
        go_c_id = result.get(UniprotNodeField.CELLULAR_COMPONENT_ID.value, {})
        assert "Q7TU21" in go_c_id
        ids = go_c_id["Q7TU21"]
        assert "0009360" in ids
        assert "0005737" in ids

    def test_go_p_id_field_created_from_go_p(self):
        # Use two semicolon-separated terms so _split_field returns a list
        # (a single-element string gets unwrapped to scalar, then character-
        # iteration would yield no GO IDs)
        data = {
            UniprotNodeField.BIOLOGICAL_PROCESS.value: {
                "Q7TU21": "DNA replication [GO:0006260];'de novo' IMP biosynthesis [GO:0006189]",
            }
        }
        result, _ = self._run(data)
        go_p_id = result.get(UniprotNodeField.BIOLOGICAL_PROCESS_ID.value, {})
        assert "0006260" in go_p_id["Q7TU21"]
        assert "0006189" in go_p_id["Q7TU21"]

    def test_go_f_id_field_created_from_go_f(self):
        # Two terms → list → IDs extracted correctly
        data = {
            UniprotNodeField.MOLECULAR_FUNCTION.value: {
                "Q7TU21": "DNA binding [GO:0003677];ATP binding [GO:0005524]",
            }
        }
        result, _ = self._run(data)
        go_f_id = result.get(UniprotNodeField.MOLECULAR_FUNCTION_ID.value, {})
        assert "0003677" in go_f_id["Q7TU21"]

    def test_terms_without_go_colon_excluded_from_id_field(self):
        # Second term has no GO id; only the first should appear in go_c_id
        data = {
            UniprotNodeField.CELLULAR_COMPONENT.value: {
                "Q7TU21": "nucleus [GO:0005634];no GO term here",
            }
        }
        result, _ = self._run(data)
        go_c_id = result[UniprotNodeField.CELLULAR_COMPONENT_ID.value]
        assert len(go_c_id["Q7TU21"]) == 1

    # ── split fields ──────────────────────────────────────────────────────────

    def test_kegg_split_on_colon(self):
        data = {
            UniprotNodeField.KEGG_IDS.value: {"Q7TU21": "pro:PMM0001"}
        }
        result, _ = self._run(data)
        assert result[UniprotNodeField.KEGG_IDS.value]["Q7TU21"] == "PMM0001"

    def test_entrez_gene_ids_takes_first(self):
        data = {
            UniprotNodeField.ENTREZ_GENE_IDS.value: {"Q7TU21": "1234567;7654321"}
        }
        result, _ = self._run(data)
        assert result[UniprotNodeField.ENTREZ_GENE_IDS.value]["Q7TU21"] == "1234567"

    def test_proteome_split_on_comma(self):
        data = {
            UniprotNodeField.PROTEOME.value: {
                "Q7TU21": "UP000000625,UP000002311"
            }
        }
        result, _ = self._run(data)
        val = result[UniprotNodeField.PROTEOME.value]["Q7TU21"]
        assert "UP000000625" in val
        assert "UP000002311" in val

    def test_single_element_list_unwrapped_for_split_fields(self):
        data = {
            UniprotNodeField.PROTEOME.value: {"Q7TU21": "UP000000625"}
        }
        result, _ = self._run(data)
        # Single entry → scalar (not a list)
        val = result[UniprotNodeField.PROTEOME.value]["Q7TU21"]
        assert val == "UP000000625"

    # ── field absent in data ─────────────────────────────────────────────────

    def test_field_not_in_data_is_skipped(self):
        # Passing a field name in `fields` list that's not in `data` → no error
        data = {UniprotNodeField.LENGTH.value: {"Q7TU21": "100"}}
        fields = [UniprotNodeField.LENGTH.value, UniprotNodeField.MASS.value]
        result, _ = preprocess_uniprot_data(data, fields)
        assert UniprotNodeField.MASS.value not in result
