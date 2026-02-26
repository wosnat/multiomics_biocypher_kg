"""Unit tests for multiomics_kg/download/download_uniprot.py.

Coverage
--------
- UniprotNodeField enum values and class methods (get_split_fields,
  get_nonuniprot_api_fields)
- DEFAULT_NODE_FIELDS composition
- download_uniprot: cache-hit skip, force-overwrite, output file creation
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.download_uniprot import (
    DEFAULT_NODE_FIELDS,
    UniprotNodeField,
    download_uniprot,
)


# ─── UniprotNodeField ─────────────────────────────────────────────────────────

class TestUniprotNodeField:
    def test_primary_gene_name_value(self):
        assert UniprotNodeField.PRIMARY_GENE_NAME.value == "gene_primary"

    def test_organism_id_value(self):
        assert UniprotNodeField.ORGANISM_ID.value == "organism_id"

    def test_go_fields_have_correct_values(self):
        assert UniprotNodeField.CELLULAR_COMPONENT.value == "go_c"
        assert UniprotNodeField.BIOLOGICAL_PROCESS.value == "go_p"
        assert UniprotNodeField.MOLECULAR_FUNCTION.value == "go_f"

    def test_go_id_fields_have_correct_values(self):
        assert UniprotNodeField.CELLULAR_COMPONENT_ID.value == "go_c_id"
        assert UniprotNodeField.BIOLOGICAL_PROCESS_ID.value == "go_p_id"
        assert UniprotNodeField.MOLECULAR_FUNCTION_ID.value == "go_f_id"

    def test_get_split_fields_returns_list(self):
        fields = UniprotNodeField.get_split_fields()
        assert isinstance(fields, list)
        assert len(fields) > 0

    def test_get_split_fields_contains_expected_fields(self):
        fields = UniprotNodeField.get_split_fields()
        assert UniprotNodeField.PROTEOME.value in fields
        assert UniprotNodeField.KEGG_IDS.value in fields
        assert UniprotNodeField.ENTREZ_GENE_IDS.value in fields
        assert UniprotNodeField.CELLULAR_COMPONENT.value in fields

    def test_get_nonuniprot_api_fields_returns_list(self):
        fields = UniprotNodeField.get_nonuniprot_api_fields()
        assert isinstance(fields, list)
        assert len(fields) > 0

    def test_get_nonuniprot_api_fields_contains_go_id_fields(self):
        fields = UniprotNodeField.get_nonuniprot_api_fields()
        assert UniprotNodeField.CELLULAR_COMPONENT_ID.value in fields
        assert UniprotNodeField.BIOLOGICAL_PROCESS_ID.value in fields
        assert UniprotNodeField.MOLECULAR_FUNCTION_ID.value in fields

    def test_get_nonuniprot_api_fields_contains_embedding_fields(self):
        fields = UniprotNodeField.get_nonuniprot_api_fields()
        assert UniprotNodeField.PROTT5_EMBEDDING.value in fields
        assert UniprotNodeField.ESM2_EMBEDDING.value in fields

    def test_split_fields_and_nonuniprot_api_fields_are_disjoint(self):
        split = set(UniprotNodeField.get_split_fields())
        non_api = set(UniprotNodeField.get_nonuniprot_api_fields())
        assert split.isdisjoint(non_api), (
            f"Fields appear in both lists: {split & non_api}"
        )


# ─── DEFAULT_NODE_FIELDS ──────────────────────────────────────────────────────

class TestDefaultNodeFields:
    def test_is_a_list(self):
        assert isinstance(DEFAULT_NODE_FIELDS, list)

    def test_contains_core_fields(self):
        assert UniprotNodeField.PRIMARY_GENE_NAME.value in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.ORGANISM_ID.value in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.PROTEIN_NAMES.value in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.LENGTH.value in DEFAULT_NODE_FIELDS

    def test_contains_go_term_fields(self):
        assert UniprotNodeField.CELLULAR_COMPONENT.value in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.BIOLOGICAL_PROCESS.value in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.MOLECULAR_FUNCTION.value in DEFAULT_NODE_FIELDS

    def test_excludes_embedding_fields(self):
        # Embeddings are handled separately by the adapter, not downloaded here
        assert UniprotNodeField.PROTT5_EMBEDDING.value not in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.ESM2_EMBEDDING.value not in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.NT_EMBEDDING.value not in DEFAULT_NODE_FIELDS

    def test_excludes_go_id_computed_fields(self):
        # go_*_id are derived fields, not fetchable from the API
        assert UniprotNodeField.CELLULAR_COMPONENT_ID.value not in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.BIOLOGICAL_PROCESS_ID.value not in DEFAULT_NODE_FIELDS
        assert UniprotNodeField.MOLECULAR_FUNCTION_ID.value not in DEFAULT_NODE_FIELDS

    def test_no_duplicates(self):
        assert len(DEFAULT_NODE_FIELDS) == len(set(DEFAULT_NODE_FIELDS))


# ─── download_uniprot ─────────────────────────────────────────────────────────

_FETCH_RAW = "multiomics_kg.download.download_uniprot.fetch_raw_uniprot"
# curl is imported lazily inside download_uniprot(); patch the source module
_PYPATH_CURL = "pypath.share.curl"


def _make_mock_curl():
    """Return a mock that supports `with curl.cache_off():`."""
    mock_curl = MagicMock()
    ctx = MagicMock()
    ctx.__enter__ = MagicMock(return_value=None)
    ctx.__exit__ = MagicMock(return_value=False)
    mock_curl.cache_off.return_value = ctx
    return mock_curl


class TestDownloadUniprot:
    def _fake_raw_data(self):
        return ({"gene_primary": {"Q7TU21": "asd"}}, {"Q7TU21"})

    def test_skips_when_cache_exists_and_no_force(self, tmp_path):
        raw_path = tmp_path / "uniprot_raw_data.json"
        raw_path.write_text("{}")

        with patch(_FETCH_RAW) as mock_fetch:
            result = download_uniprot(59919, tmp_path, force=False)

        assert result is False
        mock_fetch.assert_not_called()

    def test_downloads_when_cache_missing(self, tmp_path):
        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()) as mock_fetch,
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            result = download_uniprot(59919, tmp_path, force=False)

        assert result is True
        mock_fetch.assert_called_once()

    def test_force_redownloads_when_cache_exists(self, tmp_path):
        raw_path = tmp_path / "uniprot_raw_data.json"
        raw_path.write_text('{"old": true}')

        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()),
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            result = download_uniprot(59919, tmp_path, force=True)

        assert result is True
        saved = json.loads(raw_path.read_text())
        assert "old" not in saved

    def test_output_file_created(self, tmp_path):
        raw_path = tmp_path / "uniprot_raw_data.json"
        assert not raw_path.exists()

        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()),
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            download_uniprot(59919, tmp_path, force=False)

        assert raw_path.exists()

    def test_output_is_valid_json(self, tmp_path):
        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()),
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            download_uniprot(59919, tmp_path, force=False)

        data = json.loads((tmp_path / "uniprot_raw_data.json").read_text())
        assert isinstance(data, dict)

    def test_creates_cache_dir_if_missing(self, tmp_path):
        subdir = tmp_path / "new" / "nested" / "dir"
        assert not subdir.exists()

        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()),
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            download_uniprot(59919, subdir, force=False)

        assert subdir.exists()
        assert (subdir / "uniprot_raw_data.json").exists()

    def test_passes_node_fields_to_fetch(self, tmp_path):
        custom_fields = ["gene_primary", "length"]

        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()) as mock_fetch,
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            download_uniprot(59919, tmp_path, force=False, node_fields=custom_fields)

        called_fields = mock_fetch.call_args[0][1]
        assert called_fields == custom_fields

    def test_default_node_fields_used_when_none_given(self, tmp_path):
        with (
            patch(_FETCH_RAW, return_value=self._fake_raw_data()) as mock_fetch,
            patch(_PYPATH_CURL, _make_mock_curl()),
        ):
            download_uniprot(59919, tmp_path, force=False)

        called_fields = mock_fetch.call_args[0][1]
        assert called_fields == DEFAULT_NODE_FIELDS
