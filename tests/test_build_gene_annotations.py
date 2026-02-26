"""Unit tests for multiomics_kg/download/build_gene_annotations.py.

Coverage
--------
- Helper functions: _nonempty, _split, _coerce_to_tokens
- Named transforms: all six functions
- infer_organism_group
- load_gene_mapping, load_eggnog, load_uniprot
- AnnotationBuilder: build_wide, build_merged (all field types)
- process_strain: skip / force / output files / stats
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.build_gene_annotations import (
    AnnotationBuilder,
    _coerce_to_tokens,
    _nonempty,
    _split,
    _tx_add_go_prefix,
    _tx_extract_go_from_pipe,
    _tx_extract_pfam_ids,
    _tx_first_token_space,
    _tx_strip_function_prefix,
    _tx_strip_prefix_ko,
    infer_organism_group,
    load_eggnog,
    load_gene_mapping,
    load_uniprot,
    process_strain,
)


# ─── _nonempty ────────────────────────────────────────────────────────────────

class TestNonempty:
    def test_none_is_empty(self):
        assert _nonempty(None) is False

    def test_empty_string_is_empty(self):
        assert _nonempty("") is False

    def test_whitespace_string_is_empty(self):
        assert _nonempty("   ") is False

    def test_dash_sentinel_is_empty(self):
        assert _nonempty("-") is False

    def test_dash_with_spaces_is_empty(self):
        assert _nonempty("  -  ") is False

    def test_non_empty_string(self):
        assert _nonempty("hello") is True

    def test_zero_int_is_truthy(self):
        assert _nonempty(0) is True

    def test_false_bool_is_truthy(self):
        # Only None / empty-str / "-" are falsy; False is a real value
        assert _nonempty(False) is True

    def test_empty_list_is_empty(self):
        assert _nonempty([]) is False

    def test_list_with_empty_items_is_empty(self):
        assert _nonempty(["", "-", None]) is False

    def test_list_with_non_empty_item(self):
        assert _nonempty(["", "value"]) is True

    def test_non_empty_list(self):
        assert _nonempty(["a", "b"]) is True


# ─── _split ───────────────────────────────────────────────────────────────────

class TestSplit:
    def test_simple_comma(self):
        assert _split("a,b,c", ",") == ["a", "b", "c"]

    def test_smart_comma_does_not_split_on_comma_space(self):
        # "A, description with comma" should stay as one token
        result = _split("GO:0001,DNA replication, initiation", ",")
        assert result == ["GO:0001", "DNA replication, initiation"]

    def test_plain_delimiter(self):
        assert _split("a|b|c", "|") == ["a", "b", "c"]

    def test_strips_whitespace(self):
        # Comma NOT followed by space → splits; then tokens are stripped
        assert _split(" a,b ", ",") == ["a", "b"]

    def test_comma_space_is_not_split(self):
        # Smart comma: ", " (comma-space) is treated as internal separator, not a split point
        assert _split("a , b", ",") == ["a , b"]

    def test_skips_empty_tokens(self):
        assert _split("a,,b", ",") == ["a", "b"]

    def test_skips_dash_tokens(self):
        assert _split("a,-,b", ",") == ["a", "b"]

    def test_empty_input(self):
        assert _split("", ",") == []

    def test_non_string_input(self):
        assert _split(None, ",") == []

    def test_single_token(self):
        assert _split("hello", ",") == ["hello"]


# ─── _coerce_to_tokens ────────────────────────────────────────────────────────

class TestCoerceToTokens:
    def test_string_input(self):
        assert _coerce_to_tokens("a,b,c", ",") == ["a", "b", "c"]

    def test_list_input(self):
        assert _coerce_to_tokens(["a", "b", "c"], ",") == ["a", "b", "c"]

    def test_list_skips_empty_items(self):
        assert _coerce_to_tokens(["a", "", None, "b"], ",") == ["a", "b"]

    def test_non_string_non_list(self):
        assert _coerce_to_tokens(42, ",") == []

    def test_empty_list(self):
        assert _coerce_to_tokens([], ",") == []


# ─── Named transforms ─────────────────────────────────────────────────────────

class TestTxFirstTokenSpace:
    def test_single_word(self):
        assert _tx_first_token_space("dnaN") == "dnaN"

    def test_multiple_words(self):
        assert _tx_first_token_space("dnaN rps3") == "dnaN"

    def test_leading_whitespace(self):
        assert _tx_first_token_space("  psbA  gene") == "psbA"

    def test_empty_string(self):
        assert _tx_first_token_space("") == ""

    def test_non_string(self):
        assert _tx_first_token_space(None) == ""


class TestTxAddGoPrefix:
    def test_bare_7_digit(self):
        assert _tx_add_go_prefix("0009360") == "GO:0009360"

    def test_already_prefixed(self):
        assert _tx_add_go_prefix("GO:0009360") == "GO:0009360"

    def test_non_go_term(self):
        assert _tx_add_go_prefix("IPR001234") == "IPR001234"

    def test_dash_sentinel(self):
        assert _tx_add_go_prefix("-") == ""

    def test_empty_string(self):
        assert _tx_add_go_prefix("") == ""

    def test_strips_whitespace(self):
        assert _tx_add_go_prefix("  0009360  ") == "GO:0009360"


class TestTxStripFunctionPrefix:
    def test_strips_prefix(self):
        assert _tx_strip_function_prefix("FUNCTION: Catalyzes something") == "Catalyzes something"

    def test_case_insensitive(self):
        assert _tx_strip_function_prefix("function: Catalyzes") == "Catalyzes"

    def test_no_prefix(self):
        assert _tx_strip_function_prefix("Plain description") == "Plain description"

    def test_non_string(self):
        assert _tx_strip_function_prefix(None) == ""

    def test_empty(self):
        assert _tx_strip_function_prefix("") == ""


class TestTxStripPrefixKo:
    def test_strips_ko_prefix(self):
        assert _tx_strip_prefix_ko("ko:K02710") == "K02710"

    def test_case_insensitive(self):
        assert _tx_strip_prefix_ko("KO:K02710") == "K02710"

    def test_no_prefix(self):
        assert _tx_strip_prefix_ko("K02710") == "K02710"

    def test_empty(self):
        assert _tx_strip_prefix_ko("") == ""


class TestTxExtractGoFromPipe:
    def test_pipe_format(self):
        assert _tx_extract_go_from_pipe("DNA replication|0006260||IEA") == "GO:0006260"

    def test_bare_7_digit_fallback(self):
        assert _tx_extract_go_from_pipe("0006260") == "GO:0006260"

    def test_already_prefixed_fallback(self):
        assert _tx_extract_go_from_pipe("GO:0006260") == "GO:0006260"

    def test_pipe_with_non_digit_second_part(self):
        # Non-7-digit second segment → falls back to _tx_add_go_prefix on the *full* original string,
        # which doesn't match GO: or 7-digits → returned unchanged
        result = _tx_extract_go_from_pipe("some_term|NOTANID||IEA")
        assert result == "some_term|NOTANID||IEA"

    def test_empty(self):
        assert _tx_extract_go_from_pipe("") == ""

    def test_dash_sentinel(self):
        assert _tx_extract_go_from_pipe("-") == ""


class TestTxExtractPfamIds:
    def test_keeps_pf_tokens(self):
        assert _tx_extract_pfam_ids("TIGR00663,PF00712,IPR022634") == ["PF00712"]

    def test_multiple_pfams(self):
        assert _tx_extract_pfam_ids("PF00001,PF00002") == ["PF00001", "PF00002"]

    def test_no_pfams(self):
        assert _tx_extract_pfam_ids("TIGR00663,IPR022634") == []

    def test_non_string(self):
        assert _tx_extract_pfam_ids(None) == []

    def test_empty(self):
        assert _tx_extract_pfam_ids("") == []


# ─── infer_organism_group ─────────────────────────────────────────────────────

class TestInferOrganismGroup:
    @pytest.mark.parametrize("data_dir,expected", [
        ("cache/data/Prochlorococcus/genomes/MED4/",  "Prochlorococcus"),
        ("cache/data/Prochlorococcus/genomes/MED4",   "Prochlorococcus"),
        ("cache/data/Synechococcus/genomes/CC9311/",  "Synechococcus"),
        ("cache/data/Alteromonas/genomes/MIT1002/",   "Alteromonas"),
    ])
    def test_standard_paths(self, data_dir, expected):
        assert infer_organism_group(data_dir) == expected

    def test_absolute_path(self):
        assert infer_organism_group("/home/user/cache/data/Alteromonas/genomes/EZ55") == "Alteromonas"

    def test_no_data_segment_returns_default(self):
        assert infer_organism_group("/some/random/path") == "Prochlorococcus"


# ─── load_gene_mapping ────────────────────────────────────────────────────────

GENE_MAPPING_CSV = """\
locus_tag,protein_id,gene_names,product,product_cyanorak,start,end,strand
PMM0001,WP_011129038.1,dnaN,DNA polymerase III subunit beta,DNA polymerase III beta subunit,1,1200,+
PMM0002,WP_011129039.1,rpsT,30S ribosomal protein S20,,1250,1550,+
PMM0003,,,,hypothetical protein,1600,1900,-
"""


class TestLoadGeneMapping:
    @pytest.fixture
    def data_dir(self, tmp_path):
        p = tmp_path / "gene_mapping.csv"
        p.write_text(GENE_MAPPING_CSV)
        return str(tmp_path)

    def test_returns_dict_keyed_by_locus_tag(self, data_dir):
        result = load_gene_mapping(data_dir)
        assert set(result.keys()) == {"PMM0001", "PMM0002", "PMM0003"}

    def test_row_has_correct_fields(self, data_dir):
        result = load_gene_mapping(data_dir)
        row = result["PMM0001"]
        assert row["protein_id"] == "WP_011129038.1"
        assert row["gene_names"] == "dnaN"
        assert row["product"] == "DNA polymerase III subunit beta"

    def test_empty_protein_id_row_included(self, data_dir):
        result = load_gene_mapping(data_dir)
        assert "PMM0003" in result

    def test_raises_on_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_gene_mapping(str(tmp_path))


# ─── load_eggnog ─────────────────────────────────────────────────────────────

EGGNOG_TSV = """\
## This is a comment line — should be skipped
## Another comment
#query\tseed_ortholog\tevalue\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tPFAMs
WP_011129038.1\t592.WP001\t1e-50\tL\tDNA polymerase III beta\tdnaN\tGO:0003677,GO:0006260\t\tko:K02313\tPF00712
WP_011129039.1\t592.WP002\t2e-30\tJ\t30S ribosomal protein\trpsT\t\t\t\t
"""


class TestLoadEggnog:
    @pytest.fixture
    def data_dir(self, tmp_path):
        eggnog_dir = tmp_path / "eggnog"
        eggnog_dir.mkdir()
        (eggnog_dir / "MED4.emapper.annotations").write_text(EGGNOG_TSV)
        return str(tmp_path)

    def test_returns_dict_keyed_by_query(self, data_dir):
        result = load_eggnog(data_dir, "MED4")
        assert set(result.keys()) == {"WP_011129038.1", "WP_011129039.1"}

    def test_strips_hash_from_query_col(self, data_dir):
        result = load_eggnog(data_dir, "MED4")
        row = result["WP_011129038.1"]
        # Column name should be 'query', not '#query'
        assert "query" in row
        assert "#query" not in row

    def test_correct_field_values(self, data_dir):
        result = load_eggnog(data_dir, "MED4")
        row = result["WP_011129038.1"]
        assert row["COG_category"] == "L"
        assert row["Preferred_name"] == "dnaN"
        assert row["KEGG_ko"] == "ko:K02313"

    def test_skips_comment_lines(self, data_dir):
        # Should not have a row keyed by "## This is..."
        result = load_eggnog(data_dir, "MED4")
        for key in result:
            assert not key.startswith("#")

    def test_returns_empty_dict_when_file_missing(self, tmp_path):
        result = load_eggnog(str(tmp_path), "MISSING_STRAIN")
        assert result == {}


# ─── load_uniprot ─────────────────────────────────────────────────────────────

UNIPROT_COL_DATA = {
    "gene_primary": {"A0001": "dnaN", "A0002": "rpsT"},
    "xref_refseq": {"A0001": "WP_011129038.1", "A0002": "WP_011129039.1"},
    "reviewed": {"A0001": "reviewed", "A0002": "unreviewed"},
    "cc_function": {"A0001": "FUNCTION: Catalyzes DNA repair"},
}


class TestLoadUniprot:
    def _write_uniprot(self, path: Path, data=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data or UNIPROT_COL_DATA))

    def test_indexes_by_refseq(self, tmp_path):
        json_path = tmp_path / "uniprot_preprocess_data.json"
        self._write_uniprot(json_path)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), None, "Prochlorococcus")
        assert "WP_011129038.1" in result
        assert "WP_011129039.1" in result

    def test_correct_field_values(self, tmp_path):
        json_path = tmp_path / "uniprot_preprocess_data.json"
        self._write_uniprot(json_path)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), None, "Prochlorococcus")
        row = result["WP_011129038.1"]
        assert row["gene_primary"] == "dnaN"
        assert row["reviewed"] == "reviewed"

    def test_taxid_keyed_path_preferred(self, tmp_path):
        # Put data in taxid-keyed path
        taxid_path = tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919" / "uniprot_preprocess_data.json"
        self._write_uniprot(taxid_path)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedir", 59919, "Prochlorococcus")
        assert "WP_011129038.1" in result

    def test_falls_back_to_project_root(self, tmp_path):
        # No taxid path, but root path exists
        root_path = tmp_path / "uniprot_preprocess_data.json"
        self._write_uniprot(root_path)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedata", 99999, "Prochlorococcus")
        assert "WP_011129038.1" in result

    def test_returns_empty_when_no_file_found(self, tmp_path):
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedata", 99999, "Prochlorococcus")
        assert result == {}

    def test_handles_list_xref_refseq(self, tmp_path):
        data = {
            "xref_refseq": {"A0001": ["WP_111.1", "WP_222.2"]},
            "gene_primary": {"A0001": "test"},
        }
        json_path = tmp_path / "uniprot_preprocess_data.json"
        self._write_uniprot(json_path, data)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), None, "Prochlorococcus")
        assert "WP_111.1" in result
        assert "WP_222.2" in result

    def test_handles_semicolon_separated_xref_refseq(self, tmp_path):
        data = {
            "xref_refseq": {"A0001": "WP_111.1; WP_222.2"},
            "gene_primary": {"A0001": "test"},
        }
        json_path = tmp_path / "uniprot_preprocess_data.json"
        self._write_uniprot(json_path, data)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), None, "Prochlorococcus")
        assert "WP_111.1" in result
        assert "WP_222.2" in result


# ─── AnnotationBuilder ────────────────────────────────────────────────────────

# Minimal gene_mapping row
GM = {
    "locus_tag": "PMM0001",
    "protein_id": "WP_011129038.1",
    "gene_names": "dnaN repA",
    "gene_names_cyanorak": "dnaN_cy",
    "product": "DNA polymerase III subunit beta",
    "product_cyanorak": "DNA pol III beta (Cyanorak)",
    "start": "1",
    "end": "1200",
    "strand": "+",
    "go_process": "DNA replication|0006260||IEA",
    "go_function": "",
    "go_component": "",
    "ec_numbers": "2.7.7.7",
    "kegg": "2.7.7.7",
    "Ontology_term_ncbi": "TIGR001",
    "Ontology_term_cyanorak": "",
    "protein_domains": "TIGR00663,PF00712",
    "eggNOG": "COG1234",
}

EG = {
    "query": "WP_011129038.1",
    "COG_category": "L",
    "Description": "DNA polymerase III subunit beta",
    "Preferred_name": "dnaN",
    "GOs": "GO:0006260,GO:0006261",
    "KEGG_ko": "ko:K02313",
    "PFAMs": "PF00712",
    "EC": "2.7.7.7",
    "eggNOG_OGs": "COG1234@1|root",
    "evalue": "1e-50",
}

UP = {
    "gene_primary": "dnaN",
    "gene_names": "dnaN repA",
    "reviewed": "reviewed",
    "cc_function": "FUNCTION: Catalyzes DNA replication",
    "go_c_id": "0005737",
    "go_p_id": "0006260",
    "go_f_id": "0003677",
    "xref_pfam": "PF00712",
    "xref_refseq": "WP_011129038.1",
    "keywordid": "KW-0067;KW-0131;KW-0233",
}

MINIMAL_CONFIG = {
    "fields": {
        "locus_tag": {
            "type": "passthrough",
            "source": "gene_mapping",
            "field": "locus_tag",
        },
        "gene_name": {
            "type": "single",
            "track_source": "gene_name_source",
            "candidates": [
                {"source": "gene_mapping", "field": "gene_names_cyanorak",
                 "transform": "first_token_space", "source_label": "cyanorak"},
                {"source": "uniprot",      "field": "gene_primary",
                 "source_label": "uniprot"},
                {"source": "gene_mapping", "field": "gene_names",
                 "transform": "first_token_space", "source_label": "ncbi"},
                {"source": "eggnog",       "field": "Preferred_name",
                 "source_label": "eggnog"},
            ],
        },
        "gene_synonyms": {
            "type": "union",
            "delimiter": " ",
            "sources": [
                {"source": "gene_mapping", "field": "gene_names_cyanorak", "delimiter": " "},
                {"source": "gene_mapping", "field": "gene_names",          "delimiter": " "},
                {"source": "uniprot",      "field": "gene_names",          "delimiter": " "},
            ],
        },
        "product": {
            "type": "single",
            "track_source": "product_source",
            "candidates": [
                {"source": "gene_mapping", "field": "product_cyanorak", "source_label": "cyanorak"},
                {"source": "gene_mapping", "field": "product",          "source_label": "ncbi"},
                {"source": "eggnog",       "field": "Description",      "source_label": "eggnog"},
            ],
        },
        "start": {"type": "integer", "source": "gene_mapping", "field": "start"},
        "end":   {"type": "integer", "source": "gene_mapping", "field": "end"},
        "go_terms": {
            "type": "union",
            "filter": "^GO:",
            "sources": [
                {"source": "uniprot",      "field": "go_p_id",     "transform": "add_go_prefix"},
                {"source": "uniprot",      "field": "go_f_id",     "transform": "add_go_prefix"},
                {"source": "uniprot",      "field": "go_c_id",     "transform": "add_go_prefix"},
                {"source": "gene_mapping", "field": "go_process",  "delimiter": ",", "transform": "extract_go_from_pipe"},
                {"source": "eggnog",       "field": "GOs",         "delimiter": ","},
            ],
        },
        "ontology_terms": {
            "type": "union",
            "filter_not": "^GO:",
            "sources": [
                {"source": "gene_mapping", "field": "Ontology_term_ncbi",     "delimiter": ","},
                {"source": "gene_mapping", "field": "Ontology_term_cyanorak", "delimiter": ","},
            ],
        },
        "ec_numbers": {
            "type": "union",
            "sources": [
                {"source": "gene_mapping", "field": "ec_numbers", "delimiter": ","},
                {"source": "gene_mapping", "field": "kegg",       "delimiter": ","},
                {"source": "eggnog",       "field": "EC",         "delimiter": ","},
            ],
        },
        "kegg_ko": {
            "type": "union",
            "sources": [
                {"source": "eggnog", "field": "KEGG_ko", "delimiter": ",",
                 "transform": "strip_prefix_ko"},
            ],
        },
        "pfam_ids": {
            "type": "union",
            "sources": [
                {"source": "gene_mapping", "field": "protein_domains", "delimiter": ",",
                 "transform": "extract_pfam_ids"},
                {"source": "uniprot",      "field": "xref_pfam",       "delimiter": ";"},
                {"source": "eggnog",       "field": "PFAMs",           "delimiter": ","},
            ],
        },
        "seed_ortholog_evalue": {
            "type": "float",
            "source": "eggnog",
            "field": "evalue",
        },
        "cog_category": {
            "type": "passthrough",
            "source": "eggnog",
            "field": "COG_category",
        },
        "function_description": {
            "type": "single",
            "track_source": "function_description_source",
            "candidates": [
                {"source": "uniprot", "field": "cc_function",
                 "transform": "strip_function_prefix", "source_label": "uniprot"},
                {"source": "eggnog",  "field": "Description", "source_label": "eggnog"},
            ],
        },
        "keyword_ids": {
            "type": "passthrough_list",
            "source": "uniprot",
            "field": "keywordid",
            "delimiter": ";",
        },
    }
}


class TestAnnotationBuilderBuildWide:
    def setup_method(self):
        self.builder = AnnotationBuilder(MINIMAL_CONFIG)

    def test_prefixes_gene_mapping_fields(self):
        wide = self.builder.build_wide(GM, {}, {})
        assert "gene_mapping_locus_tag" in wide
        assert wide["gene_mapping_locus_tag"] == "PMM0001"

    def test_prefixes_eggnog_fields(self):
        wide = self.builder.build_wide({}, EG, {})
        assert "eggnog_COG_category" in wide
        assert wide["eggnog_COG_category"] == "L"

    def test_prefixes_uniprot_fields(self):
        wide = self.builder.build_wide({}, {}, UP)
        assert "uniprot_gene_primary" in wide
        assert wide["uniprot_gene_primary"] == "dnaN"

    def test_excludes_empty_values(self):
        gm = {"locus_tag": "PMM0001", "product": ""}
        wide = self.builder.build_wide(gm, {}, {})
        # Empty product should be excluded
        assert "gene_mapping_product" not in wide
        assert "gene_mapping_locus_tag" in wide

    def test_all_sources_combined(self):
        wide = self.builder.build_wide(GM, EG, UP)
        assert "gene_mapping_locus_tag" in wide
        assert "eggnog_COG_category" in wide
        assert "uniprot_gene_primary" in wide


class TestAnnotationBuilderBuildMerged:
    def setup_method(self):
        self.builder = AnnotationBuilder(MINIMAL_CONFIG)

    # ── passthrough ──────────────────────────────────────────────────────────

    def test_passthrough_basic(self):
        merged = self.builder.build_merged(GM, {}, {})
        assert merged["locus_tag"] == "PMM0001"

    def test_passthrough_missing_field_omitted(self):
        merged = self.builder.build_merged({}, {}, {})
        assert "locus_tag" not in merged

    def test_passthrough_url_decodes_value(self):
        gm = dict(GM, product_cyanorak="ATP%2FGTPase")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["product"] == "ATP/GTPase"

    # ── integer / float ──────────────────────────────────────────────────────

    def test_integer_type(self):
        merged = self.builder.build_merged(GM, {}, {})
        assert merged["start"] == 1
        assert isinstance(merged["start"], int)

    def test_integer_from_float_string(self):
        gm = dict(GM, start="1.0")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["start"] == 1

    def test_integer_invalid_omitted(self):
        gm = dict(GM, start="notanumber", end="")
        merged = self.builder.build_merged(gm, {}, {})
        assert "start" not in merged

    def test_float_type(self):
        merged = self.builder.build_merged({}, EG, {})
        assert abs(merged["seed_ortholog_evalue"] - 1e-50) < 1e-60
        assert isinstance(merged["seed_ortholog_evalue"], float)

    def test_float_invalid_omitted(self):
        eg = dict(EG, evalue="bad")
        merged = self.builder.build_merged({}, eg, {})
        assert "seed_ortholog_evalue" not in merged

    # ── single ────────────────────────────────────────────────────────────────

    def test_single_first_candidate_wins(self):
        # gene_names_cyanorak present → should win over uniprot and eggnog
        merged = self.builder.build_merged(GM, EG, UP)
        assert merged["gene_name"] == "dnaN_cy"

    def test_single_falls_through_to_next_candidate(self):
        # No cyanorak name → uniprot wins
        gm = dict(GM, gene_names_cyanorak="")
        merged = self.builder.build_merged(gm, {}, UP)
        assert merged["gene_name"] == "dnaN"
        assert merged["gene_name_source"] == "uniprot"

    def test_single_eggnog_fallback(self):
        gm = dict(GM, gene_names_cyanorak="", gene_names="")
        merged = self.builder.build_merged(gm, EG, {})
        assert merged["gene_name"] == "dnaN"
        assert merged["gene_name_source"] == "eggnog"

    def test_single_source_label_recorded(self):
        merged = self.builder.build_merged(GM, {}, {})
        assert merged.get("gene_name_source") == "cyanorak"

    def test_single_strip_function_prefix_transform(self):
        merged = self.builder.build_merged({}, {}, UP)
        assert merged.get("function_description") == "Catalyzes DNA replication"
        assert merged.get("function_description_source") == "uniprot"

    # ── union ─────────────────────────────────────────────────────────────────

    def test_union_deduplicates(self):
        # go_p_id and eggnog GOs both have GO:0006260
        merged = self.builder.build_merged(GM, EG, UP)
        go_terms = merged.get("go_terms", [])
        assert go_terms.count("GO:0006260") == 1

    def test_union_filter_keeps_only_go_terms(self):
        merged = self.builder.build_merged(GM, EG, UP)
        for term in merged.get("go_terms", []):
            assert term.startswith("GO:")

    def test_union_filter_not_excludes_go_terms(self):
        merged = self.builder.build_merged(GM, EG, UP)
        for term in merged.get("ontology_terms", []):
            assert not term.startswith("GO:")

    def test_union_add_go_prefix_transform(self):
        # UniProt go_p_id is bare 7-digit → should become GO:0006260
        merged = self.builder.build_merged({}, {}, UP)
        go_terms = merged.get("go_terms", [])
        assert "GO:0006260" in go_terms

    def test_union_strip_prefix_ko_transform(self):
        merged = self.builder.build_merged({}, EG, {})
        assert merged.get("kegg_ko") == ["K02313"]

    def test_union_extract_pfam_ids_transform(self):
        merged = self.builder.build_merged(GM, EG, UP)
        pfams = merged.get("pfam_ids", [])
        assert "PF00712" in pfams
        # TIGR00663 should be filtered out
        assert "TIGR00663" not in pfams

    def test_union_ec_numbers_deduplicated(self):
        # ec_numbers and kegg both carry 2.7.7.7; eggnog too
        merged = self.builder.build_merged(GM, EG, {})
        assert merged.get("ec_numbers", []).count("2.7.7.7") == 1

    def test_union_empty_all_sources_omitted(self):
        gm = dict(GM, Ontology_term_ncbi="", Ontology_term_cyanorak="")
        merged = self.builder.build_merged(gm, {}, {})
        assert "ontology_terms" not in merged

    # ── passthrough_list ──────────────────────────────────────────────────────

    def test_passthrough_list_splits_delimited_string(self):
        # keywordid = "KW-0067;KW-0131;KW-0233" → list of 3 tokens
        merged = self.builder.build_merged(GM, EG, UP)
        kw = merged.get("keyword_ids", [])
        assert kw == ["KW-0067", "KW-0131", "KW-0233"]

    def test_passthrough_list_pre_existing_list_passes_through(self):
        up = dict(UP, keywordid=["KW-0067", "KW-0131"])
        merged = self.builder.build_merged(GM, EG, up)
        assert merged["keyword_ids"] == ["KW-0067", "KW-0131"]

    def test_passthrough_list_missing_field_omitted(self):
        up = {k: v for k, v in UP.items() if k != "keywordid"}
        merged = self.builder.build_merged(GM, EG, up)
        assert "keyword_ids" not in merged

    # ── gene_synonyms dedup ───────────────────────────────────────────────────

    def test_gene_synonyms_excludes_canonical_name(self):
        # gene_name = "dnaN_cy"; gene_synonyms union includes "dnaN_cy" → should be removed
        merged = self.builder.build_merged(GM, {}, {})
        gene_name = merged.get("gene_name")
        synonyms = merged.get("gene_synonyms", [])
        assert gene_name not in synonyms

    def test_gene_synonyms_all_same_as_gene_name_removes_field(self):
        # Only one distinct name across all sources → synonyms removed entirely
        gm = {
            "gene_names_cyanorak": "dnaN",
            "gene_names": "dnaN",
        }
        merged = self.builder.build_merged(gm, {}, {"gene_names": "dnaN"})
        # gene_name = "dnaN"; all synonyms are "dnaN" → field should be absent
        assert "gene_synonyms" not in merged

    # ── annotation_quality ────────────────────────────────────────────────────

    def test_quality_3_when_uniprot_reviewed(self):
        merged = self.builder.build_merged(GM, EG, dict(UP, reviewed="reviewed"))
        assert merged["annotation_quality"] == 3

    def test_quality_2_when_cyanorak_product(self):
        up_no_review = {k: v for k, v in UP.items() if k != "reviewed"}
        gm = dict(GM, product_cyanorak="some product")
        merged = self.builder.build_merged(gm, {}, up_no_review)
        assert merged["annotation_quality"] == 2

    def test_quality_2_when_ncbi_product(self):
        gm = dict(GM, product_cyanorak="", product="ncbi product")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 2

    def test_quality_1_when_only_eggnog(self):
        gm = dict(GM, product_cyanorak="", product="")
        merged = self.builder.build_merged(gm, EG, {})
        assert merged["annotation_quality"] == 1

    def test_quality_0_when_nothing(self):
        gm = {k: "" for k in GM}
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0


# ─── process_strain ───────────────────────────────────────────────────────────

GENE_MAPPING_SIMPLE = """\
locus_tag,protein_id,gene_names,product,product_cyanorak,start,end,strand
PMM0001,WP_001.1,dnaN,DNA pol III beta,DNA pol III beta (Cyanorak),1,1200,+
PMM0002,WP_002.1,rpsT,,30S ribosomal S20,1300,1600,+
"""


class TestProcessStrain:
    @pytest.fixture
    def data_dir(self, tmp_path):
        d = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MED4"
        d.mkdir(parents=True)
        (d / "gene_mapping.csv").write_text(GENE_MAPPING_SIMPLE)
        return d

    @pytest.fixture
    def row(self, data_dir):
        return {
            "strain_name": "MED4",
            "data_dir": str(data_dir),
            "ncbi_taxon_id": "",  # No UniProt for this test
        }

    def test_creates_wide_and_merged_json(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        assert (data_dir / "gene_annotations_wide.json").exists()
        assert (data_dir / "gene_annotations_merged.json").exists()

    def test_merged_keyed_by_locus_tag(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        assert set(merged.keys()) == {"PMM0001", "PMM0002"}

    def test_wide_keyed_by_locus_tag(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        wide = json.loads((data_dir / "gene_annotations_wide.json").read_text())
        assert set(wide.keys()) == {"PMM0001", "PMM0002"}

    def test_skips_when_merged_exists_and_no_force(self, row, data_dir, capsys):
        merged_path = data_dir / "gene_annotations_merged.json"
        merged_path.write_text('{"PMM0001": {"locus_tag": "PMM0001"}}')
        mtime_before = merged_path.stat().st_mtime

        process_strain(row, MINIMAL_CONFIG, force=False)

        assert merged_path.stat().st_mtime == mtime_before
        out = capsys.readouterr().out
        assert "Skipping" in out

    def test_force_overwrites_existing_merged(self, row, data_dir):
        merged_path = data_dir / "gene_annotations_merged.json"
        merged_path.write_text('{"old": true}')

        process_strain(row, MINIMAL_CONFIG, force=True)

        merged = json.loads(merged_path.read_text())
        assert "old" not in merged
        assert "PMM0001" in merged

    def test_annotation_quality_in_merged(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        # PMM0001 has product_cyanorak → quality 2
        assert merged["PMM0001"]["annotation_quality"] == 2
        # PMM0002 has product_cyanorak → quality 2 as well
        assert merged["PMM0002"]["annotation_quality"] == 2

    def test_eggnog_rows_joined_via_protein_id(self, row, data_dir):
        eggnog_dir = data_dir / "eggnog"
        eggnog_dir.mkdir()
        (eggnog_dir / "MED4.emapper.annotations").write_text(
            "## comment\n"
            "#query\tCOG_category\tPreferred_name\tGOs\tKEGG_ko\tEC\tPFAMs\tevalue\n"
            "WP_001.1\tL\tdnaN\tGO:0003677\tko:K02313\t2.7.7.7\tPF00712\t1e-50\n"
        )
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        # PMM0001 has protein_id WP_001.1 which matches eggnog query
        pmm1 = merged["PMM0001"]
        assert pmm1.get("cog_category") == "L"
        assert "K02313" in pmm1.get("kegg_ko", [])

    def test_no_uniprot_load_when_taxon_id_empty(self, row, data_dir, capsys):
        # ncbi_taxon_id="" → up_data should be empty, no error
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        assert "PMM0001" in merged  # completed without error
