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
    load_eggnog,
    load_gene_mapping,
    load_uniprot,
    process_strain,
)
from multiomics_kg.download.utils.annotation_helpers import (
    _coerce_to_tokens,
    _nonempty,
    _split,
    extract_first_match_in_sources,
)
from multiomics_kg.download.utils.annotation_transforms import (
    _resolve_ec_chain,
    _tx_add_go_prefix,
    _tx_extract_go_from_pipe,
    _tx_extract_pfam_ids,
    _tx_extract_pfam_names,
    _tx_first_token_space,
    _tx_normalize_ec,
    _tx_split_cog_category,
    _tx_strip_function_prefix,
    _tx_strip_prefix_ko,
)
import multiomics_kg.download.utils.annotation_transforms as _at
from multiomics_kg.download.utils.paths import infer_organism_group


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

PROTEIN_ANNOTATIONS_DATA = {
    "A0001": {"gene_symbol": "dnaN", "refseq_ids": ["WP_011129038.1"], "is_reviewed": "reviewed", "function_description": "Catalyzes DNA repair"},
    "A0002": {"gene_symbol": "rpsT", "refseq_ids": ["WP_011129039.1"], "is_reviewed": "not reviewed"},
}


class TestLoadUniprot:
    def _write_uniprot(self, path: Path, data=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data or PROTEIN_ANNOTATIONS_DATA))

    def _taxid_path(self, tmp_path, taxid=59919):
        return tmp_path / "cache" / "data" / "Prochlorococcus" / "uniprot" / str(taxid) / "protein_annotations.json"

    def test_indexes_by_refseq(self, tmp_path):
        self._write_uniprot(self._taxid_path(tmp_path))
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), 59919, "Prochlorococcus")
        assert "WP_011129038.1" in result
        assert "WP_011129039.1" in result

    def test_correct_field_values(self, tmp_path):
        self._write_uniprot(self._taxid_path(tmp_path))
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), 59919, "Prochlorococcus")
        row = result["WP_011129038.1"]
        assert row["gene_symbol"] == "dnaN"
        assert row["is_reviewed"] == "reviewed"

    def test_taxid_keyed_path(self, tmp_path):
        self._write_uniprot(self._taxid_path(tmp_path))
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedir", 59919, "Prochlorococcus")
        assert "WP_011129038.1" in result

    def test_returns_empty_when_no_taxid(self, tmp_path):
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedata", None, "Prochlorococcus")
        assert result == {}

    def test_returns_empty_when_no_file_found(self, tmp_path):
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot("somedata", 99999, "Prochlorococcus")
        assert result == {}

    def test_handles_multiple_refseq_ids(self, tmp_path):
        data = {
            "A0001": {"gene_symbol": "test", "refseq_ids": ["WP_111.1", "WP_222.2"]},
        }
        self._write_uniprot(self._taxid_path(tmp_path), data)
        with patch("multiomics_kg.download.build_gene_annotations.PROJECT_ROOT", tmp_path):
            result = load_uniprot(str(tmp_path), 59919, "Prochlorococcus")
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
    # Scoring: 0=hypothetical-no-func, 1=hypothetical-with-func,
    #          2=real-product, 3=real-product + ≥2 structured annotations

    def test_quality_3_real_product_with_structured(self):
        # GM has real product + go_terms + ec_numbers + pfam_ids → quality 3
        merged = self.builder.build_merged(GM, EG, UP)
        assert merged["annotation_quality"] == 3

    def test_quality_2_real_product_no_structured(self):
        # Real product but no structured annotations (no GO, KEGG, EC, Pfam sources)
        gm = dict(GM, product_cyanorak="some product", product="",
                  go_process="", go_function="", go_component="",
                  ec_numbers="", kegg="", protein_domains="", eggNOG="")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 2

    def test_quality_2_real_product_one_structured(self):
        # Real product + only 1 structured field (ec_numbers) → quality 2 not 3
        gm = dict(GM, product_cyanorak="some product", product="",
                  go_process="", go_function="", go_component="",
                  kegg="", protein_domains="", eggNOG="")
        # ec_numbers="2.7.7.7" still in GM → 1 structured field
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 2

    def test_quality_1_hypothetical_with_func(self):
        # Hypothetical product but has function_description from UniProt
        gm = dict(GM, product_cyanorak="", product="hypothetical protein",
                  go_process="", go_function="", go_component="",
                  ec_numbers="", kegg="", protein_domains="", eggNOG="")
        merged = self.builder.build_merged(gm, {}, UP)
        assert merged["annotation_quality"] == 1

    def test_quality_0_hypothetical_no_func(self):
        # Hypothetical product, no function_description
        gm = dict(GM, product_cyanorak="", product="hypothetical protein",
                  go_process="", go_function="", go_component="",
                  ec_numbers="", kegg="", protein_domains="", eggNOG="")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0

    def test_quality_0_when_nothing(self):
        gm = {k: "" for k in GM}
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0

    def test_quality_0_uncharacterized_counts_as_hypothetical(self):
        gm = dict(GM, product_cyanorak="", product="Uncharacterized protein",
                  go_process="", go_function="", go_component="",
                  ec_numbers="", kegg="", protein_domains="", eggNOG="")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0

    def test_quality_0_conserved_hypothetical(self):
        gm = dict(GM, product_cyanorak="", product="conserved hypothetical protein",
                  go_process="", go_function="", go_component="",
                  ec_numbers="", kegg="", protein_domains="", eggNOG="")
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0

    def test_quality_0_hypothetical_with_structured_but_no_func(self):
        # Hypothetical product with GO + EC (structured annotations) but no
        # function_description → still quality 0, not 3.
        # Structured annotations only boost quality for non-hypothetical products.
        gm = dict(GM, product_cyanorak="", product="hypothetical protein")
        # GM has go_process + ec_numbers → structured fields resolve, but
        # product is hypothetical and no func_desc → quality 0
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["annotation_quality"] == 0

    # ── organism_strain ──────────────────────────────────────────────────────

    def test_organism_strain_set_from_param(self):
        merged = self.builder.build_merged(GM, {}, {}, organism_name="Prochlorococcus MED4")
        assert merged["organism_strain"] == "Prochlorococcus MED4"

    def test_organism_strain_absent_when_no_param(self):
        merged = self.builder.build_merged(GM, {}, {})
        assert "organism_strain" not in merged

    def test_organism_strain_absent_when_empty(self):
        merged = self.builder.build_merged(GM, {}, {}, organism_name="")
        assert "organism_strain" not in merged

    # ── gene_summary ─────────────────────────────────────────────────────────

    def test_gene_summary_full(self):
        merged = self.builder.build_merged(GM, EG, UP)
        summary = merged.get("gene_summary", "")
        assert "dnaN_cy" in summary  # gene_name (cyanorak wins)
        assert "::" in summary

    def test_gene_summary_no_function_description(self):
        # No UniProt, no eggnog → summary is gene_name :: product
        merged = self.builder.build_merged(GM, {}, {})
        summary = merged.get("gene_summary", "")
        assert "dnaN_cy" in summary
        assert "DNA pol III beta (Cyanorak)" in summary

    def test_gene_summary_uses_separator_not_pipe(self):
        merged = self.builder.build_merged(GM, EG, UP)
        summary = merged.get("gene_summary", "")
        assert "|" not in summary
        assert "::" in summary

    def test_gene_summary_deduplicates_product_and_description(self):
        # When eggnog description == product, description should not repeat
        gm = dict(GM, product_cyanorak="", product="DNA polymerase III subunit beta")
        eg = dict(EG, Description="DNA polymerase III subunit beta")
        merged = self.builder.build_merged(gm, eg, {})
        summary = merged.get("gene_summary", "")
        # "DNA polymerase III subunit beta" should appear only once
        assert summary.count("DNA polymerase III subunit beta") == 1

    def test_gene_summary_fallback_to_eggnog(self):
        # No UniProt function → falls back to eggnog description
        gm = dict(GM, product_cyanorak="", product="hypothetical protein")
        merged = self.builder.build_merged(gm, EG, {})
        summary = merged.get("gene_summary", "")
        assert "DNA polymerase III subunit beta" in summary  # eggnog Description

    def test_gene_summary_minimal_gene_only(self):
        # No product, no function → just locus_tag from gene_name
        gm = dict(GM, gene_names_cyanorak="", gene_names="", product_cyanorak="", product="")
        merged = self.builder.build_merged(gm, {}, {})
        summary = merged.get("gene_summary")
        # gene_name may be absent; if so, no summary at all
        if summary:
            assert "::" not in summary or len(summary.split("::")) <= 2

    def test_gene_summary_skips_locus_tag_gene_name(self):
        # gene_name == locus_tag → should not appear in summary
        gm = dict(GM, gene_names_cyanorak="", gene_names="PMM0001",
                  product_cyanorak="some product")
        merged = self.builder.build_merged(gm, {}, {})
        summary = merged.get("gene_summary", "")
        assert "PMM0001" not in summary
        assert "some product" in summary

    def test_gene_summary_skips_refseq_rs_gene_name(self):
        # RS-pattern gene_name → excluded
        gm = dict(GM, gene_names_cyanorak="", gene_names="ALTBGP6_RS00025",
                  locus_tag="ALTBGP6_RS00025", product_cyanorak="some toxin")
        merged = self.builder.build_merged(gm, {}, {})
        summary = merged.get("gene_summary", "")
        assert "ALTBGP6_RS00025" not in summary
        assert "some toxin" in summary

    def test_gene_summary_skips_synw_style_gene_name(self):
        # SYNW1033-style identifier → excluded
        gm = dict(GM, gene_names_cyanorak="", gene_names="SYNW1033",
                  product_cyanorak="photosystem II protein")
        merged = self.builder.build_merged(gm, {}, {})
        summary = merged.get("gene_summary", "")
        assert "SYNW1033" not in summary
        assert "photosystem II protein" in summary

    def test_gene_summary_skips_mit1002_style_gene_name(self):
        # MIT1002_00123-style identifier → excluded
        gm = dict(GM, gene_names_cyanorak="", gene_names="MIT1002_00123",
                  product_cyanorak="ABC transporter")
        merged = self.builder.build_merged(gm, {}, {})
        summary = merged.get("gene_summary", "")
        assert "MIT1002_00123" not in summary
        assert "ABC transporter" in summary

    def test_gene_summary_keeps_real_gene_name(self):
        # Real biological names like dnaN, rpsT, petB → kept
        for name in ["dnaN", "rpsT", "petB", "rlmB"]:
            gm = dict(GM, gene_names_cyanorak=name, product_cyanorak="some product")
            merged = self.builder.build_merged(gm, {}, {})
            summary = merged.get("gene_summary", "")
            assert name in summary, f"{name} should be in summary"

    def test_gene_summary_skips_duf_description(self):
        # "Protein of unknown function (DUF3464)" → excluded from summary
        gm = dict(GM, gene_names_cyanorak="", gene_names="",
                  product_cyanorak="", product="some membrane protein")
        eg = dict(EG, Description="Protein of unknown function (DUF3464)")
        merged = self.builder.build_merged(gm, eg, {})
        summary = merged.get("gene_summary", "")
        assert "unknown function" not in summary
        assert "some membrane protein" in summary

    def test_gene_summary_skips_domain_of_unknown_function(self):
        # "Domain of unknown function" variant → also excluded
        gm = dict(GM, gene_names_cyanorak="", gene_names="",
                  product_cyanorak="", product="some enzyme")
        eg = dict(EG, Description="Domain of unknown function (DUF1234)")
        merged = self.builder.build_merged(gm, eg, {})
        summary = merged.get("gene_summary", "")
        assert "unknown function" not in summary
        assert "some enzyme" in summary

    # ── all_identifiers ──────────────────────────────────────────────────────

    def test_all_identifiers_absent_with_minimal_config(self):
        # MINIMAL_CONFIG doesn't resolve protein_id/locus_tag_ncbi/gene_name_synonyms/
        # alternative_locus_tags, so all_identifiers should be absent
        merged = self.builder.build_merged(GM, EG, UP)
        assert "all_identifiers" not in merged

    def test_all_identifiers_excludes_locus_tag(self):
        merged = self.builder.build_merged(GM, EG, UP)
        ids = merged.get("all_identifiers", [])
        assert "PMM0001" not in ids  # locus_tag is scalar-indexed

    def test_all_identifiers_excludes_gene_name(self):
        merged = self.builder.build_merged(GM, EG, UP)
        ids = merged.get("all_identifiers", [])
        gene_name = merged.get("gene_name")
        assert gene_name not in ids  # gene_name is scalar-indexed

    def test_all_identifiers_sorted(self):
        merged = self.builder.build_merged(GM, EG, UP)
        ids = merged.get("all_identifiers", [])
        if ids:
            assert ids == sorted(ids)

    def test_all_identifiers_absent_when_no_alt_ids(self):
        gm = {"locus_tag": "PMM9999", "gene_names_cyanorak": "geneX"}
        merged = self.builder.build_merged(gm, {}, {})
        # gene_name = "geneX", only synonym is "geneX" (same → removed)
        # No protein_id, no locus_tag_ncbi, no other synonyms → no all_identifiers
        assert "all_identifiers" not in merged


# ─── all_identifiers with production config ──────────────────────────────────


class TestAllIdentifiersFullConfig:
    """Test all_identifiers using the production gene_annotations_config.yaml."""

    def setup_method(self):
        from multiomics_kg.download.utils.cli import load_config
        config = load_config("config/gene_annotations_config.yaml")
        self.builder = AnnotationBuilder(config)

    def test_includes_protein_id(self):
        gm = dict(GM, protein_id="WP_011129038.1", locus_tag_ncbi="TX50_RS00020")
        merged = self.builder.build_merged(gm, EG, UP)
        ids = merged.get("all_identifiers", [])
        assert "WP_011129038.1" in ids

    def test_includes_locus_tag_ncbi(self):
        gm = dict(GM, protein_id="WP_011129038.1", locus_tag_ncbi="TX50_RS00020")
        merged = self.builder.build_merged(gm, EG, UP)
        ids = merged.get("all_identifiers", [])
        assert "TX50_RS00020" in ids

    def test_includes_old_locus_tags(self):
        gm = dict(GM, old_locus_tags="PMM0001_old,PMM0001_v2")
        merged = self.builder.build_merged(gm, {}, {})
        ids = merged.get("all_identifiers", [])
        assert "PMM0001_old" in ids
        assert "PMM0001_v2" in ids

    def test_excludes_locus_tag_and_gene_name(self):
        gm = dict(GM, protein_id="WP_011129038.1", locus_tag_ncbi="TX50_RS00020")
        merged = self.builder.build_merged(gm, EG, UP)
        ids = merged.get("all_identifiers", [])
        assert merged["locus_tag"] not in ids
        assert merged["gene_name"] not in ids

    def test_is_sorted(self):
        gm = dict(GM, protein_id="WP_011129038.1", locus_tag_ncbi="TX50_RS00020",
                   old_locus_tags="ZZZ_old,AAA_old")
        merged = self.builder.build_merged(gm, EG, UP)
        ids = merged.get("all_identifiers", [])
        assert ids == sorted(ids)


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
            "preferred_name": "Prochlorococcus MED4",
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
        # PMM0001 has real product, no structured annotations → quality 2
        assert merged["PMM0001"]["annotation_quality"] == 2
        # PMM0002 has real product, no structured annotations → quality 2
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

    def test_organism_strain_from_preferred_name(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        assert merged["PMM0001"]["organism_strain"] == "Prochlorococcus MED4"
        assert merged["PMM0002"]["organism_strain"] == "Prochlorococcus MED4"

    def test_organism_strain_falls_back_to_strain_name(self, data_dir):
        row_no_pref = {
            "strain_name": "MED4",
            "data_dir": str(data_dir),
            "ncbi_taxon_id": "",
        }
        process_strain(row_no_pref, MINIMAL_CONFIG, force=True)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        assert merged["PMM0001"]["organism_strain"] == "MED4"

    def test_gene_summary_in_merged(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        summary = merged["PMM0001"].get("gene_summary", "")
        assert "dnaN" in summary
        assert "::" in summary

    def test_all_identifiers_in_merged(self, row, data_dir):
        process_strain(row, MINIMAL_CONFIG, force=False)
        merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
        pmm1 = merged["PMM0001"]
        # all_identifiers is built from gene_synonyms (minus gene_name) + other alt ID fields;
        # MINIMAL_CONFIG resolves gene_synonyms from gene_names_cyanorak + gene_names;
        # gene_names="dnaN" in GENE_MAPPING_SIMPLE → gene_name="dnaN" (first token),
        # gene_synonyms=["dnaN", "repA"...] minus gene_name
        # Depending on config fields, all_identifiers may or may not be populated
        ids = pmm1.get("all_identifiers")
        if ids is not None:
            assert isinstance(ids, list)
            assert pmm1["locus_tag"] not in ids  # locus_tag excluded


# ─── Fix D: kegg_pathway filter ───────────────────────────────────────────────

class TestKeggPathwayFilter:
    """_resolve_union filter drops map* pathway entries, keeps only ko*."""

    def _make_config(self, filter_val=None):
        cfg: dict = {
            "type": "union",
            "sources": [{"source": "eggnog", "field": "KEGG_Pathway", "delimiter": ","}],
        }
        if filter_val:
            cfg["filter"] = filter_val
        return cfg

    def setup_method(self):
        self.builder = AnnotationBuilder({"fields": {}})

    def test_without_filter_includes_map_entries(self):
        eg = {"KEGG_Pathway": "ko03030,map03030,ko00001"}
        result = self.builder._resolve_union(self._make_config(), {}, eg, {})
        assert "map03030" in result

    def test_with_ko_filter_drops_map_entries(self):
        eg = {"KEGG_Pathway": "ko03030,map03030,ko00001"}
        result = self.builder._resolve_union(self._make_config("^ko"), {}, eg, {})
        assert "ko03030" in result
        assert "ko00001" in result
        assert "map03030" not in result

    def test_empty_after_filter_returns_none(self):
        eg = {"KEGG_Pathway": "map03030,map00001"}
        result = self.builder._resolve_union(self._make_config("^ko"), {}, eg, {})
        assert result is None


# ─── Fix B: cog_category transform ────────────────────────────────────────────

class TestTxSplitCogCategory:
    def test_multi_letter(self):
        assert _tx_split_cog_category("LU") == ["L", "U"]

    def test_single_letter(self):
        assert _tx_split_cog_category("S") == ["S"]

    def test_three_letters(self):
        assert _tx_split_cog_category("EGP") == ["E", "G", "P"]

    def test_dash_sentinel(self):
        assert _tx_split_cog_category("-") == []

    def test_empty(self):
        assert _tx_split_cog_category("") == []

    def test_non_string(self):
        assert _tx_split_cog_category(None) == []

    def test_strips_whitespace(self):
        assert _tx_split_cog_category("  LU  ") == ["L", "U"]


class TestCogCategoryBuildMerged:
    """cog_category passthrough+transform produces a list."""

    def setup_method(self):
        config = {
            "fields": {
                "cog_category": {
                    "type": "passthrough",
                    "source": "eggnog",
                    "field": "COG_category",
                    "transform": "split_cog_category",
                }
            }
        }
        self.builder = AnnotationBuilder(config)

    def test_single_letter_becomes_list(self):
        merged = self.builder.build_merged({}, {"COG_category": "L"}, {})
        assert merged["cog_category"] == ["L"]

    def test_multi_letter_split(self):
        merged = self.builder.build_merged({}, {"COG_category": "LU"}, {})
        assert merged["cog_category"] == ["L", "U"]

    def test_dash_omitted(self):
        merged = self.builder.build_merged({}, {"COG_category": "-"}, {})
        assert "cog_category" not in merged

    def test_missing_field_omitted(self):
        merged = self.builder.build_merged({}, {}, {})
        assert "cog_category" not in merged


# ─── Fix C: pfam_names / pfam_ids cleanup ─────────────────────────────────────

class TestTxExtractPfamNames:
    def test_keeps_shortname_tokens(self):
        result = _tx_extract_pfam_names("PF00712,DNA_pol3_beta,DNA_pol3_beta_2")
        assert result == ["DNA_pol3_beta", "DNA_pol3_beta_2"]

    def test_no_shortnames_returns_empty(self):
        result = _tx_extract_pfam_names("PF00712,PF02767")
        assert result == []

    def test_mixed_with_tigr(self):
        result = _tx_extract_pfam_names("TIGR00663,PF00712,HATPase_c")
        assert result == ["TIGR00663", "HATPase_c"]

    def test_non_string_returns_empty(self):
        assert _tx_extract_pfam_names(None) == []

    def test_empty_string(self):
        assert _tx_extract_pfam_names("") == []


class TestPfamNamesBuildMerged:
    """pfam_names field only populated from eggnog shortnames."""

    def setup_method(self):
        config = {
            "fields": {
                "pfam_ids": {
                    "type": "union",
                    "sources": [
                        {"source": "eggnog", "field": "PFAMs", "delimiter": ",",
                         "transform": "extract_pfam_ids"},
                    ],
                },
                "pfam_names": {
                    "type": "union",
                    "sources": [
                        {"source": "eggnog", "field": "PFAMs", "delimiter": ",",
                         "transform": "extract_pfam_names"},
                    ],
                },
            }
        }
        self.builder = AnnotationBuilder(config)

    def test_pfam_ids_only_pf_tokens(self):
        eg = {"PFAMs": "PF00712,DNA_pol3_beta,PF02767"}
        merged = self.builder.build_merged({}, eg, {})
        assert merged["pfam_ids"] == ["PF00712", "PF02767"]
        assert "DNA_pol3_beta" not in merged["pfam_ids"]

    def test_pfam_names_only_shortnames(self):
        eg = {"PFAMs": "PF00712,DNA_pol3_beta,HATPase_c"}
        merged = self.builder.build_merged({}, eg, {})
        assert "DNA_pol3_beta" in merged["pfam_names"]
        assert "HATPase_c" in merged["pfam_names"]
        assert "PF00712" not in merged["pfam_names"]

    def test_pure_pf_ids_no_pfam_names_field(self):
        eg = {"PFAMs": "PF00712,PF02767"}
        merged = self.builder.build_merged({}, eg, {})
        assert "pfam_names" not in merged

    def test_no_eggnog_no_fields(self):
        merged = self.builder.build_merged({}, {}, {})
        assert "pfam_ids" not in merged
        assert "pfam_names" not in merged


# ─── Fix E: gene_synonyms / gene_name_synonyms / alternative_locus_tags ───────

LOCUS_TAG_FILTER = r"^([A-Z][A-Z0-9]*\d{4,}|.+_.+)$"
GENE_NAME_FILTER_NOT = LOCUS_TAG_FILTER


class TestGeneSynonymsDelimiterFix:
    """gene_synonyms with per-source delimiter correctly splits space-separated values."""

    def setup_method(self):
        config = {
            "fields": {
                "gene_name": {
                    "type": "single",
                    "candidates": [
                        {"source": "gene_mapping", "field": "gene_names",
                         "transform": "first_token_space"},
                    ],
                },
                "gene_synonyms": {
                    "type": "union",
                    "sources": [
                        {"source": "gene_mapping", "field": "gene_names", "delimiter": " "},
                    ],
                },
            }
        }
        self.builder = AnnotationBuilder(config)

    def test_space_delimited_splits_correctly(self):
        gm = {"gene_names": "dnaN PMM0001"}
        merged = self.builder.build_merged(gm, {}, {})
        synonyms = merged.get("gene_synonyms", [])
        # Should split into ["dnaN", "PMM0001"], then "dnaN" removed as canonical
        assert "PMM0001" in synonyms
        # dnaN is canonical gene_name → excluded
        assert "dnaN" not in synonyms

    def test_single_name_no_synonyms(self):
        gm = {"gene_names": "dnaN"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "gene_synonyms" not in merged


class TestAlternativeLocusTagsFilter:
    """alternative_locus_tags keeps only locus-tag-like tokens."""

    def setup_method(self):
        config = {
            "fields": {
                "alternative_locus_tags": {
                    "type": "union",
                    "filter": LOCUS_TAG_FILTER,
                    "sources": [
                        {"source": "gene_mapping", "field": "gene_names", "delimiter": " "},
                    ],
                },
            }
        }
        self.builder = AnnotationBuilder(config)

    def test_old_format_locus_tag_kept(self):
        # PMM0001: all-uppercase + 4 trailing digits → locus tag
        gm = {"gene_names": "dnaN PMM0001"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "PMM0001" in merged.get("alternative_locus_tags", [])

    def test_gene_name_excluded(self):
        # dnaN: has only 2 trailing digits at most → gene name, not locus tag
        gm = {"gene_names": "dnaN PMM0001"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "dnaN" not in merged.get("alternative_locus_tags", [])

    def test_rs_format_locus_tag_kept(self):
        # ALTBGP6_RS00025: contains underscore → locus tag
        gm = {"gene_names": "ALTBGP6_RS00025"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "ALTBGP6_RS00025" in merged.get("alternative_locus_tags", [])

    def test_uppercase_gene_name_not_locus_tag(self):
        # ALDH: uppercase but only 0 trailing digits → gene name, not locus tag
        gm = {"gene_names": "ALDH PMM0002"}
        merged = self.builder.build_merged(gm, {}, {})
        alt = merged.get("alternative_locus_tags", [])
        assert "ALDH" not in alt
        assert "PMM0002" in alt

    def test_all_gene_names_no_locus_tags(self):
        gm = {"gene_names": "dnaN purF rpoA"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "alternative_locus_tags" not in merged


class TestGeneNameSynonymsFilterNot:
    """gene_name_synonyms keeps only non-locus-tag tokens."""

    def setup_method(self):
        config = {
            "fields": {
                "gene_name": {
                    "type": "single",
                    "candidates": [
                        {"source": "gene_mapping", "field": "gene_names",
                         "transform": "first_token_space"},
                    ],
                },
                "gene_name_synonyms": {
                    "type": "union",
                    "filter_not": GENE_NAME_FILTER_NOT,
                    "sources": [
                        {"source": "gene_mapping", "field": "gene_names", "delimiter": " "},
                    ],
                },
            }
        }
        self.builder = AnnotationBuilder(config)

    def test_gene_name_synonym_kept(self):
        # repA is a gene name (not a locus tag)
        gm = {"gene_names": "dnaN repA PMM0001"}
        merged = self.builder.build_merged(gm, {}, {})
        synonyms = merged.get("gene_name_synonyms", [])
        assert "repA" in synonyms

    def test_locus_tag_excluded(self):
        gm = {"gene_names": "dnaN repA PMM0001"}
        merged = self.builder.build_merged(gm, {}, {})
        synonyms = merged.get("gene_name_synonyms", [])
        assert "PMM0001" not in synonyms

    def test_canonical_gene_name_excluded_from_synonyms(self):
        # dnaN = canonical gene_name → also removed from gene_name_synonyms
        gm = {"gene_names": "dnaN repA"}
        merged = self.builder.build_merged(gm, {}, {})
        synonyms = merged.get("gene_name_synonyms", [])
        assert "dnaN" not in synonyms
        assert "repA" in synonyms

    def test_all_locus_tags_no_gene_name_synonyms(self):
        gm = {"gene_names": "PMM0001 PMM0002"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "gene_name_synonyms" not in merged


# ─── extract_first_match_in_sources ──────────────────────────────────────────

_SOURCES_EG = [{"source": "eggnog", "field": "eggNOG_OGs", "delimiter": ","}]
_SOURCES_EG_THEN_GM = [
    {"source": "eggnog",       "field": "eggNOG_OGs", "delimiter": ","},
    {"source": "gene_mapping", "field": "eggNOG",     "delimiter": ","},
]

class TestExtractFirstMatchInSources:

    def test_returns_full_match_group0(self):
        eg = {"eggNOG_OGs": "COG0592@1|root,COG0592@2|Bacteria,4648R@72275|Alteromonadaceae"}
        result = extract_first_match_in_sources(
            _SOURCES_EG, {}, eg, {},
            pattern=r"^.+@72275\|Alteromonadaceae$",
            extract_group=0,
        )
        assert result == "4648R@72275|Alteromonadaceae"

    def test_extracts_cog_prefix_group1(self):
        eg = {"eggNOG_OGs": "COG0592@1|root,COG0592@2|Bacteria,4648R@72275|Alteromonadaceae"}
        result = extract_first_match_in_sources(
            _SOURCES_EG, {}, eg, {},
            pattern=r"^(COG\d+)(?:@2\|Bacteria)?$",
            extract_group=1,
        )
        # COG0592@1|root doesn't match (wrong suffix); COG0592@2|Bacteria matches → group1 = COG0592
        assert result == "COG0592"

    def test_plain_cog_from_gene_mapping_fallback(self):
        # eggnog has no matching entry; gene_mapping plain COG#### should match
        gm = {"eggNOG": "COG1234"}
        result = extract_first_match_in_sources(
            _SOURCES_EG_THEN_GM, gm, {}, {},
            pattern=r"^(COG\d+)(?:@2\|Bacteria)?$",
            extract_group=1,
        )
        assert result == "COG1234"

    def test_eggnog_source_wins_over_gene_mapping(self):
        eg = {"eggNOG_OGs": "COG9999@2|Bacteria"}
        gm = {"eggNOG": "COG1234"}
        result = extract_first_match_in_sources(
            _SOURCES_EG_THEN_GM, gm, eg, {},
            pattern=r"^(COG\d+)(?:@2\|Bacteria)?$",
            extract_group=1,
        )
        assert result == "COG9999"

    def test_no_match_returns_none(self):
        eg = {"eggNOG_OGs": "bactNOG00989,cyaNOG01040"}
        result = extract_first_match_in_sources(
            _SOURCES_EG, {}, eg, {},
            pattern=r"^COG\d+$",
        )
        assert result is None

    def test_missing_source_field_returns_none(self):
        result = extract_first_match_in_sources(
            _SOURCES_EG, {}, {}, {},
            pattern=r"^.+@72275\|Alteromonadaceae$",
        )
        assert result is None

    def test_pattern_requires_full_match(self):
        # "COG0592@1|root" should NOT match "^COG\d+$" (not a full match)
        eg = {"eggNOG_OGs": "COG0592@1|root"}
        result = extract_first_match_in_sources(
            _SOURCES_EG, {}, eg, {},
            pattern=r"^COG\d+$",
        )
        assert result is None


# ─── AnnotationBuilder: extract_first_match type ────────────────────────────

_EXTRACT_CONFIG = {
    "fields": {
        "alteromonadaceae_og": {
            "type": "extract_first_match",
            "pattern": r"^.+@72275\|Alteromonadaceae$",
            "extract_group": 0,
            "sources": [{"source": "eggnog", "field": "eggNOG_OGs", "delimiter": ","}],
        },
        "bacteria_cog_og": {
            "type": "extract_first_match",
            "pattern": r"^(COG\d+)(?:@2\|Bacteria)?$",
            "extract_group": 1,
            "sources": [
                {"source": "eggnog",       "field": "eggNOG_OGs", "delimiter": ","},
                {"source": "gene_mapping", "field": "eggNOG",     "delimiter": ","},
            ],
        },
    }
}

class TestAnnotationBuilderExtractFirstMatch:
    def setup_method(self):
        self.builder = AnnotationBuilder(_EXTRACT_CONFIG)

    def test_alteromonadaceae_og_extracted(self):
        eg = {"eggNOG_OGs": "COG0592@1|root,COG0592@2|Bacteria,4648R@72275|Alteromonadaceae"}
        merged = self.builder.build_merged({}, eg, {})
        assert merged["alteromonadaceae_og"] == "4648R@72275|Alteromonadaceae"

    def test_alteromonadaceae_og_absent_for_pro_syn(self):
        # Pro/Syn has no @72275 entry
        eg = {"eggNOG_OGs": "COG0592@1|root,COG0592@2|Bacteria"}
        merged = self.builder.build_merged({}, eg, {})
        assert "alteromonadaceae_og" not in merged

    def test_bacteria_cog_og_from_at_format(self):
        eg = {"eggNOG_OGs": "COG0592@1|root,COG0592@2|Bacteria"}
        merged = self.builder.build_merged({}, eg, {})
        assert merged["bacteria_cog_og"] == "COG0592"

    def test_bacteria_cog_og_from_plain_cog(self):
        # plain COG#### from gene_mapping (Pro/Syn style); no eggnog source
        gm = {"eggNOG": "COG5678"}
        merged = self.builder.build_merged(gm, {}, {})
        assert merged["bacteria_cog_og"] == "COG5678"

    def test_bacteria_cog_og_absent_when_no_cog(self):
        eg = {"eggNOG_OGs": "bactNOG00989,cyaNOG01040"}
        merged = self.builder.build_merged({}, eg, {})
        assert "bacteria_cog_og" not in merged

    def test_unknown_type_skipped(self):
        config = {"fields": {"x": {"type": "unknown_type", "source": "eggnog", "field": "f"}}}
        builder = AnnotationBuilder(config)
        merged = builder.build_merged({}, {"f": "val"}, {})
        assert "x" not in merged


# ─── normalize_ec transform ──────────────────────────────────────────────────

# Fake transfer map for deterministic tests (no dependency on ec_data.json)
_FAKE_EC_MAP = {
    "1.6.99.5": ["1.6.5.11"],               # single successor
    "3.6.3.14": ["7.1.2.2"],                # single successor (class 7 reclassification)
    "1.1.1.5": ["1.1.1.303", "1.1.1.304"],  # multi-successor
    "1.13.11.42": [],                         # deleted
    "3.6.1.3": [],                            # deleted
}


@pytest.fixture(autouse=False)
def _fake_ec_map(monkeypatch):
    """Inject a deterministic EC transfer map for tests."""
    monkeypatch.setattr(_at, "_EC_TRANSFER_MAP", _FAKE_EC_MAP)
    yield
    # Reset so lazy-load works again for other tests
    monkeypatch.setattr(_at, "_EC_TRANSFER_MAP", None)


class TestNormalizeEcTransform:
    """Unit tests for the _tx_normalize_ec transform function."""

    def test_current_ec_passes_through(self, _fake_ec_map):
        assert _tx_normalize_ec("2.7.7.7") == "2.7.7.7"

    def test_single_successor(self, _fake_ec_map):
        assert _tx_normalize_ec("1.6.99.5") == "1.6.5.11"

    def test_multi_successor_returns_list(self, _fake_ec_map):
        result = _tx_normalize_ec("1.1.1.5")
        assert result == ["1.1.1.303", "1.1.1.304"]

    def test_deleted_returns_empty(self, _fake_ec_map):
        assert _tx_normalize_ec("1.13.11.42") == ""

    def test_partial_ec_passes_through(self, _fake_ec_map):
        # Partial ECs (with dashes) should pass through unchanged
        assert _tx_normalize_ec("1.1.1.-") == "1.1.1.-"

    def test_two_level_partial_passes_through(self, _fake_ec_map):
        assert _tx_normalize_ec("1.1.-.-") == "1.1.-.-"

    def test_one_level_partial_passes_through(self, _fake_ec_map):
        assert _tx_normalize_ec("1.-.-.-") == "1.-.-.-"

    def test_empty_returns_empty(self, _fake_ec_map):
        assert _tx_normalize_ec("") == ""

    def test_dash_returns_empty(self, _fake_ec_map):
        assert _tx_normalize_ec("-") == ""

    def test_whitespace_stripped(self, _fake_ec_map):
        assert _tx_normalize_ec("  1.6.99.5  ") == "1.6.5.11"


class TestNormalizeEcInBuilder:
    """Integration tests: normalize_ec transform within AnnotationBuilder union resolution."""

    EC_CONFIG = {
        "fields": {
            "ec_numbers": {
                "type": "union",
                "filter": r"^\d+\.[\d\-]+\.[\d\-]+[\.\-]",
                "sources": [
                    {"source": "gene_mapping", "field": "ec_numbers", "delimiter": ",",
                     "transform": "normalize_ec"},
                    {"source": "gene_mapping", "field": "kegg", "delimiter": ",",
                     "transform": "normalize_ec"},
                    {"source": "eggnog", "field": "EC", "delimiter": ",",
                     "transform": "normalize_ec"},
                ],
            },
        },
    }

    def setup_method(self):
        self.builder = AnnotationBuilder(self.EC_CONFIG)

    def test_transferred_ec_replaced(self, _fake_ec_map):
        gm = {"ec_numbers": "1.6.99.5"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.6.5.11" in merged["ec_numbers"]
        assert "1.6.99.5" not in merged["ec_numbers"]

    def test_deleted_ec_dropped(self, _fake_ec_map):
        gm = {"ec_numbers": "1.13.11.42"}
        merged = self.builder.build_merged(gm, {}, {})
        assert merged.get("ec_numbers") is None  # only EC was deleted → no field

    def test_multi_successor_expanded(self, _fake_ec_map):
        gm = {"ec_numbers": "1.1.1.5"}
        merged = self.builder.build_merged(gm, {}, {})
        ec = merged["ec_numbers"]
        assert "1.1.1.303" in ec
        assert "1.1.1.304" in ec
        assert "1.1.1.5" not in ec

    def test_partial_ec_preserved(self, _fake_ec_map):
        gm = {"ec_numbers": "1.1.1.-"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.1.1.-" in merged["ec_numbers"]

    def test_two_level_partial_preserved(self, _fake_ec_map):
        gm = {"ec_numbers": "1.1.-.-"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.1.-.-" in merged["ec_numbers"]

    def test_one_level_partial_preserved(self, _fake_ec_map):
        gm = {"ec_numbers": "1.-.-.-"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.-.-.-" in merged["ec_numbers"]

    def test_kegg_ko_rejected_by_filter(self, _fake_ec_map):
        eg = {"EC": "K02169,2.7.7.7"}
        merged = self.builder.build_merged({}, eg, {})
        ec = merged["ec_numbers"]
        assert "2.7.7.7" in ec
        assert "K02169" not in ec

    def test_mixed_valid_transferred_deleted_partial(self, _fake_ec_map):
        gm = {"ec_numbers": "2.7.7.7,1.6.99.5,1.13.11.42,1.1.-.-"}
        merged = self.builder.build_merged(gm, {}, {})
        ec = merged["ec_numbers"]
        assert "2.7.7.7" in ec       # valid: kept
        assert "1.6.5.11" in ec       # transferred: replaced
        assert "1.13.11.42" not in ec  # deleted: dropped
        assert "1.6.99.5" not in ec    # old: gone
        assert "1.1.-.-" in ec         # partial: kept

    def test_dedup_across_sources_after_normalization(self, _fake_ec_map):
        # gene_mapping.ec_numbers and kegg both have the same obsolete EC
        gm = {"ec_numbers": "1.6.99.5", "kegg": "1.6.99.5"}
        merged = self.builder.build_merged(gm, {}, {})
        # Successor should appear exactly once
        assert merged["ec_numbers"].count("1.6.5.11") == 1


class TestResolveEcChain:
    """Unit tests for _resolve_ec_chain (transfer chain resolution)."""

    def test_current_ec_returns_itself(self):
        tmap = {"1.1.1.1": ["1.1.1.2"]}
        assert _resolve_ec_chain("2.7.7.7", tmap) == ["2.7.7.7"]

    def test_single_hop(self):
        tmap = {"1.6.99.5": ["1.6.5.11"]}
        assert _resolve_ec_chain("1.6.99.5", tmap) == ["1.6.5.11"]

    def test_chain_depth_2(self):
        # 1.6.99.5 → 1.6.5.11 → 1.6.5.9
        tmap = {"1.6.99.5": ["1.6.5.11"], "1.6.5.11": ["1.6.5.9"]}
        assert _resolve_ec_chain("1.6.99.5", tmap) == ["1.6.5.9"]

    def test_chain_ends_at_deleted(self):
        # 1.6.2.1 → 1.6.99.3 → deleted
        tmap = {"1.6.2.1": ["1.6.99.3"], "1.6.99.3": []}
        assert _resolve_ec_chain("1.6.2.1", tmap) == []

    def test_deleted_entry_returns_empty(self):
        tmap = {"1.13.11.42": []}
        assert _resolve_ec_chain("1.13.11.42", tmap) == []

    def test_multi_successor_no_chain(self):
        tmap = {"1.1.1.5": ["1.1.1.303", "1.1.1.304"]}
        assert _resolve_ec_chain("1.1.1.5", tmap) == ["1.1.1.303", "1.1.1.304"]

    def test_multi_successor_one_branch_chains(self):
        # A → [B, C]; B → D (chain); C is current
        tmap = {"A": ["B", "C"], "B": ["D"]}
        assert _resolve_ec_chain("A", tmap) == ["D", "C"]

    def test_cycle_guard(self):
        # Pathological cycle: A → B → A
        tmap = {"A": ["B"], "B": ["A"]}
        result = _resolve_ec_chain("A", tmap)
        # Should not infinite loop; returns something reasonable
        assert isinstance(result, list)


class TestNormalizeEcChainedInBuilder:
    """Integration: chained transfers resolved end-to-end in AnnotationBuilder."""

    EC_CONFIG = {
        "fields": {
            "ec_numbers": {
                "type": "union",
                "filter": r"^\d+\.[\d\-]+\.[\d\-]+[\.\-]",
                "sources": [
                    {"source": "gene_mapping", "field": "ec_numbers", "delimiter": ",",
                     "transform": "normalize_ec"},
                ],
            },
        },
    }

    # Chained fake map: already resolved (as _build_ec_transfer_map would produce)
    _CHAINED_EC_MAP = {
        "1.6.99.5": ["1.6.5.9"],   # was 1.6.99.5 → 1.6.5.11 → 1.6.5.9
        "1.6.5.11": ["1.6.5.9"],   # intermediate is also resolved
        "1.6.2.1": [],              # chain ends at deleted
        "1.6.99.3": [],             # deleted
    }

    @pytest.fixture(autouse=True)
    def _chained_ec_map(self, monkeypatch):
        monkeypatch.setattr(_at, "_EC_TRANSFER_MAP", self._CHAINED_EC_MAP)
        yield
        monkeypatch.setattr(_at, "_EC_TRANSFER_MAP", None)

    def setup_method(self):
        self.builder = AnnotationBuilder(self.EC_CONFIG)

    def test_depth2_chain_resolves_to_final(self):
        gm = {"ec_numbers": "1.6.99.5"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.6.5.9" in merged["ec_numbers"]
        assert "1.6.5.11" not in merged["ec_numbers"]
        assert "1.6.99.5" not in merged["ec_numbers"]

    def test_intermediate_also_resolves(self):
        gm = {"ec_numbers": "1.6.5.11"}
        merged = self.builder.build_merged(gm, {}, {})
        assert "1.6.5.9" in merged["ec_numbers"]
        assert "1.6.5.11" not in merged["ec_numbers"]

    def test_chain_ending_at_deleted_drops(self):
        gm = {"ec_numbers": "1.6.2.1"}
        merged = self.builder.build_merged(gm, {}, {})
        assert merged.get("ec_numbers") is None
