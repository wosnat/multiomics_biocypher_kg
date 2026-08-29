import pytest

from multiomics_kg.utils.ncbifam import parse_hmm_pgap_rows, HYPOTH_FAMILY_TYPES, parse_tigr_role_names, parse_tigr_role_link


def _row(**kw):
    base = {"#ncbi_accession": "NF000001.1", "family_type": "equivalog",
            "product_name": "widget synthase", "gene_symbol": "", "comment": "",
            "gene_synonyms": "", "ec_numbers": "", "go_terms": "", "pmids": ""}
    base.update(kw)
    return base


def test_parse_basic_and_version_strip():
    ref = parse_hmm_pgap_rows([_row()])
    assert ref["NF000001"]["name"] == "widget synthase"
    assert ref["NF000001"]["family_type"] == "equivalog"
    assert "gene_symbol" not in ref["NF000001"]        # sparse
    assert "description" not in ref["NF000001"]


def test_parse_rich_row():
    ref = parse_hmm_pgap_rows([_row(**{"#ncbi_accession": "TIGR00198.1",
        "gene_symbol": "katG", "comment": "catalase-peroxidase HPI",
        "ec_numbers": "1.11.1.21", "pmids": "9006042,9871101"})])
    e = ref["TIGR00198"]
    assert e["gene_symbol"] == "katG"
    assert e["description"] == "catalase-peroxidase HPI"
    assert e["ec_numbers"] == ["1.11.1.21"] and e["pmids"] == ["9006042", "9871101"]


def test_hypoth_types_constant():
    assert "hypoth_equivalog" in HYPOTH_FAMILY_TYPES
    assert "hypoth_equivalog_domain" in HYPOTH_FAMILY_TYPES


def test_gene_synonyms_comma_delimited_real_cell():
    # Real observed cell value from hmm_PGAP.tsv (TIGR-style row), pinning the
    # comma delimiter empirically verified against the downloaded TSV.
    ref = parse_hmm_pgap_rows([_row(gene_synonyms="cybA,dhsC")])
    assert ref["NF000001"]["gene_synonyms"] == ["cybA", "dhsC"]


def test_go_terms_stored_raw():
    ref = parse_hmm_pgap_rows([_row(go_terms="GO:0008493,GO:0015904")])
    assert ref["NF000001"]["go_terms"] == ["GO:0008493", "GO:0015904"]


def test_no_version_suffix_accession_unaffected():
    ref = parse_hmm_pgap_rows([_row(**{"#ncbi_accession": "NF000002"})])
    assert "NF000002" in ref


# ── TIGRFAMs 15.0 role archive tests ──────────────────────────────────────


ROLE_NAMES = """\
role_id:\t100\tmainrole:\tCentral intermediary metabolism
role_id:\t100\tsub1role:\tAmino sugars
role_id:\t132\tmainrole:\tDNA metabolism
role_id:\t132\tsub1role:\tDNA replication, recombination, and repair
role_id:\t719\tsub1role:\tOrphan without mainrole
""".splitlines()

ROLE_LINK = """\
TIGR00001\t132
TIGR00002\t100
TIGR00003\t719
TIGR00004\t999
""".splitlines()


def test_parse_tigr_role_names_pairs_main_and_sub():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert roles["100"] == {"mainrole": "Central intermediary metabolism", "sub1role": "Amino sugars"}
    assert roles["132"]["sub1role"] == "DNA replication, recombination, and repair"


def test_parse_tigr_role_names_drops_roles_without_mainrole():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert "719" not in roles


def test_parse_tigr_role_link_keeps_only_named_roles():
    roles = parse_tigr_role_names(ROLE_NAMES)
    link = parse_tigr_role_link(ROLE_LINK, roles)
    assert link == {"TIGR00001": "132", "TIGR00002": "100"}


def test_parse_tigr_role_link_strips_version_suffix():
    roles = parse_tigr_role_names(ROLE_NAMES)
    assert parse_tigr_role_link(["TIGR00005.1\t132"], roles) == {"TIGR00005": "132"}


def test_parsers_fail_loud_on_nonempty_garbage():
    with pytest.raises(ValueError):
        parse_tigr_role_names(["this is not a role line", "neither is this"])
    with pytest.raises(ValueError):
        parse_tigr_role_link(["garbage"], {"132": {"mainrole": "x", "sub1role": "y"}})


def test_parsers_accept_empty_input():
    assert parse_tigr_role_names([]) == {}
    assert parse_tigr_role_link([], {}) == {}
