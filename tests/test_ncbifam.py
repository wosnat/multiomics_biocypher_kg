from multiomics_kg.utils.ncbifam import parse_hmm_pgap_rows, HYPOTH_FAMILY_TYPES


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
