"""A DE table lists each gene once: enforce_one_row_per_gene.

Two rows that resolve to one gene share an expression-edge id, and BioCypher keeps
one of them. When one row names the gene more directly than the other, the weaker
row lost a tie it should never have entered: it goes to its one remaining
candidate, or stays unresolved. Fixtures mirror the cases found 2026-11-25
(backlog "Two table rows -> one gene").
"""
from multiomics_kg.utils.gene_id_utils import (
    MappingData,
    enforce_one_row_per_gene,
    resolve_row,
)


def _run(rows, name_col, id_columns, md):
    results = [resolve_row(r, name_col, id_columns, md) for r in rows]
    return enforce_one_row_per_gene(rows, results, name_col, id_columns, md)


def test_symbol_onto_gene_listed_by_locus_tag_moves_to_other_candidate():
    # Al-Hosani 2015: the table lists A9601_16521 (atpC, gamma) by locus tag and
    # also has a row "atpC", the epsilon subunit A9601_16411 (Cyanorak synonym atpC).
    md = MappingData(
        multi_lookup={"atpC": ["A9601_16411", "A9601_16521"]},
        named_lookup={"atpC": ["A9601_16521"]},
        locus_tags={"A9601_16411", "A9601_16521"},
    )
    rows = [{"Gene id": "A9601_16521"}, {"Gene id": "atpC"}]
    out = _run(rows, "Gene id", [], md)
    assert out[0] == ("A9601_16521", "locus_tag:Gene id")
    assert out[1][0] == "A9601_16411"
    assert out[1][1].startswith("reassigned:")


def test_weaker_symbol_gives_way_to_singleton_symbol():
    # Al-Hosani 2015: Cyanorak names A9601_09461 (signal peptidase II) "ispA" with
    # synonym lspA; A9601_11751 (crtE) carries synonym ispA. Rows "ispA" and "lspA".
    md = MappingData(
        multi_lookup={"ispA": ["A9601_09461", "A9601_11751"], "lspA": ["A9601_09461"]},
        named_lookup={"ispA": ["A9601_09461"]},
        locus_tags={"A9601_09461", "A9601_11751"},
    )
    rows = [{"Gene id": "ispA"}, {"Gene id": "lspA"}]
    out = _run(rows, "Gene id", [], md)
    assert out[1] == ("A9601_09461", "multi:Gene id")
    assert out[0][0] == "A9601_11751"


def test_symbol_loses_to_accession_and_has_no_other_candidate():
    # Kratzl 2024: proA.1/Q31KZ6 resolves by accession; row proA/Q31KX4 (accession
    # unknown) falls back to the symbol and lands on the same gene.
    md = MappingData(
        multi_lookup={"Q31KZ6": ["H6G84_01725"], "proA": ["H6G84_01725"],
                      "proA.1": ["H6G84_01725"]},
        tier2_tokens={"Q31KZ6"},
        locus_tags={"H6G84_01725"},
    )
    cols = [{"column": "Protein ID", "id_type": "uniprot_accession"},
            {"column": "Gene Name", "id_type": "gene_name"}]
    rows = [{"Gene Name": "proA.1", "Protein ID": "Q31KZ6"},
            {"Gene Name": "proA", "Protein ID": "Q31KX4"}]
    out = _run(rows, "Gene Name", cols, md)
    assert out[0][0] == "H6G84_01725"
    assert out[1][0] is None
    assert out[1][1].startswith("duplicate_claim:")


def test_old_alias_loses_to_row_naming_the_locus_tag():
    # He 2022: gene-PMM0236 (the KG gene itself) and gene-PMM1806 (an old tag of it).
    md = MappingData(
        specific_lookup={"gene-PMM1806": "PMM0236", "PMM1806": "PMM0236"},
        locus_tags={"PMM0236"},
    )
    rows = [{"GeneID": "gene-PMM0236"}, {"GeneID": "gene-PMM1806"}]
    out = _run(rows, "GeneID", [], md)
    assert out[0][0] == "PMM0236"
    assert out[1] == (None, "duplicate_claim:tier1:GeneID")


def test_equally_direct_rows_are_left_alone():
    # Coe 2024: NCBI and RAST models of one gene, both rows carry the locus tag.
    md = MappingData(locus_tags={"PMN2A_0530"})
    cols = [{"column": "Gene ID", "id_type": "annotation_specific"}]
    rows = [{"NCBI ID_3": "PMN2A_0530", "Gene ID": "cds-WP_011294627.1"},
            {"NCBI ID_3": "PMN2A_0530", "Gene ID": "fig|1218.94.peg.1247"}]
    out = _run(rows, "NCBI ID_3", cols, md)
    assert out == [("PMN2A_0530", "locus_tag:NCBI ID_3")] * 2


def test_unaffected_rows_pass_through():
    md = MappingData(locus_tags={"A", "B"})
    rows = [{"id": "A"}, {"id": "B"}, {"id": "zzz"}]
    out = _run(rows, "id", [], md)
    assert out == [("A", "locus_tag:id"), ("B", "locus_tag:id"), (None, "unresolved")]
