"""Tests for PDF page-chunking and per-chunk extraction merge."""
from multiomics_kg.extraction.pdf_utils import page_chunks
from multiomics_kg.extraction.topics.extract import merge_extractions


# ── page_chunks ──


def test_page_chunks_exact_multiple():
    assert page_chunks(29, 15) == [(0, 14), (15, 29)]


def test_page_chunks_remainder():
    assert page_chunks(30, 15) == [(0, 14), (15, 29), (30, 30)]


def test_page_chunks_single_chunk():
    assert page_chunks(10, 15) == [(0, 10)]


def test_page_chunks_single_page():
    assert page_chunks(0, 15) == [(0, 0)]


# ── merge_extractions ──


def _chunk(genes=None, pathways=None, conf="high", unc="none", note=""):
    return {
        "genes": genes or [],
        "pathways": pathways or [],
        "self_assessment": conf,
        "uncaptured_identifiers": unc,
        "uncaptured_identifiers_note": note,
    }


def test_merge_dedups_gene_and_unions_identifiers():
    c1 = _chunk(genes=[{"gene_name": "ntcA", "surface_form": "NtcA", "identifiers": ["PMM0246"],
                        "strain": "Prochlorococcus MED4", "prominence": "peripheral", "evidence": "a"}])
    c2 = _chunk(genes=[{"gene_name": "ntcA", "surface_form": "NtcA", "identifiers": ["PMM0246b"],
                        "strain": "Prochlorococcus MED4", "prominence": "central", "evidence": "b"}])
    merged = merge_extractions([c1, c2])
    assert len(merged["genes"]) == 1
    g = merged["genes"][0]
    assert set(g["identifiers"]) == {"PMM0246", "PMM0246b"}
    assert g["prominence"] == "central"   # central-wins


def test_merge_keeps_distinct_genes_and_strains():
    c1 = _chunk(genes=[{"gene_name": "ntcA", "surface_form": "", "identifiers": [],
                        "strain": "Prochlorococcus MED4", "prominence": "central", "evidence": ""}])
    c2 = _chunk(genes=[{"gene_name": "ntcA", "surface_form": "", "identifiers": [],
                        "strain": "Prochlorococcus MIT9313", "prominence": "central", "evidence": ""}])
    merged = merge_extractions([c1, c2])
    assert len(merged["genes"]) == 2   # same name, different strain → distinct


def test_merge_dedups_pathways_by_kegg_id():
    c1 = _chunk(pathways=[{"surface_form": "glycolysis", "kegg_pathway_id": "ko00010",
                           "prominence": "peripheral", "evidence": ""}])
    c2 = _chunk(pathways=[{"surface_form": "Glycolysis/Gluconeogenesis", "kegg_pathway_id": "ko00010",
                           "prominence": "central", "evidence": ""}])
    merged = merge_extractions([c1, c2])
    assert len(merged["pathways"]) == 1
    assert merged["pathways"][0]["prominence"] == "central"


def test_merge_confidence_min_and_uncaptured_max():
    merged = merge_extractions([
        _chunk(conf="high", unc="none"),
        _chunk(conf="low", unc="many", note="tables"),
    ])
    assert merged["self_assessment"] == "low"        # lowest confidence
    assert merged["uncaptured_identifiers"] == "many"  # highest gap
    assert "tables" in merged["uncaptured_identifiers_note"]
