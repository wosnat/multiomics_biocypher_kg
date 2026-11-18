"""Tests for multiomics_kg/extraction/pdf_utils.py.

Regression coverage for the reference-page scanner, which silently degrades:
`find_references_page` catches every exception and returns None, so a parser
failure is indistinguishable from "this paper has no References heading". The
practical cost is that the whole PDF (bibliography included) gets shipped to
the extraction model instead of just the narrative pages.
"""
from pathlib import Path

import pytest

from multiomics_kg.extraction.pdf_utils import decode_pdf_glyphs, find_references_page

PROJECT_ROOT = Path(__file__).parent.parent

# Lu 2026 (AEM). Its embedded fonts declare per-character /Widths as indirect
# references, which pypdf <= 6.7.2 stored unresolved and then tried to sum:
#   TypeError: unsupported operand type(s) for +: 'float' and 'IndirectObject'
# (pypdf/_font.py Font.text_width, via the dict-encoding branch of
# _collect_tt_t1_character_widths, which stores widths without int()/get_object()).
# Fixed upstream; guarded here so a pypdf downgrade can't silently reintroduce it.
LU_2026_PDF = (
    PROJECT_ROOT
    / "data/Prochlorococcus/papers_and_supp/Lu 2026"
    / "lu-et-al-2026-vesicle-associated-exudates-from-alteromonas-enhance-growth-"
      "and-survival-of-prochlorococcus-in-batch.pdf"
)


@pytest.mark.skipif(not LU_2026_PDF.exists(), reason="Lu 2026 PDF not present")
def test_find_references_page_handles_indirect_font_widths():
    """The scanner must locate the References heading, not swallow a parser error."""
    page = find_references_page(LU_2026_PDF)

    assert page is not None, (
        "find_references_page returned None for a PDF that does have a "
        "References section -- the scanner most likely swallowed a pypdf "
        "exception. Check the logged 'Failed scanning for references' error."
    )
    # Bibliography sits in the last third of this 21-page article.
    assert 10 < page < 21, f"References page index {page} is implausible"


@pytest.mark.skipif(not LU_2026_PDF.exists(), reason="Lu 2026 PDF not present")
def test_find_references_page_logs_nothing_on_success(caplog):
    """A successful scan must not emit the failure log line."""
    with caplog.at_level("ERROR"):
        find_references_page(LU_2026_PDF)

    assert not [r for r in caplog.records if "Failed scanning for references" in r.message], (
        "pdf_utils logged a scan failure; the PDF parser raised and the error "
        "was swallowed into a None return."
    )


@pytest.mark.parametrize("raw, decoded", [
    # Steglich 2006 title: the glyph before the dagger prints as nothing.
    ("Light Sensing in Prochlorococcus/H17188†", "Light Sensing in Prochlorococcus†"),
    # Dominguez 2017 abstract.
    ("P value of /H110210.05", "P value of <0.05"),
    # A 5-digit code running into a number: "m-2 s-1".
    ("m /H110022/H18528s/H110021", "m −2·s−1"),
    # A 4-digit code running into a number.
    ("0.1/H9262m", "0.1μm"),
    ("/H92621", "μ1"),
    # PNAS header separator never becomes the KG's "|" array delimiter.
    ("12557 /H20841PNAS", "12557  PNAS"),
])
def test_decode_pdf_glyphs(raw, decoded):
    assert decode_pdf_glyphs(raw) == decoded


def test_decode_pdf_glyphs_drops_unknown_codes():
    assert decode_pdf_glyphs("sort /H20888x and /H99991") == "sort x and 1"
