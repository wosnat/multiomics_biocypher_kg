"""Tests for multiomics_kg/extraction/pdf_utils.py.

Regression coverage for the reference-page scanner, which silently degrades:
`find_references_page` catches every exception and returns None, so a parser
failure is indistinguishable from "this paper has no References heading". The
practical cost is that the whole PDF (bibliography included) gets shipped to
the extraction model instead of just the narrative pages.
"""
from pathlib import Path

import pytest

from multiomics_kg.extraction.pdf_utils import find_references_page

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
