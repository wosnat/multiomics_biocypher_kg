# multiomics_kg/extraction/pdf_utils.py
"""Shared PDF utilities: page-to-base64 conversion and text extraction."""
import base64
import io
import logging
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

try:
    from pypdf import PdfReader, PdfWriter
except Exception:
    PdfReader = None
    PdfWriter = None


def pdf_pages_to_base64(pdf_path: Path,
                        page_range: Optional[tuple[int, int]] = None,
                        ) -> list[str]:
    """Convert PDF pages to base64-encoded single-page PDFs.

    Args:
        pdf_path: Path to the PDF file.
        page_range: (start, end) 0-indexed inclusive. None = all pages.

    Returns:
        List of base64 strings, one per page.
    """
    if PdfReader is None:
        logger.warning("pypdf not installed — cannot read %s", pdf_path)
        return []
    try:
        reader = PdfReader(pdf_path)
        total = len(reader.pages)
        start = page_range[0] if page_range else 0
        end = page_range[1] if page_range else total - 1
        start = max(0, min(start, total - 1))
        end = max(start, min(end, total - 1))

        result = []
        for i in range(start, end + 1):
            writer = PdfWriter()
            writer.add_page(reader.pages[i])
            buf = io.BytesIO()
            writer.write(buf)
            result.append(base64.standard_b64encode(buf.getvalue()).decode())
        return result
    except Exception as e:
        logger.error("Failed to read PDF %s: %s", pdf_path, e)
        return []


def extract_pdf_text(pdf_path: Path, max_pages: int = 30) -> str:
    """Extract text from PDF, stopping at References section."""
    if PdfReader is None:
        return ""
    try:
        text = ""
        with open(pdf_path, "rb") as f:
            reader = PdfReader(f)
            for i in range(min(len(reader.pages), max_pages)):
                page_text = reader.pages[i].extract_text() or ""
                for pat in [r'\nReferences\s*\n', r'\nREFERENCES\s*\n',
                            r'\nBibliography\s*\n']:
                    m = re.search(pat, page_text)
                    if m:
                        return text + page_text[:m.start()]
                text += page_text
        return text
    except Exception as e:
        logger.error("Failed to extract text from %s: %s", pdf_path, e)
        return ""


def collect_pdf_files(paper_dir: Path,
                      main_pdf_path: Optional[Path] = None,
                      max_file_size_mb: int = 10,
                      ) -> list[Path]:
    """Collect all PDF files in a paper directory.

    Returns main PDF first (if provided), then supplementary PDFs.
    Skips files larger than max_file_size_mb.
    """
    pdfs = []
    skip = {"_resolved", "_report", "cluster_extraction"}
    max_bytes = max_file_size_mb * 1024 * 1024

    if main_pdf_path and main_pdf_path.exists():
        pdfs.append(main_pdf_path)

    for p in sorted(paper_dir.rglob("*.pdf")):
        if p == main_pdf_path:
            continue
        if any(s in p.name for s in skip):
            continue
        if p.stat().st_size > max_bytes:
            continue
        pdfs.append(p)

    return pdfs
