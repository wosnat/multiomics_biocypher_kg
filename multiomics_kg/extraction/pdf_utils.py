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


def upload_pdf(client, pdf_path: Path) -> str:
    """Upload a PDF via the OpenAI Files API and return the file_id.

    `client` is an OpenAI client passed in by the caller (this module keeps no
    openai dependency of its own). Lives here so every extraction module shares
    one uploader instead of redefining it.
    """
    with open(pdf_path, "rb") as fh:
        f = client.files.create(file=fh, purpose="user_data")
    return f.id


# ── Page-range chunking (for fitting large PDFs under per-request token caps) ──


def count_pages(pdf_path: Path) -> int:
    """Number of pages in a PDF (0 if unreadable)."""
    if PdfReader is None:
        return 0
    try:
        return len(PdfReader(pdf_path).pages)
    except Exception as e:
        logger.error("Failed to count pages in %s: %s", pdf_path, e)
        return 0


def find_references_page(pdf_path: Path) -> Optional[int]:
    """Return the 0-indexed page where the References/Bibliography section starts.

    Used to skip reference pages (fewer tokens, less name-drop noise). Returns
    None if no such heading is found.
    """
    if PdfReader is None:
        return None
    patterns = [r"\nReferences\s*\n", r"\nREFERENCES\s*\n",
                r"\nBibliography\s*\n", r"\nLiterature Cited\s*\n"]
    try:
        reader = PdfReader(pdf_path)
        for i, page in enumerate(reader.pages):
            text = page.extract_text() or ""
            if any(re.search(p, text) for p in patterns):
                return i
        return None
    except Exception as e:
        logger.error("Failed scanning for references in %s: %s", pdf_path, e)
        return None


def page_chunks(last_page: int, chunk_pages: int) -> list[tuple[int, int]]:
    """Inclusive (start, end) page ranges covering 0..last_page in chunk_pages steps."""
    chunks: list[tuple[int, int]] = []
    start = 0
    while start <= last_page:
        end = min(start + chunk_pages - 1, last_page)
        chunks.append((start, end))
        start = end + 1
    return chunks


def write_page_range_pdf(pdf_path: Path, start: int, end: int, dest: Path) -> Optional[Path]:
    """Write pages [start, end] (inclusive, 0-indexed) of pdf_path to dest. Returns dest."""
    if PdfReader is None or PdfWriter is None:
        return None
    try:
        reader = PdfReader(pdf_path)
        total = len(reader.pages)
        start = max(0, start)
        end = min(end, total - 1)
        writer = PdfWriter()
        for i in range(start, end + 1):
            writer.add_page(reader.pages[i])
        with open(dest, "wb") as f:
            writer.write(f)
        return Path(dest)
    except Exception as e:
        logger.error("Failed writing page range %d-%d of %s: %s", start, end, pdf_path, e)
        return None
