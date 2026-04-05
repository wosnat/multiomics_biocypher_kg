# multiomics_kg/extraction/structured_parsers.py
"""Shared file reading utilities for XLS, HTML, TXT, DOC, DOCX."""
import logging
import subprocess
import tempfile
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class _HTMLTextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self._parts: list[str] = []

    def handle_data(self, data):
        self._parts.append(data)

    def get_text(self) -> str:
        return " ".join(self._parts)


def read_html(path: Path) -> str:
    try:
        raw = path.read_text(encoding="utf-8", errors="replace")
        parser = _HTMLTextExtractor()
        parser.feed(raw)
        return parser.get_text()
    except Exception as e:
        logger.error("Failed to read HTML %s: %s", path, e)
        return ""


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        logger.error("Failed to read text %s: %s", path, e)
        return ""


def read_xls(path: Path) -> str:
    """Read all sheets of an XLS/XLSX file as text."""
    try:
        xls = pd.ExcelFile(path)
        parts = []
        for sheet in xls.sheet_names:
            df = pd.read_excel(xls, sheet_name=sheet)
            parts.append(f"--- Sheet: {sheet} ---")
            parts.append(df.to_string(max_rows=200))
        return "\n".join(parts)
    except Exception as e:
        logger.error("Failed to read XLS %s: %s", path, e)
        return ""


def read_xls_as_df(path: Path, sheet: int = 0,
                    header: object = None) -> pd.DataFrame:
    """Read an XLS/XLSX sheet as a DataFrame."""
    return pd.read_excel(path, sheet_name=sheet, header=header)


def read_doc(path: Path) -> str:
    """Read old-format .doc file using antiword."""
    try:
        result = subprocess.run(
            ["antiword", str(path)],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            return result.stdout
        logger.error("antiword failed on %s: %s", path, result.stderr[:200])
        return ""
    except FileNotFoundError:
        logger.warning("antiword not installed — cannot read %s", path)
        return ""
    except Exception as e:
        logger.error("Failed to read DOC %s: %s", path, e)
        return ""


def read_docx(path: Path) -> str:
    """Read .docx file using python-docx."""
    try:
        import docx
        doc = docx.Document(str(path))
        return "\n".join(p.text for p in doc.paragraphs if p.text.strip())
    except Exception as e:
        logger.error("Failed to read DOCX %s: %s", path, e)
        return ""


def read_file(path: Path) -> str:
    """Read any supported file type, dispatching by extension."""
    suffix = path.suffix.lower()
    readers = {
        ".html": read_html, ".htm": read_html,
        ".txt": read_text,
        ".xls": read_xls, ".xlsx": read_xls,
        ".doc": read_doc,
        ".docx": read_docx,
    }
    reader = readers.get(suffix)
    if reader is None:
        return ""
    return reader(path)


def scan_supplementary_text(paper_dir: Path,
                            max_file_size: int = 500_000,
                            ) -> dict[str, str]:
    """Scan all readable non-PDF files in the paper directory.

    Returns {relative_path: text_content}.
    Skips generated files, huge files, and raw data directories.
    """
    skip_patterns = ["_resolved.", "_report.", "paperconfig", "cluster_extraction"]
    skip_dirs = {"Goldenspike", "Expression"}
    skip_suffixes = {".pdf", ".eps", ".csv", ".zip", ".gz", ".fasta", ".faa",
                     ".gff", ".gbk", ".gbff", ".png", ".jpg", ".gif"}

    sources: dict[str, str] = {}
    for path in sorted(paper_dir.rglob("*")):
        if not path.is_file():
            continue
        if any(p in path.name for p in skip_patterns):
            continue
        if any(d in path.parts for d in skip_dirs):
            continue
        suffix = path.suffix.lower()
        if suffix in skip_suffixes:
            continue
        if path.stat().st_size > max_file_size:
            continue
        text = read_file(path)
        if text.strip():
            sources[str(path.relative_to(paper_dir))] = text
    return sources
