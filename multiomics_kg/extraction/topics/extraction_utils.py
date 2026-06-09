"""Topic-extraction file I/O and paperconfig helpers.

Single source of truth for where the raw topics extraction lives and how the
per-paper strain set is derived. Used by extract.py (Stage 1), the resolution
step (Stage 2), and the adapter (Stage 3).
"""
import json
import logging
from pathlib import Path

from multiomics_kg.utils.paperconfig_utils import (
    get_experiments,
    get_publication,
    load_all_paperconfigs,
)

logger = logging.getLogger(__name__)

# All per-paper topic artifacts live in this subfolder of the paper directory.
TOPICS_SUBDIR = "publication_topics"
TOPICS_FILENAME = "topics.json"               # raw extraction (Stage 1)
RESOLVED_FILENAME = "topics_resolved.json"    # resolved (Stage 2)
REPORT_FILENAME = "resolution_report.txt"     # diagnostic report (Stage 2)

# Legacy flat-file locations (pre-subfolder); migrated on demand.
LEGACY_TOPICS_FILENAME = "publication_topics.json"
LEGACY_RESOLVED_FILENAME = "publication_topics_resolved.json"


def get_paper_strains(config: dict) -> list[str]:
    """Candidate strain set for a paper = organisms across its experiments.

    Includes each experiment's ``organism`` and any ``treatment_organism``
    (coculture partners like *Alteromonas* HOT1A3 are real genome strains in the
    KG and may be discussed). Order-preserving, de-duplicated.
    """
    strains: list[str] = []
    for exp in get_experiments(config).values():
        if not isinstance(exp, dict):
            continue
        for key in ("organism", "treatment_organism"):
            org = exp.get(key)
            if org:
                strains.append(org)
    return list(dict.fromkeys(strains))


def resolve_pdf_path(config: dict, paper_dir: Path) -> Path | None:
    """Resolve the main PDF path for a paper, or None if absent/missing."""
    pub = get_publication(config)
    pdf_str = pub.get("papermainpdf", "")
    if pdf_str:
        pdf_path = Path(pdf_str)
        if not pdf_path.is_absolute():
            pdf_path = Path.cwd() / pdf_path
        return pdf_path if pdf_path.exists() else None
    # Fall back to any PDF in the paper dir.
    pdfs = sorted(paper_dir.glob("*.pdf"))
    return pdfs[0] if pdfs else None


def find_all_papers() -> list[tuple[Path, dict, Path, list[str]]]:
    """Discover every paperconfig that has a usable main PDF.

    Returns list of (paper_dir, pub_config, pdf_path, strains).
    Papers without a resolvable PDF are skipped (no narrative to extract from).
    """
    papers: list[tuple[Path, dict, Path, list[str]]] = []
    for pc_path, pc in load_all_paperconfigs():
        paper_dir = pc_path.parent
        pub = get_publication(pc)
        if not pub or pub.get("skip_pdf_extraction"):
            continue
        pdf_path = resolve_pdf_path(pc, paper_dir)
        if pdf_path is None:
            logger.info("No PDF for %s — skipping topic extraction", paper_dir.name)
            continue
        papers.append((paper_dir, pub, pdf_path, get_paper_strains(pc)))
    return papers


# ── Artifact paths (subfolder) ──


def topics_dir(paper_dir: Path) -> Path:
    return Path(paper_dir) / TOPICS_SUBDIR


def topics_path(paper_dir: Path) -> Path:
    return topics_dir(paper_dir) / TOPICS_FILENAME


def resolved_path(paper_dir: Path) -> Path:
    return topics_dir(paper_dir) / RESOLVED_FILENAME


def report_path(paper_dir: Path) -> Path:
    return topics_dir(paper_dir) / REPORT_FILENAME


def _legacy_topics_path(paper_dir: Path) -> Path:
    return Path(paper_dir) / LEGACY_TOPICS_FILENAME


def _legacy_resolved_path(paper_dir: Path) -> Path:
    return Path(paper_dir) / LEGACY_RESOLVED_FILENAME


def migrate_legacy(paper_dir: Path) -> bool:
    """Move pre-subfolder flat files into the subfolder. Returns True if moved."""
    moved = False
    for legacy, new in (
        (_legacy_topics_path(paper_dir), topics_path(paper_dir)),
        (_legacy_resolved_path(paper_dir), resolved_path(paper_dir)),
    ):
        if legacy.exists() and not new.exists():
            new.parent.mkdir(parents=True, exist_ok=True)
            legacy.replace(new)
            moved = True
    return moved


def has_topics(paper_dir: Path) -> bool:
    return topics_path(paper_dir).exists() or _legacy_topics_path(paper_dir).exists()


def load_topics(paper_dir: Path) -> dict:
    """Load the raw topics extraction (subfolder, falling back to legacy path).

    Returns {} if missing/unparseable.
    """
    p = topics_path(paper_dir)
    if not p.exists():
        p = _legacy_topics_path(paper_dir)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        logger.warning("Failed to parse topics JSON: %s", p)
        return {}


def save_topics(
    paper_dir: Path,
    metadata: dict,
    genes: list[dict],
    pathways: list[dict],
) -> Path:
    """Write the raw topics extraction JSON into the subfolder. Returns the path."""
    output = {"metadata": metadata, "genes": genes, "pathways": pathways}
    p = topics_path(paper_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(output, indent=2, default=str), encoding="utf-8")
    # Clean up any legacy flat file now that the subfolder copy is authoritative.
    legacy = _legacy_topics_path(paper_dir)
    if legacy.exists():
        legacy.unlink()
    return p
