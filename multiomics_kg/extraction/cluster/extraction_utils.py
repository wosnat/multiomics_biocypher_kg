"""Extraction file I/O and data utilities.

Single source of truth for extraction file locations and JSON structure.
Used by the adapter, extract.py, and report/verify tooling.
"""
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

EXTRACTIONS_DIR = "cluster_extractions"


def load_extraction(paper_dir: Path, entry_key: str) -> dict[str, dict]:
    """Load cluster extraction data from JSON.
    Returns the clusters dict: {cluster_key: {id, name, functional_description, ...}}.
    Returns {} if file is missing or has wrong format.
    """
    json_path = Path(paper_dir) / EXTRACTIONS_DIR / f"{entry_key}.json"
    if not json_path.exists():
        return {}
    try:
        data = json.loads(json_path.read_text())
    except Exception:
        logger.warning("Failed to parse extraction JSON: %s", json_path)
        return {}
    if "clusters" not in data:
        logger.warning("Extraction JSON missing 'clusters' key: %s", json_path)
        return {}
    return data["clusters"]


def save_extraction(paper_dir: Path, entry_key: str, metadata: dict, clusters: dict[str, dict]) -> Path:
    """Write extraction results as JSON + markdown summary.
    Creates cluster_extractions/ dir if needed. Returns path to JSON file.
    """
    ext_dir = Path(paper_dir) / EXTRACTIONS_DIR
    ext_dir.mkdir(parents=True, exist_ok=True)

    output = {"metadata": metadata, "clusters": clusters}
    json_path = ext_dir / f"{entry_key}.json"
    json_path.write_text(json.dumps(output, indent=2, default=str))

    md_lines = [f"# {metadata.get('paper', '')} — {entry_key}\n"]
    for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
        c = clusters[key]
        direction = c.get("direction", "")
        assessment = c.get("self_assessment", "")
        md_lines.append(f"## Cluster {key} | {direction} | {assessment}\n")
        md_lines.append(f"**Name:** {c.get('name', '')}")
        md_lines.append(f"**Enrichment:** {c.get('enrichment_category', '')} "
                        f"(p={c.get('enrichment_pvalue', 'N/A')})")
        md_lines.append(f"**Functional:** {c.get('functional_description', '')}\n")
        md_lines.append(f"**Behavioral:** {c.get('behavioral_description', '')}\n")
        notes = c.get("confidence_notes", "")
        if notes:
            md_lines.append(f"**Notes:** {notes}\n")
        quotes = c.get("supporting_quotes", [])
        if quotes:
            md_lines.append("**Quotes:**")
            for q in quotes:
                md_lines.append(f"- [{q.get('location', '')}] {q.get('quote', '')}")
            md_lines.append("")
    md_path = ext_dir / f"{entry_key}.md"
    md_path.write_text("\n".join(md_lines))

    return json_path


def get_cluster_data(clusters: dict, key) -> dict:
    """Get data for a single cluster by key. Returns {} if not found."""
    return clusters.get(str(key), {})


def list_extraction_files(paper_dir: Path) -> list[str]:
    """List entry keys that have extraction JSONs in paper_dir/cluster_extractions/."""
    ext_dir = Path(paper_dir) / EXTRACTIONS_DIR
    if not ext_dir.is_dir():
        return []
    return sorted(p.stem for p in ext_dir.glob("*.json"))
