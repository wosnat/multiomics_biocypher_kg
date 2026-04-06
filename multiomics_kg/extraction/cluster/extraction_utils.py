"""Extraction file I/O and data utilities.

Single source of truth for extraction file locations and JSON structure.
Used by the adapter, extract.py, and report/verify tooling.
"""
import json
import logging
import re
from pathlib import Path

import pandas as pd

from multiomics_kg.utils.paperconfig_utils import (
    load_all_paperconfigs,
    iter_cluster_tables,
)

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


# ── Data loading ──


def find_all_entries() -> list[tuple[Path, str, dict, dict, dict]]:
    """Discover all gene_clusters entries across all paperconfigs.
    Returns list of (paper_dir, entry_key, table_config, pub_config, extraction_config).
    """
    entries = []
    for pc_path, pc in load_all_paperconfigs():
        paper_dir = pc_path.parent
        pub = pc.get("publication", {})
        extraction = pc.get("extraction", {})
        for entry_key, table in iter_cluster_tables(pc):
            entries.append((paper_dir, entry_key, table, pub, extraction))
    return entries


def load_cluster_summaries(table_config: dict) -> dict[str, dict]:
    """Load cluster gene counts and sample genes from CSV.
    Returns {cluster_key: {gene_count: int, sample_genes: list[str]}}.
    """
    csv_path = Path(table_config["filename"])
    gene_id_col = table_config["gene_id_col"]
    cluster_col = table_config["cluster_col"]
    skip_rows = table_config.get("skip_rows", 0)
    df = pd.read_csv(csv_path, skiprows=skip_rows if skip_rows else None)

    clusters = {}
    for cid, group in df.groupby(cluster_col):
        ckey = _cluster_val_to_str(cid)
        if not ckey:
            continue
        genes = group[gene_id_col].dropna().tolist()
        clusters[ckey] = {
            "gene_count": len(genes),
            "sample_genes": [str(g) for g in genes[:5]],
        }
    return clusters


def _cluster_val_to_str(val) -> str:
    """Convert cluster column value to clean string key."""
    if pd.isna(val):
        return ""
    if isinstance(val, float) and val == int(val):
        return str(int(val))
    return str(val).strip()


# ── Cluster key matching ──


def match_cluster_keys(
    parsed_clusters: list[dict],
    expected_keys: set[str],
) -> tuple[dict[str, dict], list[dict]]:
    """Map model output to actual cluster keys.
    Tries: name regex -> id suffix -> case-insensitive -> positional fallback.
    Returns (matched: {key: extraction_dict}, unmatched: [extraction_dict]).
    """
    remaining_keys = set(expected_keys)
    matched = {}
    unmatched_pass1 = []

    for c in parsed_clusters:
        key = _try_match_key(c, remaining_keys)
        if key:
            matched[key] = c
            remaining_keys.discard(key)
        else:
            unmatched_pass1.append(c)

    # Positional fallback: if unmatched count == remaining keys count
    if unmatched_pass1 and len(unmatched_pass1) == len(remaining_keys):
        remaining_sorted = sorted(remaining_keys, key=lambda x: (not x.isdigit(), x))
        for c, key in zip(unmatched_pass1, remaining_sorted):
            matched[key] = c
        return matched, []

    return matched, unmatched_pass1


def _try_match_key(extraction: dict, expected_keys: set[str]) -> str | None:
    """Try to match a single extraction to a cluster key."""
    name = extraction.get("name", "")
    ext_id = extraction.get("id", "")

    # Try "cluster KEY" in name
    m = re.search(r"cluster\s+(\S+)", name, re.IGNORECASE)
    if m:
        candidate = m.group(1).rstrip(",.):")
        if candidate in expected_keys:
            return candidate
        for k in expected_keys:
            if k.lower() == candidate.lower():
                return k

    # Try suffix of id
    m = re.search(r"_([^_]+)$", ext_id)
    if m:
        candidate = m.group(1)
        if candidate in expected_keys:
            return candidate
        for k in expected_keys:
            if k.lower() == candidate.lower():
                return k

    return None
