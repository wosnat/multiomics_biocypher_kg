"""Data loading and scanning for the cluster review UI."""

import json
import logging
from pathlib import Path
from typing import Optional

import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager

logger = logging.getLogger(__name__)


def scan_papers_with_clusters(paperconfig_paths: list[Path]) -> list[dict]:
    """Scan paperconfigs and return entries that have gene_clusters.

    Returns list of dicts with keys: paper_name, entry_key, paperconfig_path,
    paper_dir, organism.
    """
    results = []
    for pc_path in paperconfig_paths:
        pc_path = Path(pc_path)
        if not pc_path.exists():
            continue
        with open(pc_path) as f:
            config = yaml.safe_load(f)
        paper_name = config.get("publication", {}).get("papername", pc_path.parent.name)
        supp = config.get("publication", {}).get("supplementary_materials", {})
        for key, entry in supp.items():
            if entry.get("type") == "gene_clusters":
                results.append({
                    "paper_name": paper_name,
                    "entry_key": key,
                    "paperconfig_path": pc_path,
                    "paper_dir": pc_path.parent,
                    "organism": entry.get("organism", ""),
                })
    return results


def load_entry_summary(paper_dir: Path, entry_key: str) -> dict:
    """Load summary stats for one extraction entry.

    Returns dict with cluster_count, verdicts, review_statuses, has_run.
    """
    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.get_current_run()

    if run_dir is None:
        return {"cluster_count": 0, "verdicts": {}, "review_statuses": {}, "has_run": False}

    stage2 = rm.read_stage(run_dir, 2)
    stage3 = rm.read_stage(run_dir, 3)
    stage4 = rm.read_stage(run_dir, 4)

    cluster_keys = sorted(set(list(stage2.keys()) + list(stage3.keys())))
    verdicts: dict[str, int] = {}
    review_statuses: dict[str, int] = {}

    for ck in cluster_keys:
        v = stage3.get(ck, {}).get("verdict", "none")
        verdicts[v] = verdicts.get(v, 0) + 1
        r = stage4.get(ck, {}).get("status", "unreviewed")
        review_statuses[r] = review_statuses.get(r, 0) + 1

    return {
        "cluster_count": len(cluster_keys),
        "verdicts": verdicts,
        "review_statuses": review_statuses,
        "has_run": True,
    }


def compute_entry_status_color(paper_dir: Path, entry_key: str) -> str:
    """Compute rollup color for an entry based on review state.

    Returns: 'green', 'light_green', 'yellow', 'red'
    """
    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.get_current_run()

    if run_dir is None:
        return "red"

    stage2 = rm.read_stage(run_dir, 2)
    stage4 = rm.read_stage(run_dir, 4)
    cluster_keys = list(stage2.keys())

    if not cluster_keys:
        return "red"

    all_approved = True
    any_stale = False

    for ck in cluster_keys:
        review = stage4.get(ck, {})
        status = review.get("status", "")
        if status not in ("approve", "edit"):
            all_approved = False
        if status == "stale":
            any_stale = True

    if not all_approved:
        if any_stale:
            return "yellow"
        return "red"

    # All approved -- check if reviewed in current run
    has_all_current = all(
        stage4.get(ck, {}).get("reviewed_in_run") == run_dir.name
        for ck in cluster_keys
    )
    if has_all_current:
        return "green"
    return "light_green"


def export_issue_report(entries: list[tuple[Path, str]]) -> list[dict]:
    """Export issue report rows across all entries.

    Returns list of dicts suitable for CSV export.
    """
    rows = []
    for paper_dir, entry_key in entries:
        cache_dir = paper_dir / ".extraction_cache"
        rm = RunManager(cache_dir, entry_key)
        run_dir = rm.get_current_run()
        if run_dir is None:
            continue

        stage3 = rm.read_stage(run_dir, 3)
        stage4 = rm.read_stage(run_dir, 4)
        metadata_path = run_dir / "metadata.json"
        paper_name = entry_key
        if metadata_path.exists():
            with open(metadata_path) as f:
                meta = json.load(f)
            paper_name = meta.get("paper", entry_key)

        for ck in sorted(stage4.keys()):
            review = stage4[ck]
            if review.get("status") in ("approve",):
                continue
            verdict = stage3.get(ck, {}).get("verdict", "")
            rows.append({
                "paper": paper_name,
                "entry": entry_key,
                "cluster": ck,
                "verdict": verdict,
                "review_status": review.get("status", ""),
                "issues": ",".join(review.get("issues", [])),
                "failing_stages": ",".join(review.get("failing_stages", [])),
                "notes": review.get("notes", ""),
            })
    return rows
