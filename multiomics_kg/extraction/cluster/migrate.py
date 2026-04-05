"""One-time migration of legacy cluster_extraction_*.json into versioned run structure."""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from multiomics_kg.extraction.cluster.run_manager import RunManager

logger = logging.getLogger(__name__)


def migrate_extraction_json(
    paper_dir: Path, entry_key: str, cache_dir: Path
) -> Optional[Path]:
    """Migrate a legacy cluster_extraction_{entry_key}.json into the new run structure."""
    legacy_path = paper_dir / f"cluster_extraction_{entry_key}.json"
    if not legacy_path.exists():
        return None

    with open(legacy_path) as f:
        data = json.load(f)

    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.create_run()

    if "stage1_merged" in data:
        rm.write_stage(run_dir, 1, data["stage1_merged"])
    if "stage2_results" in data:
        rm.write_stage(run_dir, 2, data["stage2_results"])
    if "stage3_validation" in data:
        rm.write_stage(run_dir, 3, data["stage3_validation"])

    metadata = dict(data.get("metadata", {}))
    metadata["migrated_from"] = str(legacy_path)
    metadata["migrated_at"] = datetime.now().isoformat(timespec="seconds")
    with open(run_dir / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    rm.finalize_run(run_dir)
    logger.info("Migrated %s → %s", legacy_path, run_dir)
    return run_dir


def migrate_all_tolonen(project_root: Path) -> list[Path]:
    """Migrate both Tolonen 2006 extraction JSONs."""
    paper_dir = project_root / "data" / "Prochlorococcus" / "papers_and_supp" / "tolonen 2006"
    cache_dir = paper_dir / ".extraction_cache"
    results = []
    for entry_key in ["med4_kmeans_nstarvation", "mit9313_kmeans_nstarvation"]:
        run_dir = migrate_extraction_json(paper_dir, entry_key, cache_dir)
        if run_dir:
            results.append(run_dir)
    return results
