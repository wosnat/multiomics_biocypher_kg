"""Manages versioned extraction run directories and stage file I/O."""

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Optional


class RunManager:
    """Manages the .extraction_cache/{entry_key}/ directory for one paperconfig entry.

    Directory layout:
        {cache_root}/{entry_key}/
            runs/
                2026-04-04T10-30-00/
                    stage1_merged.json
                    stage2_results.json
                    stage3_validation.json
                    stage4_review.json
                ...
            current -> runs/2026-04-04T10-30-00  (symlink)

        {cache_root}/shared/
            pdf_text.json
            pages/
            embeddings.npz
            rag_chunks.json
            tables/
    """

    STAGE_FILES = {
        1: "stage1_merged.json",
        2: "stage2_results.json",
        3: "stage3_validation.json",
        4: "stage4_review.json",
    }

    def __init__(self, cache_root: Path, entry_key: str):
        self.cache_root = Path(cache_root)
        self.entry_key = entry_key
        self.entry_dir = self.cache_root / entry_key
        self.runs_dir = self.entry_dir / "runs"
        self.current_link = self.entry_dir / "current"
        self.shared_dir = self.cache_root / "shared"

    def create_run(self) -> Path:
        """Create a new timestamped run directory and update 'current' symlink."""
        timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        run_dir = self.runs_dir / timestamp
        run_dir.mkdir(parents=True, exist_ok=True)
        # Update symlink atomically
        tmp_link = self.current_link.with_suffix(".tmp")
        if tmp_link.exists() or tmp_link.is_symlink():
            tmp_link.unlink()
        tmp_link.symlink_to(os.path.relpath(run_dir, self.entry_dir))
        os.replace(tmp_link, self.current_link)
        return run_dir

    def get_current_run(self) -> Optional[Path]:
        """Return the current run directory, or None if no runs exist."""
        if self.current_link.is_symlink():
            target = self.current_link.resolve()
            if target.exists():
                return target
        return None

    def list_runs(self) -> list[Path]:
        """Return all run directories, oldest first."""
        if not self.runs_dir.exists():
            return []
        return sorted(
            [d for d in self.runs_dir.iterdir() if d.is_dir()],
            key=lambda p: p.name,
        )

    def write_stage(self, run_dir: Path, stage: int, data: dict) -> Path:
        """Write stage data to the appropriate JSON file in run_dir."""
        filename = self.STAGE_FILES[stage]
        path = run_dir / filename
        with open(path, "w") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        return path

    def read_stage(self, run_dir: Path, stage: int) -> dict:
        """Read stage data from run_dir. Returns {} if file missing."""
        filename = self.STAGE_FILES[stage]
        path = run_dir / filename
        if not path.exists():
            return {}
        with open(path) as f:
            return json.load(f)

    def read_current_stage(self, stage: int) -> dict:
        """Read stage data from the current run. Returns {} if no current run."""
        run_dir = self.get_current_run()
        if run_dir is None:
            return {}
        return self.read_stage(run_dir, stage)

    @staticmethod
    def compute_input_hash(cluster_key: str, stage1_merged: dict) -> str:
        """Compute a deterministic hash of stage1 inputs for one cluster."""
        cluster_data = stage1_merged.get(str(cluster_key), {})
        serialized = json.dumps(cluster_data, sort_keys=True, ensure_ascii=False)
        return hashlib.sha256(serialized.encode()).hexdigest()[:16]

    def copy_forward_reviews(
        self, prev_run: Path, new_run: Path, new_stage1: dict
    ) -> None:
        """Copy review decisions from prev_run to new_run.

        - If input hash matches: review carries forward as-is
        - If input hash differs: review status set to 'stale'
        - Clusters not in new_stage1 are dropped
        """
        prev_review = self.read_stage(prev_run, 4)
        if not prev_review:
            return

        new_review = {}
        for cluster_key, review_entry in prev_review.items():
            if cluster_key not in new_stage1:
                continue
            new_hash = self.compute_input_hash(cluster_key, new_stage1)
            old_hash = review_entry.get("input_hash", "")
            entry_copy = dict(review_entry)
            if new_hash != old_hash:
                entry_copy["status"] = "stale"
            entry_copy["input_hash"] = new_hash
            new_review[cluster_key] = entry_copy

        if new_review:
            self.write_stage(new_run, 4, new_review)
