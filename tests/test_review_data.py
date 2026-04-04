import json
from pathlib import Path

import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager
from multiomics_kg.review.review_data import (
    compute_entry_status_color,
    export_issue_report,
    load_entry_summary,
    scan_papers_with_clusters,
)


def _setup_paper_with_extraction(tmp_path, paper_name, entry_key, clusters):
    """Create a minimal paperconfig + extraction cache for testing."""
    paper_dir = tmp_path / paper_name
    paper_dir.mkdir(parents=True)

    config = {
        "publication": {
            "papername": paper_name,
            "papermainpdf": "",
            "supplementary_materials": {
                entry_key: {
                    "type": "gene_clusters",
                    "filename": str(paper_dir / "clusters.csv"),
                    "organism": "Test Organism",
                    "gene_id_col": "gene",
                    "cluster_col": "cluster",
                }
            },
        }
    }
    with open(paper_dir / "paperconfig.yaml", "w") as f:
        yaml.dump(config, f)

    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.create_run()

    stage1 = {
        str(k): {"enrichment_category": [{"value": f"cat_{k}"}]} for k in clusters
    }
    stage2 = {
        str(k): {"id": f"cluster_{k}", "functional_description": f"desc_{k}"}
        for k in clusters
    }
    stage3 = {str(k): {"verdict": v} for k, v in clusters.items()}

    rm.write_stage(run_dir, 1, stage1)
    rm.write_stage(run_dir, 2, stage2)
    rm.write_stage(run_dir, 3, stage3)

    return paper_dir, rm, run_dir


def test_scan_papers_finds_cluster_entries(tmp_path):
    paper_dir, _, _ = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass", 2: "warn"}
    )
    results = scan_papers_with_clusters([paper_dir / "paperconfig.yaml"])
    assert len(results) == 1
    assert results[0]["paper_name"] == "Paper1"
    assert results[0]["entry_key"] == "entry_a"


def test_load_entry_summary(tmp_path):
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass", 2: "warn", 3: "fail"}
    )
    summary = load_entry_summary(paper_dir, "entry_a")
    assert summary["cluster_count"] == 3
    assert summary["verdicts"]["pass"] == 1
    assert summary["verdicts"]["warn"] == 1
    assert summary["verdicts"]["fail"] == 1


def test_compute_entry_status_color_all_approved(tmp_path):
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    review = {
        "1": {
            "status": "approve",
            "reviewed_in_run": run_dir.name,
            "input_hash": rm.compute_input_hash("1", rm.read_stage(run_dir, 1)),
            "issues": [],
            "failing_stages": [],
            "notes": "",
            "edited_fields": {},
            "reviewed_at": "2026-04-04T12:00:00",
        },
    }
    rm.write_stage(run_dir, 4, review)
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "green"


def test_compute_entry_status_color_carried_forward(tmp_path):
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    review = {
        "1": {
            "status": "approve",
            "reviewed_in_run": "some-old-run",
            "input_hash": rm.compute_input_hash("1", rm.read_stage(run_dir, 1)),
            "issues": [],
            "failing_stages": [],
            "notes": "",
            "edited_fields": {},
            "reviewed_at": "2026-04-04T12:00:00",
        },
    }
    rm.write_stage(run_dir, 4, review)
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "light_green"


def test_compute_entry_status_color_unreviewed(tmp_path):
    paper_dir, _, _ = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "red"


def test_export_issue_report(tmp_path):
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "warn"}
    )
    review = {
        "1": {
            "status": "flag-issue",
            "issues": ["hallucinated", "cross_contamination"],
            "failing_stages": ["stage1_visual"],
            "notes": "Bad quote",
            "input_hash": "abc",
            "edited_fields": {},
            "reviewed_at": "2026-04-04T12:00:00",
            "reviewed_in_run": run_dir.name,
        },
    }
    rm.write_stage(run_dir, 4, review)
    # Write metadata for paper name lookup
    with open(run_dir / "metadata.json", "w") as f:
        json.dump({"paper": "Paper1"}, f)
    rows = export_issue_report([(paper_dir, "entry_a")])
    assert len(rows) == 1
    assert rows[0]["issues"] == "hallucinated,cross_contamination"
    assert rows[0]["failing_stages"] == "stage1_visual"
    assert rows[0]["notes"] == "Bad quote"
