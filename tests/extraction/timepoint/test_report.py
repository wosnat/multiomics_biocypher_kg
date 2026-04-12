"""Tests for report.py — aggregates extraction JSONs into a markdown report."""
import json
from pathlib import Path

from multiomics_kg.extraction.timepoint.report import (
    aggregate_reports,
    render_markdown,
)


def _seed(paper_dir: Path, paper: str, analyses: list[dict], status: str = "complete"):
    paper_dir.mkdir(parents=True, exist_ok=True)
    (paper_dir / "extractions").mkdir(exist_ok=True)
    (paper_dir / "extractions" / "timepoint.json").write_text(json.dumps({
        "metadata": {"paper": paper, "status": status, "missing_analyses": []},
        "analyses": analyses,
    }))


def test_aggregate_counts_self_assessments(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "exponential", "self_assessment": "high",
         "timepoint": "24h", "timepoint_hours": 24, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    _seed(tmp_path / "p2", "P2", [
        {"analysis_id": "b", "growth_phase": "unknown", "self_assessment": "low",
         "timepoint": "unknown", "timepoint_hours": None, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1", tmp_path / "p2"])
    assert agg["self_assessment_counts"] == {"high": 1, "low": 1}


def test_aggregate_other_slug_frequencies(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "other:heat_acclimated",
         "self_assessment": "medium", "timepoint": "48h", "timepoint_hours": 48,
         "assessment_notes": "", "supporting_quotes": [], "source_figures": []},
    ])
    _seed(tmp_path / "p2", "P2", [
        {"analysis_id": "b", "growth_phase": "other:heat_acclimated",
         "self_assessment": "high", "timepoint": "72h", "timepoint_hours": 72,
         "assessment_notes": "", "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1", tmp_path / "p2"])
    assert agg["other_slug_counts"] == {"other:heat_acclimated": 2}


def test_aggregate_collects_unknown_with_evidence(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "unknown", "self_assessment": "low",
         "timepoint": "unknown", "timepoint_hours": None,
         "assessment_notes": "Paper does not describe phase",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1"])
    assert len(agg["unknowns"]) == 1
    assert agg["unknowns"][0]["analysis_id"] == "a"
    assert "Paper does not describe phase" in agg["unknowns"][0]["notes"]


def test_aggregate_includes_partial_missing(tmp_path):
    _seed(tmp_path / "p1", "P1", [], status="partial")
    d = tmp_path / "p1" / "extractions" / "timepoint.json"
    data = json.loads(d.read_text())
    data["metadata"]["missing_analyses"] = [{"analysis_id": "x", "reason": "not_returned"}]
    d.write_text(json.dumps(data))
    agg = aggregate_reports([tmp_path / "p1"])
    assert len(agg["partial_extractions"]) == 1
    assert agg["partial_extractions"][0]["paper"] == "P1"


def test_render_markdown_has_required_sections(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "other:x", "self_assessment": "medium",
         "timepoint": "24h", "timepoint_hours": 24, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1"])
    md = render_markdown(agg)
    assert "# Timepoint Extraction Report" in md
    assert "## Self-assessment distribution" in md
    assert "## other:* slugs" in md
    assert "## Unknowns" in md
    assert "## Partial extractions" in md
