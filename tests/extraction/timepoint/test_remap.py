"""Tests for remap.py."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.remap import remap_value


def _seed_paper(paper_dir: Path, growth_phase: str) -> None:
    paper_dir.mkdir(parents=True, exist_ok=True)
    (paper_dir / "paperconfig.yaml").write_text(f"""publication:
  papername: X
  papermainpdf: data/x.pdf
  experiments:
    exp1: {{name: T, omics_type: RNASEQ}}
  supplementary_materials:
    tbl:
      type: csv
      filename: data/x.csv
      statistical_analyses:
        - id: DE_a
          experiment: exp1
          name_col: g
          logfc_col: lfc
          timepoint: "24h"
          timepoint_hours: 24
          growth_phase: {growth_phase}
""")
    (paper_dir / "extractions").mkdir(exist_ok=True)
    (paper_dir / "extractions" / "timepoint.json").write_text(json.dumps({
        "metadata": {"paper": "X", "status": "complete", "missing_analyses": []},
        "analyses": [{
            "analysis_id": "DE_a", "experiment_key": "exp1",
            "fields_requested": ["growth_phase"],
            "timepoint": "24h", "timepoint_hours": 24.0,
            "growth_phase": growth_phase,
            "self_assessment": "high", "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "p"}],
            "source_figures": [],
        }],
    }))


def test_remap_rewrites_paperconfig_and_json(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")

    remap_value(
        paper_dirs=[paper_dir],
        from_value="other:heat_acclimated",
        to_value="heat_acclimated",   # pretend validator has been updated
    )

    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "heat_acclimated"

    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "heat_acclimated"


def test_remap_appends_provenance_to_assessment_notes(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    note = ext["analyses"][0]["assessment_notes"]
    assert "Remapped from other:heat_acclimated" in note
    assert "heat_acclimated" in note


def test_remap_rejects_invalid_to_value(tmp_path):
    with pytest.raises(SystemExit):
        remap_value([tmp_path], "other:x", "malformed value with spaces")


def test_remap_accepts_other_to_other_rename(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclim")
    remap_value([paper_dir], "other:heat_acclim", "other:heat_acclimated")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "other:heat_acclimated"


def test_remap_merge_into_existing_canonical(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_shocked")
    remap_value([paper_dir], "other:heat_shocked", "acute_stress")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "acute_stress"


def test_remap_skips_analyses_without_match(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "exponential")
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated")
    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        data = yaml.load(f)
    # Unchanged
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "exponential"


def test_remap_dry_run_does_not_write(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")
    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        before = f.read()
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated", dry_run=True)
    with open(paper_dir / "paperconfig.yaml") as f:
        after = f.read()
    assert before == after
