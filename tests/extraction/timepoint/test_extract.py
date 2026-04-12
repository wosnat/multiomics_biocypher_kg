"""Tests for extract.py. LLM calls are mocked — no network."""
import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extract import (
    build_background,
    build_targets,
    extract_one_paper,
)


@pytest.fixture
def paper_dir(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Fake 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: N-starvation
      organism: Prochlorococcus MED4
      treatment_type: [nitrogen]
      omics_type: RNASEQ
  supplementary_materials:
    tbl:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_n_24h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_24h
        - id: DE_n_48h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_48h
""")
    # papermainpdf doesn't need to exist for build_background; extract_one_paper
    # will need it. Create a dummy file.
    (tmp_path.parent / "data").mkdir(parents=True, exist_ok=True)
    (tmp_path.parent / "data" / "fake.pdf").write_text("dummy")
    return tmp_path


def test_build_background_includes_experiments_block(paper_dir):
    bg = build_background(paper_dir / "paperconfig.yaml")
    assert bg["papername"] == "Fake 2026"
    assert "exp1" in bg["experiments"]
    assert bg["experiments"]["exp1"]["omics_type"] == "RNASEQ"


def test_build_targets_computes_fields_requested(paper_dir):
    analyses = [
        {"id": "DE_n_24h", "experiment": "exp1", "name_col": "gene",
         "logfc_col": "log2fc_24h"},
        {"id": "DE_n_48h", "experiment": "exp1", "name_col": "gene",
         "logfc_col": "log2fc_48h", "timepoint": "48h", "timepoint_hours": 48},
    ]
    targets = build_targets(analyses, validate=False)
    # First analysis: all three requested
    assert set(targets[0]["fields_requested"]) == {
        "timepoint", "timepoint_hours", "growth_phase",
    }
    assert targets[0]["logfc_col"] == "log2fc_24h"
    # Second analysis: only growth_phase requested (timepoint + hours set)
    assert targets[1]["fields_requested"] == ["growth_phase"]


def test_build_targets_dropped_when_all_populated(paper_dir):
    analyses = [{
        "id": "DE_done",
        "experiment": "exp1",
        "name_col": "gene", "logfc_col": "lfc",
        "timepoint": "24h", "timepoint_hours": 24, "growth_phase": "exponential",
    }]
    targets = build_targets(analyses, validate=False)
    assert targets == []


def test_extract_one_paper_happy_path(paper_dir, monkeypatch):
    fake_llm = MagicMock()
    fake_llm.call.return_value = (
        {
            "analyses": [
                {
                    "analysis_id": "DE_n_24h",
                    "timepoint": "24h",
                    "timepoint_hours": 24.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high",
                    "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "Methods"}],
                    "source_figures": [],
                },
                {
                    "analysis_id": "DE_n_48h",
                    "timepoint": "48h",
                    "timepoint_hours": 48.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high",
                    "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "Methods"}],
                    "source_figures": [],
                },
            ]
        },
        {"input_tokens": 100, "output_tokens": 200, "model": "fake-model"},
    )
    json_path = extract_one_paper(
        paper_dir,
        llm_client=fake_llm,
        validate=False,
    )
    assert json_path.exists()
    data = json.loads(json_path.read_text())
    assert data["metadata"]["status"] == "complete"
    assert data["metadata"]["missing_analyses"] == []
    assert data["metadata"]["model"] == "fake-model"
    assert len(data["analyses"]) == 2
    # experiment_key injected by extract.py
    assert all(a["experiment_key"] == "exp1" for a in data["analyses"])
    # fields_requested injected
    assert all("fields_requested" in a for a in data["analyses"])


def test_extract_one_paper_partial_marks_status(paper_dir):
    fake_llm = MagicMock()
    fake_llm.call.return_value = (
        {
            "analyses": [
                # Only 1 of 2 returned — DE_n_48h missing → partial
                {
                    "analysis_id": "DE_n_24h",
                    "timepoint": "24h", "timepoint_hours": 24.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high", "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "M"}],
                    "source_figures": [],
                },
            ]
        },
        {"input_tokens": 50, "output_tokens": 100, "model": "fake-model"},
    )
    json_path = extract_one_paper(paper_dir, llm_client=fake_llm, validate=False)
    data = json.loads(json_path.read_text())
    assert data["metadata"]["status"] == "partial"
    assert data["metadata"]["missing_analyses"] == [
        {"analysis_id": "DE_n_48h", "reason": "not_returned"},
    ]
    assert len(data["analyses"]) == 1


def test_extract_one_paper_skips_when_no_fields_requested(paper_dir):
    # Pre-populate both analyses
    pc_path = paper_dir / "paperconfig.yaml"
    yaml = YAML()
    with open(pc_path) as f:
        data = yaml.load(f)
    for a in data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"]:
        a["timepoint"] = "24h"
        a["timepoint_hours"] = 24
        a["growth_phase"] = "exponential"
    with open(pc_path, "w") as f:
        yaml.dump(data, f)

    fake_llm = MagicMock()
    result = extract_one_paper(paper_dir, llm_client=fake_llm, validate=False)
    # No LLM call
    fake_llm.call.assert_not_called()
    # No JSON written (or return None, per spec)
    assert result is None


from unittest.mock import MagicMock


def test_resume_only_extracts_missing_analyses(paper_dir):
    # Seed a partial extraction with DE_n_24h valid and DE_n_48h missing
    sig = __import__(
        "multiomics_kg.extraction.timepoint.extraction_utils", fromlist=["compute_paperconfig_signature"]
    ).compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    from multiomics_kg.extraction.timepoint.extraction_utils import save_extraction_json
    save_extraction_json(paper_dir, {
        "paper": "Fake 2026", "doi": "", "status": "partial",
        "paperconfig_signature": sig, "input_tokens": 0, "output_tokens": 0,
        "missing_analyses": [{"analysis_id": "DE_n_48h", "reason": "not_returned"}],
    }, [{
        "analysis_id": "DE_n_24h", "experiment_key": "exp1",
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
        "timepoint": "24h", "timepoint_hours": 24.0, "growth_phase": "exponential",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "M"}],
        "source_figures": [],
    }])

    fake_llm = MagicMock()
    fake_llm.call.return_value = ({
        "analyses": [{
            "analysis_id": "DE_n_48h",
            "timepoint": "48h", "timepoint_hours": 48.0,
            "growth_phase": "nutrient_limited",
            "self_assessment": "high", "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "M"}],
            "source_figures": [],
        }]
    }, {"input_tokens": 50, "output_tokens": 100, "model": "fake-model"})

    from multiomics_kg.extraction.timepoint.extract import extract_one_paper
    extract_one_paper(paper_dir, llm_client=fake_llm, validate=False, resume=True)

    # LLM should have been called with only the missing analysis in targets
    call_args = fake_llm.call.call_args
    prompt = call_args.kwargs.get("prompt") or call_args.args[0]
    assert "DE_n_48h" in prompt
    assert "DE_n_24h" not in prompt  # valid row not re-asked

    # Final JSON should now be complete and contain both rows
    import json
    final = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert final["metadata"]["status"] == "complete"
    assert {r["analysis_id"] for r in final["analyses"]} == {"DE_n_24h", "DE_n_48h"}


def test_retry_ignores_existing_json(paper_dir):
    # Seed a complete JSON
    sig = __import__(
        "multiomics_kg.extraction.timepoint.extraction_utils", fromlist=["compute_paperconfig_signature"]
    ).compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    from multiomics_kg.extraction.timepoint.extraction_utils import save_extraction_json
    save_extraction_json(paper_dir, {
        "paper": "Fake 2026", "status": "complete", "paperconfig_signature": sig,
        "missing_analyses": [], "input_tokens": 0, "output_tokens": 0,
    }, [])

    fake_llm = MagicMock()
    fake_llm.call.return_value = ({"analyses": []}, {"input_tokens": 0, "output_tokens": 0, "model": "m"})

    from multiomics_kg.extraction.timepoint.extract import extract_one_paper
    # With retry=True, the existing JSON is ignored and every analysis is requested.
    # (paperconfig has two analyses, all without timepoint fields → both requested)
    extract_one_paper(paper_dir, llm_client=fake_llm, validate=False, retry=True)
    call_args = fake_llm.call.call_args
    prompt = call_args.kwargs.get("prompt") or call_args.args[0]
    assert "DE_n_24h" in prompt
    assert "DE_n_48h" in prompt
