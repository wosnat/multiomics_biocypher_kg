"""Tests for timepoint extraction utilities."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    EXTRACTIONS_DIR,
    compute_paperconfig_signature,
    find_analyses,
    iter_paperconfigs,
    load_extraction_json,
    save_extraction_json,
)


@pytest.fixture
def paper_dir(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Test 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: Test
      organism: Prochlorococcus MED4
      omics_type: RNASEQ
  supplementary_materials:
    tbl_a:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_foo_24h
          experiment: exp1
          timepoint: "24h"
          timepoint_hours: 24
          growth_phase: exponential
          name_col: gene
          logfc_col: log2fc_24h
        - id: DE_foo_48h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_48h
""")
    return tmp_path


def test_find_analyses_returns_all(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    ids = [a["id"] for a in analyses]
    assert ids == ["DE_foo_24h", "DE_foo_48h"]


def test_find_analyses_includes_logfc_col(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    assert analyses[0]["logfc_col"] == "log2fc_24h"
    assert analyses[1]["logfc_col"] == "log2fc_48h"


def test_find_analyses_includes_experiment_key(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    assert all(a["experiment"] == "exp1" for a in analyses)


def test_save_and_load_extraction_json(paper_dir):
    metadata = {"paper": "Test 2026", "status": "complete", "missing_analyses": []}
    analyses = [
        {"analysis_id": "DE_foo_24h", "timepoint": "24h", "growth_phase": "exponential"},
    ]
    path = save_extraction_json(paper_dir, metadata, analyses)
    assert path == paper_dir / EXTRACTIONS_DIR / "timepoint.json"
    assert path.exists()

    loaded = load_extraction_json(paper_dir)
    assert loaded["metadata"]["paper"] == "Test 2026"
    assert loaded["analyses"][0]["analysis_id"] == "DE_foo_24h"


def test_load_extraction_json_returns_none_if_missing(tmp_path):
    assert load_extraction_json(tmp_path) is None


def test_compute_signature_stable_across_runs(paper_dir):
    sig1 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    sig2 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    assert sig1 == sig2
    assert len(sig1) == 40  # sha1 hex


def test_compute_signature_changes_on_edit(paper_dir):
    sig1 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    # Add a new analysis
    pc = paper_dir / "paperconfig.yaml"
    text = pc.read_text()
    text = text.replace(
        "          logfc_col: log2fc_48h\n",
        "          logfc_col: log2fc_48h\n"
        "        - id: DE_foo_72h\n"
        "          experiment: exp1\n"
        "          name_col: gene\n"
        "          logfc_col: log2fc_72h\n",
    )
    pc.write_text(text)
    sig2 = compute_paperconfig_signature(pc)
    assert sig1 != sig2


def test_iter_paperconfigs_reads_list_files(tmp_path, monkeypatch):
    list_file = tmp_path / "paperconfig_files.txt"
    p1 = tmp_path / "paper1" / "paperconfig.yaml"
    p1.parent.mkdir()
    p1.write_text("publication: {papername: p1}\n")
    list_file.write_text(f"# comment\n{p1}\n")
    paths = list(iter_paperconfigs([list_file]))
    assert paths == [p1]


def test_compute_fields_requested_all_missing():
    from multiomics_kg.extraction.timepoint.extraction_utils import (
        compute_fields_requested,
    )
    analysis = {"id": "x", "name_col": "g", "logfc_col": "lfc"}
    assert compute_fields_requested(analysis, validate=False) == [
        "timepoint", "timepoint_hours", "growth_phase",
    ]


def test_compute_fields_requested_timepoint_hours_null_is_request():
    from multiomics_kg.extraction.timepoint.extraction_utils import (
        compute_fields_requested,
    )
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": None,
        "growth_phase": "exponential",
    }
    # null is a "needs filling" signal — request timepoint_hours only
    assert compute_fields_requested(analysis, validate=False) == ["timepoint_hours"]


def test_compute_fields_requested_empty_string_is_request():
    from multiomics_kg.extraction.timepoint.extraction_utils import (
        compute_fields_requested,
    )
    analysis = {
        "id": "x", "timepoint": "", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=False) == ["timepoint"]


def test_compute_fields_requested_all_set_returns_empty():
    from multiomics_kg.extraction.timepoint.extraction_utils import (
        compute_fields_requested,
    )
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=False) == []


def test_compute_fields_requested_validate_mode_always_all_three():
    from multiomics_kg.extraction.timepoint.extraction_utils import (
        compute_fields_requested,
    )
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=True) == [
        "timepoint", "timepoint_hours", "growth_phase",
    ]


# ---------------------------------------------------------------------------
# validate_llm_payload tests
# ---------------------------------------------------------------------------

def test_validate_llm_payload_all_valid():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {
        "DE_a": ["timepoint", "timepoint_hours", "growth_phase"],
        "DE_b": ["timepoint", "timepoint_hours", "growth_phase"],
    }
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint": "24h",
                "timepoint_hours": 24,
                "growth_phase": "exponential",
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": ["quote1"],
                "source_figures": ["Table 1"],
            },
            {
                "analysis_id": "DE_b",
                "timepoint": "48h",
                "timepoint_hours": 48,
                "growth_phase": "stationary",
                "self_assessment": "medium",
                "assessment_notes": "ok",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 2
    assert missing == []


def test_validate_llm_payload_detects_missing_analysis():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {
        "DE_a": ["timepoint", "growth_phase"],
        "DE_b": ["timepoint", "growth_phase"],
    }
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint": "24h",
                "growth_phase": "exponential",
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 1
    assert {"analysis_id": "DE_b", "reason": "not_returned"} in missing


def test_validate_llm_payload_rejects_invalid_growth_phase():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["growth_phase"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "growth_phase": "maybe_stressed",
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert len(missing) == 1
    assert missing[0]["analysis_id"] == "DE_a"
    assert missing[0]["reason"] == "invalid_growth_phase: maybe_stressed"


def test_validate_llm_payload_rejects_bad_self_assessment():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["timepoint"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint": "24h",
                "self_assessment": "sure",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert len(missing) == 1
    assert missing[0]["analysis_id"] == "DE_a"
    assert "invalid_self_assessment" in missing[0]["reason"]


def test_validate_llm_payload_rejects_missing_requested_field():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["growth_phase"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                # growth_phase absent
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert {"analysis_id": "DE_a", "reason": "missing_field: growth_phase"} in missing


def test_validate_llm_payload_rejects_timepoint_hours_not_numeric():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["timepoint_hours"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint_hours": "24h",  # string instead of number
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert len(missing) == 1
    assert "timepoint_hours_not_numeric" in missing[0]["reason"]


def test_validate_llm_payload_accepts_null_timepoint_hours():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["timepoint_hours"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint_hours": None,
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 1
    assert missing == []


def test_validate_llm_payload_rejects_hallucinated_analysis_id():
    from multiomics_kg.extraction.timepoint.extraction_utils import validate_llm_payload

    requested_map = {"DE_a": ["timepoint"]}
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_ghost",
                "timepoint": "24h",
                "self_assessment": "high",
                "assessment_notes": "clear",
                "supporting_quotes": [],
                "source_figures": [],
            },
        ]
    }
    with pytest.raises(ValueError, match="DE_ghost"):
        validate_llm_payload(payload, requested_map)
