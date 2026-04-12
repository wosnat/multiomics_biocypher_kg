"""Tests for merge.py."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.merge import merge_one_paper
from multiomics_kg.extraction.timepoint.extraction_utils import (
    compute_paperconfig_signature,
    save_extraction_json,
)


def _write_paperconfig(tmp_path: Path, timepoint=None, timepoint_hours=None, growth_phase=None) -> Path:
    pc = tmp_path / "paperconfig.yaml"
    lines = [
        "publication:",
        "  papername: Fake",
        "  papermainpdf: data/fake.pdf",
        "  experiments:",
        "    exp1: {name: Test, omics_type: RNASEQ}",
        "  supplementary_materials:",
        "    tbl:",
        "      type: csv",
        "      filename: data/fake.csv",
        "      statistical_analyses:",
        "        - id: DE_a",
        "          experiment: exp1",
        "          name_col: gene",
        "          logfc_col: lfc",
    ]
    if timepoint is not None:
        lines.append(f"          timepoint: \"{timepoint}\"")
    if timepoint_hours is not None:
        lines.append(f"          timepoint_hours: {timepoint_hours}")
    if growth_phase is not None:
        lines.append(f"          growth_phase: {growth_phase}")
    pc.write_text("\n".join(lines) + "\n")
    return pc


def _valid_analysis(aid="DE_a", **overrides) -> dict:
    base = {
        "analysis_id": aid, "experiment_key": "exp1",
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
        "timepoint": "24h", "timepoint_hours": 24.0, "growth_phase": "exponential",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }
    base.update(overrides)
    return base


def test_merge_writes_three_fields_into_paperconfig(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )

    merge_one_paper(tmp_path, force=False)

    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["timepoint"] == "24h"
    assert a["timepoint_hours"] == 24.0
    assert a["growth_phase"] == "exponential"


def test_merge_refuses_partial_without_force(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "partial",
         "missing_analyses": [{"analysis_id": "DE_b", "reason": "not_returned"}],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )
    with pytest.raises(SystemExit) as exc:
        merge_one_paper(tmp_path, force=False)
    assert exc.value.code != 0


def test_merge_allows_partial_with_force(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "partial",
         "missing_analyses": [{"analysis_id": "DE_b", "reason": "not_returned"}],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )
    merge_one_paper(tmp_path, force=True)  # should not raise


def test_merge_refuses_on_signature_mismatch(tmp_path):
    pc = _write_paperconfig(tmp_path)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": "0" * 40},  # wrong sig
        [_valid_analysis()],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)


def test_merge_warns_on_ghost_analysis_id(tmp_path, caplog):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    # Ghost analysis: JSON has DE_ghost but paperconfig doesn't
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(aid="DE_ghost"), _valid_analysis(aid="DE_a")],
    )
    with caplog.at_level("WARNING"):
        merge_one_paper(tmp_path, force=False)
    assert any("DE_ghost" in r.message for r in caplog.records)


def test_merge_refuses_overwrite_of_differing_value(tmp_path):
    pc = _write_paperconfig(tmp_path, growth_phase="acclimated_steady_state")
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="exponential")],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)


def test_merge_allows_overwrite_with_force(tmp_path):
    pc = _write_paperconfig(tmp_path, growth_phase="acclimated_steady_state")
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="exponential")],
    )
    merge_one_paper(tmp_path, force=True)
    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "exponential"


def test_merge_rejects_malformed_growth_phase(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="banana")],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)
