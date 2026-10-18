"""Control timing fields (plans/per_arm_sampling_and_shared_model.md).

Validator rules 1-5 (scripts/validate_paperconfig.py::validate_control_timing) and
the omics_adapter Experiment-node emission of control_timing / design_notes.
"""
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from validate_paperconfig import VALID_CONTROL_TIMINGS, validate_control_timing  # noqa: E402

from multiomics_kg.adapters.omics_adapter import OMICSAdapter  # noqa: E402


def _exp(**extra):
    base = {"organism": "Prochlorococcus MED4", "omics_type": "RNASEQ"}
    base.update(extra)
    return base


def _pub(default="matched", **experiments):
    pub = {"experiments": experiments}
    if default is not None:
        pub["control_timing"] = default
    return pub


def _supp(*pairs):
    """pairs of (experiment_key, timepoint)."""
    return {
        "t1": {
            "type": "csv",
            "statistical_analyses": [
                {"id": f"a{i}", "experiment": k, "timepoint": tp}
                for i, (k, tp) in enumerate(pairs)
            ],
        }
    }


def test_vocab_loaded_from_controlled_vocabularies():
    assert VALID_CONTROL_TIMINGS == {"matched", "baseline", "mixed"}


# --- rule 1 ------------------------------------------------------------------

def test_rule1_paper_default_resolves():
    errors, _ = validate_control_timing(_pub(e1=_exp()), _supp(("e1", "1h")))
    assert errors == []


def test_rule1_missing_everywhere_is_error():
    errors, _ = validate_control_timing(_pub(None, e1=_exp()), _supp(("e1", "1h")))
    assert any("missing control_timing" in e for e in errors)


def test_rule1_override_without_default_ok():
    errors, _ = validate_control_timing(
        _pub(None, e1=_exp(control_timing="matched")), _supp(("e1", "1h"))
    )
    assert errors == []


def test_rule1_bad_value_is_error():
    errors, _ = validate_control_timing(
        _pub(e1=_exp(control_timing="same_time")), _supp(("e1", "1h"))
    )
    assert any("control_timing" in e and "same_time" in e for e in errors)


def test_rule1_bad_paper_default_is_error():
    errors, _ = validate_control_timing(_pub("nope", e1=_exp()), _supp(("e1", "1h")))
    assert any("nope" in e for e in errors)


def test_rule1_non_de_experiment_needs_nothing():
    errors, warnings = validate_control_timing(
        _pub(None, e1=_exp(control_timing="matched"), metab=_exp()),
        _supp(("e1", "1h")),
    )
    assert errors == [] and warnings == []


# --- rule 2 ------------------------------------------------------------------

@pytest.mark.parametrize("timing", ["baseline", "mixed"])
def test_rule2_baseline_mixed_need_design_notes(timing):
    errors, _ = validate_control_timing(
        _pub(e1=_exp(control_timing=timing)), _supp(("e1", "1h"), ("e1", "2h"))
    )
    assert any("requires design_notes" in e for e in errors)


def test_rule2_baseline_from_paper_default_needs_design_notes():
    errors, _ = validate_control_timing(_pub("baseline", e1=_exp()), _supp(("e1", "1h")))
    assert any("requires design_notes" in e for e in errors)


def test_rule2_satisfied():
    errors, _ = validate_control_timing(
        _pub("baseline", e1=_exp(design_notes="1h vs time 0 (exponential).")),
        _supp(("e1", "1h")),
    )
    assert errors == []


# --- rule 3 ------------------------------------------------------------------

def test_rule3_mixed_single_timepoint_is_error():
    errors, _ = validate_control_timing(
        _pub(e1=_exp(control_timing="mixed", design_notes="x")), _supp(("e1", "day 18"))
    )
    assert any(">= 2 distinct timepoints" in e for e in errors)


def test_rule3_mixed_two_timepoints_ok():
    errors, _ = validate_control_timing(
        _pub(e1=_exp(control_timing="mixed", design_notes="x")),
        _supp(("e1", "day 11"), ("e1", "day 18")),
    )
    assert errors == []


# --- rule 4 ------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["day 18 | day 14", "the control's day"])
def test_rule4_design_notes_forbidden_chars(bad):
    errors, _ = validate_control_timing(
        _pub(e1=_exp(design_notes=bad)), _supp(("e1", "1h"))
    )
    assert any("design_notes must not contain" in e for e in errors)


# --- rule 5 ------------------------------------------------------------------

def test_rule5_control_timing_on_non_de_experiment_warns():
    errors, warnings = validate_control_timing(
        _pub(e1=_exp(), metab=_exp(control_timing="matched")), _supp(("e1", "1h"))
    )
    assert errors == []
    assert any("control_timing" in w and "metab" in w for w in warnings)


# --- adapter -----------------------------------------------------------------

@pytest.fixture
def adapter_nodes():
    tmp = tempfile.mkdtemp()
    data_file = os.path.join(tmp, "de.csv")
    pd.DataFrame({
        "Synonym": ["PMM0001", "PMM0002"],
        "log2_fold_change": [1.5, -2.0],
        "adjusted_p_value": [0.01, 0.02],
    }).to_csv(data_file, index=False)

    def _exp_block(**extra):
        block = {
            "name": "Test", "organism": "Prochlorococcus MED4",
            "treatment_condition": "T", "control_condition": "C",
            "omics_type": "RNASEQ", "test_type": "DESeq2", "treatment_type": ["nitrogen"],
        }
        block.update(extra)
        return block

    def _analysis(aid, exp_key):
        return {"id": aid, "experiment": exp_key, "name_col": "Synonym",
                "logfc_col": "log2_fold_change", "adjusted_p_value_col": "adjusted_p_value"}

    config = {"publication": {
        "papername": "Test CT",
        "doi": "10.1234/test.ct",
        "control_timing": "matched",
        "experiments": {
            "plain": _exp_block(statistical_model_id="legacy_key"),  # ignored: concept dropped
            "base": _exp_block(control_timing="baseline",
                               design_notes="day 18 vs day 7 (exponential); it's | odd"),
            "no_de": _exp_block(),
        },
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": data_file,
            "statistical_analyses": [_analysis("a1", "plain"), _analysis("a2", "base"),
                                     ],
        }},
    }}
    config_file = os.path.join(tmp, "paperconfig.yaml")
    with open(config_file, "w") as f:
        yaml.dump(config, f)
    adapter = OMICSAdapter(config_file=config_file)
    adapter.extracted_data = {"title": "T", "authors": [], "journal": "", "year": 2024,
                              "abstract": "", "description": ""}
    return {n[0]: n[2] for n in adapter.get_nodes() if n[1] == "experiment"}


def _by_key(nodes, key):
    return next(v for k, v in nodes.items() if k.endswith("_" + key))


def test_adapter_paper_default_applied(adapter_nodes):
    props = _by_key(adapter_nodes, "plain")
    assert props["control_timing"] == "matched"
    assert "design_notes" not in props


def test_adapter_override(adapter_nodes):
    assert _by_key(adapter_nodes, "base")["control_timing"] == "baseline"


def test_adapter_never_emits_statistical_model_id(adapter_nodes):
    assert all("statistical_model_id" not in p for p in adapter_nodes.values())


def test_adapter_design_notes_cleaned(adapter_nodes):
    notes = _by_key(adapter_nodes, "base")["design_notes"]
    assert "|" not in notes and "'" not in notes
    assert notes.startswith("day 18 vs day 7")


def test_adapter_no_control_timing_without_de(adapter_nodes):
    assert "control_timing" not in _by_key(adapter_nodes, "no_de")
