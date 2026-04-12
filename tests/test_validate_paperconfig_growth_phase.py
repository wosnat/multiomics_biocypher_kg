"""Tests for the growth_phase / timepoint additions to validate_paperconfig.py."""
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).parent.parent / ".claude/skills/paperconfig"
sys.path.insert(0, str(SKILL_DIR))

from validate_paperconfig import (
    VALID_GROWTH_PHASES,
    is_valid_growth_phase,
)


def test_is_valid_growth_phase_accepts_canonical():
    for v in VALID_GROWTH_PHASES:
        assert is_valid_growth_phase(v), f"rejected canonical value: {v}"


def test_is_valid_growth_phase_accepts_other_slug():
    assert is_valid_growth_phase("other:heat_acclimated")
    assert is_valid_growth_phase("other:x")  # minimal valid


def test_is_valid_growth_phase_rejects_other_empty():
    assert not is_valid_growth_phase("other:")
    assert not is_valid_growth_phase("other")  # missing colon


def test_is_valid_growth_phase_rejects_bare_nonenum():
    assert not is_valid_growth_phase("maybe_stressed")
    assert not is_valid_growth_phase("")


def test_valid_growth_phases_includes_all_11():
    expected = {
        "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
        "infected", "recovery", "diel", "darkness", "death", "acute_stress",
        "unknown",
    }
    assert VALID_GROWTH_PHASES == expected
