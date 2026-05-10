# tests/test_tcdb_diamond.py
"""Unit tests for multiomics_kg.utils.tcdb_diamond."""
from multiomics_kg.utils import tcdb_diamond  # noqa: F401
from multiomics_kg.utils.tcdb_diamond import truncate_tcid


def test_truncate_tcid_keeps_first_n_parts():
    assert truncate_tcid("1.A.11.1.5", 3) == "1.A.11"
    assert truncate_tcid("1.A.11.1.5", 4) == "1.A.11.1"
    assert truncate_tcid("1.A.11.1.5", 5) == "1.A.11.1.5"


def test_truncate_tcid_passthrough_when_already_short():
    assert truncate_tcid("1.A.11", 5) == "1.A.11"
    assert truncate_tcid("1.A", 3) == "1.A"


def test_truncate_tcid_empty_input_returns_empty():
    assert truncate_tcid("", 3) == ""


from multiomics_kg.utils.tcdb_diamond import classify_hit


def _hit(identity, qcov, scov, length=200, evalue=1e-10):
    """Build a hit dict with sensible defaults — only override what the test cares about."""
    return {
        "identity": identity, "qcov": qcov, "scov": scov,
        "length": length, "evalue": evalue,
    }


def test_classify_tier_1_high_identity_high_qcov():
    # 80%/85% identity over 80%/75% qcov — tier 1
    assert classify_hit(_hit(80.0, 85.0, 50.0)) == 1
    assert classify_hit(_hit(70.0, 70.0, 50.0)) == 1


def test_classify_tier_2_mid_identity_mid_qcov():
    # 50% identity / 65% qcov — passes tier 2 floor, fails tier 1
    assert classify_hit(_hit(50.0, 65.0, 30.0)) == 2
    assert classify_hit(_hit(40.0, 60.0, 30.0)) == 2


def test_classify_tier_3_gblast3_floor():
    # 35% identity (no tier 1/2), qcov >= 40 — tier 3
    assert classify_hit(_hit(35.0, 45.0, 30.0)) == 3
    # qcov < 40 BUT scov >= 40 — still tier 3 (gblast3 OR rule)
    assert classify_hit(_hit(35.0, 25.0, 50.0)) == 3
    # 25% identity is fine for tier 3 (no identity floor)
    assert classify_hit(_hit(25.0, 45.0, 30.0)) == 3


def test_classify_drops_hit_below_floor():
    # qcov AND scov both < 40 — fails gblast3 OR rule
    assert classify_hit(_hit(35.0, 30.0, 30.0)) is None
    # length < 50 — fails HSP-length floor
    assert classify_hit(_hit(80.0, 90.0, 90.0, length=49)) is None
    # e-value > 0.001 — fails e-value gate
    assert classify_hit(_hit(80.0, 90.0, 90.0, evalue=0.01)) is None


def test_classify_handles_boundary_values_inclusively():
    # All thresholds are >= (not >) — boundary inputs should pass
    assert classify_hit(_hit(70.0, 70.0, 0.0)) == 1
    assert classify_hit(_hit(40.0, 60.0, 0.0)) == 2
    assert classify_hit(_hit(0.0, 40.0, 0.0)) == 3
