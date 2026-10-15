from multiomics_kg.utils.metabolite_identity.lumping import (
    SEED_LUMPING_IDS, substrate_breadth, update_lumping_pin)

KEPT = set(SEED_LUMPING_IDS) | {"2.A", "9.B.1"}
LEVELS = {i: i.count(".") for i in KEPT}
COUNTS = {i: 60 for i in SEED_LUMPING_IDS} | {"3.A.3": 47, "9.B.1": 3, "2.A": 900}


def test_first_run_seeds_literal_13_not_threshold():
    pin, rep = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    assert pin["lumping_tcdb_ids"] == sorted(SEED_LUMPING_IDS) and rep["seeded"]
    assert "3.A.3" in pin["lumping_tcdb_ids"]                 # 47 molecules, still pinned
    assert pin["frozen_kept_ids"] == sorted(KEPT)


def test_frozen_id_drift_reported_not_applied():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    pin2, rep = update_lumping_pin(pin, KEPT, LEVELS, dict(COUNTS, **{"9.B.1": 99}), "2026-11-01")
    assert "9.B.1" not in pin2["lumping_tcdb_ids"]
    assert [d["tc_id"] for d in rep["lumping_set_drift"]] == ["3.A.3", "9.B.1"]
    assert pin2["frozen_on"] == "2026-10-15"


def test_newly_kept_id_auto_classified_and_appended():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    kept = KEPT | {"2.A.99"}
    pin2, rep = update_lumping_pin(pin, kept, dict(LEVELS, **{"2.A.99": 2}), dict(COUNTS, **{"2.A.99": 300}), "x")
    assert "2.A.99" in pin2["lumping_tcdb_ids"] and rep["lumping_set_added"] == ["2.A.99"]
    assert "2.A.99" not in pin2["frozen_kept_ids"]


def test_substrate_breadth():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "d")
    assert substrate_breadth("3.A.3", 2, pin) == "lumping"
    assert substrate_breadth("9.B.1", 2, pin) == "specific"
    assert substrate_breadth("2.A", 1, pin) == "not_assessed"
