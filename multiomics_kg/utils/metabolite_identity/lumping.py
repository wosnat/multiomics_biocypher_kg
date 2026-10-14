"""Stage O: the pinned TCDB substrate-lumping set (spec §3.4a). A documented R3 exception."""
from __future__ import annotations

SEED_LUMPING_IDS = ("2.A.1", "2.A.1.1", "2.A.1.2", "2.A.1.3", "2.A.3", "2.A.6", "2.A.6.2",
                    "2.A.7", "2.A.66", "2.A.66.1", "3.A.1", "3.A.1.1", "3.A.3")
LUMPING_MIN_LEVEL = 2
LUMPING_MIN_METABOLITES = 50


def threshold_lumping(kept, levels, substrate_counts):
    return {t for t in kept if levels.get(t, 0) >= LUMPING_MIN_LEVEL
            and substrate_counts.get(t, 0) >= LUMPING_MIN_METABOLITES}


def update_lumping_pin(pin, kept, levels, substrate_counts, today):
    seeded = pin is None
    if seeded:
        pin = {"frozen_on": today, "frozen_kept_ids": sorted(kept),
               "lumping_tcdb_ids": sorted(SEED_LUMPING_IDS)}
    thr = threshold_lumping(kept, levels, substrate_counts)
    frozen = set(pin["frozen_kept_ids"])
    current = set(pin["lumping_tcdb_ids"])
    judged = frozen | (current - frozen)          # frozen ids + ids appended on earlier builds
    drift = [{"tc_id": t, "pinned": "lumping" if t in current else "specific",
              "threshold": "lumping" if t in thr else "specific",
              "metabolite_count": substrate_counts.get(t, 0)}
             for t in sorted(kept & judged) if (t in current) != (t in thr)]
    added = sorted(t for t in kept - judged if t in thr)
    pin = dict(pin, lumping_tcdb_ids=sorted(current | set(added)))
    return pin, {"lumping_set_drift": drift, "lumping_set_added": added, "seeded": seeded}


def substrate_breadth(tc_id, level, pin):
    if level < LUMPING_MIN_LEVEL:
        return "not_assessed"
    return "lumping" if tc_id in set(pin.get("lumping_tcdb_ids", [])) else "specific"
