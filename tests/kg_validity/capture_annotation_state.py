"""
Capture the Gene annotation_state / annotation_quality distribution from the
live Neo4j knowledge graph into a committed baseline artifact.

Why: bucket-movement claims across rebuilds ("the has_any_edge fix moved -420
no_evidence / +236 catch_all_only") are only reproducible when the pre-rebuild
distribution was captured as an artifact, not supplied from memory (the
2026-08-17 baseline was context, not a capture). Run this BEFORE a rebuild —
alongside the /omics-edge-snapshot baseline — then --compare after.

The artifact is committed (git history is the timeline, same convention as
snapshot_data.json). A distribution shift on rebuild is usually legitimate
(new data landed), so this is deliberately a capture/compare tool, not a
pass/fail test.

Usage:
    uv run python tests/kg_validity/capture_annotation_state.py            # print live distribution
    uv run python tests/kg_validity/capture_annotation_state.py --save     # overwrite the committed baseline
    uv run python tests/kg_validity/capture_annotation_state.py --compare  # diff live vs committed baseline
    uv run python tests/kg_validity/capture_annotation_state.py --neo4j-url bolt://localhost:7688 ...
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from neo4j import GraphDatabase

BASELINE_PATH = Path(__file__).parent / "annotation_state_baseline.json"

# annotation_state enum order (broadest bucket movement reads top-to-bottom)
STATE_ORDER = ["no_evidence", "catch_all_only", "informative_single", "informative_multi"]


def _dist(session, query: str) -> dict:
    """Run a (key, count) query and return an ordered {str(key): count} dict."""
    return {str(r[0]): r[1] for r in session.run(query).values()}


def capture(neo4j_url: str) -> dict:
    driver = GraphDatabase.driver(neo4j_url)
    with driver.session() as s:
        total = s.run("MATCH (g:Gene) RETURN count(g)").single()[0]
        state = _dist(s, """
            MATCH (g:Gene)
            RETURN coalesce(g.annotation_state, '<absent>') AS k, count(*) AS c
            ORDER BY k
        """)
        quality = _dist(s, """
            MATCH (g:Gene)
            RETURN coalesce(toString(g.annotation_quality), '<absent>') AS k, count(*) AS c
            ORDER BY k
        """)
        annotation_types = _dist(s, """
            MATCH (g:Gene)
            UNWIND coalesce(g.annotation_types, []) AS t
            RETURN t AS k, count(*) AS c
            ORDER BY k
        """)
        informative_types = _dist(s, """
            MATCH (g:Gene)
            UNWIND coalesce(g.informative_annotation_types, []) AS t
            RETURN t AS k, count(*) AS c
            ORDER BY k
        """)
        informative_size_hist = _dist(s, """
            MATCH (g:Gene)
            RETURN size(coalesce(g.informative_annotation_types, [])) AS k, count(*) AS c
            ORDER BY k
        """)
        per_organism_rows = s.run("""
            MATCH (g:Gene)
            RETURN coalesce(g.organism_name, '<absent>') AS org,
                   coalesce(g.annotation_state, '<absent>') AS st,
                   count(*) AS c
            ORDER BY org, st
        """).values()
    driver.close()

    per_organism: dict[str, dict[str, int]] = {}
    for org, st, c in per_organism_rows:
        per_organism.setdefault(org, {})[st] = c

    # Present annotation_state in enum order (unknown keys, if any, appended sorted)
    ordered_state = {k: state[k] for k in STATE_ORDER if k in state}
    ordered_state.update({k: v for k, v in sorted(state.items()) if k not in ordered_state})

    return {
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "neo4j_url": neo4j_url,
        "total_genes": total,
        "annotation_state": ordered_state,
        "annotation_quality": quality,
        "annotation_types": annotation_types,
        "informative_annotation_types": informative_types,
        "informative_annotation_types_size_histogram": informative_size_hist,
        "per_organism_annotation_state": per_organism,
    }


def print_capture(cap: dict) -> None:
    total = cap["total_genes"]
    print(f"total genes: {total}")
    for section in ("annotation_state", "annotation_quality",
                    "annotation_types", "informative_annotation_types",
                    "informative_annotation_types_size_histogram"):
        print(f"\n{section}:")
        for k, v in cap[section].items():
            print(f"  {k:24s} {v:>8d}  ({v / total:.1%})")


def compare(live: dict, baseline: dict) -> None:
    print(f"baseline captured_at: {baseline['captured_at']}")
    print(f"live     captured_at: {live['captured_at']}")
    dg = live["total_genes"] - baseline["total_genes"]
    print(f"total genes: {baseline['total_genes']} -> {live['total_genes']} ({dg:+d})")
    for section in ("annotation_state", "annotation_quality",
                    "annotation_types", "informative_annotation_types",
                    "informative_annotation_types_size_histogram"):
        old, new = baseline.get(section, {}), live.get(section, {})
        keys = list(dict.fromkeys(list(old) + sorted(k for k in new if k not in old)))
        lines = []
        for k in keys:
            o, n = old.get(k, 0), new.get(k, 0)
            if o != n:
                lines.append(f"  {k:24s} {o:>8d} -> {n:>8d}  ({n - o:+d})")
        print(f"\n{section}: {'unchanged' if not lines else ''}")
        for line in lines:
            print(line)
    # per-organism: only report organisms whose state row changed
    old_po = baseline.get("per_organism_annotation_state", {})
    new_po = live.get("per_organism_annotation_state", {})
    changed = [o for o in sorted(set(old_po) | set(new_po)) if old_po.get(o) != new_po.get(o)]
    print(f"\nper_organism_annotation_state: "
          f"{'unchanged' if not changed else f'{len(changed)} organism(s) changed'}")
    for org in changed:
        print(f"  {org}: {old_po.get(org, {})} -> {new_po.get(org, {})}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--neo4j-url", default="bolt://localhost:7687")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--save", action="store_true",
                      help="overwrite the committed baseline with the live distribution")
    mode.add_argument("--compare", action="store_true",
                      help="diff the live distribution against the committed baseline")
    args = ap.parse_args()

    live = capture(args.neo4j_url)
    if args.compare:
        if not BASELINE_PATH.exists():
            raise SystemExit(f"No baseline at {BASELINE_PATH} — run with --save first.")
        with open(BASELINE_PATH) as f:
            compare(live, json.load(f))
    elif args.save:
        with open(BASELINE_PATH, "w") as f:
            json.dump(live, f, indent=2)
            f.write("\n")
        print(f"Saved baseline to {BASELINE_PATH}")
        print_capture(live)
    else:
        print_capture(live)


if __name__ == "__main__":
    main()
