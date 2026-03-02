#!/usr/bin/env python3
"""Snapshot and compare Affects_expression_of edge counts in the Neo4j knowledge graph.

Captures per-paper edge counts before and after omics adapter changes to detect
regressions. New edges are good; lost matched edges are regressions.

Usage:
    # Capture baseline from the live graph
    uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save baseline

    # After rebuild: compare live graph to baseline
    uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare baseline

    # Compare two saved snapshots (no live Neo4j needed)
    uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before --against after

    # List saved snapshots
    uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --list
"""

import argparse
import csv
import io
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"

# Neo4j via docker (same container/credentials as gene-protein-quality)
NEO4J_CONTAINER = "deploy"
NEO4J_USER = "neo4j"
NEO4J_PASS = "neo4j"


# ---------------------------------------------------------------------------
# Neo4j helpers
# ---------------------------------------------------------------------------

def run_cypher(query: str) -> list[list[str]]:
    """Run a Cypher query via docker exec cypher-shell. Returns list of rows (header stripped)."""
    cmd = [
        "docker", "exec", NEO4J_CONTAINER,
        "cypher-shell", "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        query,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Cypher query failed:\n{result.stderr}", file=sys.stderr)
        print("Is Docker running? Try: docker compose up -d", file=sys.stderr)
        sys.exit(1)

    rows = []
    reader = csv.reader(io.StringIO(result.stdout.strip()), skipinitialspace=True)
    try:
        next(reader)  # skip header row
    except StopIteration:
        return rows
    for row in reader:
        if row:
            rows.append(row)
    return rows


def _strip_quotes(value: str) -> str:
    """cypher-shell wraps strings in double-quotes; strip them."""
    value = value.strip()
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]
    return value


# ---------------------------------------------------------------------------
# Snapshot capture
# ---------------------------------------------------------------------------

def capture_snapshot() -> dict:
    """Capture full edge metrics from the live graph."""
    snapshot: dict = {"timestamp": datetime.now(timezone.utc).isoformat()}

    # --- Total edges ---
    rows = run_cypher(
        "MATCH ()-[e:Affects_expression_of]->() RETURN count(e) AS cnt"
    )
    snapshot["total_edges"] = int(rows[0][0]) if rows else 0

    # --- Per-publication edge count (unwind publications array) ---
    rows = run_cypher("""
        MATCH ()-[e:Affects_expression_of]->()
        UNWIND e.publications AS pub
        RETURN pub, count(e) AS edge_count
        ORDER BY pub
    """)
    snapshot["per_publication"] = {
        _strip_quotes(r[0]): int(r[1]) for r in rows
    }

    # --- Per-publication by direction ---
    rows = run_cypher("""
        MATCH ()-[e:Affects_expression_of]->()
        UNWIND e.publications AS pub
        RETURN pub, e.expression_direction AS direction, count(e) AS edge_count
        ORDER BY pub, direction
    """)
    by_dir: dict = {}
    for r in rows:
        pub = _strip_quotes(r[0])
        direction = _strip_quotes(r[1]) if r[1].strip() not in ("", "null", "NULL") else "null"
        cnt = int(r[2])
        by_dir.setdefault(pub, {})[direction] = cnt
    snapshot["per_publication_by_direction"] = by_dir

    # --- By source node type ---
    rows = run_cypher("""
        MATCH (src)-[e:Affects_expression_of]->()
        RETURN labels(src)[0] AS source_type, count(e) AS edge_count
        ORDER BY source_type
    """)
    snapshot["by_source_type"] = {
        _strip_quotes(r[0]): int(r[1]) for r in rows
    }

    # --- By target gene organism (strain) ---
    rows = run_cypher("""
        MATCH ()-[e:Affects_expression_of]->(g:Gene)
        OPTIONAL MATCH (g)-[:Gene_belongs_to_organism]->(org:OrganismTaxon)
        RETURN coalesce(org.organism_name, 'unknown') AS organism, count(e) AS edge_count
        ORDER BY organism
    """)
    snapshot["by_organism"] = {
        _strip_quotes(r[0]): int(r[1]) for r in rows
    }

    # --- Per-publication gene sets (locus_tags) ---
    # Captures WHICH genes have edges, not just how many.
    # Used to detect genes that were previously matched but disappear after resolve stage.
    rows = run_cypher("""
        MATCH ()-[e:Affects_expression_of]->(g:Gene)
        UNWIND e.publications AS pub
        RETURN pub, g.locus_tag AS locus_tag
        ORDER BY pub, locus_tag
    """)
    pub_genes: dict = {}
    for r in rows:
        pub = _strip_quotes(r[0])
        locus_tag = _strip_quotes(r[1])
        pub_genes.setdefault(pub, [])
        if locus_tag not in pub_genes[pub]:
            pub_genes[pub].append(locus_tag)
    snapshot["per_publication_genes"] = pub_genes

    return snapshot


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------

def save_snapshot(snapshot: dict, name: str) -> Path:
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    path = SNAPSHOT_DIR / f"{name}.json"
    path.write_text(json.dumps(snapshot, indent=2))
    return path


def load_snapshot(name: str) -> dict:
    path = SNAPSHOT_DIR / f"{name}.json"
    if not path.exists():
        print(f"ERROR: Snapshot '{name}' not found at {path}", file=sys.stderr)
        available = [p.stem for p in sorted(SNAPSHOT_DIR.glob("*.json"))] if SNAPSHOT_DIR.exists() else []
        if available:
            print(f"Available snapshots: {', '.join(available)}", file=sys.stderr)
        sys.exit(1)
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare_snapshots(old: dict, new: dict, old_name: str, new_name: str) -> int:
    """Print a comparison report. Returns exit code (1 if regressions found)."""
    old_pubs: dict = old["per_publication"]
    new_pubs: dict = new["per_publication"]
    all_pubs = sorted(set(old_pubs) | set(new_pubs))

    regressions = []   # (pub, old_cnt, new_cnt, lost)
    improvements = []  # (pub, old_cnt, new_cnt, gained)
    new_only = []      # (pub, cnt)  — appeared in new, not in old
    missing = []       # (pub, cnt)  — in old, gone from new
    unchanged = []     # (pub, cnt)

    for pub in all_pubs:
        old_cnt = old_pubs.get(pub)
        new_cnt = new_pubs.get(pub)
        if old_cnt is None:
            new_only.append((pub, new_cnt))
        elif new_cnt is None:
            missing.append((pub, old_cnt))
        elif new_cnt < old_cnt:
            regressions.append((pub, old_cnt, new_cnt, old_cnt - new_cnt))
        elif new_cnt > old_cnt:
            improvements.append((pub, old_cnt, new_cnt, new_cnt - old_cnt))
        else:
            unchanged.append((pub, old_cnt))

    # ---- Header ----
    sep = "=" * 72
    print(f"\n{sep}")
    print("OMICS EDGE SNAPSHOT COMPARISON")
    print(f"  Before: '{old_name}'  ({old.get('timestamp', 'unknown')[:19]})")
    print(f"  After:  '{new_name}'  ({new.get('timestamp', 'unknown')[:19]})")
    print(sep)

    old_total = old.get("total_edges", 0)
    new_total = new.get("total_edges", 0)
    delta = new_total - old_total
    delta_str = f"+{delta:,}" if delta >= 0 else f"{delta:,}"
    print(f"\nTotal Affects_expression_of edges: {old_total:,} → {new_total:,}  ({delta_str})")

    # ---- By source type ----
    print("\nBy source type:")
    all_src = sorted(set(old.get("by_source_type", {})) | set(new.get("by_source_type", {})))
    for src in all_src:
        o = old.get("by_source_type", {}).get(src, 0)
        n = new.get("by_source_type", {}).get(src, 0)
        d = n - o
        d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
        print(f"  {src:<30} {o:>7,} → {n:>7,}  ({d_str})")

    # ---- By organism ----
    print("\nBy target organism:")
    all_orgs = sorted(set(old.get("by_organism", {})) | set(new.get("by_organism", {})))
    for org in all_orgs:
        o = old.get("by_organism", {}).get(org, 0)
        n = new.get("by_organism", {}).get(org, 0)
        d = n - o
        d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
        print(f"  {org:<35} {o:>7,} → {n:>7,}  ({d_str})")

    # ---- Regressions (most important) ----
    print()
    if regressions:
        bang = "!" * 72
        print(bang)
        print(f"REGRESSIONS — {len(regressions)} publication(s) LOST edges:")
        print(bang)
        for pub, old_cnt, new_cnt, lost in sorted(regressions, key=lambda x: -x[3]):
            print(f"  LOST {lost:>6,}  {pub}")
            print(f"              {old_cnt:,} → {new_cnt:,}")
    else:
        print("✓  No regressions — no publication lost edges")

    # ---- Missing publications ----
    if missing:
        print(f"\nMISSING publications — {len(missing)} disappeared entirely:")
        for pub, cnt in missing:
            print(f"  GONE          {pub}  (had {cnt:,} edges)")

    # ---- Improvements ----
    if improvements:
        print(f"\nIMPROVED — {len(improvements)} publication(s) gained edges:")
        for pub, old_cnt, new_cnt, gained in sorted(improvements, key=lambda x: -x[3]):
            print(f"  +{gained:<6,}  {pub}  ({old_cnt:,} → {new_cnt:,})")

    # ---- New publications ----
    if new_only:
        print(f"\nNEW — {len(new_only)} publication(s) added:")
        for pub, cnt in new_only:
            print(f"  NEW           {pub}  ({cnt:,} edges)")

    # ---- Unchanged ----
    print(f"\nUnchanged: {len(unchanged)} publication(s)")

    # ---- Per-direction spot-check for regressions ----
    if regressions:
        print("\nDirection breakdown for regressed publications:")
        for pub, old_cnt, new_cnt, lost in regressions:
            old_dir = old.get("per_publication_by_direction", {}).get(pub, {})
            new_dir = new.get("per_publication_by_direction", {}).get(pub, {})
            print(f"  {pub}:")
            for direction in sorted(set(old_dir) | set(new_dir)):
                o = old_dir.get(direction, 0)
                n = new_dir.get(direction, 0)
                d = n - o
                d_str = f"+{d}" if d >= 0 else str(d)
                print(f"    {direction:<6} {o:>6,} → {n:>6,}  ({d_str})")

    # ---- Gene-level diff (which locus_tags appeared / disappeared) ----
    old_genes_by_pub = old.get("per_publication_genes", {})
    new_genes_by_pub = new.get("per_publication_genes", {})
    if old_genes_by_pub or new_genes_by_pub:
        print("\n" + "-" * 72)
        print("GENE-LEVEL DIFF (per publication)")
        print("-" * 72)
        gene_problems = False
        for pub in sorted(set(old_genes_by_pub) | set(new_genes_by_pub)):
            old_set = set(old_genes_by_pub.get(pub, []))
            new_set = set(new_genes_by_pub.get(pub, []))
            lost_genes = sorted(old_set - new_set)
            gained_genes = sorted(new_set - old_set)
            if not lost_genes and not gained_genes:
                continue
            print(f"\n  {pub}")
            print(f"    Genes before: {len(old_set):,}  after: {len(new_set):,}  "
                  f"lost: {len(lost_genes):,}  gained: {len(gained_genes):,}")
            if lost_genes:
                gene_problems = True
                sample = lost_genes[:10]
                more = len(lost_genes) - len(sample)
                print(f"    LOST genes (sample): {', '.join(sample)}"
                      + (f"  ... +{more} more" if more else ""))
            if gained_genes:
                sample = gained_genes[:10]
                more = len(gained_genes) - len(sample)
                print(f"    Gained genes (sample): {', '.join(sample)}"
                      + (f"  ... +{more} more" if more else ""))
        if not gene_problems:
            print("  ✓  No genes lost across any publication")
    else:
        print("\n(Gene-level diff not available — re-capture snapshots to enable)")

    has_problems = bool(regressions or missing)
    if has_problems:
        print("\nEXIT 1: regressions or missing publications detected.")
    else:
        print("\nEXIT 0: graph looks healthy.")
    return 1 if has_problems else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Snapshot and compare Affects_expression_of edge counts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--save", metavar="NAME",
                        help="Capture snapshot from live Neo4j and save as NAME")
    parser.add_argument("--compare", metavar="OLD_NAME",
                        help="Compare OLD_NAME snapshot to live Neo4j (or use --against for two saved snapshots)")
    parser.add_argument("--against", metavar="NEW_NAME",
                        help="Compare to this saved snapshot instead of querying live Neo4j")
    parser.add_argument("--list", action="store_true",
                        help="List available saved snapshots")

    args = parser.parse_args()

    if args.list:
        if not SNAPSHOT_DIR.exists() or not list(SNAPSHOT_DIR.glob("*.json")):
            print("No snapshots found.")
            return
        snapshots = sorted(SNAPSHOT_DIR.glob("*.json"))
        print(f"{'Name':<30} {'Timestamp':<22} {'Total':>8}  Pubs")
        print("-" * 72)
        for p in snapshots:
            data = json.loads(p.read_text())
            ts = data.get("timestamp", "unknown")[:19]
            total = data.get("total_edges", "?")
            npubs = len(data.get("per_publication", {}))
            print(f"{p.stem:<30} {ts:<22} {total:>8,}  {npubs}")
        return

    if args.save:
        print(f"Capturing snapshot from live Neo4j graph...")
        snapshot = capture_snapshot()
        path = save_snapshot(snapshot, args.save)
        print(f"Saved '{args.save}' → {path}")
        print(f"  Total edges : {snapshot['total_edges']:,}")
        print(f"  Publications: {len(snapshot['per_publication'])}")
        by_src = snapshot.get("by_source_type", {})
        for src, cnt in sorted(by_src.items()):
            print(f"  {src}: {cnt:,}")
        return

    if args.compare:
        old = load_snapshot(args.compare)
        if args.against:
            new = load_snapshot(args.against)
            new_name = args.against
        else:
            print(f"Capturing live snapshot for comparison...")
            new = capture_snapshot()
            new_name = "live"
        exit_code = compare_snapshots(old, new, args.compare, new_name)
        sys.exit(exit_code)

    parser.print_help()


if __name__ == "__main__":
    main()
