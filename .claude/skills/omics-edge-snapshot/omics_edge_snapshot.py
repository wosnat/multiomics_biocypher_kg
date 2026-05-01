#!/usr/bin/env python3
"""Snapshot and compare expression edge counts in the Neo4j knowledge graph.

Captures per-paper edge counts before and after omics adapter changes to detect
regressions. Counts Changes_expression_of edges (from Experiment nodes to Gene
nodes). New edges are good; lost matched edges are regressions.

Backward compatible: when comparing against old snapshots that used
condition_edges/coculture_edges fields, the tool sums them for comparison.

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
        "MATCH ()-[r:Changes_expression_of]->() RETURN count(r) AS total"
    )
    snapshot["total_edges"] = int(rows[0][0]) if rows else 0

    # --- Per-publication edge count (via Experiment -> Publication) ---
    rows = run_cypher("""
        MATCH (pub:Publication)-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(g:Gene)
        RETURN pub.doi AS publication, count(r) AS edges
        ORDER BY publication
    """)
    snapshot["per_publication"] = {
        _strip_quotes(r[0]): int(r[1]) for r in rows
    }

    # --- Per-publication by direction ---
    rows = run_cypher("""
        MATCH (pub:Publication)-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(g:Gene)
        RETURN pub.doi AS publication, r.expression_direction AS direction, count(r) AS edges
        ORDER BY publication, direction
    """)
    by_dir: dict = {}
    for r in rows:
        pub = _strip_quotes(r[0])
        direction = _strip_quotes(r[1]) if r[1].strip() not in ("", "null", "NULL") else "null"
        cnt = int(r[2])
        by_dir.setdefault(pub, {})[direction] = cnt
    snapshot["per_publication_by_direction"] = by_dir

    # --- DM edge counts per publication + type ---
    dm_edge_types = [
        "Derived_metric_flags_gene",
        "Derived_metric_classifies_gene",
        "Derived_metric_quantifies_gene",
    ]
    snapshot["dm_edges_per_publication"] = {}
    for edge_type in dm_edge_types:
        rows = run_cypher(f"""
            MATCH (pub:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
              -[r:{edge_type}]->(g:Gene)
            RETURN pub.doi AS publication, count(r) AS edges
            ORDER BY publication
        """)
        snapshot["dm_edges_per_publication"][edge_type] = {
            _strip_quotes(r[0]): int(r[1]) for r in rows
        }
    # Total DM edge counts (derived, for quick-look reports)
    for edge_type in dm_edge_types:
        rows = run_cypher(f"MATCH ()-[r:{edge_type}]->() RETURN count(r) AS total")
        snapshot[f"total_{edge_type}"] = int(rows[0][0]) if rows else 0

    # --- By target gene organism (strain) ---
    rows = run_cypher("""
        MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
        RETURN g.organism_name AS organism, count(r) AS edges
        ORDER BY organism
    """)
    snapshot["by_organism"] = {
        _strip_quotes(r[0]): int(r[1]) for r in rows
    }

    # --- Per-publication gene sets (locus_tags) ---
    # Captures WHICH genes have edges, not just how many.
    # Used to detect genes that were previously matched but disappear after resolve stage.
    rows = run_cypher("""
        MATCH (pub:Publication)-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(g:Gene)
        RETURN pub.doi AS publication, collect(DISTINCT g.locus_tag) AS genes
        ORDER BY publication
    """)
    pub_genes: dict = {}
    for r in rows:
        pub = _strip_quotes(r[0])
        # cypher-shell returns lists as bracket-delimited strings; parse them
        genes_raw = r[1].strip()
        if genes_raw.startswith("[") and genes_raw.endswith("]"):
            genes_raw = genes_raw[1:-1]
        if genes_raw:
            genes = [_strip_quotes(g.strip()) for g in genes_raw.split(",") if g.strip()]
        else:
            genes = []
        pub_genes[pub] = sorted(genes)
    snapshot["per_publication_genes"] = pub_genes

    # --- Metabolism layer (Phase 1.2 / 1.2.1) ---
    snapshot["metabolism"] = _capture_metabolism()

    return snapshot


# Metabolism layer: node + edge counts captured by capture_snapshot().
# Reaction/Metabolite are nodes; the rest are edge labels.
METABOLISM_NODE_COUNTS = ["Reaction", "Metabolite"]
METABOLISM_EDGE_COUNTS = [
    "Gene_catalyzes_reaction",
    "Reaction_has_metabolite",
    "Reaction_in_kegg_pathway",
    "Organism_has_metabolite",
    "Metabolite_in_pathway",
]
KEGG_TERM_LEVEL_KINDS = ["category", "subcategory", "pathway", "ko"]


def _capture_metabolism() -> dict:
    """Capture metabolism-layer node + edge counts.

    Phase 1.2 introduced Reaction/Metabolite nodes + 3 metabolism edges.
    Phase 1.2.1 added the Metabolite_in_pathway edge and changed step-6 pruning
    so a few additional reaction-only pathway nodes appear.
    """
    out: dict = {"nodes": {}, "edges": {}, "kegg_term_by_level_kind": {}}
    for label in METABOLISM_NODE_COUNTS:
        rows = run_cypher(f"MATCH (n:{label}) RETURN count(n) AS n")
        out["nodes"][label] = int(rows[0][0]) if rows else 0
    for edge_label in METABOLISM_EDGE_COUNTS:
        rows = run_cypher(f"MATCH ()-[r:{edge_label}]->() RETURN count(r) AS n")
        out["edges"][edge_label] = int(rows[0][0]) if rows else 0
    rows = run_cypher(
        "MATCH (k:KeggTerm) RETURN k.level_kind AS level_kind, count(k) AS n "
        "ORDER BY level_kind"
    )
    for r in rows:
        level_kind = _strip_quotes(r[0]) if r[0].strip() not in ("", "null", "NULL") else "null"
        out["kegg_term_by_level_kind"][level_kind] = int(r[1])
    return out


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------

def save_snapshot(snapshot: dict, name: str) -> Path:
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    path = SNAPSHOT_DIR / f"{name}.json"
    path.write_text(json.dumps(snapshot, indent=2, sort_keys=True))
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

def _get_total_edges(snap: dict) -> int:
    """Extract total edge count, handling old-format snapshots.

    Old snapshots stored condition_edges + coculture_edges separately.
    New snapshots store a single total_edges field.
    """
    if "total_edges" in snap:
        return snap["total_edges"]
    # Old format: sum the two sub-counts
    return snap.get("condition_edges", 0) + snap.get("coculture_edges", 0)


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

    old_total = _get_total_edges(old)
    new_total = _get_total_edges(new)
    delta = new_total - old_total
    delta_str = f"+{delta:,}" if delta >= 0 else f"{delta:,}"
    print(f"\nTotal expression edges: {old_total:,} -> {new_total:,}  ({delta_str})")

    # Legacy per-type breakdown (only shown when comparing old-format snapshots)
    if "condition_edges" in old or "condition_edges" in new:
        for key, label in [("condition_edges", "Condition_changes_expression_of"),
                           ("coculture_edges", "Coculture_changes_expression_of")]:
            o = old.get(key, 0)
            n = new.get(key, 0)
            if o or n:
                d = n - o
                d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
                print(f"  (legacy) {label:<40} {o:>7,} -> {n:>7,}  ({d_str})")

    # ---- By organism ----
    print("\nBy target organism:")
    all_orgs = sorted(set(old.get("by_organism", {})) | set(new.get("by_organism", {})))
    for org in all_orgs:
        o = old.get("by_organism", {}).get(org, 0)
        n = new.get("by_organism", {}).get(org, 0)
        d = n - o
        d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
        print(f"  {org:<35} {o:>7,} → {n:>7,}  ({d_str})")

    # ---- DerivedMetric edge totals ----
    dm_edge_types = [
        "Derived_metric_flags_gene",
        "Derived_metric_classifies_gene",
        "Derived_metric_quantifies_gene",
    ]
    has_any_dm = any(
        old.get(f"total_{et}", 0) > 0 or new.get(f"total_{et}", 0) > 0
        for et in dm_edge_types
    )
    if has_any_dm:
        print("\nDerivedMetric edge totals:")
        for edge_type in dm_edge_types:
            total_key = f"total_{edge_type}"
            old_total = old.get(total_key, 0)
            new_total = new.get(total_key, 0)
            delta = new_total - old_total
            d_str = f"+{delta:,}" if delta >= 0 else f"{delta:,}"
            print(f"  {edge_type:<40} {old_total:>7,} -> {new_total:>7,}  ({d_str})")
            if delta < 0:
                regressions.append((f"<{edge_type} total>", old_total, new_total, -delta))

    # ---- Metabolism layer ----
    old_meta = old.get("metabolism")
    new_meta = new.get("metabolism")
    if old_meta or new_meta:
        print("\nMetabolism nodes / edges:")
        if not old_meta:
            print("  (no metabolism baseline in old snapshot)")
        if not new_meta:
            print("  (no metabolism data in new snapshot)")
        if old_meta and new_meta:
            for label in METABOLISM_NODE_COUNTS:
                o = old_meta.get("nodes", {}).get(label, 0)
                n = new_meta.get("nodes", {}).get(label, 0)
                d = n - o
                d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
                print(f"  {label + ' nodes':<35} {o:>7,} → {n:>7,}  ({d_str})")
                if d < 0:
                    regressions.append((f"<{label} nodes>", o, n, -d))
            for edge_label in METABOLISM_EDGE_COUNTS:
                o = old_meta.get("edges", {}).get(edge_label, 0)
                n = new_meta.get("edges", {}).get(edge_label, 0)
                d = n - o
                d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
                print(f"  {edge_label:<35} {o:>7,} → {n:>7,}  ({d_str})")
                # Negative delta is a regression UNLESS the edge is brand-new in 1.2.1
                # (Metabolite_in_pathway). For new edge types, baseline=0, post>0 — fine.
                if d < 0:
                    regressions.append((f"<{edge_label}>", o, n, -d))
            old_levels = old_meta.get("kegg_term_by_level_kind", {})
            new_levels = new_meta.get("kegg_term_by_level_kind", {})
            if old_levels or new_levels:
                print("  KeggTerm by level_kind:")
                for level_kind in KEGG_TERM_LEVEL_KINDS:
                    o = old_levels.get(level_kind, 0)
                    n = new_levels.get(level_kind, 0)
                    d = n - o
                    d_str = f"+{d:,}" if d >= 0 else f"{d:,}"
                    print(f"    {level_kind:<33} {o:>7,} → {n:>7,}  ({d_str})")
                    # Negative delta on KeggTerm count is a regression (data loss).
                    if d < 0:
                        regressions.append((f"<KeggTerm {level_kind}>", o, n, -d))

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
        description="Snapshot and compare Changes_expression_of edge counts (Experiment -> Gene)",
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
            total = _get_total_edges(data) if ("total_edges" in data or "condition_edges" in data) else "?"
            npubs = len(data.get("per_publication", {}))
            print(f"{p.stem:<30} {ts:<22} {total:>8,}  {npubs}")
        return

    if args.save:
        print(f"Capturing snapshot from live Neo4j graph...")
        snapshot = capture_snapshot()
        path = save_snapshot(snapshot, args.save)
        print(f"Saved '{args.save}' -> {path}")
        print(f"  Total edges             : {snapshot['total_edges']:,}")
        print(f"  Publications            : {len(snapshot['per_publication'])}")
        by_org = snapshot.get("by_organism", {})
        for org, cnt in sorted(by_org.items()):
            print(f"  {org}: {cnt:,}")
        meta = snapshot.get("metabolism", {})
        if meta:
            print("  Metabolism:")
            for label, n in meta.get("nodes", {}).items():
                print(f"    {label} nodes: {n:,}")
            for edge_label, n in meta.get("edges", {}).items():
                print(f"    {edge_label}: {n:,}")
            for level_kind, n in meta.get("kegg_term_by_level_kind", {}).items():
                print(f"    KeggTerm/{level_kind}: {n:,}")
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
