#!/usr/bin/env python
"""Gene/protein annotation quality report per organism strain.

Queries the live Neo4j knowledge graph and reports per-strain coverage of:
  proteins, function descriptions, GO terms, KEGG, eggNOG, cyanorak roles,
  product names, protein domains, and expression data.

Supports snapshot save/diff to track annotation quality over time.

Usage:
    uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py
    uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --strain MED4
    uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --save
    uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --diff
    uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --diff --save
"""

import argparse
import csv
import io
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SNAPSHOT_FILE = Path(__file__).parent / "quality_snapshot.json"

# Neo4j connection
NEO4J_CONTAINER = "deploy"
NEO4J_USER = "neo4j"
NEO4J_PASS = "neo4j"

# Ordered list of strains (genome strains only, sorted by organism then strain)
STRAIN_ORDER = [
    "MED4", "AS9601", "MIT9301", "MIT9312", "MIT9313",
    "NATL1A", "NATL2A", "RSP50",
    "CC9311", "WH8102",
    "MIT1002", "EZ55", "HOT1A3",
]

# Column display config: (key, header, description)
METRICS = [
    ("total",          "Genes",      "Total gene nodes"),
    ("has_protein",    "Protein%",   "Linked UniProt protein (Gene_encodes_protein edge)"),
    ("has_function",   "Function%",  "Non-empty function_description (UniProt, denorm onto gene)"),
    ("has_go",         "GO%",        "Any GO term: gene.go_biological_processes/Ontology_term_ncbi/Ontology_term_cyanorak OR protein.go_molecular_functions/go_cellular_components"),
    ("has_kegg",       "KEGG%",      "gene.kegg (Cyanorak/NCBI) OR protein.kegg_ids (UniProt)"),
    ("has_eggnog",     "eggNOG%",    "gene.eggNOG (Cyanorak) OR protein.eggnog_ids (UniProt)"),
    ("has_pfam",       "Pfam%",      "protein.pfam_ids (UniProt)"),
    ("has_ec",         "EC%",        "protein.ec_numbers (UniProt enzyme classification)"),
    ("has_cyanorak",   "CyanoRole%", "gene.cyanorak_Role (Cyanorak functional role)"),
    ("has_product",    "Product%",   "gene.product or gene.product_cyanorak (any product name)"),
    ("has_domains",    "Domains%",   "gene.protein_domains (InterPro via Cyanorak)"),
    ("has_expression", "Expr%",      "Has ≥1 direct Affects_expression_of edge"),
]


def run_cypher(query: str) -> list[list]:
    """Run a Cypher query via docker exec cypher-shell. Returns list of rows."""
    cmd = [
        "docker", "exec", NEO4J_CONTAINER,
        "cypher-shell", "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        query,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Cypher query failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    rows = []
    reader = csv.reader(io.StringIO(result.stdout.strip()), skipinitialspace=True)
    # First row is the header — skip it
    try:
        next(reader)
    except StopIteration:
        return rows
    for row in reader:
        if row:
            rows.append(row)
    return rows


def fetch_quality_metrics(strain_filter: str | None = None) -> dict[str, dict]:
    """Query Neo4j and return quality metrics keyed by strain name."""

    strain_clause = ""
    if strain_filter:
        strain_clause = f"AND o.strain_name = '{strain_filter}'"

    query = f"""
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.strain_name IS NOT NULL {strain_clause}
OPTIONAL MATCH (p:Protein)-[:Gene_encodes_protein]->(g)
WITH o.organism_name AS strain, g, p
OPTIONAL MATCH ()-[expr:Affects_expression_of]->(g)
WITH strain, g, p, count(expr) AS expr_count
RETURN
  strain,
  count(g)                                                                       AS total,
  count(p)                                                                       AS has_protein,
  count(CASE WHEN g.function_description IS NOT NULL
              AND g.function_description <> '' THEN 1 END)                      AS has_function,
  count(CASE WHEN g.go_biological_processes IS NOT NULL
              OR g.Ontology_term_ncbi IS NOT NULL
              OR g.Ontology_term_cyanorak IS NOT NULL
              OR p.go_molecular_functions IS NOT NULL
              OR p.go_cellular_components IS NOT NULL THEN 1 END)               AS has_go,
  count(CASE WHEN g.kegg IS NOT NULL
              OR p.kegg_ids IS NOT NULL THEN 1 END)                             AS has_kegg,
  count(CASE WHEN g.eggNOG IS NOT NULL
              OR p.eggnog_ids IS NOT NULL THEN 1 END)                           AS has_eggnog,
  count(p.pfam_ids)                                                              AS has_pfam,
  count(p.ec_numbers)                                                            AS has_ec,
  count(g.cyanorak_Role)                                                         AS has_cyanorak,
  count(CASE WHEN g.product IS NOT NULL
              OR g.product_cyanorak IS NOT NULL THEN 1 END)                      AS has_product,
  count(g.protein_domains)                                                       AS has_domains,
  count(CASE WHEN expr_count > 0 THEN 1 END)                                    AS has_expression
ORDER BY strain
"""
    rows = run_cypher(query)
    metrics = {}
    for row in rows:
        if len(row) < 13:
            continue
        strain = row[0]
        metrics[strain] = {
            "total":          int(row[1]),
            "has_protein":    int(row[2]),
            "has_function":   int(row[3]),
            "has_go":         int(row[4]),
            "has_kegg":       int(row[5]),
            "has_eggnog":     int(row[6]),
            "has_pfam":       int(row[7]),
            "has_ec":         int(row[8]),
            "has_cyanorak":   int(row[9]),
            "has_product":    int(row[10]),
            "has_domains":    int(row[11]),
            "has_expression": int(row[12]),
        }
    return metrics


def pct(count: int, total: int) -> str:
    """Format a percentage string."""
    if total == 0:
        return "  0%"
    p = 100.0 * count / total
    return f"{p:4.0f}%"


def fmt_count(count: int, total: int) -> str:
    """Format count/total for display."""
    return f"{count:>5}/{total:<5} {pct(count, total)}"


def print_report(metrics: dict[str, dict]) -> None:
    """Print a formatted quality table."""
    # Determine column order: STRAIN_ORDER first, then any extras sorted
    known = [s for s in STRAIN_ORDER if s in metrics]
    extra = sorted(s for s in metrics if s not in STRAIN_ORDER)
    strains = known + extra

    if not strains:
        print("No data found.")
        return

    # Column widths
    strain_w = max(len(s) for s in strains) + 2
    col_w = 14  # "NNNNN/NNNNN NN%"

    # Header
    header_keys = [m[1] for m in METRICS[1:]]  # skip 'total' since it's the row label
    header = f"{'Strain':<{strain_w}} {'Genes':>6}  " + "  ".join(f"{h:>{col_w}}" for h in header_keys)
    print()
    print(header)
    print("-" * len(header))

    for strain in strains:
        m = metrics[strain]
        total = m["total"]
        row = f"{strain:<{strain_w}} {total:>6}  "
        cols = []
        for key, _hdr, _desc in METRICS[1:]:
            cols.append(f"{fmt_count(m[key], total):>{col_w}}")
        row += "  ".join(cols)
        print(row)

    print()
    print("Columns: " + ", ".join(f"{h}={d}" for _, h, d in METRICS[1:]))
    print()


def load_snapshot() -> dict | None:
    """Load the saved quality snapshot if it exists."""
    if not SNAPSHOT_FILE.exists():
        return None
    with open(SNAPSHOT_FILE) as f:
        return json.load(f)


def save_snapshot(metrics: dict[str, dict]) -> None:
    """Save metrics to the snapshot file."""
    snapshot = {
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "metrics": metrics,
    }
    with open(SNAPSHOT_FILE, "w") as f:
        json.dump(snapshot, f, indent=2)
    print(f"Snapshot saved to {SNAPSHOT_FILE}")


def print_diff(current: dict[str, dict], snapshot: dict) -> None:
    """Print differences between current metrics and saved snapshot."""
    saved_at = snapshot.get("saved_at", "unknown")
    saved = snapshot.get("metrics", {})

    all_strains = sorted(set(list(current.keys()) + list(saved.keys())))

    changes = []
    for strain in all_strains:
        if strain not in saved:
            changes.append(f"  + {strain}: NEW strain (not in snapshot)")
            continue
        if strain not in current:
            changes.append(f"  - {strain}: REMOVED (not in current graph)")
            continue

        cur = current[strain]
        prev = saved[strain]
        total_cur = cur["total"]
        total_prev = prev["total"]

        strain_changes = []
        for key, hdr, _desc in METRICS:
            c_val = cur.get(key, 0)
            p_val = prev.get(key, 0)
            if c_val != p_val:
                if key == "total":
                    strain_changes.append(f"    Genes: {p_val} → {c_val} ({c_val - p_val:+d})")
                else:
                    c_pct = 100.0 * c_val / total_cur if total_cur else 0
                    p_pct = 100.0 * p_val / total_prev if total_prev else 0
                    delta_pct = c_pct - p_pct
                    delta_n = c_val - p_val
                    strain_changes.append(
                        f"    {hdr}: {p_val} ({p_pct:.0f}%) → {c_val} ({c_pct:.0f}%)  "
                        f"[{delta_n:+d}, {delta_pct:+.1f}pp]"
                    )

        if strain_changes:
            changes.append(f"  {strain}:")
            changes.extend(strain_changes)

    print(f"\nDiff vs snapshot from {saved_at}:")
    if changes:
        print("\n".join(changes))
    else:
        print("  No changes detected.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Gene/protein annotation quality report per strain."
    )
    parser.add_argument("--strain", help="Report only this strain (e.g. MED4)")
    parser.add_argument("--save", action="store_true",
                        help="Save current metrics as the quality snapshot")
    parser.add_argument("--diff", action="store_true",
                        help="Show diff vs the saved snapshot")
    args = parser.parse_args()

    print("Querying Neo4j...", end=" ", flush=True)
    metrics = fetch_quality_metrics(args.strain)
    print(f"done ({len(metrics)} strains)")

    print_report(metrics)

    if args.diff:
        snapshot = load_snapshot()
        if snapshot is None:
            print("No snapshot found. Run with --save first.")
        else:
            print_diff(metrics, snapshot)

    if args.save:
        save_snapshot(metrics)


if __name__ == "__main__":
    main()
