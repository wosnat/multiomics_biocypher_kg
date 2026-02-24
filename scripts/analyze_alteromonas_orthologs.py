"""
Analyze ortholog coverage for Alteromonas strains using shared protein_id on Gene nodes.

Genes across strains that share the same protein_id are defined as orthologs.
Reports per-strain coverage and overall ortholog group statistics.
"""

from collections import defaultdict
from neo4j import GraphDatabase

NEO4J_URL = "bolt://localhost:7687"
AUTH = ("neo4j", "")

STRAINS = {
    "insdc.gcf:GCF_901457835.2": "MIT1002",
    "insdc.gcf:GCF_901457815.2": "EZ55",
    "insdc.gcf:GCF_001578515.1": "HOT1A3",
}


def main():
    driver = GraphDatabase.driver(NEO4J_URL, auth=AUTH)
    with driver.session() as s:
        # Fetch all Alteromonas genes with protein_id
        result = s.run(
            """
            MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
            WHERE o.id IN $gcfs
            RETURN g.id AS gene_id, g.protein_id AS protein_id, o.id AS org_id
            """,
            gcfs=list(STRAINS.keys()),
        )

        # protein_id -> set of org_ids that have a gene with this protein_id
        protein_to_orgs: dict[str, set] = defaultdict(set)
        # org_id -> list of (gene_id, protein_id)
        strain_genes: dict[str, list] = defaultdict(list)

        total_rows = 0
        for rec in result:
            total_rows += 1
            org_id = rec["org_id"]
            gene_id = rec["gene_id"]
            pid = rec["protein_id"]
            strain_genes[org_id].append((gene_id, pid))
            if pid:
                protein_to_orgs[pid].add(org_id)

        # Classify proteins: cross-strain (ortholog) vs single-strain
        cross_strain_proteins = {p for p, orgs in protein_to_orgs.items() if len(orgs) >= 2}

        print(f"\n{'='*60}")
        print("Alteromonas ortholog coverage via shared protein_id")
        print(f"{'='*60}")
        print(f"\nTotal gene records fetched: {total_rows}")
        print(f"Unique protein_ids: {len(protein_to_orgs)}")
        print(f"  Cross-strain protein_ids (>=2 strains): {len(cross_strain_proteins)}")
        print(f"  Single-strain protein_ids:              {len(protein_to_orgs) - len(cross_strain_proteins)}")

        print(f"\n{'Strain':<10} {'Total':>7} {'Has protein_id':>15} {'protein_id %':>13} {'Ortholog genes':>15} {'Ortholog %':>11}")
        print("-" * 65)

        for org_id, strain in STRAINS.items():
            genes = strain_genes[org_id]
            total = len(genes)
            with_pid = sum(1 for _, pid in genes if pid)
            ortholog = sum(1 for _, pid in genes if pid and pid in cross_strain_proteins)

            pid_pct = with_pid / total * 100 if total else 0
            orth_pct = ortholog / total * 100 if total else 0

            print(
                f"{strain:<10} {total:>7} {with_pid:>15} {pid_pct:>12.1f}% {ortholog:>15} {orth_pct:>10.1f}%"
            )

        # Ortholog group size distribution
        print(f"\n{'='*60}")
        print("Ortholog group size distribution (cross-strain proteins)")
        print(f"{'='*60}")
        size_dist: dict[int, int] = defaultdict(int)
        for p in cross_strain_proteins:
            n_strains = len(protein_to_orgs[p])
            size_dist[n_strains] += 1
        for n in sorted(size_dist):
            label = "all 3 strains" if n == 3 else f"{n} strains"
            print(f"  {label}: {size_dist[n]} protein_ids")

    driver.close()


if __name__ == "__main__":
    main()
