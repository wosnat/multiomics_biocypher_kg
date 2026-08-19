#!/usr/bin/env python3
"""Derive per-strain genomic-island gene membership from Hackl 2023 mmc2.

``mmc2.xlsx`` (Table S2) holds 5,598 predicted genomic islands across 623
genomes as ``(genome_id, contig_id, start, end)`` intervals. Eleven of those
genomes are deployed strains whose assemblies match Hackl's coordinate frame
(verified by assembly-size deltas <= 236 bp against mmc1 — see
``plans/geo_paperconfig_updates.md`` item 2.4). MIT1327 is deliberately
EXCLUDED: Hackl used a different assembly (2.58 Mb / 23 contigs vs our
GCF_001632125.1), so its islands need a coordinate remap first (Tier 3).

Membership predicate is **containment** (``gene.start >= island.start AND
gene.end <= island.end``) — a gene straddling an island boundary is not an
island gene. Genes are read from each strain's committed
``gene_annotations_merged.json``; only the strain's primary contig (the one
carrying the most genes) participates, since every Hackl interval for these
strains sits on the single main contig (``<genome>_1``).

Output (one per strain, committed):
    data/Prochlorococcus/papers_and_supp/Hackl 2023/islands/<STRAIN>_island_membership.csv
    columns: locus_tag, island_id, product
    island_id = ``<contig_id>_<start>_<end>`` (contig_id starts with the
    Hackl genome_id, so identity + coordinates are both encoded — GeneCluster
    has no dedicated coordinate fields).

Plus ``islands/island_membership_report.txt`` with per-strain island/gene
counts and the MED4 hypothetical-fraction spot check (island genes are ~50%
hypothetical vs ~24% genome-wide — the flexible-genome signature that
validates the coordinate transfer biologically).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Hackl 2023")
SRC_XLSX = PAPER_DIR / "mmc2.xlsx"
OUT_DIR = PAPER_DIR / "islands"
REPORT = OUT_DIR / "island_membership_report.txt"

# Hackl genome_id -> our strain dir name (all Prochlorococcus).
# MIT1327 excluded: different assembly, coordinate frame does not transfer.
STRAINS: dict[str, str] = {
    "MED4": "MED4",
    "AS9601": "AS9601",
    "MIT9301": "MIT9301",
    "MIT9312": "MIT9312",
    "MIT9313": "MIT9313",
    "MIT9303": "MIT9303",
    "NATL1A": "NATL1A",
    "NATL2A": "NATL2A",
    "SS120": "SS120",
    "MIT1314": "MIT1314",
    "RS50": "RSP50",
}

CACHE = Path("cache/data/Prochlorococcus/genomes")


def main() -> None:
    islands = pd.read_excel(SRC_XLSX, sheet_name=0)
    OUT_DIR.mkdir(exist_ok=True)
    report_lines: list[str] = []

    for genome_id, strain in STRAINS.items():
        sub = islands[islands["genome_id"] == genome_id]
        # every wired strain's intervals sit on the single main contig; the
        # join below compares against primary-contig genes only, so fail loud
        # if a future input breaks that assumption
        assert (sub["contig_id"] == f"{genome_id}_1").all(), (
            f"{genome_id}: islands on unexpected contigs {sorted(set(sub['contig_id']))}"
        )
        anno_path = CACHE / strain / "gene_annotations_merged.json"
        genes = json.load(anno_path.open())

        # Restrict to the primary contig (most genes). Every Hackl interval
        # for these strains is on the single main contig.
        contig_counts = Counter(g.get("contig") for g in genes.values())
        primary_contig, _ = contig_counts.most_common(1)[0]
        if len(contig_counts) > 1:
            others = {c: n for c, n in contig_counts.items() if c != primary_contig}
            report_lines.append(
                f"  [note] {strain}: non-primary contigs skipped: {others}"
            )

        rows = []
        for island in sub.itertuples():
            island_id = f"{island.contig_id}_{island.start}_{island.end}"
            for lt, g in genes.items():
                if g.get("contig") != primary_contig:
                    continue
                gs, ge = g.get("start"), g.get("end")
                if gs is None or ge is None:
                    continue
                if gs >= island.start and ge <= island.end:
                    rows.append(
                        {
                            "locus_tag": lt,
                            "island_id": island_id,
                            "product": g.get("product") or "",
                        }
                    )

        df = pd.DataFrame(rows, columns=["locus_tag", "island_id", "product"])
        out = OUT_DIR / f"{strain}_island_membership.csv"
        df.to_csv(out, index=False)

        n_islands_hit = df["island_id"].nunique()
        line = (
            f"{strain:8s} islands={len(sub):3d} (with genes: {n_islands_hit:3d}) "
            f"island_genes={len(df):4d} of {len(genes):4d} "
            f"({len(df) / len(genes):5.1%} of genome)"
        )
        print(line)
        report_lines.append(line)

        if strain == "MED4":
            hyp = df["product"].str.contains("hypothetical", case=False).mean()
            genome_hyp = sum(
                1
                for g in genes.values()
                if "hypothetical" in (g.get("product") or "").lower()
            ) / len(genes)
            check = (
                f"  MED4 spot check: hypothetical fraction in islands "
                f"{hyp:.1%} vs genome-wide {genome_hyp:.1%}"
            )
            print(check)
            report_lines.append(check)

    REPORT.write_text("\n".join(report_lines) + "\n")
    print(f"report -> {REPORT}")


if __name__ == "__main__":
    main()
