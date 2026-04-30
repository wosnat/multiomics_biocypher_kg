"""Step 6 — prune-then-enrich KEGG metabolism cache.

Walks every strain's gene_annotations_merged.json to identify gene-reachable
KEGG R-numbers, expands to the C-numbers that participate in those reactions,
and enriches each with KEGG metadata (name, EC, pathways) + MNX cross-refs
(MNXR/MNXM, ChEBI, HMDB, InChIKey, formula, mass).

Output: cache/data/kegg/kegg_metabolism_xrefs.json (~1 MB).

This file is read by metabolism_adapter.py at KG build time. The 2.6 GB MNX
SQLite resolver is opened only here (prepare-data time), not at build time.

Usage:
    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]
"""
from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
OUTPUT_FILE = KEGG_CACHE_DIR / "kegg_metabolism_xrefs.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"

log = logging.getLogger(__name__)


def _collect_gene_reachable_ids(kegg_data: dict) -> tuple[set[str], set[str]]:
    """Walk all strains' gene_annotations_merged.json and KEGG link data.

    Returns (R_numbers, C_numbers): the gene-reachable subsets.
    """
    reachable_rxns: set[str] = set()
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for gene in genes.values():
            for rxn in gene.get("kegg_reactions", []) or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    reachable_rxns.add(rxn)

    rxn_to_cpds = kegg_data.get("reaction_to_compounds", {})
    reachable_cpds: set[str] = set()
    for rxn in reachable_rxns:
        for cpd in rxn_to_cpds.get(rxn, []):
            reachable_cpds.add(cpd)

    log.info(
        "Pruned to %d reactions, %d compounds",
        len(reachable_rxns),
        len(reachable_cpds),
    )
    return reachable_rxns, reachable_cpds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Overwrite existing output")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if OUTPUT_FILE.exists() and not args.force:
        log.info("Output already exists: %s (use --force to rebuild)", OUTPUT_FILE)
        return

    if not KEGG_DATA_FILE.exists():
        log.error("KEGG data not found: %s — run step 6 sub-step 1 first", KEGG_DATA_FILE)
        raise SystemExit(1)

    with open(KEGG_DATA_FILE) as f:
        kegg_data = json.load(f)

    reachable_rxns, reachable_cpds = _collect_gene_reachable_ids(kegg_data)

    # Placeholder: enrichment (MNX xrefs) will be added in Task 7
    xrefs = {
        "reactions": {r: {} for r in sorted(reachable_rxns)},
        "compounds": {c: {} for c in sorted(reachable_cpds)},
    }

    KEGG_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(xrefs, f, indent=2)
    log.info("Wrote %s", OUTPUT_FILE)


if __name__ == "__main__":
    main()
