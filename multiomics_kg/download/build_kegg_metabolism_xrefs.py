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
from multiomics_kg.utils import metabolite_utils as mu
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
OUTPUT_FILE = KEGG_CACHE_DIR / "kegg_metabolism_xrefs.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"

log = logging.getLogger(__name__)


def _collect_gene_reachable_ids(kegg_data: dict) -> tuple[set[str], set[str], set[str]]:
    """Walk all strains' gene_annotations_merged.json and KEGG link data.

    Returns (R_numbers, C_numbers, KO_pathways): gene-reachable subsets.

    `KO_pathways` is the set of pathway IDs that the kegg_anno adapter will emit
    as KeggTerm nodes — i.e. pathways referenced by at least one gene-annotated KO.
    Used to prune reaction.kegg_pathway_ids so we don't create dangling edges
    pointing to pathways with no KeggTerm node.
    """
    reachable_rxns: set[str] = set()
    reachable_kos: set[str] = set()
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for gene in genes.values():
            for rxn in gene.get("kegg_reactions", []) or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    reachable_rxns.add(rxn)
            for ko in gene.get("kegg_ko", []) or []:
                if isinstance(ko, str) and ko.startswith("K"):
                    reachable_kos.add(ko)

    rxn_to_cpds = kegg_data.get("reaction_to_compounds", {})
    reachable_cpds: set[str] = set()
    for rxn in reachable_rxns:
        for cpd in rxn_to_cpds.get(rxn, []):
            reachable_cpds.add(cpd)

    ko_to_pw = kegg_data.get("ko_to_pathways", {})
    reachable_pws: set[str] = set()
    for ko in reachable_kos:
        reachable_pws.update(ko_to_pw.get(ko, []))

    log.info(
        "Pruned to %d reactions, %d compounds, "
        "%d KO-reachable pathways (for pathway-edge pruning)",
        len(reachable_rxns),
        len(reachable_cpds),
        len(reachable_pws),
    )
    return reachable_rxns, reachable_cpds, reachable_pws


def _resolve_mnx_for_kegg_reaction(kegg_id: str, conn: sqlite3.Connection) -> str | None:
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxr_id FROM reaction_aliases WHERE source = 'kegg.reaction' AND value = ? LIMIT 1",
        (kegg_id,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _resolve_mnx_for_kegg_compound(kegg_id: str, conn: sqlite3.Connection) -> str | None:
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxm_id FROM compound_aliases WHERE source = 'kegg.compound' AND value = ? LIMIT 1",
        (kegg_id,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _query_aliases(table: str, mnx_col: str, mnx_id: str, source: str,
                   conn: sqlite3.Connection) -> list[str]:
    cur = conn.cursor()
    cur.execute(
        f"SELECT value FROM {table} WHERE {mnx_col} = ? AND source = ? ORDER BY value",
        (mnx_id, source),
    )
    return [r[0] for r in cur.fetchall()]


def _enrich_reaction(rxn_id: str, kegg_data: dict, conn: sqlite3.Connection,
                     valid_pathways: set[str]) -> dict:
    """Build the per-reaction enriched dict written to xrefs JSON."""
    raw_pws = kegg_data.get("reaction_to_pathways", {}).get(rxn_id, [])
    pruned_pws = [p for p in raw_pws if p in valid_pathways]
    out: dict = {
        "name": kegg_data.get("reaction_names", {}).get(rxn_id, ""),
        "compound_ids": kegg_data.get("reaction_to_compounds", {}).get(rxn_id, []),
        "kegg_pathway_ids": pruned_pws,
        "ec_numbers": [],
        "mnxr_id": None,
        "rhea_ids": [],
        "mass_balance": "unbalanced",
        "reaction_class": "chemical",
    }
    mnxr_id = _resolve_mnx_for_kegg_reaction(rxn_id, conn)
    if mnxr_id is None:
        return out
    out["mnxr_id"] = mnxr_id

    out["rhea_ids"] = _query_aliases("reaction_aliases", "mnxr_id", mnxr_id, "rhea", conn)

    cur = conn.cursor()
    cur.execute(
        "SELECT classifs, is_balanced, is_transport FROM reactions WHERE mnxr_id = ?",
        (mnxr_id,),
    )
    row = cur.fetchone()
    if row:
        classifs, is_balanced, is_transport = row
        if classifs:
            out["ec_numbers"] = [c.strip() for c in classifs.split(";") if c.strip()]
        out["mass_balance"] = "balanced" if (is_balanced or "").upper() == "B" else "unbalanced"
        out["reaction_class"] = "transport" if (is_transport or "").upper() == "T" else "chemical"
    return out


def _enrich_compound(cpd_id: str, kegg_data: dict, conn: sqlite3.Connection) -> dict:
    """Build the per-compound enriched dict written to xrefs JSON."""
    out: dict = {
        "name": kegg_data.get("compound_names", {}).get(cpd_id, ""),
        "mnxm_id": None,
        "formula": None,
        "mass": None,
        "inchikey": None,
        "smiles": None,
        "chebi_id": None,
        "hmdb_id": None,
    }
    mnxm_id = _resolve_mnx_for_kegg_compound(cpd_id, conn)
    if mnxm_id is None:
        return out
    out["mnxm_id"] = mnxm_id

    cur = conn.cursor()
    cur.execute(
        "SELECT formula, mass, inchikey, smiles FROM compounds WHERE mnxm_id = ?",
        (mnxm_id,),
    )
    row = cur.fetchone()
    if row:
        out["formula"] = row[0] or None
        out["mass"] = row[1] if row[1] is not None else None
        out["inchikey"] = row[2] or None
        out["smiles"] = row[3] or None

    chebi = _query_aliases("compound_aliases", "mnxm_id", mnxm_id, "chebi", conn)
    if chebi:
        out["chebi_id"] = chebi[0]
    hmdb = _query_aliases("compound_aliases", "mnxm_id", mnxm_id, "hmdb", conn)
    if hmdb:
        out["hmdb_id"] = hmdb[0]
    return out


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if OUTPUT_FILE.exists() and not force:
        log.info(f"{OUTPUT_FILE} exists; use --force to rebuild.")
        return

    if not KEGG_DATA_FILE.exists() or force:
        log.info("Loading KEGG hierarchy + metabolism endpoints (downloads from REST if not cached) ...")
        from multiomics_kg.utils.kegg_utils import download_kegg_data
        download_kegg_data(KEGG_CACHE_DIR.parent, force=force)
    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run prepare_data.sh step 0 sub-step 7 first."
        )

    log.info("Loading KEGG data ...")
    kegg_data = json.loads(KEGG_DATA_FILE.read_text())

    log.info("Pruning to gene-reachable subset ...")
    rxn_ids, cpd_ids, valid_pathways = _collect_gene_reachable_ids(kegg_data)

    log.info("Opening MNX resolver ...")
    conn = mu.open_resolver(RESOLVER_DB)

    log.info(f"Enriching {len(rxn_ids)} reactions ...")
    reactions = {
        rxn_id: _enrich_reaction(rxn_id, kegg_data, conn, valid_pathways)
        for rxn_id in sorted(rxn_ids)
    }

    log.info(f"Enriching {len(cpd_ids)} compounds ...")
    compounds = {cpd_id: _enrich_compound(cpd_id, kegg_data, conn) for cpd_id in sorted(cpd_ids)}

    conn.close()

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps({"reactions": reactions, "compounds": compounds}, indent=2))
    log.info(f"Wrote {OUTPUT_FILE} ({OUTPUT_FILE.stat().st_size:,} bytes, "
             f"{len(reactions)} reactions, {len(compounds)} compounds).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild the output cache even if it exists.")
    args = parser.parse_args()
    main(force=args.force)
