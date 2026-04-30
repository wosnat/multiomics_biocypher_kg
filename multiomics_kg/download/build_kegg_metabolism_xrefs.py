"""Step 6 — unified pruned KEGG data cache.

Walks every strain's gene_annotations_merged.json to identify gene-reachable
KEGG entities, prunes the raw KEGG dataset to that subset, enriches reactions
and compounds with MNX cross-refs, and writes a single
cache/data/kegg/kegg_data.json (~500 KB).

Both kegg_anno_adapter and metabolism_adapter read this file. No per-adapter
pruning at iteration time.

Pathway-set rule (Option B): pathway IDs in the cache are exactly those reachable
from gene-KOs or gene-reactions. compound→pathway lists are filtered to that set,
so compound-only pathways (e.g. ko05140 Leishmaniasis) are excluded.

Usage:
    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]
"""
from __future__ import annotations

import argparse
import contextlib
import json
import logging
import sqlite3
from pathlib import Path

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils import kegg_utils
from multiomics_kg.utils import metabolite_utils as mu
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"

log = logging.getLogger(__name__)


# ── Reachability ──────────────────────────────────────────────────────────────

def _gene_reachable_sets(raw: dict) -> dict[str, set[str]]:
    """Walk all strains to compute the gene-reachable {KOs, reactions, pathways, compounds}.

    Returns a dict with keys 'kos', 'rxns', 'cpds', 'pws'.
    Pathway set is KO ∪ Rxn-reachable (Option B).
    """
    kos: set[str] = set()
    rxns: set[str] = set()
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for g in genes.values():
            for ko in g.get("kegg_ko") or []:
                if isinstance(ko, str) and ko.startswith("K"):
                    kos.add(ko)
            for rxn in g.get("kegg_reactions") or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    rxns.add(rxn)

    rxn_to_cpds = raw.get("reaction_to_compounds", {})
    cpds: set[str] = set()
    for r in rxns:
        cpds.update(rxn_to_cpds.get(r, []))

    ko_to_pw = raw.get("ko_to_pathways", {})
    rxn_to_pw = raw.get("reaction_to_pathways", {})
    pws: set[str] = set()
    for k in kos:
        pws.update(ko_to_pw.get(k, []))
    for r in rxns:
        pws.update(rxn_to_pw.get(r, []))

    log.info(
        f"Gene-reachable: {len(kos)} KOs, {len(rxns)} reactions, "
        f"{len(cpds)} compounds, {len(pws)} pathways (KO∪Rxn)"
    )
    return {"kos": kos, "rxns": rxns, "cpds": cpds, "pws": pws}


# ── Hierarchy parents ─────────────────────────────────────────────────────────

def _hierarchy_parents(raw: dict, pws: set[str]) -> tuple[set[str], set[str]]:
    """Compute the subcategory + category sets that need to be emitted."""
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})
    subs = {pw_to_sub[p] for p in pws if p in pw_to_sub}
    cats = {sub_to_cat[s] for s in subs if s in sub_to_cat}
    return subs, cats


# ── MNX enrichment ────────────────────────────────────────────────────────────

def _enrich_reaction(rxn_id: str, raw: dict, conn: sqlite3.Connection,
                     allowed_pathways: set[str]) -> dict:
    out: dict = {
        "name": raw.get("reaction_names", {}).get(rxn_id, ""),
        "pathways": [p for p in raw.get("reaction_to_pathways", {}).get(rxn_id, []) if p in allowed_pathways],
        "compounds": list(raw.get("reaction_to_compounds", {}).get(rxn_id, [])),
        "ec_numbers": [],
        "mnxr_id": None,
        "rhea_ids": [],
        "mass_balance": "unbalanced",
        "reaction_class": "chemical",
    }
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxr_id FROM reaction_aliases WHERE source='kegg.reaction' AND value=? LIMIT 1",
        (rxn_id,),
    )
    row = cur.fetchone()
    if row is None:
        return out
    out["mnxr_id"] = row[0]

    cur.execute(
        "SELECT value FROM reaction_aliases WHERE mnxr_id=? AND source='rhea' ORDER BY value",
        (out["mnxr_id"],),
    )
    out["rhea_ids"] = [r[0] for r in cur.fetchall()]

    cur.execute(
        "SELECT classifs, is_balanced, is_transport FROM reactions WHERE mnxr_id=?",
        (out["mnxr_id"],),
    )
    row = cur.fetchone()
    if row:
        classifs, is_balanced, is_transport = row
        if classifs:
            out["ec_numbers"] = [c.strip() for c in classifs.split(";") if c.strip()]
        out["mass_balance"] = "balanced" if (is_balanced or "").upper() == "B" else "unbalanced"
        out["reaction_class"] = "transport" if (is_transport or "").upper() == "T" else "chemical"
    return out


def _enrich_compound(cpd_id: str, raw: dict, conn: sqlite3.Connection,
                     allowed_pathways: set[str]) -> dict:
    out: dict = {
        "name": raw.get("compound_names", {}).get(cpd_id, ""),
        "formula": None,
        "mass": None,
        "inchikey": None,
        "smiles": None,
        "mnxm_id": None,
        "chebi_id": None,
        "hmdb_id": None,
        "pathways": [p for p in raw.get("compound_to_pathways", {}).get(cpd_id, []) if p in allowed_pathways],
    }
    cur = conn.cursor()
    cur.execute(
        "SELECT mnxm_id FROM compound_aliases WHERE source='kegg.compound' AND value=? LIMIT 1",
        (cpd_id,),
    )
    row = cur.fetchone()
    if row is None:
        return out
    out["mnxm_id"] = row[0]

    cur.execute(
        "SELECT formula, mass, inchikey, smiles FROM compounds WHERE mnxm_id=?",
        (out["mnxm_id"],),
    )
    row = cur.fetchone()
    if row:
        out["formula"] = row[0] or None
        out["mass"] = row[1] if row[1] is not None else None
        out["inchikey"] = row[2] or None
        out["smiles"] = row[3] or None

    cur.execute(
        "SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='chebi' ORDER BY value LIMIT 1",
        (out["mnxm_id"],),
    )
    row = cur.fetchone()
    if row:
        out["chebi_id"] = row[0]
    cur.execute(
        "SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='hmdb' ORDER BY value LIMIT 1",
        (out["mnxm_id"],),
    )
    row = cur.fetchone()
    if row:
        out["hmdb_id"] = row[0]
    return out


# ── Build ─────────────────────────────────────────────────────────────────────

def build_pruned_kegg_data(raw: dict, conn: sqlite3.Connection, out_path: Path) -> None:
    """Build the pruned kegg_data.json from raw KEGG data + MNX resolver."""
    sets = _gene_reachable_sets(raw)
    kos, rxns, cpds, pws = sets["kos"], sets["rxns"], sets["cpds"], sets["pws"]
    subs, cats = _hierarchy_parents(raw, pws)

    ko_to_pw = raw.get("ko_to_pathways", {})
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})

    out = {
        "kos": {
            k: {
                "name": raw.get("ko_names", {}).get(k, ""),
                "pathways": [p for p in ko_to_pw.get(k, []) if p in pws],
            } for k in sorted(kos)
        },
        "pathways": {
            p: {
                "name": raw.get("pathway_names", {}).get(p, ""),
                **({"subcategory": pw_to_sub[p]} if p in pw_to_sub else {}),
            } for p in sorted(pws)
        },
        "subcategories": {
            s: {
                "name": raw.get("subcategory_names", {}).get(s, ""),
                **({"category": sub_to_cat[s]} if s in sub_to_cat else {}),
            } for s in sorted(subs)
        },
        "categories": {
            c: {"name": raw.get("category_names", {}).get(c, "")}
            for c in sorted(cats)
        },
        "reactions": {
            r: _enrich_reaction(r, raw, conn, pws) for r in sorted(rxns)
        },
        "compounds": {
            c: _enrich_compound(c, raw, conn, pws) for c in sorted(cpds)
        },
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, separators=(",", ":")))
    log.info(
        f"Wrote {out_path}: {len(out['kos'])} KOs, {len(out['pathways'])} pathways, "
        f"{len(out['reactions'])} reactions, {len(out['compounds'])} compounds"
    )


# ── Raw → in-memory dict ──────────────────────────────────────────────────────

def _parse_brite_supplements(raw_dir: Path) -> dict:
    """Parse br_ko00001.json into the BRITE-derived supplementary maps.

    Returns: pathway_to_subcategory, subcategory_names, subcategory_to_category,
    category_names, plus _brite_pathway_names (merged into pathway_names in
    _parse_raw_into_dict via {**api_pw_names, **brite_pw_names} — BRITE wins on
    overlap, matching prior download_kegg_data() semantics).
    """
    brite_json = json.loads((raw_dir / "br_ko00001.json").read_text())
    pw_to_sub, sub_names, sub_to_cat, cat_names, brite_pw_names = (
        kegg_utils._parse_brite_hierarchy(brite_json)
    )
    return {
        "pathway_to_subcategory": pw_to_sub,
        "subcategory_names": sub_names,
        "subcategory_to_category": sub_to_cat,
        "category_names": cat_names,
        "_brite_pathway_names": brite_pw_names,
    }


def _parse_raw_into_dict(cache_root: Path) -> dict:
    """Parse raw KEGG cache files into the same dict shape that download_kegg_data
    used to return (so build_pruned_kegg_data can stay agnostic of the on-disk layout).
    """
    raw_dir = cache_root / "kegg" / "raw"

    api_pw_names = kegg_utils._parse_pathway_ko_names(
        (raw_dir / "list_pathway_ko.txt").read_text()
    )
    brite_supp = _parse_brite_supplements(raw_dir)
    # BRITE pathway names win on overlap, matching download_kegg_data() behaviour.
    pathway_names = {**api_pw_names, **brite_supp.pop("_brite_pathway_names")}

    return {
        "ko_names": kegg_utils._parse_ko_names(
            (raw_dir / "list_ko.txt").read_text()
        ),
        "pathway_names": pathway_names,
        "ko_to_pathways": kegg_utils._parse_ko_to_pathways(
            (raw_dir / "link_pathway_ko.txt").read_text()
        ),
        "reaction_names": kegg_utils._parse_reaction_names(
            (raw_dir / "list_reaction.txt").read_text()
        ),
        "compound_names": kegg_utils._parse_compound_names(
            (raw_dir / "list_compound.txt").read_text()
        ),
        "reaction_to_compounds": kegg_utils._parse_reaction_to_compounds(
            (raw_dir / "link_compound_reaction.txt").read_text()
        ),
        "reaction_to_pathways": kegg_utils._parse_reaction_to_pathways(
            (raw_dir / "link_pathway_reaction.txt").read_text()
        ),
        "compound_to_pathways": kegg_utils._parse_compound_to_pathways(
            (raw_dir / "link_pathway_compound.txt").read_text()
        ),
        **brite_supp,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cache_root = PROJECT_ROOT / "cache" / "data"

    if KEGG_DATA_FILE.exists() and not force:
        log.info(f"{KEGG_DATA_FILE} exists; use --force to rebuild.")
        return

    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run prepare_data.sh step 0 sub-step 7 first."
        )

    log.info("Ensuring KEGG raw cache (downloads from KEGG REST if missing) ...")
    kegg_utils.download_kegg_raw(cache_root, force=force)

    log.info("Parsing raw KEGG into in-memory dict ...")
    raw = _parse_raw_into_dict(cache_root)

    log.info("Opening MNX resolver ...")
    with contextlib.closing(mu.open_resolver(RESOLVER_DB)) as conn:
        log.info("Building pruned kegg_data.json ...")
        build_pruned_kegg_data(raw, conn, KEGG_DATA_FILE)
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if kegg_data.json exists.")
    args = parser.parse_args()
    main(force=args.force)
