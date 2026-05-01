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

Step 6 also downloads the 3 TCDB reference TSVs (tc_classes, tc_subclasses,
families) and writes the assembled hierarchy to
cache/data/tcdb/tcdb_hierarchy.json. The TCDB lift is unconditional and not
gated on KEGG reachability.

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


# ── TCDB hierarchy ────────────────────────────────────────────────────────────

_TC_CLASS_NAMES = {
    "1": "Channels and Pores",
    "2": "Electrochemical Potential-driven Transporters",
    "3": "Primary Active Transporters",
    "4": "Group Translocators",
    "5": "Transmembrane Electron Carriers",
    "8": "Auxiliary Transport Proteins",
    "9": "Incompletely Characterized Transport Systems",
}


def _parse_tcdb_families(path: Path) -> dict[str, str]:
    """Parse families.tsv → {tc_family_id: description}."""
    out: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                out[parts[0]] = parts[1].strip()
    return out


def _parse_tcdb_superfamilies(path: Path) -> list[tuple[str, str, str, str, str]]:
    """Parse superfamilies.tsv → [(tcid, subfamily, family, abbreviation, superfamily), ...]."""
    rows: list[tuple[str, str, str, str, str]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                rows.append((parts[0], parts[1], parts[2], parts[3], parts[4]))
    return rows


def _parse_tcdb_substrates(path: Path) -> dict[str, list[str]]:
    """Parse substrates.tsv → {tcid_specificity: [substrate_name, ...]}.

    Substrate column is pipe-separated 'CHEBI:N;name|CHEBI:N;name'. We keep
    just the name (post-`;`) for substrate_classes; full CHEBI mapping is a
    Phase 1.3 enhancement.
    """
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            tcid, raw = parts
            names = []
            for entry in raw.split("|"):
                if ";" in entry:
                    _, _, name = entry.partition(";")
                    if name.strip():
                        names.append(name.strip())
            if names:
                out[tcid] = names
    return out


def build_tcdb_hierarchy(
    out_path: Path,
    families_path: Path,
    superfamilies_path: Path,
    substrates_path: Path,
) -> int:
    """Build tcdb_hierarchy.json by joining the 3 TCDB sources."""
    fam_descs = _parse_tcdb_families(families_path)
    super_rows = _parse_tcdb_superfamilies(superfamilies_path)
    substrates = _parse_tcdb_substrates(substrates_path)

    h: dict[str, dict] = {}

    # Synthesize class + subclass entries from any TCID prefix we observe
    seen_classes: set[str] = set()
    seen_subclasses: set[str] = set()

    for tcid, subfam, fam, abbr, superfam in super_rows:
        parts = tcid.split(".")
        # Need at least class.subclass.family.subfam.specificity = 5 parts
        if len(parts) < 5:
            continue
        cls = parts[0]
        subcls = ".".join(parts[:2])

        # Class (level 0)
        if cls not in seen_classes:
            h[cls] = {
                "name": _TC_CLASS_NAMES.get(cls, ""),
                "level": 0,
                "level_kind": "tc_class",
                "parent": None,
            }
            seen_classes.add(cls)

        # Subclass (level 1)
        if subcls not in seen_subclasses:
            h[subcls] = {
                "name": "",
                "level": 1,
                "level_kind": "tc_subclass",
                "parent": cls,
            }
            seen_subclasses.add(subcls)

        # Family (level 2)
        if fam not in h:
            h[fam] = {
                "name": fam_descs.get(fam, ""),
                "level": 2,
                "level_kind": "tc_family",
                "parent": subcls,
                "abbreviation": abbr,
            }

        # Subfamily (level 3)
        if subfam not in h:
            h[subfam] = {
                "name": "",
                "level": 3,
                "level_kind": "tc_subfamily",
                "parent": fam,
            }

        # Specificity (level 4)
        node: dict = {
            "name": "",
            "level": 4,
            "level_kind": "tc_specificity",
            "parent": subfam,
        }
        if tcid in substrates:
            node["substrate_classes"] = substrates[tcid]
        if superfam:
            node["superfamily"] = superfam
        h[tcid] = node

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(h, indent=2, sort_keys=True))
    log.info(f"  tcdb_hierarchy.json: {len(h)} entries")
    return len(h)


def _build_tcdb_hierarchy(cache_root: Path) -> int:
    """Parse the 3 TCDB TSVs into cache_root/tcdb/tcdb_hierarchy.json. Returns entry count."""
    tcdb_dir = cache_root / "tcdb"
    return build_tcdb_hierarchy(
        out_path=tcdb_dir / "tcdb_hierarchy.json",
        families_path=tcdb_dir / "families.tsv",
        superfamilies_path=tcdb_dir / "superfamilies.tsv",
        substrates_path=tcdb_dir / "substrates.tsv",
    )


# ── Reachability ──────────────────────────────────────────────────────────────

def _gene_reachable_sets(raw: dict) -> dict[str, set[str]]:
    """Walk all strains to compute the gene-reachable {KOs, reactions, pathways, compounds, TCDB IDs}.

    Returns a dict with keys 'kos', 'rxns', 'cpds', 'pws', 'tcdb_ids'.
    Pathway set is KO ∪ Rxn-reachable (Option B).
    """
    kos: set[str] = set()
    rxns: set[str] = set()
    tcdb_ids: set[str] = set()
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
            for tc in g.get("transporter_classification") or []:
                if isinstance(tc, str) and tc:
                    tcdb_ids.add(tc)

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
        f"{len(cpds)} compounds, {len(pws)} pathways (KO∪Rxn), "
        f"{len(tcdb_ids)} TCDB IDs"
    )
    return {"kos": kos, "rxns": rxns, "cpds": cpds, "pws": pws, "tcdb_ids": tcdb_ids}


# ── Hierarchy parents ─────────────────────────────────────────────────────────

def _hierarchy_parents(raw: dict, pws: set[str]) -> tuple[set[str], set[str]]:
    """Compute the subcategory + category sets that need to be emitted."""
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})
    subs = {pw_to_sub[p] for p in pws if p in pw_to_sub}
    cats = {sub_to_cat[s] for s in subs if s in sub_to_cat}
    return subs, cats


# ── MNX enrichment (batched) ──────────────────────────────────────────────────

# SQLite default: 999 host parameters per statement. Use 900 for safety.
_BATCH_SIZE = 900


def _batched(seq: list, n: int = _BATCH_SIZE):
    """Yield successive n-sized chunks from seq."""
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def _bulk_enrich_reactions(
    conn: sqlite3.Connection,
    kegg_ids: list[str],
    allowed_pathways: set[str],
    raw: dict,
) -> dict[str, dict]:
    """Batched enrichment for all KEGG reactions in one go.

    Returns dict mapping kegg_reaction_id → enrichment dict (same shape as the
    per-entity _enrich_reaction output). Reactions with no MNX entry get a stub
    with mnxr_id=None and default fields, so every requested ID appears in the
    output.

    Replaces ~3 queries × N reactions with ~3 batched queries total (per chunk).
    """
    cur = conn.cursor()

    # Initialize stubs for every requested reaction
    rxn_names = raw.get("reaction_names", {})
    rxn_to_pw = raw.get("reaction_to_pathways", {})
    rxn_to_cpds = raw.get("reaction_to_compounds", {})

    out: dict[str, dict] = {}
    for rxn_id in kegg_ids:
        out[rxn_id] = {
            "name": rxn_names.get(rxn_id, ""),
            "pathways": [p for p in rxn_to_pw.get(rxn_id, []) if p in allowed_pathways],
            "compounds": list(rxn_to_cpds.get(rxn_id, [])),
            "ec_numbers": [],
            "mnxr_id": None,
            "rhea_ids": [],
            "mass_balance": "unbalanced",
            "reaction_class": "chemical",
        }

    # Phase 1: bulk-resolve kegg.reaction → mnxr_id
    kegg_to_mnxr: dict[str, str] = {}
    for chunk in _batched(kegg_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxr_id FROM reaction_aliases "
            f"WHERE source='kegg.reaction' AND value IN ({placeholders})",
            chunk,
        )
        for value, mnxr_id in cur.fetchall():
            kegg_to_mnxr[value] = mnxr_id

    # Apply phase-1 results
    for rxn_id, mnxr_id in kegg_to_mnxr.items():
        out[rxn_id]["mnxr_id"] = mnxr_id

    if not kegg_to_mnxr:
        return out

    # Phase 2: bulk-fetch MNX-side properties keyed by mnxr_id
    mnxr_ids = list(set(kegg_to_mnxr.values()))

    # 2a. reactions table → ec_numbers, mass_balance, reaction_class
    rxn_props: dict[str, dict] = {}
    for chunk in _batched(mnxr_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxr_id, classifs, is_balanced, is_transport "
            f"FROM reactions WHERE mnxr_id IN ({placeholders})",
            chunk,
        )
        for mnxr_id, classifs, is_balanced, is_transport in cur.fetchall():
            rxn_props[mnxr_id] = {
                "ec_numbers": (
                    [c.strip() for c in classifs.split(";") if c.strip()]
                    if classifs else []
                ),
                "mass_balance": "balanced" if (is_balanced or "").upper() == "B" else "unbalanced",
                "reaction_class": "transport" if (is_transport or "").upper() == "T" else "chemical",
            }

    # 2b. reaction_aliases → rhea_ids
    rhea_by_mnxr: dict[str, list[str]] = {}
    for chunk in _batched(mnxr_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxr_id, value FROM reaction_aliases "
            f"WHERE source='rhea' AND mnxr_id IN ({placeholders}) ORDER BY value",
            chunk,
        )
        for mnxr_id, value in cur.fetchall():
            rhea_by_mnxr.setdefault(mnxr_id, []).append(value)

    # Stitch phase-2 results into per-reaction output
    for rxn_id, mnxr_id in kegg_to_mnxr.items():
        props = rxn_props.get(mnxr_id, {})
        out[rxn_id]["ec_numbers"] = props.get("ec_numbers", [])
        out[rxn_id]["mass_balance"] = props.get("mass_balance", "unbalanced")
        out[rxn_id]["reaction_class"] = props.get("reaction_class", "chemical")
        out[rxn_id]["rhea_ids"] = rhea_by_mnxr.get(mnxr_id, [])

    return out


def _bulk_enrich_compounds(
    conn: sqlite3.Connection,
    kegg_ids: list[str],
    allowed_pathways: set[str],
    raw: dict,
) -> dict[str, dict]:
    """Batched enrichment for all KEGG compounds.

    Same pattern as _bulk_enrich_reactions: phase 1 resolves kegg.compound →
    mnxm_id, phase 2 bulk-fetches MNX-side properties (formula/mass/inchikey/
    smiles + chebi/hmdb aliases), then assembles per-compound output.
    """
    cur = conn.cursor()

    cpd_names = raw.get("compound_names", {})
    cpd_to_pw = raw.get("compound_to_pathways", {})

    out: dict[str, dict] = {}
    for cpd_id in kegg_ids:
        out[cpd_id] = {
            "name": cpd_names.get(cpd_id, ""),
            "formula": None, "mass": None, "inchikey": None, "smiles": None,
            "mnxm_id": None, "chebi_id": None, "hmdb_id": None,
            "pathways": [p for p in cpd_to_pw.get(cpd_id, []) if p in allowed_pathways],
        }

    # Phase 1: bulk-resolve kegg.compound → mnxm_id
    kegg_to_mnxm: dict[str, str] = {}
    for chunk in _batched(kegg_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxm_id FROM compound_aliases "
            f"WHERE source='kegg.compound' AND value IN ({placeholders})",
            chunk,
        )
        for value, mnxm_id in cur.fetchall():
            kegg_to_mnxm[value] = mnxm_id

    for cpd_id, mnxm_id in kegg_to_mnxm.items():
        out[cpd_id]["mnxm_id"] = mnxm_id

    if not kegg_to_mnxm:
        return out

    mnxm_ids = list(set(kegg_to_mnxm.values()))

    # 2a. compounds table → formula, mass, inchikey, smiles
    cpd_props: dict[str, dict] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, formula, mass, inchikey, smiles "
            f"FROM compounds WHERE mnxm_id IN ({placeholders})",
            chunk,
        )
        for mnxm_id, formula, mass, inchikey, smiles in cur.fetchall():
            cpd_props[mnxm_id] = {
                "formula": formula or None,
                "mass": mass,  # numeric, no `or None` coercion needed
                "inchikey": inchikey or None,
                "smiles": smiles or None,
            }

    # 2b. compound_aliases → chebi_id (one lowest-id per mnxm)
    chebi_by_mnxm: dict[str, str] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='chebi' AND mnxm_id IN ({placeholders}) ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            # Keep first (lowest-sorted) chebi id per mnxm
            chebi_by_mnxm.setdefault(mnxm_id, value)

    # 2c. compound_aliases → hmdb_id (same pattern)
    hmdb_by_mnxm: dict[str, str] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='hmdb' AND mnxm_id IN ({placeholders}) ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            hmdb_by_mnxm.setdefault(mnxm_id, value)

    # Stitch phase-2 results into per-compound output
    for cpd_id, mnxm_id in kegg_to_mnxm.items():
        props = cpd_props.get(mnxm_id, {})
        out[cpd_id]["formula"] = props.get("formula")
        out[cpd_id]["mass"] = props.get("mass")
        out[cpd_id]["inchikey"] = props.get("inchikey")
        out[cpd_id]["smiles"] = props.get("smiles")
        out[cpd_id]["chebi_id"] = chebi_by_mnxm.get(mnxm_id)
        out[cpd_id]["hmdb_id"] = hmdb_by_mnxm.get(mnxm_id)

    return out


# ── TCDB pruning + substrate resolution ──────────────────────────────────────


def _prune_tcdb(
    hierarchy: dict,
    seed_ids: set[str],
) -> tuple[set[str], dict[str, list[str]]]:
    """Bidirectional prune of the TCDB hierarchy.

    For each seed ID present in `hierarchy`, walk **up** to the tc_class root
    and **down** to all tc_specificity leaves. Returns:

      kept: set of all TCDB IDs that survive pruning.
      leaf_subs: {leaf_tc_id: [substrate_str, ...]} for kept leaves that carry
                 substrate_classes.

    Seeds not present in the hierarchy are silently skipped.

    Cycle prevention for the recursive descent uses a separate `down_visited`
    set so the guard is correct even when a node was added to `kept` by a
    walk-up pass before any of its descendants have been visited.
    """
    kept: set[str] = set()
    leaf_subs: dict[str, list[str]] = {}
    parent_of = {tc: data.get("parent") for tc, data in hierarchy.items()}
    children_of: dict[str, list[str]] = {}
    for tc, parent in parent_of.items():
        if parent is not None:
            children_of.setdefault(parent, []).append(tc)

    def walk_up(tc: str) -> None:
        cur: str | None = tc
        while cur is not None and cur in hierarchy:
            if cur in kept:
                # Already walked up from here.
                return
            kept.add(cur)
            cur = parent_of.get(cur)

    down_visited: set[str] = set()

    def walk_down(tc: str) -> None:
        if tc in down_visited:
            return
        down_visited.add(tc)
        if tc not in hierarchy:
            return
        kept.add(tc)
        node = hierarchy[tc]
        if node.get("level_kind") == "tc_specificity":
            subs = node.get("substrate_classes")
            if subs:
                leaf_subs[tc] = list(subs)
        for child in children_of.get(tc, []):
            walk_down(child)

    for seed in seed_ids:
        if seed not in hierarchy:
            continue
        walk_up(seed)
        walk_down(seed)
    return kept, leaf_subs


def _resolve_substrates(
    leaf_subs: dict[str, list[str]],
    conn: sqlite3.Connection,
) -> tuple[dict[str, list[str]], dict[str, dict]]:
    """Resolve `CHEBI:NNNN;name` substrate strings to Metabolite primary IDs.

    Returns:
      leaf_to_primary_ids: {leaf_tcid: [primary_id, ...]}, sorted+deduped.
      compound_props: {primary_id: {name, formula, mass, inchikey, mnxm_id, chebi_id}}.
    """
    from multiomics_kg.utils.metabolite_utils import (
        resolve_metabolite, mnxm_to_primary_id,
    )

    leaf_to_primary_ids: dict[str, list[str]] = {}
    compound_props: dict[str, dict] = {}
    cur = conn.cursor()

    for leaf, subs in leaf_subs.items():
        primary_ids: list[str] = []
        for sub_str in subs:
            # Format: 'CHEBI:NNNN;name' or just 'name'
            if ":" not in sub_str:
                continue
            chebi_part, _, name_part = sub_str.partition(";")
            mnxm, _method = resolve_metabolite(chebi_part, conn)
            if mnxm is None:
                log.debug(f"TCDB substrate unresolved: {sub_str!r} (leaf {leaf})")
                continue
            primary = mnxm_to_primary_id(mnxm, conn)
            primary_ids.append(primary)

            if primary in compound_props:
                # Already collected props for this primary ID — skip the lookup.
                continue

            # Pull MNX-side properties for the additional_compounds entry
            cur.execute(
                "SELECT name, formula, mass, inchikey FROM compounds WHERE mnxm_id = ?",
                (mnxm,),
            )
            row = cur.fetchone()
            mnx_name, formula, mass, inchikey = row if row else (name_part, None, None, None)
            cur.execute(
                "SELECT value FROM compound_aliases "
                "WHERE source='chebi' AND mnxm_id = ? ORDER BY value LIMIT 1",
                (mnxm,),
            )
            chebi_row = cur.fetchone()
            chebi_id = chebi_row[0] if chebi_row else None

            compound_props[primary] = {
                "name": (mnx_name or name_part).strip(),
                "formula": formula or None,
                "mass": mass,
                "inchikey": inchikey or None,
                "mnxm_id": mnxm,
                "chebi_id": chebi_id,
            }
        leaf_to_primary_ids[leaf] = sorted(set(primary_ids))

    return leaf_to_primary_ids, compound_props


def _fold_substrates_into_kegg_data(
    kegg_data: dict,
    leaf_to_primary_ids: dict[str, list[str]],
    compound_props: dict[str, dict],
) -> None:
    """Mutate kegg_data in place: tag compounds with evidence_sources, add transport-only
    chemistry under additional_compounds.

    1. Every entry in kegg_data['compounds'] gets 'metabolism' in its evidence_sources.
    2. For each transport-substrate primary ID:
       - If kegg.compound:* and already in kegg_data['compounds']: append 'transport'.
       - Else: add to kegg_data['additional_compounds'] with evidence_sources=['transport'].
    3. Pre-existing additional_compounds entries' evidence_sources lists are unioned.
    """
    # 1. Tag existing compounds
    for cpd in kegg_data.setdefault("compounds", {}).values():
        sources = cpd.setdefault("evidence_sources", [])
        if "metabolism" not in sources:
            sources.append("metabolism")

    additional: dict[str, dict] = kegg_data.setdefault("additional_compounds", {})

    # 2. Visit every transport substrate
    for primary_ids in leaf_to_primary_ids.values():
        for primary_id in primary_ids:
            if primary_id.startswith("kegg.compound:"):
                kegg_cpd_id = primary_id[len("kegg.compound:"):]
                if kegg_cpd_id in kegg_data["compounds"]:
                    sources = kegg_data["compounds"][kegg_cpd_id]["evidence_sources"]
                    if "transport" not in sources:
                        sources.append("transport")
                    continue
            # Else: add to additional_compounds
            if primary_id in additional:
                if "transport" not in additional[primary_id].setdefault("evidence_sources", []):
                    additional[primary_id]["evidence_sources"].append("transport")
            elif primary_id in compound_props:
                additional[primary_id] = {
                    **compound_props[primary_id],
                    "evidence_sources": ["transport"],
                }


# ── Build ─────────────────────────────────────────────────────────────────────

def build_pruned_kegg_data(raw: dict, conn: sqlite3.Connection, out_path: Path) -> None:
    """Build the pruned kegg_data.json from raw KEGG data + MNX resolver."""
    sets = _gene_reachable_sets(raw)
    kos, rxns, cpds, pws = sets["kos"], sets["rxns"], sets["cpds"], sets["pws"]
    subs, cats = _hierarchy_parents(raw, pws)

    ko_to_pw = raw.get("ko_to_pathways", {})
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})

    rxn_enriched = _bulk_enrich_reactions(conn, sorted(rxns), pws, raw)
    cpd_enriched = _bulk_enrich_compounds(conn, sorted(cpds), pws, raw)

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
        "reactions": rxn_enriched,
        "compounds": cpd_enriched,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
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
    overlap).
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
    """Parse raw KEGG cache files into a flat in-memory dict shape (so
    build_pruned_kegg_data can stay agnostic of the on-disk layout)."""
    raw_dir = cache_root / "kegg" / "raw"

    api_pw_names = kegg_utils._parse_pathway_ko_names(
        (raw_dir / "list_pathway_ko.txt").read_text()
    )
    brite_supp = _parse_brite_supplements(raw_dir)
    # BRITE pathway names win on overlap.
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
            f"{RESOLVER_DB} missing — run `bash scripts/refresh_mnx.sh` first."
        )

    log.info("Ensuring KEGG raw cache (downloads from KEGG REST if missing) ...")
    kegg_utils.download_kegg_raw(cache_root, force=force)

    log.info("Ensuring TCDB reference TSVs are downloaded ...")
    from multiomics_kg.download.download_metabolism_reference import download_all
    download_all(cache_root=cache_root, force=force, sources=["tcdb"])

    log.info("Building tcdb_hierarchy.json ...")
    _build_tcdb_hierarchy(cache_root)

    log.info("Parsing raw KEGG into in-memory dict ...")
    raw = _parse_raw_into_dict(cache_root)

    log.info("Opening MNX resolver ...")
    with contextlib.closing(mu.open_resolver(RESOLVER_DB)) as conn:
        log.info("Building pruned kegg_data.json ...")
        build_pruned_kegg_data(raw, conn, KEGG_DATA_FILE)

        log.info("Pruning TCDB hierarchy + resolving transport substrates ...")
        hierarchy = json.loads(
            (cache_root / "tcdb" / "tcdb_hierarchy.json").read_text()
        )
        # Re-collect gene-reachable sets (cheap — just walks gene_annotations files)
        # to get the TCDB seed set.
        sets = _gene_reachable_sets(raw)
        seed_tcdb_ids = sets["tcdb_ids"]

        kept, leaf_subs = _prune_tcdb(hierarchy, seed_tcdb_ids)
        leaf_to_primary, compound_props = _resolve_substrates(leaf_subs, conn)

        kegg_data = json.loads(KEGG_DATA_FILE.read_text())
        _fold_substrates_into_kegg_data(kegg_data, leaf_to_primary, compound_props)
        KEGG_DATA_FILE.write_text(json.dumps(kegg_data, indent=2, sort_keys=True))

        (cache_root / "tcdb" / "tcdb_pruned.json").write_text(json.dumps({
            "kept_tcdb_ids": sorted(kept),
            "leaf_substrates": {
                leaf: ids for leaf, ids in sorted(leaf_to_primary.items())
            },
        }, indent=2, sort_keys=True))
        log.info(
            f"  TCDB: {len(kept)} kept IDs, "
            f"{len(leaf_to_primary)} leaves with substrates"
        )
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if kegg_data.json exists.")
    args = parser.parse_args()
    main(force=args.force)
