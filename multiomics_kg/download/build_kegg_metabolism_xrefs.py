"""Step 6 — unified pruned KEGG data cache.

Walks every strain's gene_annotations_merged.json to identify gene-reachable
KEGG entities, prunes the raw KEGG dataset to that subset, enriches reactions
and compounds with MNX cross-refs, and writes a single
cache/data/kegg/kegg_data.json (~3-4 MB indented for git-friendly diffs).

Both kegg_anno_adapter and metabolism_adapter read this file. No per-adapter
pruning at iteration time.

Pathway-set rule (Option B): pathway IDs in the cache are exactly those reachable
from gene-KOs, gene-reactions, OR gene-annotated TCDB substrates resolved to KEGG
compounds (transport-reachable). compound→pathway lists are filtered to that set,
so compound-only pathways (e.g. ko05140 Leishmaniasis) are excluded.

Compound-set rule: catalysis_cpds (gene-Reaction → compound) are unioned with
substrate_kegg_cpds (gene-TCDB → MNX-resolved kegg.compound:* primary). The
resulting extended cpds set rides through `_bulk_enrich_compounds` so transport-
only KEGG compounds get their pathways/MNX cross-refs natively. Non-KEGG
substrate primaries (chebi:*, mnx:*) become additional_compounds entries.

Step 6 also downloads the 3 TCDB reference TSVs (tc_classes, tc_subclasses,
families) and writes the assembled hierarchy to
cache/data/tcdb/tcdb_hierarchy.json. The TCDB lift is unconditional and not
gated on KEGG reachability.

Usage:
    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force] [--refetch-raw]

  --force         Rebuild kegg_data.json from already-cached raw inputs.
                  Does NOT re-download raw KEGG REST or TCDB TSVs — fast iteration
                  loop for paperconfig / metabolite_aliases.yaml edits.
  --refetch-raw   Re-download raw KEGG REST + TCDB TSVs, then rebuild. Implies
                  --force. Use only when KEGG / TCDB upstream releases.
"""
from __future__ import annotations

import argparse
import contextlib
import json
import logging
import sqlite3
from pathlib import Path

import pandas as pd
import yaml

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils import kegg_utils
from multiomics_kg.utils import metabolite_utils as mu
from multiomics_kg.utils.gene_id_utils import load_gene_annotations
from multiomics_kg.utils.metabolite_utils import mnxm_to_primary_id, resolve_metabolite
from multiomics_kg.download.download_metabolism_reference import download_all
from multiomics_kg.utils.paperconfig_utils import (
    iter_metabolite_assays_tables,
    load_all_paperconfigs,
)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
RESOLVER_DB = mu.get_mnx_data_dir() / "metabolite_resolver.db"

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
    """Parse substrates.tsv → {tcid_specificity: ['CHEBI:N;name', ...]}.

    Substrate column is pipe-separated 'CHEBI:N;name|CHEBI:N;name'. The full
    string is preserved (CHEBI: prefix included) so downstream MNX resolution
    in _resolve_substrates can map substrates to canonical Metabolite primary
    IDs. Earlier revisions stripped the prefix at parse time, leaving only
    the human-readable name — which silently broke transport-substrate
    resolution (every leaf returned 0 primary IDs).
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
            substrates = [entry.strip() for entry in raw.split("|") if entry.strip()]
            if substrates:
                out[tcid] = substrates
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
    """Parse the 3 TCDB TSVs into cache_root/tcdb/tcdb_hierarchy.json. Returns entry count.

    Raw TSVs live under cache_root/tcdb/raw/ (mirrors cache_root/kegg/raw/) so
    that the committed step-6 outputs (tcdb_hierarchy.json, tcdb_pruned.json)
    sit alongside but separately from the gitignored raw downloads.
    """
    tcdb_dir = cache_root / "tcdb"
    raw_dir = tcdb_dir / "raw"
    return build_tcdb_hierarchy(
        out_path=tcdb_dir / "tcdb_hierarchy.json",
        families_path=raw_dir / "families.tsv",
        superfamilies_path=raw_dir / "superfamilies.tsv",
        substrates_path=raw_dir / "substrates.tsv",
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

    skipped_seeds: list[str] = []
    for seed in seed_ids:
        if seed not in hierarchy:
            skipped_seeds.append(seed)
            continue
        walk_up(seed)
        walk_down(seed)
    if skipped_seeds:
        sample = ", ".join(sorted(skipped_seeds)[:5])
        log.info(
            f"  TCDB: {len(skipped_seeds)} gene-annotated seed(s) not in curated "
            f"hierarchy (e.g. {sample})"
        )
    return kept, leaf_subs


def _resolve_substrates(
    leaf_subs: dict[str, list[str]],
    conn: sqlite3.Connection,
) -> tuple[dict[str, list[str]], dict[str, dict]]:
    """Resolve `CHEBI:NNNN;name` substrate strings to Metabolite primary IDs.

    Returns:
      leaf_to_primary_ids: {leaf_tcid: [primary_id, ...]}, sorted+deduped.
      compound_props: {primary_id: {name, formula, mass, inchikey, mnxm_id, chebi_id}}
                      for non-KEGG primaries only (KEGG primaries flow through
                      `_bulk_enrich_compounds`, so they don't need props here).

    Batched against the MNX SQLite — one IN-clause query per phase instead of
    a per-substrate loop with 3-4 roundtrips each. Roughly 17K queries → ~6.
    """
    cur = conn.cursor()

    # Pass 1: parse all substrate strings → (leaf, chebi_value, name_part).
    # Format: 'CHEBI:NNNN;name' or just 'name'. MNX stores CHEBI aliases as
    # (source='chebi', value='<bare-number>'), so strip the prefix here.
    parsed: list[tuple[str, str, str]] = []
    for leaf, subs in leaf_subs.items():
        for sub_str in subs:
            if ":" not in sub_str:
                continue
            chebi_part, _, name_part = sub_str.partition(";")
            prefix, _, chebi_value = chebi_part.partition(":")
            if prefix.strip().upper() != "CHEBI" or not chebi_value:
                log.debug(f"TCDB substrate not in CHEBI form: {sub_str!r} (leaf {leaf})")
                continue
            parsed.append((leaf, chebi_value.strip(), name_part))

    # Always emit an entry per leaf, even if all substrates fail to parse/resolve
    leaf_to_primary_ids: dict[str, list[str]] = {leaf: [] for leaf in leaf_subs}
    compound_props: dict[str, dict] = {}
    if not parsed:
        return leaf_to_primary_ids, compound_props

    # Pass 2: chebi_value → mnxm_id (batched). First-by-mnxm_id wins on collision.
    distinct_chebi = sorted({chebi for _, chebi, _ in parsed})
    chebi_to_mnxm: dict[str, str] = {}
    for chunk in _batched(distinct_chebi):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxm_id FROM compound_aliases "
            f"WHERE source='chebi' AND value IN ({placeholders}) "
            f"ORDER BY mnxm_id",
            chunk,
        )
        for value, mnxm_id in cur.fetchall():
            chebi_to_mnxm.setdefault(value, mnxm_id)

    # Pass 3: distinct mnxm_ids → kegg.compound + chebi aliases (priority order
    # for primary-ID resolution; mirrors mnxm_to_primary_id).
    mnxm_ids = sorted(set(chebi_to_mnxm.values()))
    mnxm_to_kegg: dict[str, str] = {}
    mnxm_to_chebi: dict[str, str] = {}
    for chunk in _batched(mnxm_ids):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='kegg.compound' AND mnxm_id IN ({placeholders}) "
            f"ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            mnxm_to_kegg.setdefault(mnxm_id, value)
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='chebi' AND mnxm_id IN ({placeholders}) "
            f"ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            mnxm_to_chebi.setdefault(mnxm_id, value)

    def _primary_for(mnxm: str) -> str:
        if mnxm in mnxm_to_kegg:
            return f"kegg.compound:{mnxm_to_kegg[mnxm]}"
        if mnxm in mnxm_to_chebi:
            return f"chebi:{mnxm_to_chebi[mnxm]}"
        return f"mnx:{mnxm}"

    # Pass 4: bulk-fetch compounds props ONLY for mnxm_ids whose primary is non-KEGG
    # (KEGG primaries get their props from _bulk_enrich_compounds via the extended
    # cpds set, so we can skip them here).
    nonkegg_mnxms = [m for m in mnxm_ids if m not in mnxm_to_kegg]
    mnxm_to_props: dict[str, tuple] = {}
    for chunk in _batched(nonkegg_mnxms):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT mnxm_id, name, formula, mass, inchikey "
            f"FROM compounds WHERE mnxm_id IN ({placeholders})",
            chunk,
        )
        for mnxm_id, name, formula, mass, inchikey in cur.fetchall():
            mnxm_to_props[mnxm_id] = (name, formula, mass, inchikey)

    # Pass 5: assemble outputs in original iteration order so first-occurrence
    # name_part wins for non-KEGG primaries (matches pre-batch semantics).
    for leaf, chebi_value, name_part in parsed:
        mnxm = chebi_to_mnxm.get(chebi_value)
        if mnxm is None:
            log.debug(f"TCDB substrate unresolved: chebi:{chebi_value} (leaf {leaf})")
            continue
        primary = _primary_for(mnxm)
        leaf_to_primary_ids[leaf].append(primary)

        if primary.startswith("kegg.compound:") or primary in compound_props:
            continue
        mnx_name, formula, mass, inchikey = mnxm_to_props.get(mnxm, (None, None, None, None))
        compound_props[primary] = {
            "name": (mnx_name or name_part).strip(),
            "formula": formula or None,
            "mass": mass,
            "inchikey": inchikey or None,
            "mnxm_id": mnxm,
            "chebi_id": mnxm_to_chebi.get(mnxm),
        }

    for leaf in leaf_to_primary_ids:
        leaf_to_primary_ids[leaf] = sorted(set(leaf_to_primary_ids[leaf]))

    return leaf_to_primary_ids, compound_props


def _harvest_paper_metabolites(
    paperconfigs: list[tuple[Path, dict]],
    conn: sqlite3.Connection,
) -> dict:
    """Walk all metabolite_assays_table entries; resolve every metabolite name.

    Resolution priority per row:
      1. id_col (when set + non-empty in the row): parse per id_type → primary ID
      2. aliases_file (paper-local YAML map): name → primary ID override
      3. MNX resolver against name_col
      4. else → unresolved

    Returns a dict:
      alias_to_primary:    {name_or_alias: primary_id}
      resolution_methods:  {name_or_alias: method_str}
      unresolved:          sorted list[str]
      paper_kegg_cpds:     set[str]   — KEGG-prefix-stripped C-numbers (mirrors substrate_kegg_cpds)
      paper_non_kegg_compounds: {primary_id: {name, formula, mass, inchikey, mnxm_id, chebi_id}}
      per_paper:           {paperconfig_path: {resolved, unresolved, total}}
    """
    alias_to_primary: dict[str, str] = {}
    resolution_methods: dict[str, str] = {}
    unresolved: set[str] = set()
    paper_kegg_cpds: set[str] = set()
    paper_non_kegg_compounds: dict[str, dict] = {}
    per_paper: dict[str, dict[str, int]] = {}

    cur = conn.cursor()

    def _record(name: str, primary: str | None, method: str) -> bool:
        if primary is None:
            unresolved.add(name)
            resolution_methods.setdefault(name, "unresolved")
            return False
        alias_to_primary[name] = primary
        resolution_methods[name] = method
        if primary.startswith("kegg.compound:"):
            paper_kegg_cpds.add(primary[len("kegg.compound:"):])
        elif primary not in paper_non_kegg_compounds:
            # Pull MNX-side props for the additional_compounds entry
            mnxm = primary[len("mnx:"):] if primary.startswith("mnx:") else None
            if mnxm:
                cur.execute(
                    "SELECT name, formula, mass, inchikey FROM compounds WHERE mnxm_id = ?",
                    (mnxm,),
                )
                row = cur.fetchone()
                mnx_name, formula, mass, inchikey = row if row else (None, None, None, None)
            else:
                mnx_name, formula, mass, inchikey = None, None, None, None
            paper_non_kegg_compounds[primary] = {
                "name": (mnx_name or name).strip(),
                "formula": formula or None,
                "mass": mass,
                "inchikey": inchikey or None,
                "mnxm_id": mnxm,
                "chebi_id": primary[len("chebi:"):] if primary.startswith("chebi:") else None,
            }
        return True

    for pc_path, cfg in paperconfigs:
        for entry_key, entry in iter_metabolite_assays_tables(cfg):
            csv_path = entry.get("filename")
            if not csv_path:
                continue
            csv_full = (PROJECT_ROOT / csv_path) if not Path(csv_path).is_absolute() else Path(csv_path)
            if not csv_full.exists():
                log.warning(f"[paper-metabolites] CSV missing: {csv_full}")
                continue

            name_col = entry.get("name_col") or "compound"
            id_col = entry.get("id_col") or ""
            id_type = entry.get("id_type") or ""
            aliases_file = entry.get("aliases_file") or ""

            aliases: dict[str, str] = {}
            if aliases_file:
                aliases_path = pc_path.parent / aliases_file
                if aliases_path.exists():
                    try:
                        loaded = yaml.safe_load(aliases_path.read_text()) or {}
                        if isinstance(loaded, dict):
                            aliases = {str(k): str(v) for k, v in loaded.items()}
                    except yaml.YAMLError as e:
                        log.warning(f"[paper-metabolites] cannot parse {aliases_path}: {e}")

            try:
                skip_rows = int(entry.get("skip_rows", 0) or 0)
                df = pd.read_csv(csv_full, dtype=str, keep_default_na=False, skiprows=skip_rows or None)
            except Exception as e:
                log.warning(f"[paper-metabolites] cannot read {csv_full}: {e}")
                continue

            stats = per_paper.setdefault(str(pc_path), {"resolved": 0, "unresolved": 0, "total": 0})

            for _, row in df.iterrows():
                name = (row.get(name_col) or "").strip()
                if not name:
                    continue
                stats["total"] += 1
                if name in alias_to_primary:
                    stats["resolved"] += 1
                    continue

                # 1. id_col direct
                primary, method = None, None
                if id_col and id_type:
                    raw = (row.get(id_col) or "").strip()
                    # "0", "-", "NA", "N/A", "none" are common placeholders for
                    # "no ID assigned" in published metabolomics tables (e.g.
                    # Kujawinski 2023's KEGG column for unidentified compounds).
                    # Treating them as real IDs creates bogus compound entries
                    # like `compounds['0']`. Fall through to alias/MNX resolution.
                    if raw and raw.lower() not in {"0", "-", "na", "n/a", "none"}:
                        primary = raw if ":" in raw else f"{id_type}:{raw}"
                        method = (
                            "kegg_direct" if id_type == "kegg.compound"
                            else "chebi_direct" if id_type == "chebi"
                            else "mnx_direct" if id_type == "mnx"
                            else "id_direct"
                        )

                # 2. aliases_file override
                if primary is None and name in aliases:
                    primary, method = aliases[name], "alias_override"

                # 3. MNX resolver
                if primary is None:
                    mnxm, mnx_method = resolve_metabolite(name, conn)
                    if mnxm is not None:
                        primary = mnxm_to_primary_id(mnxm, conn)
                        method = (
                            "name_match" if mnx_method in ("xref:exact", "name:normalized")
                            else "ambiguous_multi_id" if mnx_method in ("xref:ambiguous", "ambiguous")
                            else "name_match"
                        )

                if _record(name, primary, method or "unresolved"):
                    stats["resolved"] += 1
                else:
                    stats["unresolved"] += 1

    return {
        "alias_to_primary": alias_to_primary,
        "resolution_methods": resolution_methods,
        "unresolved": sorted(unresolved),
        "paper_kegg_cpds": paper_kegg_cpds,
        "paper_non_kegg_compounds": paper_non_kegg_compounds,
        "per_paper": per_paper,
    }


def _fold_paper_metabolites_into_kegg_data(
    kegg_data: dict,
    *,
    paper_kegg_cpds: set[str],
    paper_non_kegg_compounds: dict[str, dict],
) -> None:
    """Add 'metabolomics' to evidence_sources for paper-measured compounds.

    KEGG-flavor paper compounds are already in kegg_data['compounds'] (because
    main() extended the cpds set before build_pruned_kegg_data ran). We just
    union 'metabolomics' into their evidence_sources.

    Non-KEGG paper compounds enter kegg_data['additional_compounds'] (or update
    an existing entry to add 'metabolomics' to its evidence_sources).
    """
    compounds = kegg_data.setdefault("compounds", {})
    additional = kegg_data.setdefault("additional_compounds", {})

    # paper_kegg_cpds elements are bare C-numbers (prefix-stripped by
    # _harvest_paper_metabolites). compounds is keyed by bare C-numbers too.
    # Every paper_kegg_cpd is in the extended cpds set passed to
    # build_pruned_kegg_data, so the entry must exist.
    for cpd_id in paper_kegg_cpds:
        entry = compounds[cpd_id]
        sources = list(entry.get("evidence_sources") or [])
        if "metabolomics" not in sources:
            sources.append("metabolomics")
        entry["evidence_sources"] = sources

    for primary_id, props in paper_non_kegg_compounds.items():
        if primary_id in additional:
            sources = list(additional[primary_id].get("evidence_sources") or [])
            if "metabolomics" not in sources:
                sources.append("metabolomics")
            additional[primary_id]["evidence_sources"] = sources
        else:
            additional[primary_id] = {
                **props,
                "evidence_sources": ["metabolomics"],
            }


def _write_metabolite_id_mapping(harvest: dict, out_path: Path) -> None:
    """Write metabolite_id_mapping.json. v1: only name_lookup populated.

    Schema reserves three tiers (specific_lookup, multi_lookup, name_lookup)
    + conflicts + compounds; future ID-column harvesting will fill the
    other tiers. v1 fills name_lookup only.
    """
    name_lookup: dict[str, list[str]] = {}
    compounds: dict[str, dict] = {}
    methods = harvest.get("resolution_methods") or {}
    for name, primary in (harvest.get("alias_to_primary") or {}).items():
        name_lookup.setdefault(name, []).append(primary)
        compounds.setdefault(primary, {"aliases": [], "first_method": methods.get(name, "")})
        if name not in compounds[primary]["aliases"]:
            compounds[primary]["aliases"].append(name)

    data = {
        "specific_lookup": {},
        "multi_lookup": {},
        "name_lookup": name_lookup,
        "conflicts": {},
        "compounds": compounds,
        "schema_version": 1,
        "per_paper_stats": harvest.get("per_paper") or {},
        "unresolved": harvest.get("unresolved") or [],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2, sort_keys=True))


def _fold_substrates_into_kegg_data(
    kegg_data: dict,
    *,
    catalysis_cpds: set[str],
    substrate_kegg_cpds: set[str],
    leaf_to_primary: dict[str, list[str]],
    compound_props: dict[str, dict],
) -> None:
    """Tag every kegg compound entry with evidence_sources and append non-KEGG
    transport substrates as additional_compounds entries.

    KEGG-flavor substrates (kegg.compound:C*) are already in kegg_data['compounds']
    because the caller passed an extended cpds set (catalysis_cpds ∪
    substrate_kegg_cpds) into build_pruned_kegg_data; their pathways and MNX
    cross-refs were populated by _bulk_enrich_compounds. This step only:

      1. Tags each compound entry with evidence_sources based on origin:
         - in catalysis_cpds only: ['metabolism']
         - in substrate_kegg_cpds only: ['transport']
         - in both: ['metabolism', 'transport']
      2. Adds non-KEGG substrate primaries (chebi:NNNN, mnx:MNXM*) to
         kegg_data['additional_compounds'] with evidence_sources=['transport'].
    """
    for cpd_id, entry in kegg_data.get("compounds", {}).items():
        sources = list(entry.get("evidence_sources") or [])
        if cpd_id in catalysis_cpds and "metabolism" not in sources:
            sources.append("metabolism")
        if cpd_id in substrate_kegg_cpds and "transport" not in sources:
            sources.append("transport")
        entry["evidence_sources"] = sources

    additional: dict[str, dict] = {}
    for primary_ids in leaf_to_primary.values():
        for primary_id in primary_ids:
            if primary_id.startswith("kegg.compound:"):
                continue
            if primary_id in additional or primary_id not in compound_props:
                continue
            additional[primary_id] = {
                **compound_props[primary_id],
                "evidence_sources": ["transport"],
            }
    kegg_data["additional_compounds"] = additional


# ── Build ─────────────────────────────────────────────────────────────────────

def build_pruned_kegg_data(
    raw: dict,
    conn: sqlite3.Connection,
    *,
    sets: dict[str, set[str]],
) -> dict:
    """Build the pruned kegg_data dict from raw KEGG data + MNX resolver.

    Returns the assembled dict; the caller is responsible for writing it to
    disk (after any final folds — evidence-source tagging, paper-metabolite
    union — happen in-place on the returned dict).

    `sets` matches the dict shape returned by `_gene_reachable_sets`, optionally
    extended by the caller with transport-reachable cpds + pws. Keys consumed:
    'kos', 'rxns', 'cpds', 'pws'. The caller controls extension so transport
    substrates (TCDB → MNX → kegg.compound:C*) can ride through the regular
    `_bulk_enrich_compounds` pipeline and pick up their pathway/MNX cross-refs.

    Filters `pws` to pathways with non-empty names. KEGG has a class of "summary
    of natural product biosynthesis" maps (ko010**) and drug-classification maps
    (ko07***) that have NO direct KOs/reactions and so don't appear in
    `/list/pathway/ko` (our names source). Transport-substrate compound→pathway
    extension can drag these in; dropping them avoids nameless KeggTerm nodes.
    """
    pathway_names = raw.get("pathway_names", {})
    kos, rxns, cpds = sets["kos"], sets["rxns"], sets["cpds"]
    pws_in = sets["pws"]
    pws = {p for p in pws_in if pathway_names.get(p)}
    dropped_pws = pws_in - pws
    if dropped_pws:
        log.info(
            f"Dropped {len(dropped_pws)} nameless pathways "
            f"(KEGG meta-classification maps with no KOs/reactions, e.g. ko010**, ko07***)"
        )
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

    log.info(
        f"Built kegg_data dict: {len(out['kos'])} KOs, {len(out['pathways'])} pathways, "
        f"{len(out['reactions'])} reactions, {len(out['compounds'])} compounds"
    )
    return out


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

def main(force: bool = False, refetch_raw: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cache_root = PROJECT_ROOT / "cache" / "data"

    # `refetch_raw` implies rebuild (raw inputs changed → outputs are stale).
    force = force or refetch_raw

    if KEGG_DATA_FILE.exists() and not force:
        log.info(f"{KEGG_DATA_FILE} exists; use --force to rebuild.")
        return

    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run `bash scripts/refresh_mnx.sh` first, "
            "or set MNX_DATA_DIR to a checkout that already has the resolver."
        )

    log.info("Ensuring KEGG raw cache (downloads from KEGG REST if missing) ...")
    kegg_utils.download_kegg_raw(cache_root, force=refetch_raw)

    log.info("Ensuring TCDB reference TSVs are downloaded ...")
    download_all(cache_root=cache_root, force=refetch_raw, sources=["tcdb"])

    log.info("Building tcdb_hierarchy.json ...")
    _build_tcdb_hierarchy(cache_root)

    log.info("Parsing raw KEGG into in-memory dict ...")
    raw = _parse_raw_into_dict(cache_root)

    log.info("Opening MNX resolver ...")
    with contextlib.closing(mu.open_resolver(RESOLVER_DB)) as conn:
        log.info("Computing gene-reachable sets (catalysis side) ...")
        catalysis_sets = _gene_reachable_sets(raw)
        catalysis_cpds = catalysis_sets["cpds"]
        catalysis_pws = catalysis_sets["pws"]

        log.info("Pruning TCDB hierarchy + resolving transport substrates ...")
        hierarchy = json.loads(
            (cache_root / "tcdb" / "tcdb_hierarchy.json").read_text()
        )
        kept, leaf_subs = _prune_tcdb(hierarchy, catalysis_sets["tcdb_ids"])
        leaf_to_primary, compound_props = _resolve_substrates(leaf_subs, conn)
        total_substrate_ids = sum(len(ids) for ids in leaf_to_primary.values())
        distinct_substrate_ids = len({pid for ids in leaf_to_primary.values() for pid in ids})

        # Split substrate primaries: kegg.compound:C* feed back into the regular
        # _bulk_enrich_compounds pipeline (so they get pathways + MNX cross-refs);
        # non-KEGG primaries (chebi:*, mnx:*) become additional_compounds entries.
        substrate_kegg_cpds: set[str] = set()
        substrate_non_kegg_count = 0
        for primary_ids in leaf_to_primary.values():
            for primary_id in primary_ids:
                if primary_id.startswith("kegg.compound:"):
                    substrate_kegg_cpds.add(primary_id[len("kegg.compound:"):])
                else:
                    substrate_non_kegg_count += 1

        # Transport-reachable pathways: gene-annotated TCDB → MNX-resolved KEGG
        # compound → KEGG-curated pathway list. Adds pathways like glycolysis
        # when sucrose is a substrate even if no gene catalyzes it locally.
        cpd_to_pw = raw.get("compound_to_pathways", {})
        transport_pws: set[str] = {
            p for cpd in substrate_kegg_cpds for p in cpd_to_pw.get(cpd, [])
        }

        # ── Phase 2 metabolomics: harvest paper-measured metabolites ─────────
        log.info("Harvesting paper-measured metabolites from paperconfigs ...")
        paperconfigs = load_all_paperconfigs()
        paper_harvest = _harvest_paper_metabolites(paperconfigs, conn)
        paper_kegg_cpds: set[str] = paper_harvest["paper_kegg_cpds"]
        paper_non_kegg_compounds: dict[str, dict] = paper_harvest["paper_non_kegg_compounds"]

        # Paper-measurement-reachable pathways (KEGG-primary measured compounds only)
        paper_pws: set[str] = {
            p for cpd in paper_kegg_cpds for p in cpd_to_pw.get(cpd, [])
        }

        n_paper_total = sum(s["total"] for s in paper_harvest["per_paper"].values())
        n_paper_resolved = sum(s["resolved"] for s in paper_harvest["per_paper"].values())
        n_paper_unresolved = sum(s["unresolved"] for s in paper_harvest["per_paper"].values())
        log.info(
            f"  Paper-metabolite harvest: {n_paper_resolved}/{n_paper_total} resolved "
            f"({n_paper_unresolved} unresolved); "
            f"{len(paper_kegg_cpds)} kegg.compound, "
            f"{len(paper_non_kegg_compounds)} non-KEGG primaries"
        )

        extended_sets = dict(catalysis_sets)
        extended_sets["cpds"] = catalysis_cpds | substrate_kegg_cpds | paper_kegg_cpds
        extended_sets["pws"] = catalysis_pws | transport_pws | paper_pws

        log.info(
            f"TCDB substrate primaries: {distinct_substrate_ids} distinct "
            f"({len(substrate_kegg_cpds)} kegg.compound, "
            f"{distinct_substrate_ids - len(substrate_kegg_cpds)} non-KEGG)"
        )
        log.info(
            f"Extended sets: {len(extended_sets['cpds'])} compounds "
            f"(+{len(substrate_kegg_cpds - catalysis_cpds)} transport-reachable, "
            f"+{len(paper_kegg_cpds - catalysis_cpds - substrate_kegg_cpds)} paper-measurement-reachable), "
            f"{len(extended_sets['pws'])} pathways "
            f"(+{len(transport_pws - catalysis_pws)} transport-reachable, "
            f"+{len(paper_pws - catalysis_pws - transport_pws)} paper-measurement-reachable)"
        )

        log.info("Building pruned kegg_data ...")
        kegg_data = build_pruned_kegg_data(raw, conn, sets=extended_sets)

        log.info("Tagging evidence_sources + adding non-KEGG transport substrates ...")
        _fold_substrates_into_kegg_data(
            kegg_data,
            catalysis_cpds=catalysis_cpds,
            substrate_kegg_cpds=substrate_kegg_cpds,
            leaf_to_primary=leaf_to_primary,
            compound_props=compound_props,
        )

        # Phase 2: union "metabolomics" into evidence_sources for paper compounds
        log.info("Folding paper-measured metabolites into kegg_data ...")
        _fold_paper_metabolites_into_kegg_data(
            kegg_data,
            paper_kegg_cpds=paper_kegg_cpds,
            paper_non_kegg_compounds=paper_non_kegg_compounds,
        )

        # Write metabolite_id_mapping.json (consumed by step 7)
        mapping_out = PROJECT_ROOT / "cache" / "data" / "metabolomics" / "metabolite_id_mapping.json"
        _write_metabolite_id_mapping(paper_harvest, mapping_out)
        log.info(f"  metabolite_id_mapping.json written: {mapping_out}")
        n_compounds = len(kegg_data["compounds"])
        n_additional = len(kegg_data["additional_compounds"])
        n_overlap = sum(
            1 for c in kegg_data["compounds"].values()
            if "metabolism" in c["evidence_sources"]
            and "transport" in c["evidence_sources"]
        )
        n_transport_only_kegg = sum(
            1 for c in kegg_data["compounds"].values()
            if c["evidence_sources"] == ["transport"]
        )
        KEGG_DATA_FILE.parent.mkdir(parents=True, exist_ok=True)
        KEGG_DATA_FILE.write_text(json.dumps(kegg_data, indent=2, sort_keys=True))

        (cache_root / "tcdb" / "tcdb_pruned.json").write_text(json.dumps({
            "kept_tcdb_ids": sorted(kept),
            "leaf_substrates": {
                leaf: ids for leaf, ids in sorted(leaf_to_primary.items())
            },
        }, indent=2, sort_keys=True))
        log.info(
            f"  TCDB: {len(kept)} kept IDs, "
            f"{len(leaf_to_primary)} leaves with substrates "
            f"({total_substrate_ids} total / {distinct_substrate_ids} distinct primary IDs)"
        )
        log.info(
            f"  kegg_data.json: {n_compounds} compounds "
            f"({n_overlap} metabolism+transport, "
            f"{n_transport_only_kegg} transport-only kegg), "
            f"{n_additional} additional_compounds (transport-only non-KEGG)"
        )
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild kegg_data.json from cached raw inputs "
                             "(does NOT re-download raw KEGG/TCDB).")
    parser.add_argument("--refetch-raw", action="store_true",
                        help="Also re-download raw KEGG REST files + TCDB TSVs. "
                             "Implies --force. Use when KEGG/TCDB releases.")
    args = parser.parse_args()
    main(force=args.force, refetch_raw=args.refetch_raw)
