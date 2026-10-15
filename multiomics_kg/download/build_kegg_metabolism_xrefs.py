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

Step 6 also downloads the 4 TCDB reference TSVs (families, superfamilies,
substrates, acc2tcid) and writes the assembled hierarchy to
cache/data/tcdb/tcdb_hierarchy.json. The TCDB lift is unconditional and not
gated on KEGG reachability. Node names come from the committed
cache/data/tcdb/tcdb_names.json, built separately by
`python -m multiomics_kg.download.scrape_tcdb_names` (T6) — step 6 itself
never scrapes.

TCDB substrate rollup runs in three phases:

  rollup → prune → resolve

  1. `_compute_subtree_substrate_strings` walks the FULL hierarchy bottom-up
     and unions raw `CHEBI:NNNN;name` strings from every tc_specificity
     descendant per TC ID. No MNX dependency.
  2. `_prune_tcdb` filters the hierarchy to the gene-reachable subhierarchy
     (above + below each gene-annotated seed).
  3. The rollup is filtered to kept TC IDs, the surviving distinct substrate
     strings are resolved through MNX once via `_resolve_substrate_strings`,
     and per-kept primary lists are mapped back. The adapter emits one
     Tcdb_family_transports_metabolite edge per (kept_tc_id, primary_id) pair.

Doing the rollup pre-prune ensures every kept ancestor carries substrates from
ALL its TCDB descendants, not just the gene-annotated subset (an earlier
post-prune rollup undercounted ancestor substrates by ~950 edges).

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
from collections import defaultdict
from dataclasses import dataclass, replace
from datetime import date
from pathlib import Path

import pandas as pd

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils import kegg_utils
from multiomics_kg.utils import metabolite_utils as mu
from multiomics_kg.utils.chebi_utils import (
    ChebiIndex,
    download_chebi_obo,
    index_from_cache,
    parse_chebi_obo,
)
from multiomics_kg.utils.gene_id_utils import load_gene_annotations
from multiomics_kg.utils.metabolite_utils import resolve_metabolite_candidates
from multiomics_kg.utils.metabolite_identity.canonical import load_identity_overrides
from multiomics_kg.utils.metabolite_identity.classes import COMPOUND_BRITE_TREES, download_compound_brite
from multiomics_kg.utils.metabolite_identity.outputs import write_all_atomic
from multiomics_kg.utils.metabolite_identity.paper import Override, PaperTable, harvest_paper_tables
from multiomics_kg.utils.metabolite_identity.pipeline import (
    Step6Inputs, Step6Paths, render_outputs, run_identity_stages)
from multiomics_kg.utils.metabolite_identity.roles import (
    build_role_hierarchy, build_subtree_roles, partition_role_strings, role_node_id)
from multiomics_kg.utils.paperconfig_utils import iter_metabolite_assays_tables, load_all_paperconfigs
from multiomics_kg.download.download_metabolism_reference import download_all

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
RESOLVER_DB = mu.get_mnx_data_dir() / "metabolite_resolver.db"

CHEBI_DIR = PROJECT_ROOT / "cache" / "data" / "chebi"
CHEBI_CACHE_FILE = CHEBI_DIR / "chebi_metabolite_xrefs.json"
METABOLOMICS_DIR = PROJECT_ROOT / "cache" / "data" / "metabolomics"
LOCK_FILE = METABOLOMICS_DIR / "dominant_stereo_forms.json"
IDENTITY_REPORT_FILE = METABOLOMICS_DIR / "metabolite_identity_report.json"
MAPPING_FILE = METABOLOMICS_DIR / "metabolite_id_mapping.json"
IDENTITY_OVERRIDES_FILE = PROJECT_ROOT / "config" / "metabolite_identity_overrides.yaml"

log = logging.getLogger(__name__)


# ── Stage A: ChEBI ───────────────────────────────────────────────────────────

def _load_chebi(cache_root: Path, refetch_raw: bool) -> tuple[ChebiIndex, bool]:
    """Stage A (spec §5.1, §6.4): rebuild from the raw obo whenever one is available — freshly
    downloaded, or already on disk when the download fails; reuse the committed derived cache
    only when no usable raw obo exists (TCDB outage precedent).

    Returns ``(index, from_raw)``.
    """
    raw_dir = cache_root / "chebi" / "raw"
    failure: Exception | None = None
    try:
        obo = download_chebi_obo(raw_dir, force=refetch_raw)
    except Exception as exc:
        failure = exc
        obo = raw_dir / "chebi.obo.gz"
        log.warning(f"  ChEBI download failed ({exc}); "
                    f"{'trying the raw obo on disk' if obo.exists() else 'no raw obo on disk'}")
    if failure is None or obo.exists():
        try:
            release, terms = parse_chebi_obo(obo)
            log.info(f"  ChEBI release {release}: {len(terms)} terms (from raw obo)")
            return ChebiIndex(terms, release), True
        except Exception as exc:
            failure = exc
    if not CHEBI_CACHE_FILE.exists():
        raise failure
    log.warning(f"  no usable raw ChEBI obo ({failure}); reusing the committed "
                f"{CHEBI_CACHE_FILE.name}. Ids it lacks are reported, and step 6 fails "
                f"if a canonical node needs one.")
    return index_from_cache(json.loads(CHEBI_CACHE_FILE.read_text())), False


# ── TCDB hierarchy ────────────────────────────────────────────────────────────

def _load_tcdb_names(names_path: Path | None) -> tuple[dict[str, dict], list[str]]:
    """Read the committed tcdb_names.json (built by scrape_tcdb_names — T6).

    Returns ({tc_id: {"name": str[, "description": str]}}, scraped_families).
    Missing file → empty maps (deep levels stay unnamed; the caller warns via
    _warn_unscraped_families).
    """
    if names_path is None or not names_path.exists():
        return {}, []
    data = json.loads(names_path.read_text())
    return data.get("names", {}), data.get("meta", {}).get("scraped_families", [])


def _warn_unscraped_families(kept_ids: set[str], names_path: Path) -> list[str]:
    """T6 self-healing guard: kept 4/5-part TC IDs whose family the name
    scraper has not covered — loud staleness after a strain onboarding, never
    silent. Returns the sorted missing family list."""
    needed = {".".join(tc.split(".")[:3]) for tc in kept_ids
              if len(tc.split(".")) >= 4}
    if not names_path.exists():
        if needed:
            log.warning(
                f"  cache/data/tcdb/tcdb_names.json is MISSING — {len(needed)} "
                f"families' subfamily/specificity nodes will render as bare TC "
                f"ids. Run: uv run python -m multiomics_kg.download.scrape_tcdb_names"
            )
        return sorted(needed)
    _names, scraped = _load_tcdb_names(names_path)
    missing = sorted(needed - set(scraped))
    if missing:
        log.warning(
            f"  {len(missing)} kept TCDB families not covered by tcdb_names.json "
            f"(e.g. {missing[:5]}) — their deep nodes render as bare TC ids. "
            f"Top-up: uv run python -m multiomics_kg.download.scrape_tcdb_names "
            f"--families {' '.join(missing[:5])} ..."
        )
    return missing


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


def _parse_tcdb_acc2tcid(path: Path) -> set[str]:
    """Parse acc2tcid.tsv → set of distinct 5-part TCIDs (column 1).

    acc2tcid.py is the most complete TCDB classification table (~19.5K distinct
    5-part TCIDs vs ~11.5K in superfamilies.tsv); using it as a hierarchy seed
    closes the gap for gene-annotated TCIDs that aren't curated into a
    superfamily.
    """
    out: set[str] = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1]:
                out.add(parts[1])
    return out


def build_tcdb_hierarchy(
    out_path: Path,
    families_path: Path,
    superfamilies_path: Path,
    substrates_path: Path,
    acc2tcid_path: Path,
    names_path: Path | None = None,
) -> int:
    """Build tcdb_hierarchy.json by joining the 4 TCDB sources.

    Seeds the hierarchy from the union of:
      - 3-part family IDs in families.tsv (unconditional — every documented
        family lands in the hierarchy regardless of whether it has a curated
        5-part child).
      - 5-part specificity IDs in {superfamilies.tsv col0, substrates.tsv,
        acc2tcid.tsv}. acc2tcid is the broadest catalog and previously dropped
        ~3K TCIDs that genes are annotated against.

    Class/subclass/family/subfamily ancestors are synthesized for every 5-part
    seed. When a TCID is also present in superfamilies.tsv, the family
    abbreviation and superfamily-name are pulled from there; otherwise those
    fields are blank.
    """
    fam_descs = _parse_tcdb_families(families_path)
    names, _scraped = _load_tcdb_names(names_path)
    super_rows = _parse_tcdb_superfamilies(superfamilies_path)
    substrates = _parse_tcdb_substrates(substrates_path)
    acc_tcids = _parse_tcdb_acc2tcid(acc2tcid_path)

    # Index superfamily metadata for lookup when synthesizing from other sources
    super_by_tcid: dict[str, tuple[str, str]] = {}      # 5-part → (abbr, superfam)
    fam_meta_by_id: dict[str, tuple[str, str]] = {}     # 3-part → (abbr, superfam)
    for row in super_rows:
        tcid, fam, abbr, superfam = row[0], row[2], row[3], row[4]
        super_by_tcid[tcid] = (abbr, superfam)
        fam_meta_by_id.setdefault(fam, (abbr, superfam))

    h: dict[str, dict] = {}

    def ensure_class(cls: str) -> None:
        if cls not in h:
            h[cls] = {
                "name": names.get(cls, {}).get("name", ""),
                "level": 0,
                "level_kind": "tc_class",
                "parent": None,
            }

    def ensure_subclass(subcls: str, cls: str) -> None:
        if subcls not in h:
            h[subcls] = {
                "name": names.get(subcls, {}).get("name", ""),
                "level": 1,
                "level_kind": "tc_subclass",
                "parent": cls,
            }

    def ensure_family(fam: str, subcls: str) -> None:
        if fam not in h:
            abbr, superfam = fam_meta_by_id.get(fam, ("", ""))
            entry: dict = {
                "name": fam_descs.get(fam, "") or names.get(fam, {}).get("name", ""),
                "level": 2,
                "level_kind": "tc_family",
                "parent": subcls,
            }
            if abbr:
                entry["abbreviation"] = abbr
            if superfam:
                entry["superfamily"] = superfam
            h[fam] = entry

    def ensure_subfamily(subfam: str, fam: str) -> None:
        if subfam not in h:
            h[subfam] = {
                "name": names.get(subfam, {}).get("name", ""),
                "level": 3,
                "level_kind": "tc_subfamily",
                "parent": fam,
            }

    def ensure_specificity(tcid: str, subfam: str) -> None:
        if tcid in h:
            return
        node: dict = {
            "name": names.get(tcid, {}).get("name", ""),
            "level": 4,
            "level_kind": "tc_specificity",
            "parent": subfam,
        }
        desc = names.get(tcid, {}).get("description", "")
        if desc:
            node["description"] = desc
        if tcid in substrates:
            node["substrate_classes"] = substrates[tcid]
        superfam = super_by_tcid.get(tcid, ("", ""))[1]
        if superfam:
            node["superfamily"] = superfam
        h[tcid] = node

    # Seed every documented family (3-part) — unconditional, so family-level
    # gene annotations like 1.A.11 always land in the hierarchy even without
    # curated 5-part children.
    for fam in fam_descs:
        parts = fam.split(".")
        if len(parts) != 3:
            continue
        ensure_class(parts[0])
        ensure_subclass(".".join(parts[:2]), parts[0])
        ensure_family(fam, ".".join(parts[:2]))

    # Seed every 5-part TCID we know about, synthesizing missing ancestors.
    five_part_seeds = (
        {row[0] for row in super_rows}
        | set(substrates.keys())
        | acc_tcids
    )
    for tcid in five_part_seeds:
        parts = tcid.split(".")
        if len(parts) != 5:
            continue
        cls = parts[0]
        subcls = ".".join(parts[:2])
        fam = ".".join(parts[:3])
        subfam = ".".join(parts[:4])
        ensure_class(cls)
        ensure_subclass(subcls, cls)
        ensure_family(fam, subcls)
        ensure_subfamily(subfam, fam)
        ensure_specificity(tcid, subfam)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(h, indent=2, sort_keys=True))
    log.info(f"  tcdb_hierarchy.json: {len(h)} entries")
    return len(h)


def _build_tcdb_hierarchy(cache_root: Path) -> int:
    """Parse the 4 TCDB TSVs into cache_root/tcdb/tcdb_hierarchy.json. Returns entry count.

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
        acc2tcid_path=raw_dir / "acc2tcid.tsv",
        names_path=tcdb_dir / "tcdb_names.json",
    )


# ── Reachability ──────────────────────────────────────────────────────────────

def _gene_reachable_sets(raw: dict) -> dict:
    """Walk all strains to compute the gene-reachable {KOs, reactions, pathways, compounds, TCDB IDs}.

    Returns a dict with keys 'kos', 'rxns', 'cpds', 'pws', 'tcdb_ids' (sets) plus
    'rxn_genes' ({reaction: {"{data_dir}:{locus_tag}", ...}}, the per-reaction gene keys
    stage-D catalyst counts are taken from).
    Pathway set is KO ∪ Rxn-reachable (Option B).
    """
    kos: set[str] = set()
    rxns: set[str] = set()
    tcdb_ids: set[str] = set()
    rxn_genes: dict[str, set[str]] = {}
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for gid, g in genes.items():
            for ko in g.get("kegg_ko") or []:
                if isinstance(ko, str) and ko.startswith("K"):
                    kos.add(ko)
            for rxn in g.get("kegg_reactions") or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    rxns.add(rxn)
                    rxn_genes.setdefault(rxn, set()).add(f"{row['data_dir']}:{gid}")
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
    return {"kos": kos, "rxns": rxns, "cpds": cpds, "pws": pws, "tcdb_ids": tcdb_ids,
            "rxn_genes": rxn_genes}


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
) -> tuple[set[str], dict[str, str]]:
    """Ancestor-only prune of the TCDB hierarchy.

    For each seed ID, walk **up** to the tc_class root. Returns:

      kept: set of all TCDB IDs that survive pruning.
      seed_aliases: {seed_id: anchor_id} mapping for seeds NOT in the curated
                    hierarchy. Resolved by truncating dot-segments until an
                    ancestor is found (e.g. ``3.A.1.35`` → ``3.A.1``). Callers
                    should remap the gene→TCDB attachment from `seed_id` to
                    `anchor_id` so no annotation is dropped.

    Seeds with no resolvable ancestor are dropped and logged as a warning at
    the call site.

    **No downward walk.** An earlier version also walked *down* to every
    tc_specificity leaf under each seed, so 528 gene-annotated TC IDs expanded to
    12,902 nodes — 94.5% of which had `gene_count = 0` (10,292 of them
    tc_specificity). That arm existed so substrate information on the leaves stayed
    reachable, but it duplicates work `_compute_subtree_substrate_strings` already
    does: the rollup is computed over the FULL hierarchy pre-pruning, so every kept
    ancestor already carries the union of its descendants' substrates. Keeping the
    leaves as well was a second mechanism for the same goal, and it is what forced
    the `is_promiscuous` warning flag onto consumers.

    Tradeoff: `level_kind = 'tc_specificity'` now selects the specificity nodes genes
    actually annotate, rather than every specificity node in a reachable family.

    Substrate rollup is intentionally NOT done here. The full subtree-substrate
    map is computed pre-pruning (`_compute_subtree_substrate_strings`); the caller
    filters that map to `kept` and resolves the surviving substrate strings via
    MNX. This separates structural pruning from semantic enrichment.
    """
    kept: set[str] = set()
    parent_of = {tc: data.get("parent") for tc, data in hierarchy.items()}

    def walk_up(tc: str) -> None:
        cur: str | None = tc
        while cur is not None and cur in hierarchy:
            if cur in kept:
                return
            kept.add(cur)
            cur = parent_of.get(cur)

    def find_ancestor(tc: str) -> str | None:
        """Truncate dot-segments until an ancestor is found in the hierarchy."""
        parts = tc.split(".")
        for n in range(len(parts) - 1, 0, -1):
            candidate = ".".join(parts[:n])
            if candidate in hierarchy:
                return candidate
        return None

    seed_aliases: dict[str, str] = {}
    for seed in seed_ids:
        if seed in hierarchy:
            walk_up(seed)
            continue
        anchor = find_ancestor(seed)
        if anchor is None:
            # No ancestor anywhere — caller logs a warning.
            seed_aliases[seed] = ""
            continue
        seed_aliases[seed] = anchor
        walk_up(anchor)
    return kept, seed_aliases


def _compute_subtree_substrate_strings(hierarchy: dict) -> dict[str, set[str]]:
    """Bottom-up DFS over the FULL hierarchy: for each TC ID, return the union
    of `substrate_classes` strings (raw 'CHEBI:NNNN;name') from every
    tc_specificity descendant (including self if the node is itself a leaf).

    Pure structural rollup — no MNX dependency. Caller filters the result to
    kept TC IDs after pruning, then resolves the (smaller) string union via
    `_resolve_substrate_strings`. Doing the rollup pre-prune ensures ancestor
    nodes capture substrates from ALL their TCDB descendants, not just the
    gene-annotated subset.
    """
    parent_of = {tc: e.get("parent") for tc, e in hierarchy.items()}
    children_of: dict[str, list[str]] = {}
    for tc, parent in parent_of.items():
        if parent is not None:
            children_of.setdefault(parent, []).append(tc)

    cache: dict[str, set[str]] = {}

    def dfs(node: str) -> set[str]:
        if node in cache:
            return cache[node]
        entry = hierarchy.get(node, {})
        result: set[str] = set()
        if entry.get("level_kind") == "tc_specificity":
            subs = entry.get("substrate_classes")
            if subs:
                result.update(subs)
        for child in children_of.get(node, []):
            result |= dfs(child)
        cache[node] = result
        return result

    return {tc: dfs(tc) for tc in hierarchy}


def _resolve_substrate_strings(
    substrate_strings: set[str],
    conn: sqlite3.Connection,
    chebi: ChebiIndex | None = None,
    rekeys: dict[str, str] | None = None,
    known_kegg: set[str] | None = None,
    legacy: dict[str, str] | None = None,
) -> tuple[dict[str, str], dict[str, dict]]:
    """Resolve raw 'CHEBI:NNNN;name' substrate strings to Metabolite primary IDs.

    Returns:
      string_to_primary: {substrate_str: primary_id} — only resolvable strings.
      compound_props: {primary_id: {name, formula, mass, inchikey, mnxm_id, chebi_id}}
                      for non-KEGG primaries only (KEGG primaries flow through
                      `_bulk_enrich_compounds`, so they don't need props here).

    Operates on a flat set of substrate strings (post-pruning union). MNX
    queries are batched. Iteration over `parsed` is sorted so the first-seen
    `name_part` wins for non-KEGG primaries (deterministic).

    With ``chebi`` (stage D, spec §5.2): a secondary ChEBI id is mapped to its primary
    before the MNX lookup (which tries both), and when MNX carries no KEGG alias for the
    hit, the primary's ChEBI KEGG xref supplies one (lowest C-number) — that re-key is
    recorded in ``rekeys`` as ``chebi:<primary n>`` → ``kegg.compound:C…``. A non-KEGG
    primary is then ``chebi:<primary n>``. With ``chebi=None`` the behaviour is unchanged.

    ``known_kegg`` (the raw KEGG compound list): a KEGG id becomes a node id only if KEGG
    still lists it; otherwise the next alias / xref is tried (``None`` = no filter).

    ``legacy`` (a dict, filled in place): each resolvable string's primary under the
    ``chebi=None, known_kegg=None`` rules — what the pre-identity build assigned — computed
    from the same MNX queries, so the remap baseline costs no second resolution pass.
    """
    def _listed(kid: str) -> bool:
        return known_kegg is None or kid in known_kegg

    cur = conn.cursor()
    if rekeys is None:
        rekeys = {}

    # Pass 1: parse strings → (string, chebi_value, name_part). Sort for
    # deterministic name-fallback resolution.
    parsed: list[tuple[str, str, str]] = []
    for sub_str in sorted(substrate_strings):
        if ":" not in sub_str:
            continue
        chebi_part, _, name_part = sub_str.partition(";")
        prefix, _, chebi_value = chebi_part.partition(":")
        if prefix.strip().upper() != "CHEBI" or not chebi_value:
            log.debug(f"TCDB substrate not in CHEBI form: {sub_str!r}")
            continue
        parsed.append((sub_str, chebi_value.strip(), name_part))

    if not parsed:
        return {}, {}

    # Pass 1b (ChEBI only): original value → primary ChEBI number (None if ChEBI lacks it).
    primary_num: dict[str, str | None] = {}
    if chebi is not None:
        for _, chebi_value, _ in parsed:
            if chebi_value not in primary_num:
                p = chebi.primary(chebi_value)
                primary_num[chebi_value] = p.split(":", 1)[1] if p else None

    # Pass 2: chebi_value → mnxm_id (with ChEBI: query both the original and the primary
    # number; the lowest mnxm_id hit for either wins).
    distinct_chebi = sorted({chebi_v for _, chebi_v, _ in parsed}
                            | {n for n in primary_num.values() if n})
    value_to_mnxm: dict[str, str] = {}
    for chunk in _batched(distinct_chebi):
        placeholders = ",".join("?" * len(chunk))
        cur.execute(
            f"SELECT value, mnxm_id FROM compound_aliases "
            f"WHERE source='chebi' AND value IN ({placeholders}) "
            f"ORDER BY mnxm_id",
            chunk,
        )
        for value, mnxm_id in cur.fetchall():
            value_to_mnxm.setdefault(value, mnxm_id)
    if chebi is None:
        chebi_to_mnxm = value_to_mnxm
    else:
        chebi_to_mnxm = {}
        for _, chebi_value, _ in parsed:
            hits = [m for m in (value_to_mnxm.get(chebi_value),
                                value_to_mnxm.get(primary_num.get(chebi_value) or ""))
                    if m]
            if hits:
                chebi_to_mnxm[chebi_value] = min(hits)

    # Pass 3: distinct mnxm_ids → kegg.compound + chebi aliases (legacy: the original value's
    # own hit, which may differ from the primary-number hit chosen above)
    legacy_mnxm = {cv: value_to_mnxm[cv] for _, cv, _ in parsed if cv in value_to_mnxm}
    mnxm_ids = sorted(set(chebi_to_mnxm.values())
                      | (set(legacy_mnxm.values()) if legacy is not None else set()))
    mnxm_to_kegg: dict[str, str] = {}
    mnxm_to_kegg_any: dict[str, str] = {}       # unfiltered by known_kegg (legacy rule)
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
            mnxm_to_kegg_any.setdefault(mnxm_id, value)
            if _listed(value):
                mnxm_to_kegg.setdefault(mnxm_id, value)
        cur.execute(
            f"SELECT mnxm_id, value FROM compound_aliases "
            f"WHERE source='chebi' AND mnxm_id IN ({placeholders}) "
            f"ORDER BY mnxm_id, value",
            chunk,
        )
        for mnxm_id, value in cur.fetchall():
            mnxm_to_chebi.setdefault(mnxm_id, value)

    def _primary_for(mnxm: str, pnum: str | None = None) -> str:
        if mnxm in mnxm_to_kegg:
            return f"kegg.compound:{mnxm_to_kegg[mnxm]}"
        if pnum is not None:
            kegg = sorted(k for k in (chebi.term(pnum) or {}).get("kegg", []) if _listed(k))
            if kegg:
                rekeys[f"chebi:{pnum}"] = f"kegg.compound:{kegg[0]}"
                return f"kegg.compound:{kegg[0]}"
            return f"chebi:{pnum}"
        if mnxm in mnxm_to_chebi:
            return f"chebi:{mnxm_to_chebi[mnxm]}"
        return f"mnx:{mnxm}"

    # Pass 4: bulk-fetch compounds props only for non-KEGG mnxm_ids (a ChEBI re-key may
    # still turn one of these into a KEGG primary; its props are then simply unused).
    nonkegg_mnxms = sorted({m for m in chebi_to_mnxm.values() if m not in mnxm_to_kegg})
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

    # Pass 5: assemble outputs (sorted iteration → first-seen name_part wins)
    string_to_primary: dict[str, str] = {}
    compound_props: dict[str, dict] = {}
    for sub_str, chebi_value, name_part in parsed:
        if legacy is not None and chebi_value in legacy_mnxm:
            lm = legacy_mnxm[chebi_value]
            legacy[sub_str] = (f"kegg.compound:{mnxm_to_kegg_any[lm]}" if lm in mnxm_to_kegg_any
                               else f"chebi:{mnxm_to_chebi[lm]}" if lm in mnxm_to_chebi
                               else f"mnx:{lm}")
        mnxm = chebi_to_mnxm.get(chebi_value)
        if mnxm is None:
            log.debug(f"TCDB substrate unresolved: chebi:{chebi_value} ({sub_str!r})")
            continue
        pnum = primary_num.get(chebi_value) if chebi is not None else None
        primary = _primary_for(mnxm, pnum)
        string_to_primary[sub_str] = primary

        if primary.startswith("kegg.compound:") or primary in compound_props:
            continue
        mnx_name, formula, mass, inchikey = mnxm_to_props.get(mnxm, (None, None, None, None))
        compound_props[primary] = {
            "name": (mnx_name or name_part).strip(),
            "formula": formula or None,
            "mass": mass,
            "inchikey": inchikey or None,
            "mnxm_id": mnxm,
            "chebi_id": pnum if pnum is not None else mnxm_to_chebi.get(mnxm),
        }

    return string_to_primary, compound_props


# ── Stage D: paper candidates (spec §5.2 D) ──────────────────────────────────

def _mnxm_primary(mnxm: str, conn: sqlite3.Connection, chebi: ChebiIndex,
                  rekeys: dict[str, str], known_kegg: set[str] | None = None) -> str:
    """Node id for an MNXM: its MNX KEGG alias, else the KEGG xref of one of its ChEBI
    aliases' primaries (neutral first, then lowest number; recorded in ``rekeys``), else
    that ChEBI primary, else ``mnx:<MNXM>``. With ``known_kegg`` only KEGG ids KEGG still
    lists qualify (an unlisted one falls through to the next option)."""
    for (value,) in conn.execute("SELECT value FROM compound_aliases WHERE mnxm_id=? "
                                 "AND source='kegg.compound' ORDER BY value", (mnxm,)):
        if known_kegg is None or value in known_kegg:
            return f"kegg.compound:{value}"
    prims = sorted({p for (v,) in conn.execute(
        "SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='chebi'", (mnxm,))
        if (p := chebi.primary(v))}, key=lambda c: (chebi.charge(c) != 0, int(c.split(":")[1])))
    for p in prims:
        kegg = sorted(k for k in (chebi.term(p) or {}).get("kegg", [])
                      if known_kegg is None or k in known_kegg)
        if kegg:
            rekeys[f"chebi:{p.split(':')[1]}"] = f"kegg.compound:{kegg[0]}"
            return f"kegg.compound:{kegg[0]}"
    if prims:
        return f"chebi:{prims[0].split(':')[1]}"
    return f"mnx:{mnxm}"


def _normalize_paper_id(curie: str, conn: sqlite3.Connection, chebi: ChebiIndex,
                        rekeys: dict[str, str], known_kegg: set[str] | None = None) -> str:
    """Normalize a paper-supplied CURIE to a node id: upper-case KEGG C-numbers, ChEBI
    secondary → primary (→ KEGG via its xref, recorded in ``rekeys``), MNX via
    `_mnxm_primary`. Prefixes match case-insensitively (``CHEBI:``, ``KEGG:``, ``MNX:`` …) and
    normalise to the house forms; unknown prefixes pass through unchanged.

    ``known_kegg`` filters only the ChEBI-xref / MNX routes onto KEGG. A KEGG id the paper,
    an override or the lock names is kept as-is even when KEGG no longer lists it, so stage P
    reports it (it is not a compound step 6 writes) instead of it being silently replaced."""
    prefix, _, local = curie.strip().partition(":")
    prefix, local = prefix.strip().lower(), local.strip()
    if prefix in ("kegg", "kegg.compound"):
        return f"kegg.compound:{local.strip().upper()}"
    if prefix == "chebi":
        p = chebi.primary(local)
        if p is None:
            return f"chebi:{int(local)}" if local.isdigit() else curie
        kegg = sorted(k for k in (chebi.term(p) or {}).get("kegg", [])
                      if known_kegg is None or k in known_kegg)
        if kegg:
            rekeys[f"chebi:{p.split(':')[1]}"] = f"kegg.compound:{kegg[0]}"
            return f"kegg.compound:{kegg[0]}"
        return f"chebi:{p.split(':')[1]}"
    if prefix == "mnx":
        return _mnxm_primary(local, conn, chebi, rekeys, known_kegg)
    return curie


@dataclass
class PaperCandidates:
    name_candidates: dict[str, set[str]]
    ids: set[str]
    non_kegg_props: dict[str, dict]


def _harvest_paper_candidates(
    tables: list[PaperTable],
    overrides: dict[str, dict[str, Override]],
    lock: dict,
    conn: sqlite3.Connection,
    chebi: ChebiIndex,
    rekeys: dict[str, str],
    known_kegg: set[str] | None = None,
) -> PaperCandidates:
    """Stage D, paper side (spec §5.2 D): collect CANDIDATES, never decisions. Every paper
    row's MNX name candidates, direct id, the paper's override targets and the lock targets
    are all enriched like any compound; stage J decides, stage K drops the losers.

    Rewrites every row's ``direct_id`` and every override's ``id`` IN PLACE to its normalised
    form (``_normalize_paper_id``), so stages E and J look up exactly the id harvested here —
    ``CHEBI:17234`` / ``KEGG:C00031`` / ``c00031`` included."""
    name_candidates: dict[str, set[str]] = defaultdict(set)
    ids: set[str] = set()
    for t in tables:
        for row in t.rows:
            for mnxm in resolve_metabolite_candidates(row.name, conn):
                name_candidates[row.name_key].add(_mnxm_primary(mnxm, conn, chebi, rekeys, known_kegg))
            if row.direct_id:
                row.direct_id = _normalize_paper_id(row.direct_id, conn, chebi, rekeys, known_kegg)
                ids.add(row.direct_id)
    for ovs in overrides.values():
        for k, ov in ovs.items():
            ovs[k] = replace(ov, id=_normalize_paper_id(ov.id, conn, chebi, rekeys, known_kegg))
            ids.add(ovs[k].id)
    for section in ("names", "direct_ids"):
        for key, entry in lock.get(section, {}).items():
            ids.add(_normalize_paper_id(entry["id"], conn, chebi, rekeys, known_kegg))
            if section == "direct_ids":
                ids.add(_normalize_paper_id(key, conn, chebi, rekeys, known_kegg))
    for cands in name_candidates.values():
        ids |= cands
    non_kegg: dict[str, dict] = {}
    for pid in sorted(i for i in ids if not i.startswith("kegg.compound:")):
        if pid.startswith("mnx:"):
            mnxm = pid[4:]
        else:
            row = conn.execute("SELECT mnxm_id FROM compound_aliases WHERE source='chebi' AND value=? "
                               "ORDER BY mnxm_id LIMIT 1", (pid.split(":")[1],)).fetchone()
            mnxm = row[0] if row else None
        props = conn.execute("SELECT name, formula, mass, inchikey FROM compounds WHERE mnxm_id=?",
                             (mnxm,)).fetchone() if mnxm else None
        name, formula, mass, inchikey = props if props else (None, None, None, None)
        cterm = chebi.term(pid) if pid.startswith("chebi:") else None
        non_kegg[pid] = {"name": (cterm or {}).get("name") or name, "formula": formula or None,
                         "mass": mass, "inchikey": inchikey or None, "mnxm_id": mnxm,
                         "chebi_id": pid.split(":")[1] if pid.startswith("chebi:") else None}
    return PaperCandidates(dict(name_candidates), ids, non_kegg)


def _add_paper_candidates(kegg_data: dict, non_kegg_props: dict[str, dict]) -> None:
    """Insert missing non-KEGG paper candidates into ``additional_compounds`` with
    ``evidence_sources: []``. Never tags ``metabolomics`` — stage K does that for winners."""
    additional = kegg_data.setdefault("additional_compounds", {})
    for pid, props in non_kegg_props.items():
        additional.setdefault(pid, {**props, "evidence_sources": []})


def _attach_mnx_chebi_aliases(table: dict[str, dict], conn: sqlite3.Connection) -> None:
    """Set the transient ``mnx_chebi_aliases`` field (every ChEBI alias MNX lists for the
    entry's MNXM, shortest-then-lexicographic) on each table entry with an ``mnxm_id``."""
    by_mnx: dict[str, list[str]] = defaultdict(list)
    for node, e in table.items():
        if e.get("mnxm_id"):
            by_mnx[e["mnxm_id"]].append(node)
    found: dict[str, set[str]] = defaultdict(set)
    for chunk in _batched(sorted(by_mnx)):
        q = ",".join("?" * len(chunk))
        for mnxm, value in conn.execute(
                f"SELECT mnxm_id, value FROM compound_aliases WHERE source='chebi' AND mnxm_id IN ({q})",
                chunk):
            found[mnxm].add(value)
    for mnxm, nodes in by_mnx.items():
        for n in nodes:
            table[n]["mnx_chebi_aliases"] = sorted(found.get(mnxm, ()), key=lambda v: (len(v), v))


def _to_table(kegg_data: dict) -> dict[str, dict]:
    """Flatten ``compounds`` + ``additional_compounds`` into one dict keyed by full node id
    (shallow copies — mutating an entry never touches ``kegg_data``)."""
    table = {f"kegg.compound:{k}": dict(v) for k, v in kegg_data["compounds"].items()}
    table.update({k: dict(v) for k, v in kegg_data.get("additional_compounds", {}).items()})
    return table


def _fold_substrates_into_kegg_data(
    kegg_data: dict,
    *,
    catalysis_cpds: set[str],
    substrate_kegg_cpds: set[str],
    substrate_primaries: set[str],
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

    `substrate_primaries` is the flat union of every primary ID appearing in
    any kept TC ID's subtree rollup.
    """
    for cpd_id, entry in kegg_data.get("compounds", {}).items():
        sources = list(entry.get("evidence_sources") or [])
        if cpd_id in catalysis_cpds and "metabolism" not in sources:
            sources.append("metabolism")
        if cpd_id in substrate_kegg_cpds and "transport" not in sources:
            sources.append("transport")
        entry["evidence_sources"] = sources

    additional: dict[str, dict] = {}
    for primary_id in substrate_primaries:
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


def _tcdb_bridges(cache_root: Path, hierarchy: dict, kept: set[str], previous: dict,
                  refetch_raw: bool) -> dict[str, dict]:
    """Stage N: TCDB Pfam→TC / GO→TC bridges rolled onto the kept nodes. Returns
    ``{"pfam": {...}, "go": {...}}``; ``previous`` is the committed tcdb_pruned.json (outage
    fallback). Raises ``RuntimeError`` rather than lose a committed bridge — before any
    stage-P write, so a failure writes nothing."""
    # Cross-ontology bridges. TCDB publishes curated Pfam→TC and GO→TC maps at
    # 5-part specificity depth; roll each onto the nearest surviving node
    # (same shape as subtree_substrates). Semantics: "this TC family's curated
    # reference proteins carry this domain / GO annotation" — a composition
    # statement about the family, NOT a rule for classifying the xref.
    from multiomics_kg.utils.tcdb_utils import (
        build_tcdb_bridges,
        download_tcdb_xref_map,
        load_tcdb_xref_map,
    )
    # tcdb_pruned.json is COMMITTED, so a silent degradation here is durable:
    # an empty bridge writes cleanly, drops ~8.4K edges from the next build,
    # and shows up only as a log line nobody reads. tcdb.org is a flaky host,
    # so fail loudly on both failure modes instead.
    tcdb_bridges: dict[str, dict] = {}
    for kind in ("pfam", "go"):
        dest = cache_root / "tcdb" / "raw" / f"tcdb_{kind}_map.tsv"
        try:
            download_tcdb_xref_map(kind, dest=dest, force=refetch_raw)
        except Exception as exc:
            # A cached TSV is a fine fallback. With no TSV at all we fall
            # through to the reuse-or-raise branch below, which can still
            # recover the bridge from the committed tcdb_pruned.json —
            # tcdb.org 404s its CGI endpoints during maintenance while the
            # site root stays up, and that outage must not block step 6.
            log.warning(
                f"  TCDB {kind}→TC map download failed ({exc}); "
                f"{'using cached ' + str(dest) if dest.exists() else 'no cached copy'}"
            )
        xref = load_tcdb_xref_map(kind, path=dest)
        prior = previous.get(f"{kind}_bridge", {})
        prior_pairs = sum(len(v) for v in prior.values())
        if xref:
            tcdb_bridges[kind] = build_tcdb_bridges(xref, hierarchy, kept)
            n_pairs = sum(len(v) for v in tcdb_bridges[kind].values())
            # Never trade a populated committed bridge for an empty one.
            if n_pairs == 0 and prior_pairs > 0:
                raise RuntimeError(
                    f"TCDB {kind}→TC bridge came out EMPTY but the committed "
                    f"tcdb_pruned.json holds {prior_pairs} pairs. That is data loss, not a "
                    f"rebuild — check {dest} (size "
                    f"{dest.stat().st_size if dest.exists() else 'missing'}). Delete the "
                    f"stale TSV and re-run with --refetch-raw if the map genuinely changed."
                )
        elif prior_pairs:
            # No usable map (TCDB down, nothing cached) but the committed
            # artifact still holds the bridge — reuse it rather than drop
            # every bridge edge. Filtered to THIS run's kept set, so stale
            # nodes cannot leak back in.
            tcdb_bridges[kind] = {k: v for k, v in prior.items() if k in kept}
            n_pairs = sum(len(v) for v in tcdb_bridges[kind].values())
            log.warning(
                f"  TCDB {kind}→TC map unavailable; reused {n_pairs} pairs on "
                f"{len(tcdb_bridges[kind])} kept node(s) from the committed "
                f"tcdb_pruned.json. Re-run with --refetch-raw once tcdb.org is back."
            )
        else:
            raise RuntimeError(
                f"TCDB {kind}→TC map is unavailable at {dest} and the committed "
                f"tcdb_pruned.json holds no {kind} bridge to fall back on. Refusing to "
                f"write tcdb_pruned.json without it — it is a committed artifact and the "
                f"loss would be silent. Retry when tcdb.org is reachable."
            )
        log.info(
            f"  TCDB {kind}→TC bridge: {len(xref)} published {kind} ids -> "
            f"{n_pairs} pairs on {len(tcdb_bridges[kind])} kept TC nodes"
        )
    return tcdb_bridges


def _load_compound_brite(cache_root: Path, refetch_raw: bool) -> tuple[dict, dict]:
    """Stage M inputs: the raw KEGG compound BRITE trees (downloaded, or the raw cache), and the
    committed pruned ``brite_<code>.json`` files the pipeline falls back on per tree."""
    raw = download_compound_brite(cache_root / "kegg" / "raw", force=refetch_raw)   # per-tree fallback
    committed = {code: json.loads(p.read_text()) for code in COMPOUND_BRITE_TREES
                 if (p := KEGG_CACHE_DIR / f"brite_{code}.json").exists()}
    return raw, committed


def _read_old_resolved(tables: list[PaperTable]) -> dict[tuple[str, str, int], str]:
    """The previous run's paper resolutions: ``(paper, table, row_index)`` → ``metabolite_id``
    from each existing ``<stem>_resolved.csv``. Those files are written without the source's
    preamble, so they are read with no ``skiprows`` — row ``i`` is source data row ``i``."""
    out: dict[tuple[str, str, int], str] = {}
    for t in tables:
        if t.csv_path is None:
            continue
        path = t.csv_path.with_name(t.csv_path.stem + "_resolved.csv")
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "metabolite_id" not in df.columns:
            continue
        for i, mid in enumerate(df["metabolite_id"]):
            if mid.strip():
                out[(t.paper, t.table, i)] = mid.strip()
    return out


def _missing_resolved_tables() -> list[Path]:
    """``<stem>_resolved.csv`` files step 6 owes but has not written — a newly added
    ``metabolite_assays_table``. Their presence is what lets a default (no ``--force``) run
    skip; the retired step 7 used to write them on that basis."""
    missing = []
    for _, cfg in load_all_paperconfigs():
        for _, entry in iter_metabolite_assays_tables(cfg):
            csv = entry.get("filename")
            if not csv:
                continue
            path = Path(csv) if Path(csv).is_absolute() else PROJECT_ROOT / csv
            resolved = path.with_name(path.stem + "_resolved.csv")
            if not resolved.exists():
                missing.append(resolved)
    return sorted(missing)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(force: bool = False, refetch_raw: bool = False, today: str | None = None) -> None:
    """Prepare-data step 6 (spec 2026-10-13 §5.2): stages A-D here, then
    ``run_identity_stages`` (E-M, O, P) and one atomic write of every output.

    The one write before stage P: ``cache/data/tcdb/tcdb_hierarchy.json`` is rebuilt from the
    raw TCDB TSVs (``cache/data/tcdb/raw/``) before anything else runs, outside the atomic
    write. It is deterministic given those TSVs (a run that changes nothing else rewrites it
    byte-identically); on a TCDB download / rebuild failure the committed file is reused.
    Everything else is written only after stage P validates, so a failed run writes no other
    file.
    """
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cache_root = PROJECT_ROOT / "cache" / "data"

    # `refetch_raw` implies rebuild (raw inputs changed → outputs are stale).
    force = force or refetch_raw

    if KEGG_DATA_FILE.exists() and not force:
        missing = _missing_resolved_tables()
        if not missing:
            log.info(f"{KEGG_DATA_FILE} exists; use --force to rebuild.")
            return
        log.info(f"{len(missing)} metabolite_assays_table(s) have no _resolved.csv "
                 f"(e.g. {missing[0]}); rebuilding step 6 to write them.")

    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run `bash scripts/refresh_mnx.sh` first, "
            "or set MNX_DATA_DIR to a checkout that already has the resolver."
        )

    log.info("Ensuring KEGG raw cache (downloads from KEGG REST if missing) ...")
    kegg_utils.download_kegg_raw(cache_root, force=refetch_raw)

    # TCDB's CGI endpoints go down for maintenance from time to time (they 404
    # while www.tcdb.org itself still answers). The raw TSVs are gitignored but
    # tcdb_hierarchy.json is committed, so a download failure is recoverable:
    # reuse the committed hierarchy and carry on. Only a missing hierarchy is
    # fatal.
    tcdb_hierarchy_path = cache_root / "tcdb" / "tcdb_hierarchy.json"
    try:
        log.info("Ensuring TCDB reference TSVs are downloaded ...")
        download_all(cache_root=cache_root, force=refetch_raw, sources=["tcdb"])
        log.info("Building tcdb_hierarchy.json ...")
        _build_tcdb_hierarchy(cache_root)
    except Exception as exc:
        if not tcdb_hierarchy_path.exists():
            raise
        log.warning(
            f"  TCDB reference download/rebuild failed ({exc}); reusing the "
            f"committed {tcdb_hierarchy_path.name} unchanged."
        )

    log.info("Parsing raw KEGG into in-memory dict ...")
    raw = _parse_raw_into_dict(cache_root)

    log.info("Opening MNX resolver ...")
    with contextlib.closing(mu.open_resolver(RESOLVER_DB)) as conn:
        log.info("Computing gene-reachable sets (catalysis side) ...")
        catalysis_sets = _gene_reachable_sets(raw)
        catalysis_cpds = catalysis_sets["cpds"]
        catalysis_pws = catalysis_sets["pws"]

        log.info("Loading TCDB hierarchy + computing full subtree-substrate rollup (no MNX) ...")
        hierarchy = json.loads(
            (cache_root / "tcdb" / "tcdb_hierarchy.json").read_text()
        )
        full_subtree_strings = _compute_subtree_substrate_strings(hierarchy)

        log.info("Pruning TCDB hierarchy ...")
        kept, seed_aliases = _prune_tcdb(hierarchy, catalysis_sets["tcdb_ids"])
        _warn_unscraped_families(kept, cache_root / "tcdb" / "tcdb_names.json")
        upgraded = {s: a for s, a in seed_aliases.items() if a}
        dropped = sorted(s for s, a in seed_aliases.items() if not a)
        if upgraded:
            sample = ", ".join(f"{s}→{a}" for s, a in sorted(upgraded.items())[:5])
            log.info(
                f"  TCDB: {len(upgraded)} gene-annotated seed(s) re-anchored to "
                f"nearest curated ancestor (e.g. {sample})"
            )
        if dropped:
            sample = ", ".join(dropped[:5])
            log.warning(
                f"  TCDB: {len(dropped)} gene-annotated seed(s) had no ancestor "
                f"in the curated hierarchy and were dropped (e.g. {sample}). "
                f"Likely retired TCDB IDs in eggNOG annotations."
            )

        # Filter rollup to kept TC IDs; collect distinct substrate strings to resolve
        subtree_strings: dict[str, set[str]] = {tc: full_subtree_strings.get(tc, set()) for tc in kept}
        all_substrate_strings: set[str] = set()
        for s in subtree_strings.values():
            all_substrate_strings |= s

        # A — ChEBI
        chebi, chebi_from_raw = _load_chebi(cache_root, refetch_raw)
        today = today or date.today().isoformat()

        # previous outputs: old ids for the remap, the frozen lock and pin, old paper resolutions
        prev_kegg = json.loads(KEGG_DATA_FILE.read_text()) if KEGG_DATA_FILE.exists() else {}
        previous_ids = ({f"kegg.compound:{k}" for k in prev_kegg.get("compounds", {})}
                        | set(prev_kegg.get("additional_compounds", {})))
        prev_pruned_path = cache_root / "tcdb" / "tcdb_pruned.json"
        prev_pruned: dict = {}
        if prev_pruned_path.exists():
            try:
                prev_pruned = json.loads(prev_pruned_path.read_text())
            except (OSError, json.JSONDecodeError):
                prev_pruned = {}
        lock = json.loads(LOCK_FILE.read_text()) if LOCK_FILE.exists() else {"names": {}, "direct_ids": {}}

        # C — roles out of the substrate strings, before MNX
        molecule_strings, role_map = partition_role_strings(all_substrate_strings, chebi)
        chemical_roles = build_role_hierarchy(set(role_map.values()), chebi)
        subtree_roles = build_subtree_roles(subtree_strings, role_map)
        role_of_old = {}
        with chebi.untracked():          # baseline-only lookups never shape the ChEBI cache
            for old in prev_kegg.get("additional_compounds", {}):
                if old.startswith("chebi:") and chebi.is_role(old):
                    role_of_old[old] = role_node_id(chebi.primary(old))
        log.info(f"  roles: {len(set(role_map.values()))} TCDB roles → {len(chemical_roles)} ChemicalRole nodes")

        # D — resolve molecules (secondary→primary, ChEBI KEGG re-key), harvest paper candidates
        log.info(f"Resolving {len(molecule_strings)} distinct TCDB molecule substrate strings via MNX ...")
        rekeys: dict[str, str] = {}
        # A KEGG id becomes a node id only while KEGG still lists it (no nameless nodes).
        known_kegg = set(raw.get("compound_names", {}))
        if not known_kegg:
            # An empty list would silently push every KEGG route onto ChEBI / MNX ids.
            raise RuntimeError("raw KEGG compound list is empty (list/compound parse or download "
                               "failed) — refusing to build metabolite ids without it")
        # old_string_to_primary: the previous build's string → primary map (pre-identity rules),
        # so an old id the other remap rules miss follows its string. Only molecule strings
        # matter — the remap consults it only where string_to_primary also has the string.
        old_string_to_primary: dict[str, str] = {}
        string_to_primary, compound_props = _resolve_substrate_strings(
            molecule_strings, conn, chebi=chebi, rekeys=rekeys, known_kegg=known_kegg,
            legacy=old_string_to_primary)
        subtree_substrates: dict[str, list[str]] = {}
        for tc, strings in subtree_strings.items():
            primaries = {string_to_primary[s] for s in strings if s in string_to_primary}
            if primaries:
                subtree_substrates[tc] = sorted(primaries)
        distinct_substrate_primaries = {p for ids in subtree_substrates.values() for p in ids}
        substrate_kegg_cpds = {p[len("kegg.compound:"):] for p in distinct_substrate_primaries
                               if p.startswith("kegg.compound:")}
        cpd_to_pw = raw.get("compound_to_pathways", {})
        transport_pws = {p for cpd in substrate_kegg_cpds for p in cpd_to_pw.get(cpd, [])}

        tables, paper_overrides = harvest_paper_tables(load_all_paperconfigs(), PROJECT_ROOT)
        paper = _harvest_paper_candidates(tables, paper_overrides, lock, conn, chebi, rekeys,
                                          known_kegg=known_kegg)
        # A paper / override / lock KEGG id KEGG no longer lists is not built (stage P then
        # reports its rows' target as not a compound step 6 holds; curate an override).
        paper_kegg_cpds = {i[len("kegg.compound:"):] for i in paper.ids
                           if i.startswith("kegg.compound:") and i[len("kegg.compound:"):] in known_kegg}
        paper_pws = {p for cpd in paper_kegg_cpds for p in cpd_to_pw.get(cpd, [])}

        extended_sets = dict(catalysis_sets)
        extended_sets["cpds"] = catalysis_cpds | substrate_kegg_cpds | paper_kegg_cpds
        extended_sets["pws"] = catalysis_pws | transport_pws | paper_pws
        log.info("Building pruned kegg_data ...")
        kegg_data = build_pruned_kegg_data(raw, conn, sets=extended_sets)
        _fold_substrates_into_kegg_data(kegg_data, catalysis_cpds=catalysis_cpds,
                                        substrate_kegg_cpds=substrate_kegg_cpds,
                                        substrate_primaries=distinct_substrate_primaries,
                                        compound_props=compound_props)
        _add_paper_candidates(kegg_data, paper.non_kegg_props)
        table = _to_table(kegg_data)
        _attach_mnx_chebi_aliases(table, conn)

        # N — TCDB Pfam / GO bridges (raises before any write rather than lose a bridge)
        tcdb_bridges = _tcdb_bridges(cache_root, hierarchy, kept, prev_pruned, refetch_raw)

        # M inputs — KEGG compound BRITE (raw, else committed on outage)
        brite_raw, brite_committed = _load_compound_brite(cache_root, refetch_raw)

        old_resolved = _read_old_resolved(tables)

    inp = Step6Inputs(
        table=table,
        kegg_rest={k: kegg_data[k] for k in ("kos", "pathways", "subcategories", "categories", "reactions")},
        subtree_substrates=subtree_substrates, chemical_roles=chemical_roles, subtree_roles=subtree_roles,
        rekeys=rekeys, role_of_old=role_of_old, tables=tables, paper_overrides=paper_overrides,
        name_candidates=paper.name_candidates, lock=lock, lumping_pin=prev_pruned.get("lumping_pin"),
        identity_overrides=load_identity_overrides(IDENTITY_OVERRIDES_FILE), chebi=chebi,
        chebi_from_raw=chebi_from_raw, rxn_genes=catalysis_sets["rxn_genes"],
        compound_to_pathways=cpd_to_pw, base_pathways=catalysis_pws | transport_pws,
        brite_raw=brite_raw, brite_committed=brite_committed, tcdb_kept=kept,
        tcdb_levels={t: int(hierarchy.get(t, {}).get("level", 0)) for t in kept},
        tcdb_static={"kept_tcdb_ids": sorted(kept), "pfam_bridge": tcdb_bridges["pfam"],
                     "go_bridge": tcdb_bridges["go"], "seed_aliases": dict(sorted(upgraded.items()))},
        previous_ids=previous_ids, old_resolved=old_resolved, today=today,
        previous_entries={**{f"kegg.compound:{k}": v for k, v in prev_kegg.get("compounds", {}).items()},
                          **prev_kegg.get("additional_compounds", {})},
        old_string_to_primary=old_string_to_primary, string_to_primary=string_to_primary)
    out = run_identity_stages(inp)               # raises Step6ValidationError → nothing written
    files = render_outputs(out, Step6Paths(
        kegg_data=KEGG_DATA_FILE, tcdb_pruned=prev_pruned_path, brite_dir=KEGG_CACHE_DIR,
        chebi_cache=CHEBI_CACHE_FILE, mapping=MAPPING_FILE, lock=LOCK_FILE, report=IDENTITY_REPORT_FILE))
    write_all_atomic(files)
    c = out.report["counts"]
    log.info(f"  step 6 wrote {len(files)} files: {c['metabolites']} metabolites "
             f"({c['merged_away']} merged away, {c['rekeyed']} re-keyed, {c['dropped']} dropped), "
             f"{c['roles']} ChemicalRole, {c['class_nodes']} MetaboliteClass")
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
