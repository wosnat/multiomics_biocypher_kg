"""TCDB hierarchy accessor API + curated Pfam→TC map.

Reads cache/data/tcdb/tcdb_hierarchy.json (built by `multiomics_kg.download.build_kegg_metabolism_xrefs` in prepare_data step 6).
"""
from __future__ import annotations

import json
import urllib.request
from collections import defaultdict
from pathlib import Path

DEFAULT_PATH = Path("cache/data/tcdb/tcdb_hierarchy.json")
_CACHE: dict[str, dict] | None = None

# TCDB's curated Pfam → TC association table. The source URL embeds the data in
# HTML, so `download_pfam_to_tc_map` strips non-data lines before caching.
PFAM_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/pfam.py"
PFAM_MAP_PATH = Path("cache/data/tcdb/raw/tcdb_pfam_map.tsv")


def load_tcdb() -> dict[str, dict]:
    """Load the TCDB hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_tcdb(tc_id: str) -> bool:
    """Return True if `tc_id` exists as a key in the hierarchy."""
    if not tc_id:
        return False
    return tc_id in load_tcdb()


def tcdb_ancestors(tc_id: str) -> list[str]:
    """Return root-to-parent ancestor chain for `tc_id`. Empty for unknown / root."""
    h = load_tcdb()
    if tc_id not in h:
        return []
    chain: list[str] = []
    parent = h[tc_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))


# ============================================================================
# Curated Pfam → TC map (backs the Pfam_in_tcdb_family bridge edge)
# ============================================================================
#
# TCDB publishes a curated table of which Pfam domains are associated with which
# transporter families. It is the one piece of diamond↔eggNOG cross-checking that
# is NOT derivable from the graph, so it is modelled as a bridge edge
# (Pfam -> TcdbFamily), mirroring Pfam_in_interpro_entry.
#
# Previously downloaded + consumed inside the /tcdb-diamond runner, which made the
# runner depend on downstream annotation state. It now belongs to prepare_data
# step 6 alongside the other TCDB reference TSVs.


def download_pfam_to_tc_map(dest: Path | None = None, force: bool = False) -> Path:
    """Download TCDB's curated Pfam→TC map to a clean 3-column TSV.

    Writes ``PF_id \\t TC_id \\t family_name``. Skipped when the file exists and
    ``force`` is False. Returns the destination path; on download failure the
    path is returned unwritten so the caller can decide whether to proceed.
    """
    out = Path(dest) if dest is not None else PFAM_MAP_PATH
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not force:
        return out

    with urllib.request.urlopen(PFAM_MAP_URL) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    with open(out, "w", encoding="utf-8") as fh:
        for line in text.splitlines():
            parts = line.rstrip().split("\t")
            if len(parts) < 3 or not parts[0].startswith("PF"):
                continue
            if not (len(parts[0]) >= 7 and parts[0][2:7].isdigit()):
                continue
            if parts[1].count(".") < 2:
                continue
            fh.write("\t".join(parts[:3]) + "\n")
    return out


def load_pfam_to_tc_map(path: Path | None = None) -> dict[str, set[str]]:
    """Load the curated Pfam→TC mapping into ``{pfam_id: {tc_family, ...}}``.

    TC IDs are truncated to 3-part families — the granularity at which TCDB
    curates Pfam ↔ TC associations. Returns an empty dict when the file is absent.
    """
    p = Path(path) if path is not None else PFAM_MAP_PATH
    if not p.exists():
        return {}
    result: dict[str, set[str]] = defaultdict(set)
    with open(p, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pfam_id, tc_id = parts[0], parts[1]
            if not pfam_id.startswith("PF") or tc_id.count(".") < 2:
                continue
            result[pfam_id].add(".".join(tc_id.split(".")[:3]))
    return dict(result)
