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

# TCDB's curated cross-references. Both source URLs embed the data in HTML, so
# `download_tcdb_xref_map` strips non-data lines before caching.
PFAM_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/pfam.py"
PFAM_MAP_PATH = Path("cache/data/tcdb/raw/tcdb_pfam_map.tsv")
GO_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/go.py"
GO_MAP_PATH = Path("cache/data/tcdb/raw/tcdb_go_map.tsv")

# Bridge edges are only emitted onto tc_family (level 2) or deeper. A 5-part TCID
# whose nearest surviving ancestor is a tc_class / tc_subclass would produce an
# edge like "this domain relates to Channels and Pores" — true but useless.
MIN_BRIDGE_LEVEL = 2


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


def _valid_pfam(token: str) -> bool:
    return token.startswith("PF") and len(token) >= 7 and token[2:7].isdigit()


def _valid_go(token: str) -> bool:
    return token.startswith("GO:") and token[3:].isdigit()


_XREF_KINDS = {
    "pfam": (PFAM_MAP_URL, PFAM_MAP_PATH, _valid_pfam),
    "go": (GO_MAP_URL, GO_MAP_PATH, _valid_go),
}


def download_tcdb_xref_map(kind: str, dest: Path | None = None, force: bool = False) -> Path:
    """Download a TCDB cross-reference map to a clean 3-column TSV.

    ``kind`` is ``"pfam"`` or ``"go"``. Writes ``<xref_id> \\t <TC_id> \\t <name>``.
    Skipped when the file exists and ``force`` is False.
    """
    if kind not in _XREF_KINDS:
        raise ValueError(f"unknown TCDB xref kind {kind!r}; expected one of {sorted(_XREF_KINDS)}")
    url, default_path, is_valid = _XREF_KINDS[kind]
    out = Path(dest) if dest is not None else default_path
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not force:
        return out

    with urllib.request.urlopen(url) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    with open(out, "w", encoding="utf-8") as fh:
        for line in text.splitlines():
            parts = line.rstrip().split("\t")
            if len(parts) < 3 or not is_valid(parts[0]) or parts[1].count(".") < 2:
                continue
            fh.write("\t".join(parts[:3]) + "\n")
    return out


def load_tcdb_xref_map(kind: str, path: Path | None = None) -> dict[str, set[str]]:
    """Load a TCDB cross-reference map into ``{xref_id: {full_tcid, ...}}``.

    TCIDs are kept at their published (5-part) depth — the roll-up to a surviving
    node happens in `build_tcdb_bridges`, which needs the original depth to record
    provenance. Returns an empty dict when the file is absent.
    """
    if kind not in _XREF_KINDS:
        raise ValueError(f"unknown TCDB xref kind {kind!r}; expected one of {sorted(_XREF_KINDS)}")
    _, default_path, is_valid = _XREF_KINDS[kind]
    p = Path(path) if path is not None else default_path
    if not p.exists():
        return {}
    result: dict[str, set[str]] = defaultdict(set)
    with open(p, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not is_valid(parts[0]) or parts[1].count(".") < 2:
                continue
            result[parts[0]].add(parts[1])
    return dict(result)


def build_tcdb_bridges(
    xref_to_tcids: dict[str, set[str]],
    hierarchy: dict[str, dict],
    kept: set[str],
    min_level: int = MIN_BRIDGE_LEVEL,
) -> dict[str, dict[str, list[str]]]:
    """Roll a published xref→TCID map onto the pruned node set.

    TCDB publishes these maps at 5-part `tc_specificity` depth, but the pruned
    hierarchy keeps only gene-annotated TCIDs + ancestors — so most published
    TCIDs have no node. Each is walked up to its nearest **surviving** ancestor,
    mirroring how `subtree_substrates` rolls leaf substrates onto kept ancestors.
    Targets shallower than ``min_level`` are dropped as uninformative.

    Returns ``{tc_id: {xref_id: [original_tcids, ...]}}`` — the original
    pre-roll-up TCIDs are retained so the edge can record exactly what TCDB
    curated, and roll-up precision is not silently lost.

    **Semantics of the result**: an entry means *"TCDB's curated reference
    proteins for this transport family carry this domain / this GO annotation"*.
    That is a statement about the composition of the TC family. It reads
    correctly outward from the TcdbFamily; it does **not** license the reverse
    inference that anything carrying the xref belongs to the family.
    """
    bridges: dict[str, dict[str, list[str]]] = {}
    for xref_id, tcids in xref_to_tcids.items():
        for tcid in tcids:
            anchor: str | None = tcid if tcid in kept else None
            if anchor is None:
                cur = hierarchy.get(tcid, {}).get("parent")
                while cur is not None:
                    if cur in kept:
                        anchor = cur
                        break
                    cur = hierarchy.get(cur, {}).get("parent")
            if anchor is None:
                continue
            if hierarchy.get(anchor, {}).get("level", 0) < min_level:
                continue
            bridges.setdefault(anchor, {}).setdefault(xref_id, []).append(tcid)
    return {
        tc: {x: sorted(set(v)) for x, v in sorted(xrefs.items())}
        for tc, xrefs in sorted(bridges.items())
    }
