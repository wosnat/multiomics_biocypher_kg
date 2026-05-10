# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein post-steps.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path


def classify_hit(hit: dict) -> int | None:
    """Assign a confidence tier (1/2/3) to a parsed diamond hit row.

    Returns None when the hit fails the tier-3 floor (drop it).

    Floor (gblast3-style): e-value <= 0.001, HSP length >= 50, AND
    (qcov >= 40 OR scov >= 40). Above the floor:
      - tier 1: identity >= 70 AND qcov >= 70
      - tier 2: identity >= 40 AND qcov >= 60
      - tier 3: floor only
    """
    if hit["evalue"] > 0.001:
        return None
    if hit["length"] < 50:
        return None
    if hit["qcov"] < 40.0 and hit["scov"] < 40.0:
        return None

    if hit["identity"] >= 70.0 and hit["qcov"] >= 70.0:
        return 1
    if hit["identity"] >= 40.0 and hit["qcov"] >= 60.0:
        return 2
    return 3


def truncate_tcid(tcid: str, n_parts: int) -> str:
    """Truncate a TCID to its first n_parts dot-separated segments.

    A TCID has up to 5 parts (class.subclass.family.subfamily.specificity).
    Returns the input unchanged if it already has <= n_parts segments.
    """
    if not tcid:
        return ""
    parts = tcid.split(".")
    return ".".join(parts[:n_parts])


def consensus_collapse(hits: list[dict]) -> dict | None:
    """Collapse a list of hits (one query, top-N TCDB hits) into a per-protein call.

    Returns None when fewer than 3 leading parts agree across all hits.
    The returned dict has keys: tcid (str), agreement ("5_part" | "4_part" | "3_part"),
    n (number of hits considered).

    All hits are assumed to carry 5-part TCIDs (TCDB curates only at the
    tc_specificity leaves), so dot-prefix comparison is well-defined.
    """
    if not hits:
        return None
    parts_lists = [h["tcid"].split(".") for h in hits]

    for depth in (5, 4, 3):
        prefixes = {tuple(p[:depth]) for p in parts_lists}
        if len(prefixes) == 1:
            shared = parts_lists[0][:depth]
            return {
                "tcid": ".".join(shared),
                "agreement": f"{depth}_part",
                "n": len(hits),
            }
    return None


def compute_egn_agreement(diamond_tcid: str, egn_tcid: str | None) -> str:
    """Tag the relationship between the diamond call and eggNOG's KEGG_TC.

    Returns one of: "confirms" | "refines" | "extends" | "conflicts".
    `egn_only` (eggNOG TC present, diamond absent) is not produced here —
    those proteins simply don't appear in the calls JSON in Phase 1.

    Rules:
      - confirms: identical TCIDs OR one is a strict prefix of the other
        AT family level or below (first 3 parts match)
      - refines: eggNOG TCID is a strict prefix of diamond TCID — same lineage,
        diamond went deeper. Reported separately from confirms because this
        is the headline specificity win.
      - extends: eggNOG had no TC; diamond produced one
      - conflicts: family-level (first 3 parts) disagrees
    """
    if not egn_tcid:
        return "extends"
    if diamond_tcid == egn_tcid:
        return "confirms"

    diamond_parts = diamond_tcid.split(".")
    egn_parts = egn_tcid.split(".")

    # Family-level disagreement -> conflict
    if diamond_parts[:3] != egn_parts[:3]:
        return "conflicts"

    # Same family. Diamond strictly deeper than eggNOG -> refines.
    # eggNOG strictly deeper than diamond -> confirms (rare).
    if len(diamond_parts) > len(egn_parts) and diamond_parts[: len(egn_parts)] == egn_parts:
        return "refines"
    return "confirms"


def is_class_9(tcid: str) -> bool:
    """True iff TCID is in TCDB class 9 (Incompletely Characterized Transport Systems)."""
    if not tcid:
        return False
    return tcid.split(".", 1)[0] == "9"


def parse_diamond_row(line: str) -> dict | None:
    """Parse one diamond blastp output line (--outfmt 6, 8 columns) to a dict.

    Returns None when the row is malformed (wrong column count, non-numeric
    field). Caller is responsible for further filtering / classification.
    """
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 8:
        return None
    try:
        return {
            "query_id": parts[0],
            "subject_id": parts[1],
            "identity": float(parts[2]),
            "qcov": float(parts[3]),
            "scov": float(parts[4]),
            "length": int(parts[5]),
            "evalue": float(parts[6]),
            "bitscore": float(parts[7]),
        }
    except ValueError:
        return None


_TCID_TAIL_RE = re.compile(r"-(\d+(?:\.[A-Za-z0-9]+){2,4})$")


def parse_tcdb_subject_id(subject_id: str) -> tuple[str, str] | None:
    """Extract (accession, tcid) from a TCDB FASTA-derived subject ID.

    Header format: ``[lcl|]<accession>-<TCID>`` where TCID is dot-separated
    with 3-5 parts (e.g. ``lcl|Q9I3F6-1.A.11.1.5``).

    Returns None when the subject ID does not contain a parseable TCID tail.
    Splits on the LAST dash followed by a dotted TCID — handles UniProt
    isoform accessions (e.g. ``P12345-2``) correctly.
    """
    if not subject_id:
        return None
    if subject_id.startswith("lcl|"):
        subject_id = subject_id[4:]
    match = _TCID_TAIL_RE.search(subject_id)
    if not match:
        return None
    tcid = match.group(1)
    accession = subject_id[: match.start()]
    if not accession:
        return None
    return accession, tcid


def load_eggnog_kegg_tc(annotations_path: Path) -> dict[str, str]:
    """Read an eggNOG-mapper .emapper.annotations file -> {protein_id: kegg_tc}.

    Returns an empty dict when the file is absent. Skips rows whose KEGG_TC
    column is "-", empty, or missing. When KEGG_TC carries multiple values
    (rare; comma-separated), only the first is returned (we do not attempt
    to merge here — the merge happens in compute_egn_agreement at call time).
    """
    annotations_path = Path(annotations_path)
    if not annotations_path.exists():
        return {}

    result: dict[str, str] = {}
    KEGG_TC_COL = 17  # 0-indexed, post-emapper-v2.1 column order

    with open(annotations_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue  # the #query header line
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= KEGG_TC_COL:
                continue
            protein_id = cols[0]
            kegg_tc_raw = cols[KEGG_TC_COL].strip()
            if not kegg_tc_raw or kegg_tc_raw == "-":
                continue
            # Multi-value: take first
            kegg_tc = kegg_tc_raw.split(",")[0].strip()
            if kegg_tc:
                result[protein_id] = kegg_tc
    return result


_TIER_TO_LEVEL_KIND = {
    1: "tc_specificity",
    2: "tc_subfamily",
    3: "tc_family",
}
_AGREEMENT_TO_PARTS = {"5_part": 5, "4_part": 4, "3_part": 3}


def build_strain_calls(
    tsv_path: Path,
    eggnog_annotations_path: Path,
) -> tuple[dict, dict]:
    """Run the full per-strain pipeline: parse TSV, classify, consensus, tag.

    Returns (calls, summary):
      calls: dict keyed by protein_id with the full §6.5 record shape
      summary: dict with raw counts (matches §6.6 stdout columns)
    """
    egn_lookup = load_eggnog_kegg_tc(Path(eggnog_annotations_path))

    # Group accepted (tier-classified) hits per query
    by_query: dict[str, list[dict]] = defaultdict(list)
    raw_lines = 0
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            raw_lines += 1
            row = parse_diamond_row(line)
            if row is None:
                continue
            tier = classify_hit(row)
            if tier is None:
                continue
            parsed = parse_tcdb_subject_id(row["subject_id"])
            if parsed is None:
                continue
            _, hit_tcid = parsed
            by_query[row["query_id"]].append({
                "tcid": hit_tcid,
                "tier": tier,
                "identity": row["identity"],
                "qcov": row["qcov"],
                "scov": row["scov"],
                "evalue": row["evalue"],
                "length": row["length"],
            })

    calls: dict[str, dict] = {}
    rejected = 0
    for query_id, hits in by_query.items():
        # Consensus collapse
        consensus = consensus_collapse(hits)
        if consensus is None:
            rejected += 1
            continue

        # Consensus depth tells us how many TCID parts agree across hits.
        # The effective tier is the MORE CONSERVATIVE of:
        #   - depth-implied tier (5-part=1, 4-part=2, 3-part=3)
        #   - worst (most conservative) per-hit tier
        # This ensures a single tier-2 hit is never promoted to tier-1 level
        # even when it happens to hit a unique 5-part TCID.
        n_parts = _AGREEMENT_TO_PARTS[consensus["agreement"]]
        worst_tier = max(h["tier"] for h in hits)
        depth_tier = {5: 1, 4: 2, 3: 3}[n_parts]
        effective_tier = max(worst_tier, depth_tier)
        # Truncate TCID to the parts justified by effective_tier (not consensus depth)
        effective_n_parts = {1: 5, 2: 4, 3: 3}[effective_tier]
        called_tcid = truncate_tcid(consensus["tcid"], effective_n_parts)

        # Best (highest-identity) hit drives the metadata fields
        best = max(hits, key=lambda h: h["identity"])

        egn_tcid = egn_lookup.get(query_id)
        agreement = compute_egn_agreement(called_tcid, egn_tcid)

        # Class 9 (Incompletely Characterized) hits are capped at tier 3 —
        # identity may be high but the TCDB assignment itself is uncertain.
        # The TCID depth is kept (no truncation); only the confidence tier is lowered.
        if is_class_9(called_tcid):
            effective_tier = 3

        calls[query_id] = {
            "tcid": called_tcid,
            "level_kind": _TIER_TO_LEVEL_KIND[effective_tier],
            "tier": effective_tier,
            "identity": best["identity"],
            "qcov": best["qcov"],
            "scov": best["scov"],
            "evalue": best["evalue"],
            "length": best["length"],
            "consensus_n": consensus["n"],
            "consensus_agreement": consensus["agreement"],
            "egn_agreement": agreement,
            "egn_tcid": egn_tcid,
            "incompletely_characterized": is_class_9(called_tcid),
        }

    # Build summary
    tier_dist: dict[str, int] = defaultdict(int)
    agreement_dist: dict[str, int] = defaultdict(
        int, {"confirms": 0, "refines": 0, "extends": 0, "conflicts": 0}
    )
    for c in calls.values():
        tier_dist[str(c["tier"])] += 1
        agreement_dist[c["egn_agreement"]] += 1

    summary = {
        "raw_hit_lines": raw_lines,
        "proteins_with_hits": len(by_query),
        "proteins_with_call": len(calls),
        "proteins_rejected_by_consensus": rejected,
        "tier_distribution": dict(tier_dist),
        "agreement_distribution": dict(agreement_dist),
    }
    return calls, summary
