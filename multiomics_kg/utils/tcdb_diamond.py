# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein consensus collapse.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.

**This module produces PRIMARY SEQUENCE EVIDENCE ONLY.** It derives calls from
`protein.faa` vs. the curated TCDB diamond DB and nothing else. It deliberately
does NOT read eggNOG annotations or `gene_annotations_merged.json`.

Rationale (see `docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md`
§3.2): the previous version read `gene_annotations_merged.json` — the *output* of
prepare_data step 2 — to compute Pfam-agreement tags, while step 2 is exactly where
these calls are destined to be merged. That is a cycle: re-running step 2 mutated
calls.json even when the diamond blast was byte-identical.

Cross-source reconciliation (eggNOG agreement, Pfam corroboration) now happens
downstream where both sources are already present:

  - eggNOG agreement -> derivable from `Gene_has_tcdb_family.sources` + the
    `Tcdb_family_is_a_tcdb_family` hierarchy.
  - Pfam corroboration -> the `Pfam_in_tcdb_family` bridge edge, built from
    TCDB's curated Pfam->TC map in prepare_data step 6.

There is also no `filter_action`: the old post-hoc filter chain made a candidate's
verdict depend on whether the protein happened to have *unrelated* sibling
candidates (§3.1). The tier policy in `classify_hit` is the quality gate, and it is
sibling-independent by construction.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Iterator
import re


def classify_hit(hit: dict) -> int | None:
    """Assign a confidence tier (1/2/3) to a parsed diamond hit row.

    Returns None when the hit fails the tier-3 floor (drop it).

    Floor (gblast3-style): e-value <= 0.001, HSP length >= 50, AND
    (qcov >= 40 OR scov >= 40). Above the floor:
      - tier 1: identity >= 70 AND qcov >= 70
      - tier 2: identity >= 40 AND qcov >= 60
      - tier 3: floor only

    The tier is both a confidence band AND the depth the call is truncated to
    (see `_TIER_TO_LEVEL_KIND`): weaker similarity yields a deliberately
    broader claim. This coupling is what makes tier 3 conservative rather than
    overconfident.
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


_TIER_TO_LEVEL_KIND = {
    1: "tc_specificity",
    2: "tc_subfamily",
    3: "tc_family",
}

_AGREEMENT_TO_PARTS = {"5_part": 5, "4_part": 4, "3_part": 3}
# Down-weight scores when consensus is shallower — agreement at full depth
# is the strongest evidence; agreement only at family level is weakest.
_AGREEMENT_WEIGHT = {"5_part": 1.0, "4_part": 0.85, "3_part": 0.7}


def confidence_score(identity: float, qcov: float, agreement: str) -> float:
    """Continuous 0-1 confidence summary for a TCDB call.

    ``score = (identity / 100) * (qcov / 100) * agreement_weight``

    A complement to the discrete `tier` field: tier buckets calls into
    {1, 2, 3} for policy decisions; the score lets downstream consumers do
    their own thresholding without losing the underlying gradient.

    Inputs are taken from the BEST hit's identity + qcov (since the call's
    sequence-evidence is anchored to that hit), multiplied by the agreement
    weight (1.0 / 0.85 / 0.7 for 5/4/3-part consensus).
    """
    return (identity / 100.0) * (qcov / 100.0) * _AGREEMENT_WEIGHT[agreement]


def build_strain_calls(tsv_path: Path) -> tuple[dict, dict]:
    """Run the full per-strain pipeline: parse TSV, classify, consensus collapse.

    Takes ONLY the raw diamond TSV. No eggNOG, no merged annotations, no
    external map files — see the module docstring for why.

    Returns (calls, summary):
      calls:   {protein_id: {"calls": [candidate, ...]}}, candidates sorted by
               confidence_score descending
      summary: per-strain counts for the stdout status table

    Each candidate carries only sequence-derived evidence:
      tcid, level_kind, tier, confidence_score, identity, qcov, scov, evalue,
      length, consensus_n, consensus_agreement, incompletely_characterized
    """
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
    total_candidates = 0
    for query_id, hits in by_query.items():
        # Group hits by 3-part TC family. Each family produces one candidate;
        # multi-domain proteins (RND pump + MFP adaptor + OMF, etc.) emit
        # multiple candidates instead of being rejected by global consensus.
        by_family: dict[str, list[dict]] = defaultdict(list)
        for h in hits:
            family = ".".join(h["tcid"].split(".")[:3])
            by_family[family].append(h)

        candidates: list[dict] = []
        for _family, family_hits in by_family.items():
            # All hits in this group share at least 3 parts; consensus depth
            # within the group is whichever deeper prefix (4 or 5 parts) they
            # also share — and never returns None here since the 3-part floor
            # is already guaranteed.
            consensus = consensus_collapse(family_hits)
            if consensus is None:
                continue

            n_parts = _AGREEMENT_TO_PARTS[consensus["agreement"]]
            best_tier = min(h["tier"] for h in family_hits)
            depth_tier = {5: 1, 4: 2, 3: 3}[n_parts]
            effective_tier = max(best_tier, depth_tier)
            effective_n_parts = {1: 5, 2: 4, 3: 3}[effective_tier]
            called_tcid = truncate_tcid(consensus["tcid"], effective_n_parts)

            # Best (highest-identity) hit in THIS family drives the candidate's metadata
            best = max(family_hits, key=lambda h: h["identity"])

            candidates.append({
                "tcid": called_tcid,
                "level_kind": _TIER_TO_LEVEL_KIND[effective_tier],
                "tier": effective_tier,
                "confidence_score": round(
                    confidence_score(best["identity"], best["qcov"], consensus["agreement"]),
                    4,
                ),
                "identity": best["identity"],
                "qcov": best["qcov"],
                "scov": best["scov"],
                "evalue": best["evalue"],
                "length": best["length"],
                "consensus_n": consensus["n"],
                "consensus_agreement": consensus["agreement"],
                "incompletely_characterized": is_class_9(called_tcid),
            })

        if not candidates:
            continue

        candidates.sort(key=lambda c: -c["confidence_score"])
        calls[query_id] = {"calls": candidates}
        total_candidates += len(candidates)

    # Build summary. tier_distribution + consensus_agreement_distribution count
    # CANDIDATES (not proteins): a 2-family protein contributes 2 entries.
    # Use proteins_with_call vs total_candidates to gauge how many proteins
    # are multi-family.
    tier_dist: dict[str, int] = defaultdict(int)
    consensus_dist: dict[str, int] = defaultdict(
        int, {"5_part": 0, "4_part": 0, "3_part": 0}
    )
    candidates_per_protein: dict[str, int] = defaultdict(int)
    for rec in calls.values():
        candidates_per_protein[str(len(rec["calls"]))] += 1
        for cand in rec["calls"]:
            tier_dist[str(cand["tier"])] += 1
            consensus_dist[cand["consensus_agreement"]] += 1

    summary = {
        "raw_hit_lines": raw_lines,
        "proteins_with_hits": len(by_query),
        "proteins_with_call": len(calls),
        "total_candidates": total_candidates,
        "candidates_per_protein_distribution": dict(candidates_per_protein),
        "tier_distribution": dict(tier_dist),
        "consensus_agreement_distribution": dict(consensus_dist),
    }
    return calls, summary


# ============================================================================
# calls.json consumption helpers
# ============================================================================
#
# These utilities walk the per-strain `<strain>.tcdb.calls.json` artifact
# produced by `build_strain_calls`. The KG merge uses them to:
#   - load a strain's calls
#   - iterate every candidate (nothing is pre-filtered — see module docstring)
#   - pick the top-confidence candidate per protein
#   - extract per-protein TC family / TCID lists for the annotation merge


def load_calls_json(path: Path) -> dict[str, dict]:
    """Load a single strain's `<strain>.tcdb.calls.json` into a dict.

    Raises FileNotFoundError if the file is missing — callers should
    handle this explicitly (a strain without a calls.json is a real
    "no TCDB annotation available for this strain" case, not silent).
    Returns the raw on-disk shape: `{protein_id: {"calls": [...]}}`.
    """
    path = Path(path)
    return json.loads(path.read_text())


def iter_candidates(calls: dict[str, dict]) -> Iterator[tuple[str, dict]]:
    """Yield (protein_id, candidate_dict) for every candidate.

    Nothing is filtered — the tier policy already applied the quality gate at
    build time. Consumers that want a stricter cut should filter on `tier`,
    `identity`, or `confidence_score` explicitly, so the threshold is visible
    at the point of use rather than frozen into the artifact.

    Order within a protein matches `calls[]` ordering (descending by
    `confidence_score`).
    """
    for pid, rec in calls.items():
        for cand in rec.get("calls", []):
            yield pid, cand


def best_call(rec: dict) -> dict | None:
    """Return the highest-confidence candidate for one protein's record.

    `rec` is a single value from a calls.json dict (`calls[protein_id]`).
    Since `build_strain_calls` sorts `calls[]` by confidence_score descending,
    this is just the first entry. Returns None for an empty record.
    """
    cands = rec.get("calls", [])
    return cands[0] if cands else None


def call_tc_families(rec: dict) -> list[str]:
    """Sorted-unique list of 3-part TC families across one protein's candidates.

    Truncates each candidate's `tcid` to its 3-part family — the depth at
    which eggNOG operates and TCDB curates Pfam associations.
    """
    return sorted({".".join(c["tcid"].split(".")[:3]) for c in rec.get("calls", [])})


def call_tcids(rec: dict) -> list[str]:
    """Candidates' full `tcid` values in original order (confidence descending).

    Preserves 3-part / 4-part / 5-part depth. Use when the merge wants the
    most-specific TCID per call (vs `call_tc_families`, which collapses to
    3-part). Order is meaningful — `[0]` is the highest-confidence call.
    """
    return [c["tcid"] for c in rec.get("calls", [])]
