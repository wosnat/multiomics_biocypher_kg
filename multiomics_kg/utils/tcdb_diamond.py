# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein post-steps.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator


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


_AGREEMENT_RANK = {"confirms": 0, "refines": 1, "conflicts": 2}


def _pair_agreement(diamond_tcid: str, egn_tcid: str) -> str:
    """Per-pair agreement between one diamond TCID and one eggNOG TCID.

    Returns "confirms" | "refines" | "conflicts". Caller handles the
    `extends` (no eggNOG values) case.
    """
    if diamond_tcid == egn_tcid:
        return "confirms"

    diamond_parts = diamond_tcid.split(".")
    egn_parts = egn_tcid.split(".")

    if diamond_parts[:3] != egn_parts[:3]:
        return "conflicts"

    if len(diamond_parts) > len(egn_parts) and diamond_parts[: len(egn_parts)] == egn_parts:
        return "refines"
    return "confirms"


def compute_egn_agreement(diamond_tcid: str, egn_tcids: list[str] | str | None) -> str:
    """Tag the relationship between the diamond call and eggNOG's KEGG_TC values.

    eggNOG's KEGG_TC field is multi-valued (comma-separated in the source TSV)
    — e.g. `1.A.33.1,9.B.157.1` for MreB-like proteins. The previous behavior
    of inspecting only the first value misclassified diamond calls matching any
    later value as `conflicts`. This function now considers ALL eggNOG values
    and returns the strongest match.

    Returns one of: "confirms" | "refines" | "extends" | "conflicts".
    `egn_only` (eggNOG TC present, diamond absent) is not produced here —
    those proteins simply don't appear in the calls JSON in Phase 1.

    Rules (applied per-pair, then aggregated by best match):
      - confirms: identical TCIDs OR one is a strict prefix of the other
        AT family level or below (first 3 parts match)
      - refines: eggNOG TCID is a strict prefix of diamond TCID — same lineage,
        diamond went deeper. Reported separately from confirms because this
        is the headline specificity win.
      - extends: eggNOG had no TC values; diamond produced one
      - conflicts: ALL eggNOG values disagree at family level (first 3 parts)

    `egn_tcids` accepts a list, a single string (legacy), or None.
    """
    if egn_tcids is None or egn_tcids == "":
        return "extends"
    if isinstance(egn_tcids, str):
        egn_tcids = [egn_tcids]
    if not egn_tcids:
        return "extends"

    # Best (lowest-ranked) match wins: confirms > refines > conflicts
    return min(
        (_pair_agreement(diamond_tcid, e) for e in egn_tcids),
        key=lambda a: _AGREEMENT_RANK[a],
    )


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


def load_eggnog_kegg_tc(annotations_path: Path) -> dict[str, list[str]]:
    """Read an eggNOG-mapper .emapper.annotations file -> {protein_id: [kegg_tc, ...]}.

    eggNOG's KEGG_TC field is multi-valued (comma-separated). For example, the
    rod shape-determining protein MreB carries `1.A.33.1,9.B.157.1` — both the
    legacy Hsp70-cation-channel call and the correct MreBCD-family call. All
    values are preserved so `compute_egn_agreement` can inspect each in turn.

    Returns an empty dict when the file is absent. Skips rows whose KEGG_TC
    column is "-", empty, or missing.
    """
    annotations_path = Path(annotations_path)
    if not annotations_path.exists():
        return {}

    result: dict[str, list[str]] = {}
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
            tcs = [v.strip() for v in kegg_tc_raw.split(",") if v.strip()]
            if tcs:
                result[protein_id] = tcs
    return result


def load_pfam_to_tc_map(path: Path) -> dict[str, set[str]]:
    """Load TCDB's curated Pfam → TC mapping into {pfam_id: {tc_family, ...}}.

    Source: https://www.tcdb.org/cgi-bin/projectv/public/pfam.py — a TCDB-
    maintained 3-column TSV (PF#####\\tTC_ID\\tfamily_name) listing which Pfam
    domains are associated with which TCDB transporter families. The cached
    file is built by `run_tcdb_diamond.py` alongside `tcdb.dmnd` and lives at
    `~/tools/TCDB/DB/tcdb_pfam_map.tsv` (or `$TCDB_DATA_DIR/DB/...`).

    TC IDs are truncated to 3-part families (`tc_family` level) for
    matching; this is the granularity at which TCDB curates Pfam ↔ TC
    associations and the granularity at which `compute_egn_agreement`
    decides confirms-vs-conflicts.

    Returns an empty dict when the file is absent.
    """
    path = Path(path)
    if not path.exists():
        return {}
    result: dict[str, set[str]] = defaultdict(set)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pfam_id, tc_id = parts[0], parts[1]
            if not pfam_id.startswith("PF") or tc_id.count(".") < 2:
                continue
            family = ".".join(tc_id.split(".")[:3])
            result[pfam_id].add(family)
    return dict(result)


def load_gene_pfams(annotations_path: Path) -> dict[str, list[str]]:
    """Read `gene_annotations_merged.json` -> {protein_id: [pfam_id, ...]}.

    The merged annotation file is keyed by locus_tag; this transposes it to
    a protein_id keyed lookup so `build_strain_calls` (whose calls.json is
    keyed by protein_id) can join Pfam annotations to diamond calls in a
    single dict access.

    Returns an empty dict when the file is absent. Skips genes without a
    `protein_id` (or `protein_id_refseq` fallback) — those can't appear in
    the diamond TSV's query column anyway.
    """
    annotations_path = Path(annotations_path)
    if not annotations_path.exists():
        return {}
    raw = json.loads(annotations_path.read_text())
    result: dict[str, list[str]] = {}
    for _lt, gene in raw.items():
        pid = gene.get("protein_id") or gene.get("protein_id_refseq")
        if not pid:
            continue
        pfams = gene.get("pfam_ids") or []
        if pfams:
            result[pid] = list(pfams)
    return result


def gene_pfam_implied_families(
    gene_pfams: list[str],
    pfam_to_tc_map: dict[str, set[str]],
) -> list[str]:
    """Union the 3-part TC families implied by all of a gene's Pfams.

    Returns a sorted list of unique 3-part TC families. Empty when the gene
    has no Pfams or none of them are in TCDB's curated Pfam→TC map.

    The TCDB Pfam→TC map (~1.3K Pfams, ~8.3K (Pfam, TC) pairs) only covers
    Pfams TCDB has explicitly curated against TC families; many transport-
    related Pfams (e.g. PF04193 for SWEET-family sugar transporters) are
    absent. Genes with only un-curated Pfams produce an empty result here
    and downstream `pfam_agreement` becomes `"neutral"`.
    """
    if not gene_pfams:
        return []
    implied: set[str] = set()
    for p in gene_pfams:
        implied |= pfam_to_tc_map.get(p, set())
    return sorted(implied)


def compute_pfam_agreement(
    diamond_tcid: str,
    egn_tcids: list[str],
    pfam_implied_families: list[str],
) -> str:
    """Tag the relationship between one diamond/eggNOG candidate pair and the
    gene's precomputed Pfam-implied TC families.

    Returns one of: "confirms_diamond" | "confirms_eggnog" | "confirms_both" |
    "contradicts_both" | "neutral". Independent from `egn_agreement`: a
    `conflicts` egn_agreement combined with `confirms_diamond` pfam_agreement
    is a strong "diamond wins" signal; `conflicts` + `confirms_eggnog` is a
    strong "eggNOG wins" signal. Phase 2's merge rule combines both.

    Tag semantics:
      - confirms_both: Pfam-implied family includes BOTH diamond's call's
        family AND at least one eggNOG family. Multi-domain protein, OR a
        Pfam that TCDB curates against multiple TC families.
      - confirms_diamond: matches diamond's family but no eggNOG family.
      - confirms_eggnog: matches an eggNOG family but not diamond's.
      - contradicts_both: Pfam implies one or more families, but none match
        diamond or any eggNOG. Both sequence-based calls may be wrong.
      - neutral: Pfam implies no families (empty list — gene has no Pfams,
        or its Pfams aren't in TCDB's curated map).

    `pfam_implied_families` is computed once per protein via
    `gene_pfam_implied_families`; this function is called per candidate so
    each call's diamond TCID gets its own verdict against the same Pfam set.
    """
    if not pfam_implied_families:
        return "neutral"
    implied = set(pfam_implied_families)
    diamond_family = ".".join(diamond_tcid.split(".")[:3])
    egn_families = {".".join(e.split(".")[:3]) for e in egn_tcids}
    d_match = diamond_family in implied
    e_match = bool(implied & egn_families)
    if d_match and e_match:
        return "confirms_both"
    if d_match:
        return "confirms_diamond"
    if e_match:
        return "confirms_eggnog"
    return "contradicts_both"


_PFAM_CONFIRM = {"confirms_diamond", "confirms_both"}
_EGN_CONFIRM = {"confirms", "refines"}
DEFAULT_MIN_RELATIVE_CONFIDENCE = 0.25
DEFAULT_MIN_SINGLETON_SCORE = 0.20


def annotate_candidate_filters(
    calls: dict[str, dict],
    min_relative_confidence: float = DEFAULT_MIN_RELATIVE_CONFIDENCE,
    min_singleton_score: float = DEFAULT_MIN_SINGLETON_SCORE,
) -> dict[str, int]:
    """Set per-candidate `filter_action` based on Pfam/eggNOG agreement,
    sibling-relative confidence, and an absolute single-hit floor.
    Mutates `calls` in place; returns drop-stat counts.

    Walks each candidate and applies the first matching rule:

      1. drop_pfam_contradicts (multi-candidate only) — this candidate has
         pfam_agreement 'contradicts_both' AND some sibling has
         pfam_agreement in {confirms_diamond, confirms_both}. Pfam
         disqualifies this call.
      2. drop_egn_conflicts (multi-candidate only) — this candidate has
         egn_agreement 'conflicts' AND some sibling has egn_agreement in
         {confirms, refines}. eggNOG disqualifies this call.
      3. drop_singleton_low_score — this candidate is backed by a single
         hit (consensus_n == 1) AND confidence_score < min_singleton_score.
         Applies REGARDLESS of sibling count, including single-candidate
         proteins: a lone weak single hit is the weakest evidence we have
         and shouldn't introduce a new TC family annotation downstream.
      4. drop_low_confidence (multi-candidate only) — this candidate's
         confidence_score is below min_relative_confidence × the
         protein's best candidate's score. Diamond evidence much weaker
         than the alternative.
      5. keep — none of the above.

    The pfam/egn/relative-confidence rules require siblings (multi-candidate)
    by design — they're alternative-aware. The singleton rule is absolute and
    fires on single-candidate proteins too: weak single-hit evidence is
    classified as filtered so Phase 2 doesn't add new TC family annotations
    on this basis. Nothing is deleted from `calls[]`; the rule is purely
    annotative.

    Returns a Counter-style dict keyed by filter_action.
    """
    drop_stats: dict[str, int] = defaultdict(int)
    for rec in calls.values():
        cands = rec["calls"]
        any_pfam_confirm = any(c["pfam_agreement"] in _PFAM_CONFIRM for c in cands)
        any_egn_confirm = any(c["egn_agreement"] in _EGN_CONFIRM for c in cands)
        max_conf = max(c["confidence_score"] for c in cands)
        multi = len(cands) > 1
        for c in cands:
            if multi and c["pfam_agreement"] == "contradicts_both" and any_pfam_confirm:
                action = "drop_pfam_contradicts"
            elif multi and c["egn_agreement"] == "conflicts" and any_egn_confirm:
                action = "drop_egn_conflicts"
            elif c["consensus_n"] == 1 and c["confidence_score"] < min_singleton_score:
                action = "drop_singleton_low_score"
            elif multi and c["confidence_score"] < min_relative_confidence * max_conf:
                action = "drop_low_confidence"
            else:
                action = "keep"
            c["filter_action"] = action
            drop_stats[action] += 1
    return dict(drop_stats)


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


def build_strain_calls(
    tsv_path: Path,
    eggnog_annotations_path: Path,
    gene_annotations_path: Path | None = None,
    pfam_to_tc_map_path: Path | None = None,
) -> tuple[dict, dict]:
    """Run the full per-strain pipeline: parse TSV, classify, consensus, tag.

    Returns (calls, summary):
      calls: dict keyed by protein_id with the full §6.5 record shape
      summary: dict with raw counts (matches §6.6 stdout columns)

    `gene_annotations_path` (per-strain `gene_annotations_merged.json`) and
    `pfam_to_tc_map_path` (TCDB-published Pfam→TC map at
    `~/tools/TCDB/DB/tcdb_pfam_map.tsv`) are optional. When both are provided,
    each call gains `pfam_ids`, `pfam_tc_families`, and `pfam_agreement` fields
    (and the summary gains `pfam_agreement_distribution`); when either is
    missing, those fields are still emitted with neutral values so the schema
    stays stable.
    """
    egn_lookup = load_eggnog_kegg_tc(Path(eggnog_annotations_path))
    gene_pfams_lookup = (
        load_gene_pfams(Path(gene_annotations_path))
        if gene_annotations_path is not None
        else {}
    )
    pfam_to_tc_map = (
        load_pfam_to_tc_map(Path(pfam_to_tc_map_path))
        if pfam_to_tc_map_path is not None
        else {}
    )

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

        egn_tcids = egn_lookup.get(query_id, [])
        gene_pfams = gene_pfams_lookup.get(query_id, [])
        pfam_tc_families = gene_pfam_implied_families(gene_pfams, pfam_to_tc_map)

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
                "egn_agreement": compute_egn_agreement(called_tcid, egn_tcids),
                "pfam_agreement": compute_pfam_agreement(called_tcid, egn_tcids, pfam_tc_families),
                "incompletely_characterized": is_class_9(called_tcid),
            })

        if not candidates:
            continue

        candidates.sort(key=lambda c: -c["confidence_score"])
        calls[query_id] = {
            "egn_tcids": egn_tcids,
            "pfam_ids": gene_pfams,
            "pfam_tc_families": pfam_tc_families,
            "calls": candidates,
        }
        total_candidates += len(candidates)

    # Apply the per-candidate filter to annotate each candidate with
    # `filter_action`. Does NOT remove candidates from `calls`; consumers
    # use the field to decide what to keep. See annotate_candidate_filters.
    filter_action_dist = annotate_candidate_filters(calls)

    # Build summary. tier_distribution + agreement distributions now count
    # CANDIDATES (not proteins): a 2-family protein contributes 2 entries.
    # Use proteins_with_call vs total_candidates to gauge how many proteins
    # are multi-family.
    tier_dist: dict[str, int] = defaultdict(int)
    agreement_dist: dict[str, int] = defaultdict(
        int, {"confirms": 0, "refines": 0, "extends": 0, "conflicts": 0}
    )
    pfam_agreement_dist: dict[str, int] = defaultdict(
        int,
        {"confirms_diamond": 0, "confirms_eggnog": 0, "confirms_both": 0,
         "contradicts_both": 0, "neutral": 0},
    )
    candidates_per_protein: dict[str, int] = defaultdict(int)
    for rec in calls.values():
        candidates_per_protein[str(len(rec["calls"]))] += 1
        for cand in rec["calls"]:
            tier_dist[str(cand["tier"])] += 1
            agreement_dist[cand["egn_agreement"]] += 1
            pfam_agreement_dist[cand["pfam_agreement"]] += 1

    summary = {
        "raw_hit_lines": raw_lines,
        "proteins_with_hits": len(by_query),
        "proteins_with_call": len(calls),
        "total_candidates": total_candidates,
        "candidates_per_protein_distribution": dict(candidates_per_protein),
        "tier_distribution": dict(tier_dist),
        "agreement_distribution": dict(agreement_dist),
        "pfam_agreement_distribution": dict(pfam_agreement_dist),
        "filter_action_distribution": dict(filter_action_dist),
    }
    return calls, summary


# ============================================================================
# Phase 2 — calls.json consumption helpers
# ============================================================================
#
# These utilities walk the per-strain `<strain>.tcdb.calls.json` artifact
# produced by `build_strain_calls`. Phase 2's merge code uses them to:
#   - load a strain's calls
#   - iterate only the candidates Phase 1 endorses (filter_action == "keep")
#   - pick the top-confidence kept candidate per protein
#   - extract per-protein TC family lists for annotation merge


# Public vocabulary — exported for Phase 2 consumers
KEEP_ACTION = "keep"
FILTER_ACTIONS = (
    "keep",
    "drop_pfam_contradicts",
    "drop_egn_conflicts",
    "drop_singleton_low_score",
    "drop_low_confidence",
)


def load_calls_json(path: Path) -> dict[str, dict]:
    """Load a single strain's `<strain>.tcdb.calls.json` into a dict.

    Raises FileNotFoundError if the file is missing — callers should
    handle this explicitly (a strain without a calls.json is a real
    "no TCDB annotation available for this strain" case, not silent).
    Returns the raw on-disk shape: `{protein_id: {egn_tcids, pfam_ids,
    pfam_tc_families, calls: [...]}}`.
    """
    path = Path(path)
    return json.loads(path.read_text())


def iter_kept_candidates(calls: dict[str, dict]) -> Iterator[tuple[str, dict]]:
    """Yield (protein_id, candidate_dict) for every candidate the filter
    endorsed (filter_action == "keep").

    Skips dropped candidates entirely. A protein contributes 0 to N
    candidates depending on how many of its `calls[]` were kept. Order
    within a protein matches `calls[]` ordering (descending by
    `confidence_score`).
    """
    for pid, rec in calls.items():
        for cand in rec.get("calls", []):
            if cand.get("filter_action") == KEEP_ACTION:
                yield pid, cand


def best_kept_call(rec: dict) -> dict | None:
    """Return the highest-confidence kept candidate for one protein's record,
    or None if every candidate was filtered out.

    `rec` is a single value from a calls.json dict (`calls[protein_id]`).
    Since `build_strain_calls` already sorts `calls[]` by confidence_score
    descending, this is just the first kept entry — no re-sort needed.
    """
    for cand in rec.get("calls", []):
        if cand.get("filter_action") == KEEP_ACTION:
            return cand
    return None


def kept_tc_families(rec: dict) -> list[str]:
    """Sorted-unique list of 3-part TC families across kept candidates for
    one protein. Empty when every candidate was filtered out.

    Truncates each candidate's `tcid` to its 3-part family. Useful when
    Phase 2 wants to merge a Gene-to-TC_family edge per kept family
    (the depth at which eggNOG operates and TCDB curates Pfam mappings).
    """
    fams: set[str] = set()
    for cand in rec.get("calls", []):
        if cand.get("filter_action") != KEEP_ACTION:
            continue
        fams.add(".".join(cand["tcid"].split(".")[:3]))
    return sorted(fams)


def kept_call_tcids(rec: dict) -> list[str]:
    """List of kept candidates' full `tcid` values in original order (by
    confidence_score descending). Preserves 3-part / 4-part / 5-part depth.

    Use when Phase 2 wants the most-specific TCID per kept call (vs
    `kept_tc_families` which collapses to 3-part). Order is meaningful —
    `[0]` is the highest-confidence kept call.
    """
    return [
        c["tcid"] for c in rec.get("calls", [])
        if c.get("filter_action") == KEEP_ACTION
    ]


def summarize_filter_actions(calls: dict[str, dict]) -> dict[str, int]:
    """Recompute the filter_action distribution from a calls dict.

    Useful when Phase 2 (or an ad-hoc analysis) has re-filtered the calls
    in memory and needs an updated count. For the original on-disk
    distribution, prefer the `filter_action_distribution` field in
    `<strain>.tcdb.skill_summary.json` — this function gives equivalent
    output but recounts from the calls dict.
    """
    counts: dict[str, int] = {a: 0 for a in FILTER_ACTIONS}
    for rec in calls.values():
        for cand in rec.get("calls", []):
            action = cand.get("filter_action", "")
            if action in counts:
                counts[action] += 1
    return counts
