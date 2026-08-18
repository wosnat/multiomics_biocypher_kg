# multiomics_kg/utils/merops_diamond.py
"""MEROPS-vs-Diamond per-hit tier policy + per-protein consensus collapse.

Pure Python — no filesystem or subprocess (except the calls.json consumption
helpers at the bottom, which read a path they are given). The orchestrator in
`.claude/skills/merops-diamond/run_merops_diamond.py` is responsible for I/O.

**This module produces PRIMARY SEQUENCE EVIDENCE ONLY** (tcdb-diamond
principle, see `multiomics_kg/utils/tcdb_diamond.py`'s module docstring): it
derives calls from `protein.faa` vs. the curated MEROPS scan library
(`merops_scan.lib`, one representative peptidase-unit sequence per MEROPS
identifier) and the static family→clan reference (`family.txt`). It does NOT
read eggNOG annotations or `gene_annotations_merged.json` — that would create
a cycle with prepare_data step 2, where these calls are destined to be merged.

MEROPS hierarchy: clan → family (e.g. S08) → subfamily (S08A; only some
families have them) → identifier (S08.001). Identifier tails encode entry
kind: numeric 001-899 = characterized holotype, letter+num (A33) = putative,
9xx = NON-PEPTIDASE HOMOLOG (catalytically dead). Family first letter is the
catalytic type (S/C/M/A/T/G/N/U/P) or I = inhibitor family.

Design: `docs/superpowers/specs/2026-08-17-merops-diamond-design.md`.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

# ============================================================================
# MEROPS identifier / header parsing
# ============================================================================

# >MER0000002 - chymotrypsin A (cattle-type) (Bos taurus) [S01.001]#S01A#{peptidase unit: 34-263}~source XP_608091~
_SCAN_LIB_HEADER_RE = re.compile(
    r"^>(?P<mernum>\S+) - (?P<name>.*) \[(?P<merops_id>[^\]]+)\]#(?P<subfam>[^#]*)#"
)

_MEROPS_ID_RE = re.compile(r"^(?P<family>[A-Z]\d+)(?P<subletter>[A-Z]?)\.(?P<tail>[A-Za-z0-9]+)$")


def parse_scan_lib_header(line: str) -> dict | None:
    """Parse one merops_scan.lib FASTA header.

    Returns {mernum, name, merops_id, subfamily} or None on mismatch.
    `subfamily` is the raw `#...#` token: the subfamily code (S01A) when the
    family has subfamilies, else the family code itself (S33). NOTE: the file
    is latin-1 encoded — open with encoding="latin-1".
    """
    m = _SCAN_LIB_HEADER_RE.match(line)
    if not m:
        return None
    return {
        "mernum": m.group("mernum"),
        "name": m.group("name"),
        "merops_id": m.group("merops_id"),
        "subfamily": m.group("subfam"),
    }


def id_family(merops_id: str) -> str | None:
    """Family code of a MEROPS identifier: "S01.001" → "S01". None if unparseable."""
    m = _MEROPS_ID_RE.match(merops_id or "")
    return m.group("family") if m else None


def id_kind(merops_id: str) -> str | None:
    """Entry kind from the identifier tail.

    "holotype" (numeric 001-899), "nonpeptidase_homolog" (9xx), or
    "putative" (letter-prefixed tail, e.g. A33 / B01 / UPW). None if the
    identifier does not parse.
    """
    m = _MEROPS_ID_RE.match(merops_id or "")
    if not m:
        return None
    tail = m.group("tail")
    if tail.isdigit():
        return "nonpeptidase_homolog" if tail.startswith("9") else "holotype"
    return "putative"


def parse_merops_subject_id(subject_id: str) -> dict | None:
    """Extract fields from a rewritten diamond subject ID.

    The runner rewrites scan-lib headers to ``>MERNUM|merops_id|subfam_token``
    (e.g. ``MER0000002|S01.001|S01A``). Returns {mernum, merops_id, family,
    subfamily} where `subfamily` is None when the token equals the family
    (family has no subfamilies). Returns None on malformed input.
    """
    if not subject_id:
        return None
    parts = subject_id.split("|")
    if len(parts) != 3:
        return None
    mernum, merops_id, token = parts
    family = id_family(merops_id)
    if not mernum or family is None:
        return None
    return {
        "mernum": mernum,
        "merops_id": merops_id,
        "family": family,
        "subfamily": token if token and token != family else None,
    }


# ============================================================================
# Phase-2 node vocabulary helpers (KG integration)
# ============================================================================

# MEROPS catalytic-type letters spelled out for the KG (`MeropsFamily.catalytic_type`).
# Single letters like "S" are insider jargon; full words read correctly without docs.
CATALYTIC_TYPE_WORDS = {
    "S": "serine",
    "C": "cysteine",
    "M": "metallo",
    "A": "aspartic",
    "T": "threonine",
    "G": "glutamic",
    "N": "asparagine_lyase",  # asparagine peptide lyases (self-cleaving)
    "P": "mixed",             # clans/families of mixed catalytic type
    "U": "unknown",
}

_CLAN_CODE_RE = re.compile(r"^[A-Z][A-Z]$")
_FAMILY_CODE_RE = re.compile(r"^[A-Z]\d+$")
_SUBFAMILY_CODE_RE = re.compile(r"^[A-Z]\d+[A-Z]$")


def catalytic_type_word(code: str) -> str | None:
    """Full-word catalytic type for a clan/family/subfamily code.

    None for inhibitor codes (I-prefixed) — inhibitors have no catalytic type.
    """
    if not code or code.startswith("I"):
        return None
    return CATALYTIC_TYPE_WORDS.get(code[0], "unknown")


def family_type(code: str) -> str:
    """'inhibitor' for I-prefixed codes, else 'peptidase' (node `family_type`)."""
    return "inhibitor" if code.startswith("I") else "peptidase"


def classify_code(code: str) -> tuple[int, str] | None:
    """(level, level_kind) for a clan / family / subfamily code, None if malformed.

    clan `SC` → (0, 'merops_clan'); family `S14` → (1, 'merops_family');
    subfamily `S08A` → (2, 'merops_subfamily'). MEROPS identifiers (S08.001)
    are NOT node codes — callers resolve them to subfamily/family first.
    """
    if _CLAN_CODE_RE.match(code or ""):
        return 0, "merops_clan"
    if _FAMILY_CODE_RE.match(code or ""):
        return 1, "merops_family"
    if _SUBFAMILY_CODE_RE.match(code or ""):
        return 2, "merops_subfamily"
    return None


def edge_target_code(candidate: dict) -> str | None:
    """Node code a Gene_has_merops_family edge attaches to for one candidate.

    Tier-2/3 candidates attach at their called `code` (family or subfamily).
    Tier-1 candidates call a full MEROPS identifier (S26.014) — no identifier
    nodes exist (6 such calls in the whole 42-strain corpus), so they attach
    at the claimed subfamily, else the family; the identifier is preserved on
    the edge as `best_hit_id`. Returns None on malformed input.
    """
    if candidate.get("level_kind") == "merops_id":
        target = candidate.get("subfamily") or candidate.get("family")
    else:
        target = candidate.get("code")
    return target if target and classify_code(target) else None


def call_class(candidate: dict) -> str:
    """Read-first verdict for a candidate → edge `call_class`.

    'inhibitor' when the family is an I-family; 'nonpeptidase_homolog' when
    the best (nearest) hit is a catalytically dead .9xx entry; else
    'peptidase'. Threshold-free — surfaces MEROPS's own curation as a
    top-level token (expression_status pattern).
    """
    if candidate.get("entry_type") == "inhibitor":
        return "inhibitor"
    if candidate.get("best_hit_kind") == "nonpeptidase_homolog":
        return "nonpeptidase_homolog"
    return "peptidase"


_TRAILING_ORGANISM_RE = re.compile(r"\s*\([^()]*\)\s*$")


def _clean_type_example_name(name: str) -> str:
    """Strip the trailing organism parenthetical from a scan-lib entry name.

    "aminopeptidase N (Homo sapiens)" → "aminopeptidase N";
    "chymotrypsin A (cattle-type) (Bos taurus)" → "chymotrypsin A (cattle-type)"
    (only the LAST parenthetical is the organism).
    """
    return _TRAILING_ORGANISM_RE.sub("", name).strip()


def type_example_names(scan_lib_text: str) -> dict[str, str]:
    """{family_or_subfamily_code: type-example name} from merops_scan.lib text.

    MEROPS names a family after its type example; family.txt carries that name
    for only ~27% of families, so this derives the rest from the scan library:
    per code, prefer the `.001` entry, else the lowest-numbered characterized
    holotype (< .900), else the lowest putative (letter-tail) entry. Dead
    `.9xx` homologs never name a family — a family observed only through dead
    relatives keeps its bare code. Caller passes latin-1-decoded text.
    """
    # best[code] = (rank, sort_key, name); rank 0 = holotype, 1 = putative
    best: dict[str, tuple[int, str, str]] = {}

    def offer(code: str | None, tail: str, name: str) -> None:
        if not code or not name:
            return
        if tail.isdigit():
            if tail.startswith("9"):
                return  # nonpeptidase homolog — never a namesake
            cand = (0, tail.zfill(3), name)
        else:
            cand = (1, tail, name)
        if code not in best or cand < best[code]:
            best[code] = cand

    for line in scan_lib_text.splitlines():
        if not line.startswith(">"):
            continue
        rec = parse_scan_lib_header(line)
        if not rec:
            continue
        m = _MEROPS_ID_RE.match(rec["merops_id"])
        if not m:
            continue
        family = m.group("family")
        tail = m.group("tail")
        name = _clean_type_example_name(rec["name"])
        offer(family, tail, name)
        subfam = rec["subfamily"]
        if subfam and subfam != family:
            offer(subfam, tail, name)

    return {code: name for code, (_, _, name) in best.items()}


# ============================================================================
# Reference tables (family.txt / clan.txt)
# ============================================================================

def parse_clan_txt(text: str) -> dict[str, dict]:
    """Parse MEROPS `database_files/clan.txt` into {clan: {description, family_type}}.

    Tab-separated columns (observed release 12.5): clan, subclan?, description,
    NA?, reference, type ("peptidase"|"inhibitor"), active. Codes ending in
    "-" (e.g. "A-") are the "families not assigned to a clan" sentinel rows —
    skipped (families under them carry clan=None, matching parse_family_txt).
    """
    clans: dict[str, dict] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 6 or not parts[0].strip():
            continue
        clan = parts[0].strip()
        if clan.endswith("-"):
            continue
        clans[clan] = {
            "description": parts[2].strip() or None,
            "family_type": parts[5].strip() or None,
        }
    return clans


def parse_family_txt(text: str) -> dict[str, dict]:
    """Parse MEROPS `database_files/family.txt` into {family: {clan, name, entry_type}}.

    Tab-separated columns (observed release 12.5): family, clan, subclan,
    type_example_name, ..., entry_type ("peptidase"|"inhibitor"), ...
    Clan codes ending in "-" (e.g. "A-") are MEROPS's "not assigned to a clan"
    sentinel and are normalized to None.
    """
    families: dict[str, dict] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 7 or not parts[0].strip():
            continue
        family = parts[0].strip()
        clan = parts[1].strip() or None
        if clan and clan.endswith("-"):
            clan = None
        families[family] = {
            "clan": clan,
            "name": parts[3].strip() or None,
            "entry_type": parts[6].strip() or None,
        }
    return families


# ============================================================================
# Diamond row parsing + tier policy (tcdb-diamond parity)
# ============================================================================

def parse_diamond_row(line: str) -> dict | None:
    """Parse one diamond blastp output line (--outfmt 6, 8 columns) to a dict.

    Returns None when the row is malformed (wrong column count, non-numeric
    field). Same column layout as tcdb-diamond: qseqid sseqid pident qcovhsp
    scovhsp length evalue bitscore.
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


def classify_hit(hit: dict) -> int | None:
    """Assign a confidence tier (1/2/3) to a parsed diamond hit row.

    Identical thresholds to tcdb-diamond's `classify_hit` — deliberate parity:

    Floor (gblast3-style): e-value <= 0.001, HSP length >= 50, AND
    (qcov >= 40 OR scov >= 40). Above the floor:
      - tier 1: identity >= 70 AND qcov >= 70  → identifier-level claim
      - tier 2: identity >= 40 AND qcov >= 60  → subfamily-level claim
      - tier 3: floor only                     → family-level claim

    Returns None when the hit fails the floor (drop it).
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


# ============================================================================
# Consensus collapse
# ============================================================================

_AGREEMENT_WEIGHT = {"id": 1.0, "subfamily": 0.85, "family": 0.7}
_AGREEMENT_TO_TIER = {"id": 1, "subfamily": 2, "family": 3}


def family_group_agreement(hits: list[dict]) -> str:
    """Agreement depth within one family's hit group.

    "id" when all hits share one MEROPS identifier; "subfamily" when they
    share one real subfamily (token != family); else "family". Hits are
    parsed-subject dicts carrying merops_id / family / subfamily.
    """
    ids = {h["merops_id"] for h in hits}
    if len(ids) == 1:
        return "id"
    subfams = {h["subfamily"] for h in hits}
    if len(subfams) == 1 and next(iter(subfams)) is not None:
        return "subfamily"
    return "family"


def confidence_score(identity: float, qcov: float, agreement: str) -> float:
    """Continuous 0-1 confidence summary, tcdb-diamond formula:

    ``score = (identity / 100) * (qcov / 100) * agreement_weight``
    with agreement_weight 1.0 / 0.85 / 0.7 for id / subfamily / family
    consensus.
    """
    return (identity / 100.0) * (qcov / 100.0) * _AGREEMENT_WEIGHT[agreement]


def _truncate_call(family: str, subfamily: str | None, merops_id: str,
                   effective_tier: int) -> tuple[str, str, str | None]:
    """(code, level_kind, claimed_subfamily) for a call at `effective_tier`.

    MEROPS's hierarchy is ragged — families without subfamilies make tier-2
    claims land at family level (`tier` stays 2 as the confidence band;
    `level_kind` describes the actual claim depth). This is the one deliberate
    deviation from TCDB's locked tier↔depth coupling, where every TCID has
    all 5 parts.
    """
    if effective_tier == 1:
        return merops_id, "merops_id", subfamily
    if effective_tier == 2 and subfamily is not None:
        return subfamily, "merops_subfamily", subfamily
    return family, "merops_family", None


def build_strain_calls(tsv_path: Path, family_meta: dict[str, dict]) -> tuple[dict, dict]:
    """Full per-strain pipeline: parse TSV, classify, per-family consensus collapse.

    Inputs: the raw diamond TSV + the static family reference from
    `parse_family_txt` (clan / name / entry_type per family). Nothing else.

    Returns (calls, summary):
      calls:   {protein_id: {"calls": [candidate, ...]}}, candidates sorted by
               confidence_score descending
      summary: per-strain counts for the stdout status table + QC

    Each candidate carries only sequence-derived evidence plus static family
    metadata: code, family, subfamily, clan, catalytic_type, entry_type,
    level_kind, tier, confidence_score, identity, qcov, scov, evalue, length,
    consensus_n, consensus_agreement, best_hit_id, best_hit_mernum,
    best_hit_kind, homolog_hit_fraction.
    """
    by_query: dict[str, list[dict]] = defaultdict(list)
    raw_lines = 0
    parse_failures = 0
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            raw_lines += 1
            row = parse_diamond_row(line)
            if row is None:
                parse_failures += 1
                continue
            tier = classify_hit(row)
            if tier is None:
                continue
            subject = parse_merops_subject_id(row["subject_id"])
            if subject is None:
                parse_failures += 1
                continue
            by_query[row["query_id"]].append({
                **subject,
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
        # Group hits by MEROPS family. Each family produces one candidate;
        # multi-domain proteins (e.g. a peptidase fused to an inhibitor
        # domain) emit multiple candidates instead of being rejected.
        by_family: dict[str, list[dict]] = defaultdict(list)
        for h in hits:
            by_family[h["family"]].append(h)

        candidates: list[dict] = []
        for family, family_hits in by_family.items():
            agreement = family_group_agreement(family_hits)
            best_tier = min(h["tier"] for h in family_hits)
            depth_tier = _AGREEMENT_TO_TIER[agreement]
            effective_tier = max(best_tier, depth_tier)

            # Best (highest-identity) hit drives the candidate's metadata
            best = max(family_hits, key=lambda h: h["identity"])
            code, level_kind, claimed_subfam = _truncate_call(
                family, best["subfamily"], best["merops_id"], effective_tier
            )

            n_homolog = sum(
                1 for h in family_hits
                if id_kind(h["merops_id"]) == "nonpeptidase_homolog"
            )
            meta = family_meta.get(family, {})
            is_inhibitor = family.startswith("I")

            candidates.append({
                "code": code,
                "family": family,
                "subfamily": claimed_subfam,
                "clan": meta.get("clan"),
                "catalytic_type": None if is_inhibitor else family[0],
                "entry_type": "inhibitor" if is_inhibitor else "peptidase",
                "level_kind": level_kind,
                "tier": effective_tier,
                "confidence_score": round(
                    confidence_score(best["identity"], best["qcov"], agreement), 4
                ),
                "identity": best["identity"],
                "qcov": best["qcov"],
                "scov": best["scov"],
                "evalue": best["evalue"],
                "length": best["length"],
                "consensus_n": len(family_hits),
                "consensus_agreement": agreement,
                "best_hit_id": best["merops_id"],
                "best_hit_mernum": best["mernum"],
                "best_hit_kind": id_kind(best["merops_id"]),
                "homolog_hit_fraction": round(n_homolog / len(family_hits), 3),
            })

        if not candidates:
            continue
        candidates.sort(key=lambda c: -c["confidence_score"])
        calls[query_id] = {"calls": candidates}
        total_candidates += len(candidates)

    # Distributions count CANDIDATES (not proteins), tcdb-diamond convention.
    tier_dist: dict[str, int] = defaultdict(int)
    agreement_dist: dict[str, int] = defaultdict(
        int, {"id": 0, "subfamily": 0, "family": 0}
    )
    catalytic_dist: dict[str, int] = defaultdict(int)
    entry_type_dist: dict[str, int] = defaultdict(int)
    best_hit_kind_dist: dict[str, int] = defaultdict(int)
    candidates_per_protein: dict[str, int] = defaultdict(int)
    for rec in calls.values():
        candidates_per_protein[str(len(rec["calls"]))] += 1
        for cand in rec["calls"]:
            tier_dist[str(cand["tier"])] += 1
            agreement_dist[cand["consensus_agreement"]] += 1
            catalytic_dist[cand["catalytic_type"] or "inhibitor"] += 1
            entry_type_dist[cand["entry_type"]] += 1
            best_hit_kind_dist[cand["best_hit_kind"] or "unknown"] += 1

    summary = {
        "raw_hit_lines": raw_lines,
        "parse_failures": parse_failures,
        "proteins_with_hits": len(by_query),
        "proteins_with_call": len(calls),
        "total_candidates": total_candidates,
        "candidates_per_protein_distribution": dict(candidates_per_protein),
        "tier_distribution": dict(tier_dist),
        "consensus_agreement_distribution": dict(agreement_dist),
        "catalytic_type_distribution": dict(catalytic_dist),
        "entry_type_distribution": dict(entry_type_dist),
        "best_hit_kind_distribution": dict(best_hit_kind_dist),
    }
    return calls, summary


# ============================================================================
# calls.json consumption helpers (tcdb-diamond parity)
# ============================================================================

def load_calls_json(path: Path) -> dict[str, dict]:
    """Load a single strain's `<strain>.merops.calls.json` into a dict.

    Raises FileNotFoundError if the file is missing — a strain without a
    calls.json is a real "no MEROPS annotation available" case, not silent.
    """
    return json.loads(Path(path).read_text())


def iter_candidates(calls: dict[str, dict]) -> Iterator[tuple[str, dict]]:
    """Yield (protein_id, candidate_dict) for every candidate.

    Nothing is filtered — the tier policy already applied the quality gate at
    build time. Consumers wanting a stricter cut (e.g. excluding
    `best_hit_kind == "nonpeptidase_homolog"` or `entry_type == "inhibitor"`)
    filter explicitly, so the threshold is visible at the point of use.
    """
    for pid, rec in calls.items():
        for cand in rec.get("calls", []):
            yield pid, cand


def best_call(rec: dict) -> dict | None:
    """Highest-confidence candidate for one protein's record (calls[0])."""
    cands = rec.get("calls", [])
    return cands[0] if cands else None


def call_families(rec: dict) -> list[str]:
    """Sorted-unique MEROPS family codes across one protein's candidates."""
    return sorted({c["family"] for c in rec.get("calls", [])})


def call_codes(rec: dict) -> list[str]:
    """Candidates' called `code` values in confidence-descending order."""
    return [c["code"] for c in rec.get("calls", [])]
