"""Pure-Python parsing + summarizing of InterProScan JSON output.

No filesystem, no subprocess — the orchestrator
(``.claude/skills/interproscan-run/run_interproscan.py``) handles Docker and
I/O; this module is the unit-testable core. It turns one strain's raw
InterProScan JSON (``--formats JSON``) into the Phase-1 WP_-keyed calls dict
and the per-strain QC summary.

InterProScan is a *per-protein list-of-matches* tool (Step 0 Q1): one protein
yields many matches across member databases, each integrated (or not) into an
InterPro entry. See
``docs/superpowers/specs/2026-07-22-interproscan-domains-design.md``.
"""

from __future__ import annotations

import collections
import re
from typing import Any


_REACTOME_SPECIES_RE = re.compile(r"^R-[A-Z]{3}-(\d+)$")


def normalize_interpro_type(raw: str | None) -> str:
    """House rule R1 — lowercase snake_case. InterPro ships UPPERCASE
    (``FAMILY``, ``HOMOLOGOUS_SUPERFAMILY``, …)."""
    return (raw or "").strip().lower()


def normalize_library(raw: str | None) -> str:
    """House rule R1 — lowercase snake_case member-DB name. InterProScan
    ships UPPERCASE (``PFAM``, ``PROSITE_PATTERNS``, …)."""
    return (raw or "").strip().lower()


def normalize_pathway_xref(database: str, raw_id: str) -> str:
    """Render one ``pathwayXRefs`` entry as a ``DB:id`` string.

    Reactome stable ids are *species-scoped* (``R-HSA-73817`` human,
    ``R-MMU-73817`` mouse, ``R-DME-73817`` fly, …) and InterPro lists every
    species projection of the same curated event. On a marine-bacterial
    proteome that is pure duplication — 16 species prefixes were observed in
    the MED4 smoke, inflating the artifact ~10x — and the species is
    meaningless for our organisms. We keep only the trailing numeric stable
    id, which is the part shared across all projections: ``Reactome:73817``.
    Expand back to any species form as ``R-<SPECIES>-<id>`` if ever needed.

    Non-Reactome databases (MetaCyc, and KEGG on releases that still ship it)
    pass through untouched.
    """
    if database == "Reactome":
        m = _REACTOME_SPECIES_RE.match(raw_id)
        if m:
            return f"Reactome:{m.group(1)}"
    return f"{database}:{raw_id}"


def _entry_fields(entry: dict | None) -> tuple[str | None, str | None, str | None, list[str], list[str]]:
    """Pull (interpro_accession, description, type, go_terms, pathways) from a
    match's ``signature.entry`` block. ``entry`` is None for member-DB hits
    that InterPro has not integrated into an IPR entry."""
    if not entry:
        return None, None, None, [], []
    go_terms = sorted({x["id"] for x in (entry.get("goXRefs") or []) if x.get("id")})
    pathways = sorted({
        normalize_pathway_xref(x.get("databaseName", "?"), x["id"])
        for x in (entry.get("pathwayXRefs") or [])
        if x.get("id")
    })
    return (
        entry.get("accession"),
        entry.get("description"),
        entry.get("type"),
        go_terms,
        pathways,
    )


def _extract_matches(
    matches_json: list[dict], entry_xrefs: dict[str, dict[str, list[str]]] | None = None
) -> tuple[list[dict], list[str], list[str]]:
    """Flatten one protein's ``matches`` array into one record per
    (match × location), so multi-region signatures keep their coordinates.
    Sorted by (start, evalue, signature_accession).

    Returns ``(records, go_terms, pathways)`` where the two lists are the
    protein-level unions over every integrated entry. When *entry_xrefs* is
    given it is populated in place with ``{IPR: {go_terms, pathways}}`` — the
    normalized per-entry detail that replaces the old per-match copies.
    """
    out: list[dict] = []
    go_union: set[str] = set()
    pw_union: set[str] = set()
    for m in matches_json or []:
        sig = m.get("signature") or {}
        lib_release = sig.get("signatureLibraryRelease") or {}
        library = lib_release.get("library")
        ipr_acc, ipr_desc, ipr_type, go_terms, pathways = _entry_fields(sig.get("entry"))
        sig_acc = sig.get("accession")
        sig_desc = sig.get("description") or sig.get("name")
        for loc in m.get("locations") or []:
            out.append({
                "library": library,
                "signature_accession": sig_acc,
                "signature_description": sig_desc,
                "interpro_accession": ipr_acc,
                "interpro_description": ipr_desc,
                # ipr_type is None for unintegrated member-DB hits (test_interproscan.py
                # asserts this null semantic survives R1) — only lowercase real values.
                "interpro_type": normalize_interpro_type(ipr_type) if ipr_type is not None else None,
                "start": loc.get("start"),
                "end": loc.get("end"),
                "evalue": loc.get("evalue"),
                "score": loc.get("score"),
                # GO/pathway xrefs are a property of the *InterPro entry*, not of
                # this match, so they are NOT repeated here — that duplication
                # tripled the artifact once --goterms/--pathways were enabled.
                # The protein-level union lives on the call; the per-entry detail
                # lives in the sibling entry_xrefs table, joined on
                # `interpro_accession`. Nothing is lost.
            })
        go_union.update(go_terms)
        pw_union.update(pathways)
        if ipr_acc and entry_xrefs is not None and (go_terms or pathways):
            entry_xrefs.setdefault(ipr_acc, {"go_terms": go_terms, "pathways": pathways})
    out.sort(key=lambda d: (
        d["start"] if d["start"] is not None else 0,
        d["evalue"] if d["evalue"] is not None else 0.0,
        d["signature_accession"] or "",
    ))
    return out, sorted(go_union), sorted(pw_union)


def _aggregate(matches: list[dict], go_terms: list[str], pathways: list[str]) -> dict[str, Any]:
    """Roll a protein's flattened match list into the per-protein call dict."""
    interpro_entries = sorted({m["interpro_accession"] for m in matches if m["interpro_accession"]})
    # rolled-up summary is normalized to house-rule casing (R1); the raw
    # per-match "library" token above is left as InterProScan wrote it.
    libraries = sorted({normalize_library(m["library"]) for m in matches if m["library"]})
    return {
        "match_count": len(matches),
        "interpro_entries": interpro_entries,
        "go_terms": go_terms,
        "pathways": pathways,
        "libraries": libraries,
        "matches": matches,
    }


def parse_interproscan_json(
    data: dict, entry_xrefs: dict[str, dict[str, list[str]]] | None = None
) -> dict[str, dict]:
    """Parse an InterProScan JSON document into a WP_-keyed calls dict.

    ``data`` is the loaded ``--formats JSON`` output: ``{"results": [...]}``.
    Each result has ``xref`` (list of ``{id, name}`` — InterProScan dedups
    identical sequences, so one result may back several accessions), ``md5``,
    and ``matches``. We fan every result's matches out to every xref id (the
    FASTA first-token = RefSeq WP_ accession).

    Proteins processed but matching nothing appear in the input results with an
    empty ``matches`` array and are kept here with ``match_count == 0`` — so a
    missing key means "not processed" while a zero-match entry means "no domain
    found".

    Pass *entry_xrefs* (a dict) to also collect the per-InterPro-entry GO and
    pathway detail; see :func:`parse_entry_xrefs`.
    """
    calls: dict[str, dict] = {}
    for result in data.get("results") or []:
        matches, go_terms, pathways = _extract_matches(result.get("matches") or [], entry_xrefs)
        md5 = result.get("md5")
        entry = _aggregate(matches, go_terms, pathways)
        for xref in result.get("xref") or []:
            wp = xref.get("id")
            if not wp:
                continue
            calls[wp] = {"md5": md5, **entry}
    return calls


def parse_entry_xrefs(data: dict) -> dict[str, dict[str, list[str]]]:
    """Build the ``{IPR: {go_terms, pathways}}`` side table for one strain.

    Written alongside the calls as ``<strain>.interproscan.entry_xrefs.json``.
    GO/pathway xrefs belong to the InterPro *entry*, so storing them once per
    entry (rather than on every match of that entry) keeps the artifact small
    while staying lossless — join back on a match's ``interpro_accession``.
    Only entries that actually carry at least one xref appear.
    """
    entry_xrefs: dict[str, dict[str, list[str]]] = {}
    for result in data.get("results") or []:
        _extract_matches(result.get("matches") or [], entry_xrefs)
    return dict(sorted(entry_xrefs.items()))


def summarize(
    calls: dict[str, dict],
    *,
    strain: str,
    input_proteins: int,
    tool_version: str,
    applications: str,
    image_digest: str | None = None,
    wallclock_s: float | None = None,
    parse_failures: int = 0,
    xrefs_requested: bool | None = None,
) -> dict[str, Any]:
    """Build the per-strain ``skill_summary.json`` QC dict from parsed calls."""
    calls_made = len(calls)
    proteins_no_match = sum(1 for c in calls.values() if c["match_count"] == 0)
    total_matches = sum(c["match_count"] for c in calls.values())
    interpro_integrated = sum(
        1 for c in calls.values() for m in c["matches"] if m["interpro_accession"]
    )
    lib_counter: collections.Counter = collections.Counter()
    for c in calls.values():
        for m in c["matches"]:
            if m["library"]:
                lib_counter[m["library"]] += 1

    # GO / pathway xref coverage — the QC that tells you --goterms/--pathways
    # actually took effect (all four counters are 0 when they were omitted).
    go_seen: set[str] = set()
    pw_seen: set[str] = set()
    pw_db_counter: collections.Counter = collections.Counter()
    proteins_with_go = 0
    proteins_with_pathways = 0
    for c in calls.values():
        if c["go_terms"]:
            proteins_with_go += 1
            go_seen.update(c["go_terms"])
        if c["pathways"]:
            proteins_with_pathways += 1
            pw_seen.update(c["pathways"])
            for p in c["pathways"]:
                pw_db_counter[p.split(":", 1)[0]] += 1
    summary: dict[str, Any] = {
        "strain": strain,
        "tool_version": tool_version,
        "applications": applications,
        "input_proteins": input_proteins,
        "calls_made": calls_made,
        "proteins_no_match": proteins_no_match,
        "parse_failures": parse_failures,
        "total_matches": total_matches,
        "interpro_integrated_matches": interpro_integrated,
        "distribution": dict(sorted(lib_counter.items())),
        "sentinel_rate": round(proteins_no_match / input_proteins, 4) if input_proteins else 0.0,
        "proteins_with_go_terms": proteins_with_go,
        "distinct_go_terms": len(go_seen),
        "proteins_with_pathways": proteins_with_pathways,
        "distinct_pathways": len(pw_seen),
        "pathway_databases": dict(sorted(pw_db_counter.items())),
    }
    if xrefs_requested is not None:
        summary["xrefs_requested"] = xrefs_requested
    if image_digest is not None:
        summary["image_digest"] = image_digest
    if wallclock_s is not None:
        summary["wallclock_s"] = round(wallclock_s, 1)
    return summary
