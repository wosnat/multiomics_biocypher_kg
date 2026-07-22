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
from typing import Any


def _entry_fields(entry: dict | None) -> tuple[str | None, str | None, str | None, list[str], list[str]]:
    """Pull (interpro_accession, description, type, go_terms, pathways) from a
    match's ``signature.entry`` block. ``entry`` is None for member-DB hits
    that InterPro has not integrated into an IPR entry."""
    if not entry:
        return None, None, None, [], []
    go_terms = sorted({x["id"] for x in (entry.get("goXRefs") or []) if x.get("id")})
    pathways = sorted({
        f'{x.get("databaseName", "?")}:{x["id"]}'
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


def _extract_matches(matches_json: list[dict]) -> list[dict]:
    """Flatten one protein's ``matches`` array into one record per
    (match × location), so multi-region signatures keep their coordinates.
    Sorted by (start, evalue, signature_accession)."""
    out: list[dict] = []
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
                "interpro_type": ipr_type,
                "start": loc.get("start"),
                "end": loc.get("end"),
                "evalue": loc.get("evalue"),
                "score": loc.get("score"),
                "go_terms": go_terms,
                "pathways": pathways,
            })
    out.sort(key=lambda d: (
        d["start"] if d["start"] is not None else 0,
        d["evalue"] if d["evalue"] is not None else 0.0,
        d["signature_accession"] or "",
    ))
    return out


def _aggregate(matches: list[dict]) -> dict[str, Any]:
    """Roll a protein's flattened match list into the per-protein call dict."""
    interpro_entries = sorted({m["interpro_accession"] for m in matches if m["interpro_accession"]})
    go_terms = sorted({g for m in matches for g in m["go_terms"]})
    pathways = sorted({p for m in matches for p in m["pathways"]})
    libraries = sorted({m["library"] for m in matches if m["library"]})
    return {
        "match_count": len(matches),
        "interpro_entries": interpro_entries,
        "go_terms": go_terms,
        "pathways": pathways,
        "libraries": libraries,
        "matches": matches,
    }


def parse_interproscan_json(data: dict) -> dict[str, dict]:
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
    """
    calls: dict[str, dict] = {}
    for result in data.get("results") or []:
        matches = _extract_matches(result.get("matches") or [])
        md5 = result.get("md5")
        entry = _aggregate(matches)
        for xref in result.get("xref") or []:
            wp = xref.get("id")
            if not wp:
                continue
            calls[wp] = {"md5": md5, **entry}
    return calls


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
    }
    if image_digest is not None:
        summary["image_digest"] = image_digest
    if wallclock_s is not None:
        summary["wallclock_s"] = round(wallclock_s, 1)
    return summary
