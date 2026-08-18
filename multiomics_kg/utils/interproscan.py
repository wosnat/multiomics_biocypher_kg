"""Pure-Python parsing + summarizing of InterProScan JSON output.

No filesystem, no subprocess — the orchestrator
(``.claude/skills/interproscan-run/run_interproscan.py``) handles Docker and
I/O; this module is the unit-testable core. It turns one strain's raw
InterProScan JSON (``--formats JSON``) into the multi-ontology, WP_-keyed
faceted calls dict and the per-strain QC summary.

InterProScan is a *per-protein list-of-matches* tool (Step 0 Q1): one protein
yields many matches across member databases, each optionally integrated into
an InterPro entry. This module produces three sparse facets per protein:

- ``libraries``: one list of match rows per member DB that matched (each row
  carries its own ``ipr`` attribution — no shared "unattributed" bucket).
- ``interpro_entries``: a rollup keyed by IPR accession, counting matches and
  tracking the min evalue + the library that reported it ("count, don't
  combine" — no synthesized cross-library score).
- ``go_terms``: GO id -> sorted list of the InterPro entries that carry it,
  so GO provenance is attributable back to the entry.

No pathway xrefs anywhere — pathways are out of scope for this redesign (see
``docs/superpowers/specs/2026-07-22-interproscan-domains-design.md`` and the
multi-ontology redesign SDD).

This module stores everything in the tool's *native* casing (e.g. member-DB
names like ``PFAM``, InterPro types like ``FAMILY``). These are InterPro's /
InterProScan's own controlled vocabulary terms, preserved verbatim end to
end (parser through adapter through graph) so they stay directly comparable
to the source — see ``config/controlled_vocabularies.yaml`` for the
declared value lists.
"""

from __future__ import annotations

import collections
from typing import Any


def _strip_version(acc: str | None) -> str | None:
    """Strip a trailing ``.N`` version suffix from a signature accession
    (``NF002735.2`` -> ``NF002735``). Accessions without a dot (or ``None``)
    pass through unchanged."""
    return acc.split(".")[0] if acc else acc


def parse_interproscan_json(data: dict) -> dict[str, dict]:
    """Parse an InterProScan JSON document into a WP_-keyed faceted calls dict.

    ``data`` is the loaded ``--formats JSON`` output: ``{"results": [...]}``.
    Each result has ``xref`` (list of ``{id, name}`` — InterProScan dedups
    identical sequences, so one result may back several accessions), ``md5``,
    and ``matches``. We fan every result's facets out to every xref id (the
    FASTA first-token = RefSeq WP_ accession).

    Proteins processed but matching nothing appear in the input results with
    an empty ``matches`` array and are kept here with ``match_count: 0`` and
    all three facets empty — so a missing key means "not processed" while a
    zero-match entry means "no domain found".
    """
    calls: dict[str, dict] = {}
    for result in data.get("results") or []:
        libraries: dict[str, list[dict]] = {}
        entries: dict[str, dict] = {}          # IPR -> rollup accumulator
        go_terms: dict[str, set[str]] = {}
        n_rows = 0
        for m in result.get("matches") or []:
            sig = m.get("signature") or {}
            library = (sig.get("signatureLibraryRelease") or {}).get("library")
            entry = sig.get("entry") or {}
            ipr = entry.get("accession")
            for loc in m.get("locations") or []:
                n_rows += 1
                row = {
                    "accession": _strip_version(sig.get("accession")),
                    "name": sig.get("description") or sig.get("name"),
                    "ipr": ipr,
                    "start": loc.get("start"), "end": loc.get("end"),
                    "evalue": loc.get("evalue"), "score": loc.get("score"),
                }
                if library:
                    libraries.setdefault(library, []).append(row)
                if ipr:
                    ent = entries.setdefault(ipr, {
                        "type": entry.get("type"), "libraries": set(),
                        "match_count": 0, "start": None, "end": None,
                        "evalue": None, "evalue_library": None,
                    })
                    ent["match_count"] += 1
                    if library:
                        ent["libraries"].add(library)
                    s, e = loc.get("start"), loc.get("end")
                    if s is not None and (ent["start"] is None or s < ent["start"]):
                        ent["start"] = s
                    if e is not None and (ent["end"] is None or e > ent["end"]):
                        ent["end"] = e
                    ev = loc.get("evalue")
                    if ev is not None and (ent["evalue"] is None or ev < ent["evalue"]):
                        ent["evalue"], ent["evalue_library"] = ev, library
                    for x in entry.get("goXRefs") or []:
                        if x.get("id"):
                            go_terms.setdefault(x["id"], set()).add(ipr)
        for lib in libraries:
            libraries[lib].sort(key=lambda r: (r["start"] or 0, r["accession"] or ""))
        record = {
            "md5": result.get("md5"),
            "match_count": n_rows,
            "libraries": libraries,
            "interpro_entries": {
                k: {**v, "libraries": sorted(v["libraries"])}
                for k, v in sorted(entries.items())
            },
            "go_terms": {k: sorted(v) for k, v in sorted(go_terms.items())},
        }
        for xref in result.get("xref") or []:
            if xref.get("id"):
                calls[xref["id"]] = dict(record)
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
    xrefs_requested: bool | None = None,
) -> dict[str, Any]:
    """Build the per-strain ``skill_summary.json`` QC dict from parsed calls."""
    calls_made = len(calls)
    proteins_no_match = sum(1 for c in calls.values() if c["match_count"] == 0)
    total_matches = sum(c["match_count"] for c in calls.values())

    lib_counter: collections.Counter = collections.Counter()
    interpro_integrated = 0
    for c in calls.values():
        for lib, rows in c["libraries"].items():
            lib_counter[lib] += len(rows)
            interpro_integrated += sum(1 for r in rows if r["ipr"])

    # GO xref coverage — the QC that tells you --goterms actually took effect
    # (both counters are 0 when it was omitted).
    go_seen: set[str] = set()
    proteins_with_go = 0
    for c in calls.values():
        if c["go_terms"]:
            proteins_with_go += 1
            go_seen.update(c["go_terms"])

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
    }
    if xrefs_requested is not None:
        summary["xrefs_requested"] = xrefs_requested
    if image_digest is not None:
        summary["image_digest"] = image_digest
    if wallclock_s is not None:
        summary["wallclock_s"] = round(wallclock_s, 1)
    return summary
