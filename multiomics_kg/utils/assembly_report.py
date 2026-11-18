"""NCBI assembly-quality report for a genome strain.

Step 0 (``download_genome_data._ncbi_fetch_assembly_reports``) writes one
``assembly_report.json`` per strain directory, shaped
``{"fetched_at": "YYYY-MM-DD", "report": <NCBI Datasets v2 dataset_report>}``.
"""
from __future__ import annotations

import json
import os
from typing import Any

from multiomics_kg.utils.controlled_vocab import VOCAB

ASSEMBLY_REPORT_FILENAME = "assembly_report.json"

_SUPPRESSED_FALLBACK = "removed by RefSeq staff"


def load_assembly_report(data_dir: str) -> dict | None:
    """Return the strain's wrapped report, or None when no file exists."""
    path = os.path.join(data_dir, ASSEMBLY_REPORT_FILENAME)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _suppression_reason(text: str) -> str:
    if "Reason:" in text:
        reason = text.split("Reason:", 1)[1].strip().rstrip(".").strip()
        if reason:
            return reason
    return _SUPPRESSED_FALLBACK


def assembly_properties(wrapped: dict, accession: str) -> dict[str, Any]:
    """Map a wrapped NCBI assembly report onto OrganismTaxon properties.

    Absent stays absent (never "" or None), except checkm_contamination: NCBI
    omits zero-valued fields, so a CheckM block without it reads 0.0.
    assembly_level / assembly_status outside the declared vocabulary raise.
    """
    report = wrapped.get("report") or {}
    info = report.get("assembly_info") or {}
    stats = report.get("assembly_stats") or {}
    checkm = report.get("checkm_info")
    props: dict[str, Any] = {}

    if info.get("assembly_level"):
        level = info["assembly_level"].lower().replace(" ", "_")
        props["assembly_level"] = VOCAB.check("OrganismTaxon", "assembly_level", level)
    status = info.get("assembly_status")
    if status:
        props["assembly_status"] = VOCAB.check("OrganismTaxon", "assembly_status", status)
        if status == "suppressed":
            props["assembly_suppression_reason"] = _suppression_reason(
                info.get("suppression_reason") or "")
    current = report.get("current_accession")
    if current and current != accession:
        props["assembly_current_accession"] = current
    for prop, key in (("assembly_name", "assembly_name"),
                      ("assembly_release_date", "release_date")):
        if info.get(key):
            props[prop] = info[key]
    if wrapped.get("fetched_at"):
        props["assembly_report_date"] = wrapped["fetched_at"]

    for prop, key in (("genome_length", "total_sequence_length"),
                      ("contig_count", "number_of_contigs"),
                      ("contig_n50", "contig_n50")):
        if stats.get(key) not in (None, ""):
            props[prop] = int(stats[key])
    if stats.get("gc_percent") is not None:
        props["gc_percent"] = float(stats["gc_percent"])

    if checkm:
        if checkm.get("completeness") is not None:
            props["checkm_completeness"] = float(checkm["completeness"])
        props["checkm_contamination"] = float(checkm.get("contamination", 0.0))
        for prop in ("checkm_marker_set", "checkm_marker_set_rank", "checkm_version"):
            if checkm.get(prop):
                props[prop] = checkm[prop]
    return props
