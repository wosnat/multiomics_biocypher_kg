"""Aggregate per-paper extraction JSONs into data/timepoint_extraction_report.md."""
from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path

from multiomics_kg.extraction.timepoint.extraction_utils import (
    iter_paperconfigs,
    load_extraction_json,
)

logger = logging.getLogger(__name__)

REPORT_PATH = Path("data/timepoint_extraction_report.md")


def aggregate_reports(paper_dirs: list[Path]) -> dict:
    """Walk paper_dirs, aggregate stats from each extraction JSON."""
    sa_counter: Counter = Counter()
    other_counter: Counter = Counter()
    unknowns: list[dict] = []
    partial: list[dict] = []
    per_paper: list[dict] = []

    for paper_dir in paper_dirs:
        ext = load_extraction_json(paper_dir)
        if ext is None:
            continue
        paper = ext["metadata"].get("paper", paper_dir.name)
        status = ext["metadata"].get("status")
        rows = ext.get("analyses", [])

        if status == "partial":
            partial.append({
                "paper": paper,
                "missing_analyses": ext["metadata"].get("missing_analyses", []),
            })

        for row in rows:
            sa = row.get("self_assessment")
            if sa:
                sa_counter[sa] += 1
            gp = row.get("growth_phase")
            if isinstance(gp, str) and gp.startswith("other:"):
                other_counter[gp] += 1
            if gp == "unknown":
                unknowns.append({
                    "paper": paper,
                    "analysis_id": row.get("analysis_id"),
                    "notes": row.get("assessment_notes", ""),
                    "quotes": row.get("supporting_quotes", []),
                })
        per_paper.append({"paper": paper, "n_analyses": len(rows), "status": status})

    return {
        "self_assessment_counts": dict(sa_counter),
        "other_slug_counts": dict(other_counter),
        "unknowns": unknowns,
        "partial_extractions": partial,
        "per_paper": per_paper,
    }


def render_markdown(agg: dict) -> str:
    lines = ["# Timepoint Extraction Report", ""]

    lines.append("## Self-assessment distribution")
    for sa in ("high", "medium", "low"):
        n = agg["self_assessment_counts"].get(sa, 0)
        lines.append(f"- `{sa}`: {n}")
    lines.append("")

    lines.append("## other:* slugs")
    sorted_slugs = sorted(
        agg["other_slug_counts"].items(), key=lambda kv: kv[1], reverse=True,
    )
    if not sorted_slugs:
        lines.append("_(none)_")
    else:
        for slug, n in sorted_slugs:
            lines.append(f"- `{slug}`: {n}")
    lines.append("")

    lines.append("## Unknowns")
    if not agg["unknowns"]:
        lines.append("_(none)_")
    else:
        for u in agg["unknowns"]:
            lines.append(f"- **{u['paper']}** / `{u['analysis_id']}`: {u['notes']}")
            for q in u["quotes"]:
                lines.append(f"  - _{q.get('location', '')}_: {q.get('quote', '')}")
    lines.append("")

    lines.append("## Partial extractions")
    if not agg["partial_extractions"]:
        lines.append("_(none)_")
    else:
        for p in agg["partial_extractions"]:
            lines.append(f"- **{p['paper']}**: {len(p['missing_analyses'])} missing")
            for m in p["missing_analyses"]:
                lines.append(f"  - `{m['analysis_id']}`: {m['reason']}")
    lines.append("")

    lines.append("## Per-paper summary")
    for p in agg["per_paper"]:
        lines.append(f"- {p['paper']}: {p['n_analyses']} analyses ({p['status']})")

    return "\n".join(lines) + "\n"


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="report.py")
    parser.add_argument("--output", default=str(REPORT_PATH),
                        help=f"Output path (default: {REPORT_PATH})")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]

    agg = aggregate_reports(paper_dirs)
    md = render_markdown(agg)
    Path(args.output).write_text(md)
    logger.info("Report written: %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
