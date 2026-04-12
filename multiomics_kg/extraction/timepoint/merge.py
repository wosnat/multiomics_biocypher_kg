"""Write extraction JSON fields back into paperconfig.yaml.

Rules:
- Verbatim merge: growth_phase written exactly as it appears in JSON
  (no other:<slug> → canonical transformation; that's remap.py's job).
- Staleness check (paperconfig_signature + analysis_id cross-check)
  refuses to write without --force.
- Partial JSON (metadata.status == "partial") requires --force.
- Overwrite guard: differing existing non-null value requires --force.
- Malformed values rejected even with --force.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    TARGET_FIELDS,
    _is_valid_growth_phase,
    compute_paperconfig_signature,
    load_extraction_json,
)

logger = logging.getLogger(__name__)


def merge_one_paper(paper_dir: Path, force: bool = False) -> None:
    """Merge this paper's extraction JSON into its paperconfig.yaml.

    On any refusal, prints the reason and sys.exit(1). --force overrides
    partial / signature-mismatch / overwrite-conflict, but never malformed
    values.
    """
    paperconfig_path = paper_dir / "paperconfig.yaml"
    extraction = load_extraction_json(paper_dir)
    if extraction is None:
        logger.error("No extraction JSON at %s", paper_dir)
        sys.exit(1)

    metadata = extraction.get("metadata", {})
    rows = extraction.get("analyses", [])

    # 1. Partial guard
    if metadata.get("status") == "partial" and not force:
        logger.error(
            "Refusing to merge %s: extraction is partial. Missing: %s. "
            "Pass --force to merge the valid rows only.",
            paper_dir.name, metadata.get("missing_analyses", []),
        )
        sys.exit(1)

    # 2. Signature check
    current_sig = compute_paperconfig_signature(paperconfig_path)
    stored_sig = metadata.get("paperconfig_signature")
    if stored_sig and stored_sig != current_sig and not force:
        logger.error(
            "Refusing to merge %s: paperconfig_signature mismatch. "
            "Paperconfig changed since extraction (stored=%s, now=%s). "
            "Pass --force to merge anyway.",
            paper_dir.name, stored_sig[:8], current_sig[:8],
        )
        sys.exit(1)

    # 3. Reject malformed growth_phase values (force does not bypass this)
    for row in rows:
        gp = row.get("growth_phase")
        if gp is not None and not _is_valid_growth_phase(gp):
            logger.error(
                "Refusing to merge %s: analysis %s has malformed growth_phase %r.",
                paper_dir.name, row.get("analysis_id"), gp,
            )
            sys.exit(1)

    yaml = YAML()
    yaml.preserve_quotes = True
    with open(paperconfig_path) as f:
        data = yaml.load(f)

    analyses_by_id: dict = {}
    for table in (data.get("publication", {}).get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            analyses_by_id[a["id"]] = a

    # 4. Analysis-id cross-check
    ghosts = [r["analysis_id"] for r in rows if r["analysis_id"] not in analyses_by_id]
    for g in ghosts:
        logger.warning("Ghost analysis_id in JSON (not in paperconfig): %s", g)

    coverage_gap = [aid for aid in analyses_by_id if aid not in {r["analysis_id"] for r in rows}]
    if coverage_gap:
        logger.warning(
            "Paperconfig has analyses not covered by extraction (re-run extract?): %s",
            coverage_gap,
        )

    # 5. Overwrite guard + apply
    for row in rows:
        aid = row["analysis_id"]
        if aid in ghosts:
            continue
        target = analyses_by_id[aid]
        for field in TARGET_FIELDS:
            if field not in row:
                continue
            new_value = row[field]
            existing = target.get(field)
            if existing not in (None, "") and existing != new_value and not force:
                logger.error(
                    "Refusing to overwrite %s:%s (current=%r, json=%r). Pass --force.",
                    aid, field, existing, new_value,
                )
                sys.exit(1)
            target[field] = new_value

    with open(paperconfig_path, "w") as f:
        yaml.dump(data, f)

    logger.info("Merged %s → %s", paper_dir.name, paperconfig_path)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="merge.py")
    parser.add_argument("--paper", help="Paper name (papername) to merge")
    parser.add_argument("--all", action="store_true", help="Merge all papers that have extraction JSON")
    parser.add_argument("--force", action="store_true", help="Allow overwrite / partial / signature override")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from multiomics_kg.extraction.timepoint.extract import _resolve_paper_dirs
    paper_dirs = _resolve_paper_dirs(args)
    if not paper_dirs:
        logger.error("No papers matched selection.")
        return 2

    for pd in paper_dirs:
        if (pd / "extractions" / "timepoint.json").exists():
            merge_one_paper(pd, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
