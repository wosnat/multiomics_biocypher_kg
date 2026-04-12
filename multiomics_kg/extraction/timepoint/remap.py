"""Remap growth_phase values across paperconfigs + extraction JSONs.

Usage:
    # Promote other:<slug> → canonical (after adding to VALID_GROWTH_PHASES):
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_acclimated --to heat_acclimated

    # Rename a slug:
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_acclim --to other:heat_acclimated

    # Merge into an existing canonical value:
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_shocked --to acute_stress
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    EXTRACTION_FILENAME,
    EXTRACTIONS_DIR,
    _VALID_GROWTH_PHASES,
    iter_paperconfigs,
)

logger = logging.getLogger(__name__)

# A valid remap target is either a simple slug (lowercase alphanumeric + underscores,
# no spaces) or an other:<slug> form.  We intentionally do NOT require the value to
# already be in _VALID_GROWTH_PHASES so that callers can promote a slug *before*
# updating the validator (they should add it to the frozenset afterwards).
_VALID_TARGET_RE = re.compile(r"^(?:other:)?[a-z][a-z0-9_]*$")


def _is_valid_remap_target(v: str) -> bool:
    """Accept canonical slugs and other:<slug> forms; reject anything with spaces etc."""
    return bool(_VALID_TARGET_RE.match(v))


def remap_value(
    paper_dirs: list[Path],
    from_value: str,
    to_value: str,
    dry_run: bool = False,
) -> None:
    """Rewrite growth_phase=from_value → growth_phase=to_value across
    paperconfigs and extraction JSONs in paper_dirs.
    """
    if not _is_valid_remap_target(to_value):
        logger.error(
            "Refusing to remap: --to value %r is not a valid growth_phase "
            "(must be a lowercase slug or other:<slug>, no spaces). "
            "If promoting, update VALID_GROWTH_PHASES in the validator first.",
            to_value,
        )
        sys.exit(1)

    if not to_value.startswith("other:") and to_value not in _VALID_GROWTH_PHASES:
        logger.warning(
            "--to value %r is not in VALID_GROWTH_PHASES. Merge will reject it "
            "until you add it to validate_paperconfig.py. Proceeding on the "
            "assumption that you are about to add it.",
            to_value,
        )

    yaml = YAML()
    yaml.preserve_quotes = True
    provenance = f"Remapped from {from_value} \u2192 {to_value} on {date.today().isoformat()}"

    for paper_dir in paper_dirs:
        _remap_paperconfig(paper_dir, from_value, to_value, yaml, dry_run)
        _remap_json(paper_dir, from_value, to_value, provenance, dry_run)


def _remap_paperconfig(
    paper_dir: Path, from_value: str, to_value: str, yaml: YAML, dry_run: bool,
) -> None:
    path = paper_dir / "paperconfig.yaml"
    if not path.exists():
        return
    with open(path) as f:
        data = yaml.load(f)
    changed = False
    for table in (data.get("publication", {}).get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            if a.get("growth_phase") == from_value:
                a["growth_phase"] = to_value
                changed = True
    if changed and not dry_run:
        with open(path, "w") as f:
            yaml.dump(data, f)
        logger.info("Remapped %s", path)
    elif changed:
        logger.info("[dry-run] would remap %s", path)


def _remap_json(
    paper_dir: Path, from_value: str, to_value: str, provenance: str, dry_run: bool,
) -> None:
    path = paper_dir / EXTRACTIONS_DIR / EXTRACTION_FILENAME
    if not path.exists():
        return
    data = json.loads(path.read_text())
    changed = False
    for row in data.get("analyses", []):
        if row.get("growth_phase") == from_value:
            row["growth_phase"] = to_value
            existing_note = row.get("assessment_notes", "") or ""
            row["assessment_notes"] = (
                f"{existing_note}\n{provenance}" if existing_note else provenance
            )
            changed = True
    if changed and not dry_run:
        path.write_text(json.dumps(data, indent=2, default=str))
        logger.info("Remapped %s", path)
    elif changed:
        logger.info("[dry-run] would remap %s", path)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="remap.py")
    parser.add_argument("--from", dest="from_value", required=True,
                        help="Current value (e.g. other:heat_acclimated)")
    parser.add_argument("--to", dest="to_value", required=True,
                        help="Target value (canonical enum or other:<slug>)")
    parser.add_argument("--dry-run", action="store_true",
                        help="List affected files without writing")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]
    remap_value(paper_dirs, args.from_value, args.to_value, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
