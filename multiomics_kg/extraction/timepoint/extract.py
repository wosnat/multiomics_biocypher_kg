"""Timepoint/growth-phase LLM extraction per paper.

Usage:
    uv run python -m multiomics_kg.extraction.timepoint.extract --paper "Tetu 2019"
    uv run python -m multiomics_kg.extraction.timepoint.extract --all
    uv run python -m multiomics_kg.extraction.timepoint.extract --paper "X" --dry-run
    uv run python -m multiomics_kg.extraction.timepoint.extract --all --validate
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Protocol

from multiomics_kg.extraction.timepoint.extraction_utils import (
    _load_yaml,
    compute_fields_requested,
    compute_paperconfig_signature,
    find_analyses,
    iter_paperconfigs,
    save_extraction_json,
    validate_llm_payload,
)
from multiomics_kg.extraction.timepoint.prompts import build_prompt

logger = logging.getLogger(__name__)

_DEFAULT_PDF_CACHE = (
    Path(__file__).resolve().parent.parent.parent.parent
    / "cache"
    / "pdf_extraction_cache.json"
)


class LLMClient(Protocol):
    """Minimal protocol — any caller providing `.call(...)` satisfies this."""
    def call(
        self, prompt: str, pdf_paths: list[Path], model: str | None = None,
    ) -> tuple[dict, dict]:
        """Return (payload_dict, metadata_dict).
        metadata_dict: keys `input_tokens`, `output_tokens`, `model`.
        """
        ...


def build_background(paperconfig_path: Path) -> dict:
    """Extract the background context for the LLM prompt."""
    data = _load_yaml(paperconfig_path)
    pub = data.get("publication", {})
    return {
        "papername": pub.get("papername", ""),
        "doi": pub.get("doi", ""),
        "experiments": dict(pub.get("experiments") or {}),
        "additional_pdfs": list(
            (pub.get("extraction") or {}).get("additional_pdfs") or []
        ),
        "papermainpdf": pub.get("papermainpdf"),
    }


def build_targets(analyses: list[dict], validate: bool) -> list[dict]:
    """Build the `targets` list for the prompt. Omits analyses with empty
    fields_requested unless --validate.
    """
    targets: list[dict] = []
    for a in analyses:
        requested = compute_fields_requested(a, validate=validate)
        if not requested and not validate:
            continue
        targets.append({
            "id": a["id"],
            "experiment_key": a.get("experiment", ""),
            "logfc_col": a.get("logfc_col", ""),
            "existing": {
                "timepoint": a.get("timepoint"),
                "timepoint_hours": a.get("timepoint_hours"),
                "growth_phase": a.get("growth_phase"),
            },
            "fields_requested": requested,
        })
    return targets


def _load_pdf_cache_entry(cache_path: Path, papermainpdf: str | None) -> dict | None:
    if not papermainpdf or not cache_path.exists():
        return None
    try:
        cache = json.loads(cache_path.read_text())
    except Exception:
        return None
    entry = cache.get(papermainpdf)
    if entry and "publication" in entry:
        return entry["publication"]
    return None


def extract_one_paper(
    paper_dir: Path,
    llm_client: LLMClient,
    validate: bool = False,
    pdf_cache_path: Path | None = None,
    model: str | None = None,
) -> Path | None:
    """Extract one paper's timepoint/growth_phase metadata via the LLM.

    Returns the written JSON path, or None if no extraction was needed
    (all fields populated, default mode).
    """
    paperconfig_path = paper_dir / "paperconfig.yaml"
    analyses = find_analyses(paperconfig_path)
    targets = build_targets(analyses, validate=validate)
    if not targets:
        logger.info("No fields to extract for %s — skipping LLM call.", paper_dir)
        return None

    background = build_background(paperconfig_path)

    pdf_cache_entry = _load_pdf_cache_entry(
        pdf_cache_path or _DEFAULT_PDF_CACHE,
        background.get("papermainpdf"),
    )
    prompt = build_prompt(background, targets, pdf_cache_entry)

    pdf_paths: list[Path] = []
    if pmp := background.get("papermainpdf"):
        pdf_paths.append(Path(pmp))
    pdf_paths.extend(Path(p) for p in background.get("additional_pdfs", []))

    payload, call_meta = llm_client.call(prompt, pdf_paths, model=model)

    requested_map = {t["id"]: t["fields_requested"] for t in targets}
    valid_rows, missing = validate_llm_payload(payload, requested_map)

    # Inject experiment_key and fields_requested
    for row in valid_rows:
        aid = row["analysis_id"]
        # Look up experiment_key from targets
        for t in targets:
            if t["id"] == aid:
                row["experiment_key"] = t["experiment_key"]
                row["fields_requested"] = t["fields_requested"]
                break

    metadata = {
        "paper": background.get("papername", ""),
        "doi": background.get("doi", ""),
        "model": call_meta.get("model", ""),
        "extracted_at": datetime.utcnow().isoformat(timespec="seconds"),
        "input_tokens": call_meta.get("input_tokens", 0),
        "output_tokens": call_meta.get("output_tokens", 0),
        "paperconfig_signature": compute_paperconfig_signature(paperconfig_path),
        "status": "complete" if not missing else "partial",
        "missing_analyses": missing,
    }

    json_path = save_extraction_json(paper_dir, metadata, valid_rows)
    return json_path


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="extract.py")
    parser.add_argument("--paper", help="Paper name (papername) to extract")
    parser.add_argument("--all", action="store_true", help="Extract all papers")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print fields_requested per analysis; no LLM call")
    parser.add_argument("--validate", action="store_true",
                        help="Re-examine every field, even populated ones")
    parser.add_argument("--model", default=None,
                        help="Override default model (e.g. gpt-4.1-mini)")
    # --resume and --retry are added in Task 12.
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    paper_dirs = _resolve_paper_dirs(args)
    if not paper_dirs:
        logger.error("No papers matched selection.")
        return 2

    if args.dry_run:
        for pd in paper_dirs:
            _print_dry_run(pd, validate=args.validate)
        return 0

    from multiomics_kg.extraction.timepoint.llm_client import OpenAIResponsesClient
    client = OpenAIResponsesClient()
    had_partial = False
    for pd in paper_dirs:
        try:
            path = extract_one_paper(pd, client, validate=args.validate, model=args.model)
            if path is None:
                logger.info("[%s] no fields to extract, skipped.", pd.name)
                continue
            data = json.loads(path.read_text())
            logger.info("[%s] status=%s wrote %s",
                        pd.name, data["metadata"]["status"], path)
            if data["metadata"]["status"] == "partial":
                had_partial = True
        except Exception:
            logger.exception("[%s] extraction failed", pd.name)
            had_partial = True
    return 1 if had_partial else 0


def _resolve_paper_dirs(args) -> list[Path]:
    """Map --paper / --all to paper directories. See implementation notes."""
    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]
    if args.all:
        return paper_dirs
    if args.paper:
        matches = []
        for pd in paper_dirs:
            try:
                data = _load_yaml(pd / "paperconfig.yaml")
                if data.get("publication", {}).get("papername") == args.paper:
                    matches.append(pd)
            except Exception:
                logger.warning("Could not parse %s/paperconfig.yaml", pd)
                continue
        return matches
    return []


def _print_dry_run(paper_dir: Path, validate: bool) -> None:
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    targets = build_targets(analyses, validate=validate)
    print(f"\n{paper_dir.name}:")
    if not targets:
        print("  (all fields populated)")
        return
    for t in targets:
        print(f"  {t['id']}: needs {t['fields_requested']}  (logfc_col={t['logfc_col']!r})")


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
