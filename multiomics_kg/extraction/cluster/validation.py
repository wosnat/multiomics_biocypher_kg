# multiomics_kg/extraction/cluster/validation.py
"""Stage 3: LLM-as-judge validation against the original PDF + CSV."""
import base64
import io
import json
import logging
import re
import time
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

try:
    from pypdf import PdfReader, PdfWriter
except Exception:
    PdfReader = None
    PdfWriter = None

try:
    import openai as _openai
except Exception:
    _openai = None

from multiomics_kg.extraction.pdf_utils import pdf_pages_to_base64, collect_pdf_files
from multiomics_kg.extraction.cluster.prompts import JUDGE_PROMPT


def _build_pdf_content_parts(paper_dir: Path, main_pdf_path: Path,
                              max_pages: int = 15) -> list[dict]:
    """Build PDF content parts for the validation prompt.

    Uses the same cache as visual path to avoid re-encoding PDFs.
    """
    cache_dir = paper_dir / ".extraction_cache" / "shared"
    cache_path = cache_dir / "pdf_content_parts.json"

    if cache_path.exists():
        with open(cache_path) as f:
            cached = json.load(f)
        logger.info("Validation: loaded %d cached PDF content parts", len(cached["parts"]))
        return cached["parts"]

    # Build from scratch (shouldn't normally happen — visual runs first)
    content_parts = []
    all_pdfs = collect_pdf_files(paper_dir, main_pdf_path)
    total_pages = 0
    for pdf_path in all_pdfs:
        remaining = max_pages - total_pages
        if remaining <= 0:
            break
        for b64 in pdf_pages_to_base64(pdf_path, page_range=(0, remaining - 1)):
            content_parts.append({
                "type": "file",
                "file": {"filename": pdf_path.name,
                         "file_data": f"data:application/pdf;base64,{b64}"},
            })
            total_pages += 1
            if total_pages >= max_pages:
                break

    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump({"parts": content_parts, "total_pages": total_pages}, f)

    return content_parts


def _build_csv_summary(csv_path: Path) -> str:
    """Build a text summary of the cluster CSV."""
    try:
        df = pd.read_csv(csv_path)
        summary = f"Cluster CSV ({len(df)} rows):\n"
        summary += f"Columns: {list(df.columns)}\n"
        summary += f"Sample rows:\n{df.head(5).to_string()}\n"
        summary += f"Cluster value counts:\n{df.iloc[:, 1].value_counts().to_string()}"
        return summary
    except Exception as e:
        logger.warning("Could not read CSV for validation: %s", e)
        return ""


def _validate_single_cluster(
    cluster_key: str,
    stage2_entry: dict,
    pdf_parts: list[dict],
    csv_summary: str,
    model: str,
    analysis_name: str = "",
) -> dict:
    """Validate one cluster's description against the paper."""
    description = (
        f"Cluster {cluster_key}: id={stage2_entry.get('id', '?')}, "
        f"functional={stage2_entry.get('functional_description', '(sentinel)')}, "
        f"behavioral={stage2_entry.get('behavioral_description', '(sentinel)')}"
    )

    content_parts = list(pdf_parts)  # shallow copy — PDF pages shared
    if csv_summary:
        content_parts.append({"type": "text", "text": csv_summary})
    content_parts.append({
        "type": "text",
        "text": JUDGE_PROMPT.format(descriptions=description,
                                     analysis_name=analysis_name),
    })

    client = _openai.OpenAI()
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": content_parts}],
            temperature=0,
        )
        raw = response.choices[0].message.content.strip()
    except Exception as e:
        logger.error("Validation failed for cluster %s: %s", cluster_key, e)
        return {"verdict": "fail", "explanation": f"API error: {e}"}

    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    if not raw.strip():
        logger.error("Validation returned empty response for cluster %s", cluster_key)
        return {"verdict": "fail", "explanation": "Empty response from validator"}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Validation JSON parse failed for cluster %s: %s\nRaw response: %s",
                     cluster_key, e, raw[:500])
        return {"verdict": "fail", "explanation": f"JSON parse error: {e}"}

    # Extract the cluster entry — LLM may wrap in a dict with the key
    if isinstance(parsed, dict):
        if cluster_key in parsed and isinstance(parsed[cluster_key], dict):
            return parsed[cluster_key]
        # Try to find a single entry
        entries = [v for v in parsed.values() if isinstance(v, dict)]
        if len(entries) == 1:
            return entries[0]
        # If the dict itself has verdict, it's the direct entry
        if "verdict" in parsed:
            return parsed
    return {"verdict": "fail", "explanation": "Could not parse validation response"}


def run_validation(main_pdf_path: Path,
                   paper_dir: Path,
                   csv_path: Path,
                   stage2_results: dict[str, dict],
                   cluster_keys: list[str],
                   model: str = "gpt-4o",
                   max_pages: int = 15,
                   analysis_name: str = "",
                   ) -> dict[str, dict]:
    """Validate Stage 2 descriptions against original PDF + CSV data.

    Validates one cluster at a time to prevent cross-contamination.
    """
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    # Build shared context once (PDF pages + CSV summary)
    pdf_parts = _build_pdf_content_parts(paper_dir, main_pdf_path, max_pages)
    csv_summary = _build_csv_summary(csv_path)

    results: dict[str, dict] = {}
    for key in cluster_keys:
        s2 = stage2_results.get(key, {})
        if not s2:
            logger.warning("No stage2 data for cluster %s, skipping validation", key)
            results[key] = {"verdict": "fail", "explanation": "No synthesis data"}
            continue

        logger.info("Validating cluster %s...", key)
        result = _validate_single_cluster(key, s2, pdf_parts, csv_summary, model,
                                          analysis_name=analysis_name)
        results[key] = result
        verdict = result.get("verdict", "?")
        logger.info("  Cluster %s: %s", key, verdict)

        # Small delay between calls to avoid rate limiting
        time.sleep(2)

    # Check completeness
    missing = set(cluster_keys) - set(results.keys())
    if missing:
        logger.warning("Validation missing clusters: %s", missing)
        for k in missing:
            results[k] = {"verdict": "fail", "explanation": "Validation not run"}

    return results
