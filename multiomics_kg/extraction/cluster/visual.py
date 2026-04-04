# multiomics_kg/extraction/cluster/visual.py
"""Path: visual — Send PDF pages to gpt-4o for cluster extraction."""
import json
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)

try:
    import openai as _openai
except Exception:
    _openai = None

from multiomics_kg.extraction.pdf_utils import pdf_pages_to_base64, collect_pdf_files
from multiomics_kg.extraction.cluster.prompts import VISUAL_PROMPT, EXTRACTION_FIELDS_DESCRIPTION


import time


def _build_pdf_parts(paper_dir: Path, main_pdf_path: Path,
                     max_pages: int = 15) -> tuple[list[dict], int]:
    """Build PDF content parts once, shared across per-cluster calls."""
    all_pdfs = collect_pdf_files(paper_dir, main_pdf_path)
    content_parts = []
    total_pages = 0
    for pdf_path in all_pdfs:
        remaining = max_pages - total_pages
        if remaining <= 0:
            break
        pages_b64 = pdf_pages_to_base64(pdf_path, page_range=(0, remaining - 1))
        for b64 in pages_b64:
            content_parts.append({
                "type": "file",
                "file": {
                    "filename": pdf_path.name,
                    "file_data": f"data:application/pdf;base64,{b64}",
                },
            })
            total_pages += 1
            if total_pages >= max_pages:
                break
    return content_parts, total_pages


def _extract_single_cluster(
    cluster_key: str,
    pdf_parts: list[dict],
    cluster_keys: list[str],
    cluster_method: str,
    organism: str,
    analysis_name: str,
    model: str,
) -> dict:
    """Extract visual data for one cluster."""
    system_text = VISUAL_PROMPT.format(
        n_clusters=len(cluster_keys),
        cluster_method=cluster_method,
        organism=organism,
        cluster_key=cluster_key,
        analysis_name=analysis_name,
        fields_description=EXTRACTION_FIELDS_DESCRIPTION,
    )

    content_parts = list(pdf_parts)  # shallow copy — same PDF pages
    content_parts.append({
        "type": "text",
        "text": (f"Extract data for {organism} cluster {cluster_key} ONLY. "
                 f"Return ONLY valid JSON."),
    })

    client = _openai.OpenAI()
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_text},
                {"role": "user", "content": content_parts},
            ],
            temperature=0,
        )
        raw = response.choices[0].message.content.strip()
    except Exception as e:
        logger.error("Visual extraction failed for cluster %s: %s", cluster_key, e)
        return {}

    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    if not raw.strip():
        logger.error("Visual returned empty response for cluster %s", cluster_key)
        return {}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Visual JSON parse failed for cluster %s: %s\nRaw: %s",
                     cluster_key, e, raw[:500])
        return {}

    if not isinstance(parsed, dict):
        return {}

    # Extract cluster entry — LLM may wrap in a dict with the key
    if "enrichment_category" in parsed or "direction" in parsed:
        result = parsed  # flat response
    elif cluster_key in parsed and isinstance(parsed[cluster_key], dict):
        result = parsed[cluster_key]
    else:
        entries = [v for v in parsed.values() if isinstance(v, dict)]
        result = entries[0] if len(entries) == 1 else parsed

    result["source"] = "visual"
    result["confidence"] = "high"
    return result


def run_visual(main_pdf_path: Path,
               paper_dir: Path,
               cluster_keys: list[str],
               cluster_method: str,
               organism: str,
               treatment: str,
               model: str = "gpt-4o",
               max_pages: int = 15,
               analysis_name: str = "",
               ) -> dict[str, dict]:
    """Run visual extraction — one LLM call per cluster with shared PDF pages."""
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    pdf_parts, total_pages = _build_pdf_parts(paper_dir, main_pdf_path, max_pages)
    logger.info("Path visual: prepared %d pages from PDFs", total_pages)

    if not analysis_name:
        analysis_name = f"{organism} {cluster_method}"

    results: dict[str, dict] = {}
    for key in cluster_keys:
        logger.info("  Visual extracting cluster %s...", key)
        result = _extract_single_cluster(
            key, pdf_parts, cluster_keys, cluster_method,
            organism, analysis_name, model
        )
        if result:
            results[key] = result
        time.sleep(2)  # rate limit spacing

    return results


