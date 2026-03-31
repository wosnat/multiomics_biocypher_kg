# multiomics_kg/extraction/cluster/visual.py
"""Path: visual — Send PDF pages to gpt-4o for cluster extraction."""
import json
import logging
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

try:
    import openai as _openai
except Exception:
    _openai = None

from multiomics_kg.extraction.pdf_utils import pdf_pages_to_base64, collect_pdf_files
from multiomics_kg.extraction.cluster.prompts import VISUAL_PROMPT, EXTRACTION_FIELDS_DESCRIPTION


def run_visual(main_pdf_path: Path,
               paper_dir: Path,
               cluster_keys: list[str],
               cluster_method: str,
               organism: str,
               treatment: str,
               model: str = "gpt-4o",
               max_pages: int = 15,
               ) -> dict[str, dict]:
    if _openai is None:
        logger.error("openai package not installed")
        return {}
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
    logger.info("Path visual: sent %d pages from %d PDFs", total_pages, len(all_pdfs))
    system_text = VISUAL_PROMPT.format(
        n_clusters=len(cluster_keys),
        cluster_method=cluster_method,
        organism=organism,
        treatment=treatment,
        cluster_keys=json.dumps(cluster_keys),
        fields_description=EXTRACTION_FIELDS_DESCRIPTION,
    )
    content_parts.append({
        "type": "text",
        "text": (f"Extract cluster data for {organism}, clusters "
                 f"{json.dumps(cluster_keys)}. Return ONLY valid JSON."),
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
        logger.error("Path visual API call failed: %s", e)
        return {}
    return _parse_response(raw, cluster_keys)


def _parse_response(raw: str, cluster_keys: list[str]) -> dict[str, dict]:
    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Path visual JSON parse failed: %s\nRaw: %s", e, raw[:500])
        return {}
    valid_keys = set(cluster_keys)
    results: dict[str, dict] = {}
    if isinstance(parsed, dict):
        for k, v in parsed.items():
            norm = _normalize_key(str(k), valid_keys)
            if norm and isinstance(v, dict):
                v["source"] = "visual"
                v["confidence"] = "high"
                results[norm] = v
    return results


def _normalize_key(raw_key: str, valid_keys: set[str]) -> Optional[str]:
    if raw_key in valid_keys:
        return raw_key
    num = re.search(r'\d+', raw_key)
    if num and num.group() in valid_keys:
        return num.group()
    for vk in valid_keys:
        if vk.lower() == raw_key.lower():
            return vk
    logger.warning("Could not normalize key %r", raw_key)
    return None
