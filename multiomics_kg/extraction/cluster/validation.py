# multiomics_kg/extraction/cluster/validation.py
"""Stage 3: LLM-as-judge validation against the original PDF + CSV."""
import base64
import io
import json
import logging
import re
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


def run_validation(main_pdf_path: Path,
                   paper_dir: Path,
                   csv_path: Path,
                   stage2_results: dict[str, dict],
                   cluster_keys: list[str],
                   model: str = "gpt-4o",
                   max_pages: int = 15,
                   ) -> dict[str, dict]:
    """Validate Stage 2 descriptions against original PDF + CSV data."""
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    desc_lines = []
    for key in cluster_keys:
        s2 = stage2_results.get(key, {})
        desc_lines.append(
            f"Cluster {key}: id={s2.get('id', '?')}, "
            f"functional={s2.get('functional_description', '(sentinel)')}, "
            f"behavioral={s2.get('behavioral_description', '(sentinel)')}"
        )
    descriptions = "\n".join(desc_lines)

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

    try:
        df = pd.read_csv(csv_path)
        csv_summary = f"Cluster CSV ({len(df)} rows):\n"
        csv_summary += f"Columns: {list(df.columns)}\n"
        csv_summary += f"Sample rows:\n{df.head(5).to_string()}\n"
        csv_summary += f"Cluster value counts:\n{df.iloc[:, 1].value_counts().to_string()}"
        content_parts.append({"type": "text", "text": csv_summary})
    except Exception as e:
        logger.warning("Could not read CSV for validation: %s", e)

    content_parts.append({
        "type": "text",
        "text": JUDGE_PROMPT.format(descriptions=descriptions),
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
        logger.error("Stage 3 validation failed: %s", e)
        return {}

    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Stage 3 JSON parse failed: %s", e)
        return {}

    valid_keys = set(cluster_keys)
    results: dict[str, dict] = {}
    for k, v in parsed.items():
        num = re.search(r'\d+', str(k))
        norm = num.group() if num and num.group() in valid_keys else str(k)
        if norm in valid_keys and isinstance(v, dict):
            results[norm] = v
    return results
