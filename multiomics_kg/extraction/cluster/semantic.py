# multiomics_kg/extraction/cluster/semantic.py
"""Path: semantic — Embed paper text, retrieve per-cluster, extract with LLM."""
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

from multiomics_kg.extraction.pdf_utils import extract_pdf_text, collect_pdf_files
from multiomics_kg.extraction.rag import chunk_text, embed_texts, retrieve_top_k
from multiomics_kg.extraction.structured_parsers import scan_supplementary_text
from multiomics_kg.extraction.cluster.prompts import SEMANTIC_PROMPT, EXTRACTION_FIELDS_DESCRIPTION


def build_cluster_query(cluster_key: str, seed_data: dict) -> str:
    """Build a semantic search query seeded by visual/table results."""
    parts = [f"cluster {cluster_key}"]
    if seed_data.get("enrichment_category"):
        val = seed_data["enrichment_category"]
        if isinstance(val, list):
            parts.append(str(val[0].get("value", "")))
        else:
            parts.append(str(val))
    if seed_data.get("direction"):
        val = seed_data["direction"]
        if isinstance(val, list):
            d = str(val[0].get("value", ""))
        else:
            d = str(val)
        if d:
            parts.append(f"{d}regulated")
    if seed_data.get("genes"):
        val = seed_data["genes"]
        if isinstance(val, list) and val:
            g = val[0].get("value", "")
            if isinstance(g, list):
                parts.extend(str(x) for x in g[:5])
            else:
                parts.append(str(g))
    return " ".join(parts)


def run_semantic(main_pdf_path: Path,
                 paper_dir: Path,
                 cluster_keys: list[str],
                 organism: str,
                 treatment: str,
                 seed_data: Optional[dict[str, dict]] = None,
                 model: str = "gpt-5-nano",
                 top_k: int = 8,
                 ) -> dict[str, dict]:
    """Run Path: semantic — chunk, embed, retrieve, extract per cluster."""
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    seed_data = seed_data or {}

    all_text_parts = []
    for pdf_path in collect_pdf_files(paper_dir, main_pdf_path):
        text = extract_pdf_text(pdf_path)
        if text:
            all_text_parts.append(text)

    supp = scan_supplementary_text(paper_dir)
    for source_file, text in supp.items():
        all_text_parts.append(f"[{source_file}]\n{text}")

    combined = "\n\n".join(all_text_parts)
    if not combined.strip():
        logger.warning("No text extracted for semantic path")
        return {}

    chunks = chunk_text(combined)
    logger.info("Path semantic: %d chunks from %d chars", len(chunks), len(combined))

    chunk_embs = embed_texts(chunks)
    if not chunk_embs:
        logger.error("Embedding failed")
        return {}

    results: dict[str, dict] = {}
    for key in cluster_keys:
        query = build_cluster_query(key, seed_data.get(key, {}))
        retrieved = retrieve_top_k(query, chunks, chunk_embs, top_k=top_k)
        if not retrieved:
            continue
        passages = "\n\n".join(
            f"[{i+1}] (relevance: {score:.3f}): \"{text}\""
            for i, (text, score) in enumerate(retrieved)
        )
        extraction = _extract_from_passages(key, passages, organism, treatment,
                                            retrieved, model)
        if extraction:
            results[key] = extraction

    return results


def _extract_from_passages(cluster_key: str,
                           passages: str,
                           organism: str,
                           treatment: str,
                           retrieved: list[tuple[str, float]],
                           model: str,
                           ) -> dict:
    prompt = SEMANTIC_PROMPT.format(
        cluster_key=cluster_key,
        organism=organism,
        treatment=treatment,
        fields_description=EXTRACTION_FIELDS_DESCRIPTION,
        passages=passages,
    )
    client = _openai.OpenAI()
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
        )
        raw = response.choices[0].message.content.strip()
    except Exception as e:
        logger.error("Semantic extraction failed for cluster %s: %s",
                     cluster_key, e)
        return {}
    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        logger.error("Semantic JSON parse failed for cluster %s", cluster_key)
        return {}
    if not isinstance(parsed, dict):
        return {}
    parsed["source"] = "semantic"
    parsed["confidence"] = "medium"
    return parsed
