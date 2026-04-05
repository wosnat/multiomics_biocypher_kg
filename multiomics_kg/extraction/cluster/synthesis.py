# multiomics_kg/extraction/cluster/synthesis.py
"""Stage 2: Synthesize descriptions + id from merged Stage 1 data."""
import json
import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)

try:
    import openai as _openai
except Exception:
    _openai = None

from multiomics_kg.extraction.cluster.prompts import SYNTHESIS_PROMPT


def _format_cluster_block(key: str, merged_data: dict) -> str:
    """Format one cluster's merged data for the synthesis prompt."""
    lines = [f"### Cluster {key}"]
    gene_count = merged_data.get("gene_count", "?")
    lines[0] += f" ({gene_count} genes)"

    # Fields to skip: genes (locus tags not useful for synthesis),
    # supporting_quotes (shown separately), gene_count (metadata only),
    # retrieved_passages (raw RAG data), all_enrichments (redundant)
    skip_fields = {"supporting_quotes", "gene_count", "genes",
                   "retrieved_passages", "all_enrichments"}

    for field, options in merged_data.items():
        if field in skip_fields:
            continue
        if not isinstance(options, list):
            continue
        lines.append(f"**{field}:**")
        for opt in options:
            if isinstance(opt, dict) and "value" in opt:
                lines.append(f"  - [{opt.get('source', '?')}, "
                             f"{opt.get('confidence', '?')}]: {opt['value']}")

    quotes = merged_data.get("supporting_quotes", [])
    if quotes:
        lines.append(f"**supporting_quotes:** ({len(quotes)} total)")
        for q in quotes[:10]:
            score = q.get("relevance_score", "?")
            quote_text = q.get("quote", str(q))[:200]
            lines.append(f"  - (score={score}): \"{quote_text}\"")

    return "\n".join(lines)


def _synthesize_single_cluster(
    cluster_key: str,
    merged_data: dict,
    paper_name: str,
    organism: str,
    cluster_method: str,
    model: str,
) -> dict:
    """Synthesize descriptions for one cluster."""
    cluster_block = _format_cluster_block(cluster_key, merged_data)

    prompt = SYNTHESIS_PROMPT.format(
        paper_name=paper_name,
        organism=organism,
        cluster_method=cluster_method,
        cluster_blocks=cluster_block,
    )

    client = _openai.OpenAI()
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
        )
        raw = response.choices[0].message.content.strip()
    except Exception as e:
        logger.error("Synthesis failed for cluster %s: %s", cluster_key, e)
        return {}

    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    if not raw.strip():
        logger.error("Synthesis returned empty response for cluster %s", cluster_key)
        return {}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Synthesis JSON parse failed for cluster %s: %s\nRaw: %s",
                     cluster_key, e, raw[:500])
        return {}

    if not isinstance(parsed, dict):
        return {}

    # Extract the cluster entry — LLM may wrap in a dict with the key
    if "id" in parsed:
        return parsed  # direct flat response
    if cluster_key in parsed and isinstance(parsed[cluster_key], dict):
        return parsed[cluster_key]
    entries = [v for v in parsed.values() if isinstance(v, dict)]
    if len(entries) == 1:
        return entries[0]
    return parsed


def run_synthesis(merged_results: dict[str, dict],
                  paper_name: str,
                  organism: str,
                  treatment: str,
                  cluster_method: str,
                  model: str = "gpt-5-nano",
                  ) -> dict[str, dict]:
    """Stage 2: synthesize descriptions + id from merged data.

    Runs one LLM call per cluster to prevent cross-contamination.
    """
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    results: dict[str, dict] = {}
    seen_ids: set[str] = set()

    for key in sorted(merged_results, key=lambda k: (not k.isdigit(), k)):
        logger.info("Synthesizing cluster %s...", key)
        result = _synthesize_single_cluster(
            key, merged_results[key], paper_name, organism, cluster_method, model
        )
        if result:
            # Ensure unique IDs
            cluster_id = result.get("id", "")
            if cluster_id in seen_ids:
                cluster_id = f"{cluster_id}_{key}"
                result["id"] = cluster_id
            seen_ids.add(cluster_id)
            results[key] = result
        else:
            logger.warning("No synthesis result for cluster %s", key)

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
    return None
