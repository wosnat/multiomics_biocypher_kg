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

    for field, options in merged_data.items():
        if field in ("supporting_quotes", "gene_count"):
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


def run_synthesis(merged_results: dict[str, dict],
                  paper_name: str,
                  organism: str,
                  treatment: str,
                  cluster_method: str,
                  model: str = "gpt-5-nano",
                  ) -> dict[str, dict]:
    """Stage 2: synthesize descriptions + id from merged data."""
    if _openai is None:
        logger.error("openai package not installed")
        return {}

    blocks = []
    for key in sorted(merged_results, key=lambda k: (not k.isdigit(), k)):
        blocks.append(_format_cluster_block(key, merged_results[key]))
    cluster_blocks = "\n\n".join(blocks)

    prompt = SYNTHESIS_PROMPT.format(
        paper_name=paper_name,
        organism=organism,
        treatment=treatment,
        cluster_method=cluster_method,
        cluster_blocks=cluster_blocks,
    )

    client = _openai.OpenAI()
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
        )
        raw = response.choices[0].message.content.strip()
    except Exception as e:
        logger.error("Stage 2 synthesis failed: %s", e)
        return {}

    if raw.startswith("```"):
        raw = re.sub(r'^```\w*\n?', '', raw)
        raw = re.sub(r'\n?```$', '', raw)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error("Stage 2 JSON parse failed: %s", e)
        return {}

    valid_keys = set(merged_results.keys())
    results: dict[str, dict] = {}
    seen_ids: set[str] = set()
    for k, v in parsed.items():
        norm = _normalize_key(str(k), valid_keys)
        if norm and isinstance(v, dict):
            cluster_id = v.get("id", "")
            if cluster_id in seen_ids:
                cluster_id = f"{cluster_id}_{norm}"
                v["id"] = cluster_id
            seen_ids.add(cluster_id)
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
    return None
