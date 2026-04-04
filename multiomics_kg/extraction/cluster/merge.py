# multiomics_kg/extraction/cluster/merge.py
"""Merge visual + semantic + table results into lists per field.

No winner-picking. Each field becomes a list of
{value, source, confidence} options. Stage 2 LLM sees all options.
"""

CONFIDENCE_MAP = {
    "table": "very_high",
    "visual": "high",
    "semantic": "medium",
}

CONFIDENCE_ORDER = {"very_high": 0, "high": 1, "medium": 2, "low": 3}

# Fields that get stacked as lists of options
EXTRACTION_FIELDS = [
    "enrichment_category", "enrichment_pvalue", "enrichment_significant",
    "enrichment_details", "genes", "direction", "cluster_description",
    "temporal_pattern", "peak_time", "period_description", "light_phase",
]

# Fields pooled into one combined list
POOL_FIELDS = ["supporting_quotes", "all_enrichments"]

# Fields passed through from table (ground truth scalars)
PASSTHROUGH_FIELDS = ["gene_count", "retrieved_passages"]


def merge_paths(cluster_keys: list[str],
                visual: dict[str, dict],
                semantic: dict[str, dict],
                table: dict[str, dict],
                ) -> dict[str, dict]:
    """Merge extraction results from three paths.

    Each field becomes a list of {value, source, confidence} dicts,
    sorted by confidence (highest first).
    Quotes are pooled and sorted by relevance_score.
    """
    paths = [
        ("table", table),
        ("visual", visual),
        ("semantic", semantic),
    ]

    results = {}
    for key in cluster_keys:
        merged: dict = {}

        for field in EXTRACTION_FIELDS:
            options = []
            for source_name, path_data in paths:
                cluster_data = path_data.get(key, {})
                val = cluster_data.get(field)
                if val is None or val == "" or val == []:
                    continue
                if isinstance(val, str) and val.lower() == "not described in paper":
                    continue
                conf = cluster_data.get(f"confidence_{field}",
                                        cluster_data.get("confidence",
                                                         CONFIDENCE_MAP.get(source_name, "medium")))
                options.append({
                    "value": val,
                    "source": source_name,
                    "confidence": conf,
                })
            if options:
                options.sort(key=lambda x: CONFIDENCE_ORDER.get(x["confidence"], 9))
                merged[field] = options

        for field in POOL_FIELDS:
            pooled = []
            for source_name, path_data in paths:
                cluster_data = path_data.get(key, {})
                items = cluster_data.get(field, [])
                if isinstance(items, list):
                    for item in items:
                        if isinstance(item, dict):
                            item.setdefault("source", source_name)
                        pooled.append(item)
            if pooled:
                if field == "supporting_quotes":
                    pooled.sort(key=lambda x: x.get("relevance_score", 0),
                                reverse=True)
                merged[field] = pooled

        for field in PASSTHROUGH_FIELDS:
            for _source_name, path_data in paths:
                src_data = path_data.get(key, {})
                if field in src_data:
                    merged[field] = src_data[field]
                    break

        results[key] = merged

    return results
