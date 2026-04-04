"""Dry-run: generate all LLM prompts without calling APIs.

Loads data from the most recent extraction run and builds the exact prompts
that would be sent to each LLM stage. Writes them to a directory for review.

Usage:
    uv run python scripts/dry_run_prompts.py "data/.../tolonen 2006/paperconfig.yaml"
    uv run python scripts/dry_run_prompts.py "data/.../tolonen 2006/paperconfig.yaml" --entry med4_kmeans_nstarvation
    uv run python scripts/dry_run_prompts.py "data/.../tolonen 2006/paperconfig.yaml" --cluster 1
"""

import argparse
import json
import sys
from pathlib import Path

from dotenv import load_dotenv
load_dotenv()

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager
from multiomics_kg.extraction.cluster.prompts import (
    VISUAL_PROMPT,
    SEMANTIC_PROMPT,
    SYNTHESIS_PROMPT,
    JUDGE_PROMPT,
    EXTRACTION_FIELDS_DESCRIPTION,
)
from multiomics_kg.extraction.cluster.semantic import build_cluster_query
from multiomics_kg.extraction.cluster.merge import merge_paths


def load_paperconfig(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def get_cluster_tables(config: dict) -> dict:
    supp = config.get("publication", {}).get("supplementary_materials", {})
    return {k: v for k, v in supp.items() if v.get("type") == "gene_clusters"}


def dump_visual_prompt(entry_key: str, tconfig: dict, cluster_keys: list[str],
                       cluster_key: str = None):
    """Build the visual prompt for one cluster (or all if cluster_key is None)."""
    organism = tconfig["organism"]
    treatment = tconfig.get("treatment", "")
    cluster_method = tconfig.get("cluster_method", "unknown")
    analysis_name = tconfig.get("name", entry_key)

    keys = [cluster_key] if cluster_key else cluster_keys

    prompts = {}
    for ck in keys:
        prompt = VISUAL_PROMPT.format(
            n_clusters=len(cluster_keys),
            cluster_method=cluster_method,
            organism=organism,
            treatment=treatment,
            cluster_keys=ck,  # per-cluster: just this one
            fields_description=EXTRACTION_FIELDS_DESCRIPTION,
        )
        prompts[ck] = {
            "stage": "visual",
            "model": "gpt-4o",
            "analysis": analysis_name,
            "cluster": ck,
            "content": [
                {"type": "note", "text": f"[{15} PDF pages would be prepended here]"},
                {"type": "text", "text": prompt},
            ],
        }
    return prompts


def dump_semantic_prompt(entry_key: str, tconfig: dict, stage1_merged: dict,
                         cluster_keys: list[str], cluster_key: str = None):
    """Build the semantic prompt for one cluster."""
    organism = tconfig["organism"]
    treatment = tconfig.get("treatment", "")
    analysis_name = tconfig.get("name", entry_key)

    keys = [cluster_key] if cluster_key else cluster_keys

    prompts = {}
    for ck in keys:
        seed = stage1_merged.get(ck, {})
        query = build_cluster_query(ck, seed, organism=organism)

        prompt = SEMANTIC_PROMPT.format(
            cluster_key=ck,
            organism=organism,
            treatment=treatment,
            fields_description=EXTRACTION_FIELDS_DESCRIPTION,
            passages=f"[RAG query: \"{query}\" → top 8 chunks would be here]",
        )
        prompts[ck] = {
            "stage": "semantic",
            "model": "gpt-5-nano",
            "analysis": analysis_name,
            "cluster": ck,
            "rag_query": query,
            "content": prompt,
        }
    return prompts


def dump_synthesis_prompt(entry_key: str, tconfig: dict, stage1_merged: dict,
                          cluster_keys: list[str], cluster_key: str = None):
    """Build the synthesis prompt — currently batch, showing what per-cluster would look like."""
    organism = tconfig["organism"]
    treatment = tconfig.get("treatment", "")
    cluster_method = tconfig.get("cluster_method", "unknown")
    analysis_name = tconfig.get("name", entry_key)
    paper_name = analysis_name  # close enough for dry run

    keys = [cluster_key] if cluster_key else cluster_keys

    # For per-cluster: each cluster gets its own synthesis call
    from multiomics_kg.extraction.cluster.synthesis import _format_cluster_block

    prompts = {}
    for ck in keys:
        cluster_block = _format_cluster_block(ck, stage1_merged.get(ck, {}))

        prompt = SYNTHESIS_PROMPT.format(
            paper_name=analysis_name,
            organism=organism,
            cluster_method=cluster_method,
            cluster_blocks=cluster_block,
        )
        prompts[ck] = {
            "stage": "synthesis",
            "model": "gpt-5-nano",
            "analysis": analysis_name,
            "cluster": ck,
            "content": prompt,
        }
    return prompts


def dump_validation_prompt(entry_key: str, tconfig: dict, stage2_results: dict,
                           cluster_keys: list[str], cluster_key: str = None):
    """Build the validation prompt for one cluster."""
    analysis_name = tconfig.get("name", entry_key)

    keys = [cluster_key] if cluster_key else cluster_keys

    prompts = {}
    for ck in keys:
        s2 = stage2_results.get(ck, {})
        description = (
            f"Cluster {ck}: id={s2.get('id', '?')}, "
            f"functional={s2.get('functional_description', '(sentinel)')}, "
            f"behavioral={s2.get('behavioral_description', '(sentinel)')}"
        )
        prompt = JUDGE_PROMPT.format(
            descriptions=description,
            analysis_name=analysis_name,
        )
        prompts[ck] = {
            "stage": "validation",
            "model": "gpt-4o",
            "analysis": analysis_name,
            "cluster": ck,
            "content": [
                {"type": "note", "text": f"[{15} PDF pages would be prepended here]"},
                {"type": "note", "text": "[CSV summary would be here]"},
                {"type": "text", "text": prompt},
            ],
        }
    return prompts


def main():
    parser = argparse.ArgumentParser(description="Dry-run prompt generation")
    parser.add_argument("paperconfig", type=Path)
    parser.add_argument("--entry", help="Specific entry key")
    parser.add_argument("--cluster", help="Specific cluster key")
    parser.add_argument("--stage", choices=["visual", "semantic", "synthesis", "validation"],
                        help="Only show prompts for this stage")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Write prompts to files instead of stdout")
    args = parser.parse_args()

    config = load_paperconfig(args.paperconfig)
    paper_dir = args.paperconfig.parent
    tables = get_cluster_tables(config)

    if args.entry:
        tables = {args.entry: tables[args.entry]}

    for entry_key, tconfig in tables.items():
        print(f"\n{'='*80}")
        print(f"ENTRY: {entry_key}")
        print(f"ANALYSIS: {tconfig.get('name', entry_key)}")
        print(f"ORGANISM: {tconfig['organism']}")
        print(f"{'='*80}")

        # Load latest run data
        cache_dir = paper_dir / ".extraction_cache"
        rm = RunManager(cache_dir, entry_key)
        run_dir = rm.get_current_run()

        if run_dir is None:
            # Try latest run directory directly
            runs = rm.list_runs()
            if runs:
                run_dir = runs[-1]

        stage1 = rm.read_stage(run_dir, 1) if run_dir else {}
        stage2 = rm.read_stage(run_dir, 2) if run_dir else {}
        cluster_keys = sorted(stage1.keys(), key=lambda k: (not k.isdigit(), k))

        if not cluster_keys and stage2:
            cluster_keys = sorted(stage2.keys(), key=lambda k: (not k.isdigit(), k))

        if not cluster_keys:
            print("  No cluster data found in any run.")
            continue

        print(f"CLUSTERS: {cluster_keys}")

        stages_to_show = [args.stage] if args.stage else ["visual", "semantic", "synthesis", "validation"]

        all_prompts = {}
        for stage in stages_to_show:
            if stage == "visual":
                prompts = dump_visual_prompt(entry_key, tconfig, cluster_keys, args.cluster)
            elif stage == "semantic":
                prompts = dump_semantic_prompt(entry_key, tconfig, stage1, cluster_keys, args.cluster)
            elif stage == "synthesis":
                prompts = dump_synthesis_prompt(entry_key, tconfig, stage1, cluster_keys, args.cluster)
            elif stage == "validation":
                prompts = dump_validation_prompt(entry_key, tconfig, stage2, cluster_keys, args.cluster)
            else:
                continue

            all_prompts[stage] = prompts

            for ck, prompt_data in prompts.items():
                print(f"\n{'─'*60}")
                print(f"STAGE: {stage} | CLUSTER: {ck} | MODEL: {prompt_data.get('model','?')}")
                if "rag_query" in prompt_data:
                    print(f"RAG QUERY: {prompt_data['rag_query']}")
                print(f"{'─'*60}")

                content = prompt_data["content"]
                if isinstance(content, list):
                    for part in content:
                        if part.get("type") == "note":
                            print(f"\n  >>> {part['text']}")
                        elif part.get("type") == "text":
                            print(f"\n{part['text']}")
                else:
                    print(f"\n{content}")

        # Optionally write to files
        if args.output_dir:
            out_dir = args.output_dir / entry_key
            out_dir.mkdir(parents=True, exist_ok=True)
            for stage, prompts in all_prompts.items():
                for ck, prompt_data in prompts.items():
                    path = out_dir / f"{stage}_cluster_{ck}.json"
                    with open(path, "w") as f:
                        json.dump(prompt_data, f, indent=2, ensure_ascii=False)
            print(f"\nPrompts written to {args.output_dir / entry_key}/")


if __name__ == "__main__":
    main()
