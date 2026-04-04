# multiomics_kg/extraction/cluster/pipeline.py
"""Pipeline orchestrator: gather → synthesize → validate.

CLI entry point with --stage, --from-stage, --path, --table options.
"""
import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml

logger = logging.getLogger(__name__)


def _find_project_root(start: Path) -> Path:
    p = start
    for _ in range(10):
        if (p / "create_knowledge_graph.py").exists():
            return p
        p = p.parent
    return start


def _find_main_pdf(paper_dir: Path, config: dict) -> Optional[Path]:
    pdf_rel = config.get("publication", {}).get("papermainpdf", "")
    if not pdf_rel:
        return None
    project_root = _find_project_root(paper_dir)
    pdf_path = project_root / pdf_rel
    return pdf_path if pdf_path.exists() else None


def _lookup_doi(paper_name: str, paper_dir: Path) -> Optional[str]:
    project_root = _find_project_root(paper_dir)
    cache_path = project_root / "cache" / "pdf_extraction_cache.json"
    if not cache_path.exists():
        return None
    try:
        with open(cache_path) as f:
            cache = json.load(f)
        for entry in cache.values():
            if isinstance(entry, dict) and entry.get("papername") == paper_name:
                return entry.get("doi")
    except Exception:
        pass
    return None


def run_stage1(paperconfig_path: Path, table_key: str, table_config: dict,
               path_filter: Optional[str] = None) -> dict[str, dict]:
    """Stage 1: Run extraction paths and merge."""
    from multiomics_kg.extraction.cluster.table import run_table
    from multiomics_kg.extraction.cluster.visual import run_visual
    from multiomics_kg.extraction.cluster.semantic import run_semantic
    from multiomics_kg.extraction.cluster.merge import merge_paths

    import pandas as pd

    paper_dir = paperconfig_path.parent
    project_root = _find_project_root(paper_dir)

    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    organism = table_config["organism"]
    treatment = table_config.get("treatment", "")
    cluster_method = table_config.get("cluster_method", "unknown")
    main_pdf = _find_main_pdf(paper_dir, config)

    # Read cluster keys from CSV
    csv_path = Path(table_config["filename"])
    if not csv_path.is_absolute():
        csv_path = project_root / csv_path
    df = pd.read_csv(csv_path)
    cluster_col = table_config["cluster_col"]
    cluster_keys = [str(k) for k in sorted(df[cluster_col].unique(),
                                            key=lambda x: (not str(x).isdigit(), str(x)))]

    visual_results: dict[str, dict] = {}
    semantic_results: dict[str, dict] = {}
    table_results: dict[str, dict] = {}

    # Run table (always, unless filtered out)
    if path_filter is None or path_filter == "table":
        logger.info("Running Path: table...")
        table_results = run_table(table_config, paper_dir, cluster_keys,
                                  organism_hint=organism, project_root=project_root)
        logger.info("  table: extracted data for %d clusters", len(table_results))

    # Run visual (needs main PDF)
    if (path_filter is None or path_filter == "visual") and main_pdf:
        logger.info("Running Path: visual...")
        visual_results = run_visual(main_pdf, paper_dir, cluster_keys,
                                    cluster_method, organism, treatment)
        logger.info("  visual: extracted data for %d clusters", len(visual_results))

    # Run semantic (seeded by visual + table, needs main PDF)
    if path_filter is None or path_filter == "semantic":
        if main_pdf:
            # Pre-merge visual+table to seed semantic queries
            seed = merge_paths(cluster_keys, visual_results, {}, table_results)
            logger.info("Running Path: semantic...")
            semantic_results = run_semantic(main_pdf, paper_dir, cluster_keys,
                                           organism, treatment, seed_data=seed)
            logger.info("  semantic: extracted data for %d clusters", len(semantic_results))

    # Final merge
    merged = merge_paths(cluster_keys, visual_results, semantic_results, table_results)
    logger.info("Merged: %d clusters with data", len(merged))

    return merged


def run_pipeline(paperconfig_path: Path,
                 table_key: Optional[str] = None,
                 stage: Optional[int] = None,
                 from_stage: Optional[int] = None,
                 path_filter: Optional[str] = None,
                 force_llm: bool = False,
                 force_all: bool = False,
                 ) -> list[Path]:
    """Run the full extraction pipeline.

    Output is written to versioned run directories under
    ``paper_dir/.extraction_cache/{entry_key}/runs/<timestamp>/``.

    A backward-compatible monolithic JSON is also written to
    ``paper_dir/cluster_extraction_{tkey}.json`` so that downstream
    consumers (e.g. cluster_adapter) keep working.

    Parameters
    ----------
    force_llm : bool
        When True, re-run LLM stages even if cached results exist.
        (Passed through for future stage-level caching.)
    force_all : bool
        When True, re-run everything from scratch regardless of cache.
        (Passed through for future stage-level caching.)
    """
    from multiomics_kg.extraction.cluster.synthesis import run_synthesis
    from multiomics_kg.extraction.cluster.validation import run_validation
    from multiomics_kg.extraction.cluster.run_manager import RunManager

    import pandas as pd

    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    pub = config["publication"]
    paper_name = pub["papername"]
    paper_dir = paperconfig_path.parent
    project_root = _find_project_root(paper_dir)
    supp = pub.get("supplementary_materials", {})

    # Find gene_clusters tables
    cluster_tables = {k: v for k, v in supp.items()
                      if isinstance(v, dict) and v.get("type") == "gene_clusters"}
    if table_key:
        cluster_tables = {table_key: cluster_tables[table_key]}

    if not cluster_tables:
        logger.warning("No gene_clusters entries found")
        return []

    # Determine which stages to run
    run_stages = {1, 2, 3}
    if stage is not None:
        run_stages = {stage}
    elif from_stage is not None:
        run_stages = {s for s in run_stages if s >= from_stage}

    output_paths = []
    for tkey, tconfig in cluster_tables.items():
        organism = tconfig["organism"]
        organism_short = organism.split()[-1].lower()

        # ── Set up RunManager ──
        cache_dir = paper_dir / ".extraction_cache"
        rm = RunManager(cache_dir, tkey)
        prev_run = rm.get_current_run()
        run_dir = rm.create_run()
        logger.info("Created run directory: %s", run_dir)

        # Load previous run data if skipping stages
        prev_stage1: dict = {}
        prev_stage2: dict = {}
        prev_stage3: dict = {}
        if prev_run is not None and 1 not in run_stages:
            prev_stage1 = rm.read_stage(prev_run, 1)
            prev_stage2 = rm.read_stage(prev_run, 2)
            prev_stage3 = rm.read_stage(prev_run, 3)

        # CSV path for validation
        csv_path = Path(tconfig["filename"])
        if not csv_path.is_absolute():
            csv_path = project_root / csv_path

        # Cluster keys
        df = pd.read_csv(csv_path)
        cluster_col = tconfig["cluster_col"]
        cluster_keys = [str(k) for k in sorted(df[cluster_col].unique(),
                                                key=lambda x: (not str(x).isdigit(), str(x)))]

        treatment = tconfig.get("treatment", "")
        cluster_method = tconfig.get("cluster_method", "unknown")
        main_pdf = _find_main_pdf(paper_dir, config)
        doi = _lookup_doi(paper_name, paper_dir)

        # ── Stage 1 ──
        if 1 in run_stages:
            logger.info("=== Stage 1: Gather + Merge (%s / %s) ===", paper_name, tkey)
            merged = run_stage1(paperconfig_path, tkey, tconfig, path_filter)
        else:
            merged = prev_stage1
            logger.info("Stage 1: using cached results (%d clusters)", len(merged))

        rm.write_stage(run_dir, 1, merged)

        # ── Stage 2 ──
        if 2 in run_stages:
            logger.info("=== Stage 2: Synthesis ===")
            analysis_name = tconfig.get("name", tkey)
            stage2 = run_synthesis(merged, analysis_name, organism, treatment,
                                   cluster_method)
            logger.info("Stage 2: synthesized %d clusters", len(stage2))
        else:
            stage2 = prev_stage2

        rm.write_stage(run_dir, 2, stage2)

        # ── Stage 3 ──
        if 3 in run_stages and main_pdf:
            logger.info("=== Stage 3: Validation ===")
            analysis_name = tconfig.get("name", tkey)
            stage3 = run_validation(main_pdf, paper_dir, csv_path, stage2,
                                    cluster_keys,
                                    analysis_name=analysis_name)
            verdicts = {}
            for k, v in stage3.items():
                verdicts[v.get("verdict", "?")] = verdicts.get(v.get("verdict", "?"), 0) + 1
            logger.info("Stage 3: %s", verdicts)
        else:
            stage3 = prev_stage3

        rm.write_stage(run_dir, 3, stage3)

        # ── Copy forward reviews from previous run ──
        if prev_run is not None:
            rm.copy_forward_reviews(prev_run, run_dir, merged)

        # ── Write metadata ──
        metadata = {
            "paper": paper_name,
            "doi": doi,
            "organism": organism,
            "extracted_at": datetime.now().isoformat(timespec="seconds"),
            "table_key": tkey,
            "stages_run": sorted(run_stages),
            "force_llm": force_llm,
            "force_all": force_all,
        }
        metadata_path = run_dir / "metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

        # ── Write markdown report to run dir ──
        output_for_report = {
            "metadata": metadata,
            "stage2_results": stage2,
            "stage3_validation": stage3,
        }
        report = _generate_report(output_for_report, cluster_keys)
        report_path = run_dir / "report.md"
        with open(report_path, "w") as f:
            f.write(report)
        logger.info("Wrote %s", report_path)

        # ── Write backward-compatible monolithic JSON ──
        out_path = paper_dir / f"cluster_extraction_{tkey}.json"
        output = {
            "metadata": metadata,
            "stage1_merged": merged,
            "stage2_results": stage2,
            "stage3_validation": stage3,
        }

        with open(out_path, "w") as f:
            json.dump(output, f, indent=2, ensure_ascii=False)
        logger.info("Wrote %s", out_path)

        # Write backward-compatible review report
        compat_report_path = out_path.with_suffix(".md")
        with open(compat_report_path, "w") as f:
            f.write(report)

        # ── Finalize: update current symlink now that all data is written ──
        rm.finalize_run(run_dir)

        output_paths.append(out_path)

    return output_paths


def _generate_report(output: dict, cluster_keys: list[str]) -> str:
    """Generate markdown review report."""
    meta = output["metadata"]
    s2 = output.get("stage2_results", {})
    s3 = output.get("stage3_validation", {})

    lines = [
        f"# Cluster Extraction Review: {meta['paper']}",
        f"**Organism:** {meta['organism']}  ",
        f"**Extracted:** {meta['extracted_at']}",
        "",
    ]

    for key in cluster_keys:
        desc = s2.get(key, {})
        val = s3.get(key, {})
        verdict = val.get("verdict", "no validation")
        cluster_id = desc.get("id", "?")

        if verdict == "pass":
            marker = "PASS"
        elif verdict == "warn":
            marker = "WARN"
        elif verdict == "fail":
            marker = "FAIL"
        else:
            marker = "NOT VALIDATED"

        lines.append(f"## Cluster {key}: {cluster_id} — {marker}")
        lines.append("")

        func = desc.get("functional_description", "")
        behav = desc.get("behavioral_description", "")
        if func:
            lines.append(f"- **functional:** {func}")
        if behav:
            lines.append(f"- **behavioral:** {behav}")

        if verdict in ("warn", "fail"):
            explanation = val.get("explanation", "")
            if explanation:
                lines.append(f"- **issue:** {explanation}")
            missing = val.get("missing_info", "no")
            if missing != "no":
                lines.append(f"- **missing:** {missing}")
            hallucination = val.get("hallucination", "no")
            if hallucination != "no":
                lines.append(f"- **hallucination:** {hallucination}")

        lines.append("")

    return "\n".join(lines)


def main():
    try:
        from dotenv import load_dotenv
        load_dotenv()
    except ImportError:
        pass

    parser = argparse.ArgumentParser(
        description="Extract cluster descriptions from paper PDF + supp files"
    )
    parser.add_argument("paperconfig", type=Path, help="Path to paperconfig.yaml")
    parser.add_argument("--table", default=None, help="Process only this table key")
    parser.add_argument("--stage", type=int, choices=[1, 2, 3], default=None,
                        help="Run only this stage")
    parser.add_argument("--from-stage", type=int, choices=[1, 2, 3], default=None,
                        help="Run from this stage onwards")
    parser.add_argument("--path", choices=["visual", "semantic", "table"], default=None,
                        help="Run only this path within stage 1")
    parser.add_argument("--force-llm", action="store_true",
                        help="Re-run LLM stages even if cached")
    parser.add_argument("--force-all", action="store_true",
                        help="Re-run everything from scratch")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.paperconfig.exists():
        parser.error(f"Not found: {args.paperconfig}")

    paths = run_pipeline(
        args.paperconfig,
        table_key=args.table,
        stage=args.stage,
        from_stage=args.from_stage,
        path_filter=args.path,
        force_llm=args.force_llm,
        force_all=args.force_all,
    )
    for p in paths:
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
