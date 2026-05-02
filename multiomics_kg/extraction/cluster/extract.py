"""Cluster description extraction via OpenAI Responses API.

Usage:
    uv run python -m multiomics_kg.extraction.cluster.extract             # extract all
    uv run python -m multiomics_kg.extraction.cluster.extract --report    # generate report
    uv run python -m multiomics_kg.extraction.cluster.extract --verify    # report + checks
    uv run python -m multiomics_kg.extraction.cluster.extract --dry-run   # show plan
    uv run python -m multiomics_kg.extraction.cluster.extract --paper "Tolonen 2006"
    uv run python -m multiomics_kg.extraction.cluster.extract --entry mit9313_kmeans_nstarvation --force
"""
import argparse
import json
import logging
import os
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel

from multiomics_kg.extraction.cluster.extraction_utils import (
    find_all_entries,
    get_cluster_data,
    list_extraction_files,
    load_cluster_summaries,
    load_extraction,
    match_cluster_keys,
    save_extraction,
)

logger = logging.getLogger(__name__)


# ── Pydantic schemas ──


class SupportingQuote(BaseModel):
    quote: str
    location: str


class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    temporal_pattern: str
    expression_dynamics: str
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    source_figures: list[str]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str


class AnalysisExtraction(BaseModel):
    clusters: list[ClusterExtraction]


# ── Prompt layer ──


SHARED_RULES = """\
## Rules

- For clusters the paper does NOT discuss: set description fields to "N/A". \
Do NOT say what the paper didn't do (e.g., "Specific genes are not individually named" or \
"No enrichment was reported"). Just use "N/A" and move on. If you have nothing to say \
about a field, use "N/A" — never fill space with negative statements.
- functional_description is ONLY about gene identity and pathway membership — what genes \
or pathways are enriched in this cluster. 2-3 sentences max. \
Do NOT put expression dynamics, temporal patterns, condition descriptions, or \
periodicity/cycling information here — those belong in temporal_pattern. \
BAD: "Enriched for transcripts exhibiting 24-h periodicity in coculture under diel conditions" \
BAD: "Genes show maximum mRNA abundances at intermediate irradiance" \
GOOD: "Enriched for photosystem I and II components. Includes genes psbA, psaA, and psaB."
- Do NOT include locus tags (PMM*, PMT*, P9301_*, NATL2_*, PMN2A_RS*, tll*, SY28_*, \
BSR22_*, ALT831_RS*, MIT1002_*, SYNW*, sync_*, A9601_*, WP_*, cds-*). Use gene names only.
- Do NOT include treatment conditions in descriptions — those live on the analysis node.
- Do NOT describe cluster membership statistics (gene counts, sample gene IDs) — \
those are already stored on the cluster node.
- ONLY state what the paper EXPLICITLY says about THIS SPECIFIC cluster. \
Do NOT synthesize from other clusters, other sections, or your own knowledge.
- Each field is independent: a cluster can have a functional description but \
temporal_pattern = "N/A", or vice versa. Partial descriptions are fine.
- Max 3-5 named genes per cluster, only those the paper highlights for THIS cluster.
- For enrichment: only report enrichment the PAPER performed with p-values. \
Do NOT compute your own. Max 3 decimal places or scientific notation.
- supporting_quotes: direct quotes from the paper supporting your description. REQUIRED.
- source_figures: figure/table references used as evidence (e.g., "Figure 4A"). REQUIRED.
- id format: {{organism_short}}_{{dynamics}}_{{theme}} in snake_case. \
For undiscussed clusters: {{organism_short}}_cluster_{{KEY}}.
- name format: "{{Organism}} cluster {{KEY}} ({{theme}})" — under 60 chars. \
Use the EXACT cluster key from the list below. \
For undiscussed clusters (functional and temporal both N/A): just "{{Organism}} cluster {{KEY}}" \
with no theme — do NOT invent a theme from gene annotations.
- self_assessment: your confidence in the extraction (high/medium/low). \
assessment_notes: what you're uncertain about.
"""

TYPE_RULES = {
    "time_course": """\
## Type-Specific Guidance: Time Course

- expression_dynamics: Response timing label. Examples: "early transient", \
"late sustained", "biphasic", "gradual increase".
- temporal_pattern: When the change begins, how fast it develops, whether it \
persists or reverses. Include timing numbers from the paper. 1-2 sentences. \
Do NOT add biological interpretation.

### Example — Time-course cluster with enrichment:
{{"id": "mit9313_early_transport", "name": "MIT9313 cluster 1 (transport)", \
"functional_description": "Enriched for transport and binding (p=0.04). Contains \
nitrogen transport genes urtA and the nitrite permease, and hli genes hliS and hli7.", \
"temporal_pattern": "Most rapidly upregulated cluster, with genes responding within \
the first hours of nitrogen starvation.", \
"expression_dynamics": "early transient", \
"enrichment_category": "transport and binding", \
"enrichment_pvalue": 0.04, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", "confidence_notes": "", \
"supporting_quotes": [{{"quote": "Cluster 1, the most rapidly and highly upregulated \
genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA", \
"location": "Page 3"}}], \
"source_figures": ["Figure 2"]}}

### Example — Undiscussed cluster:
{{"id": "mit9313_cluster_8", "name": "MIT9313 cluster 8 (not discussed)", \
"functional_description": "N/A", \
"temporal_pattern": "N/A", \
"expression_dynamics": "N/A", \
"enrichment_category": "", \
"enrichment_pvalue": null, "enrichment_significant": false, \
"self_assessment": "low", "assessment_notes": "Paper does not discuss this cluster.", \
"confidence_notes": "", "supporting_quotes": [], "source_figures": []}}
""",
    "diel": """\
## Type-Specific Guidance: Diel Cycling

- expression_dynamics: Peak phase label. Examples: "peaks at dawn", "peaks at dusk", \
"peaks mid-day", "peaks mid-night".
- temporal_pattern: Peak timing in hours, periodicity, phase relative to light/dark cycle. \
1-2 sentences. Do NOT add biological interpretation.

### Example — Diel cluster with enrichment:
{{"id": "med4_dawn_photosynthesis", "name": "Prochlorococcus cluster 1 (photosynthesis)", \
"functional_description": "Enriched for photosystem I and II components (p=1.5e-9). \
Includes genes psbA, psbD, and psaA.", \
"temporal_pattern": "Genes peak in expression near dawn (8.3 h) with 24h periodicity.", \
"expression_dynamics": "peaks near dawn", \
"enrichment_category": "Photosystem I and II", \
"enrichment_pvalue": 1.5e-09, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", "confidence_notes": "", \
"supporting_quotes": [{{"quote": "Expression of approximately half of photosystem (PS) II \
genes, including reaction center genes psbA and psbD, peak in abundance at mid-day", \
"location": "Page 6"}}], \
"source_figures": ["Figure 4A", "Table S4"]}}
""",
    "condition_comparison": """\
## Type-Specific Guidance: Condition Comparison

- expression_dynamics: Condition label. Examples: "up with light", "down at cold", \
"high across all", "up in coculture".
- temporal_pattern: Which conditions drive changes, direction and magnitude of response. \
1-2 sentences. Do NOT add biological interpretation.

### Example — Condition comparison cluster:
{{"id": "syn_cold_stress_C", "name": "Synechococcus cluster C (cold stress)", \
"functional_description": "Contains ribosomal protein genes rpsA, rplB and cold-shock \
protein cspA.", \
"temporal_pattern": "Strongly upregulated at Tmin relative to Topt, with the largest \
fold-changes among all clusters.", \
"expression_dynamics": "up at cold", \
"enrichment_category": "Translation", \
"enrichment_pvalue": 0.003, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", "confidence_notes": "", \
"supporting_quotes": [{{"quote": "Cluster C genes were most strongly induced at the \
minimum growth temperature", "location": "Page 8"}}], \
"source_figures": ["Figure 5"]}}
""",
    "expression_bin": """\
## Type-Specific Guidance: Expression Bin / Periodicity

- expression_dynamics: Category label. Examples: "periodic in L:D only", "not periodic", \
"periodic across all conditions", "constitutive".
- temporal_pattern: Which conditions or categories genes fall into, any distinguishing \
pattern. 1-2 sentences.
- functional_description: ONLY pathway enrichment if the paper reports it. The expression-bin \
criterion (e.g., "genes that show 24-h periodicity") is NOT functional — it belongs in \
temporal_pattern. If no pathway enrichment is reported, use "N/A".

### Example — Periodicity expression-bin cluster:
{{"id": "med4_periodic_ld_only", "name": "MED4 cluster periodic_LD (L:D periodic)", \
"functional_description": "Enriched for photosystem genes and light-harvesting complexes.", \
"temporal_pattern": "Genes show 24h periodicity in L:D cycle but lose periodicity in \
continuous light and darkness.", \
"expression_dynamics": "periodic in L:D only", \
"enrichment_category": "Photosynthesis", \
"enrichment_pvalue": 0.001, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", "confidence_notes": "", \
"supporting_quotes": [{{"quote": "A subset of genes were periodic only under light:dark \
cycling conditions", "location": "Page 5"}}], \
"source_figures": ["Figure 3"]}}
""",
}

# Map legacy cluster_type values to TYPE_RULES keys
_TYPE_ALIAS = {
    "response_pattern": "time_course",
    "diel_cycling": "diel",
    "diel_expression_pattern": "diel",
    "periodicity_classification": "expression_bin",
    "classification": "expression_bin",
    "expression_level": "condition_comparison",
    "expression_classification": "condition_comparison",
    "expression_pattern": "condition_comparison",
}

SELF_VERIFICATION = """\
## Self-Verification

Before finalizing your output, verify:
1. Each supporting_quote is a DIRECT quote from the paper (not paraphrased).
2. Each quote actually discusses THIS specific cluster, not a different one.
3. source_figures reference figures/tables that contain data about THIS cluster.
4. No field contains biological interpretation beyond what the paper states.
5. Undiscussed clusters have "N/A" for functional_description, temporal_pattern, \
and expression_dynamics — not filler text.
"""


def build_context_block(table: dict) -> str:
    """Build context block from paperconfig entry."""
    parts = [
        f"Analysis: {table.get('name', '')}",
        f"Organism: {table.get('organism', '')}",
        f"Clustering: {table.get('cluster_method', '')}",
        f"Type: {table.get('cluster_type', '')}",
        f"Treatment: {table.get('treatment', '')}",
    ]
    if table.get("treatment_type"):
        tt = table["treatment_type"]
        if isinstance(tt, list):
            parts.append(f"Treatment categories: {', '.join(tt)}")
        else:
            parts.append(f"Treatment categories: {tt}")
    if table.get("background_factors"):
        bf = table["background_factors"]
        if isinstance(bf, list):
            parts.append(f"Background factors: {', '.join(bf)}")
        else:
            parts.append(f"Background factors: {bf}")
    if table.get("experimental_context"):
        parts.append(f"Context: {table['experimental_context']}")
    if table.get("omics_type"):
        parts.append(f"Omics: {table['omics_type']}")
    if table.get("figure_hint"):
        parts.append(f"Key figures: {table['figure_hint']}")
    if table.get("time_points"):
        parts.append(f"Time points (hours): {table['time_points']}")
    if table.get("extraction_notes"):
        parts.append(f"\nExtraction guidance: {table['extraction_notes']}")
    return "\n".join(parts)


def _get_type_rules(table: dict) -> str:
    """Get type-specific rules for this cluster type."""
    ct = table.get("cluster_type", "")
    key = _TYPE_ALIAS.get(ct, ct)
    return TYPE_RULES.get(key, f"## Type: {ct}\n\nNo specific guidance available.")


def build_prompt(table_config: dict, cluster_summaries: dict[str, dict]) -> str:
    """Assemble the full developer prompt for one analysis."""
    ctx = build_context_block(table_config)
    summaries = format_cluster_summaries(cluster_summaries)
    type_rules = _get_type_rules(table_config)
    n = len(cluster_summaries)

    return (
        f"You are extracting structured descriptions of gene expression clusters "
        f"from a scientific paper.\n\n"
        f"{ctx}\n\n"
        f"For each cluster, extract all fields in the output schema.\n\n"
        f"{SHARED_RULES}\n"
        f"{type_rules}\n"
        f"{SELF_VERIFICATION}\n"
        f"CRITICAL: You MUST extract EXACTLY {n} clusters, one for each cluster key "
        f"listed below. Use the EXACT cluster keys as they appear — do NOT renumber, "
        f"skip, or merge clusters. Every key must appear exactly once in your output.\n\n"
        f"{summaries}\n"
    )


def format_cluster_summaries(clusters: dict[str, dict]) -> str:
    """Format cluster summaries for the prompt."""
    lines = []
    for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
        info = clusters[key]
        sample = ", ".join(info["sample_genes"][:3])
        lines.append(f"Cluster {key}: {info['gene_count']} genes (sample: {sample})")
    return "\n".join(lines)


# ── LLM layer ──


def upload_pdf(client, pdf_path: Path) -> str:
    """Upload PDF via Files API, return file_id."""
    f = client.files.create(file=open(pdf_path, "rb"), purpose="user_data")
    return f.id


def extract_analysis(client, file_ids, table_config, cluster_summaries, model="gpt-4.1-mini", flex=False):
    """Run extraction for one analysis. Returns (parsed, usage_dict)."""
    dev_msg = build_prompt(table_config, cluster_summaries)
    file_inputs = [{"type": "input_file", "file_id": fid} for fid in file_ids]

    kwargs = dict(
        model=model, temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": file_inputs + [
                {"type": "input_text", "text": f"Extract descriptions for all {len(cluster_summaries)} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {"input_tokens": resp.usage.input_tokens, "output_tokens": resp.usage.output_tokens, "duration_sec": elapsed}
    return resp.output[0].content[0].parsed, usage


def extract_paper(client, file_ids, tables_and_summaries, model="gpt-4.1-mini", flex=False):
    """Extract all analyses for one paper in a single call.
    Used for multi-organism papers (e.g. Tolonen) where per-analysis confuses the model.
    """
    # Build a combined prompt from all analyses
    context_parts = []
    all_summaries_parts = []
    type_rules_parts = []
    total_clusters = 0
    for table_config, cluster_summaries in tables_and_summaries:
        context_parts.append(build_context_block(table_config))
        all_summaries_parts.append(format_cluster_summaries(cluster_summaries))
        total_clusters += len(cluster_summaries)
        type_rules_parts.append(_get_type_rules(table_config))

    dev_msg = (
        f"You are extracting structured descriptions of gene expression clusters "
        f"from a scientific paper.\n\n"
        f"{chr(10).join(context_parts)}\n\n"
        f"For each cluster, extract all fields in the output schema.\n\n"
        f"{SHARED_RULES}\n"
        f"{chr(10).join(type_rules_parts)}\n"
        f"{SELF_VERIFICATION}\n"
        f"CRITICAL: You MUST extract EXACTLY {total_clusters} clusters, one for each "
        f"cluster key listed below. Use the EXACT cluster keys as they appear — do NOT "
        f"renumber, skip, or merge clusters. Every key must appear exactly once in your output.\n\n"
        f"{chr(10).join(all_summaries_parts)}\n"
    )

    file_inputs = [{"type": "input_file", "file_id": fid} for fid in file_ids]

    kwargs = dict(
        model=model, temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": file_inputs + [
                {"type": "input_text", "text": f"Extract descriptions for all {total_clusters} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {"input_tokens": resp.usage.input_tokens, "output_tokens": resp.usage.output_tokens, "duration_sec": elapsed}
    return resp.output[0].content[0].parsed, usage


# ── Report layer ──


def generate_report(entries: list[tuple[Path, str, dict, dict]]) -> str:
    """Generate diff-friendly markdown report from existing extraction JSONs."""
    lines = ["# Cluster Extraction Report\n"]

    for paper_dir, entry_key, table_config, pub in sorted(entries, key=lambda e: (pub_name(e[3]), e[1])):
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue

        paper = pub_name(pub)
        lines.append(f"## {paper} / {entry_key}\n")

        for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
            c = clusters[key]
            dynamics = c.get("expression_dynamics", "")
            assessment = c.get("self_assessment", "")
            lines.append(f"### Cluster {key} | {dynamics} | {assessment}")
            lines.append(f"**Name:** {c.get('name', '')}")
            enrich = c.get("enrichment_category", "")
            pval = c.get("enrichment_pvalue", "N/A")
            sig = c.get("enrichment_significant", "")
            lines.append(f"**Enrichment:** {enrich} (p={pval}, sig={sig})")
            lines.append(f"**Functional:** {c.get('functional_description', '')}")
            lines.append(f"**Temporal pattern:** {c.get('temporal_pattern', '')}")
            notes = c.get("confidence_notes", "")
            if notes:
                lines.append(f"**Confidence notes:** {notes}")
            assess = c.get("assessment_notes", "")
            if assess:
                lines.append(f"**Assessment notes:** {assess}")
            figs = c.get("source_figures", [])
            if figs:
                lines.append(f"**Sources:** {', '.join(figs)}")
            quotes = c.get("supporting_quotes", [])
            if quotes:
                lines.append("**Quotes:**")
                for q in quotes:
                    lines.append(f"- [{q.get('location', '')}] {q.get('quote', '')}")
            lines.append("")

    return "\n".join(lines)


def verify_quality(entries: list[tuple[Path, str, dict, dict]]) -> list[str]:
    """Run programmatic quality checks. Returns list of warning strings."""
    warnings = []

    locus_pat = re.compile(
        r"\b("
        r"PMM\d{3,}|PMT\d{3,}|P9301_\d+|tll\d{3,}|SY28_\d+|BSR22_\d+"
        r"|NATL[12]_\d+|PMN2A_RS\d+|ALT831_RS\d+|MIT1002_\d+"
        r"|SYNW\d{4}|sync_\d{4}|A9601_\d+"
        r"|WP_\d{6,}"
        r"|cds-[A-Z]{2}_\d+"
        r")\b"
    )

    filler_phrases = [
        "likely", "possibly", "not described but",
        "specific functions are not detailed",
        "not explicitly described",
    ]

    for paper_dir, entry_key, table_config, pub in entries:
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue
        paper = pub_name(pub)
        prefix = f"[{paper} / {entry_key}"

        # Check 1: locus tags in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "temporal_pattern"):
                text = c.get(field, "")
                m = locus_pat.search(text)
                if m:
                    warnings.append(
                        f"{prefix} / cluster {key}] locus tag in {field}: {m.group()}"
                    )

        # Check 2: filler on low-confidence clusters
        for key, c in clusters.items():
            if c.get("self_assessment") == "low":
                for field in ("functional_description", "temporal_pattern"):
                    text = c.get(field, "")
                    if text and text != "N/A":
                        warnings.append(
                            f"{prefix} / cluster {key}] low confidence but {field} "
                            f"is not 'N/A': {text[:60]}..."
                        )

        # Check 3: filler phrases in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "temporal_pattern"):
                text = c.get(field, "").lower()
                for phrase in filler_phrases:
                    if phrase in text:
                        warnings.append(
                            f"{prefix} / cluster {key}] filler phrase '{phrase}' in {field}"
                        )
                        break

        # Check 4: near-identical descriptions within analysis
        descs = {}
        for key, c in clusters.items():
            fd = c.get("functional_description", "")
            if fd and fd != "N/A" and len(fd) > 50:
                prefix_50 = fd[:50]
                if prefix_50 in descs:
                    warnings.append(
                        f"{prefix} / cluster {key}] near-identical functional_description "
                        f"as cluster {descs[prefix_50]}"
                    )
                else:
                    descs[prefix_50] = key

    return warnings


def pub_name(pub: dict) -> str:
    """Extract paper name from publication config."""
    return pub.get("papername", "Unknown")


# ── CLI ──


def main():
    from dotenv import load_dotenv
    load_dotenv()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Cluster extraction pipeline")
    parser.add_argument("--paper", help="Filter by paper name (partial match)")
    parser.add_argument("--entry", help="Filter by entry key (exact)")
    parser.add_argument("--model", default=os.environ.get("CLUSTER_EXTRACTION_MODEL", "gpt-4.1-mini"))
    parser.add_argument("--flex", action="store_true", help="Flex processing (50%% cheaper)")
    parser.add_argument("--force", action="store_true", help="Overwrite existing")
    parser.add_argument("--dry-run", action="store_true", help="Show what would run")
    parser.add_argument("--report", action="store_true", help="Generate report (no API calls)")
    parser.add_argument("--verify", action="store_true", help="Run quality checks (implies --report)")
    args = parser.parse_args()

    if args.verify:
        args.report = True

    entries = find_all_entries()

    if args.paper:
        entries = [(d, k, t, p, e) for d, k, t, p, e in entries
                   if args.paper.lower() in p.get("papername", "").lower()]
    if args.entry:
        entries = [(d, k, t, p, e) for d, k, t, p, e in entries if k == args.entry]

    if args.report:
        # Report/verify take 4-tuples (no extraction config needed)
        entries_4 = [(d, k, t, p) for d, k, t, p, _e in entries]
        report = generate_report(entries_4)
        if args.verify:
            warnings = verify_quality(entries_4)
            if warnings:
                report += "\n## Warnings\n\n"
                report += "\n".join(f"- {w}" for w in warnings) + "\n"
            else:
                report += "\n## Warnings\n\nNo issues found.\n"
        report_path = Path("data/cluster_extraction_report.md")
        report_path.write_text(report)
        logger.info(f"Report written to {report_path}")
        if args.verify:
            logger.info(f"{len(warnings)} warnings")
        return

    if not entries:
        logger.warning("No matching entries found")
        return

    if args.dry_run:
        print(f"\nWould process {len(entries)} entries:\n")
        for paper_dir, entry_key, table, pub, extraction in entries:
            summaries = load_cluster_summaries(table)
            exists = entry_key in list_extraction_files(paper_dir)
            status = "EXISTS" if exists else "NEW"
            scope = extraction.get("scope", "per_analysis")
            extra_pdfs = extraction.get("additional_pdfs", [])
            extra_info = ""
            if scope != "per_analysis":
                extra_info += f" scope={scope}"
            if extra_pdfs:
                extra_info += f" +{len(extra_pdfs)} PDFs"
            print(f"  {pub_name(pub)} / {entry_key}: {len(summaries)} clusters [{status}]{extra_info}")
        return

    # Run extraction
    from openai import OpenAI
    client = OpenAI()
    total_in = total_out = total_clusters = 0

    by_paper: dict[Path, list] = {}
    for paper_dir, entry_key, table, pub, extraction in entries:
        by_paper.setdefault(paper_dir, []).append((entry_key, table, pub, extraction))

    for paper_dir, group in by_paper.items():
        paper = pub_name(group[0][2])
        extraction_config = group[0][3]

        pdf_path_str = group[0][2].get("papermainpdf", "")
        if pdf_path_str:
            pdf_path = Path(pdf_path_str)
            if not pdf_path.is_absolute():
                pdf_path = Path.cwd() / pdf_path
        else:
            pdfs = list(paper_dir.glob("*.pdf"))
            pdf_path = pdfs[0] if pdfs else None

        if not pdf_path or not pdf_path.exists():
            logger.warning(f"No PDF for {paper}, skipping")
            continue

        logger.info(f"Paper: {paper} ({len(group)} entries)")
        file_ids = [upload_pdf(client, pdf_path)]

        # Upload additional PDFs from extraction config
        additional_pdfs = extraction_config.get("additional_pdfs", [])
        for extra_pdf_str in additional_pdfs:
            extra_path = Path(extra_pdf_str)
            if not extra_path.is_absolute():
                extra_path = Path.cwd() / extra_path
            if extra_path.exists():
                logger.info(f"  Uploading additional PDF: {extra_path.name}")
                file_ids.append(upload_pdf(client, extra_path))
            else:
                logger.warning(f"  Additional PDF not found: {extra_path}")

        scope = extraction_config.get("scope", "per_analysis")

        if scope == "per_paper":
            # Combine all analyses into one call
            to_extract = []
            all_summaries_map = {}
            for entry_key, table, pub, _ext in group:
                if entry_key in list_extraction_files(paper_dir) and not args.force:
                    logger.info(f"  {entry_key}: exists, skipping")
                    continue
                summaries = load_cluster_summaries(table)
                to_extract.append((entry_key, table, pub, summaries))
                all_summaries_map[entry_key] = summaries

            if to_extract:
                tables_and_summaries = [(t, s) for _k, t, _p, s in to_extract]
                try:
                    parsed, usage = extract_paper(
                        client, file_ids, tables_and_summaries,
                        model=args.model, flex=args.flex,
                    )
                    # Collect all expected keys across analyses
                    all_expected = set()
                    for _k, _t, _p, s in to_extract:
                        all_expected.update(s.keys())
                    matched, unmatched = match_cluster_keys(
                        [c.model_dump() for c in parsed.clusters], all_expected,
                    )
                    for c in unmatched:
                        logger.warning(f"    Unmatched: {c.get('name', '?')}")

                    # Distribute matched clusters back to per-entry files
                    for entry_key, table, pub, summaries in to_extract:
                        entry_matched = {k: v for k, v in matched.items() if k in summaries}
                        metadata = {
                            "paper": paper, "doi": pub.get("doi"),
                            "organism": table.get("organism", ""),
                            "entry_key": entry_key, "model": args.model,
                            "extracted_at": datetime.now().isoformat(),
                            "input_tokens": usage["input_tokens"],
                            "output_tokens": usage["output_tokens"],
                        }
                        save_extraction(paper_dir, entry_key, metadata, entry_matched)
                        total_clusters += len(entry_matched)
                        logger.info(f"    {entry_key}: {len(entry_matched)}/{len(summaries)} matched")

                    total_in += usage["input_tokens"]
                    total_out += usage["output_tokens"]

                except Exception as e:
                    logger.error(f"  {paper} (per_paper): FAILED — {e}")
        else:
            # Default: per_analysis extraction
            for entry_key, table, pub, _ext in group:
                if entry_key in list_extraction_files(paper_dir) and not args.force:
                    logger.info(f"  {entry_key}: exists, skipping")
                    continue

                summaries = load_cluster_summaries(table)
                logger.info(f"  {entry_key}: {len(summaries)} clusters")

                try:
                    parsed, usage = extract_analysis(
                        client, file_ids, table, summaries,
                        model=args.model, flex=args.flex,
                    )

                    expected_keys = set(summaries.keys())
                    matched, unmatched = match_cluster_keys(
                        [c.model_dump() for c in parsed.clusters], expected_keys,
                    )

                    for c in unmatched:
                        logger.warning(f"    Unmatched: {c.get('name', '?')}")
                    missing = expected_keys - set(matched.keys())
                    if missing:
                        logger.warning(f"    Missing: {sorted(missing)}")

                    metadata = {
                        "paper": paper, "doi": pub.get("doi"),
                        "organism": table.get("organism", ""),
                        "entry_key": entry_key, "model": args.model,
                        "extracted_at": datetime.now().isoformat(),
                        "input_tokens": usage["input_tokens"],
                        "output_tokens": usage["output_tokens"],
                    }
                    save_extraction(paper_dir, entry_key, metadata, matched)
                    total_in += usage["input_tokens"]
                    total_out += usage["output_tokens"]
                    total_clusters += len(matched)
                    logger.info(f"    {len(matched)}/{len(summaries)} matched")

                except Exception as e:
                    logger.error(f"  {entry_key}: FAILED — {e}")

                time.sleep(2)

        # Clean up uploaded files
        for fid in file_ids:
            try:
                client.files.delete(fid)
            except Exception:
                pass

    logger.info(f"Done: {total_clusters} clusters, {total_in:,} in + {total_out:,} out tokens")


if __name__ == "__main__":
    main()
