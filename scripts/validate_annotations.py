#!/usr/bin/env python3
"""
Validate functional annotations in paper CSVs against gene_annotations_merged.json.

For each paper / CSV file / strain, compares annotation/product/description columns
in the supplementary CSV against the canonical product annotation in
gene_annotations_merged.json, reporting match rates.

Uses exact + token-overlap matching for all genes, and optionally batches a sample
of non-matching pairs through an LLM for qualitative mismatch analysis.

Usage:
    uv run python scripts/validate_annotations.py
    uv run python scripts/validate_annotations.py --papers "Biller 2018" "coe 2024"
    uv run python scripts/validate_annotations.py --no-llm
    uv run python scripts/validate_annotations.py --paperconfig-list data/Prochlorococcus/papers_and_supp/paperconfig_files.txt
"""

import argparse
import json
import logging
import os
import random
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

import dotenv
import pandas as pd
import yaml

dotenv.load_dotenv()

# Project root is one level above scripts/
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from multiomics_kg.utils.gene_id_utils import (
    SKIP_PATTERNS,
    build_id_lookup,
    get_genome_dir,
    load_gene_annotations,
    map_gene_id,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

# ─── Annotation column detection ──────────────────────────────────────────────

# Keywords whose presence in a column name suggests it contains functional text
ANNOTATION_COL_KEYWORDS = [
    "product", "annotation", "description", "function", "note",
    "genbank", "rast", "role", "gene annotation", "gene name",
    "gene_name", "gene description", "definition", "inferred function",
    "protein name", "putative function",
]

# Pattern to skip clearly non-annotation columns
_SKIP_COL_RE = re.compile(
    r"(logf[cd]|log2|fold.?change|p.?value|padj|adj\.?p|fdr|lfc|"
    r"\bstart\b|\bend\b|length|position|coord|strand|"
    r"timepoint|condition|treatment|control|diffexpressed|"
    r"kegg.class|cog.class|pathway|domain|cluster|temperature|co2|"
    r"eggnog|tigrfam|pfam)",
    re.IGNORECASE,
)

# Reference annotation fields from gene_annotations_merged.json (priority order)
REF_FIELDS = ["product", "product_cyanorak"]


def detect_annotation_columns(df: pd.DataFrame, name_col: str) -> list[str]:
    """Return columns in *df* that likely contain functional descriptions (not gene symbols)."""
    exclude = {name_col}
    found = []
    for col in df.columns:
        if col in exclude:
            continue
        col_lower = col.lower().strip()
        if _SKIP_COL_RE.search(col_lower):
            continue
        if not any(kw in col_lower for kw in ANNOTATION_COL_KEYWORDS):
            continue
        # Verify the column has non-numeric string content
        vals = df[col].dropna().astype(str).str.strip()
        vals = vals[vals != ""].head(30)
        if vals.empty:
            continue
        try:
            pd.to_numeric(vals)
            continue  # purely numeric → skip
        except (ValueError, TypeError):
            pass
        # Skip columns whose values look like gene symbols (short identifiers, no spaces)
        # Gene symbols: "dnaN", "purL", "gyrA" — typically < 12 chars, rarely contain spaces
        nonempty = vals[vals.str.len() > 0]
        avg_len = nonempty.str.len().mean()
        has_spaces_frac = nonempty.str.contains(" ").mean()
        if avg_len < 12 and has_spaces_frac < 0.2:
            # Looks like gene symbols, not descriptions — skip
            continue
        found.append(col)
    return found


# ─── Text normalization & matching ────────────────────────────────────────────

_PUNCT_RE = re.compile(r"[,;:\-/()\[\]{}]")
_SPACE_RE = re.compile(r"\s+")
# Words that add little semantic meaning for comparison
_STOP = frozenset({"protein", "the", "a", "an", "of", "and", "in", "to",
                   "for", "is", "are", "with", "from", "by", "type", "like",
                   "related", "domain", "family", "subunit", "component"})


def normalize(text: str) -> str:
    """Lowercase, strip punctuation and extra whitespace."""
    if not text or (isinstance(text, float)):
        return ""
    text = str(text).strip().lower()
    text = _PUNCT_RE.sub(" ", text)
    text = _SPACE_RE.sub(" ", text)
    return text.strip()


def token_jaccard(a: str, b: str) -> float:
    """Jaccard overlap of meaningful tokens in two annotation strings."""
    ta = set(normalize(a).split()) - _STOP
    tb = set(normalize(b).split()) - _STOP
    if not ta or not tb:
        return 0.0
    return len(ta & tb) / len(ta | tb)


def is_match(csv_annot: str, ref_annot: str, threshold: float = 0.5) -> tuple[bool, bool]:
    """Return (exact_match, token_match) for a pair of annotations.

    ``token_match`` is True if Jaccard token overlap ≥ threshold.
    """
    na = normalize(csv_annot)
    nb = normalize(ref_annot)
    exact = (na == nb) and bool(na)
    token = token_jaccard(csv_annot, ref_annot) >= threshold
    return exact, token


def get_ref_annotation(entry: dict) -> str:
    """Pick the best reference annotation from a gene_annotations_merged entry."""
    for field in REF_FIELDS:
        val = entry.get(field)
        if val and isinstance(val, str) and val.strip():
            return val.strip()
    return ""


# ─── LLM semantic batch analysis ──────────────────────────────────────────────

def _init_llm(model_name: str = "gpt-4.1-nano"):
    """Initialize a LangChain chat model (falls back gracefully)."""
    try:
        from langchain.chat_models import init_chat_model
        return init_chat_model(model=model_name, temperature=0.0)
    except Exception as exc:
        logger.warning(f"LLM init failed ({exc}); skipping LLM analysis.")
        return None


def llm_analyze_mismatches(
    mismatch_pairs: list[tuple[str, str]],
    llm,
    sample_size: int = 30,
    batch_size: int = 20,
) -> dict:
    """Ask an LLM to classify a random sample of mismatching annotation pairs.

    Returns a summary dict with counts and example commentary.
    """
    if not llm or not mismatch_pairs:
        return {}

    from langchain_core.messages import HumanMessage

    sample = random.sample(mismatch_pairs, min(sample_size, len(mismatch_pairs)))

    prompt_header = (
        "You are a bioinformatics expert. Below are pairs of gene/protein annotations "
        "from two different sources for the same gene. The first is from a paper's "
        "supplementary CSV; the second is the canonical NCBI/Cyanorak annotation.\n\n"
        "Classify each pair as one of:\n"
        "  - 'equivalent': same function, just different wording or specificity\n"
        "  - 'partial': related but notably different in specificity or scope\n"
        "  - 'different': describe different functions\n"
        "  - 'generic': one entry is 'hypothetical protein' or similarly uninformative\n\n"
        "Return a JSON object with:\n"
        "  counts: {equivalent, partial, different, generic}\n"
        "  summary: one sentence explaining the main reason for mismatches\n"
        "  examples: up to 3 representative pairs with their classification\n\n"
        "Pairs:\n"
    )

    # Build numbered list
    lines = []
    for i, (csv_a, ref_a) in enumerate(sample, 1):
        lines.append(f"{i}. CSV: \"{csv_a}\" | REF: \"{ref_a}\"")

    prompt = prompt_header + "\n".join(lines)

    try:
        response = llm.invoke([HumanMessage(content=prompt)])
        text = response.content.strip()
        # Extract JSON from response
        json_match = re.search(r"\{.*\}", text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
        return {"summary": text[:300]}
    except Exception as exc:
        logger.warning(f"LLM call failed: {exc}")
        return {}


# ─── CSV loading ──────────────────────────────────────────────────────────────

def load_csv(filename: str, skip_rows: int = 0) -> Optional[pd.DataFrame]:
    """Load a CSV, skipping the first *skip_rows* rows before the real header.

    Matches the omics adapter's behaviour: ``pd.read_csv(file, skiprows=N)``
    skips the first N rows and treats the next row as the column header.
    """
    path = PROJECT_ROOT / filename
    if not path.exists():
        logger.warning(f"CSV not found: {path}")
        return None
    try:
        return pd.read_csv(path, skiprows=skip_rows or None, encoding="utf-8-sig")
    except Exception as exc:
        logger.warning(f"Could not read {path}: {exc}")
        return None


# ─── Per-CSV validation ────────────────────────────────────────────────────────

class CsvResult:
    """Holds per-column match stats for one CSV file."""

    def __init__(self, paper: str, csv_file: str, organism: str, strain: str):
        self.paper = paper
        self.csv_file = csv_file
        self.organism = organism
        self.strain = strain
        self.columns: dict[str, dict] = {}  # col_name → stats dict
        self.error: Optional[str] = None

    def add_column(self, col: str, total: int, exact: int, token: int,
                   mismatches: list, no_ref: int, no_csv: int):
        self.columns[col] = dict(
            total=total, exact=exact, token=token,
            mismatches=mismatches, no_ref=no_ref, no_csv=no_csv,
        )

    def has_data(self) -> bool:
        return bool(self.columns)


def validate_csv_file(
    paper_name: str,
    csv_filename: str,
    analyses: list[dict],
    supp_config: dict,
    use_llm: bool,
    llm,
) -> CsvResult:
    """Validate one CSV file across its annotation columns."""
    # Pull the first analysis for name_col / organism (all share the same CSV)
    first = analyses[0]
    name_col = first.get("name_col", "")
    organism = first.get("organism", "")
    skip_rows = supp_config.get("skip_rows", 0)

    result = CsvResult(paper_name, csv_filename, organism, "")

    # Resolve genome dir
    genome_dir = get_genome_dir(organism, str(PROJECT_ROOT))
    if not genome_dir:
        result.error = f"No genome dir for organism '{organism}'"
        return result

    strain = Path(genome_dir).name
    result.strain = strain

    # Load reference annotations
    annotations = load_gene_annotations(genome_dir)
    if not annotations:
        result.error = f"gene_annotations_merged.json missing for {strain}"
        return result

    # Build ID lookup
    lookup, locus_tags, supp_keys = build_id_lookup(genome_dir)
    if not lookup:
        result.error = f"Could not build ID lookup for {strain}"
        return result

    # Load CSV
    df = load_csv(csv_filename, skip_rows)
    if df is None:
        result.error = "CSV file not found or unreadable"
        return result

    # Normalise column name  (the name_col declared in the YAML may differ from actual header)
    col_map = {c.strip(): c for c in df.columns}
    actual_name_col = col_map.get(name_col, name_col)
    if actual_name_col not in df.columns:
        # Try case-insensitive match
        for c in df.columns:
            if c.strip().lower() == name_col.strip().lower():
                actual_name_col = c
                break
    if actual_name_col not in df.columns:
        result.error = f"name_col '{name_col}' not in CSV columns: {list(df.columns)[:10]}"
        return result

    # Detect annotation columns
    annot_cols = detect_annotation_columns(df, actual_name_col)
    if not annot_cols:
        result.error = "No annotation columns detected in CSV"
        return result

    # For each gene row, map ID → locus_tag, gather annotations
    gene_rows: list[dict] = []  # {locus_tag, ref_annot, csv_col_vals}
    for _, row in df.iterrows():
        raw_id = str(row[actual_name_col]).strip()
        if not raw_id or raw_id == "nan":
            continue
        if SKIP_PATTERNS.match(raw_id):
            continue

        locus_tag, _ = map_gene_id(raw_id, lookup, locus_tags, supp_keys)
        if not locus_tag:
            continue

        entry = annotations.get(locus_tag, {})
        ref_annot = get_ref_annotation(entry)

        csv_vals = {}
        for col in annot_cols:
            val = row.get(col)
            if val and str(val).strip() not in ("", "nan", "NaN"):
                csv_vals[col] = str(val).strip()
        gene_rows.append({"locus_tag": locus_tag, "ref": ref_annot, "csv": csv_vals})

    # Per-column stats
    for col in annot_cols:
        total = 0
        exact_cnt = 0
        token_cnt = 0
        no_ref = 0
        no_csv = 0
        mismatch_pairs: list[tuple[str, str]] = []

        for row in gene_rows:
            csv_a = row["csv"].get(col, "")
            ref_a = row["ref"]

            if not csv_a:
                no_csv += 1
                continue
            if not ref_a:
                no_ref += 1
                continue

            total += 1
            ex, tk = is_match(csv_a, ref_a)
            if ex:
                exact_cnt += 1
                token_cnt += 1
            elif tk:
                token_cnt += 1
            else:
                mismatch_pairs.append((csv_a, ref_a))

        # LLM analysis of a sample of mismatches
        llm_result = {}
        if use_llm and mismatch_pairs:
            llm_result = llm_analyze_mismatches(mismatch_pairs, llm)

        result.add_column(
            col=col,
            total=total,
            exact=exact_cnt,
            token=token_cnt,
            mismatches=mismatch_pairs[:5],  # keep only a small sample for the report
            no_ref=no_ref,
            no_csv=no_csv,
        )
        result.columns[col]["llm"] = llm_result

    return result


# ─── Paperconfig parsing ───────────────────────────────────────────────────────

def load_paperconfig(config_path: str) -> Optional[dict]:
    """Load a paperconfig.yaml."""
    path = PROJECT_ROOT / config_path
    if not path.exists():
        logger.warning(f"Paperconfig not found: {path}")
        return None
    with open(path) as f:
        return yaml.safe_load(f)


def iter_csvs(config_data: dict):
    """Yield (csv_filename, supp_config, analyses) for each *unique* CSV in a paperconfig.

    Multiple supp table entries may point to the same filename; we deduplicate
    so each physical CSV is validated only once.
    """
    publication = config_data.get("publication", {})
    supps = publication.get("supplementary_materials", {})
    if not isinstance(supps, dict):
        return

    seen_files: set[str] = set()
    for supp_key, supp_val in supps.items():
        if not isinstance(supp_val, dict):
            continue
        if supp_val.get("type") != "csv":
            continue
        filename = supp_val.get("filename", "")
        if not filename or filename in seen_files:
            continue
        seen_files.add(filename)
        analyses = supp_val.get("statistical_analyses", [])
        if not analyses:
            continue
        yield filename, supp_val, analyses


# ─── Report formatting ─────────────────────────────────────────────────────────

def pct(n: int, total: int) -> str:
    if total == 0:
        return "N/A"
    return f"{n/total*100:.0f}%"


def print_report(all_results: list[CsvResult]) -> None:
    """Print the validation report to stdout."""
    print("\n" + "=" * 72)
    print("ANNOTATION VALIDATION REPORT")
    print("=" * 72)

    papers = defaultdict(list)
    for r in all_results:
        papers[r.paper].append(r)

    grand_total = grand_token = grand_gene = 0

    for paper, results in sorted(papers.items()):
        print(f"\nPaper: {paper}")
        print("-" * 60)

        for res in results:
            csv_short = Path(res.csv_file).name
            print(f"  CSV:      {csv_short}")
            print(f"  Organism: {res.organism or '(unknown)'}  [strain: {res.strain or '?'}]")

            if res.error:
                print(f"  ERROR: {res.error}")
                continue

            if not res.has_data():
                print("  (no annotation columns found)")
                continue

            for col, stats in res.columns.items():
                t = stats["total"]
                ex = stats["exact"]
                tk = stats["token"]
                no_ref = stats["no_ref"]
                no_csv = stats["no_csv"]
                misses = stats["mismatches"]
                llm = stats.get("llm", {})

                broad_match = tk  # token overlap ≥ 0.5

                print(f"\n  Annotation column: \"{col}\"")
                print(f"    Genes matched to reference:  {t + no_ref + no_csv} total gene IDs resolved")
                if no_csv:
                    print(f"    ├─ No CSV annotation:        {no_csv}")
                if no_ref:
                    print(f"    ├─ No reference annotation:  {no_ref}")
                print(f"    ├─ Comparable pairs:         {t}")
                print(f"    ├─ Exact match:              ({ex}/{t}) {pct(ex, t)}")
                print(f"    └─ Broad match (token ≥50%): ({broad_match}/{t}) {pct(broad_match, t)}")
                print(f"    Summary: ({broad_match}/{t}) {pct(broad_match, t)} match")

                grand_total += t
                grand_token += broad_match
                grand_gene += t + no_ref + no_csv

                if llm:
                    counts = llm.get("counts", {})
                    summary = llm.get("summary", "")
                    examples = llm.get("examples", [])
                    print(f"\n    LLM mismatch analysis (sample of ≤30 non-matching pairs):")
                    if counts:
                        total_llm = sum(counts.values()) or 1
                        for cat, n in counts.items():
                            print(f"      {cat:12s}: {n} ({pct(n, total_llm)})")
                    if summary:
                        print(f"      Summary: {summary}")
                    if examples:
                        print("      Examples:")
                        for ex_item in examples[:3]:
                            if isinstance(ex_item, dict):
                                cls = ex_item.get("classification", ex_item.get("class", ""))
                                csv_e = ex_item.get("csv", ex_item.get("CSV", ""))
                                ref_e = ex_item.get("ref", ex_item.get("REF", ""))
                                print(f"        [{cls}] \"{csv_e}\" → \"{ref_e}\"")

                if misses and not llm:
                    print(f"\n    Sample mismatches (up to 5):")
                    for csv_a, ref_a in misses[:5]:
                        print(f"      CSV: \"{csv_a}\"")
                        print(f"      REF: \"{ref_a}\"")
                        print()

    print("\n" + "=" * 72)
    print(f"OVERALL: ({grand_token}/{grand_total}) {pct(grand_token, grand_total)} broad match")
    print(f"         ({grand_gene} gene IDs resolved across all CSVs)")
    print("=" * 72 + "\n")


# ─── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--paperconfig-list",
        default="data/Prochlorococcus/papers_and_supp/paperconfig_files.txt",
        help="Path to the file listing paperconfig.yaml paths",
    )
    parser.add_argument(
        "--papers", nargs="+",
        help="Restrict to specific paper names (substring match)",
    )
    parser.add_argument(
        "--no-llm", action="store_true",
        help="Skip LLM semantic analysis of mismatches",
    )
    parser.add_argument(
        "--llm-model", default="gpt-4.1-nano",
        help="LangChain model string for LLM analysis (default: gpt-4.1-nano)",
    )
    parser.add_argument(
        "--match-threshold", type=float, default=0.5,
        help="Jaccard token overlap threshold for broad match (default: 0.5)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
    )
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    # Load LLM if requested
    llm = None
    use_llm = not args.no_llm
    if use_llm:
        llm = _init_llm(args.llm_model)
        if llm is None:
            use_llm = False

    # Read list of paperconfigs
    list_path = PROJECT_ROOT / args.paperconfig_list
    if not list_path.exists():
        print(f"ERROR: paperconfig list not found: {list_path}", file=sys.stderr)
        sys.exit(1)

    config_paths = [
        line.strip() for line in list_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]

    if args.papers:
        filters = [p.lower() for p in args.papers]
        config_paths = [
            p for p in config_paths
            if any(f in p.lower() for f in filters)
        ]

    if not config_paths:
        print("No paperconfigs to process.", file=sys.stderr)
        sys.exit(1)

    all_results: list[CsvResult] = []

    for config_path in config_paths:
        config_data = load_paperconfig(config_path)
        if not config_data:
            continue

        paper_name = (
            config_data.get("publication", {}).get("papername", "")
            or Path(config_path).parent.name
        )
        print(f"Processing: {paper_name} ...", end="  ", flush=True)

        paper_results = []
        for csv_file, supp_config, analyses in iter_csvs(config_data):
            res = validate_csv_file(
                paper_name=paper_name,
                csv_filename=csv_file,
                analyses=analyses,
                supp_config=supp_config,
                use_llm=use_llm,
                llm=llm,
            )
            paper_results.append(res)
            all_results.append(res)

        cols_found = sum(len(r.columns) for r in paper_results)
        errors = sum(1 for r in paper_results if r.error)
        print(f"{len(paper_results)} CSV(s), {cols_found} annotation col(s), {errors} error(s)")

    print_report(all_results)


if __name__ == "__main__":
    main()
