#!/usr/bin/env python3
"""One-time migration script: convert paperconfigs from old format to new experiment-based format.

Reads all paperconfigs via paperconfig_utils, groups analyses into experiments,
writes new-format paperconfigs with an `experiments` block and streamlined analyses.

Usage:
    python scripts/migrate_paperconfigs.py --dry-run          # write to /tmp/paperconfig_migration/
    python scripts/migrate_paperconfigs.py --dry-run --output-dir /tmp/my_test/
    python scripts/migrate_paperconfigs.py                    # overwrite originals in place
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Ensure project root is on sys.path so we can import paperconfig_utils
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from multiomics_kg.utils.paperconfig_utils import (
    get_paper_name,
    get_publication,
    get_supplementary_materials,
    iter_analyses,
    load_all_paperconfigs,
    parse_timepoint_hours,
)

# ---------------------------------------------------------------------------
# YAML ordered-dict support (preserve insertion order on dump)
# ---------------------------------------------------------------------------
class _OrderedDumper(yaml.SafeDumper):
    pass

def _dict_representer(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:map", data.items())

_OrderedDumper.add_representer(OrderedDict, _dict_representer)
_OrderedDumper.add_representer(dict, _dict_representer)

# Also handle None -> null cleanly
def _none_representer(dumper, _data):
    return dumper.represent_scalar("tag:yaml.org,2002:null", "null")

_OrderedDumper.add_representer(type(None), _none_representer)


def dump_yaml(data: dict, stream=None) -> str | None:
    return yaml.dump(
        data,
        stream=stream,
        Dumper=_OrderedDumper,
        default_flow_style=False,
        sort_keys=False,
        allow_unicode=True,
        width=120,
    )


# ---------------------------------------------------------------------------
# Time-stripping for treatment_condition
# ---------------------------------------------------------------------------
_TIME_PAREN_RE = re.compile(r"\s*\(days?\s+[\d+]+\)")


def strip_time_from_condition(condition: str) -> str:
    """Remove trailing parenthesized time references like '(day 18)', '(days 60+89)'."""
    return _TIME_PAREN_RE.sub("", condition).strip()


# ---------------------------------------------------------------------------
# Slug generation
# ---------------------------------------------------------------------------
def _make_slug(treatment_type: str, treatment_stripped: str, organism: str, omics_type: str) -> str:
    """Generate a short experiment slug from key fields."""
    # organism short: last word (strain name) or full if single word
    org_parts = organism.strip().split()
    org_short = org_parts[-1].lower() if org_parts else "unknown"

    # condition short: first 4 words, underscored, lowercase
    # For coculture: strip leading "Coculture with" / "Co-culture with" to avoid
    # redundant "coculture_coculture_with_..." slugs
    cond_text = treatment_stripped
    if treatment_type and treatment_type.lower() == "coculture":
        cond_text = re.sub(r"^[Cc]o-?culture\s+with\s+", "", cond_text)
    cond_words = re.sub(r"[^\w\s]", "", cond_text).strip().split()
    cond_short = "_".join(w.lower() for w in cond_words[:4])

    tt = treatment_type.lower().replace(" ", "_") if treatment_type else "unknown"
    omics = omics_type.lower() if omics_type else "unknown"

    return f"{tt}_{cond_short}_{org_short}_{omics}"


def _context_tag(context: str) -> str:
    """Extract a short differentiator from experimental_context for slug dedup.

    Looks for keywords like 'axenic', 'coculture', 'infected', 'uninfected',
    'darkness', 'light' to distinguish experiments that differ only by context.
    """
    ctx = context.lower()
    if "axenic" in ctx:
        return "axenic"
    if "coculture" in ctx or "co-culture" in ctx:
        return "coculture"
    if "uninfected" in ctx:
        return "uninfected"
    if "infected" in ctx:
        return "infected"
    if "darkness" in ctx or "dark" in ctx:
        return "dark"
    if "light" in ctx and "inoculum" not in ctx:
        return "light"
    # Look for concentration/inoculum info (e.g., "5x10^6", "0.5x10^6")
    m = re.search(r"(\d+(?:\.\d+)?)\s*x\s*10", ctx)
    if m:
        return f"inoc_{m.group(1).replace('.', '')}e6"
    # Fallback: first distinctive word
    words = re.sub(r"[^\w\s]", "", ctx).strip().split()
    return words[0] if words else "alt"


def deduplicate_slugs(slugs: list[str], contexts: list[str]) -> list[str]:
    """Deduplicate slugs using context tags instead of numeric suffixes.

    When slugs collide, appends a short tag extracted from experimental_context
    (e.g., '_axenic' vs '_coculture') instead of '_2'.
    Falls back to numeric suffix if context tags also collide.
    """
    # Find which slugs have duplicates
    from collections import Counter
    counts = Counter(slugs)
    duplicated = {s for s, c in counts.items() if c > 1}

    # For duplicated slugs, build tagged versions
    result = []
    seen: dict[str, int] = {}
    for slug, ctx in zip(slugs, contexts):
        if slug not in duplicated:
            result.append(slug)
        else:
            tagged = f"{slug}_{_context_tag(ctx)}"
            # If tagged version also collides, add numeric suffix
            if tagged in seen:
                seen[tagged] += 1
                result.append(f"{tagged}_{seen[tagged]}")
            else:
                seen[tagged] = 1
                result.append(tagged)
    return result


# ---------------------------------------------------------------------------
# Organism short name for experiment name
# ---------------------------------------------------------------------------
def _organism_short(organism: str) -> str:
    """Short organism name for experiment display name."""
    parts = organism.strip().split()
    if len(parts) >= 2:
        # "Prochlorococcus MED4" -> "MED4", "Alteromonas macleodii HOT1A3" -> "HOT1A3"
        # "Alteromonas MIT1002" -> "MIT1002", "Marinobacter" -> "Marinobacter"
        return parts[-1]
    return organism


# ---------------------------------------------------------------------------
# Environmental conditions lookup
# ---------------------------------------------------------------------------
def _get_env_conditions(config: dict) -> dict:
    """Get environmental_conditions block from publication."""
    pub = get_publication(config)
    return pub.get("environmental_conditions") or {}


def _lookup_env_field(env_conditions: dict, condition_id: str | None, field: str) -> str | None:
    """Look up a field from an environmental condition entry."""
    if not condition_id or not env_conditions:
        return None
    entry = env_conditions.get(condition_id)
    if not entry:
        return None
    return entry.get(field)


# ---------------------------------------------------------------------------
# Fields classification
# ---------------------------------------------------------------------------
# Fields that move from analysis to experiment
FIELDS_MOVE_TO_EXPERIMENT = {
    "type",  # -> omics_type
    "name",
    "organism",
    "control_condition",
    "treatment_condition",  # stored time-stripped
    "experimental_context",
    "test_type",
    "treatment_organism",
    "treatment_taxid",
    "treatment_assembly_accession",
}

# Fields removed entirely from analysis
FIELDS_REMOVED = {
    "environmental_control_condition_id",
    "environmental_treatment_condition_id",
}

# Fields that stay on analysis (kept as-is)
FIELDS_KEEP_ON_ANALYSIS = {
    "id",
    "timepoint",
    "name_col",
    "logfc_col",
    "adjusted_p_value_col",
    "pvalue_threshold",
    "logfc_threshold",
    "significance_mode",
    "id_columns",
    "skip_rows",
    "sep",
    "prefiltered",
    "pvalue_asterisk_in_logfc",
}


# ---------------------------------------------------------------------------
# Core migration logic
# ---------------------------------------------------------------------------
def migrate_one_paper(
    path: Path, config: dict, flags: list[str]
) -> dict | None:
    """Migrate a single paperconfig. Returns new config dict, or None to skip."""
    pub = get_publication(config)
    if not pub:
        flags.append(f"  SKIP {path.name}: no publication block")
        return None

    paper_name = get_paper_name(config, path)
    env_conditions = _get_env_conditions(config)
    supp_materials = get_supplementary_materials(config)

    # Collect all analyses with their context
    analysis_records = []
    for table_key, table, analysis in iter_analyses(config):
        organism = analysis.get("organism") or table.get("organism") or ""
        organism = str(organism).strip().strip('"')
        treatment_raw = str(analysis.get("treatment_condition", "")).strip()
        treatment_stripped = strip_time_from_condition(treatment_raw)
        control = str(analysis.get("control_condition", "")).strip()
        context = str(analysis.get("experimental_context", "")).strip()
        omics_type = str(analysis.get("type", "")).strip()
        test_type = str(analysis.get("test_type", "")).strip()

        # Use env condition IDs for grouping when available (collapses
        # Weissberg-style phase names that all point to the same condition).
        # Fall back to text-based grouping when IDs are absent.
        env_treat_id = analysis.get("environmental_treatment_condition_id")
        env_ctrl_id = analysis.get("environmental_control_condition_id")
        treatment_key = env_treat_id if env_treat_id else treatment_stripped
        control_key = env_ctrl_id if env_ctrl_id else control

        grouping_key = (organism, treatment_key, control_key, context, omics_type, test_type)

        analysis_records.append({
            "table_key": table_key,
            "table": table,
            "analysis": analysis,
            "organism": organism,
            "treatment_stripped": treatment_stripped,
            "control": control,
            "context": context,
            "omics_type": omics_type,
            "test_type": test_type,
            "grouping_key": grouping_key,
        })

    if not analysis_records:
        flags.append(f"  SKIP {paper_name}: no analyses found")
        return None

    # Group analyses by key
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for rec in analysis_records:
        groups[rec["grouping_key"]].append(rec)

    # Generate experiment entries
    raw_slugs = []
    raw_contexts = []
    group_keys_ordered = list(groups.keys())
    for key in group_keys_ordered:
        recs = groups[key]
        first = recs[0]
        first_analysis = first["analysis"]
        # Determine treatment_type
        treatment_type = _determine_treatment_type(first, env_conditions)
        # Use env condition name for slug when grouped by env ID
        env_treat_id = first_analysis.get("environmental_treatment_condition_id")
        if env_treat_id and env_treat_id in env_conditions:
            slug_treatment = env_conditions[env_treat_id].get("name", first["treatment_stripped"])
        else:
            slug_treatment = first["treatment_stripped"]
        raw_slugs.append(_make_slug(treatment_type, slug_treatment, first["organism"], first["omics_type"]))
        raw_contexts.append(first["context"])

    slugs = deduplicate_slugs(raw_slugs, raw_contexts)

    experiments = OrderedDict()
    slug_by_key: dict[tuple, str] = {}

    for i, key in enumerate(group_keys_ordered):
        slug = slugs[i]
        slug_by_key[key] = slug
        recs = groups[key]
        first_analysis = recs[0]["analysis"]
        first_rec = recs[0]

        treatment_type = _determine_treatment_type(first_rec, env_conditions)
        medium = _get_medium(first_rec, env_conditions)
        temperature = _get_temperature(first_rec, env_conditions)
        light_condition = _get_light_condition(first_rec, env_conditions)
        light_intensity = _get_light_intensity(first_rec, env_conditions)

        # Determine experiment-level treatment/control labels.
        # When grouped by env condition ID, use the condition name (more general
        # than any single phase name like "Decline" or "Long-term starvation").
        env_treat_id = first_analysis.get("environmental_treatment_condition_id")
        env_ctrl_id = first_analysis.get("environmental_control_condition_id")
        if env_treat_id and env_treat_id in env_conditions:
            exp_treatment = env_conditions[env_treat_id].get("name", first_rec["treatment_stripped"])
        else:
            exp_treatment = first_rec["treatment_stripped"]
        if env_ctrl_id and env_ctrl_id in env_conditions:
            exp_control = env_conditions[env_ctrl_id].get("name", first_rec["control"])
        else:
            exp_control = first_rec["control"]

        # Build experiment name
        org_short = _organism_short(first_rec["organism"])
        exp_name = f"{org_short} {exp_treatment} vs {exp_control} ({first_rec['omics_type']})"

        exp = OrderedDict()
        exp["name"] = exp_name
        exp["organism"] = first_rec["organism"]
        exp["treatment_condition"] = exp_treatment
        exp["control_condition"] = exp_control
        exp["experimental_context"] = first_rec["context"]
        exp["omics_type"] = first_rec["omics_type"]
        exp["test_type"] = first_rec["test_type"]
        exp["treatment_type"] = treatment_type
        exp["medium"] = medium
        exp["temperature"] = temperature
        exp["light_condition"] = light_condition
        exp["light_intensity"] = light_intensity

        # Coculture fields
        t_org = first_analysis.get("treatment_organism")
        if t_org:
            exp["treatment_organism"] = t_org
            t_taxid = first_analysis.get("treatment_taxid")
            if t_taxid is not None:
                exp["treatment_taxid"] = t_taxid
            t_acc = first_analysis.get("treatment_assembly_accession")
            if t_acc:
                exp["treatment_assembly_accession"] = t_acc

        experiments[slug] = exp

        # Flag reports
        if len(recs) == 1:
            flags.append(f"  INFO {paper_name}: experiment '{slug}' has only 1 analysis (single-point)")
        if treatment_type == "NEEDS_CURATION":
            flags.append(f"  WARN {paper_name}: experiment '{slug}' has treatment_type=NEEDS_CURATION")

    # Rebuild supplementary_materials
    new_supp = OrderedDict()
    total_analyses_before = 0
    total_analyses_after = 0

    for table_key, table in supp_materials.items():
        table_type = table.get("type", "csv")

        # Preserve non-csv entries unchanged
        if table_type in ("id_translation", "annotation_gff"):
            new_supp[table_key] = table
            continue

        if table_type != "csv":
            new_supp[table_key] = table
            continue

        old_analyses = table.get("statistical_analyses", [])
        total_analyses_before += len(old_analyses)

        new_analyses = []
        for analysis in old_analyses:
            # Find the grouping key for this analysis
            organism = analysis.get("organism") or table.get("organism") or ""
            organism = str(organism).strip().strip('"')
            treatment_raw = str(analysis.get("treatment_condition", "")).strip()
            treatment_stripped = strip_time_from_condition(treatment_raw)
            control = str(analysis.get("control_condition", "")).strip()
            context = str(analysis.get("experimental_context", "")).strip()
            omics_type = str(analysis.get("type", "")).strip()
            test_type = str(analysis.get("test_type", "")).strip()

            # Must reconstruct grouping key the same way as the first pass
            env_treat_id = analysis.get("environmental_treatment_condition_id")
            env_ctrl_id = analysis.get("environmental_control_condition_id")
            treatment_key = env_treat_id if env_treat_id else treatment_stripped
            control_key = env_ctrl_id if env_ctrl_id else control

            grouping_key = (organism, treatment_key, control_key, context, omics_type, test_type)
            slug = slug_by_key.get(grouping_key)

            if not slug:
                flags.append(f"  ERROR {paper_name}: analysis '{analysis.get('id')}' has no matching experiment!")
                # Keep analysis unchanged as fallback
                new_analyses.append(analysis)
                total_analyses_after += 1
                continue

            # Build new analysis entry
            new_analysis = OrderedDict()
            new_analysis["id"] = analysis.get("id")
            new_analysis["experiment"] = slug

            # Timepoint
            timepoint = analysis.get("timepoint")
            if timepoint is not None:
                new_analysis["timepoint"] = timepoint
                tp_hours = parse_timepoint_hours(str(timepoint))
                new_analysis["timepoint_hours"] = tp_hours
                if tp_hours is None and timepoint:
                    flags.append(f"  WARN {paper_name}: unparseable timepoint_hours for '{timepoint}' (analysis {analysis.get('id')})")
            else:
                new_analysis["timepoint_hours"] = None

            # Keep remaining per-analysis fields
            for field in (
                "name_col", "logfc_col", "adjusted_p_value_col",
                "pvalue_threshold", "logfc_threshold", "significance_mode",
                "prefiltered", "pvalue_asterisk_in_logfc",
                "id_columns", "skip_rows", "sep",
            ):
                val = analysis.get(field)
                if val is not None:
                    new_analysis[field] = val

            new_analyses.append(new_analysis)
            total_analyses_after += 1

        # Build new table entry
        new_table = OrderedDict()
        new_table["type"] = "csv"
        if table.get("filename"):
            new_table["filename"] = table["filename"]
        if table.get("organism"):
            new_table["organism"] = table["organism"]
        if table.get("skip_rows") is not None:
            new_table["skip_rows"] = table["skip_rows"]
        if table.get("sep") is not None:
            new_table["sep"] = table["sep"]
        if table.get("id_columns"):
            new_table["id_columns"] = table["id_columns"]
        if table.get("product_columns"):
            new_table["product_columns"] = table["product_columns"]
        new_table["statistical_analyses"] = new_analyses

        new_supp[table_key] = new_table

    # Round-trip validation
    if total_analyses_before != total_analyses_after:
        flags.append(
            f"  ERROR {paper_name}: analysis count mismatch! before={total_analyses_before} after={total_analyses_after}"
        )

    # Build new config
    new_config = OrderedDict()
    new_pub = OrderedDict()
    new_pub["papername"] = pub.get("papername")
    if pub.get("papermainpdf"):
        new_pub["papermainpdf"] = pub["papermainpdf"]
    new_pub["experiments"] = experiments
    new_pub["supplementary_materials"] = new_supp
    new_config["publication"] = new_pub

    return new_config


def _determine_treatment_type(rec: dict, env_conditions: dict) -> str:
    """Determine treatment_type for an experiment group."""
    analysis = rec["analysis"]

    # Priority 1: coculture detection
    if analysis.get("treatment_organism"):
        return "coculture"

    # Priority 2: environmental_conditions lookup
    for ref_field in ("environmental_treatment_condition_id", "environmental_control_condition_id"):
        cond_id = analysis.get(ref_field)
        ct = _lookup_env_field(env_conditions, cond_id, "condition_type")
        if ct:
            return ct

    return "NEEDS_CURATION"


def _get_medium(rec: dict, env_conditions: dict) -> str:
    analysis = rec["analysis"]
    for ref_field in ("environmental_treatment_condition_id", "environmental_control_condition_id"):
        cond_id = analysis.get(ref_field)
        val = _lookup_env_field(env_conditions, cond_id, "medium")
        if val:
            return val
    return "NEEDS_CURATION"


def _get_temperature(rec: dict, env_conditions: dict) -> str:
    analysis = rec["analysis"]
    for ref_field in ("environmental_treatment_condition_id", "environmental_control_condition_id"):
        cond_id = analysis.get(ref_field)
        val = _lookup_env_field(env_conditions, cond_id, "temperature")
        if val:
            return val
    return "NEEDS_CURATION"


def _get_light_condition(rec: dict, env_conditions: dict) -> str:
    analysis = rec["analysis"]
    for ref_field in ("environmental_treatment_condition_id", "environmental_control_condition_id"):
        cond_id = analysis.get(ref_field)
        val = _lookup_env_field(env_conditions, cond_id, "light_condition")
        if val:
            return val
    return "NEEDS_CURATION"


def _get_light_intensity(rec: dict, env_conditions: dict) -> str:
    analysis = rec["analysis"]
    for ref_field in ("environmental_treatment_condition_id", "environmental_control_condition_id"):
        cond_id = analysis.get(ref_field)
        val = _lookup_env_field(env_conditions, cond_id, "light_intensity")
        if val:
            return val
    return "NEEDS_CURATION"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Migrate paperconfigs to new experiment-based format")
    parser.add_argument("--dry-run", action="store_true", help="Write to output dir instead of overwriting originals")
    parser.add_argument("--output-dir", type=str, default="/tmp/paperconfig_migration/", help="Output directory for dry-run mode")
    args = parser.parse_args()

    all_configs = load_all_paperconfigs()
    print(f"Loaded {len(all_configs)} paperconfigs", file=sys.stderr)

    total_experiments = 0
    total_analyses = 0
    papers_processed = 0
    papers_skipped = 0
    flags: list[str] = []

    for path, config in all_configs:
        paper_name = get_paper_name(config, path)
        result = migrate_one_paper(path, config, flags)

        if result is None:
            papers_skipped += 1
            continue

        papers_processed += 1
        new_pub = result.get("publication", {})
        experiments = new_pub.get("experiments", {})
        total_experiments += len(experiments)

        # Count analyses in output
        for _tk, table in new_pub.get("supplementary_materials", {}).items():
            if isinstance(table, dict) and table.get("type") == "csv":
                total_analyses += len(table.get("statistical_analyses", []))

        # Write output
        if args.dry_run:
            # Preserve relative path structure under output dir
            rel = path.relative_to(PROJECT_ROOT)
            out_path = Path(args.output_dir) / rel
            out_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            out_path = path

        with open(out_path, "w") as f:
            dump_yaml(result, stream=f)

        print(f"  {paper_name}: {len(experiments)} experiments -> {out_path}", file=sys.stderr)

    # Summary
    print("\n" + "=" * 60, file=sys.stderr)
    print(f"Papers processed: {papers_processed}", file=sys.stderr)
    print(f"Papers skipped:   {papers_skipped}", file=sys.stderr)
    print(f"Total experiments: {total_experiments}", file=sys.stderr)
    print(f"Total analyses:   {total_analyses}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)

    if flags:
        print("\nFlags:", file=sys.stderr)
        for flag in flags:
            print(flag, file=sys.stderr)


if __name__ == "__main__":
    main()
