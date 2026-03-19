"""Shared utilities for reading and traversing paperconfig.yaml files."""
import re
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).parent.parent.parent
PAPERCONFIG_FILES_TXT = (
    PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
)


# ─── Loading ──────────────────────────────────────────────────────────


def load_paperconfig(path: Path) -> dict:
    """Load a single paperconfig.yaml file."""
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_all_paperconfigs(
    list_file: Path = PAPERCONFIG_FILES_TXT,
) -> list[tuple[Path, dict]]:
    """Load all paperconfigs from paperconfig_files.txt.

    Returns list of (path, config_dict). Skips comments, blank lines,
    missing files (with warning).
    """
    results = []
    with open(list_file) as f:
        for line in f:
            path_str = line.strip()
            if not path_str or path_str.startswith("#"):
                continue
            path = PROJECT_ROOT / path_str
            if not path.exists():
                print(f"  [warn] paperconfig not found: {path}")
                continue
            results.append((path, load_paperconfig(path)))
    return results


# ─── Traversal helpers ───────────────────────────────────────────────


def get_publication(config: dict) -> dict:
    """Get the publication block."""
    return config.get("publication", {})


def get_paper_name(config: dict, fallback_path: Path | None = None) -> str:
    """Get paper name, with fallback to directory name."""
    name = get_publication(config).get("papername")
    if name:
        return name
    if fallback_path:
        return fallback_path.parent.name
    return "unknown"


def get_supplementary_materials(config: dict) -> dict:
    """Get supplementary_materials dict.

    Checks publication.supplementary_materials first, then falls back to
    top-level supplementary_materials (for strain-level resource configs
    like MIT9313_resources that have no publication block).
    """
    pub = get_publication(config)
    return pub.get("supplementary_materials") or config.get("supplementary_materials") or {}


def iter_csv_tables(config: dict):
    """Yield (table_key, table_config) for csv-type supplementary tables."""
    for key, table in get_supplementary_materials(config).items():
        if table.get("type", "csv") == "csv":
            yield key, table


def iter_analyses(config: dict):
    """Yield (table_key, table_config, analysis) for all analyses across all csv tables."""
    for table_key, table in iter_csv_tables(config):
        for analysis in table.get("statistical_analyses", []):
            yield table_key, table, analysis


def get_organism_for_entry(config: dict, entry: dict) -> str | None:
    """Get organism for a supplementary table entry.

    Checks in order:
    1. Direct 'organism' on entry (id_translation, annotation_gff, table-level)
    2. First analysis's experiment block organism (new format)
    3. First analysis's direct 'organism' field (old format fallback)
    """
    # Direct organism on entry (id_translation, annotation_gff, or table-level)
    org = entry.get("organism")
    if org:
        return str(org).strip().strip('"')
    # From first analysis's experiment block (new format)
    for analysis in entry.get("statistical_analyses", []):
        try:
            return get_organism_for_analysis(config, analysis)
        except (ValueError, KeyError):
            pass
        # Old format fallback: direct organism on analysis
        org = analysis.get("organism")
        if org:
            return str(org).strip().strip('"')
    return None


# ─── Experiment lookup (new format) ──────────────────────────────────


def get_experiments(config: dict) -> dict:
    """Get the experiments block from a paperconfig."""
    return get_publication(config).get("experiments", {})


def get_experiment_for_analysis(config: dict, analysis: dict) -> dict:
    """Look up the experiment definition for a given analysis entry.

    Raises ValueError if the analysis has no 'experiment' reference or
    the reference points to an unknown experiment key.
    """
    exp_key = analysis.get("experiment")
    if not exp_key:
        raise ValueError(
            f"Analysis '{analysis.get('id')}' missing 'experiment' reference"
        )
    experiments = get_experiments(config)
    if exp_key not in experiments:
        raise ValueError(
            f"Analysis '{analysis.get('id')}' references unknown "
            f"experiment '{exp_key}'"
        )
    return experiments[exp_key]


def get_organism_for_analysis(config: dict, analysis: dict) -> str:
    """Get organism for an analysis (from its experiment block)."""
    return get_experiment_for_analysis(config, analysis)["organism"]


# ─── Timepoint normalization ─────────────────────────────────────────


def parse_timepoint_hours(timepoint: str | None) -> float | None:
    """Parse a timepoint string to hours. Returns None for unparseable values.

    Handles: "4h", "0.5h", "-12h", "day 18", "Day 2", "50h (P added)",
    "0.5h post-inoculation", "1h extended darkness (36h)".
    Returns None for: "R (rescue: ...)", "days 60+89", None.
    """
    if not timepoint:
        return None
    tp = timepoint.strip()

    # "days N+M" — pooled, unparseable
    if tp.lower().startswith("days ") and "+" in tp:
        return None
    # "R (rescue: ...)" — non-numeric
    if tp.startswith("R ") or tp == "R":
        return None
    # "Nh extended darkness (Mh)" — use absolute time M
    if "extended darkness" in tp and "(" in tp:
        m = re.search(r"\((\d+(?:\.\d+)?)h\)", tp)
        return float(m.group(1)) if m else None
    # "day N" / "Day N"
    if tp.lower().startswith("day "):
        m = re.match(r"[Dd]ay\s+(\d+(?:\.\d+)?)", tp)
        return float(m.group(1)) * 24 if m else None
    # "Nh", "N h", "-Nh", "0.5h post-inoculation", "50h (P added)"
    m = re.match(r"(-?\d+(?:\.\d+)?)\s*h", tp)
    return float(m.group(1)) if m else None
