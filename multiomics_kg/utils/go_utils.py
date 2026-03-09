"""
Lightweight Gene Ontology utilities.

Instead of downloading the full GO ontology via pypath (slow, all 30K terms via QuickGO API),
this module downloads go-basic.obo once (~7 MB), parses it into a compact namespace cache
(~3-5 MB JSON), and reuses that cache on all subsequent runs.

Cache layout:
  cache_root/go_terms/go_namespace_cache.json  — compact {go_id: {name, namespace, parents}}
  go-basic.obo                                  — downloaded temporarily, deleted after parsing
"""

import json
import logging
import re
import requests
from pathlib import Path

logger = logging.getLogger(__name__)

GO_OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

# Maps OBO namespace string → BioCypher schema node label
NAMESPACE_TO_LABEL: dict[str, str] = {
    "biological_process": "biological process",
    "molecular_function": "molecular function",
    "cellular_component": "cellular component",
}

# Allowed GO-GO edge labels (must match schema_config.yaml label_in_input values)
_ALLOWED_GO_GO_EDGE_LABELS: frozenset[str] = frozenset({
    "biological_process_is_a_biological_process",
    "biological_process_part_of_biological_process",
    "biological_process_positively_regulates_biological_process",
    "biological_process_negatively_regulates_biological_process",
    "biological_process_positively_regulates_molecular_function",
    "biological_process_negatively_regulates_molecular_function",
    "molecular_function_is_a_molecular_function",
    "molecular_function_part_of_molecular_function",
    "molecular_function_positively_regulates_molecular_function",
    "molecular_function_negatively_regulates_molecular_function",
    "cellular_component_is_a_cellular_component",
    "cellular_component_part_of_cellular_component",
})

# Regex to extract the GO ID from an is_a or relationship line
_GO_ID_RE = re.compile(r"(GO:\d+)")


def _go_cache_dir(cache_root: Path) -> Path:
    d = cache_root / "go_terms"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_obo(obo_path: Path) -> None:
    """Download go-basic.obo to obo_path (streaming, with progress log)."""
    logger.info(f"Downloading go-basic.obo from {GO_OBO_URL} …")
    resp = requests.get(GO_OBO_URL, stream=True, timeout=120)
    resp.raise_for_status()
    with open(obo_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    logger.info(f"Downloaded go-basic.obo ({obo_path.stat().st_size // 1024} KB)")


def parse_obo(obo_path: Path) -> dict[str, dict]:
    """
    Parse a go-basic.obo file and return::

        {
          "GO:0005737": {
            "name": "cytoplasm",
            "namespace": "cellular_component",
            "parents": [["GO:0005622", "is_a"], ...]
          },
          ...
        }

    Only non-obsolete [Term] stanzas are included.
    Parents are collected from ``is_a:`` and ``relationship: part_of|regulates|…`` lines.
    """
    result: dict[str, dict] = {}
    current: dict = {}
    in_term = False

    with open(obo_path, encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            if line == "[Term]":
                in_term = True
                current = {"parents": []}

            elif line.startswith("[") and line != "[Term]":
                # [Typedef] or other stanza — commit current term
                _commit_term(current, in_term, result)
                in_term = False
                current = {}

            elif line == "" and in_term:
                _commit_term(current, in_term, result)
                in_term = False
                current = {}

            elif in_term:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line.startswith("is_obsolete: true"):
                    current["obsolete"] = True
                elif line.startswith("is_a: "):
                    m = _GO_ID_RE.search(line)
                    if m:
                        current["parents"].append([m.group(1), "is_a"])
                elif line.startswith("relationship: "):
                    # e.g. "relationship: part_of GO:0043231 ! ..."
                    parts = line[len("relationship: "):].split()
                    if len(parts) >= 2:
                        relation = parts[0]  # part_of / regulates / positively_regulates / …
                        m = _GO_ID_RE.search(parts[1])
                        if m:
                            current["parents"].append([m.group(1), relation])

    # Handle last stanza if file doesn't end with blank line
    _commit_term(current, in_term, result)

    logger.info(f"Parsed {len(result)} non-obsolete GO terms from {obo_path.name}")
    return result


def _commit_term(current: dict, in_term: bool, result: dict) -> None:
    if in_term and current.get("id") and not current.get("obsolete"):
        result[current["id"]] = {
            "name": current.get("name", ""),
            "namespace": current.get("namespace", ""),
            "parents": current.get("parents", []),
        }


def build_namespace_cache(cache_root: Path, force: bool = False) -> Path:
    """
    Build (or rebuild) the compact namespace cache JSON.

    Downloads go-basic.obo, parses it, writes
    ``cache_root/go_terms/go_namespace_cache.json``, deletes the OBO file,
    and returns the cache path.
    """
    cache_path = _go_cache_dir(cache_root) / "go_namespace_cache.json"
    if cache_path.exists() and not force:
        logger.info(f"GO namespace cache already exists at {cache_path}")
        return cache_path

    obo_path = _go_cache_dir(cache_root) / "go-basic.obo"
    try:
        _download_obo(obo_path)
        go_data = parse_obo(obo_path)
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(go_data, fh, separators=(",", ":"), sort_keys=True)
        logger.info(
            f"Wrote GO namespace cache: {len(go_data)} terms, "
            f"{cache_path.stat().st_size // 1024} KB → {cache_path}"
        )
    finally:
        if obo_path.exists():
            obo_path.unlink()
            logger.debug("Deleted go-basic.obo after parsing")

    return cache_path


def load_go_data(cache_root: Path, force: bool = False) -> dict[str, dict]:
    """
    Load ``{go_id: {name, namespace, parents}}`` from cache, building it first if needed.

    Args:
        cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
        force: rebuild cache even if it exists (re-downloads OBO)

    Returns:
        dict mapping GO term IDs (e.g. ``"GO:0005737"``) to
        ``{"name": ..., "namespace": ..., "parents": [[parent_id, relation], ...]}``
    """
    cache_path = build_namespace_cache(cache_root, force=force)
    with open(cache_path, encoding="utf-8") as fh:
        data = json.load(fh)
    logger.info(f"Loaded {len(data)} GO terms from namespace cache")
    return data


def compute_ancestry_closure(seed_go_ids: set[str], go_data: dict[str, dict]) -> set[str]:
    """
    BFS from *seed_go_ids* up through all ancestors, returning the full closure.

    The closure includes the seed terms themselves plus every ancestor reachable
    via any parent relationship stored in *go_data*.

    Args:
        seed_go_ids: directly referenced GO term IDs (e.g. from gene_annotations_merged.json)
        go_data: full GO data dict from :func:`load_go_data`

    Returns:
        Set of GO IDs (seed + all transitive ancestors that exist in *go_data*)
    """
    closure: set[str] = set()
    queue = list(seed_go_ids & go_data.keys())

    while queue:
        go_id = queue.pop()
        if go_id in closure:
            continue
        closure.add(go_id)
        entry = go_data.get(go_id)
        if entry is None:
            continue
        for parent_id, _relation in entry.get("parents", []):
            if parent_id not in closure and parent_id in go_data:
                queue.append(parent_id)

    return closure


def make_go_go_edge_label(
    child_ns: str,
    relation: str,
    parent_ns: str,
) -> str | None:
    """
    Build the GO-GO edge label string from namespace and relation, returning
    ``None`` if the resulting label is not in the allowed schema set.

    Example: ``("biological_process", "is_a", "biological_process")``
    → ``"biological_process_is_a_biological_process"``
    """
    label = f"{child_ns}_{relation}_{parent_ns}"
    return label if label in _ALLOWED_GO_GO_EDGE_LABELS else None
