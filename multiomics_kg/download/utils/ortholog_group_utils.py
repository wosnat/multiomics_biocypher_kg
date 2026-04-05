"""Ortholog group extraction utilities.

Extracts ortholog group memberships from gene annotation data
(cluster_number from Cyanorak, eggnog_ogs from EggNOG-mapper).

Called by build_gene_annotations.py post-merge to populate the
`ortholog_groups` field in gene_annotations_merged.json.
"""

from __future__ import annotations

# Whitelist: organism group → list of (target_taxon_id, specificity_rank).
# Each gene extracts ALL levels in the list independently.
# Bacteria-level COGs (taxon_id=2, rank=3) are always added separately.
ORGANISM_GROUP_LEVELS: dict[str, list[tuple[int, int]]] = {
    "Prochlorococcus":     [(1212, 1), (1117, 2)],    # Prochloraceae, Cyanobacteria
    "Synechococcus":       [(1129, 1), (1117, 2)],    # Synechococcus, Cyanobacteria
    "Thermosynechococcus": [(1117, 2)],                # Cyanobacteria
    "Alteromonas":         [(72275, 1), (1224, 2)],    # Alteromonadaceae, Proteobacteria
    "Shewanella":          [(1224, 2)],                 # Proteobacteria
    "Pseudomonas":         [(1224, 2)],                 # Proteobacteria
    "Ruegeria":            [(1224, 2)],                 # Proteobacteria
    "Meiothermus":         [],                          # Bacteria-level only
}


def organism_group_from_path(data_dir: str) -> str:
    """Derive organism group from data_dir path.

    E.g. 'cache/data/Prochlorococcus/genomes/MED4/' → 'Prochlorococcus'
    Note: WH8102 (Parasynechococcus) is under Synechococcus/ in data_dir.
    """
    for group in ORGANISM_GROUP_LEVELS:
        if group in data_dir:
            return group
    return "unknown"


def _parse_eggnog_ogs(eggnog_ogs: list[str]) -> dict[int, tuple[str, str]]:
    """Parse eggnog_ogs entries into {taxon_id: (og_id, level_name)}.

    Format in gene_annotations_merged.json is "OG_ID@taxon_id|level_name"
    (pipe-separated — clean_text is NOT applied to annotation JSON).
    Entries without '@' are legacy short names and are skipped.
    """
    parsed: dict[int, tuple[str, str]] = {}
    for entry in eggnog_ogs:
        if "@" not in entry:
            continue
        og_part, rest = entry.split("@", 1)
        parts = rest.split("|", 1)
        if len(parts) != 2:
            continue
        try:
            taxon_id = int(parts[0])
        except ValueError:
            continue
        parsed[taxon_id] = (og_part, parts[1])
    return parsed


def extract_ortholog_groups(gene: dict, organism_group: str) -> list[dict]:
    """Extract ortholog group memberships from a gene's annotations.

    Args:
        gene: Gene dict from gene_annotations_merged.json.
        organism_group: One of "Prochlorococcus", "Synechococcus", "Alteromonas"
            (derived from data_dir path). Used to select the correct target
            taxonomic level for lowest-level OG.

    Returns:
        List of {og_id, source, taxonomic_level, taxon_id} dicts. Deduplicated
        by og_id (a gene cannot belong to the same OG twice).
    """
    groups: list[dict] = []
    seen_ids: set[str] = set()

    # 1. Cyanorak cluster (Pro/Syn only — Alt genes won't have cluster_number)
    cluster = gene.get("cluster_number")
    if cluster:
        og_id = f"cyanorak:{cluster}"
        groups.append({
            "og_id": og_id,
            "source": "cyanorak",
            "taxonomic_level": "curated",
            "taxon_id": 0,
            "specificity_rank": 0,
        })
        seen_ids.add(og_id)

    # 2. Parse eggnog_ogs for bacteria-level and lowest-level
    eggnog_ogs = gene.get("eggnog_ogs") or []
    parsed = _parse_eggnog_ogs(eggnog_ogs)

    # 2a. Bacteria-level COG
    if 2 in parsed:
        og_part, level_name = parsed[2]
        og_id = f"eggnog:{og_part}@2"
        if og_id not in seen_ids:
            groups.append({
                "og_id": og_id,
                "source": "eggnog",
                "taxonomic_level": level_name,
                "taxon_id": 2,
                "specificity_rank": 3,
            })
            seen_ids.add(og_id)

    # 2b. Configured intermediate levels — extract ALL independently
    levels = ORGANISM_GROUP_LEVELS.get(organism_group, [])
    for tid, rank in levels:
        if tid in parsed:
            og_part, level_name = parsed[tid]
            og_id = f"eggnog:{og_part}@{tid}"
            if og_id not in seen_ids:
                groups.append({
                    "og_id": og_id,
                    "source": "eggnog",
                    "taxonomic_level": level_name,
                    "taxon_id": tid,
                    "specificity_rank": rank,
                })
                seen_ids.add(og_id)

    return groups
