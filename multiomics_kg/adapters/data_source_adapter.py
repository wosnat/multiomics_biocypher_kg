"""DataSource adapter — emits 4 metadata nodes describing the data sources
that contribute fields to Gene records. Soft-joined to Gene via
Gene.contributing_sources (no materialized edges).

Auto-generates DataSource.info_types by walking the fields block of
gene_annotations_config.yaml and inverting source -> field-name mapping.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

import yaml


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV output (avoid ' and |)."""
    return value.replace("'", "^").replace("|", "")


class DataSourceAdapter:
    """Emit one DataSource node per logical source declared in
    gene_annotations_config.yaml.

    Required YAML keys per source:
    - logical_sources: list of dicts with id, scope, provenance,
      applies_to_organisms (optional)

    Fails loudly if logical_sources is missing rather than falling back.
    """

    def __init__(self, config_path: str | Path = "config/gene_annotations_config.yaml"):
        self.config_path = Path(config_path)
        self._config: dict[str, Any] | None = None
        self._info_types_by_source: dict[str, list[str]] = {}

    def download_data(self, cache: bool = True) -> None:
        with open(self.config_path) as f:
            self._config = yaml.safe_load(f)
        self._info_types_by_source = self._derive_info_types(self._config)

    @staticmethod
    def _derive_info_types(config: dict[str, Any]) -> dict[str, list[str]]:
        """Walk fields:; for each field rule, collect contributing source(s).

        When a rule has source_label, that label is the logical source.
        When a rule has only source, the source name is used (e.g. 'gene_mapping',
        'eggnog', 'uniprot'). The get_nodes() method then unions the yaml_source_name
        bucket and the logical src_id bucket together for each node.
        """
        info_by_src: dict[str, set[str]] = defaultdict(set)
        for field_name, rule in (config.get("fields") or {}).items():
            # Skip computed fields that have no source
            if isinstance(rule, dict) and rule.get("type") == "computed":
                continue
            for src in DataSourceAdapter._sources_from_rule(rule):
                info_by_src[src].add(field_name)
        return {k: sorted(v) for k, v in info_by_src.items()}

    @staticmethod
    def _sources_from_rule(rule: Any) -> Iterator[str]:
        """A field rule may declare a single source, a list of candidates, or
        nested options. Walk recursively; surface every source / source_label
        encountered.

        Priority: when source_label is present at the same level as source,
        yield source_label (the logical name) rather than the raw source.
        When only source is present (no source_label at this level), yield source.
        """
        if isinstance(rule, dict):
            # Check if this dict has both source and source_label at the same level
            has_source = "source" in rule
            has_source_label = "source_label" in rule

            if has_source_label:
                yield rule["source_label"]
                # Still recurse into other values except source/source_label
                for k, v in rule.items():
                    if k not in ("source", "source_label"):
                        yield from DataSourceAdapter._sources_from_rule(v)
            elif has_source:
                yield rule["source"]
                # Still recurse into other values except source
                for k, v in rule.items():
                    if k != "source":
                        yield from DataSourceAdapter._sources_from_rule(v)
            else:
                # No source at this level — recurse into all values
                for v in rule.values():
                    yield from DataSourceAdapter._sources_from_rule(v)
        elif isinstance(rule, list):
            for item in rule:
                yield from DataSourceAdapter._sources_from_rule(item)

    def _all_logical_sources(self) -> Iterator[dict[str, Any]]:
        """Iterate all logical sources from all source entries."""
        sources_block = (self._config or {}).get("sources") or {}
        for src_name, src_def in sources_block.items():
            if "logical_sources" not in src_def:
                raise ValueError(
                    f"Source '{src_name}' in gene_annotations_config.yaml "
                    f"missing required 'logical_sources' key."
                )
            for ls in src_def["logical_sources"]:
                yield {"_yaml_source_name": src_name, **ls}

    @staticmethod
    def _name_for(source_id: str) -> str:
        return {
            "ncbi": "NCBI RefSeq",
            "cyanorak": "Cyanorak",
            "uniprot": "UniProt",
            "eggnog": "EggNOG-mapper",
            "psortb": "PSORTb",
        }.get(source_id, source_id.title())

    @staticmethod
    def _description_for(source_id: str) -> str:
        return {
            "ncbi": "NCBI RefSeq genome annotation (GFF + GenBank).",
            "cyanorak": "Cyanorak: curated cyanobacteria gene clusters and roles.",
            "uniprot": "UniProt proteins (reviewed + unreviewed) cross-referenced via RefSeq protein_id.",
            "eggnog": "EggNOG-mapper functional annotations (per-protein query against the eggNOG reference).",
            "psortb": "PSORTb v3.0.3 per-protein subcellular-localization predictions (Gram-negative model).",
        }.get(source_id, "")

    def get_nodes(self) -> Iterator[tuple[str, str, dict[str, Any]]]:
        if self._config is None:
            self.download_data()
        for ls in self._all_logical_sources():
            src_id = ls["id"]
            yaml_name = ls["_yaml_source_name"]
            # info_types: union of fields attributed to this YAML source name
            # (e.g. 'gene_mapping' for ncbi/cyanorak) and fields attributed
            # directly to the logical source id (e.g. 'ncbi', 'cyanorak', 'eggnog').
            yaml_bucket = set(self._info_types_by_source.get(yaml_name, []))
            id_bucket = set(self._info_types_by_source.get(src_id, []))
            info_types = sorted(yaml_bucket | id_bucket)
            node_id = f"data_source:{src_id}"
            properties = {
                "id": _clean_str(src_id),
                "name": _clean_str(self._name_for(src_id)),
                "description": _clean_str(self._description_for(src_id)),
                "version": "",   # populated case-by-case in v2
                "scope": _clean_str(ls.get("scope", "")),
                "provenance": _clean_str(ls.get("provenance", "")),
                "applies_to_organisms": [_clean_str(o) for o in (ls.get("applies_to_organisms") or [])],
                "info_types": [_clean_str(t) for t in info_types],
            }
            yield node_id, "data_source", properties

    def get_edges(self) -> Iterator[tuple]:
        return iter(())
