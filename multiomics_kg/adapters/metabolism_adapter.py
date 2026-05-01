"""Metabolism adapter — emits Reaction + Metabolite nodes and 4 edge types.

Pure file reader: consumes
- cache/data/kegg/kegg_data.json (built by build_kegg_metabolism_xrefs, step 6)
- per-strain cache/data/<organism>/genomes/<strain>/gene_annotations_merged.json

No SQLite / no KEGG REST at build time.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterator

import chemparse

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

log = logging.getLogger(__name__)


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_KEGG_DATA = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_data.json"


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV: replace single quotes → ^, strip pipes."""
    if not isinstance(value, str):
        return value
    return value.replace("'", "^").replace("|", "")


def _drop_nulls(props: dict) -> dict:
    """Sparse-output convention: drop keys whose value is None or empty list."""
    return {k: v for k, v in props.items() if v is not None and v != []}


def _parse_elements(formula: str | None) -> list[str]:
    """Sorted unique element symbols present in a Hill-notation formula.

    Empty list when formula is null/empty. Returns empty list (rather than
    raising) on any parse failure, so a malformed KEGG formula cannot break
    the build.
    """
    if not formula:
        return []
    try:
        return sorted(chemparse.parse_formula(formula).keys())
    except Exception:
        return []


class MetabolismAdapter:
    """Single-pass adapter emitting Reaction + Metabolite nodes from the kegg_data cache."""

    def __init__(self, kegg_data_path: Path | str | None = None, test_mode: bool = False):
        self.kegg_data_path = Path(kegg_data_path) if kegg_data_path else DEFAULT_KEGG_DATA
        self.test_mode = test_mode
        self._kegg_data: dict | None = None

    def _load(self) -> dict:
        """Load the pruned kegg_data.json. Strict mode: raises if missing or corrupt."""
        if self._kegg_data is None:
            if not self.kegg_data_path.exists():
                raise FileNotFoundError(
                    f"{self.kegg_data_path} missing — run "
                    f"`bash scripts/prepare_data.sh --steps 6` first."
                )
            try:
                self._kegg_data = json.loads(self.kegg_data_path.read_text())
            except json.JSONDecodeError as e:
                raise RuntimeError(
                    f"{self.kegg_data_path} is corrupted ({e}). Rebuild via "
                    f"`bash scripts/prepare_data.sh --steps 6 --force`."
                ) from e
        return self._kegg_data

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        data = self._load()
        n = 0
        for rxn_id, rxn in data.get("reactions", {}).items():
            if self.test_mode and n >= 100:
                break
            node_id = f"kegg.reaction:{rxn_id}"
            props = _drop_nulls({
                "kegg_reaction_id": rxn_id,
                "name": _clean_str(rxn.get("name", "")),
                "ec_numbers": list(rxn.get("ec_numbers", [])),
                "kegg_pathway_ids": list(rxn.get("pathways", [])),
                "mnxr_id": rxn.get("mnxr_id"),
                "rhea_ids": list(rxn.get("rhea_ids", [])),
                "mass_balance": rxn.get("mass_balance"),
                "reaction_class": rxn.get("reaction_class"),
            })
            yield node_id, "reaction", props
            n += 1

        n = 0
        for cpd_id, cpd in data.get("compounds", {}).items():
            if self.test_mode and n >= 100:
                break
            node_id = f"kegg.compound:{cpd_id}"
            props = _drop_nulls({
                "kegg_compound_id": cpd_id,
                "name": _clean_str(cpd.get("name", "")),
                "formula": _clean_str(cpd.get("formula")),
                "mass": cpd.get("mass"),
                "inchikey": _clean_str(cpd.get("inchikey")),
                "smiles": _clean_str(cpd.get("smiles")),
                "mnxm_id": cpd.get("mnxm_id"),
                "chebi_id": cpd.get("chebi_id"),
                "hmdb_id": cpd.get("hmdb_id"),
            })
            yield node_id, "metabolite", props
            n += 1

    def _reaction_metabolite_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        for rxn_id, rxn in data.get("reactions", {}).items():
            for cpd_id in rxn.get("compounds", []):
                edge_id = f"r2m:{rxn_id}:{cpd_id}"
                yield (
                    edge_id,
                    f"kegg.reaction:{rxn_id}",
                    f"kegg.compound:{cpd_id}",
                    "reaction_has_metabolite",
                    {},
                )

    def _reaction_pathway_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        for rxn_id, rxn in data.get("reactions", {}).items():
            for pw_id in rxn.get("pathways", []):
                edge_id = f"r2p:{rxn_id}:{pw_id}"
                yield (
                    edge_id,
                    f"kegg.reaction:{rxn_id}",
                    f"kegg.pathway:{pw_id}",
                    "reaction_in_kegg_pathway",
                    {},
                )

    def _metabolite_pathway_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        for cpd_id, cpd in data.get("compounds", {}).items():
            for pw_id in cpd.get("pathways", []):
                edge_id = f"m2p:{cpd_id}:{pw_id}"
                yield (
                    edge_id,
                    f"kegg.compound:{cpd_id}",
                    f"kegg.pathway:{pw_id}",
                    "metabolite_in_pathway",
                    {},
                )

    def get_edges(self) -> Iterator[tuple[str | None, str, str, str, dict]]:
        yield from self._reaction_metabolite_edges()
        yield from self._reaction_pathway_edges()
        yield from self._metabolite_pathway_edges()


class MultiMetabolismAdapter(MetabolismAdapter):
    """Multi-strain wrapper. Nodes come from the global kegg_data cache; edges
    come from per-strain gene_annotations_merged.json (see Task 11)."""

    def __init__(self, genome_config_file: str, kegg_data_path: Path | str | None = None,
                 test_mode: bool = False):
        super().__init__(kegg_data_path=kegg_data_path, test_mode=test_mode)
        self.genome_config_file = genome_config_file

    def _gene_reaction_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
        data = self._load()
        known_rxns = set(data.get("reactions", {}).keys())
        for row in load_genome_rows():
            genes = load_gene_annotations(row["data_dir"])
            if not genes:
                continue
            for locus_tag, gene in genes.items():
                for rxn_id in gene.get("kegg_reactions", []) or []:
                    if rxn_id not in known_rxns:
                        continue  # skip dangling: keeps KG free of orphan edges
                    edge_id = f"g2r:{locus_tag}:{rxn_id}"
                    yield (
                        edge_id,
                        f"ncbigene:{locus_tag}",
                        f"kegg.reaction:{rxn_id}",
                        "gene_catalyzes_reaction",
                        {},
                    )

    def get_edges(self) -> Iterator[tuple[str | None, str, str, str, dict]]:
        yield from self._gene_reaction_edges()
        yield from self._reaction_metabolite_edges()
        yield from self._reaction_pathway_edges()
        yield from self._metabolite_pathway_edges()
