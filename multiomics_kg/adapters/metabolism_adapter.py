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

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.metabolite_identity.names import parse_elements
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
    """Sorted unique element symbols (allow-list parser; see metabolite_identity.names)."""
    return parse_elements(formula)


_MET_STR = ("name", "formula", "inchikey", "smiles", "mnxm_id", "chebi_id", "hmdb_id", "structure_kind")
_MET_LIST = ("kegg_compound_ids", "chebi_ids", "mnxm_ids", "all_identifiers", "recurated_from_ids",
             "identity_basis", "evidence_sources")


def _metabolite_props(cpd: dict, kegg_id: str | None) -> dict:
    props = {"kegg_compound_id": kegg_id, "mass": cpd.get("mass"),
             "elements": cpd.get("elements") or _parse_elements(cpd.get("formula"))}
    props.update({k: _clean_str(cpd.get(k)) for k in _MET_STR})
    props.update({k: [_clean_str(v) for v in (cpd.get(k) or [])] for k in _MET_LIST})
    return _drop_nulls(props)


class MetabolismAdapter:
    """Single-pass adapter emitting Reaction + Metabolite nodes from the kegg_data cache."""

    def __init__(self, kegg_data_path: Path | str | None = None, test_mode: bool = False):
        self.kegg_data_path = Path(kegg_data_path) if kegg_data_path else DEFAULT_KEGG_DATA
        self.test_mode = test_mode
        self._kegg_data: dict | None = None
        self._emitted_metabolite_ids: set[str] = set()

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
            props = _metabolite_props(cpd, cpd_id)
            self._emitted_metabolite_ids.add(node_id)
            yield node_id, "metabolite", props
            n += 1

        n = 0
        for primary_id, cpd in data.get("additional_compounds", {}).items():
            if self.test_mode and n >= 100:
                break
            # primary_id is already in canonical form (e.g. 'chebi:9999', 'mnx:MNXM12345')
            props = _metabolite_props(cpd, None)
            self._emitted_metabolite_ids.add(primary_id)
            yield primary_id, "metabolite", props
            n += 1

    def all_metabolite_ids(self) -> set[str]:
        """Ids emitted by get_nodes (or every compound id if get_nodes hasn't run)."""
        if self._emitted_metabolite_ids:
            return set(self._emitted_metabolite_ids)
        data = self._load()
        return ({f"kegg.compound:{k}" for k in data["compounds"]}
                | set(data.get("additional_compounds", {})))

    def _metabolite_relation_edges(self):
        data = self._load()
        known = self.all_metabolite_ids()
        entries = [(f"kegg.compound:{k}", v) for k, v in data["compounds"].items()]
        entries += list(data.get("additional_compounds", {}).items())
        skipped = 0
        for src, cpd in entries:
            if src not in known:
                continue
            for tgt in cpd.get("stereo_form_of") or []:
                if tgt not in known:
                    skipped += 1
                    continue
                yield (f"m2s:{src}:{tgt}", src, tgt, "metabolite_is_stereo_form_of",
                       {"sources": ["chebi"], "evidence": "curated"})
            for tgt in cpd.get("enantiomer_of") or []:
                if tgt not in known:
                    skipped += 1
                    continue
                assert src < tgt, f"enantiomer pair stored high->low: {src} {tgt}"
                yield (f"m2e:{src}:{tgt}", src, tgt, "metabolite_is_enantiomer_of",
                       {"sources": ["chebi"], "evidence": "curated"})
        if skipped:
            log.info(f"metabolism: skipped {skipped} stereo/enantiomer edges to unemitted metabolites")

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
        yield from self._metabolite_relation_edges()


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
