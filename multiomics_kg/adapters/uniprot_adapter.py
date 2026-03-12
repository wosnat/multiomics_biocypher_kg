"""Simplified UniProt adapter.

Reads pre-built protein_annotations.json files
(output of multiomics_kg/download/build_protein_annotations.py)
and yields protein nodes and edges.  No data downloading or preprocessing
occurs at KG build time — all preprocessing is done by prepare_data.sh.

Nodes:   Protein  (uniprot:<accession>)
Edges:
  Gene_encodes_protein            (Gene → Protein, via RefSeq WP_ join)
  Protein_belongs_to_organism     (Protein → OrganismTaxon, per assembly)
  protein_catalyzes_ec_number     (Protein → EC)
  protein_located_in_cellular_component    (Protein → GO)
  protein_involved_in_biological_process   (Protein → GO)
  protein_contributes_to_molecular_function (Protein → GO)
"""
from __future__ import annotations

import csv
import json
import os
from collections import OrderedDict
from collections.abc import Generator
from pathlib import Path
from typing import Optional

from bioregistry import normalize_curie
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent


class UniprotAdapter:
    """Single-taxid adapter: reads protein_annotations.json, yields protein nodes + edges."""

    def __init__(
        self,
        organism_group: str,
        ncbi_taxon_id: int,
        assembly_info: list[dict],
        data_dirs: list[str],
        test_mode: bool = False,
    ):
        """
        Args:
            organism_group: e.g. "Prochlorococcus" or "Alteromonas"
            ncbi_taxon_id: NCBI taxonomy ID
            assembly_info: list of dicts with keys 'accession', 'strain_name',
                           'ncbi_taxon_id'; parallel to data_dirs
            data_dirs: list of data_dir paths (one per assembly), used to find
                       gene_mapping.csv for Gene→Protein edge join
            test_mode: if True, limits output to first 100 items per generator
        """
        self.organism_group = organism_group
        self.ncbi_taxon_id = ncbi_taxon_id
        self.assembly_info = assembly_info
        self.data_dirs = data_dirs
        self.test_mode = test_mode

        self.data_path = (
            PROJECT_ROOT / "cache" / "data" / organism_group
            / "uniprot" / str(ncbi_taxon_id) / "protein_annotations.json"
        )

        # provenance
        self.data_source = "uniprot"
        self.data_licence = "CC BY 4.0"
        self.data_version = "2024_03"

        self._data: dict[str, dict] = {}
        # {refseq_WP_id: [(locus_tag, ncbi_accession), ...]}
        self._refseq_to_strains: dict[str, list[tuple[str, str]]] = {}

    def _add_prefix(self, prefix: str, identifier: str) -> str:
        return normalize_curie(f"{prefix}:{identifier}")

    def _load_gene_mapping(self) -> dict[str, list[tuple[str, str]]]:
        """Build {RefSeq WP_ → [(locus_tag, ncbi_accession)]} from gene_mapping.csv files.

        One WP_ accession may appear in multiple strains (shared proteins).
        """
        refseq_to_strains: dict[str, list[tuple[str, str]]] = {}
        for data_dir, info in zip(self.data_dirs, self.assembly_info):
            gene_mapping_path = os.path.join(data_dir, "gene_mapping.csv")
            if not os.path.exists(gene_mapping_path):
                logger.warning(
                    f"[UniprotAdapter] gene_mapping.csv not found: {gene_mapping_path}"
                )
                continue
            ncbi_acc = info["accession"]
            with open(gene_mapping_path, newline="") as f:
                for row in csv.DictReader(f):
                    protein_id = (row.get("protein_id") or "").strip()
                    locus_tag = (row.get("locus_tag") or "").strip()
                    if protein_id and locus_tag:
                        refseq_to_strains.setdefault(protein_id, []).append(
                            (locus_tag, ncbi_acc)
                        )
        return refseq_to_strains

    def download_data(self, cache: bool = True, **kwargs) -> None:
        """Load protein_annotations.json and build refseq→strain mapping."""
        if not self.data_path.exists():
            logger.warning(
                f"[UniprotAdapter] protein_annotations.json not found: {self.data_path}\n"
                "  Run: uv run python multiomics_kg/download/build_protein_annotations.py --force"
            )
            return

        with open(self.data_path) as f:
            self._data = json.load(f)

        logger.info(
            f"[UniprotAdapter] Loaded {len(self._data)} proteins "
            f"for taxid {self.ncbi_taxon_id} ({self.organism_group})"
        )

        self._refseq_to_strains = self._load_gene_mapping()

        n_total = len(self._data)
        n_mapped = sum(
            1 for entry in self._data.values()
            if any(
                rs in self._refseq_to_strains
                for rs in (entry.get("refseq_ids") or [])
            )
        )
        logger.info(
            f"[UniprotAdapter] {n_mapped}/{n_total} proteins matched "
            f"to gene_mapping across {len(self.data_dirs)} strain(s)"
        )

    def _provenance(self) -> dict:
        return {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

    def get_nodes(self) -> Generator[tuple[str, str, dict]]:
        """Yield protein nodes (uniprot:<accession>)."""
        count = 0
        for uid, entry in self._data.items():
            if self.test_mode and count >= 100:
                break
            protein_id = self._add_prefix("uniprot", uid)
            props = {
                "gene_symbol":              entry.get("gene_symbol"),
                "protein_synonyms":         entry.get("protein_synonyms"),
                "locus_tag":                entry.get("locus_tag"),
                "gene_names":               entry.get("gene_names"),
                "sequence_length":          entry.get("sequence_length"),
                "molecular_mass":           entry.get("molecular_mass"),
                "refseq_ids":               entry.get("refseq_ids"),
                "proteome_ids":             entry.get("proteome_ids"),
                "ec_numbers":               entry.get("ec_numbers"),
                "go_cellular_components":   entry.get("go_cellular_components"),
                "go_biological_processes":  entry.get("go_biological_processes"),
                "go_molecular_functions":   entry.get("go_molecular_functions"),
                "function_description":     entry.get("function_description"),
                "catalytic_activities":     entry.get("catalytic_activities"),
                "cofactor_names":           entry.get("cofactor_names"),
                "pathways":                 entry.get("pathways"),
                "transmembrane_regions":    entry.get("transmembrane_regions"),
                "signal_peptide":           entry.get("signal_peptide"),
                "functional_motifs":        entry.get("functional_motifs"),
                "domain_description":       entry.get("domain_description"),
                "protein_family":           entry.get("protein_family"),
                "string_ids":               entry.get("string_ids"),
                "eggnog_ids":               entry.get("eggnog_ids"),
                "pfam_ids":                 entry.get("pfam_ids"),
                "subcellular_location":     entry.get("subcellular_location"),
                "organism_name":            entry.get("organism_name"),
                "annotation_score":         entry.get("annotation_score"),
                "is_reviewed":              entry.get("is_reviewed"),
                "keywords":                 entry.get("keywords"),
                "keyword_ids":              entry.get("keyword_ids"),
                "caution_notes":            entry.get("caution_notes"),
                "interaction_notes":        entry.get("interaction_notes"),
                **self._provenance(),
            }
            # Sparse output: drop None values
            props = {k: v for k, v in props.items() if v is not None}
            yield protein_id, "protein", props
            count += 1

    def get_edges(self) -> Generator[tuple[None, str, str, str, dict]]:
        """Yield all protein-related edges."""
        props = self._provenance()
        count = 0

        for uid, entry in self._data.items():
            if self.test_mode and count >= 100:
                break

            protein_id = self._add_prefix("uniprot", uid)

            refseq_ids = entry.get("refseq_ids") or []
            if isinstance(refseq_ids, str):
                refseq_ids = [refseq_ids]

            # Gene_encodes_protein + Protein_belongs_to_organism
            # (join via refseq_to_strains: one edge per matched locus_tag/assembly)
            seen_edges: set[tuple[str, str]] = set()
            for refseq in refseq_ids:
                for locus_tag, ncbi_acc in self._refseq_to_strains.get(refseq, []):
                    key = (locus_tag, ncbi_acc)
                    if key in seen_edges:
                        continue
                    seen_edges.add(key)
                    gene_id = self._add_prefix("ncbigene", locus_tag)
                    org_id = self._add_prefix("insdc.gcf", ncbi_acc)
                    # Direction convention: Gene → Protein (source=gene, target=protein)
                    yield None, gene_id, protein_id, "Gene_encodes_protein", props
                    yield None, protein_id, org_id, "Protein_belongs_to_organism", props

            # remove - moved to the gene edges
            # Protein → EC
            # for ec in (entry.get("ec_numbers") or []):
            #     if ec:
            #         yield None, protein_id, self._add_prefix("eccode", ec), \
            #             "protein_catalyzes_ec_number", props

            # Protein → GO (split by namespace)
            # for go in (entry.get("go_cellular_components") or []):
            #     if go:
            #         yield None, protein_id, go.lower(), \
            #             "protein_located_in_cellular_component", props
            # for go in (entry.get("go_biological_processes") or []):
            #     if go:
            #         yield None, protein_id, go.lower(), \
            #             "protein_involved_in_biological_process", props
            # for go in (entry.get("go_molecular_functions") or []):
            #     if go:
            #         yield None, protein_id, go.lower(), \
            #             "protein_contributes_to_molecular_function", props

            count += 1


class MultiUniprot:
    """Multi-taxid wrapper: creates one UniprotAdapter per unique taxid."""

    def __init__(self, config_list_file: str, test_mode: bool = False, **kwargs):
        """
        Args:
            config_list_file: Path to cyanobacteria_genomes.csv with columns:
                ncbi_accession, ncbi_taxon_id, strain_name, data_dir, ...
                Lines starting with # are treated as comments and skipped.
                One UniprotAdapter is created per unique ncbi_taxon_id.
        """
        self.adapters: list[UniprotAdapter] = []

        # Maintain insertion order; group by taxid
        taxid_to_org_group: dict[int, str] = OrderedDict()
        taxid_to_assemblies: dict[int, list[dict]] = OrderedDict()
        taxid_to_data_dirs: dict[int, list[str]] = OrderedDict()

        with open(config_list_file, newline="") as f:
            lines = [line for line in f if not line.strip().startswith("#")]
            reader = csv.DictReader(lines)
            for row in reader:
                taxon_id_str = (row.get("ncbi_taxon_id") or "").strip()
                if not taxon_id_str:
                    continue
                organism_id = int(taxon_id_str)
                data_dir = (row.get("data_dir") or "").strip()

                # Infer organism_group from data_dir path
                org_group = "Prochlorococcus"
                parts = Path(data_dir).parts
                for i, part in enumerate(parts):
                    if part == "data" and i + 1 < len(parts):
                        org_group = parts[i + 1]
                        break

                taxid_to_org_group[organism_id] = org_group
                info = {
                    "accession": row["ncbi_accession"],
                    "strain_name": row.get("strain_name") or "",
                    "ncbi_taxon_id": organism_id,
                }
                taxid_to_assemblies.setdefault(organism_id, []).append(info)
                taxid_to_data_dirs.setdefault(organism_id, []).append(data_dir)

        for organism_id, assembly_list in taxid_to_assemblies.items():
            org_group = taxid_to_org_group[organism_id]
            adapter = UniprotAdapter(
                organism_group=org_group,
                ncbi_taxon_id=organism_id,
                assembly_info=assembly_list,
                data_dirs=taxid_to_data_dirs.get(organism_id, []),
                test_mode=test_mode,
            )
            self.adapters.append(adapter)

        logger.info(
            f"[MultiUniprot] {len(self.adapters)} unique taxid(s) "
            f"({sum(len(v) for v in taxid_to_assemblies.values())} assemblies) "
            f"from {config_list_file}"
        )

    @property
    def organism_ids(self) -> list[int]:
        return [a.ncbi_taxon_id for a in self.adapters]

    def download_data(self, cache: bool = True, **kwargs) -> None:
        for adapter in self.adapters:
            adapter.download_data(cache=cache)

    def get_nodes(self) -> Generator[tuple[str, str, dict]]:
        for adapter in self.adapters:
            yield from adapter.get_nodes()

    def get_edges(self) -> Generator[tuple[None, str, str, str, dict]]:
        for adapter in self.adapters:
            yield from adapter.get_edges()
