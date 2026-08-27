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

from multiomics_kg.utils.curie_utils import normalize_curie
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent


def _clean_str(value: str) -> str:
    """Repo-wide string-property sanitiser (see CLAUDE.md)."""
    return value.replace("'", "^").replace("|", "")


def _clean_value(value):
    """Apply `_clean_str` to a property value, recursing into list elements.

    Defense in depth: `build_protein_annotations._sanitize` already scrubs both
    characters at the source, but every other adapter sanitises at the yield
    point too, and the one time this adapter did not, 239 pipe-bearing UniProt
    `catalytic_activities` values were split apart by BioCypher's `|` array
    delimiter on the way into the CSV.
    """
    if isinstance(value, str):
        return _clean_str(value)
    if isinstance(value, list):
        return [_clean_str(v) if isinstance(v, str) else v for v in value]
    return value


# A UniProt proteome accession is adopted as an assembly's "own" proteome when at
# least this many of the assembly's WP_-matched proteins carry it. Real support is
# in the hundreds-to-thousands; the threshold only has to reject strays (a UniProt
# entry listed under several proteomes, e.g. 73 A. macleodii proteins under three).
MIN_OWN_PROTEOME_SUPPORT = 20
# ...and those must account for at least this fraction of the proteome's entries.
# Measured: own proteomes 71-100%; a different annotation build of the same
# isolate ~2% (DSS-3's UP000813672); every proteome under species-level taxid
# 28108 <= 27%.
OWN_PROTEOME_MIN_FRACTION = 0.5

_LOCUS_TAG_COLUMNS = ("locus_tag", "locus_tag_ncbi", "old_locus_tags", "locus_tag_cyanorak")


def _proteome_accessions(entry: dict) -> set[str]:
    """'UP000001026: Chromosome' → 'UP000001026' for every proteome_ids element."""
    raw = entry.get("proteome_ids") or []
    if isinstance(raw, str):
        raw = [raw]
    return {str(p).split(":")[0].strip() for p in raw if p}


class UniprotAdapter:
    """Single-taxid adapter: reads protein_annotations.json, yields protein nodes + edges.

    Which proteins make it into the graph (plans/orphan_proteins.md): UniProt is
    queried per taxid, and a taxid can be species-level (28108 = every *A.
    macleodii* isolate UniProt knows) or backed by a different annotation build
    than our RefSeq download. A protein is therefore emitted only when it is
    linkable to one of OUR assemblies — by RefSeq WP_ join, by locus-tag join,
    or by belonging to the assembly's own UniProt proteome. Everything else is
    another organism's protein and is dropped rather than left as an orphan.
    """

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

        # Provenance. No `version`: nothing in the download path
        # (download_uniprot.py) captures a UniProt release, so any release string
        # here is a guess. The previous hardcoded "2024_03" was stamped on ~89K
        # edges while the underlying uniprot_raw_data.json had in fact been
        # fetched in 2026 — a false claim is worse than an absent field. To
        # restore it, capture X-UniProt-Release at download time and plumb it
        # through protein_annotations.json rather than reintroducing a constant.
        self.data_source = "uniprot"
        self.data_licence = "CC BY 4.0"

        self._data: dict[str, dict] = {}
        # {refseq_WP_id: [(locus_tag, ncbi_accession), ...]}
        self._refseq_to_strains: dict[str, list[tuple[str, str]]] = {}
        # {any locus-tag spelling (current/NCBI/old/cyanorak): [(locus_tag, ncbi_accession), ...]}
        # Fallback gene join for UniProt entries that carry no RefSeq xref but
        # name a locus tag via gene_oln / gene_names (~7.6K proteins corpus-wide).
        self._locus_to_strains: dict[str, list[tuple[str, str]]] = {}
        # {ncbi_accession: {UniProt proteome accession, ...}} — the proteome(s)
        # backing each of OUR assemblies, derived from its WP_-matched proteins.
        self._own_proteomes: dict[str, set[str]] = {}

    def _add_prefix(self, prefix: str, identifier: str) -> str:
        return normalize_curie(f"{prefix}:{identifier}")

    def _load_gene_mapping(
        self,
    ) -> tuple[dict[str, list[tuple[str, str]]], dict[str, list[tuple[str, str]]]]:
        """Build the two join indexes from each assembly's gene_mapping.csv.

        Returns ``(refseq_to_strains, locus_to_strains)``:
          - ``{RefSeq WP_ → [(locus_tag, ncbi_accession)]}`` — primary join.
            One WP_ accession may appear in multiple strains (shared proteins).
          - ``{locus-tag spelling → [(locus_tag, ncbi_accession)]}`` over the
            current, NCBI, old and Cyanorak locus-tag columns — fallback join
            for entries with no RefSeq xref.
        """
        refseq_to_strains: dict[str, list[tuple[str, str]]] = {}
        locus_to_strains: dict[str, list[tuple[str, str]]] = {}
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
                    if not locus_tag:
                        continue
                    pair = (locus_tag, ncbi_acc)
                    if protein_id:
                        refseq_to_strains.setdefault(protein_id, []).append(pair)
                    for col in _LOCUS_TAG_COLUMNS:
                        for tok in (row.get(col) or "").replace("|", ",").split(","):
                            tok = tok.strip()
                            if tok:
                                bucket = locus_to_strains.setdefault(tok, [])
                                if pair not in bucket:
                                    bucket.append(pair)
        return refseq_to_strains, locus_to_strains

    def _derive_own_proteomes(self) -> dict[str, set[str]]:
        """{ncbi_accession: {UniProt proteome accession}} backing each assembly.

        Derived rather than registered. A proteome is adopted for an assembly
        when that assembly's WP_-matched proteins account for a MAJORITY
        (>= OWN_PROTEOME_MIN_FRACTION) of the proteome's entries in this taxid's
        download, with at least MIN_OWN_PROTEOME_SUPPORT of them, AND no other
        assembly also accounts for a majority. The direction matters: measuring
        instead how many of OUR matched proteins list the proteome is fooled by
        conserved proteins (MIT1002's matched proteins list UP000063991 59% of
        the time, yet 80% of that proteome's WP_ are absent from MIT1002's
        FASTA). Measured: own proteomes sit at 71-100% coverage; a different
        annotation build of the same isolate (DSS-3's UP000813672) at 2%;
        every proteome under species-level taxid 28108 at <= 27%, so the
        proteome-only proteins there are unattributable and dropped.
        """
        proteome_total: dict[str, int] = {}
        covered: dict[str, dict[str, int]] = {}
        for entry in self._data.values():
            proteomes = _proteome_accessions(entry)
            for up in proteomes:
                proteome_total[up] = proteome_total.get(up, 0) + 1
            for acc in {acc for _, acc in self._refseq_pairs(entry)}:
                for up in proteomes:
                    covered.setdefault(acc, {})
                    covered[acc][up] = covered[acc].get(up, 0) + 1
        majority: dict[str, set[str]] = {
            acc: {
                up for up, n in ups.items()
                if n >= MIN_OWN_PROTEOME_SUPPORT
                and n / proteome_total[up] >= OWN_PROTEOME_MIN_FRACTION
            }
            for acc, ups in covered.items()
        }
        own: dict[str, set[str]] = {}
        for acc, ups in majority.items():
            exclusive = {
                up for up in ups
                if not any(up in other for a2, other in majority.items() if a2 != acc)
            }
            if exclusive:
                own[acc] = exclusive
        return own

    # ── per-entry resolution ─────────────────────────────────────────────────

    def _refseq_pairs(self, entry: dict) -> list[tuple[str, str]]:
        refseq_ids = entry.get("refseq_ids") or []
        if isinstance(refseq_ids, str):
            refseq_ids = [refseq_ids]
        return [pair for rs in refseq_ids for pair in self._refseq_to_strains.get(rs, [])]

    def _locus_pairs(self, entry: dict) -> list[tuple[str, str]]:
        tokens: list[str] = []
        if entry.get("locus_tag"):
            tokens.append(str(entry["locus_tag"]).strip())
        names = entry.get("gene_names") or []
        if isinstance(names, str):
            names = [names]
        tokens.extend(str(n).strip() for n in names if n)
        pairs: list[tuple[str, str]] = []
        for tok in tokens:
            for pair in self._locus_to_strains.get(tok, []):
                if pair not in pairs:
                    pairs.append(pair)
        return pairs

    def _resolve(self, entry: dict) -> tuple[list[tuple[str, str]], list[str]]:
        """→ (gene pairs [(locus_tag, ncbi_acc)], assemblies the protein belongs to).

        Gene pairs come from the RefSeq join, else the locus-tag join (fallback
        only — a stale gene_oln must not add a second gene edge next to a WP_
        hit). Assemblies = those from the gene pairs ∪ those whose own proteome
        the entry lists. Empty assemblies ⇒ not our protein ⇒ not emitted.
        """
        gene_pairs = self._refseq_pairs(entry) or self._locus_pairs(entry)
        assemblies: list[str] = []
        for _, acc in gene_pairs:
            if acc not in assemblies:
                assemblies.append(acc)
        proteomes = _proteome_accessions(entry)
        if proteomes:
            for acc, own in self._own_proteomes.items():
                if acc not in assemblies and own & proteomes:
                    assemblies.append(acc)
        return gene_pairs, assemblies

    def _is_kept(self, entry: dict) -> bool:
        return bool(self._resolve(entry)[1])

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

        self._refseq_to_strains, self._locus_to_strains = self._load_gene_mapping()
        self._own_proteomes = self._derive_own_proteomes()

        n_total = len(self._data)
        n_refseq = n_locus = n_proteome_only = n_dropped = 0
        for entry in self._data.values():
            if self._refseq_pairs(entry):
                n_refseq += 1
            elif self._locus_pairs(entry):
                n_locus += 1
            elif self._is_kept(entry):
                n_proteome_only += 1
            else:
                n_dropped += 1
        logger.info(
            f"[UniprotAdapter] taxid {self.ncbi_taxon_id}: {n_total} proteins — "
            f"{n_refseq} RefSeq-joined, {n_locus} locus-tag-joined, "
            f"{n_proteome_only} own-proteome only (organism edge, no gene), "
            f"{n_dropped} dropped (not linkable to any of our "
            f"{len(self.data_dirs)} assembly/ies); own proteomes: "
            + ", ".join(f"{a}={sorted(u)}" for a, u in self._own_proteomes.items())
        )

    def _provenance(self) -> dict:
        return {
            "source": self.data_source,
            "licence": self.data_licence,
        }

    def get_nodes(self) -> Generator[tuple[str, str, dict]]:
        """Yield protein nodes (uniprot:<accession>)."""
        count = 0
        for uid, entry in self._data.items():
            if self.test_mode and count >= 100:
                break
            if not self._is_kept(entry):
                continue
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
            # Sparse output: drop None values, then sanitise every string field.
            props = {
                k: _clean_value(v) for k, v in props.items() if v is not None
            }
            yield protein_id, "protein", props
            count += 1

    def get_edges(self) -> Generator[tuple[None, str, str, str, dict]]:
        """Yield all protein-related edges."""
        props = self._provenance()
        count = 0

        for uid, entry in self._data.items():
            if self.test_mode and count >= 100:
                break

            gene_pairs, assemblies = self._resolve(entry)
            if not assemblies:
                continue
            protein_id = self._add_prefix("uniprot", uid)

            # Gene_encodes_protein — one per distinct (locus_tag, assembly). A
            # WP_ accession may map to several locus tags within one assembly
            # (80 proteins do, e.g. uniprot:A8WIB5 → PMM1896 + PMM2004): that is
            # genuinely two gene edges.
            seen_genes: set[tuple[str, str]] = set()
            for locus_tag, ncbi_acc in gene_pairs:
                gene_key = (locus_tag, ncbi_acc)
                if gene_key in seen_genes:
                    continue
                seen_genes.add(gene_key)
                gene_id = self._add_prefix("ncbigene", locus_tag)
                # Direction convention: Gene → Protein (source=gene, target=protein)
                yield None, gene_id, protein_id, "Gene_encodes_protein", props

            # Protein_belongs_to_organism — one per assembly the protein belongs
            # to, independent of whether a gene edge exists (the pre-fe5c2bb
            # semantics; the refactor had tied it to the RefSeq join, which is
            # how 25K proteins lost their organism).
            for ncbi_acc in assemblies:
                org_id = self._add_prefix("insdc.gcf", ncbi_acc)
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
