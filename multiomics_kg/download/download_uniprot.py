"""
Standalone UniProt download functions.

Extracted from multiomics_kg/adapters/uniprot_adapter.py so that the download
pipeline (download_genome_data.py) can fetch UniProt data without instantiating
the full adapter.  The adapter delegates to these functions to avoid code
duplication.

Preprocessing (integer conversion, string sanitisation, GO ID extraction) is
handled downstream by build_protein_annotations.py.
"""

import json
import logging
from contextlib import ExitStack
from enum import Enum
from pathlib import Path
from time import time

from tqdm import tqdm


class UniprotNodeField(Enum):
    """UniProt API field names used by the download pipeline."""

    PRIMARY_GENE_NAME = "gene_primary"
    LENGTH = "length"
    MASS = "mass"
    ORGANISM = "organism_name"
    ORGANISM_ID = "organism_id"
    PROTEIN_NAMES = "protein_name"
    EC = "ec"
    PROTEIN_GENE_NAMES = "gene_names"
    SEQUENCE = "sequence"
    GENE_ORDERED_LOCUS = "gene_oln"
    cc_catalytic_activity = "cc_catalytic_activity"
    cc_cofactor = "cc_cofactor"
    cc_function = "cc_function"
    cc_pathway = "cc_pathway"
    annotation_score = "annotation_score"
    cc_caution = "cc_caution"
    keywordid = "keywordid"
    keyword = "keyword"
    reviewed = "reviewed"
    cc_interaction = "cc_interaction"
    CELLULAR_COMPONENT = "go_c"
    BIOLOGICAL_PROCESS = "go_p"
    MOLECULAR_FUNCTION = "go_f"
    CELLULAR_COMPONENT_ID = "go_c_id"
    BIOLOGICAL_PROCESS_ID = "go_p_id"
    MOLECULAR_FUNCTION_ID = "go_f_id"
    ft_transmem = "ft_transmem"
    ft_signal = "ft_signal"
    cc_domain = "cc_domain"
    ft_motif = "ft_motif"
    protein_families = "protein_families"
    xref_refseq = "xref_refseq"
    xref_string = "xref_string"
    xref_eggnog = "xref_eggnog"
    xref_pfam = "xref_pfam"
    PROTEOME = "xref_proteomes"
    ENTREZ_GENE_IDS = "xref_geneid"
    KEGG_IDS = "xref_kegg"
    SUBCELLULAR_LOCATION = "subcellular_location"
    PROTT5_EMBEDDING = "prott5_embedding"
    ESM2_EMBEDDING = "esm2_embedding"
    NT_EMBEDDING = "nt_embedding"

    @classmethod
    def get_split_fields(cls) -> list:
        return [
            cls.PROTEOME.value,
            cls.PROTEIN_GENE_NAMES.value,
            cls.EC.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.KEGG_IDS.value,
            cls.CELLULAR_COMPONENT.value,
            cls.MOLECULAR_FUNCTION.value,
            cls.BIOLOGICAL_PROCESS.value,
        ]

    @classmethod
    def get_nonuniprot_api_fields(cls) -> list:
        return [
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
            cls.NT_EMBEDDING.value,
            cls.CELLULAR_COMPONENT_ID.value,
            cls.MOLECULAR_FUNCTION_ID.value,
            cls.BIOLOGICAL_PROCESS_ID.value,
        ]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default node fields used by the download pipeline (rev=False run).
# Embeddings are intentionally excluded — those are adapter-only.
# ---------------------------------------------------------------------------

DEFAULT_NODE_FIELDS: list[str] = [
    UniprotNodeField.PRIMARY_GENE_NAME.value,
    UniprotNodeField.LENGTH.value,
    UniprotNodeField.MASS.value,
    UniprotNodeField.ORGANISM.value,
    UniprotNodeField.ORGANISM_ID.value,
    UniprotNodeField.PROTEIN_NAMES.value,
    UniprotNodeField.PROTEIN_GENE_NAMES.value,
    UniprotNodeField.KEGG_IDS.value,
    UniprotNodeField.PROTEOME.value,
    UniprotNodeField.SUBCELLULAR_LOCATION.value,
    UniprotNodeField.EC.value,
    UniprotNodeField.GENE_ORDERED_LOCUS.value,
    UniprotNodeField.cc_catalytic_activity.value,
    UniprotNodeField.cc_cofactor.value,
    UniprotNodeField.cc_function.value,
    UniprotNodeField.cc_pathway.value,
    UniprotNodeField.annotation_score.value,
    UniprotNodeField.cc_caution.value,
    UniprotNodeField.keywordid.value,
    UniprotNodeField.keyword.value,
    UniprotNodeField.reviewed.value,
    UniprotNodeField.cc_interaction.value,
    UniprotNodeField.CELLULAR_COMPONENT.value,
    UniprotNodeField.BIOLOGICAL_PROCESS.value,
    UniprotNodeField.MOLECULAR_FUNCTION.value,
    UniprotNodeField.ft_transmem.value,
    UniprotNodeField.ft_signal.value,
    UniprotNodeField.cc_domain.value,
    UniprotNodeField.ft_motif.value,
    UniprotNodeField.protein_families.value,
    UniprotNodeField.xref_refseq.value,
    UniprotNodeField.xref_string.value,
    UniprotNodeField.xref_eggnog.value,
    UniprotNodeField.xref_pfam.value,
]

# ---------------------------------------------------------------------------
# Core: fetch raw UniProt data via pypath
# ---------------------------------------------------------------------------

def fetch_raw_uniprot(
    organism: int,
    node_fields: list[str],
    rev: bool = False,
    test_mode: bool = False,
) -> tuple[dict, set]:
    """Fetch raw UniProt data for an organism taxid.

    Must be called inside a ``pypath.share.curl.cache_off()`` context when a
    fresh download is required (the caller is responsible for this).

    Args:
        organism: NCBI taxid (integer).
        node_fields: list of UniprotNodeField string values to fetch.
        rev: If True, fetch SwissProt-only (reviewed); False fetches all entries.
        test_mode: If True, limit to the first 100 protein IDs.

    Returns:
        Tuple of (data_dict, uniprot_ids_set).
        data_dict: {field_name: {uniprot_id: value, ...}, ...}
        uniprot_ids_set: set of UniProt accessions fetched.
    """
    from pypath.inputs import uniprot

    embedding_fields = {
        UniprotNodeField.PROTT5_EMBEDDING.value,
        UniprotNodeField.ESM2_EMBEDDING.value,
        UniprotNodeField.NT_EMBEDDING.value,
    }

    logger.info(f"Fetching UniProt IDs for taxid={organism} (rev={rev})...")
    t0 = time()

    uniprot_ids: set = set(uniprot._all_uniprots(organism, rev))
    logger.debug(f"  found {len(uniprot_ids)} UniProt IDs")

    if test_mode:
        uniprot_ids = set(list(uniprot_ids)[:100])

    data: dict = {}
    failed_fields: list[str] = []
    for query_key in tqdm(node_fields, desc="Downloading UniProt fields"):
        if query_key in embedding_fields:
            # embeddings are handled separately by the adapter
            continue

        try:
            if query_key == UniprotNodeField.SUBCELLULAR_LOCATION.value:
                data[query_key] = uniprot.uniprot_locations(organism, rev)
            else:
                data[query_key] = uniprot.uniprot_data(
                    fields=query_key,
                    organism=organism,
                    reviewed=rev,
                )
            logger.debug(f"  field {query_key!r} downloaded")
        except Exception as e:
            logger.warning(f"  field {query_key!r} FAILED: {e}")
            failed_fields.append(query_key)

    if failed_fields:
        logger.warning(f"  {len(failed_fields)} field(s) failed: {failed_fields}")

    elapsed = round((time() - t0) / 60, 2)
    logger.info(f"Acquired UniProt data in {elapsed} mins.")

    return data, uniprot_ids


# ---------------------------------------------------------------------------
# Public entry point: download + save
# ---------------------------------------------------------------------------

def download_uniprot(
    taxid: int,
    cache_dir: Path,
    force: bool = False,
    node_fields: list[str] | None = None,
    rev: bool = False,
) -> bool:
    """Download and cache raw UniProt data for a taxid.

    Output file written to ``cache_dir``:
    - ``uniprot_raw_data.json`` — raw API data

    Preprocessing (integer conversion, string sanitisation, GO ID extraction)
    is handled downstream by build_protein_annotations.py.

    Args:
        taxid: NCBI taxid.
        cache_dir: directory where the output JSON file will be written.
        force: re-download even if the cache file already exists.
        node_fields: fields to fetch; defaults to ``DEFAULT_NODE_FIELDS``.
        rev: reviewed-only flag (False = TrEMBL + SwissProt).

    Returns:
        True if download was performed; False if skipped (cache hit).
    """
    from pypath.share import curl

    if node_fields is None:
        node_fields = DEFAULT_NODE_FIELDS

    raw_path = cache_dir / "uniprot_raw_data.json"

    if not force and raw_path.exists():
        return False

    cache_dir.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        stack.enter_context(curl.cache_off())
        data, uniprot_ids = fetch_raw_uniprot(taxid, node_fields, rev=rev)

    logger.info(f"  Saving raw data ({len(uniprot_ids)} proteins) → {raw_path}")
    with open(raw_path, "w") as f:
        json.dump(data, f, indent=4, default=str)

    ids_path = cache_dir / "uniprot_ids.json"
    with open(ids_path, "w") as f:
        json.dump(sorted(uniprot_ids), f)
    logger.info(f"  Saved {len(uniprot_ids)} protein IDs → {ids_path}")

    return True
