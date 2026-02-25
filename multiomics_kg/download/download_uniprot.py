"""
Standalone UniProt download and preprocessing functions.

Extracted from multiomics_kg/adapters/uniprot_adapter.py so that the download
pipeline (download_genome_data.py) can fetch and preprocess UniProt data without
instantiating the full adapter.  The adapter delegates to these functions to
avoid code duplication.
"""

import json
import logging
from contextlib import ExitStack
from pathlib import Path
from time import time

from tqdm import tqdm

from multiomics_kg.adapters.uniprot_adapter import UniprotNodeField

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
# Helper: extract GO id from a GO term string
# ---------------------------------------------------------------------------

def _extract_go_id(go_term: str) -> str | None:
    """Extract GO id from a GO term string.

    Example input: "aspartate-semialdehyde dehydrogenase activity [GO:0004073]"
    """
    if "GO:" in go_term:
        return go_term.split("GO:")[1].split("]")[0].strip()
    return None


# ---------------------------------------------------------------------------
# Helper: split a single multi-valued UniProt field value
# ---------------------------------------------------------------------------

def _split_field(field_key: str, field_value: str):
    """Split a UniProt field with multiple entries.

    Args:
        field_key: UniprotNodeField string value (e.g. "xref_proteomes")
        field_value: raw string value from the UniProt API

    Returns:
        Processed value (list or scalar string).
    """
    if field_value is None or field_value == "":
        return field_value

    # replace sensitive elements for admin-import
    field_value = field_value.replace("|", ",").replace("'", "^").strip()

    split_dict = {
        UniprotNodeField.PROTEOME.value: ",",
        UniprotNodeField.PROTEIN_GENE_NAMES.value: " ",
    }

    split_char = split_dict.get(field_key, ";")
    field_value = [i.strip() for i in field_value.strip(split_char).split(split_char)]

    # split colons (":") in KEGG field to extract ID part
    if field_key == UniprotNodeField.KEGG_IDS.value:
        field_value = [e.split(":")[1].strip() for e in field_value]

    # take first element only for ENTREZ_GENE_IDS
    if field_key == UniprotNodeField.ENTREZ_GENE_IDS.value:
        field_value = field_value[0]

    # unwrap single-element lists to scalar
    if isinstance(field_value, list) and len(field_value) == 1:
        field_value = field_value[0]

    return field_value


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
    for query_key in tqdm(node_fields, desc="Downloading UniProt fields"):
        if query_key in embedding_fields:
            # embeddings are handled separately by the adapter
            continue

        if query_key == UniprotNodeField.SUBCELLULAR_LOCATION.value:
            data[query_key] = uniprot.uniprot_locations(organism, rev)
        else:
            data[query_key] = uniprot.uniprot_data(
                fields=query_key,
                organism=organism,
                reviewed=rev,
            )

        logger.debug(f"  field {query_key!r} downloaded")

    elapsed = round((time() - t0) / 60, 2)
    logger.info(f"Acquired UniProt data in {elapsed} mins.")

    return data, uniprot_ids


# ---------------------------------------------------------------------------
# Core: preprocess raw UniProt data
# ---------------------------------------------------------------------------

def preprocess_uniprot_data(
    data: dict,
    node_fields: list[str],
) -> tuple[dict, set]:
    """Preprocess raw UniProt data in-place.

    Applies:
    - Integer conversion for numeric fields (LENGTH, MASS, ORGANISM_ID)
    - String sanitisation (| → , and ' → ^)
    - Field splitting for multi-valued fields
    - Subcellular location extraction
    - GO ID extraction (go_c_id, go_p_id, go_f_id fields added)

    Args:
        data: raw data dict as returned by ``fetch_raw_uniprot``.
        node_fields: list of UniprotNodeField string values (same list used
            for fetching).

    Returns:
        Tuple of (processed_data, locations_set).
        locations_set collects all unique subcellular location strings.
    """
    logger.info("Preprocessing UniProt data.")

    nonuniprot_api_fields = set(UniprotNodeField.get_nonuniprot_api_fields())
    split_fields = set(UniprotNodeField.get_split_fields())
    locations: set = set()

    int_fields = {
        UniprotNodeField.LENGTH.value,
        UniprotNodeField.MASS.value,
        UniprotNodeField.ORGANISM_ID.value,
    }

    for arg in tqdm(node_fields, desc="Processing UniProt fields"):

        if arg not in data:
            continue

        if arg in nonuniprot_api_fields:
            pass  # skip embedding / derived fields

        elif arg in int_fields:
            for protein, val in data[arg].items():
                data[arg][protein] = int(str(val).replace(",", ""))

        elif arg == UniprotNodeField.SUBCELLULAR_LOCATION.value:
            for protein, val in data[arg].items():
                parsed = []
                for element in val:
                    loc = (
                        str(element.location)
                        .replace("'", "")
                        .replace("[", "")
                        .replace("]", "")
                        .strip()
                    )
                    parsed.append(loc)
                    locations.add(loc)
                data[arg][protein] = parsed

        elif arg in split_fields:
            for protein, val in data[arg].items():
                data[arg][protein] = _split_field(arg, val)

        else:
            for protein, val in data[arg].items():
                data[arg][protein] = (
                    val.replace("|", ",").replace("'", "^").strip()
                )

    # Extract GO IDs into separate fields (go_c_id, go_p_id, go_f_id)
    go_fields = [
        UniprotNodeField.CELLULAR_COMPONENT.value,
        UniprotNodeField.BIOLOGICAL_PROCESS.value,
        UniprotNodeField.MOLECULAR_FUNCTION.value,
    ]
    for go_field in go_fields:
        if go_field not in data:
            continue
        go_id_field = go_field + "_id"
        data[go_id_field] = {}
        for protein, terms in data[go_field].items():
            data[go_id_field][protein] = [
                _extract_go_id(t) for t in terms if t and "GO:" in t
            ]

    return data, locations


# ---------------------------------------------------------------------------
# Public entry point: download + preprocess + save
# ---------------------------------------------------------------------------

def download_uniprot(
    taxid: int,
    cache_dir: Path,
    force: bool = False,
    node_fields: list[str] | None = None,
    rev: bool = False,
) -> bool:
    """Download, preprocess, and cache UniProt data for a taxid.

    Output files written to ``cache_dir``:
    - ``uniprot_raw_data.json``      — raw API data
    - ``uniprot_preprocess_data.json`` — preprocessed data

    Args:
        taxid: NCBI taxid.
        cache_dir: directory where output JSON files will be written.
        force: re-download even if cache files already exist.
        node_fields: fields to fetch; defaults to ``DEFAULT_NODE_FIELDS``.
        rev: reviewed-only flag (False = TrEMBL + SwissProt).

    Returns:
        True if download was performed; False if skipped (cache hit).
    """
    from pypath.share import curl

    if node_fields is None:
        node_fields = DEFAULT_NODE_FIELDS

    raw_path = cache_dir / "uniprot_raw_data.json"
    pre_path = cache_dir / "uniprot_preprocess_data.json"

    if not force and raw_path.exists() and pre_path.exists():
        return False

    cache_dir.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        stack.enter_context(curl.cache_off())
        data, uniprot_ids = fetch_raw_uniprot(taxid, node_fields, rev=rev)

    logger.info(f"  Saving raw data ({len(uniprot_ids)} proteins) → {raw_path}")
    with open(raw_path, "w") as f:
        json.dump(data, f, indent=4, default=str)

    data, _locations = preprocess_uniprot_data(data, node_fields)

    logger.info(f"  Saving preprocessed data → {pre_path}")
    with open(pre_path, "w") as f:
        json.dump(data, f, indent=4, default=str)

    return True
