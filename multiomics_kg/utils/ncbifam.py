"""Pure-Python parsing of the NCBIfam/PGAP HMM reference (``hmm_PGAP.tsv``).

No filesystem, no network — the download + caching orchestrator lives in
``multiomics_kg/download/build_ncbifam_reference.py``. This module is the
unit-testable core: it turns NCBI's ``hmm_PGAP.tsv`` rows into the committed
reference dict consumed by ``ncbifam_adapter`` at KG-build time.

Source file (NCBI FTP, ``https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv``):
tab-separated, one header line, columns (verified 2026-08-17)::

    #ncbi_accession, source_identifier, label, sequence_cutoff, domain_cutoff,
    hmm_length, family_type, for_structural_annotation, for_naming,
    for_AMRFinder, product_name, gene_symbol, gene_synonyms, ec_numbers,
    go_terms, pmids, taxonomic_range, taxonomic_range_name,
    taxonomic_rank_name, n_refseq_protein_hits, source, name_orig, hmm_name,
    comment

``#ncbi_accession`` carries a version suffix (``NF000001.1``, ``TIGR00198.1``)
that is stripped to the bare accession used as the KG's NCBIfam node id.
``gene_synonyms``/``ec_numbers``/``go_terms``/``pmids`` are comma-delimited
multi-valued cells (empirically verified against a downloaded copy of
``hmm_PGAP.tsv`` — e.g. ``gene_synonyms`` cell ``"cybA,dhsC"``; no
space-delimited cells were observed). ``go_terms`` values may be bare GO ids
or ``GO:NNNNNNN`` — stored **raw**, consumer normalizes (interpro-reference
precedent).

Combined result: ``{unversioned_acc: {"name", "family_type"[, "gene_symbol"]
[, "description"][, "gene_synonyms"][, "pmids"][, "ec_numbers"][, "go_terms"]}}``
— all optional keys are **sparse**: omitted (not an empty string/list) when the
source cell is empty.
"""

from __future__ import annotations

from collections.abc import Iterable

# family_type values that flag a hypothetical (uncharacterized) family — the
# product_name for these is a placeholder ("hypothetical protein"-shaped),
# not a real functional annotation. Verified against the full hmm_PGAP.tsv
# family_type distribution (2026-08-17): hypoth_equivalog 434,
# hypoth_equivalog_domain 34 (of 38,394 total rows).
HYPOTH_FAMILY_TYPES = frozenset({"hypoth_equivalog", "hypoth_equivalog_domain"})

_MULTI_VALUE_DELIM = ","


def _split_multi(raw: str) -> list[str]:
    """Split a comma-delimited multi-valued TSV cell into a list of tokens."""
    return [tok.strip() for tok in raw.split(_MULTI_VALUE_DELIM) if tok.strip()]


def parse_hmm_pgap_rows(rows: Iterable[dict]) -> dict[str, dict]:
    """Parse ``hmm_PGAP.tsv`` DictReader rows into the committed reference dict.

    ``rows`` is any iterable of dicts keyed by the TSV's column names (as
    produced by ``csv.DictReader(fh, delimiter="\\t")``). The leading ``#`` in
    ``#ncbi_accession`` is exactly what ``csv.DictReader`` returns for that
    header cell, so callers must not strip it themselves.

    Returns ``{unversioned_acc: {"name", "family_type"[, "gene_symbol"]
    [, "description"][, "gene_synonyms"][, "pmids"][, "ec_numbers"]
    [, "go_terms"]}}``. Rows with a missing/empty accession are skipped.
    """
    out: dict[str, dict] = {}
    for row in rows:
        raw_acc = (row.get("#ncbi_accession") or "").strip()
        if not raw_acc:
            continue
        acc = raw_acc.split(".", 1)[0]

        entry: dict = {
            "name": (row.get("product_name") or "").strip(),
            "family_type": (row.get("family_type") or "").strip(),
        }

        gene_symbol = (row.get("gene_symbol") or "").strip()
        if gene_symbol:
            entry["gene_symbol"] = gene_symbol

        description = (row.get("comment") or "").strip()
        if description:
            entry["description"] = description

        gene_synonyms = _split_multi(row.get("gene_synonyms") or "")
        if gene_synonyms:
            entry["gene_synonyms"] = gene_synonyms

        pmids = _split_multi(row.get("pmids") or "")
        if pmids:
            entry["pmids"] = pmids

        ec_numbers = _split_multi(row.get("ec_numbers") or "")
        if ec_numbers:
            entry["ec_numbers"] = ec_numbers

        go_terms = _split_multi(row.get("go_terms") or "")
        if go_terms:
            entry["go_terms"] = go_terms

        out[acc] = entry
    return out
