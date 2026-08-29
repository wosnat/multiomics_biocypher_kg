"""Build gene_mapping.csv from NCBI GFF/GBFF and Cyanorak GFF/GBK files.

Provides standalone functions for loading and merging gene annotation data.
Used by both CyanorakNcbi adapter and the download pipeline (download_genome_data.py).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from biocypher._logger import logger
import gffpandas.gffpandas as gffpd
from urllib.parse import unquote


# ── column configuration ──────────────────────────────────────────────────────

_CYANORAK_COLS = [
    'start', 'end', 'strand',
    'ID', 'Name', 'Ontology_term', 'cluster_number',
    'cyanorak_Role', 'cyanorak_Role_description', 'eggNOG',
    'eggNOG_description', 'kegg', 'kegg_description',
    'ontology_term_description', 'product', 'protein_domains',
    'protein_domains_description', 'tIGR_Role', 'tIGR_Role_description',
    'locus_tag',
]

_FINAL_MERGE_RENAME = {
    'Name_ncbi': 'gene_names',
    'Name_cyanorak': 'gene_names_cyanorak',
    'start_ncbi': 'start', 'end_ncbi': 'end', 'strand_ncbi': 'strand',
    'start_cyanorak': 'start_cyanorak', 'end_cyanorak': 'end_cyanorak',
    'strand_cyanorak': 'strand_cyanorak',
    'product_ncbi': 'product', 'product_cyanorak': 'product_cyanorak',
    'ID': 'locus_tag_cyanorak',
}


def _get_ncbi_cols_rename_map() -> dict[str, str]:
    """Build the column rename map for NCBI GFF merged gene+CDS DataFrame."""
    ncbi_cols_to_keep = [
        'seq_id_cds',
        'Name_gene',
        'gene_gene',
        'locus_tag_cds',
        'old_locus_tag_gene',
        'source_cds',
        'start_cds',
        'end_cds',
        'strand_cds',
        'Note_cds',
        'exception_cds',
        'inference_cds',
        'product_cds',
        'protein_id_cds',
        'Ontology_term_cds',
        'go_component_cds',
        'go_function_cds',
        'go_process_cds',
        'gene_synonym_gene',
    ]
    col_rename_map = {c: c.replace('_gene', '').replace('_cds', '') for c in ncbi_cols_to_keep}
    col_rename_map['seq_id_cds'] = 'seqid'
    col_rename_map['locus_tag_cds'] = 'locus_tag_ncbi'
    col_rename_map['old_locus_tag_gene'] = 'locus_tag'
    return col_rename_map


# ── GFF/GBK parsing ───────────────────────────────────────────────────────────

def load_gff(gff_file: str) -> pd.DataFrame:
    """Load a GFF3 file and return a DataFrame with CDS and gene rows."""
    gff_annotation = gffpd.read_gff3(gff_file)
    gff_df = gff_annotation.attributes_to_columns()
    return gff_df.loc[gff_df.type.isin(['CDS', 'gene'])]


def ncbi_merge_cds_and_gene_entries(ncbi_gff_df: pd.DataFrame) -> pd.DataFrame:
    """Merge gene and CDS entries from NCBI GFF DataFrame on ID / Parent attributes."""
    gene_df = ncbi_gff_df.loc[ncbi_gff_df.type.isin(['gene'])]
    cds_df = ncbi_gff_df.loc[ncbi_gff_df.type.isin(['CDS'])]
    df = pd.merge(gene_df, cds_df, left_on='ID', right_on='Parent', suffixes=['_gene', '_cds'])

    ncbi_col_rename_map = _get_ncbi_cols_rename_map()
    potential_cols = [
        c for c in df.columns
        if (c not in ncbi_col_rename_map) and (df[c].dropna().nunique() > 10)
    ]
    logger.info(f"NCBI, potential columns : {potential_cols}")

    # Only keep columns that exist in the merged DataFrame (e.g. old_locus_tag_gene
    # may be absent for genomes that lack old_locus_tag in their GFF)
    available_cols = {k: v for k, v in ncbi_col_rename_map.items() if k in df.columns}
    df_filter = df[list(available_cols.keys())].rename(columns=available_cols)

    # Handle multiple old_locus_tag values (URL-encoded comma-separated in GFF)
    # e.g. "PMT0003%2CPMT_0003%2CRG24_RS00015" -> ["PMT0003", "PMT_0003", "RG24_RS00015"]
    if 'locus_tag' in df_filter.columns:
        df_filter['locus_tag'] = df_filter['locus_tag'].apply(
            lambda x: unquote(str(x)) if pd.notna(x) else x
        )
        # Store all old_locus_tags as a comma-separated string
        df_filter['old_locus_tags'] = df_filter['locus_tag']
    return df_filter


def _get_cyanorak_id(rec) -> str | None:
    """Extract the Cyanorak ID from a GenBank CDS feature of a Cyanorak GBK file.

    Returns None if the feature has no 'note' qualifier or no matching note.
    """
    notes = rec.qualifiers.get('note', [])
    if not notes:
        return None
    note_name = 'cyanorak ORF Id:'
    cyanorak_ids = [i for i in notes if i.startswith(note_name)]
    if not cyanorak_ids:
        return None
    if len(cyanorak_ids) > 1:
        logger.warning(f"Multiple cyanorak IDs found for feature {rec.qualifiers.get('locus_tag', ['?'])[0]}: {cyanorak_ids}")
    return cyanorak_ids[0].replace(note_name, '').strip()


def _get_cyanorak_id_map_from_gbk(gbk_file: str) -> dict[str, str]:
    """Create a mapping from Cyanorak ID to locus_tag from a Cyanorak GBK file."""
    result: dict[str, str] = {}
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            cyanorak_id = _get_cyanorak_id(feature)
            if cyanorak_id is None:
                continue
            for locus_tag in feature.qualifiers.get('locus_tag', []):
                result[cyanorak_id] = locus_tag
    return result


def _get_ec_numbers_from_gbff(gbff_file: str) -> dict[str, list[str]]:
    """Extract EC numbers from NCBI GBFF file, keyed by locus_tag.

    Args:
        gbff_file: Path to the NCBI genomic.gbff file.

    Returns:
        Dict mapping locus_tag -> list of EC number strings (e.g. ["1.1.1.1"]).
    """
    locus_tag_to_ec: dict[str, list[str]] = {}
    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            ec_numbers = feature.qualifiers.get("EC_number", [])
            if not ec_numbers:
                continue
            for locus_tag in feature.qualifiers.get("locus_tag", []):
                if locus_tag in locus_tag_to_ec:
                    locus_tag_to_ec[locus_tag].extend(ec_numbers)
                else:
                    locus_tag_to_ec[locus_tag] = list(ec_numbers)
    return locus_tag_to_ec


# ── position-based fallback merge ─────────────────────────────────────────────

# Cyanorak-specific columns to copy during position fallback merge.
# These are columns unique to Cyanorak (no _cyanorak suffix needed) or
# shared columns that got suffixed with _cyanorak during the outer merge.
_CYAN_COLUMNS_TO_COPY = [
    'ID',
    'start_cyanorak', 'end_cyanorak', 'strand_cyanorak',
    'Name_cyanorak', 'Ontology_term_cyanorak', 'product_cyanorak',
    'cluster_number',
    'cyanorak_Role', 'cyanorak_Role_description',
    'eggNOG', 'eggNOG_description',
    'kegg', 'kegg_description',
    'ontology_term_description',
    'protein_domains', 'protein_domains_description',
    'tIGR_Role', 'tIGR_Role_description',
]


def _three_prime_end(start, end, strand):
    """Genomic coordinate of the stop codon: ``end`` on +, ``start`` on -.

    Works element-wise on aligned numpy arrays / Series as well as scalars.
    """
    return np.where(np.asarray(strand) == '-', np.asarray(start, dtype=float),
                    np.asarray(end, dtype=float))


def _seqids(df: pd.DataFrame) -> pd.Series:
    """Per-row contig id ('' when the frame has no ``seqid`` column)."""
    if 'seqid' in df.columns:
        return df['seqid'].fillna('')
    return pd.Series([''] * len(df), index=df.index)


def _contig_offsets(ncbi_sourced: pd.DataFrame) -> dict:
    """Derive, per NCBI contig, the constant added to NCBI coordinates to obtain
    Cyanorak coordinates.

    Cyanorak stores a draft genome as ONE concatenated record (NCBI contigs in
    order), so contig-relative NCBI coordinates only compare after a per-contig
    shift (PAC1 contig 5: +229,662; SB contig 2: +420,116; every closed genome:
    0). The shift is measured on genes the locus_tag merge already paired, at
    the 3' end (stop codon), because start codons legitimately differ. A contig
    is usable when a strict majority of its matched genes agree on one offset
    (a contig backed by fewer than 3 matched genes is used but warned about);
    an unmatched-only contig is usable at offset 0 only when the assembly is a
    single sequence. Contigs without a derivable offset are left out — with a
    warning naming them — and the fallback never compares coordinates on them.
    """
    seqids = _seqids(ncbi_sourced)
    distinct = list(pd.unique(seqids))
    coord_cols = ['start_ncbi', 'end_ncbi', 'start_cyanorak', 'end_cyanorak']
    matched = ncbi_sourced[
        ncbi_sourced['ID'].notna() & ncbi_sourced[coord_cols].notna().all(axis=1)
    ]
    offsets: dict = {}
    for sid in distinct:
        m = matched[seqids.loc[matched.index] == sid]
        if m.empty:
            if len(distinct) == 1:
                offsets[sid] = 0
            else:
                n_unmatched = int(((seqids == sid) & ncbi_sourced['ID'].isna()).sum())
                logger.warning(
                    f"Position fallback: contig {sid!r} has no locus_tag-matched "
                    f"gene to derive a Cyanorak offset from — its {n_unmatched} "
                    f"unmatched NCBI genes are not considered"
                )
            continue
        diff = (
            _three_prime_end(m['start_cyanorak'], m['end_cyanorak'], m['strand_ncbi'])
            - _three_prime_end(m['start_ncbi'], m['end_ncbi'], m['strand_ncbi'])
        )
        vals, counts = np.unique(diff, return_counts=True)
        best = int(np.argmax(counts))
        if counts[best] * 2 > len(diff) or len(diff) == 1:
            offsets[sid] = int(vals[best])
            if len(diff) < 3:
                logger.warning(
                    f"Position fallback: contig {sid!r} offset {offsets[sid]} rests "
                    f"on only {len(diff)} matched gene(s)"
                )
        else:
            logger.warning(
                f"Position fallback: contig {sid!r} has no majority Cyanorak "
                f"offset ({len(diff)} matched genes) — skipped"
            )
    return offsets


def _position_fallback_merge(
    ncbi_sourced: pd.DataFrame,
    cyan_only: pd.DataFrame,
) -> tuple[int, set[str]]:
    """Merge unmatched Cyanorak entries into NCBI entries by shared reading frame.

    After the initial locus_tag-based merge, some genes remain unmatched because
    NCBI old_locus_tag is incomplete (e.g., MIT9313 PMT_0107 vs Cyanorak PMT0107).
    This fallback matches by genomic position.

    A bacterial CDS is identified by its stop codon and reading frame, not by
    its start: re-annotation routinely moves the start codon (PGAP vs. Cyanorak
    differ by up to ~250 bp on MIT9313), so the 5' boundary and the resulting
    overlap ratio are NOT identity evidence and are not used. Cyanorak also
    sometimes reports the CDS without its stop codon (3' end 3 bp short) or
    truncated a few codons early, still in frame.

    Criteria for a match (a Cyanorak call C against an NCBI call N):
    - Same contig, via the per-contig coordinate offset from ``_contig_offsets``
      (a draft genome's NCBI contigs are contig-relative; Cyanorak's single
      concatenated record is not — the old fallback ignored this and produced
      cross-contig false merges on PAC1 and SB)
    - Same strand
    - In frame: the 3' ends differ by a multiple of 3 (3' end = genomic
      ``end`` on the + strand, genomic ``start`` on the - strand — strand-aware;
      the old strand-blind ``end`` gate rejected every - strand start re-call)
    - C's 3' end lies inside N's interval. In frame and inside means no stop
      codon can separate the two 3' ends, so C is the same ORF as N (up to the
      start-codon choice / a missing stop codon), by construction.
    - 1:1 only, in both directions; conflicts (one NCBI gene matching multiple
      Cyanorak calls, e.g. a fusion/split disagreement, or one Cyanorak call
      sitting in frame inside two nested NCBI genes) are skipped with a warning

    No distance bound is applied inside N's interval on purpose: an in-frame
    Cyanorak ORF ending inside an NCBI span is the same locus even when N is a
    nonsense pseudogene or frameshifted CDS that Cyanorak called as its intact
    upstream ORF — pairing it with N's locus tag is the correct outcome.

    Verified 2026-08-29 by translating both coordinate sets from the Cyanorak
    genome over 21 strains: every candidate pair (~700, incl. the 336 the old
    overlap/±bp thresholds missed and the 11 stop-codon-excluded ones they
    caught) encodes the same CDS; zero false merges.

    Modifies ncbi_sourced in-place by filling Cyanorak columns for matched rows.

    Returns:
        (n_merged, consumed_cyan_locus_tags)
    """
    ncbi_unmatched = ncbi_sourced[ncbi_sourced['ID'].isna()].copy()
    if ncbi_unmatched.empty or cyan_only.empty:
        return 0, set()

    # Shift NCBI coordinates into Cyanorak's coordinate space, per contig;
    # rows on contigs with no derivable offset are dropped from consideration.
    offsets = _contig_offsets(ncbi_sourced)
    sid = _seqids(ncbi_unmatched)
    ncbi_unmatched = ncbi_unmatched[sid.isin(list(offsets))].copy()
    if ncbi_unmatched.empty:
        return 0, set()
    shift = sid.loc[ncbi_unmatched.index].map(offsets).astype(float)
    ncbi_unmatched['_start_shifted'] = ncbi_unmatched['start_ncbi'].astype(float) + shift
    ncbi_unmatched['_end_shifted'] = ncbi_unmatched['end_ncbi'].astype(float) + shift

    # Build candidate pairs keyed by NCBI locus_tag_ncbi to detect conflicts;
    # also count how many NCBI genes each Cyanorak call lands in (reverse conflicts)
    ncbi_to_cyan: dict[str, list[tuple[str, int]]] = {}
    cyan_hits: dict[str, int] = {}

    for _cyan_idx, c in cyan_only.iterrows():
        c_start = c['start_cyanorak']
        c_end = c['end_cyanorak']
        c_strand = c['strand_cyanorak']
        if pd.isna(c_start) or pd.isna(c_end) or pd.isna(c_strand):
            continue

        same_strand = ncbi_unmatched[ncbi_unmatched['strand_ncbi'] == c_strand]
        if same_strand.empty:
            continue

        # The 3' end (stop codon) is the genomic end on +, genomic start on -
        n_starts = same_strand['_start_shifted'].values
        n_ends = same_strand['_end_shifted'].values
        c_stop = float(_three_prime_end(c_start, c_end, c_strand))
        n_stop = _three_prime_end(n_starts, n_ends, c_strand)
        in_frame = (np.abs(n_stop - c_stop) % 3) == 0
        inside = (n_starts <= c_stop) & (c_stop <= n_ends)
        mask = in_frame & inside

        for i in np.where(mask)[0]:
            n_row = same_strand.iloc[i]
            n_lt_ncbi = n_row['locus_tag_ncbi']
            ncbi_to_cyan.setdefault(n_lt_ncbi, []).append(
                (c['locus_tag'], same_strand.index[i])
            )
            cyan_hits[c['locus_tag']] = cyan_hits.get(c['locus_tag'], 0) + 1

    # Resolve: skip conflicts in either direction
    consumed_cyan = set()
    n_merged = 0

    for n_lt_ncbi, candidates in ncbi_to_cyan.items():
        if len(candidates) > 1:
            logger.warning(
                f"Position fallback: skipping {n_lt_ncbi} — "
                f"matches {len(candidates)} Cyanorak entries: "
                f"{[c[0] for c in candidates]}"
            )
            continue

        cyan_lt, ncbi_idx = candidates[0]
        # Reverse conflict: one Cyanorak call in frame inside two nested NCBI
        # genes — refuse both sides, same policy as the forward case
        if cyan_hits.get(cyan_lt, 0) > 1:
            logger.warning(
                f"Position fallback: skipping {n_lt_ncbi} — Cyanorak {cyan_lt} "
                f"also lands in {cyan_hits[cyan_lt] - 1} other NCBI gene(s)"
            )
            continue

        # Copy Cyanorak columns into the NCBI row
        cyan_row = cyan_only[cyan_only['locus_tag'] == cyan_lt].iloc[0]
        for col in _CYAN_COLUMNS_TO_COPY:
            if col in ncbi_sourced.columns and col in cyan_only.columns:
                ncbi_sourced.at[ncbi_idx, col] = cyan_row[col]

        # Preserve the Cyanorak locus_tag in old_locus_tags so it remains
        # discoverable by gene_id_mapping (papers may use this form).
        if 'old_locus_tags' not in ncbi_sourced.columns:
            ncbi_sourced['old_locus_tags'] = np.nan  # GFF had no old_locus_tag at all
        existing_olt = str(ncbi_sourced.at[ncbi_idx, 'old_locus_tags'])
        if pd.notna(ncbi_sourced.at[ncbi_idx, 'old_locus_tags']):
            if cyan_lt not in existing_olt:
                ncbi_sourced.at[ncbi_idx, 'old_locus_tags'] = f"{existing_olt},{cyan_lt}"
        else:
            ncbi_sourced.at[ncbi_idx, 'old_locus_tags'] = cyan_lt

        # Add a note flagging this as a position-based merge
        # locus_tag is still NaN here for NCBI genes without old_locus_tag
        # (filled from locus_tag_ncbi later), so name the gene the same way
        target = ncbi_sourced.at[ncbi_idx, 'locus_tag']
        if pd.isna(target):
            target = ncbi_sourced.at[ncbi_idx, 'locus_tag_ncbi']
        ncbi_sourced.at[ncbi_idx, 'position_merge_note'] = (
            f"position_merge:{cyan_lt}→{target}"
        )

        consumed_cyan.add(cyan_lt)
        n_merged += 1

    return n_merged, consumed_cyan


# ── high-level loaders ────────────────────────────────────────────────────────

def load_gff_from_ncbi_only(ncbi_gff_file: str) -> pd.DataFrame:
    """Load gene data from NCBI GFF3 file only (no Cyanorak data)."""
    ncbi_df = load_gff(ncbi_gff_file)
    ncbi_df = ncbi_merge_cds_and_gene_entries(ncbi_df)

    # Prefer old_locus_tag (canonical ID), fall back to locus_tag_ncbi (RS-format)
    if 'locus_tag' in ncbi_df.columns:
        ncbi_df['locus_tag'] = ncbi_df['locus_tag'].fillna(ncbi_df['locus_tag_ncbi'])
    else:
        ncbi_df['locus_tag'] = ncbi_df['locus_tag_ncbi']

    ncbi_df = ncbi_df.rename(columns={'Name': 'gene_names'})

    n_before = len(ncbi_df)
    ncbi_df = ncbi_df.dropna(subset=['locus_tag'])
    n_dropped = n_before - len(ncbi_df)
    if n_dropped > 0:
        logger.warning(f"Dropped {n_dropped} genes with no identifiable locus_tag")

    logger.info(f"NCBI-only gene load: {len(ncbi_df)} genes from {ncbi_gff_file}")
    return ncbi_df


def load_gff_from_ncbi_and_cyanorak(
    ncbi_gff_file: str,
    cyan_gff_file: str,
    cyan_gbk_file: str,
) -> pd.DataFrame:
    """Load and merge gene data from NCBI GFF and Cyanorak GFF/GBK files."""
    cyan_df = load_gff(cyan_gff_file)
    ncbi_df = load_gff(ncbi_gff_file)
    ncbi_df = ncbi_merge_cds_and_gene_entries(ncbi_df)

    cyanID2locus_tag_map = _get_cyanorak_id_map_from_gbk(cyan_gbk_file)
    cyan_df['locus_tag'] = cyan_df['ID'].map(cyanID2locus_tag_map)

    # Only keep relevant Cyanorak columns that are present in this GFF
    cols_to_keep = [c for c in _CYANORAK_COLS if c in cyan_df.columns]
    cyan_df = cyan_df[cols_to_keep]

    if cyan_df['locus_tag'].isna().any():
        logger.warning(
            f"{cyan_df['locus_tag'].isna().sum()} Cyanorak entries have no locus_tag after mapping. "
            f"These will be dropped from the merged dataset."
        )
    cyan_df = cyan_df.dropna(subset=['locus_tag'])

    # Handle multiple old_locus_tag values: explode to one row per tag for merge.
    # If GFF lacks old_locus_tag, fall back to locus_tag_ncbi
    if 'locus_tag' not in ncbi_df.columns:
        ncbi_df['locus_tag'] = ncbi_df['locus_tag_ncbi']
    ncbi_df['locus_tag'] = ncbi_df['locus_tag'].str.split(',')
    ncbi_exploded = ncbi_df.explode('locus_tag')
    ncbi_exploded['locus_tag'] = ncbi_exploded['locus_tag'].str.strip()

    # Outer merge to include genes from either source, not just intersection
    merge_df = pd.merge(
        ncbi_exploded, cyan_df, on='locus_tag', how='outer',
        suffixes=['_ncbi', '_cyanorak'],
    )

    # --- Deduplication ---
    # NCBI genes may have multiple old_locus_tags producing multiple rows
    # after explode. Strategy: split into NCBI-sourced vs Cyanorak-only,
    # dedup each independently, then recombine.
    has_ncbi = merge_df['locus_tag_ncbi'].notna()
    ncbi_sourced = merge_df[has_ncbi].copy()
    cyan_only = merge_df[~has_ncbi].copy()

    # For NCBI genes: prefer rows that also matched Cyanorak (have ID).
    if not ncbi_sourced.empty:
        ncbi_sourced = ncbi_sourced.sort_values(
            by='ID', ascending=True, na_position='last',
        )
        dup_mask = ncbi_sourced.duplicated(subset=['locus_tag_ncbi'], keep=False)
        if dup_mask.any():
            dup_genes = ncbi_sourced.loc[dup_mask, 'locus_tag_ncbi'].unique()
            logger.warning(
                f"{len(dup_genes)} NCBI genes have duplicate rows after merge. "
                f"Keeping best match for each."
            )
        ncbi_sourced = ncbi_sourced.drop_duplicates(subset=['locus_tag_ncbi'], keep='first')

    # For Cyanorak-only genes: remove any whose locus_tag was already captured via NCBI.
    if not cyan_only.empty:
        matched_locus_tags = set(
            ncbi_sourced.loc[ncbi_sourced['ID'].notna(), 'locus_tag'].dropna()
        ) if not ncbi_sourced.empty else set()
        cyan_only = cyan_only[~cyan_only['locus_tag'].isin(matched_locus_tags)]
        cyan_only = cyan_only.drop_duplicates(subset=['locus_tag'], keep='first')

    # --- Position-based fallback for unmatched entries ---
    # When NCBI old_locus_tag doesn't include the Cyanorak locus_tag form
    # (e.g., MIT9313 PMT_0107 vs Cyanorak PMT0107), match by shared reading frame.
    n_ncbi_matched = int(ncbi_sourced['ID'].notna().sum()) if not ncbi_sourced.empty else 0
    n_pos_merged = 0
    if not ncbi_sourced.empty and not cyan_only.empty:
        n_pos_merged, consumed_cyan_lts = _position_fallback_merge(
            ncbi_sourced, cyan_only,
        )
        if consumed_cyan_lts:
            cyan_only = cyan_only[~cyan_only['locus_tag'].isin(consumed_cyan_lts)]

    n_ncbi_only = len(ncbi_sourced) - n_ncbi_matched - n_pos_merged
    n_cyan_only = len(cyan_only)

    merge_df = pd.concat([ncbi_sourced, cyan_only], ignore_index=True)

    # For NCBI genes without old_locus_tag, locus_tag is NaN — use locus_tag_ncbi.
    merge_df['locus_tag'] = merge_df['locus_tag'].fillna(merge_df['locus_tag_ncbi'])

    n_before = len(merge_df)
    merge_df = merge_df.dropna(subset=['locus_tag'])
    n_dropped = n_before - len(merge_df)
    if n_dropped > 0:
        logger.warning(f"Dropped {n_dropped} genes with no identifiable locus_tag")

    logger.info(
        f"NCBI genes {ncbi_df.shape[0]} + Cyanorak genes {cyan_df.shape[0]} → "
        f"{len(merge_df)} total genes after deduplication"
    )
    logger.info(
        f"Gene merge: {n_ncbi_matched} matched both sources, "
        f"{n_pos_merged} position-fallback, "
        f"{n_ncbi_only} NCBI-only, {n_cyan_only} Cyanorak-only, "
        f"{len(merge_df)} total"
    )

    return merge_df.rename(columns=_FINAL_MERGE_RENAME)


# ── top-level entry point ─────────────────────────────────────────────────────

def build_gene_mapping(
    ncbi_gff_file: str,
    ncbi_gbff_file: str | None = None,
    cyan_gff_file: str | None = None,
    cyan_gbk_file: str | None = None,
) -> pd.DataFrame:
    """Build the gene mapping DataFrame from GFF/GBFF and optional Cyanorak files.

    Args:
        ncbi_gff_file: Path to NCBI genomic.gff file (required).
        ncbi_gbff_file: Path to NCBI genomic.gbff file for EC number extraction.
        cyan_gff_file: Path to Cyanorak GFF file (enables merged mode).
        cyan_gbk_file: Path to Cyanorak GBK file (enables merged mode).

    Returns:
        DataFrame with one row per gene and all available annotation columns.
    """
    has_cyanorak = bool(cyan_gff_file and cyan_gbk_file)
    if has_cyanorak:
        df = load_gff_from_ncbi_and_cyanorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
    else:
        logger.info("No Cyanorak data available; loading NCBI GFF only")
        df = load_gff_from_ncbi_only(ncbi_gff_file)

    # Merge EC numbers from NCBI GBFF if available
    if ncbi_gbff_file and os.path.exists(ncbi_gbff_file):
        ec_map = _get_ec_numbers_from_gbff(ncbi_gbff_file)
        if ec_map and 'locus_tag_ncbi' in df.columns:
            df['ec_numbers'] = df['locus_tag_ncbi'].apply(
                lambda lt: ','.join(ec_map[lt]) if pd.notna(lt) and lt in ec_map else None
            )
            n_with_ec = df['ec_numbers'].notna().sum()
            logger.info(f"EC numbers merged: {n_with_ec} genes have at least one EC number")

    return df
