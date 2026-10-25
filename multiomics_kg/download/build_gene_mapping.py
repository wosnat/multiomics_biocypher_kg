"""Build gene_mapping.csv from NCBI GFF/GBFF and Cyanorak GFF/GBK files.

Provides standalone functions for loading and merging gene annotation data.
Used by both CyanorakNcbi adapter and the download pipeline (download_genome_data.py).
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

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

    # Older NCBI GFFs (annotwriter, e.g. GCF_000218705.1 HNLC1, GCF_000179235.1
    # CB0101: ID=gene0 / ID=cds0) carry locus_tag on the gene row only. Without
    # this fallback every NCBI gene was dropped for lacking locus_tag_ncbi and
    # the strain came out Cyanorak-only.
    if 'locus_tag_gene' in df.columns:
        if 'locus_tag_cds' in df.columns:
            df['locus_tag_cds'] = df['locus_tag_cds'].fillna(df['locus_tag_gene'])
        else:
            df['locus_tag_cds'] = df['locus_tag_gene']

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
    (a contig backed by fewer than 3 matched genes is used but warned about).
    Contigs without a derivable offset are left out — with a
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
            # No paired gene on this contig: the offset is unknown. (Assuming 0
            # on a single-sequence assembly paired ~330-400 genes per strain by
            # coordinate coincidence when Cyanorak's record is laid out in
            # another frame — KORDI-49/52/100, CC9616; 2026-10-26.)
            n_unmatched = int(((seqids == sid) & ncbi_sourced['ID'].isna()).sum())
            logger.warning(
                f"Position fallback: contig {sid!r} has no paired gene to derive a "
                f"Cyanorak offset from — its {n_unmatched} unmatched NCBI genes are "
                f"not considered"
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

        _absorb_cyanorak_row(ncbi_sourced, ncbi_idx, cyan_only, cyan_lt, 'position_merge')
        consumed_cyan.add(cyan_lt)
        n_merged += 1

    return n_merged, consumed_cyan


def _absorb_cyanorak_row(ncbi_sourced: pd.DataFrame, ncbi_idx, cyan_only: pd.DataFrame,
                         cyan_lt: str, how: str) -> None:
    """Copy one Cyanorak call's columns into an unmatched NCBI row (in place).

    Shared by the position and sequence fallbacks: copies ``_CYAN_COLUMNS_TO_COPY``,
    appends the Cyanorak locus_tag to ``old_locus_tags`` (papers may use that
    form, so it must stay resolvable by gene_id_mapping) and records
    ``position_merge_note = '<how>:<cyanorak lt>→<target>'``.
    """
    cyan_row = cyan_only[cyan_only['locus_tag'] == cyan_lt].iloc[0]
    for col in _CYAN_COLUMNS_TO_COPY:
        if col in ncbi_sourced.columns and col in cyan_only.columns:
            ncbi_sourced.at[ncbi_idx, col] = cyan_row[col]

    if 'old_locus_tags' not in ncbi_sourced.columns:
        ncbi_sourced['old_locus_tags'] = pd.Series(pd.NA, index=ncbi_sourced.index, dtype=object)
    elif ncbi_sourced['old_locus_tags'].dtype != object:
        ncbi_sourced['old_locus_tags'] = ncbi_sourced['old_locus_tags'].astype(object)
    existing = ncbi_sourced.at[ncbi_idx, 'old_locus_tags']
    if pd.notna(existing):
        if cyan_lt not in str(existing):
            ncbi_sourced.at[ncbi_idx, 'old_locus_tags'] = f"{existing},{cyan_lt}"
    else:
        ncbi_sourced.at[ncbi_idx, 'old_locus_tags'] = cyan_lt

    # locus_tag is still NaN here for NCBI genes without old_locus_tag
    # (filled from locus_tag_ncbi later), so name the gene the same way
    target = ncbi_sourced.at[ncbi_idx, 'locus_tag']
    if pd.isna(target):
        target = ncbi_sourced.at[ncbi_idx, 'locus_tag_ncbi']
    ncbi_sourced.at[ncbi_idx, 'position_merge_note'] = f"{how}:{cyan_lt}→{target}"


# ── protein-sequence fallback (diamond bridge) ───────────────────────────────

_PROTEIN_MATCH_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "map_img_to_ncbi_proteins.py"
_SEQ_FALLBACK_MIN_UNPAIRED = 50   # run only when the earlier stages left this many on both sides

# QC thresholds (report-only, ``[QC]`` warnings — see _merge_qc)
QC_MAX_DUPLICATION_RATIO = 1.15
QC_MIN_SYNTENY_RATE = 0.90
QC_MIN_DIAMOND_PIDENT_P5 = 90.0


def _clean_translation(seq) -> str | None:
    """Cyanorak ``translation=`` → comparable protein: stop glyphs (``*#+``)
    stripped, uppercase, first residue set to M (Cyanorak writes the start codon
    as its own residue; NCBI as M) — same normalisation as Gene.sequence; any
    other non-letter inside the sequence becomes X."""
    if not isinstance(seq, str):
        return None
    # internal stop/ambiguity glyphs (``+``, ``*``) make diamond reject the whole FASTA
    s = re.sub(r"[^A-Z]", "X", seq.strip().upper().rstrip("*#+"))
    return ("M" + s[1:]) if len(s) > 1 else None


def _read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _sequence_fallback_merge(
    ncbi_sourced: pd.DataFrame,
    cyan_only: pd.DataFrame,
    protein_faa: str,
    cyan_translations: dict[str, str],
) -> tuple[int, set[str], dict]:
    """Pair still-unpaired NCBI and Cyanorak calls by protein sequence.

    Both earlier stages assume one assembly in one coordinate frame. For some
    Cyanorak genomes neither holds (Ford 2021 batch, 2026-10-25: WH8020,
    KORDI-49/52/100, CC9616, PCC7001, CB0101, CB0205 — 0 locus-tag matches,
    0-13 % position merges), leaving each gene as an NCBI-only plus a
    Cyanorak-only row. This stage runs ``scripts/map_img_to_ncbi_proteins.py``
    (exact → subsequence → diamond blastp; the EZ55 / MIT1002 / MIT1327 bridge)
    with the unpaired Cyanorak translations as source and the unpaired NCBI
    proteins as targets, then merges **1:1 only**: an NCBI gene claimed by
    several Cyanorak calls is skipped, and a WP accession shared by several
    unpaired NCBI genes (identical paralogs) is never a target.

    Runs only when both sides keep >= ``_SEQ_FALLBACK_MIN_UNPAIRED`` unpaired
    genes, so normal strains (a few dozen unpaired) are untouched.

    Modifies ncbi_sourced in place. Returns (n_merged, consumed_cyanorak_lts, qc).
    """
    qc: dict = {"ran": False}
    if "protein_id" not in ncbi_sourced.columns or not cyan_translations:
        qc["reason"] = "no protein_id column or no Cyanorak translations"
        return 0, set(), qc
    unmatched = ncbi_sourced[ncbi_sourced["ID"].isna() & ncbi_sourced["protein_id"].notna()]
    if len(unmatched) < _SEQ_FALLBACK_MIN_UNPAIRED or len(cyan_only) < _SEQ_FALLBACK_MIN_UNPAIRED:
        qc["reason"] = f"unpaired below threshold ({len(unmatched)} NCBI / {len(cyan_only)} Cyanorak)"
        return 0, set(), qc
    if not _PROTEIN_MATCH_SCRIPT.exists():
        logger.warning(f"Sequence fallback: {_PROTEIN_MATCH_SCRIPT} missing — skipped")
        qc["reason"] = "protein match script missing"
        return 0, set(), qc

    pid_counts = unmatched["protein_id"].value_counts()
    shared = set(pid_counts[pid_counts > 1].index)
    targets = unmatched[~unmatched["protein_id"].isin(shared)]
    lt_to_idx = dict(zip(targets["locus_tag_ncbi"], targets.index))
    ncbi_prot = _read_fasta(protein_faa)

    with tempfile.TemporaryDirectory(prefix="gene_mapping_seq_") as tmp:
        src, tgt = os.path.join(tmp, "cyanorak.faa"), os.path.join(tmp, "ncbi.faa")
        gm, out = os.path.join(tmp, "gm.csv"), os.path.join(tmp, "bridge.csv")
        n_src = 0
        with open(src, "w") as fh:
            for lt in cyan_only["locus_tag"]:
                s = _clean_translation(cyan_translations.get(lt))
                if s:
                    fh.write(f">{lt}\n{s}\n")
                    n_src += 1
        with open(tgt, "w") as fh, open(gm, "w", newline="") as gfh:
            w = csv.writer(gfh)
            w.writerow(["protein_id", "locus_tag"])
            for lt, pid in zip(targets["locus_tag_ncbi"], targets["protein_id"]):
                s = ncbi_prot.get(pid)
                if s:
                    fh.write(f">{pid}\n{s}\n")
                    w.writerow([pid, lt])
        proc = subprocess.run(
            [sys.executable, str(_PROTEIN_MATCH_SCRIPT), "--img-faa", src, "--ncbi-faa", tgt,
             "--gene-mapping", gm, "--output", out, "--source-id-col", "cyanorak_locus_tag",
             "--with-method"],
            capture_output=True, text=True,
        )
        # exit 2 = "<50 % of the source mapped" warning; the CSV is still valid
        if proc.returncode not in (0, 2) or not os.path.exists(out):
            logger.warning(f"Sequence fallback: protein match failed (rc={proc.returncode}): "
                           f"{proc.stderr.strip()[-300:]}")
            qc["reason"] = f"protein match failed rc={proc.returncode}"
            return 0, set(), qc
        with open(out) as fh:
            bridge = list(csv.DictReader(fh))

    by_ncbi: dict[str, list[dict]] = {}
    for row in bridge:
        by_ncbi.setdefault(row["locus_tag"], []).append(row)
    by_cyan: dict[str, int] = {}
    for row in bridge:
        by_cyan[row["cyanorak_locus_tag"]] = by_cyan.get(row["cyanorak_locus_tag"], 0) + 1
    candidates: dict = {}
    n_conflict = 0
    for ncbi_lt, rows in by_ncbi.items():
        if len(rows) > 1 or ncbi_lt not in lt_to_idx or by_cyan[rows[0]["cyanorak_locus_tag"]] > 1:
            n_conflict += 1
            continue
        candidates[lt_to_idx[ncbi_lt]] = rows[0]

    syntenic = _syntenic_candidates(ncbi_sourced, cyan_only, candidates)
    n_non_syntenic = len(candidates) - len(syntenic)

    consumed: set[str] = set()
    methods: dict[str, int] = {}
    pidents: list[float] = []
    for idx, row in candidates.items():
        if idx not in syntenic:
            continue
        cyan_lt, how = row["cyanorak_locus_tag"], row["match_method"]
        _absorb_cyanorak_row(ncbi_sourced, idx, cyan_only, cyan_lt, f"sequence_merge_{how}")
        consumed.add(cyan_lt)
        methods[how] = methods.get(how, 0) + 1
        if how == "diamond":
            pidents.append(float(row["pident"]))
    qc.update({
        "ran": True,
        "source_cyanorak": n_src,
        "target_ncbi": len(targets),
        "excluded_shared_wp": int(len(unmatched) - len(targets)),
        "merged": len(consumed),
        "merged_by_method": methods,
        "skipped_conflicts": n_conflict,
        "rejected_non_syntenic": n_non_syntenic,
        "diamond_pident_min": min(pidents) if pidents else None,
        "diamond_pident_p5": float(np.percentile(pidents, 5)) if pidents else None,
        "diamond_pident_median": float(np.median(pidents)) if pidents else None,
    })
    return len(consumed), consumed, qc


def _syntenic_candidates(ncbi_sourced: pd.DataFrame, cyan_only: pd.DataFrame,
                         candidates: dict) -> set:
    """Keep only sequence-match candidates that sit where their neighbours say.

    A true pair maps one locus onto one locus, so its 3'-end offset between the
    NCBI and Cyanorak frames (difference; sum for a reverse-complemented record)
    equals that of one of the two nearest paired genes on either side on the
    same NCBI contig — already-paired genes and other candidates both count, so a whole
    re-laid-out genome (CB0101: 2,749 sequence pairs, all mutually consistent)
    passes while an isolated hit on a paralog elsewhere in the genome is
    rejected. Measured 2026-10-25: without this filter the stage made 1-4
    non-syntenic merges in ~23 normal strains (MIT9313, CC9311, …), every one a
    paralog with its own locus.
    """
    if not candidates:
        return set()
    cyan_by_lt = cyan_only.set_index("locus_tag")
    recs = []   # (seqid, ncbi_start, diff, sum, candidate_idx | None)
    paired = ncbi_sourced[ncbi_sourced["ID"].notna()]
    need = ["start_ncbi", "end_ncbi", "strand_ncbi", "start_cyanorak", "end_cyanorak", "strand_cyanorak"]
    if set(need) <= set(paired.columns):
        p = paired.dropna(subset=need)
        n3 = _three_prime_end(p["start_ncbi"], p["end_ncbi"], p["strand_ncbi"])
        c3 = _three_prime_end(p["start_cyanorak"], p["end_cyanorak"], p["strand_cyanorak"])
        for sid, st, d, s in zip(_seqids(p), p["start_ncbi"], c3 - n3, c3 + n3):
            recs.append((sid, float(st), d, s, None))
    seqids = _seqids(ncbi_sourced)
    for idx, row in candidates.items():
        lt = row["cyanorak_locus_tag"]
        if lt not in cyan_by_lt.index:
            continue
        c = cyan_by_lt.loc[lt]
        if isinstance(c, pd.DataFrame):
            c = c.iloc[0]
        n = ncbi_sourced.loc[idx]
        if pd.isna(n.get("start_ncbi")) or pd.isna(c.get("start_cyanorak")):
            continue
        n3 = float(_three_prime_end(n["start_ncbi"], n["end_ncbi"], n["strand_ncbi"]))
        c3 = float(_three_prime_end(c["start_cyanorak"], c["end_cyanorak"], c["strand_cyanorak"]))
        recs.append((seqids.loc[idx], float(n["start_ncbi"]), c3 - n3, c3 + n3, idx))
    keep: set = set()
    recs.sort(key=lambda r: (r[0], r[1]))
    for i, (sid, _st, d, s, idx) in enumerate(recs):
        if idx is None:
            continue
        for j in (i - 2, i - 1, i + 1, i + 2):
            if 0 <= j < len(recs) and recs[j][0] == sid and (recs[j][2] == d or recs[j][3] == s):
                keep.add(idx)
                break
    return keep


QC_MAX_PAIR_SEQ_DISAGREE = 0.05
_PAIR_TAIL_AA = 20


def _pair_sequence_check(merge_df: pd.DataFrame, ncbi_proteins: dict, cyan_translations: dict) -> dict:
    """Do paired NCBI / Cyanorak calls encode the same protein?

    For every paired row with both sequences available, the two proteins agree
    when their C-terminal 20 aa are identical (starts are legitimately re-called,
    stops are not). Reported per pairing stage; a position stage that paired
    genes by coordinate coincidence shows up as a high disagreement rate.
    """
    out: dict = {}
    if not ncbi_proteins or not cyan_translations or "ID" not in merge_df.columns:
        return out
    note = merge_df.get("position_merge_note", pd.Series(index=merge_df.index, dtype=object)).astype(str)
    for idx, row in merge_df[merge_df["ID"].notna() & merge_df.get("protein_id", pd.Series(dtype=object)).notna()].iterrows():
        n = note.loc[idx]
        if n.startswith("sequence_merge"):
            stage, cyan_lt = "sequence", n.split(":", 1)[1].split("→")[0]
        elif n.startswith("position_merge"):
            stage, cyan_lt = "position", n.split(":", 1)[1].split("→")[0]
        else:
            stage, cyan_lt = "locus_tag", row.get("locus_tag")
        a = _clean_translation(ncbi_proteins.get(row["protein_id"]))
        b = _clean_translation(cyan_translations.get(cyan_lt))
        if not a or not b:
            continue
        rec = out.setdefault(stage, {"checked": 0, "disagree": 0})
        rec["checked"] += 1
        if a[-_PAIR_TAIL_AA:] != b[-_PAIR_TAIL_AA:]:
            rec["disagree"] += 1
    for rec in out.values():
        rec["disagree_rate"] = round(rec["disagree"] / rec["checked"], 3) if rec["checked"] else None
    return out


def _merge_qc(merge_df: pd.DataFrame, n_ncbi: int, n_cyan: int, stages: dict, seq_qc: dict,
              ncbi_proteins: dict | None = None, cyan_translations: dict | None = None) -> dict:
    """Per-strain merge QC (written to ``gene_mapping_merge_qc.json`` by step 0.5).

    - duplication_ratio = genes / max(NCBI genes, Cyanorak genes): ~1.0-1.1 when
      the two annotations were paired; ~1.8 when every gene is duplicated
      (WH8020 before the sequence fallback: 5,501 / 3,115).
    - synteny_rate over sequence-merged pairs: a pair whose 3'-end offset between
      the two coordinate frames (difference, or sum for a reverse-complemented
      record) equals that of its nearest paired neighbour on the same contig is
      syntenic; an isolated offset flags a likely paralog mix-up.
    Flags breach QC_MAX_DUPLICATION_RATIO / QC_MIN_SYNTENY_RATE /
    QC_MIN_DIAMOND_PIDENT_P5 and are also logged as ``[QC]`` warnings.
    """
    qc = {"genes": int(len(merge_df)), "ncbi_genes": int(n_ncbi), "cyanorak_genes": int(n_cyan),
          "stages": stages, "sequence_fallback": seq_qc}
    qc["duplication_ratio"] = round(len(merge_df) / max(n_ncbi, n_cyan, 1), 3)

    cols = ["start_ncbi", "end_ncbi", "start_cyanorak", "end_cyanorak"]
    paired = merge_df[merge_df[cols].notna().all(axis=1)].copy() if set(cols) <= set(merge_df.columns) else merge_df.iloc[0:0]
    synteny = None
    if not paired.empty and "position_merge_note" in paired.columns:
        n3 = _three_prime_end(paired["start_ncbi"], paired["end_ncbi"], paired["strand_ncbi"])
        c3 = _three_prime_end(paired["start_cyanorak"], paired["end_cyanorak"], paired["strand_cyanorak"])
        paired["_diff"], paired["_sum"] = c3 - n3, c3 + n3
        paired["_seqid"] = _seqids(paired)
        paired = paired.sort_values(["_seqid", "start_ncbi"])
        is_seq = paired["position_merge_note"].astype(str).str.startswith("sequence_merge")
        n_seq = int(is_seq.sum())
        if n_seq:
            ok = 0
            for sid, grp in paired.groupby("_seqid", sort=False):
                d, s, flag = grp["_diff"].values, grp["_sum"].values, is_seq.loc[grp.index].values
                for i in np.where(flag)[0]:
                    nb = [j for j in (i - 1, i + 1) if 0 <= j < len(grp)]
                    if any(d[j] == d[i] or s[j] == s[i] for j in nb):
                        ok += 1
            synteny = round(ok / n_seq, 3)
    qc["synteny_rate"] = synteny
    qc["pair_sequence_check"] = _pair_sequence_check(merge_df, ncbi_proteins or {}, cyan_translations or {})

    flags = []
    for stage, rec in qc["pair_sequence_check"].items():
        if rec["checked"] >= 20 and rec["disagree_rate"] > QC_MAX_PAIR_SEQ_DISAGREE:
            flags.append(f"{stage} pairs: {rec['disagree_rate']:.1%} encode different proteins "
                         f"(> {QC_MAX_PAIR_SEQ_DISAGREE:.0%})")
    if qc["duplication_ratio"] > QC_MAX_DUPLICATION_RATIO:
        flags.append(f"duplication_ratio {qc['duplication_ratio']} > {QC_MAX_DUPLICATION_RATIO}")
    if synteny is not None and synteny < QC_MIN_SYNTENY_RATE:
        flags.append(f"synteny_rate {synteny} < {QC_MIN_SYNTENY_RATE}")
    p5 = seq_qc.get("diamond_pident_p5")
    if p5 is not None and p5 < QC_MIN_DIAMOND_PIDENT_P5:
        flags.append(f"diamond_pident_p5 {p5:.1f} < {QC_MIN_DIAMOND_PIDENT_P5}")
    qc["flags"] = flags
    for f in flags:
        logger.warning(f"[QC] gene merge: {f}")
    return qc


# ── Cyanorak-only coordinate fallback ────────────────────────────────────────

def _ncbi_contig_lengths(ncbi_gff_file: str | None) -> dict[str, int]:
    """Contig lengths from the NCBI GFF ``##sequence-region <seqid> 1 <end>``
    directives. ``{}`` when the file is absent (unit tests patch the GFF loader
    and pass a placeholder path)."""
    lengths: dict[str, int] = {}
    if not ncbi_gff_file or not os.path.exists(ncbi_gff_file):
        return lengths
    pattern = re.compile(r'##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)')
    with open(ncbi_gff_file) as fh:
        for line in fh:
            if not line.startswith('#'):
                continue
            m = pattern.match(line)
            if m:
                lengths[m.group(1)] = int(m.group(3))
    return lengths


def _cyanorak_tiles(offsets: dict, lengths: dict[str, int]) -> list[tuple]:
    """Where each NCBI contig sits inside Cyanorak's single concatenated record.

    Cyanorak lays a draft genome out as the NCBI contigs in accession order with
    a 50 bp spacer between neighbours (verified on PAC1 and SB: the per-contig
    ``_contig_offsets`` plus the GFF contig lengths tile the record with a
    50 bp gap at every boundary, and the last tile ends exactly at the record
    length). Returns ``[(first_base, last_base, seqid), ...]`` in Cyanorak's
    1-based coordinates, sorted by offset. A single contig whose length is
    unknown gets an open-ended tile; a multi-contig assembly with a missing
    length skips that contig with a warning.
    """
    tiles: list[tuple] = []
    for sid, offset in sorted(offsets.items(), key=lambda kv: kv[1]):
        length = lengths.get(sid)
        if length is None:
            if len(offsets) == 1:
                tiles.append((offset + 1, float('inf'), sid))
            else:
                logger.warning(
                    f"Coordinate fallback: no ##sequence-region length for contig "
                    f"{sid!r} — Cyanorak-only genes on it stay unplaced"
                )
            continue
        tiles.append((offset + 1, offset + length, sid))
    return tiles


def _fill_cyanorak_only_coordinates(cyan_only: pd.DataFrame, tiles: list[tuple]) -> int:
    """Give Cyanorak-only genes NCBI-space coordinates by inverting the tiling.

    For each row with Cyanorak coordinates but no NCBI call, find the tile that
    contains BOTH ends; set ``seqid`` to that contig, shift ``start``/``end``
    into contig-relative coordinates, copy the strand, and mark
    ``coordinate_source = 'cyanorak'``. A gene whose interval crosses a tile
    boundary spans the inter-contig spacer — an assembly artefact — and is
    left unplaced with a warning naming it. Modifies ``cyan_only`` in place and
    returns the number of rows placed.
    """
    if cyan_only.empty or not tiles:
        return 0
    for col in ('seqid', 'strand_ncbi', 'coordinate_source'):
        if col not in cyan_only.columns or cyan_only[col].isna().all():
            cyan_only[col] = pd.Series(pd.NA, index=cyan_only.index, dtype=object)
    for col in ('start_ncbi', 'end_ncbi'):
        if col not in cyan_only.columns:
            cyan_only[col] = np.nan
    unplaced: list[str] = []
    n_placed = 0
    for idx, row in cyan_only.iterrows():
        c_start, c_end = row.get('start_cyanorak'), row.get('end_cyanorak')
        if pd.isna(c_start) or pd.isna(c_end):
            continue
        tile = next((t for t in tiles if t[0] <= c_start and c_end <= t[1]), None)
        if tile is None:
            unplaced.append(str(row.get('locus_tag')))
            continue
        first_base, _last_base, sid = tile
        cyan_only.at[idx, 'seqid'] = sid
        cyan_only.at[idx, 'start_ncbi'] = float(c_start) - (first_base - 1)
        cyan_only.at[idx, 'end_ncbi'] = float(c_end) - (first_base - 1)
        cyan_only.at[idx, 'strand_ncbi'] = row.get('strand_cyanorak')
        cyan_only.at[idx, 'coordinate_source'] = 'cyanorak'
        n_placed += 1
    if unplaced:
        logger.warning(
            f"Coordinate fallback: {len(unplaced)} Cyanorak-only gene(s) cross a "
            f"contig boundary or fall outside every contig — left unplaced: "
            f"{unplaced[:10]}{' …' if len(unplaced) > 10 else ''}"
        )
    return n_placed


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
    ncbi_df['coordinate_source'] = 'ncbi'

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

    # Cyanorak protein translations for the sequence fallback (not carried on)
    cyan_translations = (
        dict(zip(cyan_df['locus_tag'], cyan_df['translation']))
        if 'translation' in cyan_df.columns else {}
    )

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

    # --- Position + protein-sequence fallbacks ---
    # Position pairing needs the coordinate offset between the two frames, which
    # is measured on already-paired genes. With no locus_tag pair at all the
    # frame is unknown, so the sequence fallback runs first and the position
    # fallback then measures its offsets on the sequence-proven pairs.
    n_ncbi_matched = int(ncbi_sourced['ID'].notna().sum()) if not ncbi_sourced.empty else 0
    n_pos_merged = n_seq_merged = 0
    seq_qc: dict = {"ran": False}
    protein_faa = os.path.join(os.path.dirname(ncbi_gff_file), "protein.faa")
    ncbi_proteins = _read_fasta(protein_faa) if os.path.exists(protein_faa) else {}

    def _run_position():
        nonlocal cyan_only, n_pos_merged
        if not ncbi_sourced.empty and not cyan_only.empty:
            n, consumed = _position_fallback_merge(ncbi_sourced, cyan_only)
            n_pos_merged += n
            if consumed:
                cyan_only = cyan_only[~cyan_only['locus_tag'].isin(consumed)]

    def _run_sequence():
        nonlocal cyan_only, n_seq_merged, seq_qc
        if not ncbi_sourced.empty and not cyan_only.empty and ncbi_proteins:
            n, consumed, seq_qc = _sequence_fallback_merge(
                ncbi_sourced, cyan_only, protein_faa, cyan_translations,
            )
            n_seq_merged += n
            if consumed:
                cyan_only = cyan_only[~cyan_only['locus_tag'].isin(consumed)]

    if n_ncbi_matched == 0:
        _run_sequence()
        _run_position()
    else:
        _run_position()
        _run_sequence()
    seq_qc["order"] = "sequence_first" if n_ncbi_matched == 0 else "position_first"

    n_ncbi_only = len(ncbi_sourced) - n_ncbi_matched - n_pos_merged - n_seq_merged
    n_cyan_only = len(cyan_only)

    # --- Coordinates for Cyanorak-only genes ---
    # PGAP never called these, so their NCBI coordinate columns are blank even
    # though Cyanorak places every one. Invert Cyanorak's contig tiling to give
    # them a contig + contig-relative coordinates (identity on closed genomes).
    n_coord_placed = 0
    if not ncbi_sourced.empty:
        ncbi_sourced['coordinate_source'] = 'ncbi'
        if not cyan_only.empty:
            tiles = _cyanorak_tiles(
                _contig_offsets(ncbi_sourced), _ncbi_contig_lengths(ncbi_gff_file),
            )
            n_coord_placed = _fill_cyanorak_only_coordinates(cyan_only, tiles)

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
        f"{n_pos_merged} position-fallback, {n_seq_merged} sequence-fallback, "
        f"{n_ncbi_only} NCBI-only, {n_cyan_only} Cyanorak-only "
        f"({n_coord_placed} placed on NCBI contigs via the Cyanorak tiling), "
        f"{len(merge_df)} total"
    )

    qc = _merge_qc(
        merge_df, ncbi_df.shape[0], cyan_df.shape[0],
        {"locus_tag": n_ncbi_matched, "position": n_pos_merged, "sequence": n_seq_merged,
         "ncbi_only": n_ncbi_only, "cyanorak_only": n_cyan_only},
        seq_qc, ncbi_proteins, cyan_translations,
    )
    out = merge_df.rename(columns=_FINAL_MERGE_RENAME)
    out.attrs["merge_qc"] = qc
    return out


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
