#!/usr/bin/env python3
"""
Map researcher/IMG gene IDs to NCBI locus tags via protein sequence matching.

Bridges gene identifiers between two independent genome annotations of the
same organism by aligning their predicted protein sequences. Designed for
cases where no shared identifier exists between the annotations (different
gene-calling pipelines, different assemblies).

Use cases:
  - Alteromonas macleodii EZ55: AEZ55_NNNN (researcher) → EZ55_NNNNN (NCBI)
  - Alteromonas macleodii MIT1002: fig|226.6.peg.N (RAST) → MIT1002_NNNNN (NCBI)

Method:
  1. Parse source protein FASTA → {gene_id: sequence}
  2. Parse NCBI protein FASTA → {protein_accession: sequence}
  3. Load gene_mapping.csv → {protein_accession: locus_tag}
  4. Phase 1: exact protein sequence match
  5. Phase 2: subsequence matching (annotation boundary differences)
  6. Phase 3: Diamond blastp (near-identical, handles draft fragmentation)
  7. Output: id_translation CSV for paperconfig

Usage:
  uv run python scripts/map_img_to_ncbi_proteins.py \
    --img-faa <source_protein.faa> \
    --ncbi-faa <cache/.../protein.faa> \
    --gene-mapping <cache/.../gene_mapping.csv> \
    --output <id_translation.csv> \
    --source-id-col <column_name> \
    [--img-gff <gff_for_header_remapping>]

The optional --img-gff flag is used when FASTA headers contain numeric
OIDs instead of gene IDs; the GFF provides the OID→gene_id mapping.
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_fasta(path: Path) -> dict[str, str]:
    """Parse FASTA file → {header_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                # Take first whitespace-delimited token after >
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
    return sequences


def parse_img_gff_gene_ids(path: Path) -> dict[str, str]:
    """Parse Hennon-style IMG GFF → {AEZ55_NNNN: AEZ55_NNNN}.

    Also builds numeric_oid→gene_id mapping if needed.
    The GFF has lines like: scf... gene_id AEZ55_0001
    """
    gene_ids = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                # Extract gene_id from attributes
                attrs = parts[8]
                match = re.search(r"gene_id\s+(\S+)", attrs)
                if match:
                    gid = match.group(1)
                    gene_ids[gid] = gid
    return gene_ids


def load_gene_mapping(path: Path) -> dict[str, str]:
    """Load gene_mapping.csv → {protein_id: locus_tag}."""
    protein_to_lt = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            lt = row.get("locus_tag", "").strip()
            pid = row.get("protein_id", "").strip()
            if lt and pid:
                protein_to_lt[pid] = lt
    return protein_to_lt


def map_by_exact_match(
    img_seqs: dict[str, str],
    ncbi_seqs: dict[str, str],
    protein_to_lt: dict[str, str],
) -> tuple[dict[str, str], set[str], set[str]]:
    """Map IMG gene IDs to locus tags via exact protein sequence match.

    Returns (mapping, unmatched_img, used_ncbi).
    """
    # Build reverse index: sequence → list of NCBI protein IDs
    seq_to_ncbi = defaultdict(list)
    for pid, seq in ncbi_seqs.items():
        seq_to_ncbi[seq].append(pid)

    mapping = {}  # img_gene_id → locus_tag
    unmatched = set()
    used_ncbi = set()
    ambiguous = []

    for img_id, img_seq in img_seqs.items():
        ncbi_hits = seq_to_ncbi.get(img_seq, [])
        if len(ncbi_hits) == 1:
            pid = ncbi_hits[0]
            lt = protein_to_lt.get(pid)
            if lt:
                mapping[img_id] = lt
                used_ncbi.add(pid)
            else:
                # Protein exists but not in gene_mapping (no locus_tag link)
                unmatched.add(img_id)
        elif len(ncbi_hits) > 1:
            # Multiple NCBI proteins with identical sequence (paralogs)
            # Try to find unique locus_tag among them
            lts = set()
            for pid in ncbi_hits:
                lt = protein_to_lt.get(pid)
                if lt:
                    lts.add((pid, lt))
            if len(lts) == 1:
                pid, lt = lts.pop()
                mapping[img_id] = lt
                used_ncbi.add(pid)
            else:
                ambiguous.append((img_id, ncbi_hits))
                unmatched.add(img_id)
        else:
            unmatched.add(img_id)

    if ambiguous:
        print(f"  Ambiguous matches (identical seq, multiple locus_tags): {len(ambiguous)}")
        for img_id, hits in ambiguous[:5]:
            lts = [protein_to_lt.get(p, "?") for p in hits]
            print(f"    {img_id} → {hits} → {lts}")

    return mapping, unmatched, used_ncbi


def map_by_subsequence(
    img_seqs: dict[str, str],
    ncbi_seqs: dict[str, str],
    protein_to_lt: dict[str, str],
    unmatched_img: set[str],
    used_ncbi: set[str],
    min_overlap_frac: float = 0.95,
) -> dict[str, str]:
    """Try subsequence matching for proteins with slightly different boundaries.

    Checks if one sequence is a subsequence of the other (annotation boundary
    differences — same CDS, different start/stop prediction).
    """
    extra_mapping = {}
    still_unmatched = set()

    # Build list of unused NCBI sequences
    unused_ncbi = {
        pid: seq for pid, seq in ncbi_seqs.items() if pid not in used_ncbi
    }

    for img_id in unmatched_img:
        img_seq = img_seqs[img_id]
        if len(img_seq) < 30:
            continue  # Skip very short sequences

        best_match = None
        best_overlap = 0

        for pid, ncbi_seq in unused_ncbi.items():
            lt = protein_to_lt.get(pid)
            if not lt:
                continue

            # Check containment in either direction
            shorter, longer = (img_seq, ncbi_seq) if len(img_seq) <= len(ncbi_seq) else (ncbi_seq, img_seq)
            if shorter in longer:
                overlap_frac = len(shorter) / max(len(img_seq), len(ncbi_seq))
                if overlap_frac >= min_overlap_frac and overlap_frac > best_overlap:
                    best_match = (pid, lt)
                    best_overlap = overlap_frac

        if best_match:
            pid, lt = best_match
            extra_mapping[img_id] = lt
        else:
            still_unmatched.add(img_id)

    return extra_mapping


def map_by_diamond(
    img_seqs: dict[str, str],
    ncbi_seqs: dict[str, str],
    protein_to_lt: dict[str, str],
    unmatched_img: set[str],
    used_ncbi: set[str],
    min_identity: float = 80.0,
    min_query_coverage: float = 60.0,
) -> tuple[dict[str, str], list[tuple[str, str, int]]]:
    """Use Diamond blastp to match remaining proteins by sequence similarity.

    Catches near-identical proteins with different start codons, a few
    substitutions, or minor truncations that exact/substring matching misses.

    Draft genome annotations often split one real gene into multiple smaller
    ORFs (due to frameshifts or assembly gaps). When multiple query fragments
    hit the same NCBI target, only the longest fragment is kept — its
    expression values best represent the gene since it captured the most
    RNA-seq reads. Shorter fragments are recorded as discarded.

    Only query coverage is required (not subject coverage), because a valid
    fragment will match well across its own length but may cover only part
    of the larger NCBI protein.

    Returns (mapping, discarded) where discarded is a list of
    (img_id, locus_tag, query_len) tuples for fragments that lost to a
    longer sibling.
    """
    if not shutil.which("diamond"):
        print("  WARNING: diamond not found in PATH — skipping similarity search")
        return {}, []

    with tempfile.TemporaryDirectory(prefix="diamond_protein_map_") as tmpdir:
        tmp = Path(tmpdir)
        query_faa = tmp / "query.faa"
        db_faa = tmp / "db.faa"
        db_path = tmp / "db"
        out_tsv = tmp / "hits.tsv"

        # Write unmatched researcher proteins as query
        with open(query_faa, "w") as f:
            for img_id in sorted(unmatched_img):
                seq = img_seqs.get(img_id, "")
                if len(seq) >= 20:
                    f.write(f">{img_id}\n{seq}\n")

        # Write ALL NCBI proteins that have a locus_tag as database
        # (include already-matched ones — paralogs/transposases may map many:1)
        with open(db_faa, "w") as f:
            for pid, seq in ncbi_seqs.items():
                if pid in protein_to_lt:
                    f.write(f">{pid}\n{seq}\n")

        # Check if query file has any sequences
        if query_faa.stat().st_size == 0:
            return {}, []

        # Build Diamond DB
        try:
            subprocess.run(
                ["diamond", "makedb", "--in", str(db_faa), "-d", str(db_path)],
                capture_output=True,
                check=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            print("  WARNING: diamond makedb timed out after 120s")
            return {}, []

        # Run Diamond blastp — no subject-cover filter (fragments are expected)
        try:
            subprocess.run(
                [
                    "diamond", "blastp",
                    "-q", str(query_faa),
                    "-d", str(db_path),
                    "-o", str(out_tsv),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "evalue", "bitscore",
                    "--id", str(min_identity),
                    "--query-cover", str(min_query_coverage),
                    "--max-target-seqs", "1",
                    "--very-sensitive",
                ],
                capture_output=True,
                check=True,
                timeout=600,
            )
        except subprocess.TimeoutExpired:
            print("  WARNING: diamond blastp timed out after 600s")
            return {}, []

        # Collect all hits: img_id → (locus_tag, query_length)
        all_hits: dict[str, tuple[str, int]] = {}
        with open(out_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                img_id, pid = parts[0], parts[1]
                lt = protein_to_lt.get(pid)
                if lt and img_id not in all_hits:
                    qlen = len(img_seqs.get(img_id, ""))
                    all_hits[img_id] = (lt, qlen)

    # Deduplicate: when multiple fragments hit the same locus_tag,
    # keep only the longest fragment (most RNA-seq reads).
    # Group by locus_tag
    lt_to_fragments: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for img_id, (lt, qlen) in all_hits.items():
        lt_to_fragments[lt].append((img_id, qlen))

    extra_mapping = {}
    discarded = []
    for lt, fragments in lt_to_fragments.items():
        # Sort by length descending — longest fragment wins
        fragments.sort(key=lambda x: x[1], reverse=True)
        winner_id, winner_len = fragments[0]
        extra_mapping[winner_id] = lt
        for loser_id, loser_len in fragments[1:]:
            discarded.append((loser_id, lt, loser_len))

    return extra_mapping, discarded


def remap_img_headers(
    img_seqs: dict[str, str], img_gff_path: Path | None
) -> dict[str, str]:
    """If IMG FASTA headers are numeric OIDs, remap to AEZ55_* gene IDs.

    Returns the (possibly remapped) sequences dict.
    """
    # Check if headers look like AEZ55_* already
    sample_ids = list(img_seqs.keys())[:5]
    if any(id.startswith("AEZ55_") for id in sample_ids):
        return img_seqs  # Already has gene IDs

    if img_gff_path is None:
        print(
            "WARNING: IMG FASTA headers are not AEZ55_* format and no --img-gff provided."
        )
        print(f"  Sample headers: {sample_ids}")
        print("  Cannot remap to gene IDs. Proceeding with raw headers.")
        return img_seqs

    # Parse GFF for OID → gene_id mapping
    # IMG GFF format: scf... \t source \t CDS \t ... \t ID=NNNN;locus_tag=AEZ55_NNNN
    oid_to_gene = {}
    with open(img_gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            # Try multiple attribute formats
            gene_match = re.search(r"gene_id\s+(\S+)", attrs)
            if not gene_match:
                gene_match = re.search(r"locus_tag=([^;]+)", attrs)
            id_match = re.search(r"ID=(\d+)", attrs)
            if gene_match and id_match:
                oid_to_gene[id_match.group(1)] = gene_match.group(1)

    if not oid_to_gene:
        print("WARNING: Could not extract OID→gene_id mapping from GFF.")
        return img_seqs

    remapped = {}
    unmapped_headers = 0
    for header, seq in img_seqs.items():
        gene_id = oid_to_gene.get(header)
        if gene_id:
            remapped[gene_id] = seq
        else:
            remapped[header] = seq  # Keep original
            unmapped_headers += 1

    if unmapped_headers:
        print(f"  Remapped {len(remapped) - unmapped_headers} headers via GFF; {unmapped_headers} unmapped")
    else:
        print(f"  Remapped all {len(remapped)} headers via GFF")

    return remapped


def check_de_coverage(
    mapping: dict[str, str], csv_paths: list[Path], skip_rows: int = 2
) -> None:
    """Report how many DE genes from supplementary tables are covered."""
    de_genes = set()
    for csv_path in csv_paths:
        with open(csv_path) as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                if i < skip_rows + 1:  # skip_rows header lines + column header
                    continue
                if row and row[0].strip():
                    de_genes.add(row[0].strip())

    mapped = sum(1 for g in de_genes if g in mapping)
    unmapped = sorted(g for g in de_genes if g not in mapping)
    print(f"\nDE gene coverage: {mapped}/{len(de_genes)} ({100*mapped/len(de_genes):.1f}%)")
    if unmapped:
        print(f"Unmapped DE genes ({len(unmapped)}):")
        for g in unmapped[:20]:
            print(f"  {g}")
        if len(unmapped) > 20:
            print(f"  ... and {len(unmapped) - 20} more")


def main():
    parser = argparse.ArgumentParser(
        description="Map IMG gene IDs to NCBI locus tags via protein sequence matching"
    )
    parser.add_argument(
        "--img-faa",
        type=Path,
        required=True,
        help="IMG protein FASTA (e.g., 2785510739.genes.faa)",
    )
    parser.add_argument(
        "--ncbi-faa",
        type=Path,
        required=True,
        help="NCBI protein FASTA (e.g., cache/.../protein.faa)",
    )
    parser.add_argument(
        "--gene-mapping",
        type=Path,
        required=True,
        help="gene_mapping.csv (protein_id → locus_tag)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output id_translation CSV",
    )
    parser.add_argument(
        "--img-gff",
        type=Path,
        default=None,
        help="IMG GFF for header remapping (if FASTA headers are OIDs, not gene IDs)",
    )
    parser.add_argument(
        "--source-id-col",
        type=str,
        default="source_id",
        help="Column name for source IDs in output CSV (default: source_id)",
    )
    parser.add_argument(
        "--de-csvs",
        type=Path,
        nargs="*",
        default=[],
        help="DE CSV files to check coverage against",
    )
    parser.add_argument(
        "--skip-rows",
        type=int,
        default=2,
        help="Number of header rows to skip in DE CSVs (default: 2)",
    )
    args = parser.parse_args()

    # Validate inputs
    for p in [args.img_faa, args.ncbi_faa, args.gene_mapping]:
        if not p.exists():
            print(f"ERROR: {p} not found")
            sys.exit(1)

    print(f"IMG protein FASTA: {args.img_faa}")
    print(f"NCBI protein FASTA: {args.ncbi_faa}")
    print(f"Gene mapping: {args.gene_mapping}")
    print(f"Output: {args.output}")

    # Step 1: Parse FASTA files
    print("\nParsing FASTA files...")
    img_seqs = parse_fasta(args.img_faa)
    ncbi_seqs = parse_fasta(args.ncbi_faa)
    print(f"  IMG proteins: {len(img_seqs)}")
    print(f"  NCBI proteins: {len(ncbi_seqs)}")

    # Step 1b: Remap IMG headers if needed
    img_seqs = remap_img_headers(img_seqs, args.img_gff)

    # Step 2: Load gene_mapping
    print("\nLoading gene_mapping.csv...")
    protein_to_lt = load_gene_mapping(args.gene_mapping)
    print(f"  protein_id → locus_tag entries: {len(protein_to_lt)}")

    # Step 3: Exact sequence matching
    print("\nPhase 1: Exact sequence matching...")
    mapping, unmatched, used_ncbi = map_by_exact_match(
        img_seqs, ncbi_seqs, protein_to_lt
    )
    print(f"  Exact matches: {len(mapping)}")
    print(f"  Unmatched: {len(unmatched)}")

    # Step 4: Subsequence matching for remaining
    if unmatched:
        print("\nPhase 2: Subsequence matching (boundary differences)...")
        extra = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, unmatched, used_ncbi
        )
        mapping.update(extra)
        unmatched -= set(extra.keys())
        print(f"  Subsequence matches: {len(extra)}")
        print(f"  Still unmatched: {len(unmatched)}")

    # Step 4b: Diamond similarity search for remaining
    #
    # Draft genome annotations often split one NCBI gene into multiple
    # smaller ORFs (frameshifts / assembly gaps in the draft). Diamond
    # catches these high-identity, low-subject-coverage hits. When
    # multiple draft fragments hit the same NCBI gene, only the longest
    # fragment is kept — it captured the most RNA-seq reads and thus
    # best represents the gene's expression. Shorter fragments are
    # discarded to avoid duplicate/conflicting expression edges in the KG.
    if unmatched:
        print("\nPhase 3: Diamond blastp (near-identical proteins)...")
        extra, discarded = map_by_diamond(
            img_seqs, ncbi_seqs, protein_to_lt, unmatched, used_ncbi
        )
        mapping.update(extra)
        unmatched -= set(extra.keys())
        # Also remove discarded fragments from unmatched — they had hits
        # but were superseded by a longer fragment
        discarded_ids = {d[0] for d in discarded}
        unmatched -= discarded_ids
        print(f"  Diamond matches: {len(extra)}")
        print(f"  Discarded fragments (shorter sibling of same gene): {len(discarded)}")
        if discarded:
            for img_id, lt, qlen in discarded[:10]:
                print(f"    {img_id} ({qlen} aa) → {lt} (kept longer fragment)")
            if len(discarded) > 10:
                print(f"    ... and {len(discarded) - 10} more")
        print(f"  Still unmatched: {len(unmatched)}")

    # Step 5: Summary
    print(f"\nTotal mapping: {len(mapping)} / {len(img_seqs)} ({100*len(mapping)/len(img_seqs):.1f}%)")

    # Check for duplicate locus_tag assignments
    lt_counts = defaultdict(list)
    for img_id, lt in mapping.items():
        lt_counts[lt].append(img_id)
    duplicates = {lt: ids for lt, ids in lt_counts.items() if len(ids) > 1}
    if duplicates:
        print(f"\nWARNING: {len(duplicates)} locus_tags matched by multiple IMG genes:")
        for lt, ids in list(duplicates.items())[:10]:
            print(f"  {lt} ← {ids}")

    # Step 6: Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([args.source_id_col, "locus_tag"])
        for img_id in sorted(mapping.keys()):
            writer.writerow([img_id, mapping[img_id]])
    print(f"\nWrote {len(mapping)} mappings to {args.output}")

    # Step 7: DE gene coverage
    if args.de_csvs:
        check_de_coverage(mapping, args.de_csvs, args.skip_rows)

    # Return exit code based on mapping quality
    coverage = len(mapping) / len(img_seqs) if img_seqs else 0
    if coverage < 0.5:
        print("\nWARNING: <50% mapping — check input files")
        sys.exit(2)


if __name__ == "__main__":
    main()
