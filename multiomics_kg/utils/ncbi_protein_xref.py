"""NCBI Identical Protein Group (IPG) cross-references for legacy accessions.

Old papers cite NCBI protein records by accessions that have been superseded
in modern RefSeq:

- ``NP_xxx``, ``YP_xxx``, ``XP_xxx`` — RefSeq specific-genome accessions
  (replaced by ``WP_xxx`` non-redundant entries when the same protein appears
  in multiple genomes)
- ``CAExxx``, ``KGFxxx``, ``BAA*``, ``AAA*`` — INSDC (GenBank/EMBL/DDBJ)
  protein accessions tied to the original genome submission

The IPG database groups all accessions that point to the same protein
sequence. Querying NCBI E-utilities ``efetch`` with ``rettype=ipg`` returns,
for any input protein accession, every other accession in the same group
plus the assembly each accession belongs to.

This module:

1. Provides ``fetch_ipg_xrefs(accessions, *, cache_dir, ...)`` — batched
   lookup with per-accession disk cache. Cache hits are free; misses are
   batched into ``efetch`` POST requests (200 IDs per call) and rate-limited
   to NCBI's published limits.
2. Provides ``find_assembly_match(entries, target_assembly)`` — picks the
   IPG row whose Assembly column matches a given ``GCF_xxx``/``GCA_xxx``
   accession (with .version stripped for comparison).

The downstream pipeline (``build_gene_id_mapping``) uses this to enrich
each strain's gene-ID mapping with legacy-accession → current-accession
pairs, after which the existing iterative-convergence resolver picks them
up automatically.
"""
from __future__ import annotations

import csv
import io
import os
import re
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import dotenv
import requests

# Pick up NCBI_API_KEY (and any other secrets) from .env at the project root.
# load_dotenv silently no-ops if no .env exists.
dotenv.load_dotenv()

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 200
DEFAULT_CACHE_PATH = Path("cache/data/ncbi_ipg.sqlite")

# NCBI rate limits: 3 req/sec without API key, 10 req/sec with.
RATE_DELAY_NO_KEY = 0.34   # sec between calls
RATE_DELAY_WITH_KEY = 0.11

# Heuristic regex for "looks like an NCBI protein accession we might want to
# look up via IPG". Matches RefSeq prefixed accessions (NP_/YP_/XP_/WP_) and
# INSDC accessions (3 letters + 5+ digits, optional .version). Locus-tag-like
# strings (e.g. PMM0001, TX50_RS00020) are intentionally excluded.
LEGACY_PROTEIN_ACC_RE = re.compile(
    r"^(?:[NXYW]P_\d+(?:\.\d+)?|[A-Z]{3}\d{5,8}(?:\.\d+)?)$"
)


@dataclass(frozen=True)
class IpgEntry:
    """One row of an IPG report TSV."""
    source_db: str            # "RefSeq" or "INSDC"
    nucleotide_accession: str
    start: int
    stop: int
    strand: str
    protein_accession: str    # the current accession in this row's source
    protein_name: str
    organism: str
    strain: str
    assembly: str             # "GCF_xxx.N" or "GCA_xxx.N"


def _accession_stem(acc: str) -> str:
    """Strip the .version suffix; "GCF_000011465.1" → "GCF_000011465"."""
    return acc.split(".", 1)[0]


def _parse_ipg_tsv(text: str) -> list[IpgEntry]:
    """Parse the TSV body returned by efetch?rettype=ipg.

    The header row is:
        Id  Source  Nucleotide Accession  Start  Stop  Strand  Protein
        Protein Name  Organism  Strain  Assembly

    Empty / partial rows are skipped silently.
    """
    rows: list[IpgEntry] = []
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    header = next(reader, None)
    if not header:
        return rows
    for fields in reader:
        if len(fields) < 11:
            continue
        try:
            start = int(fields[3]) if fields[3] else 0
            stop = int(fields[4]) if fields[4] else 0
        except ValueError:
            continue
        rows.append(IpgEntry(
            source_db=fields[1],
            nucleotide_accession=fields[2],
            start=start,
            stop=stop,
            strand=fields[5],
            protein_accession=fields[6],
            protein_name=fields[7],
            organism=fields[8],
            strain=fields[9],
            assembly=fields[10],
        ))
    return rows


def _open_cache_db(path: Path) -> sqlite3.Connection:
    """Open (and lazily create) the SQLite cache. One row per accession;
    body holds the raw efetch?rettype=ipg TSV (header + zero or more rows)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS ipg_cache (
            accession  TEXT PRIMARY KEY,
            body       TEXT NOT NULL,
            fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    return conn


def _read_cached(acc: str, db: sqlite3.Connection) -> list[IpgEntry] | None:
    """Return cached IPG entries for ``acc``, or None on cache miss.
    An empty list (cached header-only body) means "we asked, no results"."""
    cur = db.execute("SELECT body FROM ipg_cache WHERE accession = ?", (acc,))
    row = cur.fetchone()
    if row is None:
        return None
    return _parse_ipg_tsv(row[0])


def _write_cache(acc: str, body: str, db: sqlite3.Connection) -> None:
    db.execute(
        "INSERT OR REPLACE INTO ipg_cache (accession, body) VALUES (?, ?)",
        (acc, body),
    )
    db.commit()


def _split_returned_ipg_by_accession(
    requested: list[str],
    body: str,
) -> dict[str, str]:
    """Split a multi-accession IPG response into per-accession TSV bodies.

    Strategy: NCBI returns IPG groups in the SAME ORDER as the requested
    accessions. We split the body by IPG group ID and attribute groups to
    inputs positionally. If the input accession itself appears in the
    group (e.g. modern WP_), we also confirm the match.

    When the group count doesn't equal len(requested) — typically because
    some requested accessions are invalid/deleted/unmapped — we still
    attribute groups by position to as many inputs as possible. The
    remaining inputs get empty results (cached as header-only files).
    """
    lines = body.splitlines()
    if not lines:
        return {acc: "" for acc in requested}
    header_line = lines[0]

    groups: list[tuple[str, list[str]]] = []
    current_id: str | None = None
    current: list[str] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        ipg_id = fields[0]
        if ipg_id != current_id:
            if current_id is not None:
                groups.append((current_id, current))
            current_id = ipg_id
            current = [line]
        else:
            current.append(line)
    if current_id is not None:
        groups.append((current_id, current))

    result: dict[str, str] = {}
    for i, acc in enumerate(requested):
        if i < len(groups):
            _, group_lines = groups[i]
            result[acc] = "\n".join([header_line, *group_lines]) + "\n"
        else:
            # Cache empty (header-only) so we don't re-fetch on next run.
            result[acc] = header_line + "\n"
    return result


def _fetch_batch(
    accessions: list[str],
    *,
    api_key: str | None,
    timeout: float,
) -> dict[str, str]:
    """One efetch call for up to BATCH_SIZE accessions. Returns per-accession
    TSV bodies (header + group rows)."""
    params = {
        "db": "protein",
        "id": ",".join(accessions),
        "rettype": "ipg",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    resp = requests.post(EFETCH_URL, data=params, timeout=timeout)
    resp.raise_for_status()
    return _split_returned_ipg_by_accession(accessions, resp.text)


def fetch_ipg_xrefs(
    accessions: Iterable[str],
    *,
    cache_path: Path = DEFAULT_CACHE_PATH,
    api_key: str | None = None,
    timeout: float = 60.0,
    progress: bool = False,
) -> dict[str, list[IpgEntry]]:
    """Fetch IPG cross-references for a list of protein accessions.

    Cached results are returned without network access. Misses are batched
    and POSTed to NCBI E-utilities ``efetch?rettype=ipg``.

    Args:
        accessions: any iterable of protein accession strings; duplicates
            are de-duplicated before lookup.
        cache_path: SQLite cache file path; auto-created with the
            ``ipg_cache`` schema on first use.
        api_key: optional NCBI API key (raises rate limit from 3 → 10 req/sec).
            Defaults to ``NCBI_API_KEY`` env var when unset.
        timeout: per-batch HTTP timeout in seconds.
        progress: print one line per batch when True.

    Returns:
        ``{accession: [IpgEntry, ...]}``. An accession with no IPG hits
        is mapped to an empty list (and cached as such).
    """
    if api_key is None:
        api_key = os.environ.get("NCBI_API_KEY") or None
    delay = RATE_DELAY_WITH_KEY if api_key else RATE_DELAY_NO_KEY

    unique = sorted({a.strip() for a in accessions if a and a.strip()})
    result: dict[str, list[IpgEntry]] = {}

    db = _open_cache_db(cache_path)
    try:
        misses: list[str] = []
        for acc in unique:
            cached = _read_cached(acc, db)
            if cached is not None:
                result[acc] = cached
            else:
                misses.append(acc)

        if progress and misses:
            print(f"  [ipg] {len(unique)} accessions, {len(misses)} cache misses")

        last_call = 0.0
        for i in range(0, len(misses), BATCH_SIZE):
            batch = misses[i : i + BATCH_SIZE]
            wait = delay - (time.monotonic() - last_call)
            if wait > 0:
                time.sleep(wait)
            if progress:
                print(f"  [ipg] batch {i // BATCH_SIZE + 1}/"
                      f"{(len(misses) + BATCH_SIZE - 1) // BATCH_SIZE} "
                      f"({len(batch)} accessions)")
            bodies = _fetch_batch(batch, api_key=api_key, timeout=timeout)
            last_call = time.monotonic()
            for acc, body in bodies.items():
                _write_cache(acc, body, db)
                result[acc] = _parse_ipg_tsv(body)
    finally:
        db.close()

    return result


def find_assembly_match(
    entries: list[IpgEntry],
    target_assembly: str,
) -> IpgEntry | None:
    """Pick the IPG entry whose Assembly column matches ``target_assembly``.

    Comparison strips ``.version`` so ``GCF_000011465.1`` matches
    ``GCF_000011465.2``. RefSeq (GCF) is preferred over INSDC (GCA).
    """
    target = _accession_stem(target_assembly)
    refseq_hit = next(
        (e for e in entries if _accession_stem(e.assembly) == target and e.source_db == "RefSeq"),
        None,
    )
    if refseq_hit:
        return refseq_hit
    return next(
        (e for e in entries if _accession_stem(e.assembly) == target),
        None,
    )


def looks_like_legacy_protein_acc(value: str) -> bool:
    """True if ``value`` matches the legacy NCBI protein-accession pattern."""
    return bool(LEGACY_PROTEIN_ACC_RE.match(value.strip()))
