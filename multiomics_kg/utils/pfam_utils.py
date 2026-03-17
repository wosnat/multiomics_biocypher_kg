"""
Pfam reference data utilities.

Downloads Pfam-A.clans.tsv.gz (~494 KB, ~27.5K entries) from EBI FTP,
parses it into a compact JSON cache, and provides lookup structures for
resolving Pfam shortnames to accessions.

Cache layout:
  cache_root/pfam/pfam_reference.json  -- compact {by_accession, by_shortname, clans}
  Pfam-A.clans.tsv.gz                  -- downloaded temporarily, deleted after parsing
"""

from __future__ import annotations

import gzip
import json
import logging
from dataclasses import dataclass, asdict
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
# EggNOG v2.1.13 uses Pfam 33.1 (2020). ~1,165 shortnames were renamed between 33.1 and current.
# We download the old release to build a rename map (old shortname → current accession).
PFAM_33_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.clans.tsv.gz"


@dataclass
class PfamEntry:
    accession: str        # PF00712
    shortname: str        # DNA_pol3_beta
    description: str      # DNA polymerase III beta subunit, N-terminal domain
    clan_accession: str   # CL0060 or ""
    clan_name: str        # DNA_clamp or ""


@dataclass
class PfamData:
    by_accession: dict[str, PfamEntry]   # PF00712 -> PfamEntry
    by_shortname: dict[str, str]         # DNA_pol3_beta -> PF00712
    clans: dict[str, str]                # CL0060 -> "DNA_clamp"


def _pfam_cache_dir(cache_root: Path) -> Path:
    d = cache_root / "pfam"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_clans_tsv(gz_path: Path) -> None:
    """Download Pfam-A.clans.tsv.gz to gz_path."""
    logger.info("Downloading Pfam-A.clans.tsv.gz from %s ...", PFAM_CLANS_URL)
    resp = requests.get(PFAM_CLANS_URL, stream=True, timeout=120)
    resp.raise_for_status()
    with open(gz_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    logger.info("Downloaded Pfam-A.clans.tsv.gz (%d KB)", gz_path.stat().st_size // 1024)


def parse_clans_tsv(gz_path: Path) -> PfamData:
    """Parse Pfam-A.clans.tsv.gz and return PfamData.

    TSV format (tab-separated, no header, 5 columns):
        accession | clan_accession | clan_name | shortname | description
    Entries without a clan have empty clan_accession and clan_name fields.
    """
    by_accession: dict[str, PfamEntry] = {}
    by_shortname: dict[str, str] = {}
    clans: dict[str, str] = {}

    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue

            accession = parts[0].strip()
            clan_acc = parts[1].strip()
            clan_name = parts[2].strip()
            shortname = parts[3].strip()
            description = parts[4].strip()

            if not accession:
                continue

            entry = PfamEntry(
                accession=accession,
                shortname=shortname,
                description=description,
                clan_accession=clan_acc,
                clan_name=clan_name,
            )
            by_accession[accession] = entry

            if shortname:
                by_shortname[shortname] = accession

            if clan_acc and clan_name:
                clans[clan_acc] = clan_name

    logger.info(
        "Parsed %d Pfam entries, %d shortnames, %d clans from %s",
        len(by_accession), len(by_shortname), len(clans), gz_path.name,
    )
    return PfamData(by_accession=by_accession, by_shortname=by_shortname, clans=clans)


def _download_and_parse_old_shortnames(gz_path: Path, current_accessions: set[str]) -> dict[str, str]:
    """Download Pfam 33.1 clans and find renamed shortnames.

    Returns a dict of old_shortname → current_accession for shortnames that
    were renamed between Pfam 33.1 and the current release (same accession,
    different shortname).
    """
    logger.info("Downloading Pfam 33.1 clans for shortname rename detection ...")
    resp = requests.get(PFAM_33_CLANS_URL, stream=True, timeout=120)
    resp.raise_for_status()
    with open(gz_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)

    old_shortname_to_acc: dict[str, str] = {}
    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[3].strip():
                old_shortname_to_acc[parts[3].strip()] = parts[0].strip()

    gz_path.unlink(missing_ok=True)

    # Keep only old shortnames whose accession still exists in current Pfam
    renames = {
        old_name: acc
        for old_name, acc in old_shortname_to_acc.items()
        if acc in current_accessions
    }
    logger.info(
        "Pfam 33.1: %d shortnames, %d renamed (same accession, different shortname in current)",
        len(old_shortname_to_acc), len(renames),
    )
    return renames


def _pfam_data_to_json(pfam_data: PfamData) -> dict:
    """Serialize PfamData to a JSON-compatible dict."""
    return {
        "by_accession": {
            acc: asdict(entry) for acc, entry in pfam_data.by_accession.items()
        },
        "by_shortname": pfam_data.by_shortname,
        "clans": pfam_data.clans,
    }


def _pfam_data_from_json(data: dict) -> PfamData:
    """Deserialize PfamData from a JSON dict."""
    by_accession = {
        acc: PfamEntry(**entry_dict)
        for acc, entry_dict in data["by_accession"].items()
    }
    return PfamData(
        by_accession=by_accession,
        by_shortname=data["by_shortname"],
        clans=data["clans"],
    )


def load_pfam_data(cache_root: Path, force: bool = False) -> PfamData:
    """Load Pfam reference data, downloading and caching if needed.

    Args:
        cache_root: project-level cache directory (e.g. Path("cache/data"))
        force: rebuild cache even if it exists (re-downloads TSV)

    Returns:
        PfamData with by_accession, by_shortname, and clans lookups.
    """
    cache_path = _pfam_cache_dir(cache_root) / "pfam_reference.json"

    if cache_path.exists() and not force:
        logger.info("Loading Pfam reference from cache: %s", cache_path)
        with open(cache_path, encoding="utf-8") as fh:
            data = json.load(fh)
        pfam_data = _pfam_data_from_json(data)
        logger.info(
            "Loaded %d Pfam entries, %d shortnames, %d clans",
            len(pfam_data.by_accession),
            len(pfam_data.by_shortname),
            len(pfam_data.clans),
        )
        return pfam_data

    pfam_dir = _pfam_cache_dir(cache_root)
    gz_path = pfam_dir / "Pfam-A.clans.tsv.gz"
    try:
        _download_clans_tsv(gz_path)
        pfam_data = parse_clans_tsv(gz_path)
    finally:
        if gz_path.exists():
            gz_path.unlink()

    # Add old shortnames from Pfam 33.1 (eggNOG's version) for rename resolution
    old_gz = pfam_dir / "Pfam33.1-clans.tsv.gz"
    try:
        renames = _download_and_parse_old_shortnames(
            old_gz, set(pfam_data.by_accession.keys())
        )
        added = 0
        for old_name, acc in renames.items():
            if old_name not in pfam_data.by_shortname:
                pfam_data.by_shortname[old_name] = acc
                added += 1
        logger.info("Added %d old shortnames from Pfam 33.1 to lookup", added)
    except Exception:
        logger.warning("Failed to download Pfam 33.1 for rename detection; continuing without", exc_info=True)

    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(_pfam_data_to_json(pfam_data), fh, separators=(",", ":"), sort_keys=True)
    logger.info(
        "Wrote Pfam reference cache: %d entries, %d shortnames, %d KB -> %s",
        len(pfam_data.by_accession),
        len(pfam_data.by_shortname),
        cache_path.stat().st_size // 1024,
        cache_path,
    )

    return pfam_data
