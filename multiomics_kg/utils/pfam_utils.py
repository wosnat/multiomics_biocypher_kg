"""
Pfam reference data utilities.

Built by prepare_data step 9 (``multiomics_kg.download.build_pfam_reference``) from
the current Pfam release: ``Pfam-A.clans.tsv.gz`` (families, short names, clans),
``Pfam-A.dead.gz`` (killed families and the family each was merged into) and
``Pfam.version.gz``. Pfam 33.1's clans file adds eggNOG's old short names.

Cache layout:
  cache_root/pfam/pfam_reference.json  -- committed {release, by_accession, by_shortname, clans, dead}
  cache_root/pfam/raw/                 -- gitignored downloads (re-fetched by --refetch-raw)
"""

from __future__ import annotations

import gzip
import json
import logging
import re
from dataclasses import dataclass, asdict, field
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
PFAM_DEAD_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.dead.gz"
PFAM_VERSION_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz"
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
    # Killed families: PF01188 -> {"shortname", "forward_to": ["PF13378"], "comment"}
    dead: dict[str, dict] = field(default_factory=dict)
    release: str = ""                    # "39.0"

    def resolve_accession(self, accession: str) -> list[str]:
        """Live accession(s) for ``accession``.

        A current family maps to itself. A killed family follows Pfam's
        ``forward_to`` chain to the family it was merged into; one killed with
        no successor ("not what it seems", removed TM helix) maps to ``[]``.
        An accession in neither table is kept as-is: the reference may simply
        predate it.
        """
        if accession in self.by_accession or accession not in self.dead:
            return [accession]
        live: list[str] = []
        seen: set[str] = set()
        todo = [accession]
        while todo:
            acc = todo.pop(0)
            if acc in seen:
                continue
            seen.add(acc)
            if acc in self.by_accession:
                if acc not in live:
                    live.append(acc)
            elif acc in self.dead:
                todo.extend(self.dead[acc].get("forward_to", []))
        return live


_ACCESSION_RE = re.compile(r"^PF(\d{1,5})(?:\.\d+)?$")


def normalize_pfam_token(token: str, pfam_data: PfamData) -> str | None:
    """Turn a raw Pfam token (accession or short name) into a PF accession.

    An accession is "PF" + digits, with an optional version suffix. Cyanorak
    writes some with a digit missing ("PF7969" next to IPR013108 = PF07969), so
    a short one is zero-padded, but only when the padded accession is known.
    Anything else is a short name ("PFK", "PFOR_II" begin with PF too) and is
    looked up as one. Returns None when the token resolves to nothing.
    """
    m = _ACCESSION_RE.match(token)
    if m:
        digits = m.group(1)
        acc = f"PF{digits.zfill(5)}"
        if len(digits) == 5 or acc in pfam_data.by_accession or acc in pfam_data.dead:
            return acc
        return None
    return pfam_data.by_shortname.get(token)


def parse_dead_families(gz_path: Path) -> dict[str, dict]:
    """Parse Pfam-A.dead.gz (Stockholm-style records) into {accession: info}.

    Each record carries ``#=GF ID`` (short name), ``#=GF AC``, ``#=GF FW``
    (forward-to accessions, space-separated, often empty) and ``#=GF CC``
    comment lines.
    """
    dead: dict[str, dict] = {}
    rec: dict = {}

    def _flush() -> None:
        if rec.get("AC"):
            dead[rec["AC"]] = {
                "shortname": rec.get("ID", ""),
                "forward_to": rec.get("FW", []),
                "comment": " ".join(rec.get("CC", [])),
            }

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("//"):
                _flush()
                rec = {}
                continue
            if not line.startswith("#=GF "):
                continue
            tag, _, value = line[5:].partition(" ")
            value = value.strip()
            if tag == "ID":
                rec["ID"] = value
            elif tag == "AC":
                rec["AC"] = value.split(".")[0]
            elif tag == "FW":
                rec["FW"] = [a.split(".")[0] for a in re.split(r"[\s;,]+", value) if a.startswith("PF")]
            elif tag == "CC":
                rec.setdefault("CC", []).append(value)
    _flush()
    logger.info("Parsed %d killed Pfam families from %s", len(dead), gz_path.name)
    return dead


def _pfam_cache_dir(cache_root: Path) -> Path:
    d = cache_root / "pfam"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download(url: str, path: Path) -> None:
    """Download ``url`` to ``path`` (written via a temp file, so a failure leaves no stub)."""
    logger.info("Downloading %s ...", url)
    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()
    tmp = path.with_suffix(path.suffix + ".part")
    with open(tmp, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    with gzip.open(tmp, "rb") as fh:  # every Pfam file we fetch is gzip; an HTML error page is not
        fh.read(1)
    tmp.replace(path)
    logger.info("Downloaded %s (%d KB)", path.name, path.stat().st_size // 1024)


def parse_release(gz_path: Path) -> str:
    """Pfam release number from Pfam.version.gz ("Pfam release       : 39.0")."""
    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.lower().startswith("pfam release"):
                return line.split(":", 1)[1].strip()
    return ""


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


def _parse_old_shortnames(gz_path: Path, current_accessions: set[str]) -> dict[str, str]:
    """Find short names renamed between Pfam 33.1 (eggNOG's release) and the current one.

    Returns a dict of old_shortname → current_accession for shortnames that
    were renamed between Pfam 33.1 and the current release (same accession,
    different shortname).
    """
    old_shortname_to_acc: dict[str, str] = {}
    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[3].strip():
                old_shortname_to_acc[parts[3].strip()] = parts[0].strip()

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
        "dead": pfam_data.dead,
        "release": pfam_data.release,
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
        dead=data.get("dead", {}),
        release=data.get("release", ""),
    )


RAW_FILES = {
    "clans": ("Pfam-A.clans.tsv.gz", PFAM_CLANS_URL),
    "dead": ("Pfam-A.dead.gz", PFAM_DEAD_URL),
    "version": ("Pfam.version.gz", PFAM_VERSION_URL),
    "clans_33": ("Pfam33.1-clans.tsv.gz", PFAM_33_CLANS_URL),
}


def reference_path(cache_root: Path) -> Path:
    return _pfam_cache_dir(cache_root) / "pfam_reference.json"


def build_pfam_reference(cache_root: Path, refetch_raw: bool = False) -> PfamData:
    """Build pfam_reference.json from the raw Pfam files under ``cache_root/pfam/raw``.

    Raw files are downloaded when missing or when ``refetch_raw`` is set.
    """
    raw_dir = _pfam_cache_dir(cache_root) / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, (name, url) in RAW_FILES.items():
        paths[key] = raw_dir / name
        if refetch_raw or not paths[key].exists():
            _download(url, paths[key])

    pfam_data = parse_clans_tsv(paths["clans"])
    pfam_data.dead = parse_dead_families(paths["dead"])
    pfam_data.release = parse_release(paths["version"])

    # Old short names from Pfam 33.1 (eggNOG's version) for rename resolution
    renames = _parse_old_shortnames(paths["clans_33"], set(pfam_data.by_accession))
    added = 0
    for old_name, acc in renames.items():
        if old_name not in pfam_data.by_shortname:
            pfam_data.by_shortname[old_name] = acc
            added += 1
    logger.info("Added %d old shortnames from Pfam 33.1 to lookup", added)

    cache_path = reference_path(cache_root)
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(_pfam_data_to_json(pfam_data), fh, separators=(",", ":"), sort_keys=True)
    logger.info(
        "Wrote Pfam %s reference: %d entries, %d shortnames, %d clans, %d killed, %d KB -> %s",
        pfam_data.release, len(pfam_data.by_accession), len(pfam_data.by_shortname),
        len(pfam_data.clans), len(pfam_data.dead), cache_path.stat().st_size // 1024, cache_path,
    )
    return pfam_data


def load_pfam_data(cache_root: Path) -> PfamData:
    """Load the committed Pfam reference (built by prepare_data step 9).

    Builds it, downloading the raw files, only when it is missing.
    """
    cache_path = reference_path(cache_root)
    if not cache_path.exists():
        logger.warning("Pfam reference missing (%s); building it now", cache_path)
        return build_pfam_reference(cache_root)
    logger.info("Loading Pfam reference from cache: %s", cache_path)
    with open(cache_path, encoding="utf-8") as fh:
        pfam_data = _pfam_data_from_json(json.load(fh))
    logger.info(
        "Loaded Pfam %s: %d entries, %d shortnames, %d clans, %d killed",
        pfam_data.release or "(release unknown)", len(pfam_data.by_accession),
        len(pfam_data.by_shortname), len(pfam_data.clans), len(pfam_data.dead),
    )
    return pfam_data
