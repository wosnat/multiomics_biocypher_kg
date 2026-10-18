"""Pin ONE GO release: go-basic.obo + go-plus.json downloaded together.

Usage: uv run python -m multiomics_kg.download.refresh_go_release [--force]

Writes (under cache/data/go_terms/):
  raw/go-basic.obo, raw/go-plus.json   (gitignored)
  go_namespace_cache.json              (committed; rebuilt from go-basic)
  go_alt_ids.json                      (committed; {alt_id: primary} of live terms)
  go_release.txt                       (committed; the pinned release)
GO is refreshed only here, never during a KG build.
"""
import argparse
import json
import logging
import re
from pathlib import Path

import requests

from multiomics_kg.utils.go_utils import (
    GO_BASIC_URL,
    GO_RELEASE_FILE,
    parse_obo,
    parse_obo_release,
    read_go_release,
    write_alt_ids,
    write_namespace_cache,
)

logger = logging.getLogger(__name__)

GO_PLUS_URL = "http://purl.obolibrary.org/obo/go/extensions/go-plus.json"
DEFAULT_CACHE_ROOT = Path("cache/data")

_VERSION_INFO_RE = re.compile(
    r'"pred"\s*:\s*"http://www\.w3\.org/2002/07/owl#versionInfo"\s*,\s*"val"\s*:\s*"([^"]+)"'
)


class GoReleaseMismatch(RuntimeError):
    """go-basic and go-plus come from different GO releases."""


def _download(url: str, dest: Path) -> None:
    logger.info(f"Downloading {url} ...")
    with requests.get(url, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)


def go_plus_release(path: Path) -> str | None:
    """owl:versionInfo of graphs[0].meta in a go-plus.json file."""
    with open(path, encoding="utf-8") as fh:
        head = fh.read(1 << 20)
    m = _VERSION_INFO_RE.search(head)
    if m:
        return m.group(1)
    with open(path, encoding="utf-8") as fh:
        meta = json.load(fh)["graphs"][0].get("meta", {})
    for bpv in meta.get("basicPropertyValues", []):
        if bpv.get("pred", "").endswith("owl#versionInfo"):
            return bpv["val"]
    return None


def refresh_go_release(cache_root: Path = DEFAULT_CACHE_ROOT, force: bool = False) -> str:
    """Pin one GO release in ``cache_root/go_terms``; return the release string."""
    gdir = Path(cache_root) / "go_terms"
    raw = gdir / "raw"
    basic, plus = raw / "go-basic.obo", raw / "go-plus.json"
    cache_json = gdir / "go_namespace_cache.json"

    if not force and basic.exists() and plus.exists() and cache_json.exists():
        pinned = read_go_release(cache_root)
        if pinned and parse_obo_release(basic) == pinned and go_plus_release(plus) == pinned:
            logger.info(f"GO release {pinned} already pinned; nothing to download")
            return pinned

    raw.mkdir(parents=True, exist_ok=True)
    tmp_basic, tmp_plus = raw / "go-basic.obo.tmp", raw / "go-plus.json.tmp"
    try:
        _download(GO_BASIC_URL, tmp_basic)
        _download(GO_PLUS_URL, tmp_plus)
        r_basic, r_plus = parse_obo_release(tmp_basic), go_plus_release(tmp_plus)
        if not r_basic or r_basic != r_plus:
            raise GoReleaseMismatch(f"go-basic release {r_basic!r} != go-plus release {r_plus!r}")
        tmp_basic.replace(basic)
        tmp_plus.replace(plus)
    finally:
        for t in (tmp_basic, tmp_plus):
            if t.exists():
                t.unlink()

    write_namespace_cache(cache_json, parse_obo(basic))
    write_alt_ids(cache_root, basic)
    (gdir / GO_RELEASE_FILE).write_text(r_basic + "\n", encoding="utf-8")
    return r_basic


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO)
    print(refresh_go_release(args.cache_root, force=args.force))


if __name__ == "__main__":
    main()
