#!/usr/bin/env python3
"""Overnight orchestrator for eggNOG on 4 new strains.

Polls for each strain's protein.faa (expected cache path), then calls
run_eggnog.py --strain <name>. Writes progress to
logs/eggnog_overnight_2026-04-14.md.

This helper lives alongside the log so the background agent's writes
stay scoped to logs/ and cache/ eggnog output directories.
"""
from __future__ import annotations

import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO = Path("/home/osnat/github/multiomics_biocypher_kg")
LOG = REPO / "logs/eggnog_overnight_2026-04-14.md"
CSV = REPO / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
RUN_EGGNOG = REPO / ".claude/skills/eggnog-run/run_eggnog.py"

MAX_SECONDS = 8 * 3600  # 8 hours global cap
POLL_SECONDS = 60
CSV_EXTRA_WAIT = 600  # wait up to 10 min after FASTA for CSV row

STRAINS = [
    ("SS120",     "cache/data/Prochlorococcus/genomes/SS120"),
    ("BL107",     "cache/data/Synechococcus/genomes/BL107"),
    ("HP15",      "cache/data/Marinobacter/genomes/HP15"),
    ("AltMedDE",  "cache/data/Alteromonas/genomes/AltMedDE"),
]


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log_event(msg: str) -> None:
    with LOG.open("a") as f:
        f.write(f"- {ts()} | {msg}\n")


def update_row(strain: str, fasta: str, started: str, finished: str, status: str) -> None:
    txt = LOG.read_text()
    pattern = re.compile(rf"^\| {re.escape(strain)} \|.*$", re.MULTILINE)
    new = f"| {strain} | {fasta} | {started} | {finished} | {status} |"
    txt2, n = pattern.subn(new, txt, count=1)
    if n == 0:
        log_event(f"[{strain}] WARNING: could not find table row to update")
        return
    LOG.write_text(txt2)


def csv_has_strain(strain: str) -> bool:
    if not CSV.exists():
        return False
    for line in CSV.read_text().splitlines():
        if line.startswith("#"):
            continue
        cols = line.split(",")
        # column index 3 is strain_name
        if len(cols) > 3 and cols[3].strip() == strain:
            return True
    return False


def count_annotations(path: Path) -> int:
    n = 0
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            n += 1
    return n


def main() -> int:
    start_epoch = time.time()
    log_event(f"orchestrator started (PID {__import__('os').getpid()})")

    def elapsed() -> int:
        return int(time.time() - start_epoch)

    for strain, cache_rel in STRAINS:
        cache_dir = REPO / cache_rel
        faa = cache_dir / "protein.faa"
        out_annot = cache_dir / "eggnog" / f"{strain}.emapper.annotations"

        log_event(f"[{strain}] waiting for {faa}")
        update_row(strain, "-", "-", "-", "waiting for FASTA")

        # Poll for FASTA
        while not (faa.exists() and faa.stat().st_size > 0):
            if elapsed() > MAX_SECONDS:
                log_event(f"[{strain}] global timeout waiting for FASTA; skipping remaining")
                update_row(strain, "not found", "-", "-", "TIMEOUT (no FASTA after 8h)")
                log_event("orchestrator exiting (global timeout)")
                return 0
            time.sleep(POLL_SECONDS)

        try:
            lines = sum(1 for _ in faa.open())
        except Exception:
            lines = -1
        log_event(f"[{strain}] FASTA detected at {faa} ({lines} lines)")
        update_row(strain, str(faa), "-", "-", "FASTA found")

        # Wait briefly for CSV row (required by run_eggnog.py)
        if not csv_has_strain(strain):
            log_event(f"[{strain}] not in cyanobacteria_genomes.csv yet; waiting up to {CSV_EXTRA_WAIT}s")
            wait = 0
            while wait < CSV_EXTRA_WAIT:
                if csv_has_strain(strain):
                    log_event(f"[{strain}] CSV row appeared after {wait}s")
                    break
                time.sleep(30)
                wait += 30
            if not csv_has_strain(strain):
                log_event(f"[{strain}] CSV row never appeared; cannot run -- skipping")
                update_row(strain, str(faa), "-", "-", "FAILED: not in genomes.csv")
                continue

        start_ts = ts()
        update_row(strain, str(faa), start_ts, "-", "running")
        log_event(f"[{strain}] invoking run_eggnog.py --strain {strain} --cpu 8")

        strain_log = REPO / f"logs/eggnog_overnight_{strain}.log"
        with strain_log.open("w") as lf:
            proc = subprocess.run(
                ["uv", "run", "python", str(RUN_EGGNOG), "--strain", strain, "--cpu", "8"],
                cwd=str(REPO),
                stdout=lf,
                stderr=subprocess.STDOUT,
            )
        end_ts = ts()
        rc = proc.returncode

        if rc == 0 and out_annot.exists() and out_annot.stat().st_size > 0:
            n = count_annotations(out_annot)
            update_row(strain, str(faa), start_ts, end_ts, f"OK ({n} annotations)")
            log_event(f"[{strain}] completed OK; {n} annotations in {out_annot}")
        else:
            update_row(strain, str(faa), start_ts, end_ts, f"FAILED (rc={rc}, see {strain_log.name})")
            log_event(f"[{strain}] FAILED rc={rc}; log at {strain_log}")

    log_event(f"orchestrator finished after {elapsed()}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
