# .claude/skills/merops-diamond/run_merops_diamond.py
"""Run diamond blastp vs. the MEROPS scan library per strain. See spec
docs/superpowers/specs/2026-08-17-merops-diamond-design.md.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

# Ensure project root is on sys.path so `multiomics_kg` is importable when
# this script is run directly (uv run python .claude/skills/merops-diamond/...).
_REPO_ROOT_EARLY = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT_EARLY) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT_EARLY))

import dotenv

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.merops_diamond import (
    build_strain_calls,
    parse_family_txt,
    parse_scan_lib_header,
)
from multiomics_kg.utils.tool_calls_io import (
    calls_path,
    load_skill_summary,
    save_calls,
    save_skill_summary,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
LOGS_DIR = REPO_ROOT / "logs" / "merops"
DEFAULT_MEROPS_DATA_DIR = Path.home() / "tools" / "MEROPS"

TOOL = "merops"

dotenv.load_dotenv(REPO_ROOT / ".env")

_MEROPS_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/merops/current_release"
_SCAN_LIB_URL = f"{_MEROPS_BASE_URL}/merops_scan.lib"
_FAMILY_TXT_URL = f"{_MEROPS_BASE_URL}/database_files/family.txt"
_CLAN_TXT_URL = f"{_MEROPS_BASE_URL}/database_files/clan.txt"


def resolve_merops_data_dir() -> Path:
    """Resolve the MEROPS data dir (env override or default ~/tools/MEROPS)."""
    env = os.environ.get("MEROPS_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_MEROPS_DATA_DIR


def _rewrite_scan_lib_headers(raw_path: Path, dest_path: Path) -> int:
    """Rewrite merops_scan.lib headers to ``>MERNUM|merops_id|subfam_token``.

    Input headers: ``>MER0000002 - name (organism) [S01.001]#S01A#{...}~...~``
    Output headers: ``>MER0000002|S01.001|S01A``

    The scan library is latin-1 encoded (en-dashes in peptidase names) — NOT
    UTF-8. Returns the number of sequences written; records with unparseable
    headers are skipped.
    """
    written = 0
    with open(raw_path, encoding="latin-1") as src, open(dest_path, "w") as out:
        write_seq = False
        for line in src:
            if line.startswith(">"):
                rec = parse_scan_lib_header(line)
                if rec:
                    out.write(f">{rec['mernum']}|{rec['merops_id']}|{rec['subfamily']}\n")
                    written += 1
                    write_seq = True
                else:
                    write_seq = False
            elif write_seq:
                out.write(line)
    return written


def _download(url: str, dest: Path, label: str) -> None:
    print(f"Downloading {label} → {dest} ...")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as exc:
        print(f"ERROR: {label} download failed: {exc}", file=sys.stderr)
        sys.exit(1)


def build_merops_db(data_dir: Path, force: bool) -> tuple[Path, dict[str, dict]]:
    """Ensure ~/tools/MEROPS/DB/{merops.dmnd,family.txt,clan.txt} exist.

    Downloads the MEROPS scan library (1.9 MB, 5,009 representative
    peptidase-unit sequences — one per MEROPS identifier) plus the
    family→clan reference tables, rewrites the FASTA headers, and runs
    ``diamond makedb``. Skipped when everything already exists and ``force``
    is False. Returns (dmnd_path, family_meta).

    Requires: ``diamond`` on PATH.
    """
    db_dir = data_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)
    dmnd = db_dir / "merops.dmnd"
    family_txt = db_dir / "family.txt"
    clan_txt = db_dir / "clan.txt"

    if not (dmnd.exists() and family_txt.exists() and not force):
        if not shutil.which("diamond"):
            print("ERROR: diamond not on PATH", file=sys.stderr)
            sys.exit(1)

        raw_lib = db_dir / "merops_scan.raw.lib"
        clean_faa = db_dir / "merops.faa"

        _download(_SCAN_LIB_URL, raw_lib, "MEROPS scan library")
        _download(_FAMILY_TXT_URL, family_txt, "MEROPS family.txt")
        _download(_CLAN_TXT_URL, clan_txt, "MEROPS clan.txt")

        print(f"Rewriting headers → {clean_faa} ...")
        n_seqs = _rewrite_scan_lib_headers(raw_lib, clean_faa)
        print(f"  wrote {n_seqs} sequences")
        if n_seqs == 0:
            print("ERROR: scan library had no parseable sequences", file=sys.stderr)
            sys.exit(1)

        print(f"Building diamond DB → {dmnd} ...")
        cmd = ["diamond", "makedb", "--in", str(clean_faa), "-d", str(db_dir / "merops")]
        result = subprocess.run(cmd, cwd=str(REPO_ROOT))
        if result.returncode != 0 or not dmnd.exists():
            print("ERROR: diamond makedb failed", file=sys.stderr)
            sys.exit(1)

    family_meta = parse_family_txt(family_txt.read_text(encoding="latin-1"))
    return dmnd, family_meta


def diamond_version() -> str:
    """`diamond version` one-liner for provenance in skill_summary."""
    try:
        out = subprocess.run(["diamond", "version"], capture_output=True, text=True)
        return out.stdout.strip() or "unknown"
    except OSError:
        return "unknown"


def run_diamond(faa: Path, dmnd: Path, out_tsv: Path, threads: int, log_path: Path) -> bool:
    """Run diamond blastp for one strain. Returns True on success.

    Floor-only filtering at the diamond step: --evalue 0.001 only.
    Identity / coverage tiering happens in build_strain_calls (Python).
    Same parameters as tcdb-diamond — deliberate parity.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "diamond", "blastp",
        "-q", str(faa),
        "-d", str(dmnd),
        "-o", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident",
        "qcovhsp", "scovhsp", "length", "evalue", "bitscore",
        "--evalue", "0.001",
        "--max-target-seqs", "25",
        "--more-sensitive",
        "--threads", str(threads),
    ]
    print(f"\n>>> diamond blastp {faa.name} → {out_tsv} (log: {log_path.relative_to(REPO_ROOT)})")
    with open(log_path, "w") as logf:
        logf.write(f"$ {' '.join(cmd)}\n\n")
        logf.flush()
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
    return result.returncode == 0


def count_fasta_records(faa: Path) -> int:
    with open(faa) as f:
        return sum(1 for line in f if line.startswith(">"))


def truncate_faa(faa: Path, n_proteins: int, dest: Path) -> Path:
    """Copy the first N sequences of a FASTA to `dest` (--limit smoke tests)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(faa) as src, open(dest, "w") as out:
        seen = 0
        for line in src:
            if line.startswith(">"):
                seen += 1
                if seen > n_proteins:
                    break
            out.write(line)
    return dest


def process_strain(
    strain: str,
    data_dir_genome: Path,
    dmnd: Path,
    family_meta: dict[str, dict],
    version: str,
    threads: int,
    force: bool,
    limit: int | None = None,
) -> tuple[str, str, dict | None]:
    """Run the per-strain pipeline. Returns (strain, status, summary_or_None).

    Status values: OK / SKIP_NO_FAA / SKIP_EXISTS / FAILED_DIAMOND
    """
    faa = data_dir_genome / "protein.faa"
    if not faa.exists():
        return strain, "SKIP_NO_FAA", None

    suffix = f".limited_{limit}" if limit else ""
    out_tsv = data_dir_genome / TOOL / f"{strain}.{TOOL}{suffix}.tsv"

    if calls_path(data_dir_genome, TOOL, strain, limited=limit).exists() and not force:
        return strain, "SKIP_EXISTS", load_skill_summary(
            data_dir_genome, TOOL, strain, limited=limit
        )

    if limit:
        truncated = Path("/tmp") / f"merops_diamond_{strain}_first{limit}.faa"
        faa = truncate_faa(faa, limit, truncated)

    t0 = time.monotonic()
    log_path = LOGS_DIR / f"{strain}{suffix}.log"
    if not run_diamond(faa, dmnd, out_tsv, threads, log_path):
        print(f"  see log: {log_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        return strain, "FAILED_DIAMOND", None

    calls, summary = build_strain_calls(out_tsv, family_meta)
    summary = {
        "strain": strain,
        "tool_version": version,
        "input_proteins": count_fasta_records(faa),
        "wallclock_s": round(time.monotonic() - t0, 1),
        **summary,
    }

    save_calls(data_dir_genome, TOOL, strain, calls, limited=limit)
    save_skill_summary(data_dir_genome, TOOL, strain, summary, limited=limit)
    return strain, "OK", summary


def main():
    parser = argparse.ArgumentParser(
        description="Run diamond vs. the MEROPS scan library per strain (Phase 1)."
    )
    parser.add_argument("--strains", nargs="+", metavar="STRAIN",
                        help="Run only these strains (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if calls.json exists")
    parser.add_argument("--refresh-merops", action="store_true",
                        help="Re-download the MEROPS scan library + reference "
                             "tables and rebuild the diamond DB")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Diamond threads (default: os.cpu_count() or 4)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: run on first N proteins of each strain "
                             "only. Outputs go to <strain>.merops.limited_<N>.* "
                             "alongside (not replacing) full-run artifacts.")
    args = parser.parse_args()

    data_dir = resolve_merops_data_dir()
    print(f"MEROPS data dir: {data_dir}")
    dmnd, family_meta = build_merops_db(data_dir, force=args.refresh_merops)
    print(f"diamond DB: {dmnd} ({len(family_meta)} families in reference)")
    version = diamond_version()

    genomes = load_genome_rows(args.strains)

    results: list[tuple[str, str, dict | None]] = []
    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        results.append(
            process_strain(strain, data_dir_genome, dmnd, family_meta, version,
                           args.threads, args.force, limit=args.limit)
        )

    # Status table. Prot = proteins with >=1 call; Cand = total candidates.
    # T1/T2/T3, pep/inh, and hom (best_hit_kind = nonpeptidase_homolog) count
    # CANDIDATES, not proteins.
    print(f"\n{'='*104}")
    cols = ("Strain", "Status", "Prot", "Cand", "T1", "T2", "T3", "pep", "inh", "hom")
    print("{:<12} {:<14} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}".format(*cols))
    print("-" * 104)
    for strain, status, summary in results:
        if summary is None:
            print(f"{strain:<12} {status:<14}")
            continue
        td = summary["tier_distribution"]
        ed = summary.get("entry_type_distribution", {})
        bd = summary.get("best_hit_kind_distribution", {})
        print(
            "{:<12} {:<14} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}".format(
                strain, status, summary["proteins_with_call"],
                summary.get("total_candidates", 0),
                td.get("1", 0), td.get("2", 0), td.get("3", 0),
                ed.get("peptidase", 0), ed.get("inhibitor", 0),
                bd.get("nonpeptidase_homolog", 0),
            )
        )
    print("=" * 104)

    failed = [s for s, st, _ in results if st.startswith("FAILED")]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
