# .claude/skills/tcdb-diamond/run_tcdb_diamond.py
"""Run diamond blastp vs. TCDB FASTA per strain. See spec
docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

# Ensure project root is on sys.path so `multiomics_kg` is importable when
# this script is run directly (uv run python .claude/skills/tcdb-diamond/...).
# Python adds the *script's* directory to sys.path, not the project root.
_REPO_ROOT_EARLY = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT_EARLY) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT_EARLY))

import dotenv

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.tcdb_diamond import build_strain_calls

REPO_ROOT = Path(__file__).resolve().parents[3]
LOGS_DIR = REPO_ROOT / "logs"
DEFAULT_TCDB_DATA_DIR = Path.home() / "tools" / "TCDB"

dotenv.load_dotenv(REPO_ROOT / ".env")


def resolve_tcdb_data_dir() -> Path:
    """Resolve the TCDB data dir (env override or default)."""
    env = os.environ.get("TCDB_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_TCDB_DATA_DIR


_TCDB_DOWNLOAD_URL = "https://tcdb.org/public/tcdb"
_TCDB_PFAM_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/pfam.py"
_TCDB_HEADER_RE = re.compile(r"^>gnl\|TC-DB\|([^|\s]+)\|(\S+)")


def _rewrite_tcdb_headers(raw_path: Path, dest_path: Path) -> int:
    """Rewrite raw TCDB FASTA headers to ``>lcl|<acc>-<tcid>`` form.

    Input headers (TCDB-canonical): ``>gnl|TC-DB|<accession>|<tcid> <description>``
    Output headers (parser-friendly): ``>lcl|<accession>-<tcid>``

    Returns the number of sequences written. Skips records with malformed
    headers (no parseable acc/tcid pair) — they would only cause spurious
    diamond hits anyway.
    """
    written = 0
    with open(raw_path) as src, open(dest_path, "w") as out:
        write_seq = False
        for line in src:
            if line.startswith(">"):
                m = _TCDB_HEADER_RE.match(line)
                if m:
                    acc, tcid = m.group(1), m.group(2)
                    out.write(f">lcl|{acc}-{tcid}\n")
                    written += 1
                    write_seq = True
                else:
                    write_seq = False
            elif write_seq:
                out.write(line)
    return written


def build_pfam_to_tc_map(data_dir: Path, force: bool) -> Path:
    """Ensure ~/tools/TCDB/DB/tcdb_pfam_map.tsv exists. Returns the path.

    Downloads TCDB's curated Pfam → TC mapping from
    https://www.tcdb.org/cgi-bin/projectv/public/pfam.py and writes a clean
    3-column TSV (PF_id, TC_id, family_name). The source URL embeds the table
    in HTML; we strip non-data lines to keep the cache parseable by
    `load_pfam_to_tc_map`. Skipped when the file already exists and ``force``
    is False. Failures are non-fatal — the diamond pipeline still runs without
    the Pfam tiebreaker, just with `pfam_agreement = 'neutral'` everywhere.
    """
    db_dir = data_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)
    out = db_dir / "tcdb_pfam_map.tsv"
    if out.exists() and not force:
        return out

    print(f"Downloading TCDB Pfam→TC mapping → {out} ...")
    try:
        with urllib.request.urlopen(_TCDB_PFAM_MAP_URL) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"WARN: Pfam map download failed: {exc} — continuing without "
              f"Pfam tiebreaker", file=sys.stderr)
        return out

    n = 0
    with open(out, "w") as fh:
        for line in text.splitlines():
            parts = line.rstrip().split("\t")
            if len(parts) < 3 or not parts[0].startswith("PF"):
                continue
            if not (len(parts[0]) >= 7 and parts[0][2:7].isdigit()):
                continue
            if parts[1].count(".") < 2:
                continue
            fh.write("\t".join(parts[:3]) + "\n")
            n += 1
    print(f"  wrote {n} Pfam→TC pairs")
    return out


def build_tcdb_db(data_dir: Path, force: bool) -> Path:
    """Ensure ~/tools/TCDB/DB/tcdb.dmnd exists. Returns the path.

    Downloads the curated TCDB FASTA from https://tcdb.org/public/tcdb,
    rewrites its headers from ``>gnl|TC-DB|<acc>|<tcid>`` to ``>lcl|<acc>-<tcid>``
    (the format consumed by parse_tcdb_subject_id), then runs ``diamond makedb``
    to build the diamond-format database. Skipped when the DB already exists
    and ``force`` is False.

    Note: Saier Lab's TCDBtools/bin/extractTCDB.pl performs the same operation
    but contains a bash-only redirect (``>&/dev/null``) that fails on dash-based
    /bin/sh systems (Ubuntu default). We replicate its behavior in pure Python
    to sidestep the shell-portability issue.

    Requires: ``diamond`` on PATH.
    """
    db_dir = data_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)
    dmnd = db_dir / "tcdb.dmnd"

    if dmnd.exists() and not force:
        return dmnd

    if not shutil.which("diamond"):
        print("ERROR: diamond not on PATH", file=sys.stderr)
        sys.exit(1)

    raw_faa = db_dir / "tcdb.raw.faa"
    clean_faa = db_dir / "tcdb.faa"

    print(f"Downloading TCDB FASTA → {raw_faa} ...")
    try:
        urllib.request.urlretrieve(_TCDB_DOWNLOAD_URL, raw_faa)
    except Exception as exc:
        print(f"ERROR: TCDB download failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Rewriting headers → {clean_faa} ...")
    n_seqs = _rewrite_tcdb_headers(raw_faa, clean_faa)
    print(f"  wrote {n_seqs} sequences")
    if n_seqs == 0:
        print("ERROR: TCDB FASTA had no parseable sequences", file=sys.stderr)
        sys.exit(1)

    print(f"Building diamond DB → {dmnd} ...")
    cmd = ["diamond", "makedb", "--in", str(clean_faa), "-d", str(db_dir / "tcdb")]
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    if result.returncode != 0:
        print(f"ERROR: diamond makedb returned {result.returncode}", file=sys.stderr)
        sys.exit(1)

    if not dmnd.exists():
        print(f"ERROR: diamond makedb completed but {dmnd} was not created", file=sys.stderr)
        sys.exit(1)

    return dmnd


def run_diamond(faa: Path, dmnd: Path, out_tsv: Path, threads: int, log_path: Path) -> bool:
    """Run diamond blastp for one strain. Returns True on success.

    Floor-only filtering at the diamond step: --evalue 0.001 only.
    Identity / coverage tiering happens in build_strain_calls (Python).

    Diamond's stdout + stderr are captured to `log_path` (overwritten on
    each run). The terminal only sees one progress line per strain.
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


def truncate_faa(faa: Path, n_proteins: int, dest: Path) -> Path:
    """Copy the first N sequences of a FASTA to `dest`. Returns dest.

    Used by the --limit flag for fast end-to-end smoke tests against a small
    subset of one strain's proteome (~10-30s instead of ~5min).
    """
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
    threads: int,
    force: bool,
    pfam_map: Path | None = None,
    limit: int | None = None,
) -> tuple[str, str, dict | None]:
    """Run the per-strain pipeline. Returns (strain, status, summary_or_None).

    When `limit` is set, runs diamond against the first `limit` proteins from
    protein.faa (truncated copy in /tmp), and writes outputs to
    `<strain>.tcdb.limited_<N>.{tsv,calls.json,skill_summary.json}` to avoid
    overwriting full-run artifacts. Skip-if-exists logic also keys on the
    limit-suffixed name so re-running with a different limit is a fresh run.

    Status values:
      OK / SKIP_NO_FAA / SKIP_EXISTS / FAILED_DIAMOND / FAILED_NO_EGGNOG
    """
    faa = data_dir_genome / "protein.faa"
    if not faa.exists():
        return strain, "SKIP_NO_FAA", None

    out_dir = data_dir_genome / "tcdb"
    suffix = f".limited_{limit}" if limit else ""
    out_tsv = out_dir / f"{strain}.tcdb{suffix}.tsv"
    out_calls = out_dir / f"{strain}.tcdb{suffix}.calls.json"
    out_summary = out_dir / f"{strain}.tcdb{suffix}.skill_summary.json"

    if out_calls.exists() and not force:
        with open(out_summary) as f:
            return strain, "SKIP_EXISTS", json.load(f)

    if limit:
        truncated = Path("/tmp") / f"tcdb_diamond_{strain}_first{limit}.faa"
        faa = truncate_faa(faa, limit, truncated)

    log_path = LOGS_DIR / f"tcdb_diamond_{strain}{suffix}.log"
    if not run_diamond(faa, dmnd, out_tsv, threads, log_path):
        print(f"  see log: {log_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        return strain, "FAILED_DIAMOND", None

    eggnog_path = data_dir_genome / "eggnog" / f"{strain}.emapper.annotations"
    if not eggnog_path.exists():
        # eggNOG missing → still proceed but agreement column will all be "extends"
        print(f"  WARN: {eggnog_path} missing — egn_agreement will all be 'extends'")

    gene_ann_path = data_dir_genome / "gene_annotations_merged.json"
    if not gene_ann_path.exists():
        # Gene annotations missing → pfam_agreement will all be 'neutral'
        print(f"  WARN: {gene_ann_path} missing — pfam_agreement will all be 'neutral'")
        gene_ann_arg = None
    else:
        gene_ann_arg = gene_ann_path

    calls, summary = build_strain_calls(out_tsv, eggnog_path, gene_ann_arg, pfam_map)

    with open(out_calls, "w") as f:
        json.dump(calls, f, indent=2, sort_keys=True)
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    return strain, "OK", summary


def main():
    parser = argparse.ArgumentParser(
        description="Run diamond vs. TCDB FASTA per strain (Phase 1)."
    )
    parser.add_argument("--strains", nargs="+", metavar="STRAIN",
                        help="Run only these strains (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if calls.json exists")
    parser.add_argument("--refresh-tcdb", action="store_true",
                        help="Re-download TCDB FASTA + diamond DB even if cached")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Diamond threads (default: os.cpu_count() or 4)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: run on first N proteins of each strain "
                             "only. Outputs go to <strain>.tcdb.limited_<N>.* "
                             "alongside (not replacing) full-run artifacts.")
    args = parser.parse_args()

    data_dir = resolve_tcdb_data_dir()
    print(f"TCDB data dir: {data_dir}")
    dmnd = build_tcdb_db(data_dir, force=args.refresh_tcdb)
    print(f"diamond DB: {dmnd}")
    pfam_map = build_pfam_to_tc_map(data_dir, force=args.refresh_tcdb)
    print(f"Pfam→TC map: {pfam_map}{' (missing)' if not pfam_map.exists() else ''}")

    # load_genome_rows exits with an error if --strains names a strain not in
    # the registry; no need to handle that here.
    genomes = load_genome_rows(args.strains)

    results: list[tuple[str, str, dict | None]] = []
    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        results.append(
            process_strain(strain, data_dir_genome, dmnd, args.threads, args.force,
                           pfam_map=pfam_map if pfam_map.exists() else None,
                           limit=args.limit)
        )

    # Status table. Prot = proteins with >=1 call; Cand = total candidates
    # (proteins with hits in multiple families contribute >1 candidate).
    # T1/T2/T3 and agreement columns count CANDIDATES, not proteins.
    # Kept = candidates passing the filter; drp_* = candidates dropped by
    # each filter rule.
    print(f"\n{'='*160}")
    cols = ("Strain", "Status", "Prot", "Cand", "Kept",
            "T1", "T2", "T3",
            "confirms", "refines", "extends", "conflicts",
            "pfam_d", "pfam_e", "pfam_b",
            "drp_pf", "drp_eg", "drp_sng", "drp_lc")
    print("{:<12} {:<14} {:>5} {:>5} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7}".format(*cols))
    print("-" * 168)
    for strain, status, summary in results:
        if summary is None:
            print(f"{strain:<12} {status:<14}")
            continue
        td = summary["tier_distribution"]
        ad = summary["agreement_distribution"]
        pad = summary.get("pfam_agreement_distribution", {})
        fad = summary.get("filter_action_distribution", {})
        print(
            "{:<12} {:<14} {:>5} {:>5} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7}".format(
                strain, status, summary["proteins_with_call"],
                summary.get("total_candidates", 0),
                fad.get("keep", 0),
                td.get("1", 0), td.get("2", 0), td.get("3", 0),
                ad.get("confirms", 0), ad.get("refines", 0),
                ad.get("extends", 0), ad.get("conflicts", 0),
                pad.get("confirms_diamond", 0), pad.get("confirms_eggnog", 0),
                pad.get("confirms_both", 0),
                fad.get("drop_pfam_contradicts", 0),
                fad.get("drop_egn_conflicts", 0),
                fad.get("drop_singleton_low_score", 0),
                fad.get("drop_low_confidence", 0),
            )
        )
    print("=" * 168)

    failed = [s for s, st, _ in results if st.startswith("FAILED")]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
