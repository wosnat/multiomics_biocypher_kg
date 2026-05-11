# .claude/skills/psortb-run/run_psortb.py
"""Run PSORTb v3.0.3 (Gram-negative) per strain via Docker. See spec
docs/superpowers/specs/2026-05-10-psortb-localization-design.md.

All KG strains are Gram-negative, so --negative is hardcoded. Phase 1 emits
inspectable per-strain JSON artifacts; KG integration is deferred to Phase 2.
"""
from __future__ import annotations

import argparse
import datetime
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Ensure project root is on sys.path so dotenv works the same as in other skills.
_REPO_ROOT_EARLY = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT_EARLY) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT_EARLY))

import dotenv

from multiomics_kg.download.utils.cli import load_genome_rows

REPO_ROOT = Path(__file__).resolve().parents[3]
LOGS_DIR = REPO_ROOT / "logs"

PSORTB_IMAGE = "brinkmanlab/psortb_commandline:1.0.2"
# psort wrapper path inside the brinkmanlab image. Verified by inspecting the
# image's PATH at /usr/local/psortb/bin/psort.
PSORTB_BINARY = "/usr/local/psortb/bin/psort"

VALID_LOCALIZATIONS = {
    "Cytoplasmic",
    "CytoplasmicMembrane",
    "Periplasmic",
    "OuterMembrane",
    "Extracellular",
    "Unknown",
}

dotenv.load_dotenv(REPO_ROOT / ".env")


def check_docker_on_path() -> None:
    """Exit 2 with install hint if docker is not on PATH."""
    if not shutil.which("docker"):
        print(
            "ERROR: docker is not on PATH. Install Docker Engine "
            "(https://docs.docker.com/engine/install/) and re-run.",
            file=sys.stderr,
        )
        sys.exit(2)


def image_present_locally() -> bool:
    """Return True iff PSORTb image is already in local Docker daemon storage."""
    result = subprocess.run(
        ["docker", "image", "inspect", PSORTB_IMAGE],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def get_image_digest() -> str | None:
    """Return the RepoDigest for the locally-cached PSORTb image, or None."""
    result = subprocess.run(
        ["docker", "inspect", "--format={{json .RepoDigests}}", PSORTB_IMAGE],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    try:
        digests = json.loads(result.stdout.strip())
    except json.JSONDecodeError:
        return None
    if digests:
        return digests[0]
    return None


def prepare_image(force_refresh: bool = False) -> int:
    """Pull the PSORTb image into Docker's daemon storage, print digest."""
    if image_present_locally() and not force_refresh:
        digest = get_image_digest()
        print(f"PSORTb image already present: {PSORTB_IMAGE}")
        if digest:
            print(f"  digest: {digest}")
        print("Re-run with --refresh-image to force a re-pull.")
        return 0

    action = "Re-pulling" if force_refresh else "Pulling"
    print(f"{action} {PSORTB_IMAGE} (~2 GB)...")
    result = subprocess.run(["docker", "pull", PSORTB_IMAGE])
    if result.returncode != 0:
        print(f"ERROR: docker pull failed (exit {result.returncode}).", file=sys.stderr)
        return 1

    digest = get_image_digest()
    print(f"\nImage ready: {PSORTB_IMAGE}")
    if digest:
        print(f"  digest: {digest}")
    return 0


def truncate_faa(faa: Path, n_proteins: int, dest: Path) -> Path:
    """Copy the first N sequences of a FASTA to `dest`. Returns dest."""
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


def parse_terse_line(line: str) -> tuple[str, str, float, str | None, float | None, bool] | None:
    """Parse one line of PSORTb terse TSV output.

    Returns (seq_id, primary_loc, primary_score, secondary_loc, secondary_score,
    is_multi) or None for malformed/header lines.
    """
    line = line.rstrip("\n")
    if not line:
        return None
    parts = line.split("\t")
    if len(parts) < 3:
        return None
    seq_id_raw, loc_field, score_field = parts[0], parts[1], parts[2]
    # Skip header line if PSORTb emits one.
    if seq_id_raw.lower() in {"seqid", "seq_id", "sequence_id"}:
        return None

    # Strip lcl| prefix if present and take the first whitespace-delimited token.
    seq_id = seq_id_raw.split()[0]
    if seq_id.startswith("lcl|"):
        seq_id = seq_id[len("lcl|"):]

    loc_field = loc_field.strip()
    score_field = score_field.strip()

    loc_parts = [p for p in loc_field.split("/") if p]
    score_parts = [p for p in score_field.split("/") if p]

    if not loc_parts:
        return None

    primary_loc = loc_parts[0]
    try:
        primary_score = float(score_parts[0]) if score_parts else 0.0
    except ValueError:
        primary_score = 0.0

    secondary_loc: str | None = None
    secondary_score: float | None = None
    is_multi = False
    if len(loc_parts) >= 2:
        is_multi = True
        # Order primary/secondary by descending score to match the spec.
        try:
            paired = [
                (loc_parts[i], float(score_parts[i]) if i < len(score_parts) else 0.0)
                for i in range(len(loc_parts))
            ]
        except ValueError:
            paired = [(loc_parts[i], 0.0) for i in range(len(loc_parts))]
        paired.sort(key=lambda x: x[1], reverse=True)
        primary_loc, primary_score = paired[0]
        secondary_loc, secondary_score = paired[1]

    return seq_id, primary_loc, primary_score, secondary_loc, secondary_score, is_multi


def parse_terse_tsv(tsv_path: Path) -> dict[str, dict]:
    """Parse the full terse TSV into per-protein records."""
    calls: dict[str, dict] = {}
    with open(tsv_path) as f:
        for line in f:
            parsed = parse_terse_line(line)
            if parsed is None:
                continue
            seq_id, prim_loc, prim_score, sec_loc, sec_score, is_multi = parsed
            calls[seq_id] = {
                "localization": prim_loc,
                "score": prim_score,
                "secondary_localization": sec_loc,
                "secondary_score": sec_score,
                "is_multi_localized": is_multi,
                "is_unknown": prim_loc == "Unknown",
            }
    return calls


def summarize_calls(calls: dict[str, dict]) -> dict:
    """Build per-strain summary stats from a calls dict."""
    counts: dict[str, int] = {}
    n_unknown = 0
    n_multi = 0
    for rec in calls.values():
        counts[rec["localization"]] = counts.get(rec["localization"], 0) + 1
        if rec["is_unknown"]:
            n_unknown += 1
        if rec["is_multi_localized"]:
            n_multi += 1
    total = len(calls)
    return {
        "protein_count": total,
        "localization_counts": dict(sorted(counts.items())),
        "unknown_rate": (n_unknown / total) if total else 0.0,
        "multi_localized_rate": (n_multi / total) if total else 0.0,
    }


def run_psortb_container(
    data_dir_genome: Path,
    fasta_basename: str,
    out_tsv: Path,
    log_path: Path,
) -> bool:
    """Run the brinkmanlab PSORTb container on one strain's FASTA.

    Mounts `data_dir_genome` at /tmp/psortb in the container, runs psort on
    the named FASTA file, captures stdout (terse TSV) to `out_tsv` and stderr
    to `log_path`. Returns True on exit-0.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{data_dir_genome.resolve()}:/tmp/psortb",
        PSORTB_IMAGE,
        PSORTB_BINARY,
        "--negative",
        "--output", "terse",
        f"/tmp/psortb/{fasta_basename}",
    ]
    print(f"\n>>> psortb {fasta_basename} → {out_tsv.relative_to(REPO_ROOT)} "
          f"(log: {log_path.relative_to(REPO_ROOT)})")
    with open(out_tsv, "w") as outf, open(log_path, "w") as logf:
        logf.write(f"$ {' '.join(cmd)}\n\n")
        logf.flush()
        result = subprocess.run(cmd, stdout=outf, stderr=logf)
    return result.returncode == 0


def process_strain(
    strain: str,
    data_dir_genome: Path,
    image_digest: str | None,
    force: bool,
    limit: int | None,
) -> tuple[str, str, dict | None]:
    """Per-strain pipeline. Returns (strain, status, summary_or_None).

    Status values:
      OK / SKIP_NO_FAA / SKIP_EXISTS / FAILED_DOCKER / FAILED_NO_OUTPUT
    """
    faa = data_dir_genome / "protein.faa"
    if not faa.exists():
        return strain, "MISSING_INPUT", None

    out_dir = data_dir_genome / "psortb"
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = f".limited_{limit}" if limit else ""
    out_tsv = out_dir / f"{strain}.psortb{suffix}.terse.tsv"
    out_calls = out_dir / f"{strain}.psortb{suffix}.calls.json"
    out_summary = out_dir / f"{strain}.psortb{suffix}.skill_summary.json"

    if out_calls.exists() and not force:
        try:
            with open(out_summary) as f:
                return strain, "SKIP_EXISTS", json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return strain, "SKIP_EXISTS", None

    # Decide which FASTA to feed the container. For --limit, write a truncated
    # copy alongside protein.faa so the volume mount picks it up.
    if limit:
        truncated_name = f"protein.limited_{limit}.faa"
        truncated_path = data_dir_genome / truncated_name
        truncate_faa(faa, limit, truncated_path)
        fasta_basename = truncated_name
    else:
        fasta_basename = "protein.faa"

    log_path = LOGS_DIR / f"psortb_{strain}{suffix}.log"

    started_at = datetime.datetime.now(datetime.timezone.utc)
    t0 = time.monotonic()
    ok = run_psortb_container(data_dir_genome, fasta_basename, out_tsv, log_path)
    elapsed = time.monotonic() - t0
    finished_at = datetime.datetime.now(datetime.timezone.utc)

    # Clean up the truncated FASTA copy used for --limit smoke tests.
    if limit:
        try:
            (data_dir_genome / fasta_basename).unlink()
        except FileNotFoundError:
            pass

    if not ok:
        print(f"  see log: {log_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        return strain, "FAILED_DOCKER", None

    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return strain, "FAILED_NO_OUTPUT", None

    calls = parse_terse_tsv(out_tsv)
    summary = summarize_calls(calls)
    summary.update({
        "strain": strain,
        "wallclock_seconds": round(elapsed, 2),
        "image_digest": image_digest,
        "psortb_args": ["--negative", "--output", "terse"],
        "started_at": started_at.isoformat(timespec="seconds"),
        "finished_at": finished_at.isoformat(timespec="seconds"),
    })

    with open(out_calls, "w") as f:
        json.dump(calls, f, indent=2, sort_keys=True)
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    return strain, "OK", summary


def print_status_table(results: list[tuple[str, str, dict | None]]) -> None:
    """Print a per-strain status table summarizing the run."""
    print(f"\n{'=' * 110}")
    cols = ("Strain", "Status", "N", "Cyto", "CytoMem", "Peri", "OutMem", "Extra", "Unk", "Multi")
    print("{:<12} {:<16} {:>6} {:>6} {:>8} {:>6} {:>7} {:>6} {:>5} {:>6}".format(*cols))
    print("-" * 110)
    for strain, status, summary in results:
        if summary is None:
            print(f"{strain:<12} {status:<16}")
            continue
        counts = summary.get("localization_counts", {})
        n_multi = round(summary.get("multi_localized_rate", 0.0) * summary.get("protein_count", 0))
        print(
            "{:<12} {:<16} {:>6} {:>6} {:>8} {:>6} {:>7} {:>6} {:>5} {:>6}".format(
                strain, status,
                summary.get("protein_count", 0),
                counts.get("Cytoplasmic", 0),
                counts.get("CytoplasmicMembrane", 0),
                counts.get("Periplasmic", 0),
                counts.get("OuterMembrane", 0),
                counts.get("Extracellular", 0),
                counts.get("Unknown", 0),
                n_multi,
            )
        )
    print("=" * 110)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run PSORTb v3.0.3 per strain via Docker (Phase 1)."
    )
    parser.add_argument("--strain", help="Run only this strain (e.g. MED4)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if calls.json exists")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: run on first N proteins of each "
                             "strain only. Outputs go to "
                             "<strain>.psortb.limited_<N>.* alongside "
                             "(not replacing) full-run artifacts.")
    img_group = parser.add_mutually_exclusive_group()
    img_group.add_argument("--prepare-image", action="store_true",
                           help="Pull the PSORTb Docker image (~2 GB) and exit.")
    img_group.add_argument("--refresh-image", action="store_true",
                           help="Re-pull the PSORTb Docker image even if present, then exit.")
    args = parser.parse_args()

    check_docker_on_path()

    if args.prepare_image:
        return prepare_image(force_refresh=False)
    if args.refresh_image:
        return prepare_image(force_refresh=True)

    if not image_present_locally():
        print(
            f"ERROR: PSORTb image {PSORTB_IMAGE} is not present locally.\n"
            "Run --prepare-image first (one-time ~2 GB download):\n"
            "  uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image",
            file=sys.stderr,
        )
        return 2

    image_digest = get_image_digest()

    # Run on any row with a strain_name + data_dir, regardless of organism_type.
    # protein.faa presence is checked per-strain inside process_strain().
    genomes = load_genome_rows([args.strain] if args.strain else None)
    if not genomes:
        print(
            f"No genome rows found"
            f"{f' for strain {args.strain}' if args.strain else ''}.",
            file=sys.stderr,
        )
        return 1

    print(f"PSORTb image: {PSORTB_IMAGE}")
    if image_digest:
        print(f"  digest: {image_digest}")
    print(f"Strains to process: {len(genomes)}")

    results: list[tuple[str, str, dict | None]] = []
    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        results.append(
            process_strain(strain, data_dir_genome, image_digest, args.force, args.limit)
        )

    print_status_table(results)

    failed = [s for s, st, _ in results if st.startswith("FAILED")]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
