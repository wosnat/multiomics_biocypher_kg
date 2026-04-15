#!/usr/bin/env bash
# Overnight orchestrator for eggNOG on 4 new strains.
# Polls for each strain's protein.faa, then calls run_eggnog.py --strain <name>.
# Writes progress to logs/eggnog_overnight_2026-04-14.md.
# This helper script lives alongside the log so the background agent's
# writes stay scoped to logs/ and cache/ eggnog output directories.

set -u
REPO=/home/osnat/github/multiomics_biocypher_kg
LOG="$REPO/logs/eggnog_overnight_2026-04-14.md"
START_EPOCH=$(date +%s)
MAX_SECONDS=$((8 * 3600))  # 8 hours

# strain:org_dir:cache_rel
STRAINS=(
  "SS120:Prochlorococcus:cache/data/Prochlorococcus/genomes/SS120"
  "BL107:Synechococcus:cache/data/Synechococcus/genomes/BL107"
  "HP15:Marinobacter:cache/data/Marinobacter/genomes/HP15"
  "AltMedDE:Alteromonas:cache/data/Alteromonas/genomes/AltMedDE"
)

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

append_event() {
  echo "- $(timestamp) | $1" >> "$LOG"
}

update_row() {
  local strain="$1" fasta="$2" started="$3" finished="$4" status="$5"
  python3 - "$LOG" "$strain" "$fasta" "$started" "$finished" "$status" <<'PY'
import sys, re
path, strain, fasta, started, finished, status = sys.argv[1:]
txt = open(path).read()
pattern = re.compile(rf"^\| {re.escape(strain)} \|.*$", re.MULTILINE)
new = f"| {strain} | {fasta} | {started} | {finished} | {status} |"
txt2 = pattern.sub(new, txt, count=1)
open(path, "w").write(txt2)
PY
}

elapsed() { echo $(( $(date +%s) - START_EPOCH )); }

append_event "orchestrator started (PID $$)"

for entry in "${STRAINS[@]}"; do
  IFS=: read -r STRAIN ORG CACHE_REL <<< "$entry"
  CACHE_DIR="$REPO/$CACHE_REL"
  FAA="$CACHE_DIR/protein.faa"
  OUT_ANNOT="$CACHE_DIR/eggnog/${STRAIN}.emapper.annotations"

  append_event "[$STRAIN] waiting for $FAA"
  update_row "$STRAIN" "-" "-" "-" "waiting for FASTA"

  # Poll for FASTA
  while [ ! -s "$FAA" ]; do
    if [ "$(elapsed)" -gt "$MAX_SECONDS" ]; then
      append_event "[$STRAIN] global timeout waiting for FASTA; skipping remaining"
      update_row "$STRAIN" "not found" "-" "-" "TIMEOUT (no FASTA after 8h)"
      continue 2
    fi
    sleep 60
  done

  FASTA_TS="$(timestamp)"
  LINES=$(wc -l < "$FAA" 2>/dev/null || echo "?")
  append_event "[$STRAIN] FASTA detected at $FAA ($LINES lines)"
  update_row "$STRAIN" "$FAA" "-" "-" "FASTA found"

  # run_eggnog.py reads cyanobacteria_genomes.csv; wait briefly if row missing
  CSV="$REPO/data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
  if ! grep -E "^[^#]*,${STRAIN}," "$CSV" >/dev/null 2>&1; then
    append_event "[$STRAIN] not in cyanobacteria_genomes.csv yet; waiting up to 10 min"
    CSV_WAIT=0
    while [ $CSV_WAIT -lt 600 ]; do
      if grep -E "^[^#]*,${STRAIN}," "$CSV" >/dev/null 2>&1; then
        append_event "[$STRAIN] CSV row appeared after ${CSV_WAIT}s"
        break
      fi
      sleep 30
      CSV_WAIT=$((CSV_WAIT + 30))
    done
    if ! grep -E "^[^#]*,${STRAIN}," "$CSV" >/dev/null 2>&1; then
      append_event "[$STRAIN] CSV row never appeared; run_eggnog.py cannot run -- skipping"
      update_row "$STRAIN" "$FAA" "-" "-" "FAILED: not in genomes.csv"
      continue
    fi
  fi

  START_TS="$(timestamp)"
  update_row "$STRAIN" "$FAA" "$START_TS" "-" "running"
  append_event "[$STRAIN] invoking run_eggnog.py --strain $STRAIN --cpu 8"

  STRAIN_LOG="$REPO/logs/eggnog_overnight_${STRAIN}.log"
  (cd "$REPO" && uv run python .claude/skills/eggnog-run/run_eggnog.py --strain "$STRAIN" --cpu 8) \
    > "$STRAIN_LOG" 2>&1
  RC=$?

  END_TS="$(timestamp)"
  if [ $RC -eq 0 ] && [ -s "$OUT_ANNOT" ]; then
    N=$(grep -cv '^#' "$OUT_ANNOT" 2>/dev/null || echo "?")
    update_row "$STRAIN" "$FAA" "$START_TS" "$END_TS" "OK ($N annotations)"
    append_event "[$STRAIN] completed OK; $N annotations in $OUT_ANNOT"
  else
    update_row "$STRAIN" "$FAA" "$START_TS" "$END_TS" "FAILED (rc=$RC, see $STRAIN_LOG)"
    append_event "[$STRAIN] FAILED rc=$RC; log at $STRAIN_LOG"
  fi
done

append_event "orchestrator finished after $(elapsed)s"
