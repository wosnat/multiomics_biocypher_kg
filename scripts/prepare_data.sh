#!/usr/bin/env bash
# Orchestration: download + preprocess genome annotation data.
#
# Step 0 — Download genome data (NCBI GFF, Cyanorak GFF/GBK, UniProt, gene_mapping.csv)
#           calls: multiomics_kg/download/download_genome_data.py --steps 1 2 3 5
#           (eggNOG step 4 is intentionally skipped — run /eggnog-run skill separately)
#           Use --skip-cyanorak to skip Cyanorak downloads (step 2) when server is slow
# Step 1 — Build per-strain gene annotation tables (gene_annotations_merged.json)
#           calls: multiomics_kg/download/build_gene_annotations.py
# Step 2 — Build per-taxid protein annotation tables (protein_annotations.json)
#           calls: multiomics_kg/download/build_protein_annotations.py
#           Requires step 0 (UniProt data must be cached first)
#
# Logs: logs/prepare_data_step0.log, logs/prepare_data_step1.log, logs/prepare_data_step2.log
#       Monitor with: tail -f logs/prepare_data_step0.log
#
# Usage:
#   ./scripts/prepare_data.sh
#   ./scripts/prepare_data.sh --force
#   ./scripts/prepare_data.sh --skip-cyanorak --force   # skip slow Cyanorak server
#   ./scripts/prepare_data.sh --strains MED4 MIT9313
#   ./scripts/prepare_data.sh --steps 0 1 2
#   ./scripts/prepare_data.sh --steps 0 --force
#   ./scripts/prepare_data.sh --steps 2 --force         # rebuild protein_annotations.json only
#   ./scripts/prepare_data.sh --steps 1 --strains MED4 --force

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$LOG_DIR"

# ── parse args ────────────────────────────────────────────────────────────────

FORCE=""
STEPS="0 1 2"
STRAINS=()
SKIP_CYANORAK=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)          FORCE="--force"; shift ;;
        --skip-cyanorak)  SKIP_CYANORAK=1; shift ;;
        --steps)
            STEPS=""
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STEPS="$STEPS $1"; shift
            done
            ;;
        --strains)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STRAINS+=("$1"); shift
            done
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

STRAINS_ARG=""
if [[ ${#STRAINS[@]} -gt 0 ]]; then
    STRAINS_ARG="--strains ${STRAINS[*]}"
fi

# ── helpers ───────────────────────────────────────────────────────────────────

run_step() {
    local step_num="$1"
    local label="$2"
    local log="$3"
    shift 3
    local cmd=("$@")

    echo ""
    echo "══════════════════════════════════════════════════════════════════════"
    echo "  Step $step_num — $label"
    echo "  Log: $log"
    echo "  Cmd: ${cmd[*]}"
    echo "══════════════════════════════════════════════════════════════════════"
    echo ""

    # tee: show output live AND write to log
    if ! "${cmd[@]}" 2>&1 | tee "$log"; then
        echo ""
        echo "ERROR: Step $step_num failed — check $log" >&2
        exit 1
    fi
    echo ""
    echo "  Step $step_num complete. Full log: $log"
}

# ── steps ─────────────────────────────────────────────────────────────────────

cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

echo "prepare_data.sh: steps=[${STEPS}]${STRAINS_ARG:+ strains=[${STRAINS[*]}]}${FORCE:+ (force)}${SKIP_CYANORAK:+ (skip-cyanorak)}"
echo "(step 1 = gene annotations, step 2 = protein annotations)"
echo "Project root: $PROJECT_ROOT"
echo "Logs dir:     $LOG_DIR"

for step in $STEPS; do
    case "$step" in
        0)
            if [[ $SKIP_CYANORAK -eq 1 ]]; then
                DOWNLOAD_SUBSTEPS="1 3 5"
                STEP0_LABEL="Download genome data (NCBI + UniProt + gene_mapping; Cyanorak SKIPPED)"
            else
                DOWNLOAD_SUBSTEPS="1 2 3 5"
                STEP0_LABEL="Download genome data (NCBI + Cyanorak + UniProt + gene_mapping)"
            fi
            run_step 0 \
                "$STEP0_LABEL" \
                "$LOG_DIR/prepare_data_step0.log" \
                uv run python multiomics_kg/download/download_genome_data.py \
                    --steps $DOWNLOAD_SUBSTEPS \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        1)
            run_step 1 \
                "Build gene annotation tables (gene_annotations_merged.json)" \
                "$LOG_DIR/prepare_data_step1.log" \
                uv run python multiomics_kg/download/build_gene_annotations.py \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        2)
            run_step 2 \
                "Build protein annotation tables (protein_annotations.json)" \
                "$LOG_DIR/prepare_data_step2.log" \
                uv run python multiomics_kg/download/build_protein_annotations.py \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        *)
            echo "Unknown step: $step (valid: 0 1 2)" >&2
            exit 1
            ;;
    esac
done

echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo "  All steps complete."
echo "══════════════════════════════════════════════════════════════════════"
