#!/usr/bin/env bash
# Refresh MNX (MetaNetX) cross-reference cache + SQLite resolver.
#
# Run this ONLY when MNX releases a new version. The resolver is heavy (~2.5 GB
# SQLite, ~30 minutes wall time) and does not depend on which strains are in
# the project — it's a global metabolism cross-reference lookup.
#
# After this completes, run `bash scripts/prepare_data.sh --steps 6 --force`
# to regenerate cache/data/kegg/kegg_data.json with the new MNX cross-refs.
#
# Outputs land at $MNX_DATA_DIR (set in .env, e.g. ~/tools/mnx) — defaults to
# cache/data/mnx/ when unset. Sharing one MNX_DATA_DIR across checkouts avoids
# the ~4 GB duplication and the ~30 min rebuild cost:
#   $MNX_DATA_DIR/{chem_prop,chem_xref,reac_prop,reac_xref}.tsv  (~1.5 GB raw)
#   $MNX_DATA_DIR/metabolite_resolver.db                         (~2.5 GB SQLite)
#
# Usage:
#   bash scripts/refresh_mnx.sh           # download + build (skip if cached)
#   bash scripts/refresh_mnx.sh --force   # re-download + rebuild

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$LOG_DIR"

FORCE=""
if [[ "${1:-}" == "--force" ]]; then
    FORCE="--force"
fi

cd "$PROJECT_ROOT"

echo "══════════════════════════════════════════════════════════════════════"
echo "  refresh_mnx.sh — MNX TSV download + SQLite resolver build"
echo "  Log: $LOG_DIR/refresh_mnx.log"
echo "══════════════════════════════════════════════════════════════════════"

{
    echo ""
    echo "── 1. Download MNX TSVs (4 files, ~1.5 GB) ──"
    uv run python -m multiomics_kg.download.download_metabolism_reference \
        --sources mnx $FORCE

    echo ""
    echo "── 2. Build metabolite resolver SQLite (~30 min, ~2.5 GB output) ──"
    uv run python -m multiomics_kg.download.build_mnx_resolver $FORCE
} 2>&1 | tee "$LOG_DIR/refresh_mnx.log"

echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo "  refresh_mnx.sh complete. Next: bash scripts/prepare_data.sh --steps 6 --force"
echo "══════════════════════════════════════════════════════════════════════"
