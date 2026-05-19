#!/bin/bash
# Build the SignalP 6.0 (fast) Docker image as `signalp6`.
#
# Usage:
#   bash .claude/skills/signalp-run/build.sh
#
# Reads the licensed tarball (download from
# https://services.healthtech.dtu.dk/services/SignalP-6.0/; official install
# steps at https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md),
# stages the Docker build context under ~/tools/signalp-6.0i/, builds the
# image, and tags it `signalp6`. The image lives in Docker's daemon storage —
# nothing is written into the repo.
#
# Environment overrides:
#   SIGNALP_TARBALL    Path to the licensed tar.gz
#                      (default: ~/tools/signalp-6.0i.fast.tar.gz)
#   SIGNALP_BUILD_DIR  Build-context directory
#                      (default: ~/tools/signalp-6.0i)
#   SIGNALP_IMAGE_TAG  Docker tag for the built image (default: signalp6)
set -euo pipefail

TARBALL="${SIGNALP_TARBALL:-$HOME/tools/signalp-6.0i.fast.tar.gz}"
BUILD_DIR="${SIGNALP_BUILD_DIR:-$HOME/tools/signalp-6.0i}"
IMAGE_TAG="${SIGNALP_IMAGE_TAG:-signalp6}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [ ! -f "$TARBALL" ]; then
    echo "ERROR: SignalP tarball not found: $TARBALL"
    echo ""
    echo "To get it:"
    echo "  1. Register at https://services.healthtech.dtu.dk/services/SignalP-6.0/"
    echo "  2. Accept the academic license, download 'SignalP-6.0 fast'"
    echo "  3. Move the file to $TARBALL"
    echo "     (or set SIGNALP_TARBALL to your chosen path)"
    exit 1
fi

mkdir -p "$BUILD_DIR"
cp "$TARBALL" "$BUILD_DIR/signalp.tar.gz"
cp "$SCRIPT_DIR/Dockerfile" "$BUILD_DIR/Dockerfile"

echo "Building $IMAGE_TAG"
echo "  tarball:    $TARBALL"
echo "  context:    $BUILD_DIR"
echo "  Dockerfile: $SCRIPT_DIR/Dockerfile (copied into context)"
echo ""
echo "First build downloads Python 3.10 base + PyTorch <2.0 (~1 GB)."
echo ""

docker build -t "$IMAGE_TAG" "$BUILD_DIR"

echo ""
echo "Done! Image tagged '$IMAGE_TAG' (in Docker daemon storage; not in repo)."
echo ""
echo "Verify:"
echo "  docker run --rm $IMAGE_TAG --help"
echo ""
echo "Run on all strains:"
echo "  uv run python .claude/skills/signalp-run/run_signalp.py"
