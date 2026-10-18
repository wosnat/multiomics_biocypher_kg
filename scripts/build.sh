#!/bin/bash
set -euo pipefail

cd /usr/app/
# Skip biocypher-out/: local (non-Docker) build outputs, tens of GB, never read
# here (the Docker build writes to /usr/app/data/build2neo). Copying it filled
# the disk via the container's writable layer.
find /src -mindepth 1 -maxdepth 1 ! -name '.*' ! -name biocypher-out -exec cp -r {} . \;
cp /src/.env .
cp config/biocypher_docker_config.yaml config/biocypher_config.yaml
# Use persistent cache volume so pip/uv skip downloads on subsequent builds
pip install --cache-dir /cache/pip uv
export UV_CACHE_DIR=/cache/uv
uv sync
# add clean build (set via env variable in docker-compose.yml)
if [ "${BUILD2NEO_CLEANUP:-no}" == 'yes' ]; then
        echo "BUILD2NEO_CLEANUP is yes"
        echo "rm -rf /usr/app/data/build2neo"
        rm -rf /usr/app/data/build2neo
fi
EXTRA_ARGS=""
if [ "${BUILD_GO:-no}" == 'yes' ]; then
        EXTRA_ARGS="$EXTRA_ARGS --go"
fi
if [ "${BUILD_EC:-no}" == 'yes' ]; then
        EXTRA_ARGS="$EXTRA_ARGS --ec"
fi
uv run python3 create_knowledge_graph.py $EXTRA_ARGS
chmod -R 777 biocypher-log