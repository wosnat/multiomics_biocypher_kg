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
# Git provenance of the checkout this build ran from, for Schema_info.git_* on a
# local build (/release-kg sets KG_GIT_* itself, which take precedence in
# post-import.sh). Dirty = modified tracked files; untracked scratch dirs ignored.
if git config --global --add safe.directory /src 2>/dev/null && git -C /src rev-parse HEAD >/dev/null 2>&1; then
        {
                echo "git_sha=$(git -C /src rev-parse HEAD)"
                echo "git_sha_short=$(git -C /src rev-parse --short HEAD)"
                echo "git_branch=$(git -C /src rev-parse --abbrev-ref HEAD)"
                if [ -n "$(git -C /src status --porcelain --untracked-files=no)" ]; then
                        echo "git_dirty=true"
                else
                        echo "git_dirty=false"
                fi
        } > data/build2neo/git_provenance.txt
        cat data/build2neo/git_provenance.txt
else
        echo "WARNING: /src is not a git checkout; Schema_info.git_* will read 'unknown'"
fi
chmod -R 777 biocypher-log