#!/bin/bash -c
cd /usr/app/
cp -r /src/* .
cp /src/.env .
cp config/biocypher_docker_config.yaml config/biocypher_config.yaml
pip install uv
uv sync
# add clean build (set via env variable in docker-compose.yml)
if [ "${BUILD2NEO_CLEANUP}" == 'yes' ]; then
        echo "BUILD2NEO_CLEANUP is yes"
        echo "rm -rf /usr/app/data/build2neo"
        rm -rf /usr/app/data/build2neo
fi
uv run python3 create_knowledge_graph.py
chmod -R 777 biocypher-log