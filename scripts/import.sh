#!/bin/bash
set -euo pipefail

sleep 2
if [ ! -f /data/build2neo/neo4j-admin-import-call.sh ]; then
  echo "ERROR: /data/build2neo/neo4j-admin-import-call.sh not found. Build step may have failed." >&2
  exit 1
fi
chmod +x /data/build2neo/neo4j-admin-import-call.sh
/data/build2neo/neo4j-admin-import-call.sh
neo4j start
sleep 10
neo4j stop
