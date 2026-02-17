#!/bin/bash
set -euo pipefail

sleep 2
if [ ! -f /data/build2neo/neo4j-admin-import-call.sh ]; then
  echo "ERROR: /data/build2neo/neo4j-admin-import-call.sh not found. Build step may have failed." >&2
  exit 1
fi
chmod +x /data/build2neo/neo4j-admin-import-call.sh
/data/build2neo/neo4j-admin-import-call.sh
# Copy import report to accessible location
cp /var/lib/neo4j/import.report /data/build2neo/import.report 2>/dev/null || true
cp /var/lib/neo4j/import.report /output/import.report 2>/dev/null || true
neo4j start
sleep 10
neo4j stop
