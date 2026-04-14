#!/bin/bash
set -euo pipefail

OUTPUT_DIR=/output
REPORT_SRC=/var/lib/neo4j/import.report
VOLUME_REPORT=/data/build2neo/import.report

mkdir -p "$OUTPUT_DIR"
# The host-side ./output/ mount may be root-owned from prior runs; make it
# writable so cp below can always succeed.
chmod 777 "$OUTPUT_DIR" || true

IMPORT_EXIT=1
write_status() {
  local trap_exit=$?
  echo "exit_code=${IMPORT_EXIT}" > "$OUTPUT_DIR/import.status" || true
  if [ -f "$REPORT_SRC" ]; then
    cp "$REPORT_SRC" "$OUTPUT_DIR/import.report" || true
    cp "$REPORT_SRC" "$VOLUME_REPORT" || true
  else
    : > "$OUTPUT_DIR/import.report" || true
  fi
  exit "$trap_exit"
}
trap write_status EXIT

sleep 2
if [ ! -f /data/build2neo/neo4j-admin-import-call.sh ]; then
  echo "ERROR: /data/build2neo/neo4j-admin-import-call.sh not found. Build step may have failed." >&2
  exit 1
fi
chmod +x /data/build2neo/neo4j-admin-import-call.sh

# Run the import without letting `set -e` abort before we capture the exit code.
set +e
/data/build2neo/neo4j-admin-import-call.sh
IMPORT_EXIT=$?
set -e

if [ "$IMPORT_EXIT" -ne 0 ]; then
  echo "neo4j-admin import failed with exit code $IMPORT_EXIT" >&2
  exit "$IMPORT_EXIT"
fi

neo4j start
sleep 10
neo4j stop
