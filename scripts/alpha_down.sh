#!/usr/bin/env bash
# Operator wrapper: tear down the alpha stack.
#
# Usage:
#   ./scripts/alpha_down.sh         # stop alpha-deploy, KEEP volumes (the
#                                   # active-color data volume is the only
#                                   # copy of the live KG — do not casually
#                                   # delete it)
#   ./scripts/alpha_down.sh -v      # also delete the project's anonymous
#                                   # volumes (does NOT touch the named
#                                   # kg-alpha-blue / kg-alpha-green
#                                   # external volumes — those persist by
#                                   # name across this command)
#
# To delete an alpha color volume explicitly (e.g. recovering disk), use:
#   docker volume rm kg-alpha-blue   # or kg-alpha-green
# but only AFTER confirming alpha-deploy is down and the other color is the
# one in the active-color marker.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ ! -f .env.alpha ]]; then
  echo "FATAL: .env.alpha is missing." >&2
  exit 1
fi

set -a; . ./.env.alpha; set +a

# Default to blue if the marker is missing — `down` is destructive against
# state, so we want to read whichever color compose registered with the
# kg-alpha project. The marker is the source of truth post-release.
COLOR_MARKER=".alpha_active_color"
ACTIVE_COLOR=$(cat "$COLOR_MARKER" 2>/dev/null || echo "blue")

export ALPHA_DATA_VOLUME="kg-alpha-${ACTIVE_COLOR}"
export KG_DEPLOY_HTTP_BIND="${ALPHA_BIND_IP}:17474"
export KG_DEPLOY_BOLT_BIND="${ALPHA_BIND_IP}:17687"

exec docker compose -p kg-alpha \
  -f docker-compose.yml \
  -f docker-compose.alpha.yml \
  --env-file .env.alpha down "$@"
