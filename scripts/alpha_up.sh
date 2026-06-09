#!/usr/bin/env bash
# Operator wrapper: bring up / inspect the alpha stack on the lab box.
#
# Usage:
#   ./scripts/alpha_up.sh up -d           # bring alpha-deploy up on the active color
#   ./scripts/alpha_up.sh logs -f deploy  # tail the live alpha-deploy logs
#   ./scripts/alpha_up.sh ps              # list alpha containers
#
# DO NOT use this to perform a RELEASE — use `/release-kg <version> --target
# local`. This wrapper is for ad-hoc operator inspection / restart after the
# release flow has already established an active color.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ ! -f .env.alpha ]]; then
  echo "FATAL: .env.alpha is missing." >&2
  echo "Copy .env.alpha.example to .env.alpha and fill in the placeholders." >&2
  exit 1
fi

# Resolve active color from the state marker — required so we don't silently
# come back on a stale color. /release-kg --target local writes this file at
# the end of every successful release.
COLOR_MARKER=".alpha_active_color"
if [[ ! -f $COLOR_MARKER ]]; then
  echo "FATAL: $COLOR_MARKER missing — no release has been cut against this box yet." >&2
  echo "Run \`/release-kg <version> --target local\` first." >&2
  exit 1
fi
ACTIVE_COLOR=$(cat "$COLOR_MARKER")

# Source .env.alpha so ALPHA_BIND_IP is in this shell — we need it to build
# the host:port strings for KG_DEPLOY_*_BIND. (--env-file feeds compose, not
# this wrapper.)
set -a; . ./.env.alpha; set +a

export ALPHA_DATA_VOLUME="kg-alpha-${ACTIVE_COLOR}"
export KG_DEPLOY_HTTP_BIND="${ALPHA_BIND_IP}:17474"
export KG_DEPLOY_BOLT_BIND="${ALPHA_BIND_IP}:17687"

exec docker compose -p kg-alpha \
  -f docker-compose.yml \
  -f docker-compose.alpha.yml \
  --env-file .env.alpha "$@"
