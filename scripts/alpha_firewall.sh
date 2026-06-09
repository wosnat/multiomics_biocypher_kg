#!/usr/bin/env bash
# Operator helper: apply / inspect / remove the DOCKER-USER firewall allowlist
# that restricts the alpha Neo4j ports to a confirmed lab subnet.
#
# WHY: the alpha deploy binds to the campus-routable IP 132.75.249.47, so the
# published ports 17474/17687 are reachable by anything on campus that can route
# there. Docker publishes ports by inserting rules into the DOCKER-USER iptables
# chain, which is consulted BEFORE ufw's INPUT chain — so a plain `ufw allow/deny`
# does NOT filter Docker-published ports. We filter in DOCKER-USER instead.
# See plans/alpha_release.md §2.6.
#
# Usage (needs root — run with sudo):
#   sudo ./scripts/alpha_firewall.sh status                 # show current DOCKER-USER rules for the alpha ports
#   sudo ./scripts/alpha_firewall.sh apply  <LAB_SUBNET_CIDR> [--persist]
#   sudo ./scripts/alpha_firewall.sh remove <LAB_SUBNET_CIDR>
#
# Examples:
#   sudo ./scripts/alpha_firewall.sh apply 132.75.249.32/27 --persist
#   sudo ./scripts/alpha_firewall.sh status
#
# The rule installed per port is:
#   iptables -I DOCKER-USER -p tcp --dport <PORT> ! -s <CIDR> -j DROP
# i.e. "drop traffic to <PORT> that is NOT from <CIDR>" — an allowlist for the
# lab subnet. `apply` is idempotent (checks with -C before inserting with -I).
#
# Confirm <LAB_SUBNET_CIDR> with IT before trusting it. The box's own on-link
# subnet is 132.75.249.32/27 (interface enp6s0f1 = 132.75.249.47/27); use the
# subnet your TESTERS sit on. Do NOT allowlist the whole campus /16.

set -euo pipefail

# Alpha published ports (keep in sync with ALPHA_PUBLISHED_{HTTP,BOLT}_PORT in
# .claude/skills/release-kg/release_kg.py).
ALPHA_PORTS=(17687 17474)
CHAIN="DOCKER-USER"

die() { echo "FATAL: $*" >&2; exit 1; }

require_root() {
  [[ "${EUID:-$(id -u)}" -eq 0 ]] || die "must run as root (use sudo)."
}

validate_cidr() {
  local cidr="$1"
  [[ "$cidr" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}/[0-9]{1,2}$ ]] \
    || die "CIDR '$cidr' is malformed (expected A.B.C.D/NN, e.g. 132.75.249.32/27)."
  local prefix="${cidr##*/}"
  # Refuse dangerously broad allowlists — a /16 here would expose the alpha ports
  # to the entire campus, defeating the point.
  if (( prefix < 24 )); then
    die "refusing CIDR '$cidr': prefix /$prefix is too broad (< /24). Do not allowlist the campus /16; use the tester subnet (the box's own is 132.75.249.32/27)."
  fi
}

rule_args() {  # $1 = port
  echo "-p tcp --dport $1 ! -s $CIDR -j DROP"
}

cmd_status() {
  require_root
  echo "=== current $CHAIN chain ==="
  iptables -L "$CHAIN" -n --line-numbers
  echo
  for port in "${ALPHA_PORTS[@]}"; do
    if iptables -L "$CHAIN" -n | grep -q "dpt:$port"; then
      echo "port $port: rule PRESENT ✓"
    else
      echo "port $port: NO rule ✗ (port is unrestricted)"
    fi
  done
}

cmd_apply() {
  require_root
  [[ -n "${CIDR:-}" ]] || die "apply requires a CIDR. Usage: $0 apply <LAB_SUBNET_CIDR> [--persist]"
  validate_cidr "$CIDR"
  for port in "${ALPHA_PORTS[@]}"; do
    # shellcheck disable=SC2046
    if iptables -C "$CHAIN" $(rule_args "$port") 2>/dev/null; then
      echo "port $port: allowlist rule already present for $CIDR — skipping"
    else
      # shellcheck disable=SC2046
      iptables -I "$CHAIN" $(rule_args "$port")
      echo "port $port: inserted DROP for sources outside $CIDR ✓"
    fi
  done
  if [[ "${PERSIST:-0}" == "1" ]]; then
    echo "=== persisting via netfilter-persistent ==="
    if ! command -v netfilter-persistent >/dev/null 2>&1; then
      apt-get install -y iptables-persistent
    fi
    netfilter-persistent save
    echo "persisted ✓"
  else
    echo
    echo "NOTE: rules are NOT persisted across reboot. Re-run with --persist, or:"
    echo "      apt-get install -y iptables-persistent && netfilter-persistent save"
  fi
  echo
  cmd_status
}

cmd_remove() {
  require_root
  [[ -n "${CIDR:-}" ]] || die "remove requires the same CIDR used in apply. Usage: $0 remove <LAB_SUBNET_CIDR>"
  validate_cidr "$CIDR"
  for port in "${ALPHA_PORTS[@]}"; do
    # shellcheck disable=SC2046
    while iptables -C "$CHAIN" $(rule_args "$port") 2>/dev/null; do
      # shellcheck disable=SC2046
      iptables -D "$CHAIN" $(rule_args "$port")
      echo "port $port: removed DROP rule for $CIDR"
    done
  done
  echo "done. (If you persisted earlier, re-run netfilter-persistent save.)"
}

ACTION="${1:-}"
CIDR="${2:-}"
PERSIST=0
[[ "${3:-}" == "--persist" ]] && PERSIST=1

case "$ACTION" in
  status) cmd_status ;;
  apply)  cmd_apply  ;;
  remove) cmd_remove ;;
  *) cat >&2 <<EOF
Usage:
  sudo $0 status
  sudo $0 apply  <LAB_SUBNET_CIDR> [--persist]
  sudo $0 remove <LAB_SUBNET_CIDR>

Confirm <LAB_SUBNET_CIDR> with IT. The box's own /27 is 132.75.249.32/27.
EOF
     exit 2 ;;
esac
