#!/usr/bin/env python3
"""Thin wrapper around .claude/skills/paperconfig/validate_paperconfig.py that
injects the 4 strains added in the 2026-04-14 CSV-ready-papers batch into the
CANONICAL_GENOMIC_ORGANISMS set before delegating to the real validator.

Use this while the upstream validator's hardcoded organism list is stale.
Remove once the upstream validator is updated (tracked in the follow-up task).

Usage:
    uv run python scripts/validate_paperconfig_with_new_strains.py <paperconfig.yaml>
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VALIDATOR_PATH = REPO / ".claude" / "skills" / "paperconfig" / "validate_paperconfig.py"

NEW_STRAINS = {
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
    "Synechococcus sp. BL107",
    "Marinobacter adhaerens DSM 23420 / HP15",
    "Alteromonas mediterranea DE",
}


def main() -> int:
    spec = importlib.util.spec_from_file_location("validate_paperconfig", VALIDATOR_PATH)
    module = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    assert spec and spec.loader, "Failed to load validator spec"
    spec.loader.exec_module(module)  # type: ignore[union-attr]

    # Inject the new strains into the canonical organism set
    module.CANONICAL_GENOMIC_ORGANISMS |= NEW_STRAINS

    if len(sys.argv) == 2 and sys.argv[1] == "--all":
        all_configs = module.load_all_paperconfigs()
        if not all_configs:
            print("No paperconfigs found in paperconfig_files.txt")
            return 1
        total = len(all_configs)
        passed = 0
        failed = 0
        for path, _ in all_configs:
            print(f"\n{'='*60}")
            print(f"Validating: {path}\n")
            ok = module.validate(str(path))
            passed += int(bool(ok))
            failed += int(not ok)
        print(f"\n{'='*60}")
        print(f"Results: {passed}/{total} passed, {failed}/{total} failed")
        return 0 if failed == 0 else 1

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_paperconfig.yaml>")
        return 2

    config_path = sys.argv[1]
    print(f"Validating: {config_path}\n")
    ok = module.validate(config_path)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main() or 0)
