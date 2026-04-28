"""Step 0 sub-step 7 — Build resolver + hierarchy caches.

Phase 1.1A skeleton — main() raises NotImplementedError. Real implementation
lands in Phase 1.1B once the actual MNX/TCDB/CAZy file shapes are confirmed.

Outputs (when implemented):
- cache/data/mnx/metabolite_resolver.db        (SQLite)
- cache/data/tcdb/tcdb_hierarchy.json
- cache/data/cazy/cazy_hierarchy.json
- cache/data/mnx/metabolite_id_mapping_report.json
"""
from __future__ import annotations

import argparse


def main(force: bool = False) -> None:
    raise NotImplementedError("Phase 1.1B")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild caches even if they exist.")
    args = parser.parse_args()
    main(force=args.force)
