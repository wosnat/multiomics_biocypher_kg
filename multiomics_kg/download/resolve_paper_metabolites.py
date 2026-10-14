"""Retired (2026-10): step 7 is folded into step 6.

Step 6 (build_kegg_metabolism_xrefs) decides every paper-metabolite resolution and writes
each <stem>_resolved.csv from the same in-memory decisions, so the graph's metabolite set
and the paper tables can't drift (spec 2026-10-13 §5.3). Re-resolve a paper with
`bash scripts/prepare_data.sh --steps 6 --force` (about a minute, no network).
"""
import sys


def main(argv=None) -> int:
    print("step 7 is folded into step 6 (build_kegg_metabolism_xrefs); nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
