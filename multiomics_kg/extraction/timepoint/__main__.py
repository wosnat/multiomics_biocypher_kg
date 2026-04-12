"""Dispatch `python -m multiomics_kg.extraction.timepoint <subcommand>`.

Subcommands: extract, merge, remap, report.
"""
import sys


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python -m multiomics_kg.extraction.timepoint "
              "{extract|merge|remap|report} [args...]", file=sys.stderr)
        return 2

    sub = sys.argv[1]
    argv = sys.argv[2:]

    if sub == "extract":
        from multiomics_kg.extraction.timepoint.extract import _main as extract_main
        return extract_main(argv)
    if sub == "merge":
        from multiomics_kg.extraction.timepoint.merge import _main as merge_main
        return merge_main(argv)
    if sub == "remap":
        from multiomics_kg.extraction.timepoint.remap import _main as remap_main
        return remap_main(argv)
    if sub == "report":
        from multiomics_kg.extraction.timepoint.report import _main as report_main
        return report_main(argv)

    print(f"Unknown subcommand: {sub}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
