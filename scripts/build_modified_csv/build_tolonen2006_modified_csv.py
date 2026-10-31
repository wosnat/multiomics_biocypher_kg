"""Tolonen 2006 MED4 expression table: drop the garbled PMM1118 (hli4) row.

In `medExpression.1.csv` (and the original xls) the PMM1118 row holds 15 values under
8 value columns: two probesets (MED4_ARR_1112_at and _s_at) are interleaved, so every
column after t=0h carries the wrong probeset's number (t=3h = the _s_at t=0h value,
t=6h = the _at t=3h value, ...). Rebuilding the row needs the Goldenspike per-timepoint
files, which do not cover t=12h or cyanate, so the row is dropped.

Writes `medExpression.1_modified.csv` next to the original; everything else is copied
byte for byte.

Run: uv run python scripts/build_modified_csv/build_tolonen2006_modified_csv.py
"""

from pathlib import Path

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/tolonen 2006")
SRC = PAPER_DIR / "medExpression.1.csv"
DST = SRC.with_name(SRC.stem + "_modified.csv")
DROP = {"PMM1118"}


def main() -> None:
    lines = SRC.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    kept = [ln for ln in lines if ln.split(",", 1)[0] not in DROP]
    if len(lines) - len(kept) != len(DROP):
        raise SystemExit(f"expected to drop {len(DROP)} rows, dropped {len(lines) - len(kept)}")
    DST.write_text("".join(kept), encoding="utf-8")
    print(f"wrote {DST} ({len(kept) - 1} genes; dropped {sorted(DROP)})")


if __name__ == "__main__":
    main()
