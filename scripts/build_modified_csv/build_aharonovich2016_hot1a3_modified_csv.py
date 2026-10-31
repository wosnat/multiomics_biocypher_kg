"""Aharonovich 2016 HOT1A3 tables: re-derive log2 fold change with a +1 pseudocount.

The two Alteromonas HOT1A3 tables publish a linear "fold change" computed as
(a + 0.01) / (b + 0.01) on the Rockhopper expression values (the +.01 columns are
in the CSV). When one arm is zero this gives 100 x count + 1: 1 read vs 0 reads is
101-fold (log2 6.7), 20 vs 0 is 2001-fold. Those zero-count genes take every
one of the top-100 |FC| rows, so rank_by_effect would rank sequencing depth.

This script writes `<table>_modified.csv` with

    log2_fold_change = log2((a + 1) / (b + 1))     (treatment / control)

which bounds a 1-vs-0 gene at log2 1. Rows with both arms zero (no measurement)
and rows the paper gives no fold change for get a blank log2_fold_change, so the
adapter builds no edge for them. q-values are kept as published.

Before writing, it checks that the published ratio equals (a + 0.01) / (b + 0.01)
for every row that has one, so a misread column order fails loudly.

Run: uv run python scripts/build_modified_csv/build_aharonovich2016_hot1a3_modified_csv.py
"""

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Prochlorococcus/papers_and_supp/Aharonovich 2016")
# (file, treatment-arm label, control-arm label); positional columns in both files:
# 0 Synonym, 1 Product, 2 treatment expression, 3 control expression,
# 4/5 the same + 0.01, 6 qValue, 7 published fold change (treatment / control)
TABLES = [
    ("HOT1A3 coculture with 9313 ver coculture with 9313 diluted.csv", "9313+1A3", "9313+1A3 dil"),
    ("HOT1A3 coculture with MED4 vs coculture with MIT9313.csv", "9313+1A3", "med4+1A3"),
]
PSEUDOCOUNT = 1.0


def build(name: str, treat: str, ctrl: str) -> None:
    src = PAPER_DIR / name
    dst = src.with_name(src.stem + "_modified.csv")
    raw = pd.read_csv(src, header=None, skiprows=1, encoding="utf-8-sig")
    header = pd.read_csv(src, nrows=0, encoding="utf-8-sig").columns.tolist()

    a = pd.to_numeric(raw[2], errors="coerce")
    b = pd.to_numeric(raw[3], errors="coerce")
    fc = pd.to_numeric(raw[7], errors="coerce")

    has_fc = fc.notna()
    expected = (a + 0.01) / (b + 0.01)
    rel = ((fc - expected).abs() / expected)[has_fc]
    if (rel > 1e-3).any():
        raise SystemExit(f"{name}: published fold change != (a+0.01)/(b+0.01) on "
                         f"{int((rel > 1e-3).sum())} rows — column order misread?")

    both_zero = (a == 0) & (b == 0)
    keep = has_fc & ~both_zero & a.notna() & b.notna()
    log2fc = [round(math.log2((x + PSEUDOCOUNT) / (y + PSEUDOCOUNT)), 6) if k else None
              for x, y, k in zip(a, b, keep)]

    out = pd.DataFrame({
        "Synonym": raw[0],
        "Product": raw[1],
        f"Expression {treat}": a,
        f"Expression {ctrl}": b,
        header[6].strip(): raw[6],
        f"published fold change (+0.01 pseudocount) {header[7].strip()}": fc,
        "log2_fold_change": log2fc,
    })
    out.to_csv(dst, index=False)
    lf = pd.Series([v for v in log2fc if v is not None])
    print(f"wrote {dst.name}: {int(keep.sum())} rows with log2_fold_change "
          f"({int((has_fc & both_zero).sum())} both-zero rows blanked); "
          f"log2 range {lf.min():.2f} .. {lf.max():.2f} "
          f"(published {math.log2(fc[keep].min()):.2f} .. {math.log2(fc[keep].max()):.2f})")


def main() -> None:
    for name, treat, ctrl in TABLES:
        build(name, treat, ctrl)


if __name__ == "__main__":
    main()
