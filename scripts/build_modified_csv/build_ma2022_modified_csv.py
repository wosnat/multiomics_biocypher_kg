"""Ma 2022 Table S4 (iTRAQ proteomics): sign-marked ratios -> log2_fold_change.

The "Mean Ratio" column writes a down-regulated protein's ratio with a minus sign
(-0.55 means coculture/axenic = 0.55) and an up-regulated one as the plain ratio
(1.31). The methods use "fold change > 1.2" as the cutoff, and 1/1.2 = 0.83 is
exactly the edge of the down values (-0.55 ... -0.83), so |x| is the ratio:

    log2_fold_change = log2(|Mean Ratio|)

Every negative value must have |x| < 1 and every positive one x > 1; anything
else means the reading above is wrong, so the script fails rather than guess.

Writes `Supplementary Table S4 Selected differentially expressed proteins_modified.csv`
next to the original (same 2 title rows, so `skip_rows: 2` still applies).

Run: uv run python scripts/build_modified_csv/build_ma2022_modified_csv.py
"""

import math
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Ma 2022")
SRC = PAPER_DIR / "Supplementary Table S4 Selected differentially expressed proteins.csv"
DST = SRC.with_name(SRC.stem + "_modified.csv")
SKIP_ROWS = 2


def main() -> None:
    title = SRC.read_text(encoding="utf-8-sig").splitlines()[:SKIP_ROWS]
    df = pd.read_csv(SRC, skiprows=SKIP_ROWS, encoding="utf-8-sig")
    ratio = pd.to_numeric(df["Mean Ratio"], errors="raise")

    bad = df[((ratio < 0) & (ratio.abs() >= 1)) | ((ratio > 0) & (ratio <= 1)) | (ratio == 0)]
    if not bad.empty:
        raise SystemExit(f"Mean Ratio values outside the sign-marked-ratio reading:\n{bad}")

    df["log2_fold_change"] = [round(math.log2(abs(x)), 6) for x in ratio]

    with DST.open("w", encoding="utf-8", newline="") as f:
        for line in title:
            f.write(line + "\n")
        df.to_csv(f, index=False)
    n_down = int((ratio < 0).sum())
    print(f"wrote {DST} ({len(df)} proteins: {len(df) - n_down} up, {n_down} down; "
          f"log2 range {df['log2_fold_change'].min():.3f} .. {df['log2_fold_change'].max():.3f})")


if __name__ == "__main__":
    main()
