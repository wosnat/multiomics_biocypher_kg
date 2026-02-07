"""Extract locus_tag IDs from Thompson 2011 supplementary tables.

The Gene column contains composite IDs like 'PMED4_00651 (PMM0063)' or
'P9313_01131 (PMT0107)'. This script extracts the parenthesized locus tag
(e.g., PMM0063, PMT0107) into a new 'locus_tag' column.
"""

import pandas as pd
from pathlib import Path

dir_path = Path(__file__).parent

for csv_file in ["supp table 1.csv", "supp table 2.csv"]:
    df = pd.read_csv(dir_path / csv_file)
    df["primary_id"] = df["Gene"].str.extract(r"^([^(]+)")[0].str.strip()
    df["locus_tag"] = df["Gene"].str.extract(r"\((\w+)\)")
    out_file = csv_file.replace(".csv", "_with_locus_tag.csv")
    df.to_csv(dir_path / out_file, index=False)
    matched = df["locus_tag"].notna().sum()
    total = len(df)
    print(f"{csv_file}: {matched}/{total} IDs extracted -> {out_file}")
