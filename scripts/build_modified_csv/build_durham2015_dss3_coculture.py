"""Build the R. pomeroyi DSS-3 coculture DE table for Durham 2015 (PNAS).

Source: ``pnas.1413137112.st01.docx`` (Table S1): DSS-3 genes differentially
expressed after 8 h of co-culture with Thalassiosira pseudonana CCMP1335
(n = 4) vs DSS-3 alone in f/2 (+Si, -B12) (n = 2). baySeq, three independent
subsampled analyses (Likelihood.1-3 / FDR.1-3); a gene is listed when
FDR < 0.05 in at least two of the three analyses -> significant-only table.

Fold change is LINEAR and signed: positive = X-fold up in co-culture,
negative = X-fold down (values run -52 .. +373.5, none in (-1, 1); the paper
text says "up-regulated 56- to 374-fold"). Converted to log2:
``log2fc = sign(FC) * log2(|FC|)``.

``fdr_median`` = median of FDR.1-3. Because a gene is listed only when >= 2 of
3 FDRs are < 0.05, the median is < 0.05 for every listed gene, so it doubles
as the KG significance column.

The docx is parsed with the standard library (zipfile + ElementTree) so no
python-docx dependency is needed.

Skipped: ST02 (T. pseudonana genes; diatom not a KG genome), ST03 (ortholog
context table), ST04.

Idempotent; never edits the original docx.
"""
from __future__ import annotations

import math
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import median

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/Durham 2015"
SRC = PAPER / "pnas.1413137112.st01.docx"
OUT = PAPER / "durham2015_st01_dss3_coculture.csv"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def read_table(path: Path) -> list[list[str]]:
    root = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    tbl = next(root.iter(W + "tbl"))
    return [
        ["".join(t.text or "" for t in c.iter(W + "t")).strip() for c in r.findall(W + "tc")]
        for r in tbl.findall(W + "tr")
    ]


def main() -> None:
    rows = read_table(SRC)
    hdr_i = next(i for i, r in enumerate(rows) if r and r[0] == "Locus Tag")
    header = rows[hdr_i]
    recs = []
    for r in rows[hdr_i + 1:]:
        if len(r) != len(header) or not r[0]:
            continue
        d = dict(zip(header, r))
        fc = float(d["Fold Change"])
        if abs(fc) < 1:
            raise ValueError(f"unexpected |FC| < 1 for {d['Locus Tag']}: {fc}")
        fdrs = [float(d[f"FDR.{k}"]) for k in (1, 2, 3)]
        recs.append({
            "locus_tag": d["Locus Tag"],
            "product": d["Description"].replace("|", ",").replace("'", "^"),
            "fold_change_linear": fc,
            "log2fc_coculture_vs_alone": math.copysign(math.log2(abs(fc)), fc),
            "fdr_1": fdrs[0], "fdr_2": fdrs[1], "fdr_3": fdrs[2],
            "fdr_median": median(fdrs),
            "n_fdr_below_0.05": sum(f < 0.05 for f in fdrs),
        })
    out = pd.DataFrame(recs)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(out)} rows; "
          f"max median FDR {out['fdr_median'].max():.4f}; "
          f"rows with <2 FDR<0.05: {(out['n_fdr_below_0.05'] < 2).sum()}")


if __name__ == "__main__":
    main()
