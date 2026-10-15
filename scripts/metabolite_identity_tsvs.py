"""Write the metabolite-identity release TSVs from step 6's identity report.

Reads  cache/data/metabolomics/metabolite_identity_report.json
Writes docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv
         old_id, new_id, basis, role_id       (report["remap"])
       docs/kg-changes/metabolite-identity/metabolite_assay_repoints.tsv
         paper, table, paper_name, old_id, new_id, stereo_assignment, reason   (report["repoints"])

Both are sorted, tab-separated and newline-terminated. The report's remap / repoints are
diagnostics against its baseline (the previous kegg_data.json), so run this on the report of
the step-6 run that made the change being released, not on a later no-change rerun.

    uv run python scripts/metabolite_identity_tsvs.py [--report PATH] [--out-dir PATH]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "cache" / "data" / "metabolomics" / "metabolite_identity_report.json"
OUT_DIR = ROOT / "docs" / "kg-changes" / "metabolite-identity"

REMAP_COLS = ("old_id", "new_id", "basis", "role_id")
REPOINT_COLS = ("paper", "table", "paper_name", "old_id", "new_id", "stereo_assignment", "reason")


def _cell(v) -> str:
    s = "" if v is None else str(v)
    if "\t" in s or "\n" in s:
        raise ValueError(f"tab/newline inside a TSV cell: {s!r}")
    return s


def render_tsv(rows: list[dict], cols: tuple[str, ...]) -> str:
    lines = sorted("\t".join(_cell(r.get(c)) for c in cols) for r in rows)
    return "\n".join(["\t".join(cols), *lines]) + "\n"


def write_tsvs(report: dict, out_dir: Path) -> dict[Path, str]:
    files = {out_dir / "metabolite_id_remap.tsv": render_tsv(report.get("remap", []), REMAP_COLS),
             out_dir / "metabolite_assay_repoints.tsv": render_tsv(report.get("repoints", []), REPOINT_COLS)}
    out_dir.mkdir(parents=True, exist_ok=True)
    for path, text in files.items():
        path.write_text(text)
    return files


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", type=Path, default=REPORT)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    report = json.loads(args.report.read_text())
    for path, text in write_tsvs(report, args.out_dir).items():
        print(f"wrote {path} ({text.count(chr(10)) - 1} rows)")


if __name__ == "__main__":
    main()
