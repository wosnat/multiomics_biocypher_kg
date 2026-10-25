"""Build Ford 2021 per-strain SBP class tables + the paperconfig.

Ford et al. 2021 (Biochem Soc Trans 49:2465, doi 10.1042/BST20200244) is a
review whose Supplementary Data S1 (`bst-2020-0244c_supp.xlsx`) lists the 1,257
ABC-transporter substrate-binding proteins (SBPs) found by text search of
Cyanorak v2.1 across 97 picocyanobacterial genomes, one row per gene, with the
gene's Cyanorak cluster (CLOG).

The paper's functional calls are per cluster; each gene inherits its
cluster's call:
  * Figure 2 groups every cluster into a major nutrient class
    (Carbon / Nitrogen / Phosphorus / Sulphur / Trace metals / Unknown).
  * Figure 1 labels most clusters with a predicted substrate; Table 1 and the
    text add a few more.
Both are transcribed in NUTRIENT_CLASS / SUBSTRATE below with their source.
Nothing is inferred beyond what the paper states: a cluster with no stated
substrate is `not_stated`.

They enter the KG as two categorical DerivedMetrics per strain
(`sbp_nutrient_class`, `sbp_predicted_substrate`), each strain's pair bound to
one annotation Experiment (`omics_type: BIOINFORMATIC`,
`treatment_type: [genomic_analysis]`, `background_factors: []` — computed from
genome annotation, no sample, no culture context, no DE).

Outputs (idempotent):
  sbp_classes/<strain>_sbp_classes.csv      one row per SBP gene of that strain
  paperconfig.yaml                          one experiment + derived_metrics_table per strain

Strains are matched to the KG by Cyanorak organism name
(`cyanobacteria_genomes.csv` `cyanorak_organism`). Strains with no registered
genome are skipped and reported (UH18301: no public assembly).

Run:  uv run python scripts/build_modified_csv/build_ford2021_sbp_classes.py
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "data/Synechococcus/papers_and_supp/Ford 2021"
PAPER_REL = "data/Synechococcus/papers_and_supp/Ford 2021"
SUPP = PAPER_DIR / "bst-2020-0244c_supp.xlsx"
GENOMES_CSV = ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"

# ── Figure 2: cluster -> major nutrient class ────────────────────────────────
# Column groups of Figure 2. The "Hypothetical" column is the caption's list
# (CK_264, CK_2249, CK_3177, CK_4819, CK_4821, CK_8828, CK_23352, CK_31409,
# CK_32631, CK_33184, CK_34148, CK_37737, CK_44411, CK_44797, CK_46634,
# CK_46636); CK_31409 also has its own column under Unknown.
# CK_48211 (1 gene, Cyanobium NS01) is in S1 but named nowhere in Figure 2;
# every S1 cluster without a named nutrient column is in the hypothetical
# group, so it is "unknown" (the caption's CK_4821 matches no S1 cluster).
# CK_2249, CK_3177, CK_8828 and CK_37737 carry a Figure 1 substrate label yet
# sit in the Figure 2 hypothetical group; both calls are kept as the paper
# states them (nutrient class unknown, substrate from Figure 1).
_FIG2 = {
    "carbon": [1455, 1342],
    "nitrogen": [7046, 1489, 1014, 1944, 3072, 2165, 76],
    "phosphorus": [23, 43821, 860, 6203, 56876, 1829],
    "sulfur": [37708, 57276],
    "trace_metals": [40935, 919, 2462, 1600, 53501, 68, 57079, 51448],
    "unknown": [31409, 264, 2249, 3177, 4819, 8828, 23352, 32631, 33184,
                34148, 37737, 44411, 44797, 46634, 46636, 48211],
}
NUTRIENT_CLASS = {ck: cls for cls, cks in _FIG2.items() for ck in cks}

# ── Figure 1 / Table 1 / text: cluster -> predicted substrate ────────────────
# (value, source). Values are KG-minted snake_case (vocabulary rule R1).
SUBSTRATE = {
    23: ("phosphate", "Fig. 1"),
    1829: ("phosphate", "Fig. 1"),
    43821: ("phosphate", "Fig. 1 (PstS)"),
    860: ("phosphonate", "Fig. 1 (PhnD1; text: binds phosphite)"),
    6203: ("phosphonate", "Fig. 1 (PhnD2)"),
    56876: ("phosphonate", "text (PhnD3)"),
    8828: ("phosphonate", "Fig. 1"),
    37737: ("hypothetical_phosphonate", "Fig. 1"),
    2165: ("cyanate", "Fig. 1"),
    2249: ("nitrate", "Fig. 1"),
    3072: ("nitrate_nitrite", "Table 1 (NrtA)"),
    1489: ("amino_acids", "Fig. 1 (text: acidic or polar amino acids)"),
    1014: ("peptide", "Fig. 1"),
    1944: ("glycine_betaine", "Fig. 1"),
    76: ("urea", "Fig. 1 (UrtA)"),
    1342: ("sugar", "Fig. 1 (Table 1: predicted chitobiose, ChiE)"),
    1455: ("glucose_glycerol_trehalose", "Fig. 1"),
    57079: ("iron", "Fig. 1 (IdiA1)"),
    68: ("iron", "Fig. 1 (IdiA2)"),
    51448: ("iron_siderophore", "Fig. 1 (text: Fe siderophore/vitamin B12)"),
    919: ("manganese", "Fig. 1"),
    3177: ("manganese", "Fig. 1"),
    2462: ("zinc", "Fig. 1"),
    1600: ("manganese_zinc", "Fig. 1"),
    40935: ("molybdenum", "text"),
    53501: ("cobalt_nickel", "text"),
}
NOT_STATED = "not_stated"
CLASS_VALUES = list(_FIG2)
SUBSTRATE_VALUES = sorted({v for v, _ in SUBSTRATE.values()}) + [NOT_STATED]


def _ck_num(cluster: str) -> int:
    return int(re.sub(r"^CK_0*", "", cluster))


def _slug(strain: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", strain.lower()).strip("_")


def _clean(s) -> str:
    """Supp descriptions use '_' for spaces; drop the KG-unsafe characters."""
    if not isinstance(s, str):
        return ""
    return s.replace("_", " ").replace("'", "^").replace("|", "").strip()


def load_genomes() -> dict[str, dict]:
    with open(GENOMES_CSV) as fh:
        rows = csv.DictReader(l for l in fh if l.strip() and not l.startswith("#"))
        return {r["cyanorak_organism"]: r for r in rows if r["cyanorak_organism"]}


def main() -> None:
    df = pd.read_excel(SUPP)
    df["ck"] = df["Cluster"].map(_ck_num)
    missing = sorted(set(df["ck"]) - set(NUTRIENT_CLASS))
    assert not missing, f"clusters without a Figure 2 class: {missing}"
    df["cyanorak_organism"] = df["Cyanorak ID"].str.extract(r"^CK_(.+)_\d+$")[0]

    genomes = load_genomes()
    out_dir = PAPER_DIR / "sbp_classes"
    out_dir.mkdir(exist_ok=True)

    entries: list[tuple[str, dict]] = []
    skipped = Counter()
    for org, sub in sorted(df.groupby("cyanorak_organism"), key=lambda kv: kv[0]):
        g = genomes.get(org)
        if g is None:
            skipped[org] = len(sub)
            continue
        strain = g["strain_name"]
        csv_name = f"{strain}_sbp_classes.csv"
        rows = pd.DataFrame({
            "Cyanorak ID": sub["Cyanorak ID"],
            "cyanorak_cluster": sub["Cluster"],
            "gene_name": sub["Gene"].fillna(""),
            "product": sub["Description"].map(_clean),
            "nutrient_class": sub["ck"].map(NUTRIENT_CLASS),
            "predicted_substrate": sub["ck"].map(lambda c: SUBSTRATE.get(c, (NOT_STATED,))[0]),
        })
        rows.to_csv(out_dir / csv_name, index=False)
        entries.append((_slug(strain), {
            "organism": g["preferred_name"], "csv": f"{PAPER_REL}/sbp_classes/{csv_name}",
            "n": len(rows),
        }))

    write_paperconfig(entries)
    print(f"{len(entries)} strains, {sum(e['n'] for _, e in entries)} SBP genes written")
    if skipped:
        print("skipped (no registered genome):", dict(skipped))


def _yaml_list(values: list[str], indent: str) -> list[str]:
    return [f'{indent}- "{v}"' for v in values]


def write_paperconfig(entries: list[tuple[str, dict]]) -> None:
    class_desc = (
        "Major nutrient class of the ABC substrate-binding protein (SBP), from Ford 2021 Figure 2: "
        "carbon, nitrogen, phosphorus, sulfur, trace_metals, or unknown (hypothetical/unknown-function "
        "clusters). Called per Cyanorak cluster and inherited by each member gene. Prediction from "
        "genome annotation, not experimentally validated."
    )
    subst_desc = (
        "Substrate the paper predicts for the Cyanorak cluster of the SBP (Ford 2021 Figure 1, Table 1 and "
        "text), inherited by each member gene; not_stated when the paper names none. Some clusters "
        "with a stated substrate sit in the unknown nutrient class (Figure 2 hypothetical group). "
        "Prediction from genome annotation, not experimentally validated."
    )
    lines = [
        "# Generated by scripts/build_modified_csv/build_ford2021_sbp_classes.py — edit the script, not this file.",
        "# Ford 2021 is a review: no experiments, no DE. Its Supplementary Data S1 lists the 1,257",
        "# ABC substrate-binding proteins (SBPs) Cyanorak v2.1 annotates across 97 picocyanobacterial",
        "# genomes. Each strain gets one annotation experiment (omics_type BIOINFORMATIC, treatment_type",
        "# [genomic_analysis], background_factors [] — no sample, no culture, no DE) and a",
        "# derived_metrics_table with two categorical metrics: the paper's Figure 2 nutrient class and",
        "# its predicted substrate, both called per Cyanorak cluster.",
        "# 01whole.pdf is the first author's 2017 MRes thesis on Synechococcus SBPs (background only).",
        "publication:",
        '  papername: "Ford 2021"',
        '  doi: "10.1042/BST20200244"',
        f'  papermainpdf: "{PAPER_REL}/bst-2020-0244c.pdf"',
        "",
        "  experiments:",
    ]
    for key, e in entries:
        lines += [
            f"    {key}_sbp_annotation:",
            f'      name: "{e["organism"]} ABC substrate-binding protein annotation survey (Ford 2021)"',
            f'      organism: "{e["organism"]}"',
            "      omics_type: BIOINFORMATIC",
            "      treatment_type: [genomic_analysis]",
            "      background_factors: []   # computed from genome annotation: no sample, no culture",
            '      treatment_condition: "Text search of Cyanorak v2.1 annotations for ABC-transporter substrate-binding proteins (no experimental treatment)"',
            f'      experimental_context: "{e["n"]} SBP genes in this genome (Ford 2021 Supplementary Data S1); nutrient class and predicted substrate are the per-Cyanorak-cluster calls of the paper (Figures 1-2, Table 1, text), not experimentally validated."',
        ]
    lines += ["", "  supplementary_materials:"]
    for key, e in entries:
        lines += [
            f"    {key}_sbp:",
            "      type: derived_metrics_table",
            f'      filename: "{e["csv"]}"',
            f'      organism: "{e["organism"]}"',
            f"      experiment: {key}_sbp_annotation",
            '      name_col: "Cyanorak ID"',
            "      id_columns:",
            '        - column: "Cyanorak ID"',
            "          id_type: locus_tag_cyanorak",
            "      product_columns:",
            '        - column: "product"',
            "      metrics:",
            "        - metric_type: sbp_nutrient_class",
            '          name: "Predicted SBP nutrient class (Ford 2021)"',
            "          value_kind: categorical",
            "          value_col: nutrient_class",
            "          allowed_categories:",
            *_yaml_list(CLASS_VALUES, "            "),
            f'          field_description: "{class_desc}"',
            "        - metric_type: sbp_predicted_substrate",
            '          name: "Predicted SBP substrate (Ford 2021)"',
            "          value_kind: categorical",
            "          value_col: predicted_substrate",
            "          allowed_categories:",
            *_yaml_list(SUBSTRATE_VALUES, "            "),
            f'          field_description: "{subst_desc}"',
            "",
        ]
    (PAPER_DIR / "paperconfig.yaml").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
