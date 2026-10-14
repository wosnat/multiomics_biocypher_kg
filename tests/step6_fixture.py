"""Synthetic step-6 inputs for the pure identity pipeline (stages E-P).

The scenario, end to end:
- Nitrate: KEGG ``C00244`` (catalysis) and the TCDB substrate ``chebi:17632`` are ChEBI
  conjugates, so they merge (``identity_basis: [chebi_conjugate]``) onto the KEGG node.
- Alanine stereo family, linked by ChEBI ``is_a`` / enantiomer with a shared InChIKey first
  block: L-Ala ``C00041`` (2 catalyst genes via R1), D-Ala ``C00133`` (1 gene via R2), the
  unspecified DL ``C01401`` (0 genes, transport only). The paper row "Alanine" matches both
  L and D, one family, so the dominant-form rule picks L (2 >= 2 x 1) and freezes it.
- Tryptophan: L-Trp ``C00078`` (1 gene via R3) and D-Trp ``C00525``, a measurement-only
  candidate (``evidence_sources: []``) that only exists because a previous build resolved
  paper row 1 to it. The unqualified "Tryptophan" now resolves to L-Trp, so stage K drops
  D-Trp and records it in L-Trp's ``recurated_from_ids``.
- Pathways: ko00250 / ko00380 are catalysis pathways (``base_pathways``); ko00999 (own
  subcategory and category) is reachable only through D-Trp, so it is pruned with it.
- One role substrate (``chebi:10004`` → ``chebi.role:33229``) was partitioned out in stage C.
- One ChEBI → KEGG re-key from stage D (``chebi:16977`` → ``C00041``).
- KEGG BRITE br08001 puts ``C00041`` under Peptides / Amino acids; every other compound tree
  holds only ``C99999`` (not ours), so it prunes to empty.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from multiomics_kg.utils.metabolite_identity.canonical import IdentityOverrides
from multiomics_kg.utils.metabolite_identity.classes import COMPOUND_BRITE_TREES
from multiomics_kg.utils.metabolite_identity.paper import PaperRow, PaperTable
from multiomics_kg.utils.metabolite_identity.pipeline import Step6Inputs
from tests.metabolite_identity_helpers import cpd, make_chebi, term

PAPER = "data/Prochlorococcus/papers_and_supp/Fixture 2026/paperconfig.yaml"
TABLE = "fixture_metabolites"
K = "kegg.compound:"


def _chebi():
    return make_chebi(
        CHEBI_48107=term("nitric acid", charge=0, kegg=["C00244"], conj=["CHEBI:17632"],
                         inchikey="GRYLNZFGIOXLOG-UHFFFAOYSA-N"),
        CHEBI_17632=term("nitrate", charge=-1, conj=["CHEBI:48107"],
                         inchikey="NHNBFGGVMKEFGY-UHFFFAOYSA-N"),
        CHEBI_16977=term("L-alanine", charge=0, kegg=["C00041"], is_a=["CHEBI:16449"],
                         enan=["CHEBI:15570"], inchikey="QNAYBMKLOCPYGJ-REOHCLBHSA-N"),
        CHEBI_15570=term("D-alanine", charge=0, kegg=["C00133"], is_a=["CHEBI:16449"],
                         enan=["CHEBI:16977"], inchikey="QNAYBMKLOCPYGJ-UWTATZPHSA-N"),
        CHEBI_16449=term("alanine", charge=0, kegg=["C01401"],
                         inchikey="QNAYBMKLOCPYGJ-UHFFFAOYSA-N"),
        CHEBI_16828=term("L-tryptophan", charge=0, kegg=["C00078"], enan=["CHEBI:16296"],
                         inchikey="QIVBCDIJIAJPQS-VIFPVBQESA-N"),
        CHEBI_16296=term("D-tryptophan", charge=0, kegg=["C00525"], enan=["CHEBI:16828"],
                         inchikey="QIVBCDIJIAJPQS-SECBINFHSA-N"),
    )


def _table():
    return {
        f"{K}C00244": cpd("Nitrate", formula="NO3", inchikey="NHNBFGGVMKEFGY-UHFFFAOYSA-N",
                          chebi_id="48107", mnxm_id="MNXM732398"),
        "chebi:17632": cpd("nitrate", formula="NO3", chebi_id="17632", ev=("transport",)),
        f"{K}C00041": cpd("L-Alanine", formula="C3H7NO2", inchikey="QNAYBMKLOCPYGJ-REOHCLBHSA-N",
                          chebi_id="16977", mnxm_id="MNXM32", pathways=("ko00250",)),
        f"{K}C00133": cpd("D-Alanine", formula="C3H7NO2", inchikey="QNAYBMKLOCPYGJ-UWTATZPHSA-N",
                          chebi_id="15570"),
        f"{K}C01401": cpd("Alanine", formula="C3H7NO2", inchikey="QNAYBMKLOCPYGJ-UHFFFAOYSA-N",
                          chebi_id="16449", ev=("transport",)),
        f"{K}C00078": cpd("L-Tryptophan", formula="C11H12N2O2",
                          inchikey="QIVBCDIJIAJPQS-VIFPVBQESA-N", chebi_id="16828",
                          pathways=("ko00380",)),
        # Measurement-only candidate: enriched at stage D only because a paper row pointed at it.
        f"{K}C00525": cpd("D-Tryptophan", formula="C11H12N2O2",
                          inchikey="QIVBCDIJIAJPQS-SECBINFHSA-N", chebi_id="16296", ev=(),
                          pathways=("ko00999",)),
    }


def _empty_tree(code):
    return {"name": code, "children": [{"name": "X", "children": [{"name": "C99999  none"}]}]}


def make_inputs(tmp_path: Path, extra_rows=(), extra_candidates=None) -> Step6Inputs:
    names = [("Alanine", None), ("Tryptophan", None), *extra_rows]
    csv = tmp_path / "paper" / "fixture_metabolites.csv"
    csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"compound": [n for n, _ in names],
                  "value": [str(i + 1) for i in range(len(names))]}).to_csv(csv, index=False)
    rows = [PaperRow(PAPER, TABLE, i, n, d, "kegg_direct" if d else None)
            for i, (n, d) in enumerate(names)]
    tables = [PaperTable(PAPER, TABLE, csv, 0, rows)]

    name_candidates = {"alanine": {f"{K}C00133", f"{K}C00041"},
                       "tryptophan": {f"{K}C00525", f"{K}C00078"}}
    for k, v in (extra_candidates or {}).items():
        name_candidates[k] = set(name_candidates.get(k, set())) | set(v)

    table = _table()
    reactions = {"R1": {"compounds": ["C00041", "C00244"], "pathways": []},
                 "R2": {"compounds": ["C00133"], "pathways": []},
                 "R3": {"compounds": ["C00078"], "pathways": ["ko00380", "ko00999"]}}
    brite_raw = {code: _empty_tree(code) for code in COMPOUND_BRITE_TREES}
    brite_raw["br08001"] = {"name": "br08001", "children": [
        {"name": "Peptides", "children": [
            {"name": "Amino acids", "children": [{"name": "C00041  L-Alanine"}]}]}]}
    return Step6Inputs(
        table=table,
        kegg_rest={"kos": {"K00001": {"name": "k", "pathways": ["ko00250", "ko00999"]}},
                   "pathways": {"ko00250": {"name": "Alanine metabolism", "subcategory": "09105"},
                                "ko00380": {"name": "Tryptophan metabolism", "subcategory": "09105"},
                                "ko00999": {"name": "Paper-only pathway", "subcategory": "09199"}},
                   "subcategories": {"09105": {"name": "Amino acid metabolism", "category": "09100"},
                                     "09199": {"name": "Paper-only subcategory", "category": "09190"}},
                   "categories": {"09100": {"name": "Metabolism"}, "09190": {"name": "Paper-only"}},
                   "reactions": reactions},
        subtree_substrates={"2.A.1": ["chebi:17632", f"{K}C00244", f"{K}C01401"],
                            "3.A.3": [f"{K}C00244"]},
        chemical_roles={"chebi.role:33229": {
            "name": "vitamin (role)", "parents": [], "level": 0, "level_kind": "chebi_role",
            "all_identifiers": ["chebi.role:33229", "chebi:10004", "chebi:33229"]}},
        subtree_roles={"2.A.1": ["chebi.role:33229"]},
        rekeys={"chebi:16977": f"{K}C00041"},
        role_of_old={"chebi:10004": "chebi.role:33229"},
        tables=tables,
        paper_overrides={},
        name_candidates=name_candidates,
        lock={},
        lumping_pin=None,
        identity_overrides=IdentityOverrides([], set()),
        chebi=_chebi(),
        chebi_from_raw=True,
        rxn_genes={"R1": {"s:G1", "s:G2"}, "R2": {"s:G3"}, "R3": {"s:G4"}},
        # ko00999 is reachable only through D-Trp, a paper candidate that stage K drops.
        compound_to_pathways={"C00525": ["ko00999"], "C00078": ["ko00380"], "C00041": ["ko00250"]},
        base_pathways={"ko00250", "ko00380"},
        brite_raw=brite_raw,
        brite_committed={},
        tcdb_kept={"2.A.1", "3.A.3"},
        tcdb_levels={"2.A.1": 2, "3.A.3": 2},
        tcdb_static={"kept_tcdb_ids": ["2.A.1", "3.A.3"], "pfam_bridge": {}, "go_bridge": {},
                     "seed_aliases": {}},
        previous_ids=set(table) | {"chebi:10004", "chebi:16977"},
        old_resolved={(PAPER, TABLE, 0): f"{K}C00133", (PAPER, TABLE, 1): f"{K}C00525"},
        today="2026-10-15",
    )
