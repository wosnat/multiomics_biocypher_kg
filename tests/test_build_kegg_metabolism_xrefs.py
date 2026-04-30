"""Unit tests for step 6: unified pruned kegg_data.json builder."""
from __future__ import annotations
import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def _make_resolver(tmp_path):
    """Minimal MNX resolver with synthetic compounds and reactions."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
            reference TEXT, classifs TEXT, is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'I=', 'WQZ-...', 'OC[CH]1...', 0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO reactions VALUES
            ('MNXR101234', '...', 'kegg.reaction:R00200', '2.7.1.40', 'B', NULL);
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
    """)
    conn.commit()
    return conn


@pytest.fixture
def synthetic_kegg_raw(tmp_path):
    """Return a kegg_data dict matching the parser output, suitable as input to step 6."""
    return {
        "ko_names":               {"K02338": "DNA polymerase III subunit beta",
                                   "K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00710": "Carbon fixation",
                                   "ko03030": "DNA replication",
                                   "ko00500": "Starch and sucrose metabolism"},
        "ko_to_pathways":         {"K02338": ["ko03030"], "K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00710": "09102",
                                   "ko03030": "09124", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism",
                                   "09102": "Energy metabolism",
                                   "09124": "Replication and repair"},
        "subcategory_to_category":{"09101": "09100", "09102": "09100", "09124": "09120"},
        "category_names":         {"09100": "Metabolism", "09120": "Genetic Information Processing"},
        "reaction_names":         {"R00200": "ATP:pyruvate ..."},
        "compound_names":         {"C00031": "D-glucose", "C00074": "PEP",
                                   "C99999": "compound-only-pathway-trigger"},
        "reaction_to_compounds":  {"R00200": ["C00031", "C00074"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"],
                                   "C99999": ["ko99999"]},   # ko99999 is compound-only
    }


def test_build_pruned_kegg_data(tmp_path, monkeypatch, synthetic_kegg_raw):
    """End-to-end: 1 strain, 2 genes, ensure pruned structure is correct."""
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
        "gene2": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])

    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, out_path)

    data = json.loads(out_path.read_text())

    # KOs: only gene-annotated KOs survive
    assert set(data["kos"].keys()) == {"K02338", "K01006"}
    assert data["kos"]["K02338"]["name"] == "DNA polymerase III subunit beta"
    assert data["kos"]["K02338"]["pathways"] == ["ko03030"]

    # Pathways: KO ∪ Rxn-reachable. K02338→ko03030, K01006→ko00010, R00200→ko00010
    # ko00710 not present — no KO/Rxn maps there.
    assert set(data["pathways"].keys()) == {"ko03030", "ko00010"}

    # Reactions: only gene-reachable
    assert set(data["reactions"].keys()) == {"R00200"}
    rxn = data["reactions"]["R00200"]
    assert rxn["name"].startswith("ATP:pyruvate")
    assert rxn["pathways"] == ["ko00010"]
    assert rxn["compounds"] == ["C00031", "C00074"]
    assert rxn["mnxr_id"] == "MNXR101234"
    assert rxn["ec_numbers"] == ["2.7.1.40"]
    assert rxn["mass_balance"] == "balanced"
    assert rxn["reaction_class"] == "chemical"

    # Compounds: only reaction-reachable
    assert set(data["compounds"].keys()) == {"C00031", "C00074"}
    glucose = data["compounds"]["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"

    # Option B: compound.pathways pruned to KO ∪ Rxn-reachable set.
    # ko00500 is NOT in the set (no KO or Rxn maps there in the fixture).
    # ko00010 IS in the set.
    assert glucose["pathways"] == ["ko00010"]

    # Subcategories/categories present iff parent of an emitted pathway
    assert set(data["subcategories"].keys()) == {"09101", "09124"}
    assert data["subcategories"]["09101"]["category"] == "09100"
    assert set(data["categories"].keys()) == {"09100", "09120"}


def test_compound_only_pathways_dropped(tmp_path, monkeypatch, synthetic_kegg_raw):
    """The C99999→ko99999 compound-only pathway should not appear anywhere.

    Critical for Option B: even if a compound were in the pruned set, its pathways
    should be filtered to KO ∪ Rxn-reachable. C99999 isn't even in the pruned set
    (no reaction reaches it), so this is a double-defensive check.
    """
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, out_path)
    data = json.loads(out_path.read_text())
    assert "C99999" not in data["compounds"]
    assert "ko99999" not in data["pathways"]


def test_compound_pathway_filter_drops_unevidenced(tmp_path, monkeypatch):
    """Compound pathways are filtered to KO ∪ Rxn-reachable.

    Force a case where a compound IS in the gene-reachable set (via R00200) but has
    a pathway (ko00500) that's not reachable from any gene-KO or gene-reaction.
    Expected: compound emitted, but ko00500 NOT in compound.pathways.
    """
    raw = {
        "ko_names":               {"K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00500": "Starch and sucrose"},
        "ko_to_pathways":         {"K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism"},
        "subcategory_to_category":{"09101": "09100"},
        "category_names":         {"09100": "Metabolism"},
        "reaction_names":         {"R00200": "rxn"},
        "compound_names":         {"C00031": "D-glucose"},
        "reaction_to_compounds":  {"R00200": ["C00031"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"]},  # ko00500 unevidenced
    }
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "g1": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(raw, conn, out_path)
    data = json.loads(out_path.read_text())

    # ko00500 must NOT be in pathways top-level (Option B: only KO ∪ Rxn-reachable)
    assert "ko00500" not in data["pathways"]
    # And compound's pathway list filtered
    assert data["compounds"]["C00031"]["pathways"] == ["ko00010"]
