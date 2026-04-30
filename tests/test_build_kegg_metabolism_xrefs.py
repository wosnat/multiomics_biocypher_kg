"""Unit tests for step 6: build_kegg_metabolism_xrefs."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def test_collect_gene_reachable_kegg_ids(tmp_path, monkeypatch):
    # Two synthetic strains
    s1 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MED4"
    s2 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MIT9301"
    _write_strain_annotations(s1, {
        "PMM0001": {"kegg_reactions": ["R00200", "R00010"], "kegg_ko": ["K00001"]},
        "PMM0002": {"kegg_reactions": ["R00010"]},
        "PMM0003": {},  # no reactions, no KOs
    })
    _write_strain_annotations(s2, {
        "P9301_001": {"kegg_reactions": ["R12345"]},
    })

    # Stub: load_genome_rows returns rows pointing at the two strains
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [
        {"data_dir": str(s1)},
        {"data_dir": str(s2)},
    ])

    # Stub kegg_data with reaction_to_compounds + ko_to_pathways maps
    kegg_data = {
        "reaction_names": {
            "R00200": "pyruvate kinase reaction",
            "R00010": "phosphofructokinase",
            "R12345": "test rxn",
        },
        "compound_names": {
            "C00074": "phosphoenolpyruvate",
            "C00008": "ADP",
            "C00031": "D-glucose",
            "C00002": "ATP",
            "C99999": "test cpd",
        },
        "reaction_to_compounds": {
            "R00200": ["C00074", "C00008"],
            "R00010": ["C00031", "C00002"],
            "R12345": ["C99999"],
            "R_unused": ["C_should_not_appear"],
        },
        "reaction_to_pathways": {
            "R00200": ["ko00010"],
        },
        "ko_to_pathways": {
            "K00001": ["ko00010", "ko00710"],
            "K99999": ["ko99999"],  # KO not in any gene → pathway pruned
        },
    }

    rxns, cpds, pws = bx._collect_gene_reachable_ids(kegg_data)

    # Only gene-reachable R-numbers
    assert rxns == {"R00200", "R00010", "R12345"}
    # C-numbers reachable from those reactions
    assert cpds == {"C00074", "C00008", "C00031", "C00002", "C99999"}
    # R_unused / C_should_not_appear are pruned
    assert "R_unused" not in rxns
    assert "C_should_not_appear" not in cpds
    # KO-reachable pathways: K00001 maps to ko00010 + ko00710; K99999 not in any gene
    assert pws == {"ko00010", "ko00710"}
    assert "ko99999" not in pws


def test_collect_handles_missing_kegg_reactions_field(tmp_path, monkeypatch):
    s = tmp_path / "strain"
    _write_strain_annotations(s, {
        "g1": {},  # no kegg_reactions key, no kegg_ko key
        "g2": {"kegg_reactions": []},  # empty list
        "g3": {"kegg_reactions": ["R00001"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(s)}])
    kegg_data = {"reaction_to_compounds": {"R00001": ["C00001"]}, "ko_to_pathways": {}}
    rxns, cpds, pws = bx._collect_gene_reachable_ids(kegg_data)
    assert rxns == {"R00001"}
    assert cpds == {"C00001"}
    assert pws == set()


import sqlite3


def test_enrich_reaction(tmp_path):
    """Reaction enrichment: KEGG name + pathway + MNX cross-refs."""
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

        INSERT INTO reactions VALUES
            ('MNXR101234', '1 C00074 + 1 C00008 = 1 C00022 + 1 C00002',
             'kegg.reaction:R00200', '2.7.1.40', 'B', NULL);
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
        INSERT INTO reaction_aliases VALUES ('rhea', '10828', 'MNXR101234');

        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'InChI=1/C6H12O6',
             'WQZGKKKJIJFFOK-GASJEMHNSA-N', 'OC[CH]1OC(O)[CH](O)[CH](O)[CH]1O',
             0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('hmdb', 'HMDB0000122', 'MNXM41');
    """)
    conn.commit()

    kegg_data = {
        "reaction_names": {"R00200": "pyruvate kinase reaction"},
        "compound_names": {"C00031": "D-glucose"},
        "reaction_to_compounds": {"R00200": ["C00074", "C00008"]},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
    }

    # Both pathways are valid — all survive
    rxn = bx._enrich_reaction("R00200", kegg_data, conn, valid_pathways={"ko00010", "ko00710"})
    assert rxn["name"] == "pyruvate kinase reaction"
    assert rxn["compound_ids"] == ["C00074", "C00008"]
    assert rxn["kegg_pathway_ids"] == ["ko00010", "ko00710"]
    assert rxn["mnxr_id"] == "MNXR101234"
    assert rxn["rhea_ids"] == ["10828"]
    assert rxn["mass_balance"] == "balanced"
    assert rxn["reaction_class"] == "chemical"
    assert rxn["ec_numbers"] == ["2.7.1.40"]


def test_enrich_reaction_prunes_dangling_pathways(tmp_path):
    """Pathway IDs not in valid_pathways are pruned to avoid dangling edges."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
            reference TEXT, classifs TEXT, is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
    """)
    conn.commit()

    kegg_data = {
        "reaction_names": {"R00200": "pyruvate kinase reaction"},
        "compound_names": {},
        "reaction_to_compounds": {"R00200": []},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
    }

    # Only ko00010 is KO-reachable; ko00710 should be pruned
    rxn = bx._enrich_reaction("R00200", kegg_data, conn, valid_pathways={"ko00010"})
    assert rxn["kegg_pathway_ids"] == ["ko00010"]
    assert "ko00710" not in rxn["kegg_pathway_ids"]


def test_enrich_compound(tmp_path):
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'InChI=...',
             'WQZGKKKJIJFFOK-GASJEMHNSA-N', 'OC[CH]1...', 0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('hmdb', 'HMDB0000122', 'MNXM41');
    """)
    conn.commit()

    kegg_data = {"compound_names": {"C00031": "D-glucose"}}
    cpd = bx._enrich_compound("C00031", kegg_data, conn)
    assert cpd["name"] == "D-glucose"
    assert cpd["formula"] == "C6H12O6"
    assert cpd["mass"] == 180.063
    assert cpd["inchikey"] == "WQZGKKKJIJFFOK-GASJEMHNSA-N"
    assert cpd["mnxm_id"] == "MNXM41"
    assert cpd["chebi_id"] == "17234"
    assert cpd["hmdb_id"] == "HMDB0000122"


def test_enrich_compound_orphan_no_mnx_match(tmp_path):
    """C-number with no MNX entry: enrichment returns just KEGG name."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
    """)
    conn.commit()

    kegg_data = {"compound_names": {"C99999": "obscure compound"}}
    cpd = bx._enrich_compound("C99999", kegg_data, conn)
    assert cpd["name"] == "obscure compound"
    assert cpd.get("mnxm_id") is None
    assert cpd.get("formula") is None
