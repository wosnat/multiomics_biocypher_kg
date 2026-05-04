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
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    data = bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, sets=sets)

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
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    data = bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, sets=sets)
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
    sets = bx._gene_reachable_sets(raw)
    data = bx.build_pruned_kegg_data(raw, conn, sets=sets)

    # ko00500 must NOT be in pathways top-level (Option B: only KO ∪ Rxn-reachable)
    assert "ko00500" not in data["pathways"]
    # And compound's pathway list filtered
    assert data["compounds"]["C00031"]["pathways"] == ["ko00010"]


def test_bulk_enrich_reactions_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_reactions(conn, kegg_ids, allowed_pathways, raw) returns
    a dict mapping kegg_reaction_id → enrichment dict, with the same fields
    as the per-entity _enrich_reaction.
    """
    conn = _make_resolver(tmp_path)
    raw = {
        "reaction_names": {"R00200": "ATP:pyruvate ..."},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
        "reaction_to_compounds": {"R00200": ["C00031"]},
    }
    allowed = {"ko00010", "ko00710"}

    result = bx._bulk_enrich_reactions(conn, ["R00200", "R99999"], allowed, raw)

    # Every requested ID is in the result (R99999 has no MNX entry but still gets a stub)
    assert set(result.keys()) == {"R00200", "R99999"}

    r = result["R00200"]
    assert r["name"].startswith("ATP:pyruvate")
    assert r["mnxr_id"] == "MNXR101234"
    assert r["ec_numbers"] == ["2.7.1.40"]
    assert r["mass_balance"] == "balanced"
    assert r["reaction_class"] == "chemical"
    assert r["pathways"] == ["ko00010", "ko00710"]
    assert r["compounds"] == ["C00031"]

    r99 = result["R99999"]
    assert r99["mnxr_id"] is None
    assert r99["ec_numbers"] == []
    assert r99["mass_balance"] == "unbalanced"  # default
    assert r99["reaction_class"] == "chemical"  # default


def test_step6_builds_tcdb_hierarchy(tmp_path):
    """Step 6 must produce cache/data/tcdb/tcdb_hierarchy.json (was previously
    built by build_metabolite_resolver.py)."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    # Stage minimal TCDB TSVs the way download_metabolism_reference.download_all would:
    cache_root = tmp_path / "cache" / "data"
    tcdb_dir = cache_root / "tcdb"
    raw_dir = tcdb_dir / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "families.tsv").write_text("1.A.1\tThe Voltage-gated Ion Channel\n")
    (raw_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (raw_dir / "substrates.tsv").write_text(
        "1.A.1.5.2\tCHEBI:9314;sucrose|CHEBI:3308;calcium(2+)\n"
    )
    (raw_dir / "acc2tcid.tsv").write_text("P12345\t1.A.1.5.2\n")

    out = tcdb_dir / "tcdb_hierarchy.json"
    mod._build_tcdb_hierarchy(cache_root=cache_root)

    assert out.exists()
    import json
    h = json.loads(out.read_text())
    # 5 levels for 1.A.1.5.2 → expect class+subclass+family+subfamily+specificity
    assert "1" in h
    assert "1.A" in h
    assert "1.A.1" in h
    assert "1.A.1.5" in h
    assert "1.A.1.5.2" in h
    assert h["1.A.1.5.2"]["level_kind"] == "tc_specificity"
    # substrate_classes preserves the full 'CHEBI:N;name' string so downstream
    # MNX resolution can map to a canonical Metabolite primary ID.
    assert h["1.A.1.5.2"]["substrate_classes"] == [
        "CHEBI:9314;sucrose",
        "CHEBI:3308;calcium(2+)",
    ]


def test_build_tcdb_hierarchy_seeds_3part_families_without_5part_children(tmp_path):
    """Family-level entries in families.tsv must land in the hierarchy even
    when superfamilies.tsv has no 5-part child for them. Regression for the
    1.A.11 (Ammonium Transporter Channel) class of dropped seeds."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    raw_dir = cache_root / "tcdb" / "raw"
    raw_dir.mkdir(parents=True)
    # 1.A.11 has a description but no curated 5-part children
    (raw_dir / "families.tsv").write_text(
        "1.A.1\tThe Voltage-gated Ion Channel\n"
        "1.A.11\tThe Ammonium Transporter Channel (Amt) Family\n"
    )
    (raw_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (raw_dir / "substrates.tsv").write_text("")
    (raw_dir / "acc2tcid.tsv").write_text("")

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    assert "1.A.11" in h
    assert h["1.A.11"]["level_kind"] == "tc_family"
    assert h["1.A.11"]["name"] == "The Ammonium Transporter Channel (Amt) Family"
    assert h["1.A.11"]["parent"] == "1.A"


def test_build_tcdb_hierarchy_seeds_5part_from_acc2tcid(tmp_path):
    """5-part TCIDs present only in acc2tcid.tsv (not in superfamilies.tsv)
    must land in the hierarchy with synthesized ancestors. Regression for
    1.A.12.2.2 / 1.A.30.2.1 / 2.A.23.1.1 class of dropped seeds."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    raw_dir = cache_root / "tcdb" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "families.tsv").write_text("1.A.12\tThe Intracellular Chloride Channel\n")
    (raw_dir / "superfamilies.tsv").write_text("")
    # The TCID has a substrate even though no superfamilies row exists
    (raw_dir / "substrates.tsv").write_text("1.A.12.2.2\tCHEBI:3731;chloride\n")
    (raw_dir / "acc2tcid.tsv").write_text("Q9X9X9\t1.A.12.2.2\n")

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    # Full 5-level chain synthesized
    assert {"1", "1.A", "1.A.12", "1.A.12.2", "1.A.12.2.2"} <= set(h.keys())
    assert h["1.A.12.2.2"]["level_kind"] == "tc_specificity"
    assert h["1.A.12.2.2"]["substrate_classes"] == ["CHEBI:3731;chloride"]
    assert h["1.A.12.2"]["level_kind"] == "tc_subfamily"
    assert h["1.A.12"]["name"] == "The Intracellular Chloride Channel"


def test_bulk_enrich_compounds_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_compounds returns a dict mapping kegg_compound_id → enrichment."""
    conn = _make_resolver(tmp_path)
    raw = {
        "compound_names": {"C00031": "D-glucose", "C99999": "obscure"},
        "compound_to_pathways": {
            "C00031": ["ko00010", "ko00500"],
            "C99999": [],
        },
    }
    allowed = {"ko00010"}

    result = bx._bulk_enrich_compounds(conn, ["C00031", "C99999"], allowed, raw)

    assert set(result.keys()) == {"C00031", "C99999"}

    glucose = result["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["pathways"] == ["ko00010"]  # ko00500 filtered out by allowed

    obscure = result["C99999"]
    assert obscure["mnxm_id"] is None
    assert obscure["chebi_id"] is None
    assert obscure["formula"] is None


# ── Step 6 part 2: TCDB substrate resolution + evidence_sources tagging ──────


def test_compounds_get_metabolism_evidence_source():
    """A compound in catalysis_cpds with no transport overlap gets
    evidence_sources=['metabolism']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {
        "compounds": {
            "C00031": {"name": "D-Glucose", "formula": "C6H12O6"},
            "C00002": {"name": "ATP", "formula": "C10H16N5O13P3"},
        },
        "reactions": {},
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds={"C00031", "C00002"},
        substrate_kegg_cpds=set(),
        leaf_to_primary={},
        compound_props={},
    )
    for cpd in kegg_data["compounds"].values():
        assert cpd["evidence_sources"] == ["metabolism"]


def test_transport_only_compounds_land_in_additional_compounds():
    """A non-KEGG TCDB substrate primary becomes an additional_compounds entry
    with evidence_sources=['transport']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {"compounds": {}, "reactions": {}}
    leaf_to_primary = {"1.A.1.5.2": ["chebi:9999"]}
    compound_props = {
        "chebi:9999": {
            "name": "tetracycline", "formula": "C22H24N2O8",
            "mass": 444.43, "inchikey": None,
            "mnxm_id": "MNXM00099", "chebi_id": "9999",
        }
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds=set(),
        substrate_kegg_cpds=set(),
        leaf_to_primary=leaf_to_primary,
        compound_props=compound_props,
    )
    assert "chebi:9999" in kegg_data["additional_compounds"]
    entry = kegg_data["additional_compounds"]["chebi:9999"]
    assert entry["evidence_sources"] == ["transport"]
    assert entry["chebi_id"] == "9999"


def test_overlap_compound_gets_both_evidence_sources():
    """A compound in BOTH catalysis_cpds AND substrate_kegg_cpds gets
    evidence_sources=['metabolism', 'transport']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {
        "compounds": {"C00089": {"name": "Sucrose", "formula": "C12H22O11"}},
        "reactions": {},
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds={"C00089"},
        substrate_kegg_cpds={"C00089"},
        leaf_to_primary={"1.A.1.5.2": ["kegg.compound:C00089"]},
        compound_props={},
    )
    # KEGG-flavor substrate stays out of additional_compounds — it's already in compounds.
    assert "kegg.compound:C00089" not in kegg_data["additional_compounds"]
    assert "C00089" not in kegg_data["additional_compounds"]
    assert kegg_data["compounds"]["C00089"]["evidence_sources"] == ["metabolism", "transport"]


def test_transport_only_kegg_compound_in_compounds_with_pathways(tmp_path, monkeypatch):
    """Refactor invariant: when a kegg.compound primary is in substrate_kegg_cpds
    but NOT in catalysis_cpds, build_pruned_kegg_data (called with the extended
    cpds set) puts the entry in compounds (not additional_compounds) AND populates
    its pathways field with the gene-reachable pathways (extended pws set).

    Setup: the substrate compound C00089 participates in 2 KEGG pathways. Only
    one (ko00010) is in the extended pws; ko00500 is excluded. Verify the entry
    lands in compounds with pathways=['ko00010'].
    """
    raw = {
        "ko_names":               {"K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00500": "Starch and sucrose"},
        "ko_to_pathways":         {"K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism"},
        "subcategory_to_category": {"09101": "09100"},
        "category_names":         {"09100": "Metabolism"},
        "reaction_names":         {},
        "compound_names":         {"C00089": "Sucrose"},
        "reaction_to_compounds":  {},
        "reaction_to_pathways":   {},
        "compound_to_pathways":   {"C00089": ["ko00010", "ko00500"]},
    }
    # Catalysis side: 1 KO (ko00010), no reactions, no compounds.
    catalysis_sets = {
        "kos": {"K01006"},
        "rxns": set(),
        "cpds": set(),
        "pws": {"ko00010"},
        "tcdb_ids": set(),
    }
    # Caller-extended: C00089 is a transport-reachable kegg compound. Curate
    # extended_pws = {"ko00010"} so ko00500 stays excluded; we want to confirm
    # the compound's pathways field is filtered to that extended set.
    extended_sets = {
        **catalysis_sets,
        "cpds": {"C00089"},
        "pws": {"ko00010"},
    }

    conn = sqlite3.connect(":memory:")
    conn.executescript("""
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, formula TEXT, mass REAL,
            inchikey TEXT, smiles TEXT);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, classifs TEXT,
            is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
    """)
    data = bx.build_pruned_kegg_data(raw, conn, sets=extended_sets)
    # Entry lands in compounds (not additional_compounds — the latter is created
    # by _fold_substrates_into_kegg_data, not build_pruned_kegg_data).
    assert "C00089" in data["compounds"]
    # Pathways filtered to extended pws — ko00500 dropped.
    assert data["compounds"]["C00089"]["pathways"] == ["ko00010"]
    # ko00500 does not appear in top-level pathways either.
    assert "ko00500" not in data["pathways"]


def test_prune_tcdb_walks_up_and_down():
    """_prune_tcdb walks up to tc_class AND down to tc_specificity for each seed."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"name": "Channels", "level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"name": "VIC", "level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"name": "", "level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["calcium(2+)"]},
        "1.A.1.5.3": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["sodium(1+)"]},
        # An unrelated branch that should NOT be kept
        "2": {"name": "ECP", "level": 0, "level_kind": "tc_class", "parent": None},
        "2.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "2"},
    }
    # Seed at family level — walk up to 1, walk down to BOTH leaves
    kept, leaf_subs, seed_aliases = _prune_tcdb(hierarchy, {"1.A.1"})
    assert kept == {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2", "1.A.1.5.3"}
    assert leaf_subs == {
        "1.A.1.5.2": ["calcium(2+)"],
        "1.A.1.5.3": ["sodium(1+)"],
    }
    assert seed_aliases == {}
    # Branch 2 untouched
    assert "2" not in kept


def test_prune_tcdb_seed_at_leaf():
    """Seeding at a leaf still walks up to root."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["x"]},
    }
    kept, leaf_subs, seed_aliases = _prune_tcdb(hierarchy, {"1.A.1.5.2"})
    assert kept == {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2"}
    assert leaf_subs == {"1.A.1.5.2": ["x"]}
    assert seed_aliases == {}


def test_prune_tcdb_unknown_seed_walks_up_to_ancestor():
    """Seeds not in the hierarchy are remapped to the nearest curated ancestor
    (so retired eggNOG TCIDs collapse to family-level, not lost)."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "3": {"level": 0, "level_kind": "tc_class", "parent": None},
        "3.A": {"level": 1, "level_kind": "tc_subclass", "parent": "3"},
        "3.A.1": {"level": 2, "level_kind": "tc_family", "parent": "3.A"},
    }
    # 3.A.1.35 doesn't exist — should anchor onto 3.A.1.
    kept, leaf_subs, seed_aliases = _prune_tcdb(hierarchy, {"3.A.1.35"})
    assert seed_aliases == {"3.A.1.35": "3.A.1"}
    assert kept == {"3", "3.A", "3.A.1"}
    assert leaf_subs == {}


def test_prune_tcdb_unknown_seed_with_no_ancestor_marked_dropped():
    """When no truncation hits an ancestor, the seed is recorded with empty
    alias so the caller can warn and drop it."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
    }
    kept, leaf_subs, seed_aliases = _prune_tcdb(hierarchy, {"99.X.Y.Z.Q"})
    assert kept == set()
    assert leaf_subs == {}
    assert seed_aliases == {"99.X.Y.Z.Q": ""}


def test_resolve_substrates_strips_chebi_prefix(tmp_path):
    """Regression: substrate strings 'CHEBI:NNNN;name' must resolve through the
    MNX SQLite, which stores chebi aliases as (source='chebi', value='<NNNN>') —
    NOT 'CHEBI:NNNN'. Earlier resolver implementation passed the prefixed form
    to a value-only query, returning 0 hits for every substrate."""
    import sqlite3

    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrates

    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compound_aliases (
            source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id)
        );
        CREATE TABLE compounds (
            mnxm_id TEXT PRIMARY KEY,
            name TEXT, formula TEXT, mass REAL, inchikey TEXT
        );
    """)
    # Bare value, source='chebi' — matches MNX's storage convention.
    cur.execute("INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000')")
    cur.execute("INSERT INTO compounds VALUES ('MNXM50000', 'sucrose', 'C12H22O11', 342.30, 'CZMRCDWAGMRECN-UGDNZRGBSA-N')")
    conn.commit()

    leaf_to_primary, compound_props = _resolve_substrates(
        {"1.A.1.5.2": ["CHEBI:9314;sucrose"]},
        conn,
    )

    assert leaf_to_primary["1.A.1.5.2"] == ["chebi:9314"], (
        f"Expected ['chebi:9314'], got {leaf_to_primary['1.A.1.5.2']!r}. "
        "If empty, _resolve_substrates failed to strip the CHEBI: prefix."
    )
    entry = compound_props["chebi:9314"]
    assert entry["mnxm_id"] == "MNXM50000"
    assert entry["chebi_id"] == "9314"
    assert entry["name"] == "sucrose"
    assert entry["formula"] == "C12H22O11"


def _resolve_subs_db(tmp_path):
    """Shared resolver fixture for batched _resolve_substrates tests."""
    import sqlite3
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compound_aliases (
            source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id)
        );
        CREATE TABLE compounds (
            mnxm_id TEXT PRIMARY KEY,
            name TEXT, formula TEXT, mass REAL, inchikey TEXT
        );
    """)
    return conn


def test_resolve_substrates_dedupes_compound_props_across_leaves(tmp_path):
    """Two leaves citing the same CHEBI substrate yield one compound_props entry,
    and first-occurrence's name_part wins for the fallback name."""
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.execute("INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000')")
    # No `compounds` row → mnx name is None, so fallback = name_part from input
    conn.commit()

    leaf_to_primary, compound_props = bx._resolve_substrates(
        {
            "1.A.1.5.2": ["CHEBI:9314;sucrose"],
            "2.A.1.1.1": ["CHEBI:9314;table sugar"],
        },
        conn,
    )

    assert leaf_to_primary["1.A.1.5.2"] == ["chebi:9314"]
    assert leaf_to_primary["2.A.1.1.1"] == ["chebi:9314"]
    assert list(compound_props.keys()) == ["chebi:9314"]
    # First leaf in iteration order wins for the name fallback.
    assert compound_props["chebi:9314"]["name"] == "sucrose"


def test_resolve_substrates_kegg_primary_skips_compound_props(tmp_path):
    """When CHEBI resolves to an mnxm with a kegg.compound alias, the primary is
    `kegg.compound:Cnnnnn` and compound_props has NO entry — KEGG primaries flow
    through `_bulk_enrich_compounds` via the extended cpds set."""
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.executescript("""
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', 'C6H12O6', 180.06, 'WQZ-');
    """)
    conn.commit()

    leaf_to_primary, compound_props = bx._resolve_substrates(
        {"3.A.1.1.1": ["CHEBI:17234;glucose"]},
        conn,
    )

    assert leaf_to_primary["3.A.1.1.1"] == ["kegg.compound:C00031"]
    assert "kegg.compound:C00031" not in compound_props
    assert compound_props == {}


def test_resolve_substrates_unresolved_substrate_keeps_leaf_with_empty_list(tmp_path):
    """Unresolvable substrates are silently dropped, but every leaf passed in
    must still appear in the output (with possibly empty list) — downstream code
    iterates leaf_to_primary.values() expecting all leaves."""
    conn = _resolve_subs_db(tmp_path)
    # No aliases inserted — every substrate is unresolved.

    leaf_to_primary, compound_props = bx._resolve_substrates(
        {
            "1.X.X.X.X": ["CHEBI:99999;mystery"],
            "2.Y.Y.Y.Y": ["malformed_no_colon"],
        },
        conn,
    )

    assert leaf_to_primary == {"1.X.X.X.X": [], "2.Y.Y.Y.Y": []}
    assert compound_props == {}


def test_resolve_substrates_mixed_kegg_and_chebi_primaries(tmp_path):
    """Mixed kegg-priority + chebi-only primaries resolve correctly in one batch.

    (The `mnx:*` fallback in _primary_for is by-construction unreachable from
    _resolve_substrates: every inbound CHEBI match guarantees the mnxm has at
    least one chebi alias, so mnxm_to_chebi is always populated for nonkegg
    mnxm_ids. Defensive code, no test needed.)
    """
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.executescript("""
        -- glucose: kegg + chebi present → kegg.compound:* primary
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', 'C6H12O6', 180.06, '');

        -- sucrose: chebi only → chebi:* primary
        INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000');
        INSERT INTO compounds VALUES ('MNXM50000', 'sucrose', 'C12H22O11', 342.30, '');
    """)
    conn.commit()

    leaf_to_primary, compound_props = bx._resolve_substrates(
        {
            "1.A": ["CHEBI:17234;glucose"],
            "1.B": ["CHEBI:9314;sucrose"],
        },
        conn,
    )

    assert leaf_to_primary["1.A"] == ["kegg.compound:C00031"]
    assert leaf_to_primary["1.B"] == ["chebi:9314"]
    # Only non-KEGG primaries have compound_props
    assert set(compound_props.keys()) == {"chebi:9314"}
    assert compound_props["chebi:9314"]["mnxm_id"] == "MNXM50000"
    assert compound_props["chebi:9314"]["chebi_id"] == "9314"
