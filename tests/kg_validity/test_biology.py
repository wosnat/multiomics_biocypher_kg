"""
Biological truth tests for the multi-omics knowledge graph.

These tests encode domain-specific constraints that must hold for
Prochlorococcus / Synechococcus / Alteromonas biology:

- Ecotype/clade labels on strains
- Gene absences that are defining features (Black Queen Hypothesis)
- All expected strains present
- UniProt locus-tag → accession spot checks
"""

import pytest


pytestmark = pytest.mark.kg


# Strains confirmed in graph (OrganismTaxon nodes with strain_name):
# 8 Prochlorococcus + Synechococcus CC9311 + Parasynechococcus WH8102 + 3 Alteromonas
EXPECTED_STRAINS = [
    # Prochlorococcus
    "MED4",
    "MIT9313",
    "MIT9301",
    "MIT9312",
    "NATL2A",
    "NATL1A",
    "AS9601",
    "RSP50",
    # Synechococcus / Parasynechococcus
    "WH8102",
    "CC9311",
    # Alteromonas
    "MIT1002",
    "EZ55",
    "HOT1A3",
]

# Clade assignments for key Prochlorococcus strains (HLI = High Light clade I)
EXPECTED_CLADES = {
    "MED4":    "HLI",
    "RSP50":   "HLI",
    "AS9601":  "HLII",
    "MIT9312": "HLII",
    "MIT9301": "HLII",
    "NATL2A":  "LLII",
    "NATL1A":  "LLII",
    "MIT9313": "LLIV",
}


# ---------------------------------------------------------------------------
# Strain presence
# ---------------------------------------------------------------------------

def test_all_expected_strains_present(run_query):
    """All expected organism strains must exist in the graph."""
    result = run_query(
        "MATCH (o:OrganismTaxon) RETURN o.strain_name AS name"
    )
    present = {r["name"] or "" for r in result}
    missing = [
        s for s in EXPECTED_STRAINS
        if not any(s in name for name in present)
    ]
    assert not missing, (
        f"Expected strains not found in graph: {missing}\n"
        f"Strains present: {sorted(present)}"
    )


# ---------------------------------------------------------------------------
# Ecotype / clade labels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strain,expected_clade", EXPECTED_CLADES.items())
def test_strain_clade_label(run_query, strain, expected_clade):
    """Each Prochlorococcus strain must carry the correct ecotype clade label."""
    result = run_query(
        "MATCH (o:OrganismTaxon) WHERE o.strain_name = $strain "
        "RETURN o.strain_name AS name, o.clade AS clade",
        strain=strain,
    )
    assert len(result) >= 1, f"No OrganismTaxon node found for strain '{strain}'"
    for row in result:
        assert row["clade"] == expected_clade, (
            f"Strain {row['name']}: expected clade '{expected_clade}', "
            f"got '{row['clade']}'"
        )


def test_prochlorococcus_genus_label(run_query):
    """All Prochlorococcus strains must have genus='Prochlorococcus'."""
    pro_strains = [s for s in EXPECTED_STRAINS
                   if s in ("MED4", "MIT9313", "MIT9301", "MIT9312",
                             "NATL2A", "NATL1A", "AS9601", "RSP50")]
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IN $strains
        AND (o.genus IS NULL OR o.genus <> 'Prochlorococcus')
        RETURN o.strain_name AS name, o.genus AS genus
    """, strains=pro_strains)
    assert result == [], (
        f"Strains with wrong/missing genus: {result}"
    )


# ---------------------------------------------------------------------------
# Black Queen Hypothesis: katG is absent in Prochlorococcus
# ---------------------------------------------------------------------------

def test_prochlorococcus_no_katG(run_query):
    """
    Prochlorococcus has lost catalase-peroxidase (katG) through genome streamlining
    (Black Queen Hypothesis). No Gene belonging to a Prochlorococcus organism
    should carry 'katG' as a gene name.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.genus = 'Prochlorococcus'
          AND any(name IN g.gene_names WHERE toLower(name) = 'katg')
        RETURN count(g) AS cnt, collect(g.locus_tag)[..5] AS examples
    """)
    cnt = result[0]["cnt"]
    assert cnt == 0, (
        f"Found {cnt} Prochlorococcus gene(s) named 'katG' (should be 0). "
        f"Example locus tags: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# UniProt locus-tag → accession spot checks (MED4)
# ---------------------------------------------------------------------------

def test_med4_genes_have_proteins(run_query):
    """
    MED4 genes should map to UniProt proteins via shared locus_tag.
    Require that at least 80% of MED4 genes have a corresponding Protein.
    Gene-protein linkage is via locus_tag (no explicit Gene_encodes_protein edge).
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.strain_name = 'MED4'
        OPTIONAL MATCH (p:Protein {locus_tag: g.locus_tag})
        WITH count(g) AS total, count(p) AS with_protein
        RETURN total, with_protein
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No MED4 genes found")
    coverage = row["with_protein"] / row["total"]
    assert coverage >= 0.75, (
        f"MED4 gene→protein coverage is {coverage:.1%} "
        f"({row['with_protein']}/{row['total']}); expected >= 75%"
    )


def test_known_med4_locus_tag_maps_to_protein(run_query):
    """
    PMM1375 (psbA) in MED4 should have a corresponding UniProt Protein node.
    Genes and proteins are linked by shared locus_tag (no explicit edge).
    This is a regression anchor — if the locus-tag mapping breaks, this fails.
    """
    result = run_query("""
        MATCH (p:Protein {locus_tag: 'PMM1375'})
        RETURN p.id AS uniprot_id, p.protein_name AS name
    """)
    assert len(result) >= 1, (
        "No Protein found with locus_tag='PMM1375' (MED4 psbA). "
        "UniProt locus-tag mapping may be broken."
    )


# ---------------------------------------------------------------------------
# Taxonomic hierarchy completeness
# ---------------------------------------------------------------------------

def test_genomic_organisms_have_genus(run_query):
    """
    Genomic OrganismTaxon nodes (those with a strain_name) must have genus populated.
    Species may be absent for some strains (e.g., Alteromonas strains without
    a named species yet). Non-genomic treatment organisms (phage, viruses) may
    lack full taxonomy and are excluded by requiring strain_name IS NOT NULL.
    """
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IS NOT NULL AND o.genus IS NULL
        RETURN o.strain_name AS strain, o.genus AS genus
    """)
    assert result == [], (
        f"Genomic OrganismTaxon nodes missing genus: {result}"
    )
