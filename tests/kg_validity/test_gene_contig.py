"""KG validity: every Gene must have a non-null contig.

F4 ship 4.1 requires the contig (NCBI GFF seqid) on every Gene node so
genomic-neighbor lookups have a valid discriminator.
"""

import pytest


pytestmark = pytest.mark.kg


def test_every_gene_has_contig(run_query):
    """Every Gene node must have a non-null contig property."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.contig IS NULL OR g.contig = ''
        RETURN count(*) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} Gene nodes have null/empty contig"


def test_contig_consistent_within_organism(run_query):
    """All genes from a single organism should map to at least one contig.

    Single-chromosome organisms have 1 contig; multi-replicon organisms have
    more, but every organism must have at least one.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WITH o.preferred_name AS organism, collect(DISTINCT g.contig) AS contigs
        WHERE size(contigs) = 0
        RETURN organism, contigs
    """)
    assert result == [], f"organisms with no distinct contig: {result}"


def test_med4_has_known_chromosome(run_query):
    """Spot-check: MED4 genes should land on its known NCBI accession."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.preferred_name = 'Prochlorococcus MED4'
        RETURN DISTINCT g.contig AS contig
    """)
    contigs = {row["contig"] for row in result}
    # MED4 single-replicon — expect one chromosome accession
    assert len(contigs) == 1, f"MED4 expected 1 contig, got {len(contigs)}: {contigs}"
