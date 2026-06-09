"""Publication "discusses" topic extraction.

Extracts the genes and KEGG pathways a paper substantively discusses in its
narrative (not just the supplementary DE tables), for the
``Publication_discusses_gene`` / ``Publication_discusses_kegg_pathway`` edges.

Stage 1 of the extract -> resolve -> adapter pipeline. See
``plans/publication_discusses_edges.md``.
"""
