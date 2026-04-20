// Biller 2018 retrofit-removal baseline capture.
//
// Publication: Biller et al. 2018 (mSystems), DOI 10.1128/mSystems.00040-18
// Purpose: snapshot KG state BEFORE removing the three retrofitted
// ClusteringAnalysis nodes (natl2a_periodicity, mit1002_periodicity,
// natl2a_darkness_survival) so Plan 2 can assert that the retrofit removal
// deleted exactly these cluster edges and replaced them with DerivedMetric
// edges.
//
// Note on relationship labels: the clustering-hierarchy edges in the live KG
// use PascalCase without underscores (`PublicationHasClusteringAnalysis`,
// `ClusteringAnalysisHasGeneCluster`, `ExperimentHasClusteringAnalysis`),
// while `Gene_in_gene_cluster` uses underscore form. This script uses the
// labels as they appear in the live graph (verified via `type(r)` lookup).
//
// To re-run:
//   cypher-shell -a bolt://localhost:7687 --format plain \
//     -P "pub_id => \"doi:10.1128/mSystems.00040-18\"" \
//     < docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher

:param pub_id => "doi:10.1128/mSystems.00040-18";

MATCH (p:Publication {id: $pub_id})
RETURN 'publication_exists' AS metric, count(p) AS value
UNION ALL
MATCH (p:Publication {id: $pub_id})-[:Has_experiment]->(e:Experiment)
RETURN 'experiment_count' AS metric, count(e) AS value
UNION ALL
MATCH (p:Publication {id: $pub_id})-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(:Gene)
RETURN 'changes_expression_of_edges' AS metric, count(r) AS value
UNION ALL
MATCH (p:Publication {id: $pub_id})-[:PublicationHasClusteringAnalysis]->(ca:ClusteringAnalysis)
RETURN 'clustering_analysis_count' AS metric, count(ca) AS value
UNION ALL
MATCH (ca:ClusteringAnalysis)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_clustering_analysis_count' AS metric, count(ca) AS value
UNION ALL
MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_gene_cluster_count' AS metric, count(gc) AS value
UNION ALL
MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(:GeneCluster)-[r:Gene_in_gene_cluster]->(:Gene)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_gene_in_gene_cluster_edges' AS metric, count(r) AS value
UNION ALL
MATCH (:Experiment)-[r:ExperimentHasClusteringAnalysis]->(ca:ClusteringAnalysis)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_experiment_has_ca_edges' AS metric, count(r) AS value;
