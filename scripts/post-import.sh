#!/bin/bash
set -euo pipefail

echo "=== Post-process: Starting Neo4j ==="
neo4j start
sleep 15

echo "=== Post-process: Running post-import.cypher ==="
cypher-shell -f /scripts/post-import.cypher

echo "=== Post-process: Edge counts ==="
cypher-shell "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS homolog_edges;"
cypher-shell "MATCH ()-[r:Affects_expression_of_homolog]->() RETURN count(r) AS homolog_expression_edges;"

echo "=== Post-process complete ==="
neo4j stop
