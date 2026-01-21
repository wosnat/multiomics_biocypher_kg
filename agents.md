# Goal of this file

This file contain instructions for AI agents for creating and structuring the code of this repo

# Goal of the repo
This repo contains code for building a knowledge graph that describes the known genomic information on the cyanobacteria Prochlorococcus and marine bacteria Alteromonas. and the interaction between them as measured in coculture experiments. These include information on genes, proteins and metabolites from public databases (uniprot, kegg, ec, cyanorak, and more). As well as results of RNASEQ, affymetrix arrays, proteomics and metablomics from this study and other publications.

The resulting knowledge graph will be fed to LLM agents and used to analyze fresh omics data.
We are using biocypher infrastructure which can create neo4j graphs.

# directory structure of the repo:
* create_knowledge_graph.py - the main script that build the graph.
* config/schema_config.yaml - the graph schema
* multiomics_kg/adapters - adapters used to build the graph
* data/Prochlorococcus/papers_and_supp - existing publications. each subfolder contains the PDF of a publication as well as relevant supplemental informantion files and tables.
* data/Prochlorococcus/genomes/MED4 - genomic information on Prochlorococcus MED4 - the first organism to insert into the graph.

# implementation strategy
Start with a graph for a single Prochlorococcus strain - Prochlorococcus MED4. Reuse and tailor adapters from existing biocypher implementations where possible. 
Leverage unit tests whereever possible to check continuous correction as more data is added.


Use a combination of LLM based entity extraction and summarization, automatic downloads, and manual coding to convert the files to the required format. cache intermediate results in json file to avoid LLM costs and time consuming reruns.


# for the omics result adapter
Use a config file per publication, located in the folder of the paper (under data\Prochlorococcus\papers_and_supp)
YAML format

List all of the tables that need to be loaded. 
For each table, what is the entities column and what it represent (e.g., genes, proteins, pathways) and which statistical tests are included. 
For each of these tests, add to the config all fields as listed in the schema: config\schema_config.yaml. Connect supp table columns to statistical test fields.

Use this config in the adapter to create the relevant nodes and edges.
