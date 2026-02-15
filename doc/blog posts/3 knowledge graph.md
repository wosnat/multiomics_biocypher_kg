How do you integrate Omics data into a GenAI LLM?
In previous posts, I described our mission: extending GenAI models with scientific omics data. I also explained our first design choice—working at the level of individual fold changes.

In this post, I will go over the different options for presenting this data to an LLM, their pros and cons, and the motivation behind our final choice: the Knowledge Graph.

Common integration options: Fine-tuning and RAG
There are two methods often used to teach an LLM new data: fine-tuning and RAG.

Fine-tuning involves retraining an existing model using a custom training set that includes knowledge and tasks targeting a specific application. This can be done for foundation models like GPT-4 or open-source models like Llama. Since the model is already pre-trained, we are simply "tweaking" it rather than starting from scratch. The result is a model with relevant knowledge built directly into its weights.

Retrieval-Augmented Generation (RAG) is a technique where customization happens in the prompt instead of changing the model itself. The user prompt is automatically augmented with relevant data injected from an external source, accompanied by an instruction to answer based on that data.

RAG development is often based on a vector database containing relevant text and data. This database is queried using keywords or semantic searches. Each model invocation involves a vector DB search, injection into the prompt, and feeding the revised prompt into the LLM.

Pros and Cons
Fine-tuning makes the model highly customized. It becomes intimately familiar with your data. However, if your data is dynamic, you must frequently retrain the model. This process is expensive, requiring significant expertise and computing resources.

RAG is more flexible and supports data changes by design. It is lightweight and easier to implement. However, the quality of the results is highly dependent on your ability to select the most relevant sections from the database. It also requires a larger context window and can run into the model’s token limits.

The Complexity of Multi-omics Data
The main elements in a cell include genomic data, enzymes, proteins, and metabolites. Extensive research on these exists in online databases like UniProt, NCBI, and ChEMBL. These elements are all interconnected in a complex web of interactions.

For example, a specific gene has an associated protein enzyme. That protein may be part of multiple metabolic pathways, have a specific cellular location, and be involved in various biological processes. Enzymes catalyze chemical reactions involving metabolites as both reactants and products.

Multi-omics assays study all of these elements:

RNA sequencing looks at cell adaptation by tracking genes transcribed into proteins.

Proteomics examines protein enzymes.

Metabolomics measures the concentration of metabolites.

DNA resequencing studies genotypes induced by mutations.

Because this data is a web of complicated relationships, our integration method must reflect that connectivity.

The Knowledge Graph Approach
There is a third approach to integration: instead of changing the model or the prompt, we can give the model a tool (similar to how models use web search) to acquire information.

Since omics data is highly interconnected, using a graph to represent it is an intriguing and powerful approach. This is the Knowledge Graph (KG)—an organized representation of real-world data where semantic interpretations are made explicit.

Take, for example, these Cypher queries (the language used to query graphs like Neo4j):

(Nutrient starvation)-[affects_expression_of {direction: 'up'}]->(gene)

(protein)-[involved_in_biological_process]->(pathway)

Complex questions become much easier to answer: "What are the genes in this pathway?" or "Are there stress-related genes that are upregulated in co-culture, but not under starvation?"

Knowledge Graph as a Tool
Each Knowledge Graph has a well-defined schema, which an LLM agent can use to form precise queries.

While "Graph RAG" (injecting graph results into a prompt) is an option, we found it less flexible because a single user prompt may require multiple, iterative graph queries. Instead, we elected to use a sub-agent.

This sub-agent receives a text prompt, translates it into one or more Cypher queries, and returns the results formatted as LLM-friendly text.

Benefits of this modular pattern:

Context Management: The sub-agent has its own context, keeping messy graph query details out of the main agent's window.

Specialization: The sub-agent prompt can include targeted instructions and examples of useful Cypher queries specific to our schema.

Maintainability: The graph query module is self-contained and can be evaluated or improved separately.

Our Technical Stack
Biolink: An ontology for biology-related knowledge graphs.

BioCypher: A framework for building biomedical knowledge graphs.

PyPath: A Python module for processing molecular biology data.

Neo4j: The graph database used to host and query the data.

LangChain/LangGraph: The LLM orchestration framework.

In the next post, we will walk through the first steps of building a Knowledge Graph, including Python code examples.

