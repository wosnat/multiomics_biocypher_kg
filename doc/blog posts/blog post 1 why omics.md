
Title: Stop Doing Simplistic Omics. Start Using AI.
Subtitle: We have decades of data at our fingertips. It’s time we used AI to understand the cell as a whole, not just as a list of upregulated genes.

Why use AI to research omics?

We use AI for everything these days, right? So why not use it to research omics? Let me share some of my experience. I have just finished my PhD in biology, where I had to analyze omics results (transcriptomics and proteomics), and I was a little bit frustrated by how simplistic these analysis methods are. We don't really know everything that's going on in the cell, so we end up doing basic statistical tests. For example, we look for enriched pathways, or we pick a few differentially expressed (DE) genes and from there we come up with a hypothesis of what is actually going on.

It's really hard to understand exactly what is going on in cellular biology because the cell machinery is so complex. However, in the last few decades, there have been so many studies of transcriptomics, omics, and genetics in general. The wealth of data at our fingertips is immense. Unfortunately, it's beyond human capacity to assimilate all of that into coherent research. Dealing with large datasets is one of AI's strong suits, making it ideal for this purpose.

AI in biomedical research

Several AI models have been developed in the fields of genetics and proteomics. The most famous of these is AlphaFold, which uses AI to predict the structure of protein folding. Protein Language Models (PLM) use a list of amino acids as a vocabulary, and genomic Language Models (gLM) use a list of nucleotides (A, G, C, T). These models take existing genes and proteins and generate an encoding—which you can think of as a kind of dimensionality reduction.

This encoding can be used for downstream applications. For example, you can use it for clustering; if a novel protein is close to an existing encoding, then that gives us a hint of what that protein's function might be. Other applications include identifying characteristics, such as signal peptides, motifs, antimicrobial peptides, and more.

In the field of omics analysis, there have been several models proposed for studying scRNA-seq (single-cell RNA sequencing). This analysis returns transcription data for hundreds of thousands of cells. An AI model encoding this transcription pattern can be used to cluster cells of the same tissues together and identify their specific activation patterns. This is also promising.

AI for omics analysis

All these models aim to decipher the language of the cell by detecting hidden patterns. My goal, however, is slightly different: I want to use AI to analyze current omics results within the context of the vast ocean of published data and public genomic databases.

The truth is, even with AI, much of our current work remains simplistic. We aren't fully leveraging the hundreds of sequencing studies published in recent years. I started asking myself: what if we took all those existing studies, fed them to an AI model, and asked it to help us understand a new study in the context of everything we already know? By moving beyond basic enrichment tests and statistical correlations, we can build a "knowledge center" that provides a truly holistic approach. That is where we will find real biological insights.

In the next post, I'll talk about how to actually do this. Stay tuned.
