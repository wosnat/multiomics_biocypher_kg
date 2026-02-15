# Title: Building an AI for Omics
# Subtitle: why Fold Change is the "Sweet Spot" for Training Omics AI


# tags:
Bioinformatics (Very active community)

Artificial Intelligence (The broadest reach)

Data Science (For those interested in data "sweet spots")

Genomics (To capture the omics/biology crowd)

Biotech (To reach the industry/startup audience)




# Toward an AI That Understands Omics

We want to create an AI system that can understand omics and answer complex questions about the data. But how do we actually go about it? The first question to answer is: what level of detail should we use in the integration?

# The Flow of RNA-seq Analysis

As a driving example, let’s look at RNA-seq analysis. We start with an experiment involving treatments and controls. RNA is extracted and sequenced, resulting in FASTQ files. These files contain millions of short reads from the extracted mRNA.

The next step is to map those reads to a reference genome, which basically gives us a list of reads per gene. We can then count how many reads are mapped to each specific gene. However, this isn't enough because every sample and every sequencing run is slightly different. We need to normalize the data, that's the next step. We normalize for the length of the gene and the total number of reads in the sample, giving us FPKM (Fragments Per Kilobase of transcript per Million mapped reads).

But even that isn't enough; a number in isolation is hard to interpret. The next step is to compare the treatment to the control. How different is this number from the control? This is the Fold Change (FC). Ideally, you sequence multiple replicates to gain statistical power. Finally, we can identify "differentially expressed genes", those with a statistically valid fold change. The final step is to look into which pathways these genes belong to and what they actually do.

There are other things you can do besides counting genes. You can look for additional features; for example, some reads may differ from your reference genome. This mismatch allows you to detect mutations in your sequenced cells, such as SNPs (Single Nucleotide Polymorphisms). You can also look for antisense reads, which are reads in the opposite direction of the gene. These are often regulatory mechanisms; an antisense RNA can adhere to the mRNA transcript and stop it from being translated. There is a wealth of knowledge beyond just fold change in these sequencing results.

# The Details in the publications

In final publications, you usually only get the "bottom line." The authors describe the experiment, analyze the results, and present a hypothesis. In a good publication, the authors select a very small subset of the differentially expressed genes or pathways to build a physiological narrative.

You cannot realistically go over every low-level result in a paper; it’s not human-readable, and frankly, it would be a bit boring. If we only look at the publication text, we only see a tiny subset of the study's findings. This isn't necessarily "cherry-picking" - it’s because listing everything as-is has low scientific value if you can't make sense of it. However, when integrating this into an AI knowledge system, relying only on the genes mentioned in the text leaves out a massive amount of useful data.

# The Level of Detail for the AI Engine

When deciding what level of detail to feed into the AI model, we considered a few options:

Publications: Too high-level. There is too much data loss.

Normalized Read Counts: Too low-level. The AI gets lost in the details, as this data is hard to interpret without specific analytical tools.

Fold Change Statistics: The "Sweet Spot".

We selected the level of fold change: which genes are up- or down-regulated, by how much, under which conditions, and at what p-value. A single publication might have 10–20 different fold change statistics (e.g., across different timepoints and different treatments). These results are often tucked away in supplementary tables that can be used almost as-is.

Overall, fold change is the ideal level. Unlike raw numbers, it has a clear biological interpretation and a P-value to evaluate validity. It captures all differentially expressed genes, even those the authors didn't discuss, ensuring important data isn't lost.

Fold change statistics are ubiquitous and a good starting point. However, some papers offer additional analysis, like periodicity, transcriptomic dynamics, or co-expression.
In the world of omics, there are additional assays, such as proteomics and metabolomics. Each needs to be added to the knowledge base in a way that allows for holistic reasoning, which is beyond the scope of this post.

Currently, our multi-omics knowledge base contains 110,000 fold change statistics across 50 different treatments, uploaded from 19 publications. This provides a powerful starting point for our AI.

# Coming Next: Connecting the Dots
Now that we’ve found the "sweet spot" for our data, the next big challenge is figuring out how to organize it so an AI can actually reason across it. 110,000 data points are great, but they don't mean much if they’re just sitting in a pile. In Post #3, I’ll dive into why we’re using a Knowledge Graph to map these complex biological relationships and turn them into a system that truly "thinks".

# About the Project
We are building a next-generation Multi-Omics Knowledge Base designed for holistic biological reasoning. By integrating high-fidelity data like fold-change statistics from thousands of studies, we are creating an AI system that doesn't just read papers - it understands the underlying science. Currently, our engine hosts over 110,000 data points from 50+ treatments, serving as a scalable foundation for the future of automated biological discovery.

Enjoyed this deep dive? This is the second post in my series on building a brain for omics. Follow me here on Medium to stay updated as we scale our knowledge base and tackle the next challenge: holistic reasoning across different sequencing assays.


Pro-Tip for Medium:
At the bottom of your post, you can also add a link to your first post. Something like:

Missed the beginning? Read [Post #1 Title] here.


Pro-Tip for your Medium Series:
Since this is a multi-part project, I highly recommend creating a "Medium Series" or a "List" titled something like "Building an AI for Omics." You can then link to that list in every post, making it easy for new readers to binge-read from the beginning!
