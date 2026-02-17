# The Omics Deluge: Why I’m Training an AI to Read 30 Papers So I Don’t Have To

# The Problem: Data is Cheap, Insight is Expensive

Omics has essentially turned biology into a massive data science project. Thirty years ago, sequencing the human genome cost billions. Today? I can plug a sequencer the size of a stapler into my laptop and sequence seawater in my kitchen if I want to.

We’ve solved the "how to get data" problem. But we’ve created a massive "what does this mean" bottleneck.

When I run a sequencing experiment on Prochlorococcus—a tiny microbe that basically breathes for the planet by producing 10% of our oxygen—I don’t get a simple answer. I get a spreadsheet. Usually, it’s a list of about 3,000 upregulated genes.

# The 3,000 Gene Headache

The standard move here is to group genes into "pathways." If you see 50 genes for nitrogen transport, you assume the cell is hungry for nitrogen.

But biology is rarely that polite.

Many genes are moonlighters; they do three different things depending on the day. When you're looking at 3,000 active genes, how do you know what the cell is actually prioritizing? In my recent paper in Nature Microbiology, we modeled four specific ways these microbes might be interacting—things like "Exoenzyme" recycling or "ROS" detoxification.

The math gave us a likely winner, but math is just a simulation. To prove it, I need to see if the RNA and protein levels in the lab actually back up the theory.

# Enter Spock

This is where it gets frustrating. As a scientist, I can’t manually cross-reference 3,000 genes against every paper ever written on marine microbes. I’d be doing nothing else for a decade.

So, I’m building Spock.

It’s not just a chatbot. It’s a system designed to be a "Logic Officer." Right now, we’re training it on a curated set of about 30 high-quality studies—the core knowledge of our specific niche. The goal isn't to have the AI write my paper for me. It’s to have the AI "read" my raw data and tell me: "Based on these 30 studies, your current gene expression looks like a 90% match for the Exoenzyme model you predicted."

# What’s Next?

In the next post, I’ll stop talking about the "why" and show you the "how." I’ll be breaking down why a standard AI isn’t enough and why we’re building a Knowledge Graph to give Spock a brain that actually understands biological context.