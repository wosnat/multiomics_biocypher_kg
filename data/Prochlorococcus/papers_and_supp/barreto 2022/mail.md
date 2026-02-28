---------- Forwarded message ---------
From: Malisano Barreto Filho, Marcelo <barretof@uab.edu>
Date: Thu, 26 Feb 2026 at 19:24
Subject: Re: Translation tables for cyanobacterial transcriptomes
To: Daniel Sher <dsher@univ.haifa.ac.il>
Cc: Morris, James Jeffrey <evolve@uab.edu>, Osnat Weissberg <osnat.weissberg@gmail.com>


Hi Daniel and Osnat,It’s great to hear from you! Thanks for your message and for sharing the details of your project. It sounds like a really exciting approach.I’m happy to help. I’ve attached the data spreadsheets, which contain the translation information for the gene symbols (where available), gene products, UniProt IDs, EC numbers, and other details you’re looking for. These are quite valuable; we originally obtained them from the CYORF (Cyanobacteria Open Reading Frame) database, which compiled information on cyanobacterial genes and genomes with a focus on community-driven annotation. CYORF was developed through a collaboration between the Kazusa DNA Research Institute, the Cyanobacteria DNA Chip Consortium, and Kyoto University (GenomeNet, which maintains the KEGG database). The information in the spreadsheets is based on the KEGG IDs for MIT9312, CC9311, and WH8102. Unfortunately, the CYORF project ended in 2022, possibly due to funding issues (http://cyano.genome.jp/).We also created annotation packages that you can use in your code to add any information to your count table using our pipeline. You can find them in our OSF repository, along with all our code and data using the link below

https://osf.io/dzfht/overview?view_only=b59419aedcf040d7a9afff0cc072e00b 

.For Alteromonas, we did not have a specific kegg organismal id, so we had to use a different approach using ko numbers (from KAAS) which allows us to conduct our pathway analysis. I am also attaching the annotation.gtf file we used for feature counts. I think that featurecounts require .gtf to work, can't remember now. But that is how you can obtain the ids we used, and they are easily interconvertable.If you need any additional files or clarification, just let me know.I’m looking forward to seeing how your integration progresses. Osnat’s blog sounds great, so please do share it when it’s ready. Also, over the past year, I’ve gained experience in metabolic reconstruction, curation, FBA/FVA, and community metabolic models. If you’re interested in collaborating on any related projects in the future, please let me know.Best wishes,
Marcelo