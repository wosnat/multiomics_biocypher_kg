This file is README.txt for Supplemental Material relevant to GSEA (Gene Set Enrichment Analysis) 
for 

"Global gene expression of Prochlorococcus ecotypes in response to changes in nitrogen availability"

by Andrew C. Tolonen, John Aach, Debbie Lindell, Zackary I. Johnson, Trent Rector, Robert Steen, 
George M. Church, Sallie W. Chisholm


For general information on GSEA, see Subramanian et al. (2005) PNAS 102:15545


For questions on this material, please contact John Aach at aach@receptor.med.harvard.edu


The following are descriptions of the files provided in the GSEA directory:

GSEA_goldenspike_15_15.ppt

GSEA_goldenspike_20_20.ppt

GSEA_goldenspike_ortho_15.ppt

GSEA_for_NtcA_on_repressed_genes.ppt

------------------------------------

These are PowerPoint files containing graphs relevant to the use of GSEA to detect overrepresentation 
of genes with NtcA motifs in their promoters among significantly upregulated and downregulated genes 
at time points in the Nitrogen starvation experiments described in the article text.  

The three files GSEA_goldenspike_15_15.ppt, GSEA_goldenspike_20_20.ppt, and GSEA_goldenspike_ortho_15.ppt
look for overrepresentation of genes with NtcA motifs in their promotors among significantly
upregulated genes.  The three files therefore present evidence relevant to whether NtcA may be 
acting as an activator under N starvation.  The three files differ in the stringency and selection 
criteria used to identify genes with NtcA sites. 

We used the NtcA site score rankings of Su et al. (2005) Nucleic Acids Research 33: 5156 to 
identify genes with NtcA sites in their promoters. As explained in our article (see Methods), 
the method of Su et al. (2005) gives bonus points to genes whose NtcA sites are apparently 
conserved across other cyanobacterial genomes, and therefore tends to penalize sites for genes specific 
to any one species.  We therefore considered both ranks of NtcA sites among genes with orthologs,
and without orthologs, separately in identifying sets of genes to be used with GSEA.

GSEA_goldenspike_15_15.ppt defined 'genes with NtcA sites' as those whose NtcA site score ranks 
were among the top 15 for genes with orthologs, and also among the top 15 for genes without 
orthologs.  This is the way we identified genes with NtcA sites in our article.

GSEA_goldenspike_20_20.ppt and GSEA_goldenspike_ortho_15.ppt were alternative analyses that used 
different definitions of 'genes with NtcA sites'.  Although we did not feature them in our 
article, we took them into account there in correcting for multiple hypotheses in our Table 2 
and provide them here for completeness.  The first (_20_20) defined 'genes with NtcA sites' as 
those whose NtcA site score ranks were among the top 20 for genes with orthologs, and also among 
the top 20 for genes without orthologs.  The second (_ortho_15) defined 'genes with NtcA sites' as 
only those genes with orthologs whose NtcA site score ranks were among the top 15 for genes with 
orthologs.

In the text of our manuscript, we also discuss reports of NtcA functioning like a repressor
vs. an activator for some genes in some circumstances.  We looked for evidence of this function
by using GSEA to look for overrepresentation of genes with NtcA motifs in their promotors among
significantly downregulated genes. The file GSEA_for_NtcA_on_repressed_genes.ppt describes several
classes of genes for which we conducted GSEA tests that are distinguished by orthology and by 
the relationship of putative NtcA sites to putative -10 boxes identified by Su et.al. (2005), along
with the GSEA results.

Each PowerPoint file presents two graphs for each GSEA test.  The first is a running Enrichment 
Score (ES) graph of the sort described in Subramanian et al. (2005) that displays graphically 
the ES at each position in the complete set of genes presented to GSEA (which is ordered by degree 
of significant upregulation), and which identifies by a dashed vertical line the position with the 
maximum ES, which is used to define the GSEA statistic.  The second graph is a histogram of GSEA 
statistics based on the same gene expression data in the same order but with 5000 random reassignments 
of the NtcA sites among the genes, and which is used to estimate the P value of the GSEA score for 
the actual data.  The title of each graph identifies the Prochlorococcus species (MED4 or MIT9313) 
and time point (T0, T3, T6, T12, T24, T48) in the Nitrogen Starvation series described in our article.  
Additional information presented with the graphs includes the maximum ES value (ES), the P value 
based on the 5000 random reassignments of NtcA sites (P), a Chebyschev upper bound approximation 
of this P value based on the variance of the GSEA scores for the 5000 random reassignments (Pc; this 
estimate proved not to be useful), and the size of the GSEA 'leading edge' (le).  The P value (P)
gives the probability that the degree of enrichment of NtcA sites among the genes that is represented 
by the GSEA score could arise by chance alone. Among GSEA tests looking for evidence of NtcA 
function as an activator, the _15_15 and _20_20 analyses were performed for all time points for 
both Prochlorococcus species, while the _ortho_15 analysis was performed for only a subset of time 
points.  For the GSEA tests lookinf for evidence of NtcA repressor function, only the 6hr and 
48hr time points were considered.


gsea_pvalues.xls

----------------

This Excel spreadsheet contains the P values for enrichment of genes with NtcA sites among genes 
ordered by significant upregulation for the complete set of analyses described above.  It is a 
superset of the results presented in Table 2 of our article.  The results of GSEA tests for 
NtcA repressor function are not included in this spreadsheet.


NtcA_gsea_tests.tar.gz

----------------------

This file is a compressed archive that contains the MatLab script used to perform the GSEA analyses,
and all of the numerical data generated by all of the analyses mentioned above.  To access the
individual GSEA test information, one must decompress the archive with gunzip and then untar it 
with tar -xvf.  The individual files are:

gsea2.m: The MatLab function used to perform the analysis. Comments at the beginning of this function 
  describe parameters and their usage, and are also available by typing "help gsea2" from the MatLab 
  command line.

Su_MED4_ranks.txt and Su_MIT9313_ranks.txt: Extracts of supplemental material from Su et al. (2005)
  that are augmented with ortholog and non-ortholog rank information, and which were used for input to 
  gsea2.m to analyze NtcA site enrichment.

gsea_pvalues.txt : A text file version of gsea_pvalues.xls (described above).


files of the form gsea_(species)_(timepoint)_(NtcA_set).txt (30 files): These files contain the output
  from gsea2.m for each of the GSEA analyses for NtcA inducer function.  Each file corresponds to
  a pair of GSEA test figures in one of GSEA_goldenspike_15_15.ppt, GSEA_goldenspike_20_20.ppt, and 

  GSEA_goldenspike_ortho_15.ppt.  

  The conventions encoded in the file names are illustrated as follows: gsea_MED4_T12_15_15.txt 
  is the output for the analysis of time point T12 for MED4 using the _15_15 definition of NtcA sites 
  described above.  

  Each of these files summarizes the parameters passed to gsea2 to perform the analysis, the 
  statistics generated for the analysis (most of which are presented in the PowerPoints), and the 
  running list of ES scores generated for genes and identified NtcA sites that are computed as part 
  of the analysis.

files of the form gsea_(species)_(timepoint)_repress_(NtcA_set).txt (10 files): These files contain
  the output from gsea2.m for each of the GSEA analyses for NtcA repressor function.  The NtcA_set
  keywords are either "_15_15", "_15_15_no10box", and "_10pos_lt_15", which mean, respectively:

  _15_15: these are identical to the NtcA activator GSEA tests except that they look for enrichment
    among repressed vs. activator genes.  These tests were negative controls.

  _15_15_no10box: These consider the top 15 NtcA sites identified by Su et.al. (2005) for genes that 
    did 
not have a putative -10 box downstream of the putative NtcA site, where genes with orthologs
    were ranked separately from genes without orthologs.

  _10pos_lt_15: These consider all NtcA sites identified by Su et.al. (2005) for genes that had 
    a putative NtcA site that was less than 15 bp away from (upstream of) a putative -10 box.

  Further details on these options are given in GSEA_for_NtcA_on_repressed_genes.ppt.

