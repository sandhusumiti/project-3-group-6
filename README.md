# project-3-group-6


## biologist.Rmd 
Creates a heatmap based on the normalized counts matrix, filters genes for enrichment analysis

## run_limma.R
This file uses microarray results and the tox-group information as input to perform a differential expression analysis using limma. The output is three CSV files that contain information about the differential expression of each probe.

## make_histograms_significant_genes.R
This file uses the three CSV files from the run_limma.R script and a CSV file of the probes mapping to genes as input. This code generates histograms of the log fold change values of the significant genes for each chemical. This code also generates volcano plots of the significant genes for each chemical.

## part_6.R
This file uses the same input files as make_histograms_significant_genes.R and the DESeq files for each chemical. This code maps the probes in the limma output to the refseqids in the DESeq output and calculates the concordance between the two files. This concordance is then plotted as a function of the number of differentially expressed genes in each analysis. The concordance values for above- and below-median expression genes was also calculated. These results were plotted on a bar graph that displayed the overall, above- and below-median concordance values for each chemical. 
