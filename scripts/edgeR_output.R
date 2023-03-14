#!! INCLUDE SUBSETTING!!!

#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		19-APR-2013
### Modified:		19-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description:	Generates a number of output tables and plots of each individual (single or pairwise) edgeR analysis conducted with edgeR_dge_disp_splice.R in which all relevant data objects are stored in a list 'objects'; output files are written to the output folders defined in edgeR_dge_disp_splice.R 	  
### Arguments: 		1. Path containing the saved objects resulting from edgeR_dge_disp_splice.R [STRING | CHARACTER]; 2. Glob-style pattern for file selection [STRING | CHARACTER] 
### Output: 		MDS-plot, heatmap
### Usage example:	Rscript edgeR_output.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- args[1]
pattern <- args[2]
suppressMessages(library(edgeR))
###

### B. Import data
# Import file table from input directory
files <- dir(path=in_folder, pattern=glob2rx(pattern), full.names=TRUE, recursive=TRUE)
## Traverse file list row by row; returns list of ls 
ls <- lapply(files, function(filename) {
	# Load object 'objects' saved in file 'filename'
	load(filename)		

	### C. Write tables
	## Print tables of log2 fold changes, counts per million, P values and false discovery rates
	# All exons
	if (exists("objects$all_exon", inherits=FALSE)) write.table(objects$all_exon, file=paste(objects$prefix, "fc_cpm_p_fdr_all_exons.tab", sep="_"), quote=FALSE, sep="\t")
	# Differentially expressed exons
	if (exists("objects$diff_exon", inherits=FALSE)) write.table(objects$diff_exon, file=paste(objects$prefix, "fc_cpm_p_fdr_diff_expr_exons.tab", sep="_"), quote=FALSE, sep="\t")
	# All genes
	if (exists("objects$all_gene", inherits=FALSE)) write.table(objects$all_gene, file=paste(objects$prefix, "fc_cpm_p_fdr_all_genes.tab", sep="_"), quote=FALSE, sep="\t")
	# Differentially expressed genes
	if (exists("objects$diff_gene", inherits=FALSE)) write.table(objects$diff_gene, file=paste(objects$prefix, "fc_cpm_p_fdr_diff_expr_genes.tab", sep="_"), quote=FALSE, sep="\t")
	## Alternative splicing analysis
	# Print table of log2 fold changes of each exon relative to the fold change of the corresponding gene (i.e. all exons)
	if (exists("objects$all_exon_spl", inherits=FALSE)) write.table(objects$all_exon_spl, file=paste(objects$prefix, "splice_anal_rel_fc_all_exons.tab", sep="_"), quote=FALSE, sep="\t")
	# Print P values and false discoveries for evidence suggesting differential splicing -> all genes
	if (exists("objects$all_gene_spl", inherits=FALSE)) write.table(objects$all_gene_spl, file=paste(objects$prefix, "splice_anal_cpm_p_fdr_all_genes.tab", sep="_"), quote=FALSE, sep="\t")
	# Print P values and false discoveries for evidence suggesting differential splicing -> differentially spliced genes
	if (exists("objects$diff_gene_spl", inherits=FALSE)) write.table(objects$diff_gene_spl, file=paste(objects$prefix, "splice_anal_cpm_p_fdr_diff_spliced_genes.tab", sep="_"), quote=FALSE, sep="\t")
	###
	
	### D. Generate plots
	
})

					### E. Clean-up & save session image
					# Remove unused/temp variables
					rm(df, df_files, input_df, ls, prefix)
					# Save workspace image/session
					save.image(file=paste(prefix, "_image.Rdata", sep=""))
					###
