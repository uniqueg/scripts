#######
### GENERAL:
### --------
### Author: 		Christina Herrmann
### Created: 		30-APR-2013
### Modified:		30-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	From objects created by edgeR_???? script,a MDS plot is generated for all samples, regardless of treatment groups. Finally, a heatmap of tagwise dispersion values is generated (x = groups, y = features)  
### Arguments: 		1. Tab-separated file containing no header and the following column data: I. PATH to or FILE NAME of count table(s) of format  [FEATURE | STRING | must be unique!] [TAB] [COUNT | INTEGER] with header line (output of find_overlaps_bam_within_bed.R and derived), II. Group name, III. (only required if I. does not indicate absolute file path!) Glob-style pattern indicating which files to select from indicated PATH in I. 
### Output: 		MDS-plot, heatmap
### Usage example:	Rscript edgeR_MDS_all_NINA.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- args[1]
pattern <- args[2]
out_prefix <- args[3]
suppressMessages(library(edgeR))
###

### B. Import data
# Import file table from input directory
files <- dir(path=in_folder, pattern=glob2rx(pattern), full.names=TRUE, recursive=TRUE)

# remove all pairwise comparisons files (those containing "vs" in file name)
# !!!!!!!!!!!!! removes all files if no file contains "vs" in file name
#vs <- grep("*_vs_*",files)
#files <- files[-vs]

# vector for group names of loaded objects, used for naming list elements later
groups <- NULL

## Traverse file list row by row; returns list pseudo_ls, each entry=1 tissue, contains list of matrices for plots
pseudo_ls <- lapply(files, function(filename) {
			# Print status message
			cat("Processing file '", filename , "'...\n", sep="")
			# load file for current tissue
			load(filename)
			# Retrieve group of current object 
			group <- ifelse(length(objects$row)>4,paste(objects$row[1], objects$row[5], sep="_vs_"),objects$row[1])
			
			# genes: extract relevant data from object for current tissue
			gene <- lapply(objects$DATA, function (set) {
						mt <- set$gene$dge$pseudo.counts
						dimnames(mt) <- list(rownames(mt),rep(group,dim(mt)[2]))
						return(mt)
					})
			# exons: extract relevant data from object for current tissue
			exon <- lapply(objects$DATA, function (set) {
						mt <- set$exon$dge$pseudo.counts
						dimnames(mt) <- list(rownames(mt),rep(group,dim(mt)[2]))
						return(mt)
					})
			groups <<- c(groups, group)
			# Print status message
			cat("Processed file '", filename , "'!\n", sep="")
			return(list(gene=gene,exon=exon))
		})

# Set names of pseudo_ls entries to corresponding group
names(pseudo_ls) <- groups

# Combine pseudo.counts for all groups for all gene subsets into 
# one matrix for each gene subset
#all_gene <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$gene$all)))
ref_all_gene <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$gene$refseq_all)))
#ref_spl_gene <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$gene$refseq_splice)))
#ref_3_gene <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$gene$refseq_3._end)))

#all_exon <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$exon$all)))
ref_all_exon <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$exon$refseq_all)))
ref_spl_exon <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$exon$refseq_splice)))
ref_3_exon <- na.omit(do.call(cbind, lapply(pseudo_ls, function(group) group$exon$refseq_3._end)))



### C. MDS PLOT ALL GROUPS

## Preparations
# Set colours to be used		
palette(c("yellow2","orchid","orange","tomato","tomato4","pink","peachpuff","turquoise","plum","steelblue1","royalblue4","lightblue","palegreen","brown","deeppink2","red","tomato2","yellow4","seagreen","springgreen2","darkgreen","wheat4","grey","black","lightgrey"))

## MDS Plots

# MDS: function for plotMDS, input: mt:matrices with pseudo.counts for all groups, out:output prefix
MDS <- function(mt,out){
	# Construct output file name
	out_name <- paste(out_prefix, out, "_MDS.pdf", sep="")
	
	# Save Data matrix to table
	write.table(mt, file=paste(out_prefix, out, "_pseudo.counts_matrix.tab", sep=""), quote=FALSE, sep="\t")
	
	pdf(out_name, height=6, width=8)
	
	par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE)
	plotMDS(mt, gene.selection="common",labels=NULL, cex=0.8, col=as.numeric(as.factor(colnames(mt))) , main=paste("Multidimensional Scaling Plot\n",out))
	
	legend("topright", inset=c(-0.45,0), legend=levels(as.factor(colnames(mt))), fill=1:length(levels(as.factor(colnames(mt)))))
	
	dev.off()
}

#plotMDS genes
MDS(ref_all_gene, "Refseq_all_Gene")
#MDS(ref_spl_gene, "Refseq_splice_Gene")
#MDS(ref_3_gene, "Refseq_3UTR_Gene")
#plotMDS exons
MDS(ref_all_exon, "Refseq_all_Exon")
MDS(ref_spl_exon, "Refseq_splice_Exon")
MDS(ref_3_exon, "Refseq_3UTR_Exon")


