#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		19-APR-2013
### Modified:		19-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	From objects created by edgeR_???? script,a MDS plot is generated for all samples, regardless of treatment groups. Finally, a heatmap of tagwise dispersion values is generated (x = groups, y = features)  
### Arguments: 		1. Tab-separated file containing no header and the following column data: I. PATH to or FILE NAME of count table(s) of format  [FEATURE | STRING | must be unique!] [TAB] [COUNT | INTEGER] with header line (output of find_overlaps_bam_within_bed.R and derived), II. Group name, III. (only required if I. does not indicate absolute file path!) Glob-style pattern indicating which files to select from indicated PATH in I. 
### Output: 		Heatmap
### Usage example:	Rscript edgeR_heatmap.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- "." # args[1]
pattern <- "*objects.Rdata.2"  # args[2]
out_prefix <- "./TEST_" # args[3]
suppressMessages(library(edgeR))
###

### B. Import/process data
# Import file table from input directory
files <- dir(path=in_folder, pattern=glob2rx(pattern), full.names=TRUE, recursive=TRUE)
# Remove all non-pairwise comparisons files (those NOT! containing "vs" in file name)
files <- files[grep("*_vs_*",files)]
# Vector for group names of loaded objects, used for naming list elements later
groups <- NULL
# Vector for FDR cutoff, used later
fdr <- NULL
## Traverse file list row by row; returns list 'fc_ls'  containing LOG2 FOLD CHANGES
fc_ls <- lapply(files, function(filename) {
	# Print status message
	cat("Processing file '", filename , "'...\n", sep="")
	## Load file for current tissue
	load(filename)
	# Extract group/comparison name
	group <- ifelse(length(objects$row)>4,paste(objects$row[1], objects$row[5], sep="_vs_"),objects$row[1])
	## GENES: Extract fold changes vector as matrix
	gene <- lapply(objects$DATA, function (set) {
		matrix(set$gene$tt$table[,1], dimnames=list(rownames(set$gene$tt$table), group))	
	})
	## EXONS: Extract fold changes vector as matrix
	exon <- lapply(objects$DATA, function (set) {
		matrix(set$exon$tt$table[,1], dimnames=list(rownames(set$exon$tt$table), group))	
	})
	## SPLICE: Extract fold changes vector as matrix
	splice <- lapply(objects$DATA, function (set) {
		matrix(set$splice$exon_df[,1], dimnames=list(rownames(set$splice$exon_df), group))
	})
	# Add current group to global character vector 'groups'
	groups <<- c(groups, group)
	# Extract FDR cutoff
	if (is.null(fdr)) fdr <<- objects$p_cutoff
	# Print status message
	cat("Processed file '", filename , "'!\n", sep="")
	# Return list of data
	return(list(gene=gene, exon=exon, splice=splice))
})
## Traverse file list row by row; returns list 'fdr_ls' containing FALSE DISCOVERY RATES
fdr_ls <- lapply(files, function(filename) {
	# Print status message
	cat("Processing file '", filename , "'...\n", sep="")
	## Load file for current tissue
	load(filename)
	# Extract group/comparison name
	group <- ifelse(length(objects$row)>4,paste(objects$row[1], objects$row[5], sep="_vs_"),objects$row[1])
	## GENES: Extract fold changes vector as matrix
	gene <- lapply(objects$DATA, function (set) {
		matrix(set$gene$tt$table[,4], dimnames=list(rownames(set$gene$tt$table), group))	
	})
	## EXONS: Extract fold changes vector as matrix
	exon <- lapply(objects$DATA, function (set) {
		matrix(set$exon$tt$table[,4], dimnames=list(rownames(set$exon$tt$table), group))	
	})
	## SPLICE: Extract fold changes vector as matrix
	splice <- lapply(objects$DATA, function (set) {
		a_ids <- sub('\\.\\d+$', "", rownames(set$splice$exon_df))
		fdr <- set$splice$gene_tt$table[match(a_ids, set$splice$gene_tt$table[,1]),5]
		matrix(fdr, dimnames=list(rownames(set$splice$exon_df), group))
	})
	# Print status message
	cat("Processed file '", filename , "'!\n", sep="")
	# Return list of data
	return(list(gene=gene, exon=exon, splice=splice))
})
## Rename list names with entries in vector 'groups'
names(fc_ls) <- groups
names(fdr_ls) <- groups
###

### C. Combine/prepare data
## FOLD CHANGES
all_gene_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$gene$all)))
ref_all_gene_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$gene$refseq_all)))
ref_spl_gene_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$gene$refseq_splice)))
ref_3_gene_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$gene$refseq_3._end)))

all_exon_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$exon$all)))
ref_all_exon_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$exon$refseq_all)))
ref_spl_exon_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$exon$refseq_splice)))
ref_3_exon_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$exon$refseq_3._end)))

all_splice_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$splice$all)))
ref_all_splice_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$splice$refseq_all)))
ref_spl_splice_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$splice$refseq_splice)))
ref_3_splice_fc <- na.omit(do.call(cbind, lapply(fc_ls, function(group) group$splice$refseq_3._end)))

## FALSE DISCOVERY RATES
all_gene_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$gene$all)))
ref_all_gene_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$gene$refseq_all)))
ref_spl_gene_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$gene$refseq_splice)))
ref_3_gene_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$gene$refseq_3._end)))

all_exon_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$exon$all)))
ref_all_exon_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$exon$refseq_all)))
ref_spl_exon_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$exon$refseq_splice)))
ref_3_exon_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$exon$refseq_3._end)))

all_splice_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$splice$all)))
ref_all_splice_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$splice$refseq_all)))
ref_spl_splice_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$splice$refseq_splice)))
ref_3_splice_fdr <- na.omit(do.call(cbind, lapply(fdr_ls, function(group) group$splice$refseq_3._end)))

## FOLD CHANGES FDR < CUTOFF
all_gene_fc_cut_ls <- sapply(colnames(all_gene_fc), function(col_name) all_gene_fc[all_gene_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_all_gene_fc_cut_ls <- sapply(colnames(ref_all_gene_fc), function(col_name) ref_all_gene_fc[ref_all_gene_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_spl_gene_fc_cut_ls <- sapply(colnames(ref_spl_gene_fc), function(col_name) ref_spl_gene_fc[ref_spl_gene_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_3_gene_fc_cut_ls <- sapply(colnames(ref_3_gene_fc), function(col_name) ref_3_gene_fc[ref_3_gene_fdr[ ,col_name] < fdr, col_name, drop=FALSE])

all_exon_fc_cut_ls <- sapply(colnames(all_exon_fc), function(col_name) all_exon_fc[all_exon_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_all_exon_fc_cut_ls <- sapply(colnames(ref_all_exon_fc), function(col_name) ref_all_exon_fc[ref_all_exon_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_spl_exon_fc_cut_ls <- sapply(colnames(ref_spl_exon_fc), function(col_name) ref_spl_exon_fc[ref_spl_exon_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_3_exon_fc_cut_ls <- sapply(colnames(ref_3_exon_fc), function(col_name) ref_3_exon_fc[ref_3_exon_fdr[ ,col_name] < fdr, col_name, drop=FALSE])

all_splice_fc_cut_ls <- sapply(colnames(all_splice_fc), function(col_name) all_splice_fc[all_splice_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_all_splice_fc_cut_ls <- sapply(colnames(ref_all_splice_fc), function(col_name) ref_all_splice_fc[ref_all_splice_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_spl_splice_fc_cut_ls <- sapply(colnames(ref_spl_splice_fc), function(col_name) ref_spl_splice_fc[ref_spl_splice_fdr[ ,col_name] < fdr, col_name, drop=FALSE])
ref_3_splice_fc_cut_ls <- sapply(colnames(ref_3_splice_fc), function(col_name) ref_3_splice_fc[ref_3_splice_fdr[ ,col_name] < fdr, col_name, drop=FALSE])

## DENSITIES
all_gene_den_ls <- apply(all_gene_fc, 2, function(col) density(col, from=-10, to=10))
ref_all_gene_den_ls <- apply(ref_all_gene_fc, 2, function(col) density(col, from=-10, to=10))
ref_spl_gene_den_ls <- apply(ref_spl_gene_fc, 2, function(col) density(col, from=-10, to=10))
ref_3_gene_den_ls <- apply(ref_3_gene_fc, 2, function(col) density(col, from=-10, to=10))

all_exon_den_ls <- apply(all_exon_fc, 2, function(col) density(col, from=-10, to=10))
ref_all_exon_den_ls <- apply(ref_all_exon_fc, 2, function(col) density(col, from=-10, to=10))
ref_spl_exon_den_ls <- apply(ref_spl_exon_fc, 2, function(col) density(col, from=-10, to=10))
ref_3_exon_den_ls <- apply(ref_3_exon_fc, 2, function(col) density(col, from=-10, to=10))

all_splice_den_ls <- apply(all_splice_fc, 2, function(col) density(col, from=-10, to=10))
ref_all_splice_den_ls <- apply(ref_all_splice_fc, 2, function(col) density(col, from=-10, to=10))
ref_spl_splice_den_ls <- apply(ref_spl_splice_fc, 2, function(col) density(col, from=-10, to=10))
ref_3_splice_den_ls <- apply(ref_3_splice_fc, 2, function(col) density(col, from=-10, to=10))

## DENSITIES FDR < CUTOFF
all_gene_den_cut_ls <- lapply(all_gene_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_all_gene_den_cut_ls <- lapply(ref_all_gene_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_spl_gene_den_cut_ls <- lapply(ref_spl_gene_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_3_gene_den_cut_ls <- lapply(ref_3_gene_fc_cut_ls, function(elem) density(elem, from=-10, to=10))

all_exon_den_cut_ls <- lapply(all_exon_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_all_exon_den_cut_ls <- lapply(ref_all_exon_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_spl_exon_den_cut_ls <- lapply(ref_spl_exon_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_3_exon_den_cut_ls <- lapply(ref_3_exon_fc_cut_ls, function(elem) density(elem, from=-10, to=10))

all_splice_den_cut_ls <- lapply(all_splice_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_all_splice_den_cut_ls <- lapply(ref_all_splice_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_spl_splice_den_cut_ls <- lapply(ref_spl_splice_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
ref_3_splice_den_cut_ls <- lapply(ref_3_splice_fc_cut_ls, function(elem) density(elem, from=-10, to=10))
###

### D. Save session image
save.image(objects, file=paste(out_prefix, "logFC_dens_FDR_objects.Rdata", sep=""))
###


### Tables
ls_all_gene <- sapply(colnames(all_gene_fc), function(col) {
	diff <- all_gene_fc[all_gene_fdr[,col] < fdr,col]	
	len_all <- length(all_gene_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_all_gene, file="ls_all_gene.tab", quote=FALSE, sep="\t")

ls_ref_all_gene <- sapply(colnames(ref_all_gene_fc), function(col) {
	diff <- ref_all_gene_fc[ref_all_gene_fdr[,col] < fdr,col]	
	len_all <- length(ref_all_gene_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_all_gene, file="ls_ref_all_gene.tab", quote=FALSE, sep="\t")

ls_ref_spl_gene <- sapply(colnames(ref_spl_gene_fc), function(col) {
	diff <- ref_spl_gene_fc[ref_spl_gene_fdr[,col] < fdr,col]	
	len_all <- length(ref_spl_gene_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_spl_gene, file="ls_ref_spl_gene.tab", quote=FALSE, sep="\t")

ls_ref_3_gene <- sapply(colnames(ref_3_gene_fc), function(col) {
	diff <- ref_3_gene_fc[ref_3_gene_fdr[,col] < fdr,col]	
	len_all <- length(ref_3_gene_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_3_gene, file="ls_ref_3_gene.tab", quote=FALSE, sep="\t")

ls_all_exon <- sapply(colnames(all_exon_fc), function(col) {
	diff <- all_exon_fc[all_exon_fdr[,col] < fdr,col]	
	len_all <- length(all_exon_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_all_exon, file="ls_all_exon.tab", quote=FALSE, sep="\t")

ls_ref_all_exon <- sapply(colnames(ref_all_exon_fc), function(col) {
	diff <- ref_all_exon_fc[ref_all_exon_fdr[,col] < fdr,col]	
	len_all <- length(ref_all_exon_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_all_exon, file="ls_ref_all_exon.tab", quote=FALSE, sep="\t")

ls_ref_spl_exon <- sapply(colnames(ref_spl_exon_fc), function(col) {
	diff <- ref_spl_exon_fc[ref_spl_exon_fdr[,col] < fdr,col]	
	len_all <- length(ref_spl_exon_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_spl_exon, file="ls_ref_spl_exon.tab", quote=FALSE, sep="\t")

ls_ref_3_exon <- sapply(colnames(ref_3_exon_fc), function(col) {
	diff <- ref_3_exon_fc[ref_3_exon_fdr[,col] < fdr,col]	
	len_all <- length(ref_3_exon_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_3_exon, file="ls_ref_3_exon.tab", quote=FALSE, sep="\t")

ls_all_splice <- sapply(colnames(all_splice_fc), function(col) {
	diff <- all_splice_fc[all_splice_fdr[,col] < fdr,col]	
	len_all <- length(all_splice_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_all_splice, file="ls_all_splice.tab", quote=FALSE, sep="\t")

ls_ref_all_splice <- sapply(colnames(ref_all_splice_fc), function(col) {
	diff <- ref_all_splice_fc[ref_all_splice_fdr[,col] < fdr,col]	
	len_all <- length(ref_all_splice_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_all_splice, file="ls_ref_all_splice.tab", quote=FALSE, sep="\t")

ls_ref_spl_splice <- sapply(colnames(ref_spl_splice_fc), function(col) {
	diff <- ref_spl_splice_fc[ref_spl_splice_fdr[,col] < fdr,col]	
	len_all <- length(ref_spl_splice_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_spl_splice, file="ls_ref_spl_splice.tab", quote=FALSE, sep="\t")

ls_ref_3_splice <- sapply(colnames(ref_3_splice_fc), function(col) {
	diff <- ref_3_splice_fc[ref_3_splice_fdr[,col] < fdr,col]	
	len_all <- length(ref_3_splice_fc[,col])	
	len_diff <- length(diff)
	len_pc <- len_diff/len_all * 100
	len_gr1_up <- length(diff[diff < 0])
	len_gr2_up <- length(diff[diff > 0])
	skew <- len_gr1_up/len_gr2_up
	if (skew < 1) skew <- 1/skew
	return(data.frame(len_all=len_all, len_diff=len_diff, len_pc=len_pc, len_gr1_up=len_gr1_up, len_gr2_up=len_gr2_up, skew=skew))
}, USE.NAMES=TRUE)
write.table(ls_ref_3_splice, file="ls_ref_3_splice.tab", quote=FALSE, sep="\t")



ls_ref_all_gene_fc <- sapply(colnames(all_gene_fc), function(col) all_gene_fc[all_gene_fdr[,col] < fdr,col, drop=FALSE], USE.NAMES=TRUE)
ls_ref_spl_gene_fc <- sapply(colnames(all_gene_fc), function(col) all_gene_fc[all_gene_fdr[,col] < fdr,col, drop=FALSE], USE.NAMES=TRUE)
ls_ref_3_gene_fc <- sapply(colnames(all_gene_fc), function(col) all_gene_fc[all_gene_fdr[,col] < fdr,col, drop=FALSE], USE.NAMES=TRUE)


###################################### ALL FEATS
################ GENE
# ref splice factors
pdf("gene_splice.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_spl_gene_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_spl_gene_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_spl_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_spl_gene_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_spl_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_spl_gene_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref 3' end processing
pdf("gene_3.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_3_gene_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_3_gene_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_3_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_3_gene_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_3_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_3_gene_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref all genes
pdf("gene_all.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_all_gene_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_all_gene_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_all_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_all_gene_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_all_gene_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_all_gene_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

################### SPLICE
# ref splice factors
pdf("splice_splice.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_spl_splice_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_spl_splice_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_spl_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_spl_splice_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_spl_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_spl_splice_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref 3' end processing
pdf("splice_3.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_3_splice_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_3_splice_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_3_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_3_splice_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_3_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_3_splice_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref all splices
pdf("splice_all.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_all_splice_den_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_all_splice_den_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_all_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_all_splice_den_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_all_splice_den_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_all_splice_den_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

###################################### DIFF
################ GENE
# ref splice factors
pdf("gene_splice_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_spl_gene_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_spl_gene_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_spl_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_spl_gene_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_spl_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_spl_gene_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref 3' end processing
pdf("gene_3_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_3_gene_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_3_gene_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_3_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_3_gene_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_3_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_3_gene_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref all genes
pdf("gene_all_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_all_gene_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_all_gene_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_all_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_all_gene_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_all_gene_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_all_gene_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

################### SPLICE
# ref splice factors
pdf("splice_splice_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_spl_splice_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_spl_splice_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_spl_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_spl_splice_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_spl_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_spl_splice_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref 3' end processing
pdf("splice_3_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_3_splice_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_3_splice_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_3_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_3_splice_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_3_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_3_splice_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()

# ref all splices
pdf("splice_all_diff.pdf")
plot(-1, xlab="log2 fold change", ylab="density", xlim=c(-10, 10), ylim=c(0, 1.5))
ctrl <- sapply(names(ref_all_splice_den_cut_ls), function(group) {
	if (!(any(grep("*iPSC*", group)) && any(grep("*parental*", group)))) {
		lines(ref_all_splice_den_cut_ls[[group]])
		group
	}
})
smpl <- sapply(names(ref_all_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && !any(grep("*fam*", group))) {
		lines(ref_all_splice_den_cut_ls[[group]], col="red")
		group
	}
})
smpl_fam <- sapply(names(ref_all_splice_den_cut_ls), function(group) {
	if (any(grep("*iPSC*", group)) && any(grep("*parental*", group)) && any(grep("*fam*", group))) {
		lines(ref_all_splice_den_cut_ls[[group]], col="red", lty=2)
		group
	}
})
dev.off()



## work on legend with grep
legend("topright", legend=c(unlist(ctrl), unlist(smpl), unlist(smpl_fam)))


ref_spl_gene_den_cut_ls[[group]]

	### E. Plot density functions
	dens <- function(mt,out){
		# Construct output file name
		out_name <- paste(out_prefix, out, "_heatmap.pdf", sep="")
		# Save Data matrix to table
		write.table(mt, file=paste(out_prefix, out, "_tagwise_dispersion_matrix.tab", sep=""), quote=FALSE, sep="\t")

	}


	### FUNCTION CALLS!!!
	heat(ref_spl_gene, "ref_spl_gene")


	# Compute density functions of all and differentially expressed features
	dens_all <- density(fc, from=min(fc[is.finite(fc)]), to=max(fc[is.finite(fc)]))
		#dens_diff <- density(fc, from=min(fc[is.finite(fc)]), to=max(fc[is.finite(fc)]))
	# Plot
	pdf(file=paste(prefix, "_fc_densities.pdf", sep="_"))
	plot(dens_all, xlab="log2 fold change")
		#lines(dens_diff, lty=2, col="red")
	dump <- dev.off()
