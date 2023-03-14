#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		03-APR-2013
### Modified:		03-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR
### Description: 	Differential gene expression analysis of a DGE list object
### Arguments: 		1. Path to DGE list object; 2. Output file prefix 
### Output: 		BCV and smear plot; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage:			Rscript edgeR_analysis_dge.R ./dge_list_object.R prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
# Load library
library(edgeR)
## Pass arguments
dge_l_file <- args[1]
prefix <- args[2]
###


### B. Import data
# Load DGE list object
dge_l_file <- load(dge_l_file)
# Assign to variable
cts <- get(dge_l_file)
###


### C. Calculation / analysis
# Get count summary
summ_cts <- summary(cts$counts) 
# Get number of unique counts
no_cts <- dim(cts)[1]
# Calculate normalization factors
cts_norm_fact <- calcNormFactors(cts)
# Calculate common dispersion 
cts_comm_disp <- estimateCommonDisp(cts_norm_fact)
# Calculate tagwise dispersion
cts_tag_wise_disp <- estimateTagwiseDisp(cts_comm_disp)
# Exact negative binomial tagwise tests
cts_exact_test <- exactTest(cts_tag_wise_disp)
# Calculate differentially expressed

prefix <- "exons"
summ_de <- summary(decideTestsDGE(cts_exact_test))
# Subset top tags (FDR < 0.05)
tags_de <- topTags(cts_exact_test, n=sum(summ_de[c(1,3)]))
# Get count table normalized to counts per million
cpm_de <- cpm(cts_tag_wise_disp)[rownames(tags_de),]			## Subscript out of bounds
###


			### E. Write log and tables
			## E1. Write log
			writeLines("### Sample information\n")
			cts_norm_fact$samples
			writeLines("\n\n### Count summary\n")
			summ_cts
					writeLines("\n\n### Count summary (counts per million)\n")
					summary(cpm_de)
			writeLines("\n\n### Number of unique counts\n")
			no_cts
			writeLines("\n\n### Common dispersion\n")
			cts_comm_disp$common.dispersion
			writeLines("\n\n### Pseudo/normalized library size\n")
			cts_comm_disp$pseudo.lib.size
			writeLines("\n\n### Sample comparison\n")
			cts_exact_test$comparison
			writeLines("\n\n### Differentially expressed (FDR < 0.05)\n")
			sum(summ_de[c(1,3)])
			writeLines("\n\n### Downregulated, unchanged, upregulated\n")
			summ_de
			writeLines("\n\n### Session info\n")
			print.default(sessionInfo())
			## E2. Write tables
			write.table(cts$counts, file=paste(prefix, "counts_raw.tsv", sep="_"), quote=FALSE, sep="\t")
			write.table(cts_comm_disp$pseudo.counts, file=paste(prefix, "counts_norm.tsv", sep="_"), quote=FALSE, sep="\t")
					write.table(cpm_de, file=paste(prefix, "counts_norm_cpm.tsv", sep="_"), quote=FALSE, sep="\t")
			write.table(cts_exact_test$table, file=paste(prefix, "diff_exp_all.tsv", sep="_"), quote=FALSE, sep="\t")
			write.table(tags_de$table, file=paste(prefix, "diff_exp_fdr_cutoff.tsv", sep="_"), quote=FALSE, sep="\t")
			###

### D. Plots
## D1. Tagwise dispersion vs log2(cpm)
pdf(file=paste(prefix, "BCV_plot.pdf", sep="_"), width = 6, height = 6)
plotBCV(cts_tag_wise_disp, cex=0.4)
dev.off()					
## D2. Tagwise log2(FC) vs log2(cpm) (~MA plot)
pdf(file=paste(prefix, "smear_plot.pdf", sep="_"), width = 6, height = 6)
detags <- rownames(cts_tag_wise_disp)[as.logical(decideTestsDGE(cts_exact_test))]
plotSmear(cts_exact_test, de.tags=detags)
# Blue lines indicate log2(FC) > 1 and < -1 
abline(h = c(-1, 1), col = "blue")
dev.off()
###	


### E. Write log and tables
## E1. Write log
writeLines("### Sample information\n")
cts_norm_fact$samples
writeLines("\n\n### Count summary\n")
summ_cts
writeLines("\n\n### Count summary (counts per million)\n")
summary(cpm_de)
writeLines("\n\n### Number of unique counts\n")
no_cts
writeLines("\n\n### Common dispersion\n")
cts_comm_disp$common.dispersion
writeLines("\n\n### Pseudo/normalized library size\n")
cts_comm_disp$pseudo.lib.size
writeLines("\n\n### Sample comparison\n")
cts_exact_test$comparison
writeLines("\n\n### Differentially expressed (FDR < 0.05)\n")
sum(summ_de[c(1,3)])
writeLines("\n\n### Downregulated, unchanged, upregulated\n")
summ_de
writeLines("\n\n### Session info\n")
print.default(sessionInfo())
## E2. Write tables
write.table(cts$counts, file=paste(prefix, "counts_raw.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cts_comm_disp$pseudo.counts, file=paste(prefix, "counts_norm.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cpm_de, file=paste(prefix, "counts_norm_cpm.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cts_exact_test$table, file=paste(prefix, "diff_exp_all.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(tags_de$table, file=paste(prefix, "diff_exp_fdr_cutoff.tsv", sep="_"), quote=FALSE, sep="\t")
###


### F. Clean-up and save
# Remove unused/temp files
rm(files, group, detags)
# Save workspace image
save.image(file=paste(prefix, "image.Rdata", sep="_"))
###






















### Load data objects, rename and set names (if necessary)
# Load GenomicRanges library
library(GenomicRanges)
## GENES
load("summarized_overlaps_hg19_genes_abyzov_ipsc_and_parental.R")
# Rename object
genes <- overlaps
## EXONS
# Load exon data object
load("summarized_overlaps_hg19_exons_abyzov_ipsc_and_parental.R")
# Rename object
exons <- overlaps
# Set names attribute
names(rowData(exons)) <- rowData(exons)$name
###

### Define sample groups
parental <- c(1,5,9,10,14,18,22)
iPSC <- c(2:4,6:8,11:13,15:17,19:21,23:25)
###

### Prepare count tables
## Extract counts
counts_genes_parental <- assays(genes)$counts[,parental]
counts_genes_iPSC <- assays(genes)$counts[,iPSC]
counts_exons_parental <- assays(exons)$counts[,parental]
counts_exons_iPSC <- assays(exons)$counts[,iPSC]
## Shorten column names
colnames(counts_genes_parental) <- sub("segemehl_", "", basename(colnames(counts_genes_parental)))
colnames(counts_genes_iPSC) <- sub("segemehl_", "", basename(colnames(counts_genes_iPSC)))
colnames(counts_exons_parental) <- sub("segemehl_", "", basename(colnames(counts_exons_parental)))
colnames(counts_exons_iPSC) <- sub("segemehl_", "", basename(colnames(counts_exons_iPSC)))
###

# Remove unused objects
rm(overlaps, exons, genes, parental, iPSC)

## Save count tables
save(counts_genes_parental, counts_genes_iPSC, counts_exons_parental, counts_exons_iPSC, file="count_tables_abyzov.R")
write.table(counts_genes_parental, file = "counts_genes_parental", quote=FALSE, sep = "\t")
write.table(counts_genes_iPSC, file = "counts_genes_iPSC", quote=FALSE, sep = "\t")
write.table(counts_exons_parental, file = "counts_exons_parental", quote=FALSE, sep = "\t")
write.table(counts_exons_iPSC, file = "counts_exons_iPSC", quote=FALSE, sep = "\t")

## Load count tables
load("count_tables_abyzov.R")
counts_genes_parental <- as.matrix(read.table("counts_genes_parental"))
counts_genes_iPSC <- as.matrix(read.table("counts_genes_iPSC"))
counts_exons_parental <- as.matrix(read.table("counts_exons_parental"))
counts_exons_iPSC <- as.matrix(read.table("counts_exons_iPSC"))

### Prepare DGE list objects and save
# Load edgeR library
library(edgeR)
## Genes
counts_genes <- cbind(counts_genes_parental, counts_genes_iPSC)
genes_dge_l <- DGEList(counts=counts_genes, group=c(rep("parental",7),rep("iPSC",18)))
save(genes_dge_l, file="genes_dge_list_abyzov.R")
## Exons
counts_exons <- cbind(counts_exons_parental, counts_exons_iPSC)
exons_dge_l <- DGEList(counts=counts_exons, group=c(rep("parental",7),rep("iPSC",18)))
save(exons_dge_l, file="exons_dge_list_abyzov.R")
###

### STARTING POINT: Load library and objects
library(edgeR)
load("genes_dge_list_abyzov.R")
load("exons_dge_list_abyzov.R")
###

### Differential gene expression analysis (edgeR)

## GENES
genes_norm_fact <- calcNormFactors(genes_dge_l)
# Calculate common dispersion and 
genes_comm_disp <- estimateCommonDisp(genes_norm_fact)
# Calculate tagwise dispersion
genes_tag_wise_disp <- estimateTagwiseDisp(genes_comm_disp)
# Exact negative binomial tagwise tests
genes_exact_test <- exactTest(genes_tag_wise_disp)
# Calculate differentially expressed
summ_de_genes <- summary(decideTestsDGE(genes_exact_test))
# Subset top tags (FDR < 0.05)
tags_de_genes <- topTags(genes_exact_test, n=sum(summ_de_genes[c(1,3)]))
# Get count table normalized to counts per million
cpm_de_genes <- cpm(genes_tag_wise_disp)[rownames(tags_de_genes),]
## Write tables
write.table(genes_dge_l$counts, file="counts_raw_genes.tsv", quote=FALSE, sep="\t")
write.table(genes_comm_disp$pseudo.counts, file="counts_norm_genes.tsv", quote=FALSE, sep="\t")
write.table(cpm_de_genes, file="counts_norm_cpm_genes.tsv", quote=FALSE, sep="\t")
write.table(genes_exact_test$table, file="diff_expr_all_genes.tsv", quote=FALSE, sep="\t")
write.table(tags_de_genes$table, file="diff_exp_fdr_cutoff_genes.tsv", quote=FALSE, sep="\t")
## Tagwise dispersion vs log2(cpm)
pdf(file="BCV_plot_genes.pdf", width = 6, height = 6)
plotBCV(genes_tag_wise_disp, cex=0.4)
dev.off()					
## Tagwise log2(FC) vs log2(cpm) (~MA plot)
pdf(file="smear_plot_genes.pdf", width = 6, height = 6)
detags <- rownames(genes_tag_wise_disp)[as.logical(decideTestsDGE(genes_exact_test))]
plotSmear(genes_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")
dev.off()
## EXONS
exons_norm_fact <- calcNormFactors(exons_dge_l)
# Calculate common dispersion and 
exons_comm_disp <- estimateCommonDisp(exons_norm_fact)
# Calculate tagwise dispersion
exons_tag_wise_disp <- estimateTagwiseDisp(exons_comm_disp)
# Exact negative binomial tagwise tests
exons_exact_test <- exactTest(exons_tag_wise_disp)
# Calculate differentially expressed
summ_de_exons <- summary(decideTestsDGE(exons_exact_test))
# Subset top tags (FDR < 0.05)
tags_de_exons <- topTags(exons_exact_test, n=sum(summ_de_exons[c(1,3)]))
# Get count table normalized to counts per million
cpm_de_exons <- cpm(exons_tag_wise_disp)[rownames(tags_de_exons),]
## Sample comparison
pdf(file="MDS_plot_exons.pdf", width = 6, height = 6)
plotMDS(exons_comm_disp)
dev.off()
## Write tables
write.table(exons_dge_l$counts, file="counts_raw_exons.tsv", quote=FALSE, sep="\t")
write.table(exons_comm_disp$pseudo.counts, file="counts_norm_exons.tsv", quote=FALSE, sep="\t")
write.table(cpm_de_exons, file="counts_norm_cpm_exons.tsv", quote=FALSE, sep="\t")
write.table(exons_exact_test$table, file="diff_expr_all_exons.tsv", quote=FALSE, sep="\t")
write.table(tags_de_exons$table, file="diff_exp_fdr_cutoff_exons.tsv", quote=FALSE, sep="\t")
## Tagwise dispersion vs log2(cpm)
pdf(file="BCV_plot_exons.pdf", width = 6, height = 6)
plotBCV(exons_tag_wise_disp, cex=0.4)
dev.off()					
## Tagwise log2(FC) vs log2(cpm) (~MA plot)
pdf(file="smear_plot_exons.pdf", width = 6, height = 6)
detags <- rownames(exons_tag_wise_disp)[as.logical(decideTestsDGE(exons_exact_test))]
plotSmear(exibs_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")
dev.off()
###


#### LEFTOVERS ####

## Sample comparison
pdf(file="MDS_plot_genes.pdf", width = 6, height = 6)
plotMDS(genes_comm_disp)
dev.off()

mean_counts_genes_parental <- rowMeans(counts_genes_parental)
mean_counts_genes_iPSC <- rowMeans(counts_genes_iPSC)
mean_counts_exons_parental <- rowMeans(counts_exons_parental)
mean_counts_exons_iPSC <- rowMeans(counts_exons_iPSC)

sd_counts_genes_parental <- apply(counts_genes_parental, 1, sd)
sd_counts_genes_iPSC <- apply(counts_genes_iPSC, 1, sd)
sd_counts_exons_parental <- apply(counts_exons_parental, 1, sd)
sd_counts_exons_iPSC <- apply(counts_exons_iPSC, 1, sd)

####



















