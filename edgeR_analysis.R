#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		05-APR-2013
### Modified:		05-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Differential gene expression analysis of two sample groups using the R/Bioconductor package edgeR
### Arguments: 		1./4. Paths to files of each group; 2./5. Pattern for file selection for each group; 3./6. names for each group; 7. output file prefix (MUST exist; may include file path) 
### Output: 		BCV, MDS and smear plots; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage example:	Rscript edgeR_analysis.R ./wt *.tab$ 2 wt ./ko *.tab$ 2 ko path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
path1 <- args[1]
pattern1 <- args[2]
group1 <- args[3]
path2 <- args[4]
pattern2 <- args[5]
group2 <- args[6]
prefix <- args[7]
# Load library
suppressMessages(library(edgeR))
###

### B. Import data
## Read indicated files from indicated directories
files1 <- dir(path=path1, pattern=pattern1, full.names=TRUE)
files2 <- dir(path=path2, pattern=pattern2, full.names=TRUE)
files <- c(files1, files2)
# Assign indicated group names and sample numbers
group <- factor( c( rep.int(group1, length(files1)) , rep.int(group2, length(files2)) ) )
# Read files into DGEList object 'cts'
cts <- readDGE(files, group=group, labels=basename(files), row.names=NULL)
# Get number of unique counts
no_cts <- dim(cts)[1]
###

### C. Calculation / analysis
# Calculate normalization factors
cts_norm_fact <- calcNormFactors(cts)
# Calculate common dispersion and 
cts_comm_disp <- estimateCommonDisp(cts_norm_fact)
# Calculate tagwise dispersion
cts_tag_wise_disp <- estimateTagwiseDisp(cts_comm_disp)
# Exact negative binomial tagwise tests
cts_exact_test <- exactTest(cts_tag_wise_disp)
# Calculate differentially expressed
summ_de <- summary(decideTestsDGE(cts_exact_test))
# Subset top tags (FDR < 0.05)
tags_de <- topTags(cts_exact_test, n=sum(summ_de[c(1,3)]))
# Get count table normalized to counts per million
cpm_de <- cpm(cts_tag_wise_disp)[rownames(tags_de),]
###

### D. Write tables
write.table(cts$counts, file=paste(prefix, "counts_raw.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cts_comm_disp$pseudo.counts, file=paste(prefix, "counts_norm.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cpm_de, file=paste(prefix, "counts_norm_cpm.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(cts_exact_test$table, file=paste(prefix, "diff_exp_all.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(tags_de$table, file=paste(prefix, "diff_exp_fdr_cutoff.tsv", sep="_"), quote=FALSE, sep="\t")
###

## E. Plots
## E1. Smear plot (~MA; tagwise log2 FC vs log2 cpm)
pdf(file=paste(prefix, "smear_plot.pdf", sep="_"), width = 6, height = 6)
detags <- rownames(cts_tag_wise_disp)[as.logical(decideTestsDGE(cts_exact_test))]
plotSmear(cts_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
dev.off()
## E2. Biological covariance plot (BCV; tagwise dispersion vs log2 cpm)
pdf(file=paste(prefix, "BCV_plot.pdf", sep="_"), width = 6, height = 6)
plotBCV(cts_tag_wise_disp, cex=0.4)
dev.off()
## E3. Multidimensional scaling plot (MDS; sample relations)
pdf(file=paste(prefix, "MDS_plot.pdf", sep="_"), width = 6, height = 6)
plotMDS(cts_comm_disp, labels=NULL, col=ifelse(cts$samples[,2] == group1, "blue", "red"))
dev.off()
###

### F. Write lig, clean up and save image
## Write log
cat("\n### Files read for group 1:\n")
files1
cat("# Total number:\n")
length(files1)
cat("##\n\n### Files read for group 2:\n")
files2
cat("# Total number:\n")
length(files2)
cat("##\n\n### Sample information:\n")
cts$samples
cat("##\n\n### Count summary:\n")
summary(cts$counts)
cat("##\n\n### Count summary (counts per million):\n")
summary(cpm_de)
cat("##\n\n### Number of unique counts:\n")
no_cts
cat("##\n\n### Common dispersion:\n")
cts_comm_disp$common.dispersion
cat("##\n\n### Pseudo/normalized library size:\n")
cts_comm_disp$pseudo.lib.size
cat("##\n\n### Sample comparison:\n")
cts_exact_test$comparison
cat("##\n\n### Differentially expressed (FDR < 0.05):\n")
sum(summ_de[c(1,3)])
cat("##\n\n### Downregulated, unchanged, upregulated:\n")
summ_de
cat("##\n\n### Session info:\n")
print.default(sessionInfo())
cat("##\n\n")
# Remove unused/temp files
rm(files, group, detags)
# Save workspace image/session
save.image(file=paste(prefix, "image.Rdata", sep="_"))
###
