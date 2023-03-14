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
### Arguments: 		1./4. Paths to files of each group; 2./5. Pattern for file selection for each group; 3./6. names for each group; 7. Value for common dispersion; 8. Output file prefix (MUST exist; may include file path) 
### Output: 		Smear plots; tables of all and differentially expressed genes (FDR < 0.5); log (to STDOUT) 
### Usage example:	Rscript edgeR_analysis_no_repl.R ./wt *.tab$ 2 wt ./ko *.tab$ 2 ko 0.4 path/to/out/files/prefix
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
disp <- as.numerical(args[7])
prefix <- args[8]
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
# Read files into DGEList object 'dge'
dge <- readDGE(files, group=group, labels=basename(files), row.names=NULL)
###

### C. Calculation / analysis
# Calculate normalization factors
dge <- calcNormFactors(dge)
# Exact negative binomial tagwise tests
dgex <- exactTest(dge, dispersion=disp)
# Calculate differentially expressed
de <- decideTestsDGE(dgex)
# Summary de object
summ_de <- summary(de)
# Subset top tags (FDR < 0.05)
tags_de <- topTags(dgex, n=sum(summ_de[c(1,3)]))
###

### D. Write tables
write.table(dgex$table, file=paste(prefix, "diff_exp_all.tsv", sep="_"), quote=FALSE, sep="\t")
write.table(tags_de$table, file=paste(prefix, "diff_exp_fdr_cutoff.tsv", sep="_"), quote=FALSE, sep="\t")
###

## E. Plots
## E1. Smear plot (~MA; tagwise log2 FC vs log2 cpm)
pdf(file=paste(prefix, "smear_plot.pdf", sep="_"))
detags <- rownames(dge)[as.logical(de)]
plotSmear(dgex, de.tags=detags)
abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
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
dge$samples
cat("##\n\n### Count summary:\n")
summary(dge$counts)
cat("##\n\n### Count summary (counts per million):\n")
summary(cpm_de)
cat("##\n\n### Number of unique counts:\n")
dim(dge)[1]
cat("##\n\n### Sample comparison:\n")
dgex$comparison
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