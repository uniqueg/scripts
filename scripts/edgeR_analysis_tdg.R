#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		24-JAN-2013
### Modified:		25-JAN-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR
### Description: 	Differential gene expression analysis of two sample groups using the R/Bioconductor package edgeR
### Arguments: 		1./5. Paths to files of each group; 2./6. Pattern for file selection for each group; 3./7. number of samples for each group; 4./8. names for each group; 9. output file path (must exist!) 
### Output: 		BCV and smear plot; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage:			Rscript edgeR_analysis.R ./wt *.tab$ 2 wt ./ko *.tab$ 2 ko ./out
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
# Load library
library(edgeR)
###

### B. Import data
# Read indicated files from indicated directories
files <- c( dir(path=args[1], pattern=args[2], full.names=TRUE), dir(path=args[5], pattern=args[6], full.names=TRUE) )
files
# Assign indicated group names and sample numbers
group <- factor( c( rep.int(args[4], as.integer(args[3])) , rep.int(args[8], as.integer(args[7])) ) )
# Read files into DGEList object 'cts'
cts <- readDGE(files, group=group, labels=basename(files))
# Get count summary
summ_cts <- summary(cts$counts) 
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

## D. Plots
## D1. Tagwise dispersion vs log2(cpm)
pdf(file=paste(args[9], sep="/", "BCV_plot.pdf"), width = 6, height = 6)
plotBCV(cts_tag_wise_disp, cex=0.4)
dev.off()					
## D2. Tagwise log2(FC) vs log2(cpm) (~MA plot)
pdf(file=paste(args[9], sep="/", "smear_plot.pdf"), width = 6, height = 6)
detags <- rownames(cts_tag_wise_disp)[as.logical(decideTestsDGE(cts_exact_test))]
plotSmear(cts_exact_test, de.tags=detags)
# Blue lines indicate log2(FC) > 1 and < -1 
abline(h = c(-1, 1), col = "blue")
dev.off()
###	

### E. Write, save & clean-up
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
write.table(cts$counts, file=paste(args[9], sep="/", "counts_raw.tsv"), quote=FALSE, sep="\t")
write.table(cts_comm_disp$pseudo.counts, file=paste(args[9], sep="/", "counts_norm.tsv"), quote=FALSE, sep="\t")
write.table(cpm_de, file=paste(args[9], sep="/", "counts_norm_cpm.tsv"), quote=FALSE, sep="\t")
write.table(cts_exact_test$table, file=paste(args[9], sep="/", "diff_exp_all.tsv"), quote=FALSE, sep="\t")
write.table(tags_de$table, file=paste(args[9], sep="/", "diff_exp_fdr_cutoff.tsv"), quote=FALSE, sep="\t")
## E3. Remove unused/temp files & save workspace image/session
rm(files, group, detags)
save.image(file=paste(args[9], sep="/", "image.Rdata"))
###
