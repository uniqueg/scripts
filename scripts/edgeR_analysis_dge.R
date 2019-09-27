#######
### GENERAL:
### --------
### Author: Alexander Kanitz
### Created: 03-APR-2013
### Modified: 03-APR-2013
### Language: R
### Version: 2.15.2
### Requirements: Bioconductor_2.11, edgeR
### Description: Differential gene expression analysis of a DGE list object
### Arguments: 1. Path to DGE list object; 2. Output file prefix 
### Output: BCV and smear plot; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage: Rscript edgeR_analysis_dge.R ./dge_list_object.R prefix
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
summ_de <- summary(decideTestsDGE(cts_exact_test))
# Subset top tags (FDR < 0.05)
tags_de <- topTags(cts_exact_test, n=sum(summ_de[c(1,3)]))
# Get count table normalized to counts per million
cpm_de <- cpm(cts_tag_wise_disp)[rownames(tags_de),]
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
