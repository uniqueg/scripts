#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		21-DEC-2012
### Modified:		21-DEC-2012
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR_3.0.7, limma_3.14.3
#######

#######
### DESCRIPTION:
### ------------
### Differential gene expression analysis of murine embryonic stem cells with a heterozygous (+/-) or homozygous (-/-) knockout for TDG. The differential analysis was performed using the R/Bioconductor package edgeR.
#######

#######
### OUTPUT:
### -------
### BCV and smear plot; table of differentially expressed genes; various other count tables (see section E2.); log (to STDOUT)
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Import data
### C. Calculation / analysis
### D. Plots
### E. Write, save & clean-up
#######

### A. Pre-requisites
# Load library
library(edgeR)
###

### B. Import data
## Read .tab files from current working directory and merge into DGEList object 'cts'
files <- dir(pattern="*.tab$")
group <- factor(c("TDG+/-", "TDG+/-", "TDG-/-", "TDG-/-"))
labels <- c("TDG+/-_1", "TDG+/-_2", "TDG-/-_1", "TDG-/-_2")
cts <- readDGE(files, group=group, labels=labels)
summ_cts <- summary(cts$counts) 
no_cts <- dim(cts)[1]
###

### C. Calculation / analysis
## C1. Calculate normalization factors, common/tagwise dispersion and do exact negative binomial tagwise tests
cts_norm_fact <- calcNormFactors(cts)
cts_comm_disp <- estimateCommonDisp(cts_norm_fact)
cts_tag_wise_disp <- estimateTagwiseDisp(cts_comm_disp)
cts_exact_test <- exactTest(cts_tag_wise_disp)
## C2. Extract/subset differentially expressed genes (FDR < 0.05)
summ_de <- summary(decideTestsDGE(cts_exact_test))
tags_de <- topTags(cts_exact_test, n=sum(summ_de[c(1,3)]))
cpm_de <- cpm(cts_tag_wise_disp)[rownames(tags_de),]
###

## D. Plots
## D1. Tagwise dispersion vs log2(cpm)
pdf("BCV_plot.pdf", width = 6, height = 6)
plotBCV(cts_tag_wise_disp, cex=0.4)
dev.off()					
## D2. Tagwise log2(FC) vs log2(cpm) (~MA plot)
pdf("smear_plot.pdf", width = 6, height = 6)
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
write.table(cts$counts, file="counts_raw.tsv", quote=FALSE)
write.table(cts_comm_disp$pseudo.counts, file="counts_norm.tsv", quote=FALSE)
write.table(cpm_de, file="counts_norm_cpm.tsv", quote=FALSE)
write.table(cts_exact_test$table, file="diff_exp_all.tsv", quote=FALSE)
write.table(tags_de$table, file="diff_exp_fdr_cutoff.tsv", quote=FALSE)
## E3. Remove unused/temp files & save workspace image/session
rm(files, group, labels, detags)
save.image(file = "image.Rdata")
###
