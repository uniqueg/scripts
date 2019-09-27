#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		05-APR-2013
### Modified:		28-MAY-2013 Christina Herrmann
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Differential gene expression analysis of sample groups using the R/Bioconductor package edgeR
###			Input from single count table for all samples
### Arguments: 		1. absolute file name for count table; 2. output file prefix (MUST exist; may include file path) 
### Output: 		BCV, MDS and smear plots; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage example:	Rscript edgeR_analysis_abyzov_unique.R ./count_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
counts <- as.character(args[1])
prefix <- args[2]
# Load library
suppressMessages(library(edgeR))
###

### B. Import data
## load count table
counts <- read.table(counts)
# Assign indicated group names and sample numbers
group <- factor( c( "fam03_parental", rep("fam03_iPSC", 3), "fam03_parental", rep("fam03_iPSC", 3), "fam03_parental", "fam03_parental", rep("fam03_iPSC", 3), "fams1123_parental", rep("fams1123_iPSC", 3), "fams1123_parental", rep("fams1123_iPSC", 3),"fams1123_parental", rep("fams1123_iPSC", 3), "H1_iPSC" ))
# Read count table into DGEList object 'dge'
dge <- DGEList(counts, group=group)
###

### C. Calculation / analysis
# Calculate normalization factors
dge <- calcNormFactors(dge)
# Calculate common dispersion and 
dge <- estimateCommonDisp(dge)
# Calculate tagwise dispersion
dge <- estimateTagwiseDisp(dge)
# Exact negative binomial tagwise tests, pairs to compare have to be specified
fam03_exact_test <- exactTest(dge, pair=c("fam03_parental","fam03_iPSC"))
fams1123_exact_test <- exactTest(dge, pair=c("fams1123_parental","fams1123_iPSC"))
#parental_exact_test <- exactTest(dge, pair=c("fam03_parental","fams1123_parental"))
#ipsc_exact_test <- exactTest(dge, pair=c("fam03_iPSC","fams1123_iPSC"))
# iPSC all vs. parental all
all <- dge
all$samples$group <- as.factor(sub("fam03_", "", dge$samples$group))
all$samples$group <- as.factor(sub("fams1123_", "", dge$samples$group))
all_exact_test <- exactTest(all, pair=c("parental","iPSC"))
# Calculate differentially expressed
fam03_summ_de <- summary(decideTestsDGE(fam03_exact_test))
fams1123_summ_de <- summary(decideTestsDGE(fams1123_exact_test))
#parental_summ_de <- summary(decideTestsDGE(parental_exact_test))
#ipsc_summ_de <- summary(decideTestsDGE(ipsc_exact_test))
all_summ_de <- summary(decideTestsDGE(all_exact_test))
# Subset top tags (FDR < 0.05)
fam03_tags_de <- topTags(fam03_exact_test, n=sum(fam03_summ_de[c(1,3)]))
fams1123_tags_de <- topTags(fams1123_exact_test, n=sum(fams1123_summ_de[c(1,3)]))
all_tags_de <- topTags(all_exact_test, n=sum(all_summ_de[c(1,3)]))
# Get count table normalized to counts per million
#cpm_de <- cpm(dge)[rownames(tags_de),]
###

### D. Write tables
write.table(all_exact_test$table, file=paste(prefix, "edgeR_out_all.tab", sep="_"), quote=FALSE, sep="\t")
write.table(fam03_exact_test$table, file=paste(prefix, "edgeR_out_fam03.tab", sep="_"), quote=FALSE, sep="\t")
write.table(fams1123_exact_test$table, file=paste(prefix, "edgeR_out_fams1123.tab", sep="_"), quote=FALSE, sep="\t")

write.table(fam03_tags_de$table, file=paste(prefix, "diff_exp_fam03.tab", sep="_"), quote=FALSE, sep="\t")
write.table(fams1123_tags_de$table, file=paste(prefix, "diff_exp_fams1123.tab", sep="_"), quote=FALSE, sep="\t")
write.table(all_tags_de$table, file=paste(prefix, "diff_exp_all.tab", sep="_"), quote=FALSE, sep="\t")
###

## E. Plots
## E1. Smear plot (~MA; tagwise log2 FC vs log2 cpm)
# all
pdf(file=paste(prefix, "smear_plot_all.pdf", sep="_"), width = 6, height = 6)
detags <- rownames(all)[as.logical(decideTestsDGE(all_exact_test))]
plotSmear(all_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
dev.off()
# fam03
pdf(file=paste(prefix, "smear_plot_fam03.pdf", sep="_"), width = 6, height = 6)
detags <- rownames(dge)[as.logical(decideTestsDGE(fam03_exact_test))]
plotSmear(fam03_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
dev.off()
#fams1123
pdf(file=paste(prefix, "smear_plot_fams1123.pdf", sep="_"), width = 6, height = 6)
detags <- rownames(dge)[as.logical(decideTestsDGE(fams1123_exact_test))]
plotSmear(fams1123_exact_test, de.tags=detags)
abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
dev.off()


## E2. Biological covariance plot (BCV; tagwise dispersion vs log2 cpm)
#pdf(file=paste(prefix, "BCV_plot.pdf", sep="_"), width = 6, height = 6)
#plotBCV(cts_tag_wise_disp, cex=0.4)
#dev.off()

## E3. Multidimensional scaling plot (MDS; sample relations)
palette(c("yellow2","orchid","orange","tomato","tomato4","pink","peachpuff","turquoise","plum","steelblue1","royalblue4","lightblue","palegreen","brown","deeppink2","red","tomato2","yellow4","seagreen","springgreen2","darkgreen","wheat4","grey","black","lightgrey"))
pdf(file=paste(prefix, "MDS_plot_all.pdf", sep="_"), width = 6, height = 6)
plotMDS(dge, labels=NULL, col=as.numeric(as.factor(group)))
legend("topright", inset=c(-0.45,0), legend=levels(as.factor(group)), fill=1:length(levels(as.factor(group))))
dev.off()
###

### F. Write log, clean up and save image
## Write log
#cat("\n### Files read for group 1:\n")
#files1
#cat("# Total number:\n")
#length(files1)
#cat("##\n\n### Files read for group 2:\n")
#files2
#cat("# Total number:\n")
#length(files2)
#cat("##\n\n### Sample information:\n")
#cts$samples
#cat("##\n\n### Count summary:\n")
#summary(cts$counts)
#cat("##\n\n### Count summary (counts per million):\n")
#summary(cpm_de)
#cat("##\n\n### Number of unique counts:\n")
#dim(cts)[1]
#cat("##\n\n### Common dispersion:\n")
#cts_comm_disp$common.dispersion
#cat("##\n\n### Pseudo/normalized library size:\n")
#cts_comm_disp$pseudo.lib.size
#cat("##\n\n### Sample comparison:\n")
#cts_exact_test$comparison
#cat("##\n\n### Differentially expressed (FDR < 0.05):\n")
#sum(summ_de[c(1,3)])
#cat("##\n\n### Downregulated, unchanged, upregulated:\n")
#summ_de
#cat("##\n\n### Session info:\n")
#print.default(sessionInfo())
#cat("##\n\n")
# Remove unused/temp files
#rm(files, group, detags)
# Save workspace image/session
save.image(file=paste(prefix, "image.Rdata", sep="_"))
###
