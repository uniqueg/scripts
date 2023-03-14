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




















