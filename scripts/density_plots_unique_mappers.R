#######
### GENERAL:
### --------
### Author: 		Christina Herrmann
### Created: 		30-MAY-2013
### Modified:		20-MAY-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Density plots for fold change output from edgeR analysis,
###					all genes and subsets containing only splicing factors are used
### Arguments: 		-
### Output: 		Density plots 
### Usage example:	Rscript ./density_plots_unique_mappers
#######

# Load objects from edgeR analysis of unique mappers (fam03, fams1123, all ...)
load("abyzov_for_plots_image.Rdata")

# Read in list of ENSEMBL gene_ids of splice factors and 3' end processing factors
# manual list
splice_names <- read.delim("/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/output/unique_mappers/plots/splice_factors_man_ENSEMBL_gene_names.tab",header=FALSE)
splice_names_char <- as.character(splice_names[,1])
# GO term list
splice_GO_names <- read.delim("/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/output/unique_mappers/plots/splice_factors_GO_ENSEMBL_gene_names.tab",header=FALSE)
splice_GO_names_char <- as.character(splice_GO_names[,1])
# 3' end factors
utr_names <- read.delim("3_end_names.tab",header=FALSE)
utr_names_char <- as.character(utr_names[,1])

# Subset splice factors from all_topTags
splice <- na.omit(all_topTags$table[splice_names_char,])
utr <- na.omit(all_topTags$table[utr_names_char,])
splice_GO <- na.omit(all_topTags$table[splice_GO_names_char,])


# Subset splice factors from all_tags_de
splice_de <- na.omit(all_tags_de$table[splice_names_char,])
utr_de <- na.omit(all_tags_de$table[utr_names_char,])
splice_GO_de <- na.omit(all_tags_de$table[splice_GO_names_char,])

# Density plots
# all genes topTags vs. all genes differentially expressed

all_FC <- all_topTags$table[,1]
all_de_FC <- all_tags_de$table[,1]

all_dens <- density(all_FC, from=-12, to=12)
all_de_dens <- density(all_de_FC, from=-12, to=12)

png("density_all_VS_all_de.png")
plot(all_dens, xlab="log2 fold change", ylab="density",lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All genes\nvs\ndifferentially expressed genes")
lines(all_de_dens)
dev.off()

# splice factors topTags vs. splice factors differentially expressed

splice_FC <- splice[,1]
splice_de_FC <- splice_de[,1]

splice_dens <- density(splice_FC, from=-12, to=12)
splice_de_dens <- density(splice_de_FC, from=-12, to=12)

png("density_splice_VS_splice_de.png")
plot(splice_dens, xlab="log2 fold change", ylab="density", col="red", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="Splice factors\nvs\ndifferentially expressed splice factors")
lines(splice_de_dens, col="red")
dev.off()

# all genes topTags vs. splice factors topTags

png("density_all_VS_splice.png")
plot(all_dens, xlab="log2 fold change", ylab="density", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All genes\nvs\nsplice factors")
lines(splice_dens, col="red", lty=2)
dev.off()

# all genes differentially expressed vs. splice factors differentially expressed

png("density_all_de_VS_splice_de.png")
plot(all_de_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All differentially expressed genes\nvs\ndifferentially expressed splice factors")
lines(splice_de_dens, col="red")
dev.off()

# all vs all de vs. splice vs splice de

png("density_all_VS_all_de_VS_splice_VS_splice_de.png")
plot(all_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All differentially expressed genes\nvs\ndifferentially expressed splice factors", lty=2)
lines(all_de_dens)
lines(splice_dens, col="red", lty=2)
lines(splice_de_dens, col="red")
dev.off()

# UTR topTags vs. UTR differentially expressed

utr_FC <- utr[,1]
utr_de_FC <- utr_de[,1]

utr_dens <- density(utr_FC, from=-12, to=12)
utr_de_dens <- density(utr_de_FC, from=-12, to=12)

png("density_utr_VS_utr_de.png")
plot(utr_dens, xlab="log2 fold change", ylab="density", col="blue", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="3'end factors\nvs\ndifferentially expressed 3'end factors")
lines(utr_de_dens, col="blue")
dev.off()

# all genes topTags vs. UTR topTags
png("density_all_VS_utr.png")
plot(all_dens, xlab="log2 fold change", ylab="density", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All genes\nvs\n3'end factors")
lines(utr_dens, col="blue", lty=2)
dev.off()

# all genes differentially expressed vs. UTR differentially expressed
png("density_all_de_VS_utr_de.png")
plot(all_de_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All differentially expressed genes\nvs\ndifferentially expressed 3'end factors")
lines(utr_de_dens, col="blue")
dev.off()

# all vs all de vs. splice vs splice de vs utr vs utr de

png("density_all_VS_all_de_VS_splice_VS_splice_de_VS_utr_VS_utr_de.png")
plot(all_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All genes\nvs\n splice factors\nvs\n3' end factors", lty=2)
lines(all_de_dens)
lines(splice_dens, col="red", lty=2)
lines(splice_de_dens, col="red")
lines(utr_dens, col="blue", lty=2)
lines(utr_de_dens, col="blue")
dev.off()

#all de vs. splice de vs. utr de

png("density_all_de_VS_splice_de_VS_utr_de.png")
plot(all_de_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All DE genes\nvs\nDE splice factors vs DE 3' end factors")
lines(splice_de_dens, col="red")
lines(utr_de_dens, col="blue")
dev.off()

# splice factors (GO) topTags vs. splice factors (GO) differentially expressed

splice_GO_FC <- splice_GO[,1]
splice_GO_de_FC <- splice_GO_de[,1]

splice_GO_dens <- density(splice_GO_FC, from=-12, to=12)
splice_GO_de_dens <- density(splice_GO_de_FC, from=-12, to=12)

png("density_splice_GO_VS_splice_GO_de.png")
plot(splice_GO_dens, xlab="log2 fold change", ylab="density", col="green", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="Splice factors (GO)\nvs\ndifferentially expressed splice factors (GO)")
lines(splice_GO_de_dens, col="green")
dev.off()

# all genes topTags vs. splice factors (GO) topTags

png("density_all_VS_splice_GO.png")
plot(all_dens, xlab="log2 fold change", ylab="density", lty=2, xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All genes\nvs\nsplice factors (GO)")
lines(splice_GO_dens, col="green", lty=2)
dev.off()

# all genes differentially expressed vs. splice factors (GO) differentially expressed

png("density_all_de_VS_splice_GO_de.png")
plot(all_de_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2),
		main="All differentially expressed genes\nvs\ndifferentially expressed splice factors (GO)")
lines(splice_GO_de_dens, col="green")
dev.off()

# Splice man vs splice GO
png("density_splice_man_VS_splice_GO_de.png")
plot(splice_dens, xlab="log2 fold change", ylab="density", xlim=c(-12, 12), ylim=c(0, 1.2), col="red", lty=2,
		main="Splice factors (manual list)\nvs\nsplice factors (GO)")
lines(splice_de_dens, col="red")
lines(splice_GO_dens, col="green", lty=2)
lines(splice_GO_de_dens, col="green")
dev.off()

