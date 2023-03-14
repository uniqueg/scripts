## For each tissue / cell line / clinical sample (n > 1) / iPSC & parent family / iPSC individual
## For all cell lines / clinical samples / iPSC & parent
# Plot MDS for each
# Plot MDS for all
# Write table of common dispersions for all (alternative: first value of tagwise dispersion table is common dispersion
# Heat map
# Run script for ALL transcripts, splice factors & 3'end processing subsets

##! What to implement
# Differential expression analysis GENES (edgeR; decide exact tests)
# Differential expression analysis EXONS (edgeR; decide exact tests)
# Normalized fold changes per EXON (Difference log fold changes EXONS-GENES)
# Confidence of alternative splicing per GENE (edgeR::spliceVariants)

##? What output and how to plot/visualize?
#? For each individual analysis?
#? For all analyses together? 

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
df_files <- args[1]
prefix <- args[2]
# Load library
suppressMessages(library(edgeR))
###

### B. Import data
# Import file info table from file
#input_df <- read.table(df_files, header=FALSE, sep="\t", colClasses=c("character", "factor", "factor", "integer", "character"), na.strings=c("NA",""), fill=TRUE)
input_df <- read.table("df_test", header=FALSE, sep="\t", colClasses=c("character", "factor", "character"), na.strings=c("NA",""), fill=TRUE)
## Traverse file info data frame row by row; returns list of data frames of (absolute) filenames and groups
ls <- apply(input_df, 1, function(row) {
	## If no pattern for file selection is given
	if (is.na(row[3])) {
		# Assign given absolut file location to files
		files <- row[1]
	} else {
		# Else derive one or more files from indicated path and pattern
		files <- dir(path=row[1], pattern=glob2rx(row[3]), full.names=TRUE)
	}
	# Construct and return data frame
	data.frame(file=files, group=row[2], row.names=basename(files))
})
###

### C. MDS PLOT ALL GROUPS
# Build data frame containg all groups from list of data frames
df <- do.call(rbind,ls)
# Read all files into a DGEList object and append to list of DGEList objects
dge_all <- readDGE(as.character(df$file), group=df$group, labels=basename(as.character(df$file)), columns=c(1,2), row.names=NULL)
## Plot MDS
pdf(paste(prefix, "MDS_all_groups", sep="_"), height=6, width=6)
plotMDS(dge_all, labels="o", col=as.numeric(dge_all$samples$group), main="Multidimensional scaling plot\nAll samples")
#legend("topright", inset=c(-0.2,0), legend=levels(dge_all$samples$group), fill=1:length(levels(dge_all$samples$group)))
dev.off()
###

### D. GROUPWISE MDS PLOTS AND DISPERSION TABLES
# Read files of individual groups into a list of DGEList objects
dge_ls <- lapply(ls, function(group) dge <- readDGE(as.character(group$file), group=group$group, labels=basename(as.character(group$file)), columns=c(1,2), row.names=NULL))
## Apply over each DGEList object
df_ls <- lapply(dge_ls, function(dge) {
	## Plot MDS
	#pdf("test.pdf", height=6, width=6)
	#plotMDS(dge, labels="o", col=as.numeric(dge$samples$group), main=paste("Multidimensional scaling plot\nGroup: ",levels(dge$samples$group)))
	#legend("topright", inset=c(-0.2,0), legend=levels(dge$samples$group), fill=1:length(levels(dge$samples$group)))
	# Normalize library sizes
	dge <- calcNormFactors(dge)
	## Calculate common and tagwise dispersions
	dge <- estimateCommonDisp(dge)
	dge <- estimateTagwiseDisp(dge)
	# Construct data frame of common and tagwise dispersions
	df <- data.frame(c(dge$common.dispersion, dge$tagwise.dispersion), row.names=c("common dispersion", row.names(dge$counts)))
	# Write to file
	write.table(df, file=paste(prefix, "_dispersion_", levels(dge$samples$group), sep="_") , quote=FALSE, sep="\t")
	# Return data frame of 
	return(df)
})

### E. HEATMAP ALL GROUPS
# Combine dispersion value data frames for all libraries into one common matrix
mt_all <- as.matrix(do.call(cbind, df_ls))
# Set group names as column names
colnames(mt_all) <- levels(dge_all$samples$group)
## Plot heatmap



# Read files of class 'exon' into DGEList object 'dge_exons'
#dge_exons <- readDGE(as.character(ls$exons$file), group=ls$exons$group, labels=basename(as.character(ls$exons$file)), columns=c(1,2), row.names=NULL)
# Read files of class 'exon' into DGEList object 'dge_genes'
	## Update library sizes
	#dge_exons$samples$lib.size <- ls$exons$size
	#dge_genes$samples$lib.size <- ls$genes$size
# Remove temporary objects
rm(input_df, df, ls)
###

# for genes: DGE

design_genes <- model.matrix(~dge_genes$samples$group, data=dge_genes$samples)
dge_genes <- calcNormFactors(dge_genes)
dge_genes <- estimateGLMCommonDisp(dge_genes, design_genes)
dge_genes <- estimateGLMTrendedDisp(dge_genes, design_genes)
dge_genes <- estimateGLMTagwiseDisp(dge_genes, design_genes)
fit_genes <- glmFit(dge_genes, design_genes)
lrt_genes <- glmLRT(fit_genes, coef=2:dim(design_genes)[2])


# for exons: DGE + spliceVariants


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

### C. Calculation / analysis
# Calculate normalization factors
cts_norm_fact <- calcNormFactors(cts)
# Calculate common dispersion and 
cts_comm_disp <- estimateCommonDisp(cts_norm_fact)
# Calculate tagwise dispersion
cts_tag_wise_disp <- estimateTagwiseDisp(cts_comm_disp)
# Calculate alternatively spliced genes
cts_exact_test <- spliceVariants(cts_tag_wise_disp, feature_names, dispersion=NULL, group=group, estimate.genewise.disp=TRUE, trace=TRUE)
# DOES NOT GIVE MEANINGFUL RESULT FOR SPLICEVARIANT! // Calculate differentially expressed
#summ_de <- summary(decideTestsDGE(cts_exact_test))
# Subset top tags (FDR < 0.05)
tags_de <- topTags(cts_exact_test, n=sum(cts_exact_test$table$PValue < 0.05))
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
dim(cts)[1]
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







#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		05-APR-2013
### Modified:		05-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Compute relative log2 fold changes between two edgeR analyses (query and subject); query feature names are tested for suffixes added by base::make.unique ([name].[serial_number]) before looking up corresponding values in subject to allow for multiple queries matching to a single subject (e.g. exon to gene comparison) 
### Arguments: 		1. Query file; 2. Subject file; 3. Output file prefix (may contain path!); 1./2. are tables of differentially expressed features generated by "edgeR_analysis.R" (column 1 = features [rownames], column 2 = log2 FC); possible duplicate query feature names are allowed in principal, but have to have beeen made unique with base::make.unique BEFORE using this script!
### Output: 		Table with relative log2 fold changes for each query feature 
### Usage example:	Rscript compare_edgeR.R ./query.tab ./subject.tab path/to/output/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
query <- args[1]
subject <- args[2]
prefix <- args[3]
###

### B. Load data
# Load query data frame
query_df <- read.table(query)
# Load subject data frame
subject_df <- read.table(subject)
###

### C. Generate table of relative fold changes
# Create (ambiguous) query feature names from unique identifiers
feature_names <- sub('\\.\\d+$', "", rownames(query_df))
# Create vector of subject indices for each query feature name
lookup_indices <- match(feature_names, rownames(subject_df))
# Calculate relative log2 fold changes between query and subject (differences of log2 query and log2 subject)
rel_fc_query <- query_df[,1] - subject_df[lookup_indices,1]
# Bind (ambiguous) identifiers and relative log2 fold changes to data.frame
rel_fc_query_df <- data.frame("Features"=feature_names, "Log2FC"=rel_fc_query, row.names=row.names(query_df), stringsAsFactors=FALSE, check.rows=TRUE)
###

### D. Write output table and plot results
# Write table of relative fold changes
write.table(rel_fc_query_df, file=paste(prefix, "relative_fold_changes.tsv", sep="_"), quote=FALSE, sep="\t")
# Write summary to <STDOUT>
quantile(rel_fc_query_df[,2], c(.01, .02, .05, .1, .15, .25, .5, .75, .85, .90, .95, .98, .99))
# Plot histogram
pdf(file=paste(prefix, "relative_fold_changes.pdf", sep="_"), width = 6, height = 6)
breaks <- seq(from=floor(min(rel_fc_query_df[,2])), to=ceiling(max(rel_fc_query_df[,2])), by=1)
hist(rel_fc_query_df[,2], breaks=breaks, freq=FALSE, xlab="Relative log2 fold change")
dev.off()
###

### E. Write output & save session image
# Remove unused/temp variables
rm(feature_names, lookup_indices, rel_fc_query, breaks)
# Save workspace image/session
save.image(file=paste(prefix, "image.Rdata", sep="_"))
###