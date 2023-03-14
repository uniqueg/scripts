#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		22-APR-2013
### Modified:		22-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Differential gene expression analysis of two sample groups using the R/Bioconductor package edgeR
### Arguments: 		1. Tab-separated file (no header!) of format: Group 1 (1) name, (2) path, (3) pattern, (4) recursive [TRUE or FALSE]; group 2 (5) name, (6) path, (7) pattern, (8) recursive [FILE | TAB/TSV]; 2. P value cutoff (e.g. 0.05) for writing tables of differentially expressed/alternatively spliced genes [NUMERIC]; 3. Character representation of logical ["TRUE" | "FALSE"] indicating whether a splice variant analysis should be performed (counts need to be supplied individually for each exon; in this implementation, exon identifiers need to be derived from corresponding gene names with limma::makeUnique, i.e. of the form 'gene', 'gene.1', 'gene.2' etc. for consecutive exons); 4. Base path to the directory in which output files are written (a separate folder is produced here for each individual pairwise comparison!); 5. Prefix that is prepended to each output filename [CHARACTER] 
### Output: 		BCV, MDS and smear plots; table of all features (fold change, counts per million, P value, FDR; if argument 'splice' == "TRUE" then also P value and FDR with regards to alternative splicing evidence); table of differentially expressed (and - if argument 'splice' == "TRUE" - alternatively spliced) genes (P value < indicated cutoff); log (to STDOUT) 
### Usage example:	Rscript edgeR_pairwose.R ./input_files.tab 0.05 TRUE ./output first_run_
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
df_files <- args[1]
p_cutoff <- as.numeric(args[2])
splice <- as.logical(args[3])
out_folder <- args[4]
out_prefix <- args[5]
# Load library
suppressMessages(library(edgeR))
###

### B. Import data
# Import file info table from file
input_df <- read.table(df_files, header=FALSE, sep="\t", colClasses=c("character", "factor", "character"), na.strings=c("NA",""), fill=TRUE)
## Run pairwise edgeR differential gene expression analysis for each row in 'input_df' data frame
ls <- apply(input_df, 1, function(row) {
	# Generate list of objects to be saved for later use
	objects <- list()
	## Add important objects to 'objects' list
	objects$row <- row
	objects$p_cutoff <- p_cutoff
	objects$splice <- splice
	## For each comparison group, read files specified by path and pattern
	objects$files1 <- dir(path=row[2], pattern=glob2rx(row[3]), recursive=row[4], full.names=TRUE)
	objects$files2 <- dir(path=row[6], pattern=glob2rx(row[7]), recursive=row[8], full.names=TRUE)
	# Assign indicated group names and sample numbers
	objects$group <- factor(c(rep.int(row[1], length(objects$files1)), rep.int(row[5], length(objects$files2))))
	# Read files into DGEList object
	objects$dge <- readDGE(c(objects$files1, objects$files2), group=objects$group, labels=basename(c(objects$files1, objects$files2)), row.names=NULL)	
	###
	
	### C. DGE analysis
	# Normalize library sizes
	objects$dge <- calcNormFactors(objects$dge)
	## Calculate common and tagwise dispersions 
	objects$dge <- estimateCommonDisp(objects$dge)
	objects$dge <- estimateTagwiseDisp(objects$dge)
	# Exact tagwise negative binomial tests
	objects$dgex <- exactTest(objects$dge)
	# Decide tests and calculate FDR
	objects$de <- decideTestsDGE(objects$dgex, p.value=p_cutoff)
	objects$tt <- topTags(objects$dgex, n=dim(objects$dgex)[1])
	###
	
	### D. Splice variants
	if (splice == "TRUE") {
		# Extract ambiguous names from exon names made unique with limma::makeUnique
		feature_names <- sub('\\.\\d+$', "", rownames(objects$dge$counts))
		# Run edgeR::spliceVariants
		objects$dgex_spl <- spliceVariants(objects$dge, feature_names, estimate.genewise.disp=TRUE)
		# Calculate FDR
		objects$tt_spl <- topTags(objects$dgex_spl, n=dim(objects$dgex_spl)[1])
	}
	###
	
	### E. Generate output (tables, plots and log file)
	# Generate comparison identifier (of format 'group1_vs_group2')
	comp <- paste(row[1], "_vs_", row[5], sep="")
	# Generate output folder name for current comparison
	out_path <- paste(out_folder, "/", comp, "/", sep="")
	# Create output folder
	dir.create(out_path, mode = "0775")
	# Save comparison information, DGE and TopTags objects
	save(objects, file=paste(out_path, out_prefix, comp, "_objects.Rdata", sep=""))
	## Write tables
	# Extract fold change, expression (counts per million), P value, FDR table
	objects$all <- objects$tt$table
	# Subset differentially expressed
	objects$diff <- objects$all[abs(objects$de) == 1, ]
	# Print table (all features) 
	write.table(objects$all, file=paste(out_path, out_prefix, comp, "_fc_cpm_p_fdr_all_features.tab", sep=""), quote=FALSE, sep="\t")
	# Print table (differentially expressed)
	write.table(objects$diff, file=paste(out_path, out_prefix, comp, "_fc_cpm_p_fdr_diff_expr.tab", sep=""), quote=FALSE, sep="\t")
	# P values and FDRs for alternative splicing (only if 'splice' is "TRUE"); all features
	if (splice == "TRUE") write.table(objects$tt_spl$table, file=paste(out_path, out_prefix, comp, "_p_fdr_splice_all_features.tab", sep=""), quote=FALSE, sep="\t")
	# Subset with P values lesser than specified value
	if (splice == "TRUE") write.table(objects$tt_spl$table[objects$tt_spl$table$PValue <= p_cutoff, ], file=paste(out_path, out_prefix, comp, "_p_fdr_splice_alt_spliced.tab", sep=""), quote=FALSE, sep="\t")
	## Multidimensional scaling plot (MDS; sample relations)
	pdf(file=paste(out_path, out_prefix, comp, "_MDS_plot.pdf", sep=""))
	plotMDS(objects$dge, labels=NULL, col=as.numeric(objects$dge$samples$group))
	dump <- dev.off()
	## Biological covariance plot (BCV; tagwise dispersion vs log2 cpm)
	pdf(file=paste(out_path, out_prefix, comp, "_BCV_plot.pdf", sep=""))
	plotBCV(objects$dge, cex=0.4)
	dump <- dev.off()
	## Smear plot (~MA; tagwise log2 FC vs log2 cpm)
	pdf(file=paste(out_path, out_prefix, comp, "_smear_plot.pdf", sep=""))
	plotSmear(objects$dgex, de.tags=rownames(objects$dge)[as.logical(objects$de)])
	abline(h = c(-1, 1), col = "blue")		# Blue lines indicate log2(FC) > 1 and < -1 
	dump <- dev.off()
	## Plot log2 fold change density functions
	# Find non-infinite fold change limits
	fc_min <- min(objects$all[is.finite(objects$all[,1]),1])
	fc_max <- max(objects$all[is.finite(objects$all[,1]),1])
	# Compute density functions of all and differentially expressed features
	dens_all <- density(objects$all[,1], from=fc_min, to=fc_max)
	dens_diff <- density(objects$diff[,1], from=fc_min, to=fc_max)
	# Plot
	pdf(file=paste(out_path, out_prefix, comp, "_fc_densities.pdf", sep=""))
	plot(dens_all, xlab="log2 fold change")
	lines(dens_diff, lty=2, col="red")
	dump <- dev.off()
	## P value ~ fold change plots
	# Only fold changes logarithmic
	pdf(file=paste(out_path, out_prefix, comp, "_p_vs_log_fc.pdf", sep=""))
	plot(objects$all[,1], objects$all[,3], xlab="log2 fold change", ylab="P value")
	points(objects$diff[,1], objects$diff[,3], pch="x", col="red")
	dump <- dev.off()
	# Log10 P values
	pdf(file=paste(out_path, out_prefix, comp, "_log_p_vs_log_fc.pdf", sep=""))
	plot(objects$all[,1], log10(objects$all[,3]), xlab="log2 fold change", ylab="log10 P value")
	points(objects$diff[,1], log10(objects$diff[,3]), pch="x", col="red")
	dump <- dev.off()	
	###

	### F. Write log
	# Generate log file name
	log <- paste(out_path, out_prefix, comp, ".log", sep="")
	sink(log)
	cat("##\n\n### Sample information:\n")
	objects$dge$samples
	cat("##\n\n### Number of unique features:\n")
	dim(objects$dge)[1]
	cat("##\n\n### Common dispersion:\n")
	objects$dge$common.dispersion
	cat("##\n\n### Pseudo/normalized library size:\n")
	objects$dge$pseudo.lib.size
	cat("##\n\n### Sample comparison:\n")
	objects$dgex$comparison
	cat("##\n\n### P value cutoff:\n")
	objects$p_cutoff
	cat("##\n\n### Differentially expressed features:\n")
	sum(summary(objects$de)[c(1,3)])
	cat("##\n\n### Downregulated, unchanged, upregulated:\n")
	summary(objects$de)
	cat("##\n\n### Quantiles all features:\n")
	quantile(objects$all[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
	cat("##\n\n### Quantiles differentially expressed features:\n")
	quantile(objects$diff[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
	cat("##\n\n### Session info:\n")
	print(sessionInfo())
	cat("##\n\n")
	sink()
})
###


QUERY = exons
SUBJECT = genes

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