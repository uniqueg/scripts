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
out_folder <- args[3]
out_prefix <- args[4]
# Load library
suppressMessages(library(edgeR))
###

### B. Import/prepare data
## B1. Read input data table and apply over each row
# Import file info table from file
input_df <- read.table(df_files, header=FALSE, sep="\t", colClasses=c("character", "factor", "character"), na.strings=c("NA",""), fill=TRUE)
## Apply over each comparison
ls <- apply(input_df, 1, function(row) {
	
	## B2. Generate list of objects to be saved for later use; add important objects
	objects <- list()
	objects$row <- row
	objects$p_cutoff <- p_cutoff
	
	## B3. Load exon count tables
	## For each comparison group, read files specified by path and pattern
	objects$files1 <- dir(path=row[2], pattern=glob2rx(row[3]), recursive=row[4], full.names=TRUE)
	objects$files2 <- dir(path=row[6], pattern=glob2rx(row[7]), recursive=row[8], full.names=TRUE)
	# Assign indicated group names and sample numbers
	objects$group <- factor(c(rep.int(row[1], length(objects$files1)), rep.int(row[5], length(objects$files2))))
	# Read files into DGEList object
	objects$dge_exon <- readDGE(c(objects$files1, objects$files2), group=objects$group, labels=basename(c(objects$files1, objects$files2)), row.names=NULL)	

	## B4. Derive gene counts from exon counts
	# Extract ambiguous names from exon names made unique with limma::makeUnique
	feature_names <- sub('\\.\\d+$', "", rownames(objects$dge$counts))
	# Assemble gene count table	
	gene_counts <- aggregate(objects$dge$counts, list(feature_names), sum)
	# Assign rownames
	rownames(gene_counts) <- gene_counts[,1]
	# Remove feature name column
	gene_counts <- gene_counts[,2:dim(gene_counts)[2]]
	# Generate DGEList object
	objects$dge_gene <- DGEList(gene_counts, objects$group)
	###
	
	### C. edgeR DGE analysis
	## Normalize library sizes
	objects$dge_exon <- calcNormFactors(objects$dge_exon)
	objects$dge_gene <- calcNormFactors(objects$dge_gene)
	## Calculate common dispersions 
	objects$dge_exon <- estimateCommonDisp(objects$dge_exon)
	objects$dge_gene <- estimateCommonDisp(objects$dge_gene)
	## Calculate tagwise dispersions
	objects$dge_exon <- estimateTagwiseDisp(objects$dge_exon)
	objects$dge_gene <- estimateTagwiseDisp(objects$dge_gene)
	## Exact tagwise negative binomial tests
	objects$dgex_exon <- exactTest(objects$dge_exon)
	objects$dgex_gene <- exactTest(objects$dge_gene)
	## Decide tests
	objects$de_exon <- decideTestsDGE(objects$dgex_exon, p.value=p_cutoff)
	objects$de_gene <- decideTestsDGE(objects$dgex_gene, p.value=p_cutoff)
	## Calculate FDR
	objects$tt_exon <- topTags(objects$dgex_exon, n=dim(objects$dgex_exon)[1])
	objects$tt_gene <- topTags(objects$dgex_gene, n=dim(objects$dgex_gene)[1]) 
	## Extract fold change, expression (counts per million), P value, FDR tables
	objects$all_exon <- objects$tt_exon$table
	objects$all_gene <- objects$tt_gene$table
	# Subset differentially expressed
	objects$diff_exon <- objects$all_exon[abs(objects$de_exon) == 1, ]
	objects$diff_gene <- objects$all_gene[abs(objects$de_gene) == 1, ]
	###
	
	### D. Splice variant analysis (exons only!)
	## D1. edgeR: Calculate P values and FDRs for alternative splicing per gene
	# Run spliceVariants()
	objects$dgex_gene_spl <- spliceVariants(objects$dge_exon, feature_names, estimate.genewise.disp=TRUE)
	# Calculate FDR
	objects$tt_gene_spl <- topTags(objects$dgex_gene_spl, n=dim(objects$dgex_gene_spl)[1])
	# Extract fold change, expression (counts per million), P value, FDR tables; Note: fold change and expression will be empty!
	objects$all_gene_spl <- objects$tt_gene_spl$table
	# Subset differentially expressed
	objects$diff_gene_spl <- objects$all_gene_spl[objects$all_gene_spl$PValue < p_cutoff, ]
	
	## D2. Calculate relative fold changes per exon
	# Create vector of subject indices for each exon feature name
	lookup_indices <- match(feature_names, rownames(objects$all_gene))
	# Calculate relative log2 fold changes between exon and gene
	rel_fc <- objects$all_exon[,1] - objects$all_gene[lookup_indices,1]
	# Convert to data frame
	objects$all_exon_spl <- data.frame("Log2FC"=rel_fc, row.names=row.names(objects$all_exon))
	###
	
	### E. Generate output (tables, plots and log file)
	## E1. Prepare folders and build reusable output file name fragments
	# Generate comparison identifier (of format 'group1_vs_group2')
	comp <- paste(row[1], "_vs_", row[5], sep="")
	# Generate output folder name for current comparison
	out_path <- paste(out_folder, "/", comp, "/", sep="")
	# Create output folder
	dir.create(out_path, mode = "0775")
	
	## E2. Save 'objects' object
	# Save comparison information, DGE and TopTags etc. objects
	save(objects, file=paste(out_path, out_prefix, comp, "_objects.Rdata", sep=""))
	
	## E3. Write tables
	# All exons 
	write.table(objects$all_exon, file=paste(out_path, out_prefix, comp, "_table_all_exons.tab", sep=""), quote=FALSE, sep="\t")
	# Differentially expressed exons
	write.table(objects$diff_exon, file=paste(out_path, out_prefix, comp, "_table_diff_expr_exons.tab", sep=""), quote=FALSE, sep="\t")
	# All genes 
	write.table(objects$all_gene, file=paste(out_path, out_prefix, comp, "_table_all_genes.tab", sep=""), quote=FALSE, sep="\t")
	# Differentially expressed genes
	write.table(objects$diff_gene, file=paste(out_path, out_prefix, comp, "_table_diff_expr_genes.tab", sep=""), quote=FALSE, sep="\t")
	# Alternatively spliced exons (fold changes only!)
	write.table(objects$all_exon_spl, file=paste(out_path, out_prefix, comp, "_table_all_spl_exons.tab", sep=""), quote=FALSE, sep="\t")
	# All genes splicing (statistics only!)
	write.table(objects$all_gene_spl, file=paste(out_path, out_prefix, comp, "_table_all_spl_genes.tab", sep=""), quote=FALSE, sep="\t")
	# Alternatively spliced genes (statistics only!)
	write.table(objects$diff_gene_spl, file=paste(out_path, out_prefix, comp, "_table_diff_spl_genes.tab", sep=""), quote=FALSE, sep="\t")
	
	## E4. Plots
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
	cat("##\n\n### P value cutoff:\n")
	objects$p_cutoff

	
	### FOR EXONS, GENES & (PARTLY) SPLICING
	cat("##\n\n### Number of unique features:\n")
	dim(objects$dge)[1]
	cat("##\n\n### Common dispersion:\n")
	objects$dge$common.dispersion
	cat("##\n\n### Pseudo/normalized library size:\n")
	objects$dge$pseudo.lib.size
	cat("##\n\n### Sample comparison:\n")
	objects$dgex$comparison
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


