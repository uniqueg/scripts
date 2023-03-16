#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		25-APR-2013
### Modified:		25-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Differential gene expression, dispersion and alternative splicing analysis using the R/Bioconductor package edgeR
### Arguments: 		1. Tab-separated file (no header!) indicating folders that contain exon count tables; FORMAT: Group 1 (1) name, (2) path, (3) pattern for file selection, (4) recursive [TRUE or FALSE]; OPTIONAL: group 2 (5) name, (6) path, (7) pattern for file selection, (8) recursive [FILE | TAB/TSV]; 2. P value cutoff (e.g. 0.05) for defining differentially expressed/alternatively spliced genes/exons [NUMERIC]; 3. Base path to the directory in which output files are written (a separate folder is produced here for each individual pairwise comparison!); 4. Prefix that is prepended to each output filename [CHARACTER] 
### Output: 		BCV, MDS and smear plots; table of all features (fold change, counts per million, P value, FDR; if argument 'splice' == "TRUE" then also P value and FDR with regards to alternative splicing evidence); table of differentially expressed (and - if argument 'splice' == "TRUE" - alternatively spliced) genes (P value < indicated cutoff); log (to STDOUT) 
### Usage example:	Rscript edgeR_dge_disp_splice.R ./input_files.tab 0.05 ./output test_run_
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
input_df <- read.table(df_files, header=FALSE, sep="\t", colClasses=c("factor", "character", "character", "logical"), na.strings="", fill=TRUE)
## Apply over each comparison
ls <- apply(input_df, 1, function(row) {
	
	## B2. Generate list of objects to be saved for later use; add important objects
	objects <- list()
	objects$row <- row
	objects$p_cutoff <- p_cutoff
	
	## B3. Load exon count tables
	## For each comparison group, read files specified by path and pattern
	objects$files1 <- dir(path=objects$row[2], pattern=glob2rx(objects$row[3]), recursive=objects$row[4], full.names=TRUE)
	if (is.na(objects$row[8])) {
		objects$files2 <- NULL
	} else {
		objects$files2 <- dir(path=objects$row[6], pattern=glob2rx(objects$row[7]), recursive=objects$row[8], full.names=TRUE)		
	}
	# Assign indicated group names and sample numbers
	objects$group <- factor(c(rep.int(objects$row[1], length(objects$files1)), rep.int(objects$row[5], length(objects$files2))))
	# Read files into DGEList object
	objects$dge_exon <- readDGE(c(objects$files1, objects$files2), group=objects$group, labels=basename(c(objects$files1, objects$files2)), row.names=NULL)	

	## B4. Derive gene counts from exon counts
	# Extract ambiguous names from exon names made unique with limma::makeUnique
	feature_names <- sub('\\.\\d+$', "", rownames(objects$dge$counts))
	# Assemble gene count table	
	gene_counts <- aggregate(objects$dge$counts, list(feature_names), sum)
	# Assign rownames
	rownames(gene_counts) <- gene_counts[ ,1]
	# Remove feature name column
	gene_counts <- gene_counts[ ,-1]
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

	## Check whether one (dispersion analysis only!) or two sample groups are given...
	if (length(levels(objects$group)) == 2) {
		
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
	}

	### E. Generate output (tables, plots and log file)
	## E1. Prepare folders and build reusable output file name fragments
	## Generate comparison identifier (of format 'group1(_vs_group2)') depending on whether single group or pairwise analysis
	if (length(levels(objects$group)) == 2) {
		comp <- paste(objects$row[1], "_vs_", objects$row[5], sep="")
	} else {
		comp <- objects$row[1]
	}
	# Generate output folder name for current comparison
	out_path <- paste(out_folder, "/", comp, "/", sep="")
	# Create output folder
	dir.create(out_path, mode = "0775")
	# Generate output file path/prefix for current comparison
	objects$prefix <- paste(out_path, out_prefix, comp, sep="")
	
	## E2. Save 'objects' object
	# Save comparison information, DGE and TopTags etc. objects
	save(objects, file=paste(objects$prefix, "objects.Rdata", sep="_"))
	
	### F. Write log
	# Generate log file name
	log <- paste(objects$prefix, "log", sep=".")
	## Write output to sink
	sink(log)
	## Print sample information and significance cutoff
	cat("##\n\n### Sample information:\n")
	objects$dge$samples
	cat("##\n\n### Sample comparison:\n")
	cat(objects$dgex$comparison, "\n")
	cat("##\n\n### P value cutoff:\n")
	objects$p_cutoff
	## FOR EXONS, GENES & (PARTLY) SPLICING
	cat("##\n\n### Number of unique features:\n")
	cat("Exons: ", dim(objects$dge_exon)[1], "\n")
	cat("Genes: ", dim(objects$dge_gene)[1], "\n")
	cat("Genes (splice analysis): ", dim(objects$dge_gene_spl)[1], "\n")
	cat("##\n\n### Common dispersion:\n")
	cat("Exons: ", objects$dge_exon$common.dispersion, "\n")
	cat("Genes: ", objects$dge_gene$common.dispersion, "\n")
	cat("##\n\n### Pseudo/normalized library size:\n")
	cat("Exons: ", objects$dge_exon$pseudo.lib.size, "\n")
	cat("Genes: ", objects$dge_gene$pseudo.lib.size, "\n")
	##  Only required if pairwise comparison...
	if (length(levels(objects$group)) == 2) {	
		# all_exon all_gene diff_exon diff_gene
		cat("##\n\n### Differentially expressed features:\n")
		cat("Exons: ", dim(objects$diff_exon)[1], "\n")
		cat("Genes: ", dim(objects$diff_gene)[1], "\n")
		cat("Genes (splice analysis): ", dim(objects$diff_gene_spl)[1], "\n")	
		cat("##\n\n### Downregulated, unchanged, upregulated:\n")
		cat("Exons:\n")
		summary(objects$de_exon)
		cat("Genes:\n")
		summary(objects$de_gene)
		# all_exon all_gene diff_exon diff_gene all_exon_spl
		cat("##\n\n### Quantiles all exons:\n")
		quantile(objects$all_exon[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
		cat("##\n\n### Quantiles all genes:\n")
		quantile(objects$all_gene[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
		cat("##\n\n### Quantiles all exons (splice analysis):\n")
		quantile(objects$all_exon_spl[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
		cat("##\n\n### Quantiles differentially expressed exons:\n")
		quantile(objects$diff_exon[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
		cat("##\n\n### Quantiles differentially expressed genes:\n")
		quantile(objects$diff_gene[,1], c(.01, .02, .05, .1, .15, .20, .25, .333, .5, .667, .75, .85, .90, .95, .98, .99))
	}
	## Print session info	
	cat("##\n\n### Session info:\n")
	print(sessionInfo())
	cat("##\n\n")
	sink()
})
###