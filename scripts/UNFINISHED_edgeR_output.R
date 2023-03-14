### Improve: 
# Remove $all_genes, $diff_genes objects from edgeR_dge_disp_splice.R
# Add parameters (files1, files2, row, p_cutoff etc. to objects$ATTRIBUTES)
###

#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		19-APR-2013
### Modified:		19-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description:	Generates a number of output tables for one or more (single or pairwise) edgeR analyses conducted with edgeR_dge_disp_splice.R. Accepts a table with subsets of transcript IDs (one subset per column; transcript IDs must be a subset of those available in the original analysis) for which output is computed individually (header of each column is treated as subset name). The list object 'objects' that is used to store all relevant data is modified/updated and overwritten(!); output files are written to the output folders defined in edgeR_dge_disp_splice.R. 
### Arguments: 		1. Path containing the saved objects resulting from edgeR_dge_disp_splice.R [STRING | CHARACTER]; 2. Glob-style pattern for file selection [STRING | CHARACTER]; 3. Table containing transcript identifiers for subsetting (one column = one subset; tab-separated; requires header: subset name) [FILE | TAB]
### Output: 		Updates objects list 'object' (overwrites original file!); various tables
### Usage example:	Rscript edgeR_tables.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- "."					#args[1]
pattern <- "*_objects.Rdata"		#args[2]
subsets_file <- "subsets"			#args[3]
suppressMessages(library(edgeR))
###

### B. Import data
# Import file table from input directory
files <- dir(path=in_folder, pattern=glob2rx(pattern), full.names=TRUE, recursive=TRUE)
# Import gene subsets table
subsets_gene <- as.list(read.delim(subsets_file, stringsAsFactors=FALSE, na.strings=c("", NA)))
# Remove NA & convert to list
subsets_gene <- lapply(subsets_gene, FUN=na.omit)
# Load 'objects' object from first file
load(files[1])
## Apply over each element of gene subset list...
subsets_exon <- lapply(subsets_gene, function(set) {
			## ...and extract exon identifier subsets 
			as.character(unlist(sapply(set, function(gene) {
										rownames(objects$dge_exon[which(sub('\\.\\d+$', "", rownames(objects$dge_exon)) == gene), ])
									})))
		})
## Add all_genes and all_exons to subsets_gene and subsets_exon respectively
subsets_gene$all <- rownames(objects$dge_gene)
subsets_exon$all <- rownames(objects$dge_exon)
# Remove loaded 'objects' object
rm(objects, subsets_file, in_folder, pattern)
###

### C. Handle input files one by one
## Traverse file list row by row; returns list of ls 
ls <- lapply(files, function(filename) {
	# Load object 'objects' saved in file 'filename'
	load(filename)	
	# Add column 'dispersion' to DGEX object objects$dgex_gene_spl
	if (!is.null(objects$dgex_gene)) objects$dgex_gene_spl$table$dispersion <- objects$dgex_gene_spl$dispersion
	## Add subset lists to 'objects' objects
	objects$subsets_gene <- subsets_gene
	objects$subsets_exon <- subsets_exon
	# Create empty list 'DATA' within object 'objects'
	objects$DATA <- list()
	
	## C1. Write common dispersion to file
	if (!is.null(objects$dge_gene$common.dispersion)) writeLines(as.character(objects$dge_gene$common.dispersion), con=paste(objects$prefix, "gene", "common_dispersion", sep="_"))
	if (!is.null(objects$dge_exon$common.dispersion)) writeLines(as.character(objects$dge_exon$common.dispersion), con=paste(objects$prefix, "exon", "common_dispersion", sep="_"))
	
	## C2. Subset gene objects
	for (set in names(objects$subsets_gene)) {
		# Create empty list for subsets
		objects$DATA[[set]] <- list()
		## Create emtpy list for gene, exon and splice objects
		objects$DATA[[set]][["gene"]] <- list()
		objects$DATA[[set]][["exon"]] <- list()
		objects$DATA[[set]][["splice"]] <- list()
		## Subset gene objects
		if (!is.null(objects$dge_gene)) objects$DATA[[set]][["gene"]]$dge <- objects$dge_gene[match(objects$subsets_gene[[set]], rownames(objects$dge_gene)), ]
		if (!is.null(objects$dgex_gene)) objects$DATA[[set]][["gene"]]$dgex <- objects$dgex_gene[match(objects$subsets_gene[[set]], rownames(objects$dgex_gene$table)), ]
		if (!is.null(objects$de_gene)) objects$DATA[[set]][["gene"]]$de <- objects$de_gene[match(objects$subsets_gene[[set]], rownames(objects$dgex_gene$table)), ]
		if (!is.null(objects$tt_gene)) objects$DATA[[set]][["gene"]]$tt <- objects$tt_gene[match(objects$subsets_gene[[set]], rownames(objects$tt_gene$table)), ]
		if (!is.null(objects$dgex_gene_spl)) objects$DATA[[set]][["splice"]]$gene_dgex <- objects$dgex_gene_spl[match(objects$subsets_gene[[set]], rownames(objects$dgex_gene_spl$table)), ]
		if (!is.null(objects$tt_gene_spl)) objects$DATA[[set]][["splice"]]$gene_tt <- objects$tt_gene_spl[match(objects$subsets_gene[[set]], as.character(objects$tt_gene_spl$table[ ,1])), ]
		## Subset exon objects
		if (!is.null(objects$dge_exon)) objects$DATA[[set]][["exon"]]$dge <- objects$dge_exon[match(objects$subsets_exon[[set]], rownames(objects$dge_exon)), ]
		if (!is.null(objects$dgex_exon)) objects$DATA[[set]][["exon"]]$dgex <- objects$dgex_exon[match(objects$subsets_exon[[set]], rownames(objects$dgex_exon$table)), ]
		if (!is.null(objects$de_exon)) objects$DATA[[set]][["exon"]]$de <- objects$de_exon[match(objects$subsets_exon[[set]], rownames(objects$dgex_exon$table)), ]
		if (!is.null(objects$tt_exon)) objects$DATA[[set]][["exon"]]$tt <- objects$tt_exon[match(objects$subsets_exon[[set]], rownames(objects$tt_exon$table)), ]
		if (!is.null(objects$all_exon_spl)) objects$DATA[[set]][["splice"]]$exon_df <- data.frame(Log2FC=objects$all_exon_spl[match(objects$subsets_exon[[set]], rownames(objects$all_exon_spl)), ], row.names=objects$subsets_exon[[set]])
		## Add output file prefixes
		objects$DATA[[set]][["gene"]]$prefix <- paste(objects$prefix, set, "gene", sep="_")
		objects$DATA[[set]][["exon"]]$prefix <- paste(objects$prefix, set, "exon", sep="_")
		objects$DATA[[set]][["splice"]]$prefix <- paste(objects$prefix, set, "splice", sep="_")
		
		## C3. Write tables for gene analysis objects
		## Table (FC, CPM, stats...)
		if (!is.null(objects$DATA[[set]][["gene"]]$tt$table)) {
			write.table(objects$DATA[[set]][["gene"]]$tt$table, file=paste(objects$DATA[[set]][["gene"]]$prefix, "logFC_logCPM_P_FDR.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Counts
		if (!is.null(objects$DATA[[set]][["gene"]]$dge$counts)) {
			write.table(objects$DATA[[set]][["gene"]]$dge$counts, file=paste(objects$DATA[[set]][["gene"]]$prefix, "counts.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Pseudo-counts
		if (!is.null(objects$DATA[[set]][["gene"]]$dge$pseudo.counts)) {
			write.table(objects$DATA[[set]][["gene"]]$dge$pseudo.counts, file=paste(objects$DATA[[set]][["gene"]]$prefix, "pseudo_counts.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Tagwise dispersion
		if (!is.null(objects$DATA[[set]][["gene"]]$dge$tagwise.dispersion)) {
			temp_df <- data.frame(tagwise.dispersion=objects$DATA[[set]][["gene"]]$dge$tagwise.dispersion, row.names=rownames(objects$DATA[[set]][["gene"]]$dge))
			write.table(temp_df, file=paste(objects$DATA[[set]][["gene"]]$prefix, "tagwise_dispersion.tab", sep="_"), quote=FALSE, sep="\t")
		}		
		
		## C4. Write tables for exon analysis objects
		## Table (FC, CPM, stats...)
		if (!is.null(objects$DATA[[set]][["exon"]]$tt$table)) {
			write.table(objects$DATA[[set]][["exon"]]$tt$table, file=paste(objects$DATA[[set]][["exon"]]$prefix, "logFC_logCPM_P_FDR.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Counts
		if (!is.null(objects$DATA[[set]][["exon"]]$dge$counts)) {
			write.table(objects$DATA[[set]][["exon"]]$dge$counts, file=paste(objects$DATA[[set]][["exon"]]$prefix, "counts.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Pseudo-counts
		if (!is.null(objects$DATA[[set]][["exon"]]$dge$pseudo.counts)) {
			write.table(objects$DATA[[set]][["exon"]]$dge$pseudo.counts, file=paste(objects$DATA[[set]][["exon"]]$prefix, "pseudo_counts.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Tagwise dispersion
		if (!is.null(objects$DATA[[set]][["exon"]]$dge$tagwise.dispersion)) {
			temp_df <- data.frame(tagwise.dispersion=objects$DATA[[set]][["exon"]]$dge$tagwise.dispersion, row.names=rownames(objects$DATA[[set]][["exon"]]$dge))
			write.table(temp_df, file=paste(objects$DATA[[set]][["exon"]]$prefix, "tagwise_dispersion.tab", sep="_"), quote=FALSE, sep="\t")
		}		
		
		## C5. Write tables for splice analysis objects	
		## Splice: Table (per gene)
		if (!is.null(objects$DATA[[set]][["splice"]]$gene_tt$table)) {
			temp_df <- objects$DATA[[set]][["splice"]]$gene_tt$table[complete.cases(objects$DATA[[set]][["splice"]]$gene_tt$table[,1]), ]
			rownames(temp_df) <- temp_df[ ,1]
			temp_df <- temp_df[ ,-1]	
			write.table(temp_df, file=paste(objects$DATA[[set]][["splice"]]$prefix, "P_FDR_per_gene.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Splice: Dispersion (per gene)
		if (!is.null(objects$DATA[[set]][["splice"]]$gene_dgex$table)) {
			write.table(subset(objects$DATA[[set]][["splice"]]$gene_dgex$table, select=5, drop=FALSE), file=paste(objects$DATA[[set]][["splice"]]$prefix, "dispersion_per_gene.tab", sep="_"), quote=FALSE, sep="\t")
		}
		## Splice: Relative log2 fold change (per exon)
		if (!is.null(objects$DATA[[set]][["splice"]]$exon_df)) {
			write.table(objects$DATA[[set]][["splice"]]$exon_df, file=paste(objects$DATA[[set]][["splice"]]$prefix, "rel_logFC_per_exon.tab", sep="_"), quote=FALSE, sep="\t")
		}
	}
	
	## C4. Remove superfluous objects from 'objects' object
	objects$dge_gene <- NULL
	objects$dgex_gene <- NULL
	objects$de_gene <- NULL
	objects$tt_gene <- NULL
	objects$all_gene <- NULL		#*
	objects$diff_gene <- NULL		#*
	
	objects$dge_exon <- NULL
	objects$dgex_exon <- NULL
	objects$de_exon <- NULL
	objects$tt_exon <- NULL
	objects$all_exon <- NULL		#*
	objects$diff_exon <- NULL		#*
	
	objects$dgex_gene_spl <- NULL
	objects$tt_gene_spl <- NULL
	objects$all_gene_spl <- NULL		#*
	objects$diff_gene_spl <- NULL		#*
	objects$all_exon_spl <- NULL		
	
	# Save modified 'objects' object
	save(objects, file=paste(filename, "2", sep="."))
	# Remove 'objects' object
	rm(objects)
	# Run garbage collector
	gc()
})