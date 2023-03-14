#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		19-APR-2013
### Modified:		19-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	From objects created by edgeR_???? script,a MDS plot is generated for all samples, regardless of treatment groups. Finally, a heatmap of tagwise dispersion values is generated (x = groups, y = features)  
### Arguments: 		1. Tab-separated file containing no header and the following column data: I. PATH to or FILE NAME of count table(s) of format  [FEATURE | STRING | must be unique!] [TAB] [COUNT | INTEGER] with header line (output of find_overlaps_bam_within_bed.R and derived), II. Group name, III. (only required if I. does not indicate absolute file path!) Glob-style pattern indicating which files to select from indicated PATH in I. 
### Output: 		Heatmap
### Usage example:	Rscript edgeR_heatmap.R in_folder pattern path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- args[1]
pattern <- args[2]
out_prefix <- args[3]
suppressMessages(library(edgeR))
###

### B. Import data
# Import file table from input directory
files <- dir(path=in_folder, pattern=glob2rx(pattern), full.names=TRUE, recursive=TRUE)

# vector for group names of loaded objects, used for naming list elements later
groups <- NULL

## Traverse file list row by row; returns list tagw_ls, each entry=1 tissue, contains list of matrices for plots
tagw_ls <- lapply(files, function(filename) {
	# load file for current tissue
	load(filename)
	# Print status message
	cat("Processing file '", filename , "'...\n", sep="")
	# Retrieve group (or group comparison)of current object 
	group <- ifelse(length(objects$row)>4,paste(objects$row[1], objects$row[5], sep="_vs_"),objects$row[1])
	
	# genes: extract relevant data from object for current tissue
	gene <- lapply(objects$DATA, function (set) {
				matrix(set$gene$dge$tagwise.dispersion, dimnames=list(rownames(set$gene$dge), group))	
			})
	# exons: extract relevant data from object for current tissue
	exon <- lapply(objects$DATA, function (set) {
				matrix(set$exon$dge$tagwise.dispersion, dimnames=list(rownames(set$exon$dge), group))	
			})
	groups <<- c(groups, group)
	# Print status message
	cat("Processed file '", filename , "'!\n", sep="")
	return(list(gene=gene,exon=exon))
})

# Set names of tagw_ls entries to corresponding group
names(tagw_ls) <- groups

# Combine tagwise dispersions for all groups for all gene subsets into 
# one matrix for each gene subset
#all_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$all)))
ref_all_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_all)))
ref_spl_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_splice)))
ref_3_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_3._end)))

#all_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$all)))
ref_all_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_all)))
ref_spl_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_splice)))
ref_3_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_3._end)))



### C. HEATMAP ALL GROUPS

# heat: heatmap function. Input: matrix of interest, output prefix
heat <- function(mt,out){
	# Construct output file name
	out_name <- paste(out_prefix, out, "_heatmap.pdf", sep="")
	
	# Save Data matrix to table
	write.table(mt, file=paste(out_prefix, out, "_tagwise_dispersion_matrix.tab", sep=""), quote=FALSE, sep="\t")
	
	## Plot heatmap
	library(gplots)
	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	
	pdf(out_name, height=20, width=6)
	
	par(mar = c(0,0,0,0))
	heatmap.2 (mt,
			
			# dendrogram control
			Rowv = TRUE,
			Colv = FALSE,
			distfun = dist,
			hclustfun = hclust,
			dendrogram = c("none"),
			symm = FALSE,
			
			# data scaling
			scale = c("none"), # Default is "row" if x is not symm
			na.rm=TRUE,
			
			# image plot
			# revC = identical("Rowv"),
			# add.expr,
			
			# mapping data to colors
			# breaks=c(seq(from = 0, to = 0.05, by = 0.001),seq(from = 0.05, to = 0.1, by = 0.005),seq(from = 0.1, to = 2, by = 0.5)), # numeric vector specifying breakpoints, or integer for number of breakpoints
			symbreaks= F, # min(x < 0, na.rm=TRUE) || scale!="none",
			
			# colors
			col="heat.colors",
			
			# block separation
	#		colsep,
	#		rowsep,
	#		sepcolor="white",
	#		sepwidth=c(0.05,0.05),
			
			# cell labeling
			# cellnote, matrix of char, symbols to place in each cell, e.g. p-values
	#		notecex=1.0,
	#		notecol="cyan",
	#		na.color=par("bg"),
			
			# level trace
			trace=c("none"),
			tracecol="cyan",
			hline=NULL,  # median(breaks),
			vline=NULL, # median(breaks),
			#linecol="cyan",
			
			# Row/Column Labeling
			margins = c(5, 5),
			# ? ColSideColors,
			# ? RowSideColors,
			cexRow = 0.2, #+ 1/log10(1000), # nr: number of rows: 1000
			cexCol = 2, # + 1/log10(9), # nc: number of cols: 9
			labRow =F,
			labCol =colnames(mt),
			
			# color key + density info
			key = TRUE,
			keysize = 1.2,
			density.info=c("density"), # also possible:"density","none"),
			denscol="black",
			symkey = F , #min(x < 0, na.rm=TRUE) || symbreaks,
			densadj = 0.25,
			
			# plot labels
			main = "All Dispersions",
			xlab = NULL,
			ylab = NULL,
			
			# plot layout
			#lmat=rbind( c(0, 3), c(2,1), c(0,4) ),
			#lhei=c(0.25, 4, 0.25 ),
			lmat = matrix(c(4,4,2,3,0,1,0,1,0,0), 5, 2, byrow = TRUE),
			lhei = c(2,1,10,10,2),
			lwid = c(1,10),
			
			# extras
			
	)
	
	
	dev.off()
	par(def.par)  #- reset to default
	###
}

# Call heat for all matrices obtained in ###B.
# Heatmap will be created and matrices will be saved

#heat(all_gene, "all_gene")
heat(ref_all_gene, "ref_all_gene")
heat(ref_spl_gene, "ref_spl_gene")
heat(ref_3_gene, "ref_3_gene")

#heat(all_exon, "all_exon")
heat(ref_all_exon, "ref_all_exon")
heat(ref_spl_exon, "ref_spl_exon")
heat(ref_3_exon, "ref_3_exon")


