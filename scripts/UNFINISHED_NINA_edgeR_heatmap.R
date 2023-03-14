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
### Output: 		MDS-plot, heatmap
### Usage example:	Rscript edgeR_heatmap.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_folder <- "/home/herrmanch00/Documents/ALTSPLICE/test_heatmap/" # args[1]
pattern <- "*objects.Rdata.2"  # args[2]
out_folder <- "/home/herrmanch00/Documents/ALTSPLICE/test_heatmap" # args[3]
out_prefix <- "test_intra_"
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
	cat("Processing file '", filename , "'...\n", sep="")
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

names(tagw_ls) <- groups

# 
all_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$all)))
ref_all_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_all)))
ref_spl_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_splice)))
ref_3_gene <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$gene$refseq_3._end)))

all_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$all)))
ref_all_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_all)))
ref_spl_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_splice)))
ref_3_exon <- na.omit(do.call(cbind, lapply(tagw_ls, function(group) group$exon$refseq_3._end)))



###
	# extract relevant tables
			#	ls <- list()
			#	ls$group <- levels(objects$group)
			#	ls$exon_pseudo.counts <- objects$dge_exon$pseudo.counts
			#	colnames(ls$exon_pseudo.counts)=objects$group
			#	ls$gene_pseudo.counts <- objects$dge_gene$pseudo.counts
			#	colnames(ls$gene_pseudo.counts)=objects$group
			#	ls$exon_tagwise.dispersion <- matrix(objects$dge_exon$tagwise.dispersion, dimnames=list(rownames(objects$dge_exon$counts), makeUnique(levels(objects$group))))
			#	ls$gene_tagwise.dispersion <- matrix(objects$dge_gene$tagwise.dispersion, dimnames=list(rownames(objects$dge_gene$counts), makeUnique(levels(objects$group))))
			#			
			#	rm(objects)
			#	return(ls)		
		# Construct and return data frame
		#data.frame(file=files, group=row[2], row.names=basename(files))






#!! Set output directories
#	ls$out_path <- paste(out_folder, "/common/", sep="")


### C. MDS PLOT ALL GROUPS

		# Build data frame containg all groups from list of data frames
		#df <- do.call(rbind,ls)
		# Read all files into a DGEList object and append to list of DGEList objects
		#dge_all <- readDGE(as.character(df$file), group=df$group, labels=basename(as.character(df$file)), columns=c(1,2), row.names=NULL)

## Preparations
# Set colours to be used		
palette(c("yellow2","peachpuff", "orange","tomato4", "deeppink2","plum","turquoise","royalblue4","lightblue","palegreen", "seagreen", "springgreen2", "darkgreen", "brown", "grey", "black"))

# Combine data from all objects into one matrix for each plot
# For exon counts
exon_mds <- ls[[1]]$exon_pseudo.counts
for (i in 2:length(ls)){
	exon_mds <- cbind(exon_mds, ls[[i]]$exon_pseudo.counts)}

# For gene counts
gene_mds <- ls[[1]]$gene_pseudo.counts
for (i in 2:length(ls)){
	gene_mds <- cbind(gene_mds, ls[[i]]$gene_pseudo.counts)}

# Make vector for legend and labels
groups <- NULL
for (i in 1:length(ls)){
	groups <- c(groups, ls[[i]]$group)}

## MDS Plots
# Plot MDS for exon counts
pdf(paste(out_prefix, "MDS_plot_exons_all_groups.pdf", sep=""), height=6, width=6)

par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plotMDS(exon_mds, gene.selection="common",labels=NULL, col=as.numeric(as.factor(colnames(exon_mds))) , main="Multidimensional scaling plot\nAll samples")

legend("topright", inset=c(-0.31,0), legend=groups, fill=1:length(groups))

dev.off()

# Plot MDS for gene counts
pdf(paste(out_prefix, "MDS_plot_genes_all_groups.pdf", sep=""), height=6, width=6)

par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plotMDS(gene_mds, gene.selection="common",labels=NULL, col=as.numeric(as.factor(colnames(gene_mds))) , main="Multidimensional scaling plot\nAll samples")

legend("topright", inset=c(-0.31,0), legend=groups, fill=1:length(groups))

dev.off()
###

		#### D. GROUPWISE DISPERSION TABLES, MDS & BCV PLOTS
		## Read files of individual groups into a list of DGEList objects
		#dge_ls <- lapply(ls, function(group) dge <- readDGE(as.character(group$file), group=group$group, labels=basename(as.character(group$file)), columns=c(1,2), row.names=NULL))
		### Apply over each DGEList object
		#df_ls <- lapply(dge_ls, function(dge) {
		#	## Plot MDS if sample number > 2
		#	if (dim(dge)[2] > 2) {
		#		pdf(paste(prefix, "MDS_plot_", levels(dge$samples$group), ".pdf", sep=""), height=6, width=6)
		#		plotMDS(dge, labels=NULL, col=as.numeric(dge$samples$group), main=paste("Multidimensional scaling plot\nGroup:", levels(dge$samples$group), sep=" "))
		#		#legend("topright", inset=c(-0.2,0), legend=levels(dge$samples$group), fill=1:length(levels(dge$samples$group)))
		#		dev.off()
		#	}
		#	# Normalize library sizes
		#	dge <- calcNormFactors(dge)
		#	## Calculate common and tagwise dispersions
		#	dge <- estimateCommonDisp(dge)
		#	dge <- estimateTagwiseDisp(dge)
		#	## Plot BCV
		#	pdf(paste(prefix, "BCV_plot_", levels(dge$samples$group), ".pdf", sep=""), height=6, width=6)
		#	plotBCV(dge)
		#	dev.off()
		#	# Construct data frame of common and tagwise dispersions
		#	df <- data.frame(c(dge$common.dispersion, dge$tagwise.dispersion), row.names=c("common dispersion", row.names(dge$counts)))
		#	# Set group names as column names
		#	colnames(df) <- levels(dge$samples$group)
		#	# Write to file
		#	write.table(df, file=paste(prefix, "dispersion_", levels(dge$samples$group), ".tab", sep="") , quote=FALSE, sep="\t")
		#	# Return data frame of 
		#	return(df)
		#})

### E. HEATMAP ALL GROUPS
# Combine dispersion value matrices for all libraries into one common matrix
#mt_all <- as.matrix(do.call(cbind, df_ls))

exon_disp_all <- ls[[1]]$exon_tagwise.dispersion
for (i in 2:length(ls)){
	exon_disp_all <- cbind(exon_disp_all, ls[[i]]$exon_tagwise.dispersion)}

# Set group names as column names
#colnames(mt_all) <- levels(dge_all$samples$group)

## Plot heatmap
library(gplots)
def.par <- par(no.readonly = TRUE) # save default, for resetting...

pdf(paste(out_prefix, "heat_map_all_groups.pdf", sep=""), height=20, width=6)

par(mar = c(0,0,0,0))
heatmap.2 (exon_disp_all,
		
		# dendrogram control
		Rowv = TRUE,
		Colv = TRUE,
		distfun = dist,
		hclustfun = hclust,
		dendrogram = c("both"),
		symm = FALSE,
		
		# data scaling
		scale = c("none"), # Default is "row" if x is not symm
		na.rm=TRUE,
		
		# image plot
		# revC = identical("Rowv"),
		# add.expr,
		
		# mapping data to colors
		# breaks=c(seq(from = 0, to = 0.05, by = 0.001),seq(from = 0.05, to = 0.1, by = 0.005),seq(from = 0.1, to = 2, by = 0.5)), # numeric vector specifying breakpoints, or integer for number of breakpoints
		symbreaks= TRUE, # min(x < 0, na.rm=TRUE) || scale!="none",
		
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
		trace=c("column"),
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
		labRow =rownames(mt_all),
		labCol =colnames(mt_all),
		
		# color key + density info
		key = TRUE,
		keysize = 1.5,
		density.info=c("histogram"), # also possible:"density","none"),
		denscol="cyan",
		symkey = TRUE , #min(x < 0, na.rm=TRUE) || symbreaks,
		densadj = 0.25,
		
		# plot labels
		main = "All Dispersions",
		xlab = NULL,
		ylab = NULL,
		
		# plot layout
		#lmat=rbind( c(0, 3), c(2,1), c(0,4) ),
		#lhei=c(0.25, 4, 0.25 ),
		lmat = matrix(c(4,4,0,3,2,1,2,1), 4, 2, byrow = TRUE),
		lhei = c(2,1,10,10),
		lwid = NULL,
		
		# extras
		
)


dev.off()
par(def.par)  #- reset to default
###

### F. Clean-up & save session image
# Remove unused/temp variables
rm(df, df_files, input_df, ls, prefix)
# Save workspace image/session
save.image(file=paste(prefix, "_image.Rdata", sep=""))
###
