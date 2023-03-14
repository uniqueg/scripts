#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		19-APR-2013
### Modified:		19-APR-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	Analyzes the variances/dispersion of samples within groups. The following output is generated separately for each group: 1. Table of common and tagwise dispersions, 2. MDS plot (if sample number > 2), 3. BCV plot. A further MDS plot is generated for all samples, regardless of treatment groups. Finally, a heatmap of tagwise dispersion values is generated (x = groups, y = features)  
### Arguments: 		1. Tab-separated file containing no header and the following column data: I. PATH to or FILE NAME of count table(s) of format  [FEATURE | STRING | must be unique!] [TAB] [COUNT | INTEGER] with header line (output of find_overlaps_bam_within_bed.R and derived), II. Group name, III. (only required if I. does not indicate absolute file path!) Glob-style pattern indicating which files to select from indicated PATH in I. 
### Output: 		BCV, MDS and smear plots; table of differentially expressed genes (FDR < 0.5); various other count tables (see section E2.); log (to STDOUT) 
### Usage example:	Rscript edgeR_dispersion.R input_table path/to/out/files/prefix
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
df_files <- "df_test" #args[1]
prefix <- "./temp/test_"	#args[2]
# Load library
suppressMessages(library(edgeR))
###

### B. Import data
# Import file info table from file
input_df <- read.table(df_files, header=FALSE, sep="\t", colClasses=c("character", "factor", "character"), na.strings=c("NA",""), fill=TRUE)
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
pdf(paste(prefix, "MDS_plot_all_groups.pdf", sep=""), height=6, width=6)
plotMDS(dge_all, labels="o", col=as.numeric(dge_all$samples$group), main="Multidimensional scaling plot\nAll samples")
#legend("topright", inset=c(-0.2,0), legend=levels(dge_all$samples$group), fill=1:length(levels(dge_all$samples$group)))
dev.off()
###

### D. GROUPWISE DISPERSION TABLES, MDS & BCV PLOTS
# Read files of individual groups into a list of DGEList objects
dge_ls <- lapply(ls, function(group) dge <- readDGE(as.character(group$file), group=group$group, labels=basename(as.character(group$file)), columns=c(1,2), row.names=NULL))
## Apply over each DGEList object
df_ls <- lapply(dge_ls, function(dge) {
	## Plot MDS if sample number > 2
	if (dim(dge)[2] > 2) {
		pdf(paste(prefix, "MDS_plot_", levels(dge$samples$group), ".pdf", sep=""), height=6, width=6)
		plotMDS(dge, labels=NULL, col=as.numeric(dge$samples$group), main=paste("Multidimensional scaling plot\nGroup:", levels(dge$samples$group), sep=" "))
		#legend("topright", inset=c(-0.2,0), legend=levels(dge$samples$group), fill=1:length(levels(dge$samples$group)))
		dev.off()
	}
	# Normalize library sizes
	dge <- calcNormFactors(dge)
	## Calculate common and tagwise dispersions
	dge <- estimateCommonDisp(dge)
	dge <- estimateTagwiseDisp(dge)
	## Plot BCV
	pdf(paste(prefix, "BCV_plot_", levels(dge$samples$group), ".pdf", sep=""), height=6, width=6)
	plotBCV(dge)
	dev.off()
	# Construct data frame of common and tagwise dispersions
	df <- data.frame(c(dge$common.dispersion, dge$tagwise.dispersion), row.names=c("common dispersion", row.names(dge$counts)))
	# Set group names as column names
	colnames(df) <- levels(dge$samples$group)
	# Write to file
	write.table(df, file=paste(prefix, "dispersion_", levels(dge$samples$group), ".tab", sep="") , quote=FALSE, sep="\t")
	# Return data frame of 
	return(df)
})

### E. HEATMAP ALL GROUPS
# Combine dispersion value data frames for all libraries into one common matrix
mt_all <- as.matrix(do.call(cbind, df_ls))
# Set group names as column names
colnames(mt_all) <- levels(dge_all$samples$group)
## Plot heatmap
pdf(paste(prefix, "heat_map_all_groups.pdf", sep=""), height=6, width=6)
heatmap(mt_all, labRow=NA, col=heat.colors(256), scale="none")
dev.off()
###

### F. Clean-up & save session image
# Remove unused/temp variables
rm(df, df_files, input_df, ls, prefix)
# Save workspace image/session
save.image(file=paste(prefix, "_image.Rdata", sep=""))
###