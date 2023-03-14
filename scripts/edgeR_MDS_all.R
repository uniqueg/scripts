#######
### GENERAL:
### --------
### Author: 		Christina Herrmann
### Created: 		07-MAY-2013
### Modified:		07-MAY-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, edgeR + dependencies
### Description: 	MDS plot is created from count tables for all samples.
### Arguments: 		1. PATH to tab delimited input file containing [1] PATH to files [2] group [3] pattern for files ###			to select from PATH to files 
###			2. PATH to output files 
### Output: 		MDS-plot, heatmap
### Usage example:	Rscript edgeR_MDS_all.R path/to/in/file/ path/to/out/files/prefix
#######

### A. Pre-requisites
## Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
in_file <- args[1] 
out_prefix <- args[2] 
suppressMessages(library(edgeR))
###

### B. Import data
## Import file table from input directory
# data frame
input_df <- read.table(in_file, header=FALSE, sep="\t", 
		colClasses=c("character", "factor", "character"), 
		na.strings="", fill=TRUE)

groups <- NULL
merged_df_from_files <- function(files, colnames) {
# Description: Merges TAB files of equal length into one data frame by rownames (rownames are derived from first column, behavior potentially unstable if header not present)  
# Accepts: 1. Filenames pointing to TAB files with equal line numbers; 2. column names
# Returns: Dataframe with specified colnames
	df <- do.call(cbind,lapply(files, function(tab) read.table(tab)))
	colnames(df) <- makeUnique(rep(colnames,ncol(df)))
	groups <<- c(groups,rep(colnames,ncol(df)))
	df
}

counts_df <- do.call(cbind, apply(input_df,1,function(row){
	merged_df_from_files(dir(path=row[1], pattern=glob2rx(row[3]),
		full.names=TRUE, recursive=TRUE), row[2])
}))


# Create DGEList object

dge <- DGEList(counts_df, group=groups)

#!!! Does filtering for low count reads have to be done?
keep <- rowSums(cpm(dge)>100) >= 4
dge_no_low <- dge[keep,]
dge_no_low$samples$lib.size <- colSums(dge_no_low$counts)


### C. MDS plot

## Calculate normalisation factors
dge <- calcNormFactors(dge)
dge_no_low <- calcNormFactors(dge_no_low)

## MDS Plot for all groups
# Set colours to be used		
palette(c("yellow2","orchid","orange","tomato","tomato4","pink","peachpuff","turquoise","plum","steelblue1","royalblue4","lightblue","palegreen","brown","deeppink2","red","tomato2","yellow4","seagreen","springgreen2","darkgreen","wheat4","grey","black","lightgrey"))

# MDS: function for plotMDS, input: DGEList object of all groups, groups (vector), out_name:output prefix
MDS <- function(dge,groups,out_name){
	
	pdf(out_name, height=6, width=8)
	
	par(mar=c(5.1, 4.1, 4.1, 12), xpd=TRUE)
	plotMDS(dge, gene.selection="common",labels=NULL, cex=0.8, col=as.numeric(as.factor(groups)) , main=paste("Multidimensional Scaling Plot\n",sub("_MDS.pdf","",basename(out_name))))
	
	legend("topright", inset=c(-0.45,0), legend=levels(as.factor(groups)), fill=1:length(levels(as.factor(groups))))
	
	dev.off()
}

# Construct output file name
out_name <- paste(out_prefix, sub("input_MDS_","", basename(in_file)), "_samples", sep="")
# Plot MDS
MDS(dge,groups,paste(out_name,"_MDS.pdf",sep=""))
MDS(dge_no_low,groups, paste(out_name,"_no_low_counts","_MDS.pdf",sep=""))

