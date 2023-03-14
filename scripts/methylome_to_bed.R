#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		20-FEB-2013
### Modified:		20-FEB-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	N/A
### Description: 	Convert methylome data sets from Lister et al., Nature, 2009 (Human DNA methylomes at base resolution show widespread epigenomic differences) to BED files; only CG dinucleotides are considered
### Arguments: 		1. directory containing input files; 2. directory for writing output files 
### Output: 		BED files indicating chromosome, start/end position, name (generated from coordinates), score (ratios of methylated Cs and total Cs), strand 
### Usage:			Rscript ./methylome_to_bed.R /path/to/methylome/data/files /path/to/output/files
#######

# Pass command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Suppress scientific notation
options(scipen=999)

# Create file list
files_ls <- list.files(args[1], "^mc_")

# Traverse through each file in file list
lapply(files_ls, function (file) {
	# Load data
	df <- read.table(paste(args[1], "/", file, sep=""), header=1, as.is=TRUE)
	# Exclude all but CG dinucleotides
	df <- df[df$class == "CG", ]
	# Generate chromosome name string
	df$chr <- paste("chr", df$assembly, sep="")
	# Compute end position
	df$end <- df$position + 1
	# Generate name string from coordinates
	df$name <- paste(df$chr, ":", df$position, "-", df$end, ":", df$strand, sep="")
	# Calculate methylated frequency by dividing methylated C by total base calls
	df$score <- df$mc / df$h
	# Generate .bed file table
	bed <- cbind(df$chr, df$position, df$end, df$name, df$score, df$strand)
	# Write to file
	write.table(bed, paste(args[2], "/", file, "_hg18.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	## Write screen output
	writeLines(paste("File '", args[1], "/", file, "' processed.", sep=""))
	writeLines(paste("Output written to file '", args[2], "/", file, "_hg18.bed'.", sep=""))
})

# Reset scientific notation option
options(scipen=0)

## Write to log
writeLines("Input file directory:")
args[1]
writeLines("Files found:")
files_ls
writeLines("Output file directory:")
args[2]
writeLines("\nSession info")
print.default(sessionInfo())
###
