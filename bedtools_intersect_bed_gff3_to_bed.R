#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Dec 17, 2013
### Modified: Dec 18, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Requirements: optparse
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#

#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> BUILD DESCRIPTION <---#
version <- "Version: 1.0"
created <- "Created: 17-DEC-2013"
modified <- "Modified: 18-DEC-2013"
author <- "Author: Alexander Kanitz"
affiliation <- "Affiliation: Zavolan Group, Biozentrum, University of Basel"
description <- "Description:\nSummarizes the results of 'bedtools intersect -wo (-s) -a [ALIGNMENTS|BED] -b [FEATURES|GFF3]' between a BED read/alignment file and a GFF3 feature file. Produces a BED file indicating coordinates for each alignment and the name of the overlapped feature. Identical alignments are collapsed via the integer in the score column. Each overlapped feature is reported (i.e. alignments can be counted multiple times)."
details <- paste("", version, created, modified, author, affiliation, "" , description, sep="\n")

#---> PARSE COMMAND-LINE <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--in"), action="store", type="character", default="", help="Input filename (output of 'bedtools intersect -wo (-s) -a [ALIGNMENTS|BED] -b [FEATURES|GFF3]'; required)", metavar="file"),
		make_option(c("-o", "--out"), action="store", type="character", default="", help="Output file name (BED format; required)", metavar="file"),
		make_option(c("-t", "--type"), action="store", type="character", default="", help="Filter GFF3 entries of the indicated type (third column, e.g. 'exon').", metavar="string"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --in [FILE] --out [FILE]", option_list = option_list, add_help_option=FALSE, description=details)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE ARGUMENTS <---#
## Die if any required arguments are missing...
if 	( opt$`in`	== ""	||	opt$out	== "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
}

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#

#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> DEFINE USEFUL CHARACTER VECTORS <---#
# Print status message
if ( opt$verbose ) cat("Preparing to process input/output file formats...\n", sep="")
# Define column names
col_names <- c("rseq_r", "start_r", "stop_r", "name_r", "score_r", "str_r", "rseq_ft", "src_ft", "type_ft", "start_ft", "stop_ft", "score_ft", "str_ft", "phase_ft", "attr_ft", "overlap")
# Define output format order
out_order <- c(2,3,4,1,6,5)

#---> LOAD DATA <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", opt$`in`, "'...\n", sep="")
# Load data
df <- read.delim(opt$`in`, header=FALSE, stringsAsFactors=FALSE, col.names=col_names)

#---> SUBSETTING DATA (if applicable) <---#
if ( opt$type != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Subsetting data...\n", sep="")
	# Subset data
	df <- df[df$type == opt$type, ]
}

#---> AGGREGATE DATA <---#
# Print status message
if ( opt$verbose ) cat("Summarizing and converting data...\n", sep="")
# Aggregate data into output format
aggr <- aggregate(df$str_r, by=list(df$attr_ft, df$rseq_r, df$start_r, df$stop_r, df$str_r), length)[,out_order]

#---> WRITE OUTPUT FILE <---#
# Print status message
if ( opt$verbose ) cat("Writing output to file '", opt$out, "'...\n", sep="")
# Write data to output file
write.table(aggr, opt$out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")

#================#
#    MAIN END    #
#================#
