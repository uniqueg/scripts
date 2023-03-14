#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 5, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Adapted from: n/a
#==================#
#    HEADER END    #
#==================#


#=================#
#   IDEAS START   #
#=================#
### Comma-separated list of integers for column selection for each of the input files
### If counts and length files contain feature names (e.g. transcript ID), then allow matching based on them (equal lines not necessary, only that all factors in counts files are present in len file)
### If counts file contains a column for library names, then match based on them (equal lines not necessary, only that all factors in counts files are present in lib file)
#=================#
#    IDEAS END    #
#=================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option("--feature-counts", action="store", type="character", default="", help="REQUIRED: File with counts for each feature", metavar="file"),
		make_option("--library-sizes", action="store", type="character", default="", help="REQUIRED: File with library sizes for each feature", metavar="file"),
		make_option("--prefix", action="store", type="character", default="", help="REQUIRED: Output file prefix", metavar="file"),
		make_option("--feature-lengths", action="store", type="character", default="", help="File with lengths of each feature (in nt)", metavar="string"),
		make_option("--rownames", action="store_true", default=FALSE, help="Indicate if input files contain row names (separated from values by TAB)"), 
		make_option("--header", action="store_true", default=FALSE, help="Indicate if input files contain headers"),
		make_option("--force-cpm", action="store_true", default=FALSE, help="Output a file each for CPM and RPKM counts (Default: Only write RPKM if argument to --feature-lengths is provided)"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

## Description
description <- "Calculates CPM and/or RPKM values from feature counts, library sizes and lengths.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0 (07-MAR-2014)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --feature-counts [FILE] --library-sizes [FILE] --prefix [STRING]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$`feature-counts`	== ""	||	opt$`library-sizes`	== "" || opt$prefix == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
}

## Set defaults
if (opt$rownames) col = 2 else col = 1

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> CREATE FILE CONTAINER <---#
files <- list()

#---> IMPORT FEATURE COUNTS <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$`feature-counts`), "'...\n", sep="")
# Load table
files$cts <- read.table(opt$`feature-counts`, header=opt$header, stringsAsFactors=FALSE)

#---> IMPORT LIBRARY SIZES <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$`library-sizes`), "'...\n", sep="")
# Load table
files$lib <- read.table(opt$`library-sizes`, header=opt$header, stringsAsFactors=FALSE)

#---> IMPORT FEATURE LENGTHS <---#
if ( opt$`feature-lengths` != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Reading input file '", basename(opt$`feature-lengths`), "'...\n", sep="")
	# Load table
	files$len <- read.table(opt$`feature-lengths`, header=opt$header, stringsAsFactors=FALSE)
}

#---> VALIDATE INPUT FILES <---#
# Print status message
if ( opt$verbose ) cat("Validating input files...\n", sep="")
# Check for equal row numbers
if ( diff(range(lapply(files, nrow))) ) {
	write("[ERROR] Input files have differing row numbers.\nExecution aborted!\n", stderr())	
	stop(print_help(opt_parser))
}
# Check if enough columns are present in all files
if ( any(sapply(files, ncol) < col) ) {
	write("[ERROR] No values detected. Check presence of values and/or rownames.\nExecution aborted!\n", stderr())	
	stop(print_help(opt_parser))
}
# Check if all values are integers
if ( all(sapply(files, function(file) all(class(file[,col]) == "integer")) == FALSE) ) {
	write("[ERROR] Not all values in first (non-rowname) column are integers.\nExecution aborted!\n", stderr())	
	stop(print_help(opt_parser))
}

#---> CALCULATE CPM <---#
# Print status message
if ( opt$verbose ) cat("Calculating CPM values...\n", sep="")
# Calculate
cpm <- files$cts[,col] / files$lib[,col] * 1000000
	
#---> CALCULATE RPKM <---#
if ( opt$`feature-lengths` != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Calculating RPKM values...\n", sep="")
	# Calculate
	rpkm <- cpm / files$len[,col] * 1000
}

#---> WRITE OUTPUT FILE(S) <---#
# Print status message
if ( opt$verbose ) cat("Write output files...\n", sep="")
# 
if ( opt$`feature-lengths` == "" || opt$`force-cpm` ) write.table(cpm, file=paste0(opt$prefix, "cpm.tab"), quote=FALSE, col.names=FALSE, row.names=FALSE)
if ( opt$`feature-lengths` != "" ) write.table(rpkm, file=paste0(opt$prefix, "rpkm.tab"), quote=FALSE, col.names=FALSE, row.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#