#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 11, 2013
### Modified: Nov 11, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: rtracklayer, optparse
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--input"), action="store", type="character", default="", help="REQUIRED: Input file in BED format", metavar="file"),
		make_option(c("-o", "--output"), action="store", type="character", default="", help="REQUIRED: Output file in BED format", metavar="file"),
		make_option(c("-c", "--chain"), action="store", type="character", default="", help="REQUIRED: Appropriate liftOver CHAIN file (obtain from UCSC)", metavar="file"),
		make_option(c("-l", "--label"), action="store", type="character", default="", help="REQUIRED: Arbitraty label for the lifted coordinates (e.g. mm10, dm3 etc.)", metavar="string"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gff2 [FILE] --bed [FILE]", option_list = option_list, add_help_option=FALSE, description="\nLifts over genomic coordinates from one system to another using a UCSC chain file. Accepts and returns files in BED6 format. One input range may map to multiple output ranges (range names are made unique via appending a serial number separated by a dot).")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$input	== ""	||	opt$output	== "" || opt$chain == "" ||	opt$label	== "" ) { 
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

#---> IMPORT BED <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$input), "'...\n", sep="")
# Use rtracklayer::import method to import BED input file to GRanges object 
original <- import.bed(opt$input, asRangedData=FALSE)

#---> IMPORT CHAIN <---#
# Print status message
if ( opt$verbose ) cat("Reading chain file '", basename(opt$chain), "'...\n", sep="")
# Use rtracklayer::import method to import CHAIN file
chain <- import.chain(opt$chain)

#---> LIFT COORDINATES <---#
# Print status message
if ( opt$verbose ) cat("Lifting over coordinates...\n")
# Use rtracklayer::liftOver method to lift over coordinates to new system
lifted <- liftOver(original, chain)

#---> REFORMAT OUTPUT <---#
# Print status message
if ( opt$verbose ) cat("Preparing export...\n")
# Unlist GRangesList object
lifted <- unlist(lifted)
# Generate unique names including specified label and a serial number if more than one fragment
lifted$name <- make.unique(paste(lifted$name, opt$label, sep="_"))
# Remove GRanges names to avoid duplicate rownames (causes problems during export...)
names(lifted) <- NULL

#---> EXPORT OUTPUT <---#
# Print status message
if ( opt$verbose ) cat("Writing results to file 'opt$output'...\n")
# Use rtracklayer::export method to export GRanges object as BED file
export.bed(lifted, opt$output)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
