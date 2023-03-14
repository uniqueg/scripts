#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 4, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#    HEADER END    #
#==================#

#===================#
#   OPTIONS START   #
#===================#
#---> LOAD OPTIONS PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Generates a table of RPKM-normalized read counts from a table of raw read counts and a table of gene lengths.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 04-NOV-2014"
version <- "Version: 1.0 (04-NOV-2014)"
requirements <- "Requires: optparse"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option(c("-i", "--raw"), action="store", type="character", default="", help="REQUIRED: Raw gene counts in TAB format.", metavar="file"),
		make_option(c("-l", "--lengths"), action="store", type="character", default="", help="REQUIRED: Gene lengths file in TAB format.", metavar="file"),
		make_option(c("-o", "--rpkm"), action="store", type="character", default="", help="REQUIRED: Normalized (RPKM) counts in TAB format (output filename).", metavar="file"),
		make_option(c("-c", "--comprehensive"), action="store_true", default=FALSE, help="Print comprehensive output table including CPM (counts per million reads) and gene lengths."),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die."),
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages.")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --raw [FILE] --lengths [FILE] --rpkm [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if 	( opt$raw	== ""	||	opt$lengths	== "" || opt$rpkm == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
#===================#
#    OPTIONS END    #
#===================#

#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n\n", sep="")

#---> READ TABLE OF RAW COUNTS <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$raw), "'...\n", sep="")
# Read table
raw <- read.table(opt$raw, sep="\t", header=FALSE, row.names=1, col.names=c("id", "raw_count"))

#---> READ TABLE OF GENE LENGTHS <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$lengths), "'...\n", sep="")
# Read table
len <- read.table(opt$lengths, sep="\t", header=FALSE, row.names=1, col.names=c("id", "length"))

#---> NORMALIZE BY LIBRARY SIZE <---#
# Print status message
if ( opt$verbose ) cat("Compute counts per million reads (CPM)...\n", sep="")
# Compute CPM
cpm <- cbind(raw, cpm=raw$raw_count / sum(raw$raw_count) * 1000000)

#---> NORMALIZE BY GENE LENGTHS <---#
# Print status message
if ( opt$verbose ) cat("Compute reads per kilobase per million reads (CPM)...\n", sep="")
# Merge CPM table with gene lengths
cpm <- merge(cpm, len, by=0)
# Warn if no size was found for every gene with raw counts
if (nrow(cpm) < nrow(raw)) warning("[WARNING] Gene lengths were not supplied for all genes in the raw count table. Some genes may be missing from the output.")
# Compute RPKM
rpkm <- cbind(cpm, rpkm=cpm$cpm / cpm$length * 1000)

#---> SUBSET TABLE <---#
if ( ! opt$comprehensive ) {
	# Print status message
	if ( opt$verbose ) cat("Subsetting output table...\n", sep="")
	# Subset table
	rpkm <- rpkm[ ,c("Row.names", "rpkm")]
}

#---> WRITE OUTPUT <---#
# Print status message
if ( opt$verbose ) cat("Writing output to TAB file '", opt$rpkm, "'...\n", sep="")
# Write output file
write.table(rpkm, file=opt$rpkm, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#