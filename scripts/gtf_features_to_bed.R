#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 5, 2013
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
description <- "From a GTF gene set annotation file, generates a BED6 file of the features of the specified type.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 05-NOV-2013"
version <- "Version: 1.2 (06-AUG-2014)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename (file can also be gzipped)", metavar="file"),
		make_option(c("-o", "--bed"), action="store", type="character", default="", help="REQUIRED: BED output filename", metavar="file"),
		make_option(c("-s", "--split"), action="store", type="character", default="", help="REQUIRED: Attribute to be used for the 'name' (fourth) column in the BED output file (e.g. 'transcript_ID' or 'Name'; check GFF2 file)", metavar="string"),
		make_option(c("-t", "--type"), action="store", type="character", default="", help="Select only GTF entries of the indicated type (third column, e.g. 'exon')", metavar="string"),
		make_option(c("-c", "--score"), action="store", type="numeric", default=0, help="Value to be written in 'score' (fifth) column in the BED output file (default: 0)", metavar="float"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --bed [FILE]", option_list = option_list, add_help_option=FALSE, description="\nSelects features from a GTF file and exports them as a BED file.")
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if 	( opt$gff2	== ""	||	opt$bed	== "" || opt$split == "" ) { 
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
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> LOAD LIBRARIES <---#
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
# Load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gff2), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file as GRanges object 
gr <- import(con=opt$gff2, format="gff2", asRangedData=FALSE)

#---> SELECT FEATURES <---#
# Print status message
if ( opt$verbose ) cat("Filtering features...\n")
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
if ( opt$type != "" ) {
	gr <- gr[values(gr)[["type"]] == opt$type]
	if ( ! is.null(nrow(gr)) ) stop("No entries of the indicated type in input file! Check file and/or argument to --type.\nExecution halted.\n")
}

#---> PROCESS FEATURES <---#
# Print status message
if ( opt$verbose ) cat("Renaming and sorting features...\n")
# Add name column
gr$name <- mcols(gr)[[opt$split]]
# Set scores
gr$score <- opt$score
# Sort ranges by reference sequence, strand, start and width
gr <- gr[order(width(gr))]
gr <- gr[order(start(gr))]
gr <- gr[order(strand(gr))]
gr <- gr[order(seqnames(gr))]
	
#---> EXPORT BED <---#
# Print status message
if ( opt$verbose ) cat("Writing output...\n")
# Write output file / adds some track descriptor line to the top that is undesired!
export(object=gr, con=opt$bed, format="bed")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
