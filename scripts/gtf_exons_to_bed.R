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
description <- "From a GTF gene set annotation file, generates a BED6 file of the features of type 'exon'.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 05-NOV-2013"
version <- "Version: 1.1 (06-AUG-2014)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename.", metavar="file"),
		make_option(c("-o", "--bed"), action="store", type="character", default="", help="REQUIRED: BED output filename.", metavar="file"),
		make_option(c("-n", "--name"), action="store", type="character", default="transcript_id", help="Attribute to be used for the 'name' (fourth) column in the BED output file (check GTF file for available options; default: 'transcript_id').", metavar="string"),
		make_option(c("-s", "--score"), action="store", type="numeric", default=0, help="Score that should be set in each row of the output BED file (default: 0).", metavar="num"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die."),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]."),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --gtf [FILE] --bed [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if      ( opt$gtf       == ""   ||      opt$bed == "" || opt$name == "" ) {
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
# Print status message
if ( opt$verbose ) cat("Loading required libraries...\n")
# Load packages
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> FILTER EXONS <---#
# Print status message
if ( opt$verbose ) cat("Filtering exons...\n")
# Subset exonS (discard all other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Test if at least one region is returned
if ( ! is.null(nrow(gr)) ) stop("No entries of type 'exon' in input file! Check file.\nExecution halted.\n")

#---> GROUP ENTRIES BY NAME FIELD <---"
# Print status message
if ( opt$verbose ) cat("Grouping exons by the '", opt$name, "' attribute...\n", sep="")
# Test if specified argument to --name exists
if (! opt$name %in% names(mcols(gr)) ) stop("'", opt$name, "' is not an attribute of the input GTF file. Check the file and spelling.\nExecution aborted.")
# Split exons GRanges into GRangesList by indicated identifier
grl <- split(gr, mcols(gr)[[opt$name]])
# Unlist
gr <- unlist(grl)

#---> ADD NAME METADATA COLUMN <---#
# Print status message
if ( opt$verbose ) cat("Naming exons...\n")
# Add 'name' column
gr$name <- mcols(gr)[[opt$name]]
# Resetting range names
names(gr) <- NULL

#---> SET SCORES <---#
# Print status message
if ( opt$verbose ) cat("Setting scores...\n")
# Add 'name' column
gr$score <- opt$score

#---> EXPORT BED <---#
# Print status message
if ( opt$verbose ) cat("Writing exons to BED file '", opt$bed , "'...\n", sep="")
# Export BED file
export(object=gr, con=opt$bed, format="bed")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
