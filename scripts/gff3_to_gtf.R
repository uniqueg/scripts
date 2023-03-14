#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Aug 06, 2014
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
description <- "Converts a gene set annotation file from GFF3 to GTF (GFF2) format.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 06-AUG-2014"
version <- "Version: 1.0 (06-AUG-2014)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option(c("-i", "--gff3"), action="store", type="character", default="", help="REQUIRED: Input file in GFF3 format. Can be gzipped.", metavar="file"),
		make_option(c("-o", "--gtf"), action="store", type="character", default="", help="REQUIRED: Output filename.", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die!"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die!"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [default: TRUE]."),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gff3 [FILE] --gtf [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if  ( opt$gff3 == "" || opt$gtf == "" ) { 
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
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
# Load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading GFF3 file '", basename(opt$gff3), "'...\n", sep="")
# Import GTF
gff3 <- import(opt$gff3, format="gff3", asRangedData=FALSE)

#---> EXPORT GFF3 <---#
# Print status message
if ( opt$verbose ) cat("Exporting as GTF file '", opt$gtf, "'...\n", sep="")
# Load table
export(gff3, opt$gtf, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#