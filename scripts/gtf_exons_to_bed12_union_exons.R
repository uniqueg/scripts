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


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTIONS PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="file"),
		make_option(c("-o", "--bed12"), action="store", type="character", default="", help="REQUIRED: BED12 output filename", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)

## Description
description <- "Starting from a GTF gene set annotation file, a BED12 file containing, for each gene, the ranges resulting from the union of all of its exons.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.1 (26-MAY-2014)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, version, requirements, sep="\n")

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --gtf [FILE] --bed12 [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$gtf	== ""	||	opt$bed12	== "" ) { 
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

#---> COMPILE LIST OF GENES <---#
# Print status message
if ( opt$verbose ) cat("Compiling transcript list...\n")
# Subset EXONS (discards all other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Split exons GRanges into GRangesList by 'gene_id' 
grl <- split(gr, gr$gene_id)

#---> TAKE UNION OF EXONS <---#
# Print status message
if ( opt$verbose ) cat("Generating pseudoexons...\n")
# Make union of exons to generate pseudoexons
grl <- reduce(grl, min.gapwidth=1L)
	
#---> EXPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Writing output to file '", opt$bed12, "'...\n", sep="")
# Write output file
export(object=grl, con=opt$bed12, format="bed15")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
