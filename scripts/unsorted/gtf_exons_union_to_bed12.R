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
description <- "From a GTF gene set annotation file, generates a BED12+3 file containing the 'pseudoexons' resulting from the union of all exons of all transcripts of each gene (one row per gene).\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 05-NOV-2013"
version <- "Version: 1.2 (06-AUG-2014)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="file"),
		make_option(c("-o", "--bed12"), action="store", type="character", default="", help="REQUIRED: BED12 output filename", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --gtf [FILE] --bed12 [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if 	( opt$gtf	== ""	||	opt$bed12	== "" ) { 
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

#---> LOAD LIBRARIES <---#
# Print status message
if ( opt$verbose ) cat("Loading required libraries...\n")
# Load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> SUBSET EXONS <---#
# Print status message
if ( opt$verbose ) cat("Subsetting exons...\n")
# Subset EXONS (discards all other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Test if at least one region is returned
if ( length(gr) == 0 ) stop("No entries of type 'exon' in input file! Check file.\nExecution halted.\n")

#---> COMPILE LIST OF GENES <---#
# Print status message
if ( opt$verbose ) cat("Grouping exons by genes...\n")
# Split exons GRanges into GRangesList by 'gene_id' 
grl <- split(gr, gr$gene_id)

#---> TAKE UNION OF EXONS <---#
# Print status message
if ( opt$verbose ) cat("Merging exons...\n")
# Make union of exons to generate pseudoexons
grl <- reduce(grl, min.gapwidth=1L)
	
#---> EXPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Writing output to BED12 file '", opt$bed12, "'...\n", sep="")
# Write output file
export(object=grl, con=opt$bed12, format="bed15")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
