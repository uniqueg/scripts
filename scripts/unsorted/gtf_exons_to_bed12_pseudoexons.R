#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 5, 2013
### Modified: Nov 5, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: rtracklayer, tools, optparse
#==================#
### Description: Converts the exon entries of a GTF file to a BED12+3 file with one line per transcript
### Output: A bed 12+3 file with one line per transcript and block starts/sizes for each exon
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
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="file"),
		make_option(c("-o", "--bed12"), action="store", type="character", default="", help="REQUIRED: BED12 output filename", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --bed12 [FILE]", option_list = option_list, add_help_option=FALSE, description="\nConverts the exon entries of a GTF file to a BED12+3 file with one line per transcript")
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
if ( opt$verbose ) cat("Starting '", script, "'...\n\n", sep="")

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> COMPILE LIST OF GENES <---#
# Print status message
if ( opt$verbose ) cat("Compiling transcript list...\n")
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Split exons GRanges into GRangesList by 'gene_id' 
grl <- split(gr, gr$gene_id)

#---> MAKE UNION OF EXONS <---#
# Print status message
if ( opt$verbose ) cat("Generating pseudoexons...\n")
# Make union of exons to generate pseudoexons
grl <- reduce(grl, drop.empty.ranges=FALSE, min.gapwidth=1L)
	
#---> EXPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Writing output...\n")
# Write output file / adds some track descriptor line to the top that is undesired!
export(object=grl, con=opt$bed12, format="bed15")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
