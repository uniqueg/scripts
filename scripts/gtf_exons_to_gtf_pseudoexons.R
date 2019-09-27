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
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--infile"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="GTF file"),
		make_option(c("-o", "--outfile"), action="store", type="character", default="", help="REQUIRED: GTF output filename", metavar="GTF fle"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --infile [FILE] --outfile [FILE]", option_list = option_list, add_help_option=FALSE, description="\nConverts the exon entries of a GTF file to pseudoexons (i.e. union of all overlapping exons) and writes out a GTF file with one line per gene.")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$infile == "" || opt$outfile == "" ) { 
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

#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("rtracklayer"))

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$infile), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$infile, format="gtf", asRangedData=FALSE)

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
export(object=grl, con=opt$outfile, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
