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
### Requirements: rtracklayer, optparse
#==================#
### Description: Converts the exon entries of a GFF2 file to a BED file with one line per transcript
### Output: A BED file with one line per transcript and block starts/sizes for each exon
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
		make_option(c("-i", "--gff2"), action="store", type="character", default="", help="REQUIRED: GFF2 input filename", metavar="file"),
		make_option(c("-o", "--bed"), action="store", type="character", default="", help="REQUIRED: BED output filename", metavar="file"),
		make_option(c("-s", "--split"), action="store", type="character", default="", help="REQUIRED: Attribute to be used for the 'name' (fourth) column in the BED output file (e.g. 'transcript_ID' or 'Name'; check GFF2 file)", metavar="string"),
		make_option(c("-t", "--type"), action="store", type="character", default="", help="Filter GFF2 entries with the indicated type (third column, e.g. 'exon')", metavar="string"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gff2 [FILE] --bed [FILE]", option_list = option_list, add_help_option=FALSE, description="\nConverts the exon entries of a GFF2 file to a BED file with one line per transcript")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$gff2	== ""	||	opt$bed	== "" || opt$split == "" ) { 
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

#---> IMPORT GFF2 <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gff2), "'...\n", sep="")
# Use rtracklayer::import method to import GFF2 file to GRanges object 
gr <- import(con=opt$gff2, format="gff2", asRangedData=FALSE)

#---> COMPILE LIST OF TRANSCRIPTS <---#
# Print status message
if ( opt$verbose ) cat("Compiling transcript list...\n")
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
if ( opt$type != "" ) {
	gr <- gr[values(gr)[["type"]] == opt$type]
	if ( ! is.null(nrow(gr)) ) stop("No entries of the indicated type in input file! Check file and/or argument to --type.\nExecution halted.\n")
}
# Split exons GRanges into GRangesList by indicated identifier
grl <- split(gr, mcols(gr)[[opt$split]])
	
#---> EXPORT BED <---#
# Print status message
if ( opt$verbose ) cat("Writing output...\n")
# Write output file / adds some track descriptor line to the top that is undesired!
export(object=grl, con=opt$bed, format="bed")
# Read output file
df <- read.table(file=opt$bed, header=FALSE, sep="\t")
# Subset first six columns
df <- df[,1:6]
# Write output file
write.table(df, file=opt$bed, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#