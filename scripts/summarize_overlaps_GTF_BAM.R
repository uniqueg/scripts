#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Mar 22, 2013
### Modified: Jun 26, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.1
### Adapted from: n/a
### Requirements: rtracklayer, optparse, GenomicRanges
#==================#
### Description: Calculates overlaps of reads in BAM files with regions specified in a GTF file using GenomicRanges::summarizeOverlaps. Each read is counted at most once.
### Output: Count table of format: Feature name/ID TAB count
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#

#---> LOAD GETOPT LIBRARY <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-g", "--gtf"), action="store", type="character", default="", help="REQUIRED: Feature filename", metavar="gtf file"),
		make_option(c("-b", "--bam"), action="store", type="character", default="", help="REQUIRED: Read filename", metavar="bam file"),
		make_option(c("-o", "--out"), action="store", type="character", default="", help="REQUIRED: Output filename", metavar="tab file"),
		make_option(c("-i", "--index"), action="store", type="character", default=FALSE, help="Filename for BAM index (MUST be present and MUST end with '.bai', although it may be omitted from the specified name [DEFAULT: BAM filename and path]", metavar="bai file"),
		make_option(c("-r", "--genome"), action="store", type="character", default=NA, help="Reference genome identifier (e.g. 'hg19', 'mm10'...); may be unstable if input file contains sequence names that are not available in the indicated genome, therefore leave as/set to NA if problems occur [DEFAULT: NA]", metavar="string"),
		make_option(c("-m", "--mode"), action="store", type="character", default="IntersectionStrict", help="Mode of assigning reads to overlapping features or other ambiguities; allowed values: 'Union', 'IntersectionStrict' [DEFAULT], 'IntersectionNotEmpty'", metavar="string"),
		make_option(c("-s", "--ignore-strand"), action="store_true", default=FALSE, help="Ignore strand information [DEFAULT: FALSE]"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --bam [GTF] --out [FILE]", option_list = option_list, add_help_option=FALSE, description="\nDescription: Calculates overlaps of reads in BAM files with regions specified in a GTF file using GenomicRanges::summarizeOverlaps. Each read is counted at most once.")
opt <- parse_args(opt_parser)

# Die if any required arguments are missing...
if ( opt$bam == "" || opt$gtf == "" || opt$out == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
}

## Die if specified mode is not allowed
if (! ( opt$mode == "Union" || opt$mode == "IntersectionNotEmpty" || opt$mode == "IntersectionStrict") ) {
	write(paste0("[ERROR] The specified overlap mode (--mode ", opt$mode, ") is not allowed!\n\n"), stderr())	
	stop(print_help(opt_parser))
}
	
## Build complete index filename from default if no --index argument was provided
if (opt$index == FALSE) opt$index <- opt$bam else opt$index <- sub("\\.bai$", "", opt$index, ignore.case=TRUE) 

#---> LOAD MORE LIBRARIES <---#
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("Rsamtools"))

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#

#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n\n", sep="")

#---> IMPORT AND PREPARE GTF FEATURE FILE <---#
# Print status message
if ( opt$verbose ) cat("Reading GTF feature file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gtf <- import(con=opt$gtf, genome=opt$genome, format="gtf", asRangedData=FALSE) 
# Split features GRanges by gene_id --> GRangesList
grl <- split(gtf, gtf$gene_id)

#---> IMPORT BAM READ FILE <---#
# Print status message
if ( opt$verbose ) cat("Reading BAM read file '", basename(opt$bam), "' and associated BAI index file '", basename(opt$index), "'...\n", sep="")
# Read BAM file
bam <- BamFile(file=opt$bam, index=opt$index)
# Add to BamFileList object
bfl <- BamFileList(bam)

#---> FIND/COUNT OVERLAPS <---#
# Print status message
if ( opt$verbose ) cat("Counting overlaps...\n", sep="")
if ( opt$verbose ) cat("ignore_strand set to", opt$`ignore-strand`, "\n", sep="")
# Count overlaps
overlaps <- summarizeOverlaps(features=grl, reads=bfl, mode=opt$mode, ignore.strand=opt$`ignore-strand`)

#---> WRITE COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writinig output...\n")
# Write output file
write.table(assays(overlaps)$counts, file=opt$out, quote=FALSE, sep="\t")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")

#================#
#    MAIN END    #
#================#