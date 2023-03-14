#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 8, 2013
### Modified: Nov 8, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Requirements: rtracklayer, tools, optparse
#==================#
### Description: From a GTF file containing exons, writes out a GTF file of exons only (i.e. filters out other types of entries) and/or a GTF file of introns; script was written/optimized for gene set annotation files obtained from ENSEMBL; other GTF files may not be compatible (in particular, an attribute field "transcript_id" is required!".
### Output: A GTF file of exons and/or introns
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
		make_option(c("-g", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="file"),
		make_option(c("-e", "--exon"), action="store", type="character", default="", help="Filename for exon GTF output file (one of --exon and --intron is required!)", metavar="file"),
		make_option(c("-i", "--intron"), action="store", type="character", default="", help="Filename for intron GTF output file (one of --exon and --intron is required!)", metavar="file"),		
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --bed12 [FILE]", option_list = option_list, add_help_option=FALSE, description="\nFrom a GTF file containing exons, writes out a GTF file of exons only (i.e. filters out other types of entries) and/or a GTF file of introns; script was written/optimized for gene set annotation files obtained from ENSEMBL; other GTF files may not be compatible (in particular, an attribute field 'transcript_id' is required!).")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$gtf	== "" || ( opt$exon	== "" && opt$intron == "" ) ) { 
	write("[ERROR] Required argument(s) missing!\n", stderr())	
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
# Split exons GRanges into GRangesList by 'transcript_id' 
exons_grl <- split(gr, gr$transcript_id)

#---> EXPORT EXONS <---#
if ( opt$exon != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Writing exons...\n")
	# Convert GRangesList back to GRanges
	exons_gr <- unlist(exons_grl)
	# Remove metadata
	mcols(exons_gr) <- NULL
	# Set name field in metadata
	exons_gr$name <- names(exons_gr)
	# Remove range names to avoid duplicates (will throw error!)
	names(exons_gr) <- NULL
	# Write output file / adds some track descriptor line to the top that is undesired!
	export(object=exons_gr, con=opt$exon, format="bed")	
}

#---> DETERMINE INTRONS <---#
if (opt$intron != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Determining intron ranges...\n")
	# Create a list of transcripts with one range each, starting from the first base of the first exon to the last base of the last
	genes_grl <- range(exons_grl)
	# Calculate the intron ranges by taking the differences between the one-range and the exon transcript lists
	introns_grl <- psetdiff(genes_grl, exons_grl)

	#---> EXPORT INTRONS <---#
	# Print status message
	if ( opt$verbose ) cat("Writing introns...\n")
	# Convert GRangesList back to GRanges
	introns_gr <- unlist(introns_grl)
	# Set name field in metadata
	introns_gr$name <- names(introns_gr)
	# Remove range names to avoid duplicates (will throw error!)
	names(introns_gr) <- NULL
	# Write output file / adds some track descriptor line to the top that is undesired!
	export(object=introns_gr, con=opt$intron, format="bed")
}

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#