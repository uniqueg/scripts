#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Jan 28, 2014
### Modified: Jan 28, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: rtracklayer, optparse
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER LIBRARY <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-g", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF file (tested on ENSEMBL gene set annotation files)", metavar="file"),
		make_option(c("-p", "--prefix"), action="store", type="character", default="", help="REQUIRED: Prefix for output files", metavar="string"),
		make_option(c("-n", "--merge"), action="store_true", default=FALSE, help="Merge ranges in output file; names will then be derived from coordinates of merged ranges rather than transcript IDs!"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

## Parse command-line arguments
description <- "\nBased on a GTF file in ENSEMBL gene set annotation format writes out the coordinates of CDS and 5'- and 3'-untranslated regions in BED format."
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --input [FILE] --prefix [FILE]", option_list = option_list, add_help_option=FALSE, description=description)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$gtf	== "" || opt$prefix	== "" ) { 
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

#---> LOAD LIBRARIES <---#
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
suppressPackageStartupMessages(library("rtracklayer"))

#---> READ DATA <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gtf), "' (may take long!)...\n", sep="")
# Import GTF file as GRanges object
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> EXTRACT PROTEIN-CODING TRANSCRIPTS <---#
# Print status message
if ( opt$verbose ) cat("Extracting protein-coding transcripts...\n", sep="")
# Extract identifiers of transcripts containing a CDS
trx_id <- unique(subset(gr, type == "CDS")$transcript_id)
# Subset all ranges of coding transcripts
cdng_trx <- gr[gr$transcript_id %in% trx_id,]
## Subset exons from coding transcripts and split GRanges object by transcript ID
ex <- subset(cdng_trx, type == "exon")
ex_l <- split(ex, ex$transcript_id)
## Subset CDSs from coding transcripts and split GRanges object by transcript ID
cds <- subset(cdng_trx, type == "CDS")
names(cds) <- cds$transcript_id
cds_l <- split(cds, cds$transcript_id)

#---> DERIVE UTRS AND UTR ANNOTATIONS <---#
# Print status message
if ( opt$verbose ) cat("Deriving and annotating untranslated regions (may take long!)...\n", sep="")
# Calculate UTRs (difference of exons and CDSs)
utr_l <- psetdiff(ex_l, cds_l)
## Iterate over UTRs and CDSs...
utr_types <- mapply(function(utr,cds) {
			## Extract genomic start and end coordinates of CDS
			cds_start <- min(start(cds))
			cds_end <- max(end(cds))
			# Initialize container for UTR annotations	
			type <- NULL
			## Verify that a UTR range exists...
			if ( length(utr) > 0 ) {
				## Iterate over each UTR range...
				for ( j in 1:length(utr) ) {
					## Depending on the orientation of the UTR range relative to the CDS, add "5UTR" or "3UTR" to the annotation container
					if ( ( start(utr[j]) < cds_start && as.character(strand(utr[j])) == "+" ) || ( start(utr[j]) > cds_end && as.character(strand(utr[j])) == "-" ) ) { 
						type <- c(type, "5UTR")
					} else {
						type <- c(type, "3UTR")
					}
				}
			}
			# Annotate UTR ranges of transcript	
			return(type)
}, utr_l, cds_l)
## Unlist and combine
utr <- unlist(utr_l)
utr$type <- unlist(utr_types)
## Split by orientation
utr5 <- subset(utr, type == "5UTR")
utr3 <- subset(utr, type == "3UTR")

#---> PREPARE DATA FOR EXPORT <---#
# Print status message
if ( opt$verbose ) cat("Preparing data for export...\n", sep="")
# Remove all metadata
mcols(cds) <- NULL
mcols(utr5) <- NULL
mcols(utr3) <- NULL
## Prepare data depending on whether ranges are supposed to be merged
if ( opt$merge ) {
	## Merge ranges
	cds <- reduce(cds)
	utr5 <- reduce(utr5)
	utr3 <- reduce(utr3)
	## Set names
	cds$name <- paste(as.character(seqnames(cds)), start(cds), end(cds), as.character(strand(cds)), sep="_")
	utr5$name <- paste(as.character(seqnames(utr5)), start(utr5), end(utr5), as.character(strand(utr5)), sep="_")
	utr3$name <- paste(as.character(seqnames(utr3)), start(utr3), end(utr3), as.character(strand(utr3)), sep="_")
} else {
	## Set "name" field in metadata columns
	cds$name <- names(cds)
	utr5$name <- names(utr5)
	utr3$name <- names(utr3)
	## Delete "names" attribute (necessary for BED export)
	names(cds) <- NULL
	names(utr5) <- NULL
	names(utr3) <- NULL
}

#---> EXPORT CDS AND UTR REGIONS AS BED FILES <---#
# Print status message
if ( opt$verbose ) cat("Exporting ranges to BED files...\n", sep="")
## Export BED files
export(object=cds, con=paste0(opt$prefix, "cds.bed"), format="bed")
export(object=utr5, con=paste0(opt$prefix, "utr5.bed"), format="bed")
export(object=utr3, con=paste0(opt$prefix, "utr3.bed"), format="bed")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#