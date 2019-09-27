#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Aug 6, 2014
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
description <- "From a GTF gene set annotation file, generates a BED12+3 file containing the exons of the longest transcript isoform of each gene (one gene/transcript per row).\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 06-AUG-2014"
version <- "Version: 1.0 (06-AUG-2014)"
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
feat_gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> SUBSET EXONS <---#
# Print status message
if ( opt$verbose ) cat("Subsetting exons...\n")
# Subset EXONS (discards all other categories, e.g. CDS, start_codon etc.)
ex_gr <- feat_gr[values(feat_gr)[["type"]] == "exon"]
# Test if at least one region is returned
if ( length(ex_gr) == 0 ) stop("No entries of type 'exon' in input file! Check file.\nExecution halted.\n")

#---> DETERMINE LONGEST TRANSCRIPTS FOR EACH GENE <---#
# Print status message
if ( opt$verbose ) cat("Finding longest transcripts...\n")
# Construct data frame of exons, denoting their gene IDs, transcript IDs and widths
ex_wd <- data.frame(gene=ex_gr$gene_id, transcript=ex_gr$transcript_id, width=width(ex_gr))
# Aggregate exon widths to transcript widths
trx_wd <- aggregate(ex_wd$width, by=list(ex_wd$gene, ex_wd$transcript), sum)
# Select longest transcripts per gene
gene_max_wd <- aggregate(trx_wd$x, by=list(trx_wd$Group.1) , max)
# Add transcript IDs
gene_max_wd_trx <- merge(gene_max_wd, trx_wd)
# Remove duplicates and extract transcript IDs
trx_keep <- gene_max_wd_trx[!duplicated(gene_max_wd_trx[,1:2]),][,3]

#---> SUBSET LONGEST TRANSCRIPTS <---#
# Print status message
if ( opt$verbose ) cat("Subsetting exons of longest transcripts...\n")
# Subset longest transcripts
ex_filt_gr <- subset(ex_gr, ex_gr$transcript_id %in% trx_keep)

#---> GROUP TRANSCRIPTS <---#
# Print status message
if ( opt$verbose ) cat("Grouping exons by transcript...\n")
# Combine gene and transcript IDs
ex_filt_gr$gene_trx_id <- paste(ex_filt_gr$gene_id, ex_filt_gr$transcript_id, sep="$")
# Group exons by combined gene/transcript IDs
trx_filt_grl <- split(ex_filt_gr, ex_filt_gr$gene_trx_id)

#---> EXPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Writing output to BED12 file '", opt$bed12, "'...\n", sep="")
# Write output file
export(object=trx_filt_grl, con=opt$bed12, format="bed15")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
