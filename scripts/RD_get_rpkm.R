#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Jul 16, 2013
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
		make_option(c("-i", "--intermediate"), action="store", type="character", default="", help="REQUIRED: Intermediate (step 1) RD output file.", metavar="file"),
		make_option(c("-f", "--fraction"), action="store", type="character", default="", help="REQUIRED: Table of transcripts and fractions (obtained by parsing final (step 2) RD output file (format: ENSEMBL transcript ID \tab\ estimated fraction of expressioin of corresponding gene.", metavar="file"),
		make_option(c("-l", "--length"), action="store", type="character", default="", help="REQUIRED: Table of transcript lengths (format: ENSEMBL transcript ID \tab\ size in nt.", metavar="file"),
		make_option(c("-o", "--output"), action="store", type="character", default="", help="REQUIRED: Output filename", metavar="file"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)

## Description
description <- "Calculates RPKM values from RD's intermediate and final output. Requires a table of transcript sizes. Only works with ENSEMBL gene and transcript IDs of mouse and human. No file checks performed.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0 (16-JULY-2014)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --intermediate <FILE> --fraction <FILE> --length <FILE> --output <FILE>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$intermediate == "" || opt$fraction == "" || opt$length == "" || opt$output == "" ) { 
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

#---> LOAD RD INTERMEDIATE DATA <---#
# Print status message
if ( opt$verbose ) cat("Reading RD intermediate (step 1) from file '", basename(opt$intermediate), "'...\n", sep="")
# Load RD step 1 data
step_1 <- read.table(opt$intermediate, colClasses=c("factor", "factor", "factor", "integer", "integer", "integer", "character"))

#---> PROCESS RD INTERMEDIATE DATA <---#
# Print status message
if ( opt$verbose ) cat("Processing RD intermediate (step 1) data...\n", sep="")
# Subset transcript IDs per gene
ensg_enst <- step_1[grep("ENS(MUS)?T", step_1$V7),]
# Build list of transcript -> gene lookup tables (one list element per gene)
enst_ensg_ls <- apply(ensg_enst, 1, function(row) {
			enst <- unlist(strsplit(row[7], ","))
			data.frame(enst=enst, ensg=rep(row[1], length(enst)))
		})
# Build transcript -> gene lookup table
enst_ensg <- do.call(rbind, enst_ensg_ls)
# Subset read counts per fragment
df_frag_cnt <- step_1[grep("ENS(MUS)?T", step_1$V7, invert=TRUE),]
# Build gene -> read count lookup table
ensg_cnt <- aggregate(df_frag_cnt$V6, by=list(df_frag_cnt$V1), sum)
# Rename columns
colnames(ensg_cnt) <- c("ensg", "gene_count")
# Build gene -> transcript -> gene count lookup table
enst_gene_cnt <- merge(enst_ensg, ensg_cnt)

#---> LOAD TRANSCRIPT FRACTIONS <---#
# Print status message
if ( opt$verbose ) cat("Reading fractions of gene expression per transcript from file '", basename(opt$fraction), "'...\n", sep="")
# Load RD transcript fractions
enst_frct <- read.table(opt$fraction, colClasses=c("factor", "numeric"), col.names=c("enst", "fraction"))

#---> CALCULATE CPM <---#
# Print status message
if ( opt$verbose ) cat("Calculating counts per million reads (CPM)...\n", sep="")
# Merge counts and fraction tables
enst_cpm <- merge(enst_frct, enst_gene_cnt)
# Remove lines with NAs
enst_cpm <- na.omit(enst_cpm)
# Multiply fractions and counts with each other
enst_cpm <- cbind(enst_cpm, trx_count=enst_cpm$fraction * enst_cpm$gene_count)
# Calculate total library size
lib_size <- sum(enst_cpm$trx_count)
# Divide by library size and multiply with 10^6
enst_cpm <- cbind(enst_cpm, cpm=enst_cpm$trx_count / lib_size * 1000000)

#---> LOAD TRANSCRIPT LENGTHS <---#
# Print status message
if ( opt$verbose ) cat("Reading transcript lengths from file '", basename(opt$length), "'...\n", sep="")
# Load transcript length lookup table
enst_len <- read.table(opt$length, colClasses=c("factor", "integer"), col.names=c("enst", "length"))

#---> CALCULATE RPKM <---#
# Print status message
if ( opt$verbose ) cat("Calculating reads per kilobase per million reads (RPKM)...\n", sep="")
# Merge counts/fractions and length tables
enst_rpkm <- merge(enst_cpm, enst_len)
# Divide by transcript length and multiply with 10^3
enst_rpkm <- cbind(enst_rpkm, rpkm=enst_cpm$cpm / enst_rpkm$length * 1000)

#---> EXPORT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writing output to file '", opt$output, "'...\n", sep="")
# Write results to file
write.table(enst_rpkm, opt$output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#