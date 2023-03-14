#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Aug 7, 2014
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
description <- "Counts the overlaps between a BAM read file and a BED12 feature file.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 07-AUG-2014"
version <- "Version: 1.0 (07-AUG-2014)"
requirements <- "Requires: optparse, rtracklayer, GenomicAlignments"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> DEFINE COMMAND-LINE OPTIONS <---#
option_list <- list(
		make_option("--bed12", action="store", type="character", default="", help="REQUIRED: BED12 feature file.", metavar="file"),
		make_option("--bam", action="store", type="character", default="", help="REQUIRED: BAM read alignments file. A BAI index file must be present, either (1) with the same path and filename plus a '.bai' extension or (2) supplied via '--bai'.", metavar="file"),
		make_option("--bai", action="store", type="character", default="", help="BAM index file (default: argument to '--bam' + '.bai'", metavar="file"),		
		make_option("--counts", action="store", type="character", default="", help="REQUIRED: Output file; table of raw and normalized counts for each feature.", metavar="file"),
		make_option("--mode", action="store", type="character", default="Union", help="Counting mode (one of 'Union', 'IntersectionStrict' or 'IntersectionNotEmpty'; default: 'Union'); refer to the Bioconductor GenomicAlignments manual for details.", metavar="string"),
		make_option("--ignore-strand", action="store_true", default=FALSE, help="Do not consider the strand when counting read/feature overlaps."),
		make_option("--inter-feature", action="store_true", default=FALSE, help="Do not count reads intersecting multiple features."),
		make_option("--paired-end", action="store_false", default=TRUE, help="Specify when counting paired-end reads."),
		make_option("--include-singletons", action="store_true", default=FALSE, help="Specify if singletons should be included when counting paired-end reads."),		
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages (default)"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

#---> PARSE COMMAND-LINE OPTIONS <---#
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --bed12 [FILE] --bam [FILE] --counts [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> VALIDATE COMMAND-LINE OPTIONS <---#
## Die if any required arguments are missing...
if 	( opt$bed12 == "" || opt$bam == "" || opt$counts == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
## Set default BAM index filename if not supplied by user
if ( opt$bai == "" ) {
	opt$bai <- paste(opt$bam, "bai", sep=".")
}
## Die if illegal counting mode is specified
if ( ! ( opt$mode == "Union" || opt$mode == "IntersectionStrict" || opt$mode == "IntersectionNotEmpty" ) ) {
	write("[ERROR] Illegal argument to '--mode' option! Please use one of 'Union', 'IntersectionStrict' or 'IntersectionNotEmpty'.\n\n", stderr())
	stop(print_help(opt_parser))	
}
## Set '--include-singletons' to FALSE if reads are not paired-ended
if ( opt$`include-singletons` && ! opt$`paired-end`) {
	write("[WARNING] Option '--include-singletons' does not apply when reads are single-ended. Option ignored!\n", stderr())
	opt$`include-singletons` <- FALSE
}
#===================#
#    OPTIONS END    #
#===================#

#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> LOAD LIBRARIES <---#
# Print status message
if ( opt$verbose ) cat("Loading required libraries...\n")
# Load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }
#if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicAlignments"))) == FALSE ) { stop("Package 'GenomicAlignments' required!\nExecution aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("Rsamtools"))) == FALSE ) { stop("Package 'Rsamtools' required!\nExecution aborted.") }

#---> IMPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Importing features from BED12 file '", opt$bed12, "'...\n", sep="")
# Use rtracklayer::import method to import file as GRanges object 
gr <- import(con=opt$bed12, format="bed", asRangedData=FALSE)
# Use gene names for range names
names(gr) <- gr$name

#---> ACCESS BAM <---#
# Print status message
if ( opt$verbose ) cat("Accessing BAM file '", opt$bam, "'...\n", sep="")
# Use Rsamtools::BamFile and Rsamtools::BamFileList methods to access the file 
bfl <- BamFileList(BamFile(opt$bam, index=opt$bai))

#---> COUNTING OVERLAPS <---#
# Print status message
if ( opt$verbose ) cat("Counting overlaps...\n", sep="")
# Use Rsamtools::BamFile and Rsamtools::BamFileList methods to access the file 
se <- summarizeOverlaps(features=gr, reads=bfl, mode=opt$mode, ignore.strand=opt$`ignore-strand`, inter.feature=opt$`inter-feature`, singleEnd=opt$`paired-end`, fragments=opt$`include-singletons`)

#---> NORMALIZE READ COUNTS <---#
# Print status message
if ( opt$verbose ) cat("Normalizing raw counts...\n", sep="")
# Extract/calculate data for table columns
names=rowData(se)$name
raw_counts=as.integer(assays(se)$counts)
lib_size=sum(as.integer(assays(se)$counts))
widths=sum(width(rowData(se)$blocks))
# Build output data table
df <- data.frame(id=names, raw_counts=raw_counts, cpm=raw_counts/lib_size*1000000, size=widths, rpkm=raw_counts/lib_size/widths*1000000000)

#---> WRITE COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writing output to table '", opt$counts, "'...\n", sep="")
# Write to output file
write.table(df, file=opt$counts, quote=FALSE, sep="\t", row.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n\n")

#---> PRINT SESSION INFO <---#
if ( opt$verbose ) {
	cat("++++++++++++\nSESSION INFO\n++++++++++++\n\n")	
	writeLines(capture.output(sessionInfo()))
}
#================#
#    MAIN END    #
#================#
