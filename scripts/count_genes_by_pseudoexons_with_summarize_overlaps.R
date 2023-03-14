#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 21, 2013
### Modified: Nov 21, 2013
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
#---> LOAD OPTION PARSER <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> SET SCRIPT DESCRIPTION & USAGE <---#
script_description <- "\nIntersects sorted, indexed BAM files with ENSEMBL annotations and generates gene counts based on pseudoexons (i.e. overlapping exons within a gene are merged). Different overlap modes can be chosen, but in all cases reads may increase the counts of at most one gene each. Strand sensitivity may be turned off. Due to sequence/chromosome naming conventions, it is recommended to map against genome versions provided by ENSEMBL, but genomes from UCSC (and potentially others) may be compatible for the main chromosomes if the corresponding option is set."
usage <- "Usage: %prog [OPTIONS] --gtf [GTF] --bam [PATH] --prefix [TAB]"

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-g", "--gtf"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="gtf"),
		make_option(c("-b", "--bam"), action="store", type="character", default="", help="REQUIRED: Path to BAM input files", metavar="path"),
		make_option(c("-p", "--prefix"), action="store", type="character", default="", help="REQUIRED: Output filename", metavar="tab"),
		make_option(c("-m", "--mode"), action="store", type="character", default="IntersectionStrict", help="Overlap mode; options: 'IntersectionStrict' (default), 'Union', 'IntersectionNotEmpty' (see GenomicRanges::summarizeOverlaps for details", metavar="STRING"),
		make_option(c("-i", "--ignore-strand"), action="store_true", default=FALSE, help="Ignore strand for overlap counts"),
		make_option(c("-r", "--rename-ucsc"), action="store_true", default=FALSE, help="Prepend chromosome names with 'chr' if mappings were done against a UCSC version of the genome"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage=usage, option_list = option_list, add_help_option=FALSE, description=script_description)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$gtf	== ""	||	opt$bam	== "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
}

# Remove trailing slash "/" from BAM path name if present
if ( substr(opt$bam, nchar(opt$bam), nchar(opt$bam)) == "/" ) { opt$bam <- substr(opt$bam, 1, nchar(opt$bam)-1) }

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
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("edgeR"))

#---> GENERATE FILE LIST OF BAM FILES <---#
# Print status message
if ( opt$verbose ) cat("Searching BAM files in folder '", basename(opt$gtf), "'...\n", sep="")
# Generate a character vector of absolute filenames of BAM files in the specified folder
files <- list.files(opt$bam, ".bam$", full=TRUE)
# Create a BamFileList object of BAM files
bfl <- BamFileList(files)
# Print status message
if ( opt$verbose ) cat("The following files were found:", files, sep="\n")

#---> IMPORT GTF ANNOTATION FILE <---#
# Print status message
if ( opt$verbose ) cat("Reading annotation file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> SYNCHRONIZE CHROMOSOME NAMES <---#
if ( opt$`rename-ucsc` ) {
	# Print status message
	if ( opt$verbose ) cat("Synchronizing chromosome names...\n")
	# Create Seqinfo object with new sequence names prepended by "chr" (this is not safe for chromosome fragments etc., which will be lost downstream)
	seqinfo <- Seqinfo(paste0("chr", seqlevels(gr)))
	# Replace old with new Seqinfo object
	seqinfo(gr, new2old=as.integer(1:length(seqinfo)), force=FALSE) <- seqinfo
}
	
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

#---> COUNT OVERLAPS <---#
# Print status message
if ( opt$verbose ) cat("Counting overlaps...\n")
# Use GenomicRanges::summarizeOverlaps method to count overlaps between reads and features
overlaps <- summarizeOverlaps(grl, bfl, mode=opt$mode, ignore.strand=opt$`ignore-strand`)

#---> EXTRACT & PROCESS COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Extracting and processing count table...\n")
# Extract count table
counts <- assays(overlaps)$counts
# Rename row names
colnames(counts) <- basename(colnames(counts))

#---> EXTRACT & PROCESS COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Extracting and processing count table...\n")
# Calculate counts per million reads (CPM)
counts_cpm <- cpm(counts)
# Calculate size of each gene (i.e. total number of nucleotides in all pseudogenes)
widths <- sum(width(grl))
# Re-arrange size vector in the same order as the counts_cpm table
widths <- widths[row.names(counts_cpm)]
# Calculate reads per kilobase per million reads (RPKM)
counts_rpkm <- cpm(counts) / widths * 1000

#---> SAVE COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writing output...\n")
# Write output file
save(grl, bfl, overlaps, counts, counts_cpm, counts_rpkm, file=paste(opt$prefix, "annotation_overlaps.R", sep=""))
## Write count tables to tab files
write.table(counts, file=paste(opt$prefix, "raw_counts.tab", sep=""), quote=FALSE, sep="\t")
write.table(counts_cpm, file=paste(opt$prefix, "CPM.tab", sep=""), quote=FALSE, sep="\t")
write.table(counts_rpkm, file=paste(opt$prefix, "RPKM.tab", sep=""), quote=FALSE, sep="\t")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#