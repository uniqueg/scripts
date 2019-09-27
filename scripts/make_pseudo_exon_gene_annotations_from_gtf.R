#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Jun 25, 2013
### Modified: Jun 25, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: rtracklayer, optparse
#==================#
### Description: Merges exons in a GTF file into pseudo-exons and writes them back into a GTF file 
### Output: A GTF file of pseudoexons with the following metadata columns: type: "pseudoexon"; pseudo_exon_number; gene_id
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
		make_option(c("-i", "--in_file"), action="store", type="character", default="", help="REQUIRED: GTF input filename", metavar="file"),
		make_option(c("-o", "--out_file"), action="store", type="character", default="", help="REQUIRED: GTF output filename", metavar="file"),
		make_option(c("-g", "--genome"), action="store", type="character", default=NA, help="Genome identifier (e.g. 'hg19', 'mm10'...); may be unstable if input file contains sequence names that are not available in the indicated genome, therefore leave as/set to NA if problems occur [DEFAULT: NA]", metavar="string"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [PATH] --out-dir [STRING]", option_list = option_list, add_help_option=FALSE, description="\nDescription: Merges exons in a GTF file into pseudo-exons and writes them back into a GTF file")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$in_file	== ""	||	opt$out_file	== "" ) { 
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
if ( opt$verbose ) cat("Reading input file '", basename(opt$in_file), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gtf <- import(con=opt$in_file, genome=opt$genome, format="gtf", asRangedData=FALSE) 

#---> COMPILE LIST OF GENES/EXONS <---#
# Print status message
if ( opt$verbose ) cat("Compiling exon list...\n")
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
exons <- gtf[values(gtf)[["type"]] == "exon"]
# Split exons GRanges into GRangesList by 'gene_id' 
exons_by_gene <- split(exons, exons$gene_id)

#---> GENERATE PSEUDO-EXONS <---#
# Print status message
if ( opt$verbose ) cat("Merging exons...\n")
# For each gene (i.e. GRangesList element), merge overlapping exons into pseudo exons
pseudo_exons_by_gene <- endoapply(exons_by_gene, function(gr) reduce(gr)) 
# Print status message
if ( opt$verbose ) cat("Adding metadata...\n")
## Traverse over each gene, ...
pseudo_exons_by_gene <- endoapply(pseudo_exons_by_gene, function(gr) {
			# ... add type "pseudoexon"...
			gr$type <- "pseudoexon"
			# ... and add pseudoexon serial number for each gene
			gr$pseudo_exon_number <- 1:length(gr)
			# Return object
			return(gr)
})
# Unlist GRangesList to GRanges
pseudo_exons <- unlist(pseudo_exons_by_gene)
# Put names in meta column "gene_id" (can cause problems in 'names' attribute)
pseudo_exons$gene_id <- names(pseudo_exons)
# Reset names attribute
names(pseudo_exons) <- NULL

#---> EXPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Writinig output...\n")
# Write output file
export(object=pseudo_exons, con=opt$out_file, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
