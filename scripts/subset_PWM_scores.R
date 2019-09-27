#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 12, 2013
### Modified: Nov 12, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: optparse
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--input"), action="store", type="character", default="", help="REQUIRED: Input file (output of 'ag-scan-with-PWM.pl' by Andreas R. Gruber)", metavar="file"),
		make_option(c("-p", "--prefix"), action="store", type="character", default="", help="REQUIRED: Prefix for output files (e.g. name of protein whose PWM was analyzed)", metavar="string"),
		make_option(c("-n", "--number"), action="store", type="integer", default="0", help="Number of sites written to output (one and only one of --number or --cutoff is required)", metavar="integer"),
		make_option(c("-c", "--cutoff"), action="store", type="numeric", default="-1", help="Minimum score cutoff (between 0 and 1) for sites to be written to output (one and only one of --number or --cutoff is required)", metavar="float"),
		make_option(c("-s", "--species"), action="store", type="character", default="", help="Species filter in full name as it appears in the miRNA fasta file (e.g. 'Homo sapiens')", metavar="string"),		
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --input [FILE] --prefix [FILE]", option_list = option_list, add_help_option=FALSE, description="\nSubsets top sites from the output of the script 'ag-scan-with-PWM.pl' by Andreas R. Gruber.")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$input	== ""	||	opt$prefix	== "" || (opt$number == 0 && opt$cutoff == -1 ) ) { 
	write("[ERROR] Required argument(s) missing!\n", stderr())	
	stop(print_help(opt_parser))
}

## Die if both --number and --cutoff is indicated
if 	( opt$number != 0 && opt$cutoff != -1 ) { 
	write("[ERROR] Too many arguments! Choose only one of '--number' and '--cutoff'.\n", stderr())	
	stop(print_help(opt_parser))
}

## Set switch for analysis type
if ( opt$number != 0 ) {
	number <- TRUE
} else {
	number <- FALSE
}

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> READ DATA <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$input), "'...\n", sep="")
# Read data
df <- read.delim(opt$input, header=FALSE, col.names=c("fasta_id", "start", "stop", "score", "seq"))

#---> REARRANGE DATA <---#
# Print status message
if ( opt$verbose ) cat("Rearranging data...\n", sep="")
## Deconvolute FASTA identifier
fasta_id <- strsplit(as.character(df$fasta_id), " ")
short_name <- sapply(fasta_id, "[", 1)
long_name <- sapply(fasta_id, function(element) {paste(element[3:length(element)], collapse=" ")})
hairpin_id <- sapply(fasta_id, "[", 2)
# Recompile dataframe
df <- cbind(short_name, long_name, hairpin_id, df[,-1])
# Sort by score column
df <- df[rev(order(df$score)),]
# Initialize output string
out_string <- opt$prefix

#---> SUBSET BY SPECIES (IF DESIRED) <---#
# Check if '--species' is set
if ( opt$species != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Subsetting data by indicated species...\n", sep="")
	# Subset data
	df <- df[grep(opt$species, df$long_name),]
	# Replace spaces with underscores in species name
	species <- gsub(" ", "_", opt$species)
	# Modify output string
	out_string <- paste(out_string, "species", species, sep="_")
}

#---> SELECT TOP SITES BY NUMBER OR SCORE CUTOFF <---#
# Print status message
if ( opt$verbose ) cat("Selecting top sites...\n", sep="")
# Check switch
if ( number ) {
	# Subset data
	df <- head(df, opt$number)
	# Modify output string
	out_string <- paste(out_string, "top", opt$number, sep="_")
} else {
	# Subset data
	df <- df[df$score >= opt$cutoff,]
	# Modify output string
	out_string <- paste(out_string, "score_cutoff", opt$cutoff, sep="_")
}

#---> WRITE DATA <---#
# Print status message
if ( opt$verbose ) cat("Writing data...\n", sep="")
# Write dataframe to file
write.table(df, paste0(out_string, ".tab"), quote=FALSE, sep="\t", row.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
