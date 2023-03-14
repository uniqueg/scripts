#!/usr/bin/Rscript

### A. PRE-REQUISITES
## Load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("SRAdb"))
###

### B. PARSE ARGUMENTS
## List of allowed/recognized parameters
option_list <- list(
		make_option(c("-d", "--sql-dir"), action="store", type="character", default="", help="REQUIRED: Location of SQLite database file 'SRAmetadb.sqlite'", metavar="path"),
		make_option(c("-q", "--query"), action="store", type="character", default="", help="REQUIRED: Query; identifier of type SR[A|P|X|S|R] (NCBI Sequence Read Archive) or DR[A|P|X|S|R] (DNA Data Bank of Japan Read Archive); all available data sets associated with the query and of the indicated type will be downloaded", metavar="file"),
		make_option(c("-o", "--out-file"), action="store", type="character", default="", help="REQUIRED: Filename for output table", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
		make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!"),
		make_option(c("-i", "--sql-db-file"), action="store", type="character", default="SRAmetadb.sqlite", help="Name of SQLite database file [DEFAULT: 'SRAmetadb.sqlite']", metavar="file"),
		make_option(c("-t", "--type"), action="store", type="character", default="fastq", help="File type of data to download; options: 'fastq', 'sra', 'litesra' [DEFAULT: 'fastq']", metavar="string")
)

## Parse command-line arguments
# Initiate OptionParser
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --sql-dir [PATH] --query [STRING] --out-file [FILE]", option_list = option_list, add_help_option=FALSE)
# Parse options and write to list object 'opt'
opt <- parse_args(opt_parser)

## IF any required parameter is missing...
if 	(
		opt$`sql-dir`	== ""	||
		opt$query	== ""		||
		opt$`out-file` == ""
		) { 
	## ...write error message and die
	write("\n[ERROR] Required argument(s) missing!\n", stderr())	
	print_help(opt_parser)
	quit(status=2)
}
###


### C. MAIN
## Status message
if ( opt$verbose ) cat("\nStarting 'retrieve_file_info_from_SRA.R'...\n")

## If present, remove trailing '/' from folders
opt$`sql-dir` <- sub("/$", "", opt$`sql-dir`)

# Generate absolute filename for SQLite database file
sql_file <- paste(opt$`sql-dir`, opt$`sql-db-file`, sep="/")

# Connect to SRA
sra_con <- dbConnect(SQLite(), sql_file)

# Return with error if query ID is not available
if ( ! dim(sraConvert(opt$query, out_type=c("sra"), sra_con))[1] ) {
	## ...write error message and die
	write("\n[ERROR] Query accession not found in database!\n", stderr())	
	quit(status=3)
}

# Obtain all SRA type entities associated with the query ID
runs <- suppressWarnings(getSRAinfo(opt$query, sra_con, sraType = opt$type))

# Generate data.frame for output table
out_df <- data.frame(ftp=runs$ftp, study=runs$study, sample=runs$sample, experiment=runs$experiment, run=runs$run, file_type=rep(opt$type, nrow(runs)))

# Write data frame to file
write.table(out_df, file=opt$`out-file`, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# Disconnet from SRA
invisible(dbDisconnect(sra_con))

# Status message
if ( opt$verbose ) cat("\nDone.\n\n")

# Quit
if ( nrow(out_df) == 0 ) quit(status=1) else quit(status=0)
###
