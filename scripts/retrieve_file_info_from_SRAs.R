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

## Status message
if ( opt$verbose ) cat("Starting 'retrieve_file_info_from_SRA.R'...\n\n")

## IF any required parameter is missing...
if 	(
		opt$`sql-dir`	== ""	||
		opt$query	== ""		||
		opt$`out-file` == ""
		) { 
	## ...write error message and die
	cat("[ERROR] Required argument(s) missing!\n", file=stderr())	
	print_help(opt_parser)
	quit(status=1)
}

# Split --query
queries <- unlist(strsplit(opt$query, "\\|"))

# Global variable for accessions that were not found
not_found <- NULL

# Global error code
error_code <- 0
###


### C. MAIN
## If present, remove trailing '/' from folders
opt$`sql-dir` <- sub("/$", "", opt$`sql-dir`)

# Generate absolute filename for SQLite database file
sql_file <- paste(opt$`sql-dir`, opt$`sql-db-file`, sep="/")

# Connect to SRA
sra_con <- dbConnect(SQLite(), sql_file)

## Traverse over all queries
runs <- lapply(queries, function(query) {
	# IF query ID is not available
	if ( ! dim(sraConvert(query, out_type=c("sra"), sra_con))[1] ) {
		## ...add query to dedicated vector and return NULL
		not_found <<- c(not_found, query)
		return(NULL)

	} else {
		# Else obtain all SRA type entities associated with the query ID
		df <- suppressWarnings(getSRAinfo(query, sra_con, sraType = opt$type))
		df$accession = rep(query, nrow(df))
		return(df)
	}
})

# Disconnet from SRA
invisible(dbDisconnect(sra_con))

## Combine dataframes for each query into one
all_runs <- do.call(rbind, runs)

## IF the resulting dataframe has at least one row (i.e. at least one SRA accession was found)...
if ( ! is.null(all_runs) ) {
	# Generate data.frame for output table 
	out_df <- data.frame(ftp=all_runs$ftp, accession=all_runs$accession, study=all_runs$study, sample=all_runs$sample, experiment=all_runs$experiment, run=all_runs$run, file_type=rep(opt$type, nrow(all_runs)))
	# Write data frame to file
	write.table(out_df, file=opt$`out-file`, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
## ELSE set error code for non-existing output file
} else {
	error_code <- 2
}
	
## Write unavailable queries to dedicated files, give warning and modify error code
if ( ! is.null(not_found) ) {
	not_found_filename <- paste0(opt$`out-file`, "_not_found")
	cat("[WARNING] One or more query accessions were not found in database! File '", not_found_filename, "'with unavailable accessions written.\n", file=stderr(), sep="")	
	write.table(not_found, file=not_found_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
	error_code <- error_code + 3;
}

## Status message
if ( opt$verbose ) cat("\nDone.\n")

## Quit with error code
quit(status=error_code)
###