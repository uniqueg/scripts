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
	make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
	make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!"),
	make_option(c("-i", "--sql-db-file"), action="store", type="character", default="SRAmetadb.sqlite", help="Name for SQLite database file [DEFAULT: 'SRAmetadb.sqlite']", metavar="file"),
	make_option(c("-f", "--force"), action="store_true", default=FALSE, help="Force updating of SRA SQLite database file to indicated location, regardless of existence or time of last update")
)

## Parse command-line arguments
# Initiate OptionParser
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --sql-dir [PATH] --query [STRING]", option_list = option_list, add_help_option=FALSE)
# Parse options and write to list object 'opt'
opt <- parse_args(opt_parser)

## IF any required parameter is missing...
if(opt$`sql-dir` == "") { 
	## ...write error message and die
	write("\n[ERROR] Required argument(s) missing!\n", stderr())	
	stop(print_help(opt_parser))
}
###


### C. MAIN
## Status message
if ( opt$verbose ) cat("\nStarting program 'update_SRA_SQLite_db_file'...\n")

# If present, remove trailing '/' from folder
opt$`sql-dir` <- sub("/$", "", opt$`sql-dir`)

# Generate absolute filename for SQLite database file
sql_file <- paste(opt$`sql-dir`, opt$`sql-db-file`, sep="/")

## Check IF folder specified in '--sql-dir' exists
if(file.exists(opt$`sql-dir`)) {
	# Print status message
	if ( opt$verbose ) cat("\nSpecified directory '", opt$`sql-dir`, "' exists.\n", sep="")
## ELSE create directory
} else {
	# Print status message
	if ( opt$verbose ) cat("\nCreating directory: ", opt$`sql-dir`, "\n", sep="")
	# Make directory
	dir.create(opt$`sql-dir`, showWarnings = TRUE, recursive=TRUE)
}

## Check IF option '--force' is set
if(opt$force) {
	# Print status message
	if ( opt$verbose ) cat("\nForce updating/downloading SRA SQLite database file '", opt$`sql-db-file`, "'...\nDepending on the connection bandwith, this may take long!\n\n", sep="")
	# Download/update SRA SQLite database file
	sql_file <- getSRAdbFile(destdir=opt$`sql-dir`)
## ELSE test IF file does not exist
} else if (! file.exists(sql_file)) {
	# Print status message
	if ( opt$verbose ) cat("\nSRA SQLite database file '", opt$`sql-db-file`, "' does not exist at specified location.\nDownloading...\nDepending on the connection bandwith, this may take long!\n\n", sep="")
	# Download/update SRA SQLite database file
	sql_file <- getSRAdbFile(destdir=opt$`sql-dir`)	
## ELSE test IF file is outdated
} else if (!(Sys.Date() - as.Date(file.info(sql_file)$mtime) < 30)) {
	# Print status message
	if ( opt$verbose ) cat("\nSRA SQLite database file '", opt$`sql-db-file`, "' is outdated (older than 30 days).\nUpdating...\nDepending on the connection bandwith, this may take long!\n\n", sep="")
	# Download/update SRA SQLite database file
	sql_file <- getSRAdbFile(destdir=opt$`sql-dir`)	
## ELSE do nothing
} else {
	# Print status message
	if ( opt$verbose ) cat("\nSRA SQLite database file '", opt$`sql-db-file`, "' exists, is up to date (less than 30 days old) and the '--force' option has not been set. No download.\n", sep="")	
}

## Status message
if ( opt$verbose ) cat("\nDone.\n\n")
###
