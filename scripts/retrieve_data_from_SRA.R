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
	make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [DEFAULT]"),
	make_option(c("-s", "--silent"), action="store_false", dest="verbose", help="Shut up!"),
	make_option(c("-i", "--sql-db-file"), action="store", type="character", default="SRAmetadb.sqlite", help="Name of SQLite database file [DEFAULT: 'SRAmetadb.sqlite']", metavar="file"),
	make_option(c("-t", "--type"), action="store", type="character", default="best", help="File type of data to download; options: 'fastq', 'sra', 'litesra' and 'best' [DEFAULT]; 'best' chooses the type that has the most files associated with the query; in case of tie the order is fastq > sra > litesra", metavar="string"),
	make_option(c("-o", "--out-dir"), action="store", type="character", default="./", help="Output folder for retrieved files [DEFAULT: '.']", metavar="path"),
	make_option(c("-m", "--meta"), action="store_true", default=FALSE, help="Download only metadata")
)

## Parse command-line arguments
# Initiate OptionParser
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --sql-dir [PATH] --query [STRING]", option_list = option_list, add_help_option=FALSE)
# Parse options and write to list object 'opt'
opt <- parse_args(opt_parser)

## IF any required parameter is missing...
if 	(
	opt$`sql-dir`	== ""	||
	opt$query	== ""
	) { 
	## ...write error message and die
	write("\n[ERROR] Required argument(s) missing!\n", stderr())	
	stop(print_help(opt_parser))
}
###


### C. MAIN
## Status message
if ( opt$verbose ) cat("\nStarting 'retrieve_data_from_SRA.R'...\n")

## If present, remove trailing '/' from folders
opt$`sql-dir` <- sub("/$", "", opt$`sql-dir`)
opt$`out-dir` <- sub("/$", "", opt$`out-dir`)

# Generate absolute filename for SQLite database file
sql_file <- paste(opt$`sql-dir`, opt$`sql-db-file`, sep="/")

# Connect to SRA
sra_con <- dbConnect(SQLite(), sql_file)

# Obtain all SRA type entities associated with the query ID
assoc_entities <- sraConvert(opt$query, sra_con = sra_con)

# Subset all SRA entities of type 'run'
assoc_runs <- assoc_entities$run

## Check number of available files for each type
file_no_fq = length(listSRAfile(assoc_runs, sra_con, fileType='fastq')$ftp)
file_no_sra = length(listSRAfile(assoc_runs, sra_con, fileType='sra')$ftp)
file_no_lsra = length(listSRAfile(assoc_runs, sra_con, fileType='litesra')$ftp)

## IF --type option 'best' is set, determine which file type has the most files available, then set opt$type to the corresponding file type; in case of ties, the order 'fastq' > 'sra' > 'litesra' applies
if(opt$type == "best") {
	if(file_no_fq  >= file_no_sra && file_no_fq >= file_no_lsra) {
		opt$type <- "fastq"
	} else {
		if (file_no_sra >= file_no_lsra) {
			opt$type <- "sra"
		} else {
			opt$type <- "litesra"
		}
	}
}

# Determine files to download
runs_to_get <- listSRAfile(assoc_runs, sra_con, fileType=opt$type)$run

## IF number of available files for the specified/determined file type is smaller than the number of files associated with the query
if(length(runs_to_get) < length(assoc_runs)) {
	## ...issue a warning and list entries associated with query and available entries
	cat("\n[WARNING] Not all run data associated with the query are available for download.\n")
	cat("Runs associated with query: \n", assoc_runs, "\n")
	cat("Runs to be retrieved: \n", runs_to_get, "\n")
}

# Filter only available entries from associated entities object
assoc_entities <- subset(assoc_entities, subset=assoc_entities$run %in% runs_to_get)

# Create base output folder if not present
dir.create(opt$`out-dir`, showWarnings=FALSE, recursive=TRUE)

## Get metadata
# Extract unique submission accessions associated with query
subm_acc <- unique(assoc_entities$submission)
## For each submission accession, generate one list of meta data
meta_ls <- lapply(subm_acc, function(acc) {
	# Initialize empty list object
	ls <- list()
	## Obtain metadata from different tables (sra, submission, study, experiment, sample)		
	ls$sra <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM sra WHERE submission_accession = ", shQuote(acc))))
	ls$submission <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM submission WHERE submission_accession = ", shQuote(acc))))
	ls$study <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM study WHERE submission_accession = ", shQuote(acc))))
	ls$experiment <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM experiment WHERE submission_accession = ", shQuote(acc))))
	ls$sample <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM sample WHERE submission_accession = ", shQuote(acc))))
	## For each data frame in list object 'ls'...
	ls <- lapply(ls, function(df) {
		# ...remove all rows with only NA values,
		df <- df[ ,which(unlist(lapply(df, function(x)!all(is.na(x)))))]
		# ...transpose data frame and
		df <- t(df)
		# ...sort data frame
		df <- df[order(row.names(df)),]
	})
	# Return list
	return(ls)
})

## Write metadata to files
## Iterate over metadata list (one list element for each submission accession)
for (i in 1:length(subm_acc)) {
	## Apply over the names of the data frames present in the resulting list
	invisible(lapply(names(meta_ls[[i]]), function(df_name) {
		# ...and write data frame to file
		write.table(meta_ls[[i]][[df_name]], file=paste0(opt$`out-dir`, "/", subm_acc[i], "_meta_", df_name, ".tab"), col.names=FALSE, quote=FALSE, sep="\t", eol="\n")
	}))
}

## Apply over files to retrieve
invisible(lapply(runs_to_get, function(run) {
	## Folder name and generation
	# Generate output folder name
	out_folder <- paste(opt$`out-dir`, run, sep="/")
	# Create output directory
	dir.create(out_folder, showWarnings=FALSE, recursive=TRUE)
	
	## Run metadata
	# Get metadata for run
	sra <- as.data.frame(dbGetQuery(sra_con, paste0("SELECT * FROM run WHERE run_accession = ", shQuote(run))))
	# Remove all rows with only NA values from data frame
	sra <- sra[ ,which(unlist(lapply(sra, function(x)!all(is.na(x)))))]
	# Transpose data frame
	sra <- t(sra)
	# Sort data frame
	sra <- sra[order(row.names(sra)),]
	# Write data frame to file
	write.table(sra, file=paste0(out_folder, "/", run, "_meta_run.tab"), col.names=FALSE, quote=FALSE, sep="\t", eol="\n")
	
	## Download experiment data
	# Download file of type $type
	if (!opt$meta) {getSRAfile(run, sra_con, destDir=out_folder, fileType=opt$type)}
}))

# Disconnet from SRA
invisible(dbDisconnect(sra_con))

## Status message
if ( opt$verbose ) cat("\nDone.\n\n")
###