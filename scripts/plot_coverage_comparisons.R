## Alexander Kanitz
## 18-DEC-2013

############ Pass ###############
# 1. Number of files in group 1 #
# 2. Number of files in group 2 #
# 3. Prefix for output files    #
# 4+. Input filenames group 1/2 #
#################################

## Get/parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)
no_grp1 <- as.integer(args[1])
no_grp2 <- as.integer(args[2])
prefix <- args[3]
files <- args[-1:-3]

## Validate command-line arguments
if ( ! is.integer(no_grp1) ) {
	stop("First argument must be an integer!\nExecution aborted.")
}
if ( ! is.integer(no_grp2) ) {
	stop("Second argument must be an integer!\nExecution aborted.")
}
if ( length(files) != no_grp1 + no_grp2 ) {
	stop("Insufficient/incompatible number of arguments!\nExecution aborted.")
}
if ( length(files) == 0 ) {
	stop("No files indicated!\nExecution aborted.")
}

## Define column names
col_names <- c("rseq", "source", "type", "start", "stop", "score", "str", "phase", "attr", "pos", "cov")

## Extract library numbers
file_base <- basename(files)
lib_no <- gsub(".bed", "", file_base, fixed=TRUE)
lib_no <- gsub("miRs_", "", lib_no, fixed=TRUE)
lib_no <- gsub("of_interest_", "", lib_no, fixed=TRUE)

# Load first data set
df <- read.delim(files[1], header=FALSE, stringsAsFactors=FALSE, col.names=col_names)

## Extract names and IDs and parents
df$name <- sub("Name=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 3)))
df$id <- sub("ID=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 1)))
df$parent <- sub("Derives_from=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 4)))

## Subset data by type
mat <- unique(df[df$type == "miRNA", c("start", "stop", "name", "id", "parent")])
pre <- unique(df[df$type == "miRNA_primary_transcript", c("rseq", "source", "type", "start", "stop", "score", "str", "phase", "attr", "name", "id")])

# Convert precursor dataframe to list and add coverage
pre_ls <- lapply(split(pre, pre$name), as.list)

## Add mature data to precursor list
pre_ls <- lapply(pre_ls, function(pre) {
	positions <- which(mat$parent == pre$id)
	pre$mature <- list()
	for (pos in positions) {
		pre$mature[[mat$name[pos]]] <- list(name=mat$name[pos], id=mat$id[pos], start=mat$start[pos], stop=mat$stop[pos])
	}
	return(pre)
})

# Initialize container for coverage list element names for different coverage files
cov_names <- NULL

## Load data files and add coverages to pre_ls
for (file in 1:length(files)) {

	# Load all data sets
	df <- read.delim(files[file], header=FALSE, stringsAsFactors=FALSE, col.names=col_names)

	# Extract names
	df$name <- sub("Name=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 3)))

	# Subset coverage for precursors
	cov <- df[df$type == "miRNA_primary_transcript", c("name", "cov")]

	# Convert precursor dataframe to list and add coverage
	cov_ls <- split(cov, cov$name)
	
	# Convert coverage counts to fractions
	cov_ls <- lapply(cov_ls, function(pre) {
		tot_cov <- sum(pre$cov)
		if ( tot_cov != 0 ) pre$cov <- pre$cov / tot_cov
		return(pre)
	})

	# Generate coverage element name
	cov_name <- paste0("cov", file)
	
	# Add coverage element name to dedicated list element
	cov_names <- c(cov_names, cov_name)

	# Combine lists
	for (name in names(cov_ls)) pre_ls[[name]][[cov_name]] <- cov_ls[[name]]$cov

}

## Plot
invisible(lapply(pre_ls, function(pre) { 
		
	# Get x coordinates
	x = seq(pre$start, pre$stop, 1)

	# Set labels
	main = paste0(pre$name, " / ", pre$id);
	xlab = pre$rseq
	ylab = "Coverage [nt]"

	# Open graphics device
	pdf(paste0(prefix, "coverage_", pre$id, ".pdf"))
	
	# Plot coverage
	plot(x, rep(0, length(x)), type="n", main=main, xlab=xlab, ylab=ylab, xaxt="n", ylim=c(-0.005,0.08))

	# Add custom x axis
	axis(1, at=c(pre$start, (pre$start + pre$stop) / 2, pre$stop))
	
	# Initialize line type and group size variables
	line_type = 1
	group_size = no_grp1
	group1_correction = 0
	
	# For each coverage profile...
	for (i in 1:length(cov_names)) {
		# Change line type and group size if group 2 is reached 
		if ( i > no_grp1 ) { 
			line_type = 6
			group_size = no_grp2
			group1_correction = no_grp1
		}
		
		# Plot line
		lines(x, pre[[cov_names[i]]], lty=line_type, col=(i-group1_correction+1))
		
	}

	## Add mature miRNAs
	invisible(lapply(pre$mature, function(mat) {
		x_mat <- seq(mat$start, mat$stop, 1)
		y_mat <- rep(0, length(x_mat))
		lines(x_mat, y_mat, col=1, lwd=2, lty=2)
		text((mat$start + mat$stop) / 2, 0, labels=mat$name, pos=1, offset=0.5, col=1)
	}))

	# Plot legend
	legend("topright", legend=lib_no, col=c(2:(no_grp1+1), 2:(no_grp2+1)), lty=c(rep(1,no_grp1), rep(6, no_grp2)))
	
	# Close graphics device
	dev.off()

}))
