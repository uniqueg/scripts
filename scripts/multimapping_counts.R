#######
### GENERAL:
### --------
### Author: 	Alexander Kanitz
### Created:	07-DEC-2012
### Modified:	12-DEC-2012	
#######

#######
### FUNCTION:
### ---------
### The script analyses one or more '.multi' result files of script "compare_tsv.pl" for the benchmarking of short read mappers. The script sums up and computes the ratios of found alignments and 'real' matches. Results are grouped by number of 'real' matches, different sequence lengths and a combination of the two.
#######

#######
### ARGUMENTS:
### ----------
### 1. A folder containing the result files from "compare_tsv.pl". Note that NO '/' is accepted at the end of the folder.
### 2. A pattern (mind R usage of patterns!) indicating which files to process. Note that the files that are to be processed MUST have the file extension '.tsv.multi' which is automatically added and therefore should NOT be included in the pattern. If all '.tsv.multi' files in the folder are to be processed (preferred, as the pattern matching was not tested extensively!), indicate "".
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each processed file, an output file with the same base name and the suffix '_count_mm.tab' is created in the input file folder.
#######

#######
### USAGE:
### ------
### Rscript /path/to/mapping_counts.R path/to/files <pattern>
### Examples:
### 1. Rscript ~/mapping_counts.R ~/multi_files "" (processes all files with extension '.tsv.multi' in the folder ~/multi_files; script in home directory)
### 2. Rscript ~/mapping_counts.R . ^example (processes file 'example.tsv.multi' in the current folder; script in home directory) 
#######

#######
### OVERVIEW:
### ---------
### A. PRE-REQUISITES
### B. LOAD DATA
### C. PREPARE GENERAL RESULTS MATRIX
### D. AGGREGATE RESULTS BY NUMBER OF MATCHES AND/OR SIZE
### E. COMBINE MATRICES
#######

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
###

### B. LOAD DATA
# Use base::list.files method to create a list of files to be processed
ls_fls <- as.list(list.files(args[1], paste(args[2],".tsv.multi$", sep=""), full=TRUE))
#ls_fls <- as.list(list.files(., paste("",".tsv.multi$", sep=""), full=TRUE))
print(ls_fls)
# Use lapply to traverse through file list
discard <- lapply(ls_fls, function(u) {
	# Read data table	
	df <- read.table(u, colClasses = c("NULL", "integer", "integer", "integer"), col.names = c("id", "matches", "size", "count"))
	#df <- read.table("test_mp_map.tsv.multi", colClasses = c("NULL", "integer", "integer", "integer"), col.names = c("id", "matches", "size", "count"))
	###

	### C. PREPARE GENERAL RESULTS MATRIX
	# Group by sequence size and 'real' number of matches and count occurences (i.e. number of sequences) per group
	mat <- as.matrix(aggregate(x = df$matches, by = list(df$size, df$matches), FUN = length)) 
	# Sum total found alignments per group
	cnt <- as.matrix(aggregate(x = df$count, by = list(df$size, df$matches), FUN = sum))
	# Calculate 'real' number of alignments (occurences * match)
	mt <- cbind(mat, mat[,2] * mat[,3], cnt[,3])
	# Assign column names
	colnames(mt) <- c("Size", "Possible matches per sequence", "Number of sequences", "Total number of possible matches", "Total number of observed matches")
	###
	
	### D. AGGREGATE
	
	## D1. BY SIZE AND NUMBER OF MATCHES
	# Subset general results matrix
	mt_mat_sz <- mt[,3:5]
	# Assign row names
	rownames(mt_mat_sz) <- paste(mt[, 1],":MULTIMAPPER",mt[,2],"TIMES", sep="")
	# Add "TOTALS" row with column sums
	mt_mat_sz <- rbind(mt_mat_sz, colSums(mt_mat_sz, na.rm = TRUE))
	rownames(mt_mat_sz)[length(rownames(mt_mat_sz))] <- "TOTAL BY NUMBER OF MATCHES AND SIZE"
	# Calculate fractions (observed/real alignments) in percent
	mt_mat_sz <- cbind(mt_mat_sz, mt_mat_sz[,3] / mt_mat_sz[,2] * 100)
	colnames(mt_mat_sz)[length(colnames(mt_mat_sz))] <- "Fraction %"	
	##
	
	## D2. BY NUMBER OF MATCHES
	# Aggregate general results matrix	
	tmp_mat <- aggregate(mt[,c(3:5)],by = list(mt[,2]), FUN = sum)
	# Subset
	mt_mat <- tmp_mat[,2:4]
	# Assign row names
	rownames(mt_mat) <- paste(tmp_mat[,1],"TIMES", sep="_")	
	# Add "TOTALS" row with column sums
	mt_mat <- rbind(mt_mat, colSums(mt_mat, na.rm = TRUE))
	rownames(mt_mat)[length(rownames(mt_mat))] <- "TOTAL BY NUMBER OF MATCHES"
	# Calculate fractions (observed/real alignments) in percent
	mt_mat <- cbind(mt_mat, mt_mat[,3] / mt_mat[,2] * 100)
	colnames(mt_mat)[length(colnames(mt_mat))] <- "Fraction %"		
	##
		
	## D3. BY SIZE
	# Subset and aggregate general results matrix	
	tmp_sz <- aggregate(mt[,3:5], by = list(mt[,1]), FUN = sum)
	# Subset
	mt_sz <- tmp_sz[,2:4]
	# Assign row names
	rownames(mt_sz) <- paste(tmp_sz[, 1],"NT", sep="_")
	# Add "TOTALS" row with column sums
	mt_sz <- rbind(mt_sz, colSums(mt_sz, na.rm = TRUE))	
	rownames(mt_sz)[length(rownames(mt_sz))] <- "TOTAL BY SIZE"
	# Calculate fractions (observed/real alignments) in percent
	mt_sz <- cbind(mt_sz, mt_sz[,3] / mt_sz[,2] * 100)
	colnames(mt_sz)[length(colnames(mt_sz))] <- "Fraction %"	
	###

	## E. COMBINE MATRICES
	mt <- rbind(mt_mat, mt_sz, mt_mat_sz)
	write.table(mt, file=paste(args[1],"/",sub('.tsv.multi$', '_count_mm.tab', basename(u)), sep=""), sep="\t", quote = FALSE)
	return(mt)
	##
})
