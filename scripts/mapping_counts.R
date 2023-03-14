#######
### GENERAL:
### --------
### Author: 	Alexander Kanitz
### Created:	05-DEC-2012
### Modified:	07-DEC-2012	
#######

#######
### FUNCTION:
### ---------
### The script analyses one or more '.log' result files of script "compare_tsv.pl" for the benchmarking of short read mappers. The script sums up and computes the fractions of "TRUE", "FALSE" and "UNMAPPED" alignments for each type of sequence, different sequence lengths and a combination of the two.
#######

#######
### ARGUMENTS:
### ----------
### 1. A folder containing the result files from "compare_tsv.pl". Note that NO '/' is accepted at the end of the folder.
### 2. A pattern (mind R usage of patterns!) indicating which files to process. Note that the files that are to be processed MUST have the file extension '.tsv.log' which is automatically added and therefore should NOT be included in the pattern. If all '.tsv.log' files in the folder are to be processed (preferred, as the pattern matching was not tested extensively!), indicate "".
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each processed file, an output file with the same base name and the suffix '_count.tab' is created in the input file folder.
#######

#######
### USAGE:
### ------
### Rscript /path/to/mapping_counts.R path/to/files <pattern>
### Examples:
### 1. Rscript ~/mapping_counts.R ~/log_files "" (processes all files with extension '.tsv.log' in the folder ~/log_files; script in home directory)
### 2. Rscript ~/mapping_counts.R . ^example (processes file 'example.tsv.log' in the current folder; script in home directory) 
#######

#######
### OVERVIEW:
### ---------
### A. PRE-REQUISITES
### B. ANALYSIS
#######

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
###

### B. ANALYSIS

# Use base::list.files method to create a list of files to be processed
ls_fls <- as.list(list.files(args[1], paste(args[2],".tsv.log$", sep=""), full=TRUE))
# Use lapply to traverse through file list
discard <- lapply(ls_fls, function(u) {
	
	## B0. LOAD DATA	
	df <- read.table(u, colClasses = c("NULL", "factor", "integer", "factor"), col.names = c("id", "type", "size", "ident"))
	##
	
	## B1. ANALYSE BY TYPE
	aggr_tp <- aggregate(x = df$type, by = list(df$type, df$ident), FUN = length)
	aggr_tp <- merge(merge(aggr_tp[aggr_tp$Group.2 == "TRUE", c(1,3)], aggr_tp[aggr_tp$Group.2 == "FALSE", c(1,3)], by = 1, all = TRUE), aggr_tp[aggr_tp$Group.2 == "UNMAPPED", c(1,3)], by = 1, all = TRUE)
	mt_tp <- as.matrix(aggr_tp[,2:4])
	mt_tp <- rbind(mt_tp, colSums(mt_tp, na.rm = TRUE))
	rownames(mt_tp) <- c(as.vector(aggr_tp[,1]), "TOTAL TYPE")
	mt_tp <- cbind(mt_tp, rowSums(mt_tp, na.rm = TRUE)) 
	mt_tp <- cbind(mt_tp, mt_tp[,1:4] / mt_tp[,4] * 100)
	colnames(mt_tp) <- c("True", "False", "Unmapped", "TOTAL", "True %", "False %", "Unmapped %", "TOTAL %")
	##
	
	## B2. ANALYSE BY SIZE
	aggr_sz <- aggregate(x = df$size, by = list(df$size, df$ident), FUN = length)
	aggr_sz <- merge(merge(aggr_sz[aggr_sz$Group.2 == "TRUE", c(1,3)], aggr_sz[aggr_sz$Group.2 == "FALSE", c(1,3)], by = 1, all = TRUE), aggr_sz[aggr_sz$Group.2 == "UNMAPPED", c(1,3)], by = 1, all = TRUE)
	mt_sz <- as.matrix(aggr_sz[,2:4])
	mt_sz <- rbind(mt_sz, colSums(mt_sz, na.rm = TRUE))
	rownames(mt_sz) <- c(as.vector(aggr_sz[,1]), "TOTAL SIZE")
	mt_sz <- cbind(mt_sz, rowSums(mt_sz, na.rm = TRUE)) 
	mt_sz <- cbind(mt_sz, mt_sz[,1:4] / mt_sz[,4] * 100)
	colnames(mt_sz) <- c("True", "False", "Unmapped", "TOTAL", "True %", "False %", "Unmapped %", "TOTAL %")
	##
	
	## B3. ANALYSE BY TYPE & SIZE
	aggr_tp_sz <- aggregate(x = df$size, by = list(df$type, df$size, df$ident), FUN = length)
	aggr_tp_sz$Group.1 = paste(aggr_tp_sz$Group.2, aggr_tp_sz$Group.1, sep=":")
	aggr_tp_sz <- merge(merge(aggr_tp_sz[aggr_tp_sz$Group.3 == "TRUE", c(1,4)], aggr_tp_sz[aggr_tp_sz$Group.3 == "FALSE", c(1,4)], by = 1, all = TRUE), aggr_tp_sz[aggr_tp_sz$Group.3 == "UNMAPPED", c(1,4)], by = 1, all = TRUE)
	mt_tp_sz <- as.matrix(aggr_tp_sz[,2:4])
	mt_tp_sz <- rbind(mt_tp_sz, colSums(mt_tp_sz, na.rm = TRUE))
	rownames(mt_tp_sz) <- c(as.vector(aggr_tp_sz[,1]), "TOTAL SIZE:TYPE")
	mt_tp_sz <- cbind(mt_tp_sz, rowSums(mt_tp_sz, na.rm = TRUE)) 
	mt_tp_sz <- cbind(mt_tp_sz, mt_tp_sz[,1:4] / mt_tp_sz[,4] * 100)
	colnames(mt_tp_sz) <- c("True", "False", "Unmapped", "TOTAL", "True %", "False %", "Unmapped %", "TOTAL %")
	##
	
	## B4. COMBINE
	mt <- rbind(mt_tp, mt_sz, mt_tp_sz)
	write.table(mt, file=paste(args[1],"/",sub('.tsv.log$', '_count.tab', basename(u)), sep=""), sep="\t", quote = FALSE)
	return(mt)
	##
})
