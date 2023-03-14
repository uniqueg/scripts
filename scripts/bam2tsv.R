#######
### GENERAL:
### --------
### Author: 	Alexander Kanitz
### Created:	30-NOV-2012
### Modified:	05-DEC-2012	
#######

#######
### FUNCTION:
### ---------
### The script transforms one or more bam files to tables indicating the following values for each alignment (tab-separated): 1. identifier, 2. chromosome, 3. strand, 4. start position, 5. end position, 6. sequence. It accepts the following arguments (see USAGE): 1. a folder containing INDEXED bam files, 2. a pattern (mind R usage of patterns!) indicating which files to process.
#######

#######
### ARGUMENTS:
### ----------
### 1. A folder containing INDEXED bam files. Note that NO '/' is accepted at the end of the folder and that the corresponding '.bam.bai' files for each bam file MUST be present in the same folder too.
### 2. A pattern (mind R usage of patterns!) indicating which files to process. Note that the files that are to be processed MUST have the file extension '.bam' which is automatically added and therefore should NOT be included in the pattern. If all bam files in the folder are to be processed (preferred, as the pattern matching was not tested extensively!), indicate "".
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each processed bam file, an output file with the same base name and the file extension '.tsv' is created in the input file folder.
#######

#######
### USAGE:
### ------
### Rscript /path/to/bam2tsv.R /path/to/files <pattern>
### Examples:
### 1. Rscript ~/bam2tsv.R ~/bam_files "" (processes all files with extension 'bam' in the folder ~/bam_files; script in home directory)
### 2. Rscript ~/bam2tsv.R . ^example (processes file 'example.bam' in the current folder; script in home directory) 
#######

#######
### REQUIREMENTS:
### -------------
### Package "Rsamtools" and dependencies
#######

#######
### OVERVIEW:
### ---------
### A. PRE-REQUISITES
### B. IMPORT BAM FILES AS BAMFILELIST
### C. EXTRACT PARAMETERS FOR EACH BAM FILE
### D. WRITE DATA FRAMES AS TABLES
#######

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
# Load libraries
library(Rsamtools)
###

### B. IMPORT BAM FILES AS BAMFILELIST
# Use base::list.files method to create a char vector of bam files
cv_fls <- list.files(args[1], paste(args[2], ".bam$", sep=""), full=TRUE)
# Use Rsamtools::BamFileList method to generate file list from char vector
l_bam <- BamFileList(cv_fls)
###

### C. CREATE DATA FRAMES FOR EACH BAM FILE
# Create param object defining the parameters to extract
param <- ScanBamParam(what=c("qname", "rname", "strand", "pos", "qwidth", "seq")) # other possible values: "flag", "mapq", "cigar", "mrnm", "mpos", "isize", "qual", different tags (see Rsamtools manual for more info)
# Create list of lists (level1: list of bam files; level2: list of parameters)
ll <- lapply(l_bam, function(u) scanBam(u, param=param))
# Generate list of data frames (level1: list bam files; level2: data frame of parameters)
ldf <- lapply(ll, function(u) {
  u <- as.data.frame(do.call("DataFrame", u))
  # Calculate alignment "end" position and substitute for "qwidth
  u$qwidth <- u$pos + u$qwidth - 1
  u
})
###

### D. WRITE DATA FRAMES AS TABLES
mapply(function(u,v) write.table(u, file=paste(args[1],"/",sub('bam$', 'tsv', basename(v)), sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE), ldf, names(ldf))
###
