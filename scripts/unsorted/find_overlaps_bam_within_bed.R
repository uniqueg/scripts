### Author: Alexander Kanitz
### Created: 04-APR-2013
### Modified: 04-APR-2013
### Description: Calculates overlaps of reads in BAM files with features specified in a BED file using GenomicRanges::findOverlaps 
### Arguments: 1. BED file [FILE]; 2. Path containing sorted BAM (extension MUST be '.bam'; corresponding '.bam.bai' files with the same basename MUST be present!) [PATH]; 3. Genome (STRING | e.g. "hg19"); 4. Output file prefix, may contain path [STRING]
### Output: Count table specifying the number of reads for each feature
### Usage example: ./Rscript ./find_overlaps_bam_within_bed.R ./feature_file.bed /path/to/BAM/files/ hg19 /path/to/OUT/files/count_table  

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
bed <- args[1]
bam_path <- args[2]
gen <- args[3]
out_prefix <- args[4]
## Load libraries
library(Rsamtools)
library(rtracklayer)
###

### B. IMPORT FEATURES AND READ LIST OF BAM FILES
# Import specified BED feature file
features_gr <- import.bed(bed, genome=gen, asRangedData=FALSE)
# Create file list of BAM read files in indicated path
read_files <- list.files(bam_path,".bam$",full=TRUE)
### 

### C. LOOP OVER EACH BAM FILE
dump <- lapply(read_files, function(current_bam) {
	
	# Extract basename of current BAM file (removes path and '.bam' file extension)
	bam_basename <- sub(".bam", "", basename(current_bam))
	# Generate name for output file
	out_name <- paste(out_prefix, "_", bam_basename, ".tab", sep="")
	
	# Use Rsamtools::readBamGappedAlignments to read BAM file as Rsamtools::GappedAlignments object
	bam <- readBamGappedAlignments(current_bam)
	# Coerce to GenomicRanges::GRanges object
	reads_gr <- as(bam, "GRanges")
	# Set appropriate genome 
	genome(reads_gr) <- gen
	
	# Use GenomicRanges::findOverlaps to find reads in 'reads_gr' that overlap completely (overlap type "within"!) with one or more features in 'features_gr' 
	# Reads have to be used as 'query' in order to allow the desired usage of overlap type "within"; this requires transposition of the resulting IRanges::Hits object later on (see below)
	hits <- findOverlaps(reads_gr, features_gr, type="within")

	# Create count table of hits per feature from IRanges::Hits object with IRanges::as.table and convert to matrix
	# For IRanges::as.table to summarize reads per feature, it is necessary to swap query and subject by transposing the Hits object first with IRanges::t
	count_table <- as.matrix(as.table(t(hits)))
	## Add column label and (unique) row names
	dimnames(count_table) <- list(make.unique(features_gr$name), bam_basename)
	# Save count table
	write.table(count_table, file=out_name, quote=FALSE, sep="\t")
	# Print status message to <STDOUT>
	cat("BAM file '", current_bam, "' processed. Output written to '", out_name, "'.\n", sep = "")

})
###
