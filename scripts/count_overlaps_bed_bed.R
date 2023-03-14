#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		24-JAN-2013
### Modified:		24-JAN-2013
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, rtracklayer, GenomicRanges
### Description: 	Generate a count table from one BED file each of features and reads
### Arguments: 		1. BED file of reads; 2. BED file of features; 3. Minimum overlap (in nucleotides); 4. output file 
### Output: 		Tab-delimited count table with sample name from feature file and number of overlapping reads 
### Usage:			Rscript count_overlaps_bed_bed.R reads.bed features.bed 25 count_table.tab
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Load libraries
library(rtracklayer, GenomicRanges)
###

### B. Import data
# Use rtracklayer::import.bed() to import reads as GRanges object
reads <- import.bed(args[1], asRangedData = FALSE, genome = 'mm9')
# Use rtracklayer::import.bed() to import features of interest as GRanges object
feat <- import.bed(args[2], asRangedData = FALSE, genome = 'mm9')
###

### C. Calculate overlaps
# GenomicRanges::countOverlaps()
overlap <- countOverlaps(feat, reads, minoverlap=as.integer(args[3]), type='any')
# Write count table to output file
write.table(overlap, file = args[4], quote = FALSE, sep = '\t', row.names = mcols(feat)$name, col.names = FALSE)
###

### D. Write log
writeLines("Read file:")
args[1]
writeLines("Feature file:")
args[2]
writeLines("Minimum overlap:")
args[3]
writeLines("Output file:")
args[4]
writeLines("Summary count table: ")
summary(overlap)
writeLines("\nSession info")
print.default(sessionInfo())
###