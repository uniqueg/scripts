#######
### GENERAL:
### --------
### Author: 		Alexander Kanitz
### Created: 		20-JAN-2013
### Modified:		20-FEB-2013
### Adapted from:	subset_overlaps_bed_bed_overlap_no.R
### Language: 		R
### Version:		2.15.2
### Requirements:	Bioconductor_2.11, rtracklayer, GenomicRanges
### Description: 	Generate an overlap table from two BED files (query and subject)
### Arguments: 		1. query BED file with features of interest; 2. subject BED file with features of interest; 3. Genome (hg19, mm9, etc.); 4. Minimum overlap (in nucleotides); 5. score threshold below which intervals are discarded; 6. output file 
### Output: 		Tab-delimited count table with sample name from query file and average score as score
### Usage:			Rscript count_overlaps_bed_bed.R query.bed subject.bed hg19 25 0.5 overlaps_out.bed
### Notes:			Works/makes sense only if the number of overlaps between query and subject and subject and query are equal
#######

### A. Pre-requisites
# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Load libraries
library(rtracklayer)
###

### B. Import data
# Use rtracklayer::import.bed() to import features of interest as GRanges object
query <- import.bed(args[1], asRangedData = FALSE, genome = args[3])
# Use rtracklayer::import.bed() to import features of interest as GRanges object
subject <- import.bed(args[2], asRangedData = FALSE, genome = args[3])
###

### C. Subset data
# Subset query according to specified score threshold
query <- query[query$score >= as.numeric(args[5])]
# Subset subject according to specified score threshold
subject <- subject[subject$score >= as.numeric(args[5])]

### D. Calculate overlaps
# GenomicRanges::subsetByOverlaps()
overlap <- subsetByOverlaps(query, subject, minoverlap=as.integer(args[4]), type='any')
# Get corresponding lookup values for each overlap pair
lookup <- findOverlaps(query, subject, minoverlap=as.integer(args[4]), type='any')
# Calculate average score
mcols(overlap)$score <- (query[lookup@queryHits]$score + subject[lookup@subjectHits]$score) / 2

# Write overlap table to output file
export.bed(overlap, args[6])
###

### D. Write log
writeLines("Query file:")
args[1]
writeLines("Subject file:")
args[2]
writeLines("Genome:")
args[3]
writeLines("Minimum overlap:")
args[4]
writeLines("Score threshold:")
args[5]
writeLines("Output file:")
args[6]
writeLines("\nSession info")
print.default(sessionInfo())
###