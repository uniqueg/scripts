### Author: Alexander Kanitz
### Created: 22-MAR-2013
### Modified: 22-MAR-2013
### Description: Calculates overlaps of reads in BAM files with regions specified in a BED file using GenomicRanges::summarizeOverlaps 
### Arguments: 1. Path containing sorted BAM (BAI files MUST be present!) [PATH]; 2. BED file [FILE]; 3. Genome (STRING | e.g. "hg19"); 4. Overlap mode (STRING | e.g. "Union); 5. Output file [FILE]
### Output: SummarizedExperiment R object containing the count table 
### Usage: Rscript ./summarize_overlaps.R /path/to/BAM/files/ ./annotation_file.bed hg19 Union ./overlaps.R 

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Load libraries
library(Rsamtools)
library(rtracklayer)
###

### B. LOAD DATA
## List of read files
read_files <- list.files(args[1],".bam$",full=TRUE)
read_files_ls <- BamFileList(read_files)
# Import BED annotation
anno <- import.bed(args[2], genome=args[3], asRangedData=FALSE)
###

###. C. OVERLAPS
overlaps <- summarizeOverlaps(anno, read_files_ls, mode = args[4], ignore.strand=TRUE)
save(overlaps, file=args[5])
###
