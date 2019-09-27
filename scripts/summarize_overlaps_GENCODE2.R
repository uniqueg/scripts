### Author: Alexander Kanitz
### Created: 22-MAR-2013
### Modified: 06-MAY-2013
### Description: Calculates overlaps of reads in BAM files with regions specified in a BED file using GenomicRanges::summarizeOverlaps 
### Arguments: 1. Path containing sorted BAM (BAI files MUST be present!) [PATH]; 2. BED file [FILE]; 3. Genome (STRING | e.g. "hg19"); 4. Overlap mode (STRING | e.g. "Union); 5. Ignore strands? [LOGICAL | "TRUE" or "FALSE"]; 6. Output file prefix [STRING | may include PATH]
### Output: Count table [TAB], relevant objects [R] 
### Usage: Rscript ./summarize_overlaps.R /path/to/BAM/files/ ./annotation_file.bed hg19 Union FALSE ./overlaps 

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
bam_path <- as.character(args[1])
feat_file <-as.character(args[2])
genome <- as.character(args[3])
mode <- as.character(args[4])
ignore_strand <- as.logical(args[5])
out_prefix <- as.character(args[6])
## Load libraries
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))

### B. LOAD DATA
## Generate BamFileList from read file
read_files <- list.files(bam_path,".bam$",full=TRUE)
bam_files <- BamFileList(read_files)
# Import BED annotation
features <- import.bed(feat_file, genome=genome, asRangedData=FALSE)
# Split features by gene_id --> GRangesList by gene_id
features <- split(features, features$name)

### C. FIND OVERLAPS
overlaps <- summarizeOverlaps(features, bam_files, mode=mode, ignore.strand=ignore_strand)
###

### D. SAVE DATA
# Save important objects
save(bam_files, features, overlaps, file=paste(out_prefix, "objects.R", sep=""))
# Write features to tab file
write.table(assays(overlaps)$counts, file=paste(out_prefix, "count_table.tab", sep=""), quote=FALSE, sep="\t")
###


