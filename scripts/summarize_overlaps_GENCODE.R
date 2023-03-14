### Author: Alexander Kanitz
### Created: 22-MAR-2013
### Modified: 24-MAY-2013 by Christina Herrmann
### Description: Calculates overlaps of reads in BAM files with regions specified in an annotation file (GRangesList) using GenomicRanges::summarizeOverlaps 
### Arguments: 1. Path containing sorted BAM (BAI files MUST be present!) [PATH]; 2. annotation file [FILE]; 3. Genome (STRING | e.g. "hg19"); 4. Overlap mode (STRING | e.g. "Union", "IntersectionStrict", "IntersectionNotEmpty"); 5. Ignore strands? [LOGICAL | "TRUE" or "FALSE"]; 6. Output file prefix [STRING | may include PATH]
### Output: Count table [TAB], relevant objects [R] 
### Usage: Rscript ./summarize_overlaps_GENCODE.R /path/to/BAM/files/ ./annotation_file.Rdata hg19 IntersectionStrict TRUE ./overlaps 

### BAM files "stemcells unique_mappers": /import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/input/bam/unique_mappers_bam/TEST
### output folder: /import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/output/pseudo_exons/
### 	/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/input/bam/unique_mappers_bam/TEST
### GENCODE "pseudo-exons" annotation file: /import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/info/hs71_gene_list_w_pseudo_exons.Rdata

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
bam_path <- as.character("/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/TEST_TEMP") #as.character(args[1])
feat_file <- as.character("/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/info/hs71_gene_list_w_pseudo_exons.Rdata")#as.character(args[2])
genome <- as.character("hg19") #as.character(args[3])
mode <- as.character("IntersectionStrict") #as.character(args[4])
ignore_strand <- as.logical("TRUE") #as.logical(args[5])
out_prefix <- as.character("/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/TEST_TEMP/") # as.character(args[6])
## Load libraries
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))

### B. LOAD DATA
## List of read files
read_files <- list.files(bam_path,".bam$",full=TRUE)
# Import annotation/feature file
features <- get(load(feat_file))
###

### C. FIND OVERLAPS
## Create single count table for each BAM file with GenomicRanges::summarizeOverlaps
dump <- lapply(read_files, function(current_bam) {
	
	# Generate name for output file
	out_name <- paste(out_prefix, "count_table_", sub(".bam", "", basename(current_bam)), ".tab", sep="")
	
	# Use Rsamtools::readBamGappedAlignments to read BAM file as Rsamtools::GappedAlignments object
	bam <- readBamGappedAlignments(current_bam)
	
	# Find Overlaps
	overlaps <- summarizeOverlaps(features, bam, mode=mode, ignore.strand=ignore_strand)
	colnames(overlaps) <- sub(".bam", "", basename(current_bam))
		
	# Save count table
	write.table(assays(overlaps)$counts, file=out_name, quote=FALSE, sep="\t")
	# Print status message to <STDOUT>
	cat("BAM file '", basename(current_bam), "' processed. Output written to '", out_name, "'.\n", sep = "")

})
###
