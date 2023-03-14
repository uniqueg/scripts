### Author: Alexander Kanitz
### Created: 06-MAY-2013
### Modified: 06-MAY-2013
### Description: Calculates overlaps of reads in BAM files with subranges of regions specified in a BED file; uses GenomicRanges::summarizeOverlaps; subrange size can be specified 
### Arguments: 1. Path containing sorted BAM (BAI files MUST be present!) [PATH]; 2. BED file [FILE]; 3. Genome (STRING | e.g. "hg19"); 4. Bin size (INTEGER | default: 100); 5. Overlap mode (STRING | e.g. "Union); 6. Ignore strands? [LOGICAL | "TRUE" or "FALSE"]; 7. Output file prefix [STRING | may include PATH]
### Output: Count table [TAB], subranges [BED], relevant objects [R] 
### Usage: Rscript ./summarize_overlaps_bins.R /path/to/BAM/files/ ./annotation_file.bed hg19 200 Union FALSE ./overlaps.R 

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
bam_path <- as.character(args[1])
feat_file <- as.character(args[2])
genome <- as.character(args[3])
bin_size <- as.numeric(args[4])
mode <- as.character(args[5])
ignore_strand <- as.logical(args[6])
out_prefix <- as.character(args[7])
## Load libraries
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
###

### B. FUNCTIONS 
## >> subGRanges <<
## Description: Generate subranges of specified size of GRanges object
## Accepts: 1. GRanges object; 2. Bin size (integer; default = 100); 3. Shall last subrange be shortened? (logical; default=TRUE)
## Returns: GRanges object containing all subranges of each individual range of input GRanges object
## Details: If "shorten" is TRUE (default), the last subrange of each input range is shortened to the actual end of the range. Unique names are derived from the coordinates (original names are not considered).
subGRanges <- function(gr, binSize=100, shorten=TRUE) {
	# Traverse through all ranges in GRanges objects 'gr', generate subranges and return as GRangesList object
	grl <- lapply(seq_along(gr), function(i) {
				# Set start coordinates of subranges
				head <- as.numeric(seq(start(gr)[i], end(gr)[i], binSize))
				# Set end coordinates of subranges
				tail <- head + binSize - 1
				# Shorten end coordinate of last subrange to actual end of range
				if (shorten) tail[length(tail)] <- end(gr)[i]
				# Construct GRanges object of subranges
				GRanges(seqnames=as.character(seqnames(gr)[i]), IRanges(unlist(head), unlist(tail)), strand=strand(gr)[i], seqlengths = seqlengths(gr))
			})
	# Concatetnate GRangesList to GRanges object
	gr <- do.call(c, grl)
	# Set unique names derived from coordinates
	names(gr) <- paste(seqnames(gr), ":", start(gr), "-", end(gr), ":", strand(gr), sep="")
	# Return GRanges object
	gr
}
###

### C. PREPARE DATA
## List of read files
read_files <- list.files(bam_path, ".bam$", full=TRUE)
read_files_ls <- BamFileList(read_files)
## Import BED and generate subranges
features <- import.bed(feat_file, genome=genome, asRangedData=FALSE)
features <- subGRanges(features, binSize=bin_size)
###

### D. FIND OVERLAPS
overlaps <- summarizeOverlaps(features, read_files_ls, mode=mode, ignore.strand=ignore_strand)
###

### E. SAVE
# Save important objects
save(read_files_ls, features, overlaps, file=paste(out_prefix, "objects.R", sep=""))
# Write features to tab file
write.table(assays(overlaps)$counts, file=paste(out_prefix, "count_table", sep=""), quote=FALSE, sep="\t")
# Write subranges to BED file
export(features, file=paste(out_prefix, "features", sep=""), "bed")
###