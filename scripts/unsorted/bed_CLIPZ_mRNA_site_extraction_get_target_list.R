#!/usr/bin/env Rscript

# Alexander Kanitz, Biozentrum, University of Basel
# 05-AUG-2015

# DESCRIPTION
# Build CLIP target list from multiple replicates based on the output of the CLIPZ tool "mRNA site 
# extraction", converted to BED-like format. Reads input files from the command line. An output 
# table ("sites.csv") and a Venn diagram ("venn.tiff") depicting the overlap of sites between 
# replicates are written to the specified output directory.

# USAGE
# Usage: bed_CLIPZ_mRNA_site_extraction_get_target_list.R <BED> (<BED_2> ... <BED_N>) <OUTPUT_DIRECTORY>

# DEPENDENCIES
# The following packages are required:
# - VennDiagram
# - Bioconductor
# - GenomicRanges


###################################
##  OPTIONS & INITIALIZATION //  ##
###################################

# Set parameters
col_selector <- c(4,5,7)
col_selector_merge <- 1:4

# Parse command-line arguments
write("Parsing command-line arguments...", stdout())
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 2 ) {
    write("[ERROR] Not enough arguments! Execution aborted.", stderr())
    write("Usage: bed_CLIPZ_mRNA_site_extraction_get_target_list.R <BED> (<BED_2> ... <BED_N>) <OUTPUT_DIRECTORY>", stderr())
    quit(save="no", status=1)
}
inFiles <- args[1:length(args)-1]
outDir <- args[length(args)]

# Load packages
write("Loading packages...", stdout())
if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicRanges"))) == FALSE  ) { stop("Package 'GenomicRanges' required!\nExecution aborted.")  }
if ( suppressWarnings(suppressPackageStartupMessages(require("VennDiagram"))) == FALSE  ) { stop("Package 'VennDiagram' required!\nExecution aborted.")  }

###################################
##  // OPTIONS & INITIALIZATION  ##
###################################


####################
##  FUNCTIONS //  ##
####################

# Functions
Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

####################
##  // FUNCTIONS  ##
####################


###############
##  MAIN //  ##
###############

# Load files
write("Loading files...", stdout())
df_ls <- lapply(inFiles, function(file) {
    df <- read.table(file, header=FALSE, sep="\t", colClasses=c("character", "integer", "integer", "character", "numeric", "character", "numeric"), col.names=c("refSeq", "start", "end", "siteID", "foldEnrichment", "strand", "logPosterior"))
})

# Build replicate names
rep_names <- paste("rep", 1:length(inFiles), sep="_")

# Concatenate and merge overlapping sites
write("Concatenating and merging overlapping sites...", stdout())
merge_df <- do.call(rbind, df_ls)
merge_gr <- with(merge_df, GRanges(refSeq, IRanges(start+1, end), strand=strand, id=siteID, score=foldEnrichment, logPosterior=logPosterior))
merge_gr <- reduce(merge_gr)

# Convert site/region tables to GRanges objects
gr_ls <- lapply(df_ls, function(file) {
    query <- with(file, GRanges(refSeq, IRanges(start+1, end), strand=strand, id=siteID, score=foldEnrichment, logPosterior=logPosterior))
})

# Calculate overlaps
write("Calculating overlaps between replicates...", stdout())
ol_ls <- lapply(gr_ls, findOverlaps, subject=merge_gr, type="within")

# Evaluate overlaps
write("Evaluating overlaps...", stdout())
ol_filled_ls <- mapply(function(hits, df) {
    subset_vec <- match(1:subjectLength(hits), subjectHits(hits))
    df[subset_vec, col_selector]
}, ol_ls, df_ls, SIMPLIFY=FALSE)
ol_df <- do.call(cbind, ol_filled_ls)
colnames(ol_df) <- paste(colnames(ol_df), rep_names, sep=".")
rownames(ol_df) <- NULL
mergedRegions <- as.data.frame(merge_gr)[, col_selector_merge]

# Calculate summary statistics
write("Calculating summary statistics...", stdout())
replicateHits <- length(df_ls) - rowSums(is.na(ol_df)) / length(col_selector)
foldEnrichmentMean <- rowMeans(ol_df[, seq(2, ncol(ol_df), length(col_selector)), drop=FALSE], na.rm=TRUE)
foldEnrichmentStDev <- apply(ol_df[, seq(2, ncol(ol_df), length(col_selector)), drop=FALSE], 1, sd, na.rm=TRUE)
logPost_ls <- split(ol_df[, seq(3, ncol(ol_df), length(col_selector)), drop=FALSE], rownames(ol_df))
logPost_ls <- logPost_ls[order(as.numeric(names(logPost_ls)))]
chiSquare <- sapply(logPost_ls, function(logPost) -2 * sum(log(-logPost[! is.na(logPost)])))
pValue <- mapply(function(logPost, chiSq) pchisq(chiSq, df=2 * length(logPost), lower.tail=FALSE), logPost_ls, chiSquare)

# Compile & sort output table
write("Compiling and sorting output table...", stdout())
result_df <- cbind(mergedRegions, replicateHits, foldEnrichmentMean, foldEnrichmentStDev, chiSquare, pValue, ol_df)
result_df <- result_df[order(result_df$pValue),]

# Extract information for Venn diagram & plot
write("Plotting Venn diagram...", stdout())
venn_df <- ol_df[,seq(1, ncol(ol_df), 3), drop=FALSE]
colnames(venn_df) <- rep_names
venn_ls <- lapply(as.list(as.data.frame(ifelse(is.na(venn_df), NA, 1:nrow(venn_df)))), na.omit)
venn.diagram(venn_ls, file.path(outDir, "venn.tiff"))

# Write table
write("Writing output table...", stdout())
write.table(result_df, file.path(outDir, "sites.csv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Save session
write("Saving session...", stdout())
save.image(file.path(outDir, "session.R"))

# Print status message & exit
write("Done.", stdout())
quit(save="no", status=0)

###############
##  // MAIN  ##
###############
