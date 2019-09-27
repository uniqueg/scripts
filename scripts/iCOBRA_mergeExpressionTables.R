#!/usr/bin/env Rscript

## Alexander Kanitz, Biozentrum, University of Basel
## alexander.kanitz@alumni.ethz.ch
## 27-OCT-2015

## DESCRIPTION
# Merges expression table files of the format '<feature ID> TAB <expression value>' (exactly two 
# columns!) present in a specified directory.

## NOTES
# Since this script was written for a specific purpose, several assumptions are made that influence 
# the behavior in a way that is likely not desired for generic execution. In particular, the scripts 
# assumes that all files in the specified directory are expression tables of the described type. 
# Moreover, files containing 'Ground_truth' or 'A-seq-2' in their respective filenames are excluded. 
# The output file contains the feature names in the first column, and the expression values 
# extracted from each file in the remaining columns. A header line is included, consisting of the
# literal 'featureID' for the first column and the basenames of each file up to the first dot for 
# the remaining columns.

## USAGE
# mergeExpressionTables.R <INPUT_DIRECTORY> <OUTPUT_FILE>

## COMMAND-LINE ARGUMENTS ##
args <- commandArgs(trailingOnly = TRUE)
inDir <- args[1]
outFile <- args[2]

## FUNCTIONS ##
merge_df_from_ls <- function(ls, by=0, all=TRUE) {
        mdf <- ls[[1]]
        for (df in 2:length(ls)) {
                mdf <- merge(mdf, ls[[df]], by=by, all=all)
        }
        return(mdf)
}

## MAIN ##
# Load files
files <- dir(inDir, full.names=TRUE)
files <- files[grep("Ground_truth", files, invert=TRUE)]
files <- files[grep("A-seq-2", files, invert=TRUE)]
df_ls <- lapply(files, function(file) {
    read.table(file, col.names=c("featureID", unlist(strsplit(basename(file), ".", fixed=TRUE))[1]), sep="\t", header=FALSE)
})
# Merge tables
df <- merge_df_from_ls(df_ls, by=1)
# Write output
write.table(df, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
