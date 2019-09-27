#!/usr/bin/env Rscript

## Alexander Kanitz, Biozentrum, University of Basel
## alexander.kanitz@alumni.ethz.ch
## 27-OCT-2015

## DESCRIPTION
# Given an expression table of the format '<feature ID> TAB <expression value>' (exactly two 
# columns!) and a precompiled category subset R object, builds a matrix of features (rows) and 
# category subsets (columns), indicating for each field either the expression value from the 
# expression table (if the feature is included in the category subset), or NA (if the feature is 
# *not* included in the category subset).

## NOTES
# Since this script was written for a specific purpose, several assumptions are made that influence 
# the behavior in a way that is likely not desired for generic execution. In particular, the 
# category subset object is a binary R object assumed to be a list of lists of character vectors. 
# The main container is assumed to be called 'subsets', and each element of the main container 
# represents a feature type (e.g. genes, transcripts etc.). Each feature type list element 
# represents a sublist of category subset vectors, each of which indicating the IDs of features 
# belonging to the subset. The desired feature type (i.e. the actual name of the main container 
# element) has to be specified. Note that it is assumed that all entries (i.e. feature IDs) in all 
# subset vectors of the given feature type are present in the <feature ID> column of the specified 
# expressioin table. In addition to these format requirements, note that the first element of each 
# feature type sublist (i.e. the first vector) is excluded from processing. If this behavior is not 
# desired, comment out the line 'subsets <- subsets[-1]'. Finally, complex assumptions are made 
# about the names of each category subset vector. If you obtain unexpected category names in the 
# headers of output files, modify the 'formatNames' functions according to your needs.

## USAGE
# buildTruthTableWithCategories.R <TRUTH_TABLE> <CATEGORY_OBJECT> <FEATURE_TYPE> <OUT_FILE>

## COMMAND-LINE ARGUMENTS ##
args <- commandArgs(trailingOnly = TRUE)
truth=args[1]
categories=args[2]
type=args[3]
outFile <- args[4]

## FUNCTIONS ##
formatNames <- function(names) {
        names <- gsub("--op--A--cl--", "(A)", names)
	names <- gsub("_--op--.*--cl--", "", names)
	names <- gsub("--col--", "", names)
	names <- gsub("_to_", "-", names)
	names <- gsub("--dot--", ".", names)
	names <- gsub("_", " ", names)
	return(names)
}

## MAIN ##
# Load truth table
df <- read.table(truth, header=FALSE, col.names=c("featureID", "allFeatures"), row.names=1)
vec <- setNames(df[,1], row.names(df))
# Load category subset object
load(categories)
subsets <- subsets[[type]]
subsets <- subsets[-1]
names(subsets) <- formatNames(names(subsets))
# Get list of expression value/NA vectors for each category/subset
ls <- lapply(subsets, function(subset) {
	ifelse(names(vec) %in% subset, vec, NA)
})
# Build table
df <- as.data.frame(do.call(cbind, ls))
df <- cbind(featureID=names(vec), df)
# Write output
write.table(df, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
