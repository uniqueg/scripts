#!/bin/bash

#==================#
#   HEADER START   #
#==================#
### Name: sam_cat_sort.sh
### Created: Aug 28, 2013
### Modified: Aug 28, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: Concatenates a number of SAM files, then writes out the resulting file sorted by read name (QNAME). Header lines are NOT ignored, so for the most stable/safe experience, strip headers before usage.
### Notes:
### - Requires a recent version of GNU sort (tested with version 8.9)
### - Sorting is performed in "version" mode (refer to GNU sort manual) via the first field (tab-separated)
### - Default resources allocated for sorting: 16G buffer size, 8 threads (change in PRE-REQUISITES section)
### - Default temporary directory: $TMP (change below) (change in PRE-REQUISITES section)
### Output: Writes concatenated, QNAME-sorted SAM file to STDOUT
### Usage: bash ./sam_cat_sort.sh [FILE/S]
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#

###########################################################
# Define variables that change the behavior of the sorting
buffer="16G"
threads=8
tmp_dir=$TMP
###########################################################

# Define variable containing usage information
usage="bash ./$0 input_sam_file/s [FILE/S]"

## Die if wrong number of arguments
if [ $# = 0 ]; then
	echo -e "[ERROR] No input files!\n\nUsage: $usage\n" >&2
	exit 1
fi

## Die if input files do not exist
for file in $*; do
	if [ ! -f $file ]; then
		echo -e "[ERROR] File '$file' not found!\n\nUsage: $usage\n" >&2
		exit 1
	fi
done

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#

# Concatenate and sort input files
cat $* | sort --sort=version --field-separator=$'\t' --key=1 --stable --buffer-size=$buffer --parallel=$threads --temporary-directory=$tmp_dir

# Exit script
exit 0

#================#
#    MAIN END    #
#================#
