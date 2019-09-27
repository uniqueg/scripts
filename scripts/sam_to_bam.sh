#!/bin/bash

#==================#
#   HEADER START   #
#==================#
### Name: sam2bam.sh
### Created: Nov 27, 2012
### Modified: Mar 27, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.2
### Adapted from: n/a
### Requirements: SAMTools
#==================#
### Description: Uses SAMTools package to generate sorted, indexed BAM files from one or more SAM files 
### Arguments: 1. Filename(s) with suffix(es) "sam", wildcards accepted
### Output: For each processed SAM file, a sorted BAM file and the corresponding ".bam.bai" index file are generated
### Usage: bash ./sam2bam.sh [SAM FILE 1] [SAM FILE 2] ...
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
## CHECK FOR PRESENCE OF COMMAND LINE ARGUMENT(S)
if [ $# == 0 ]; then 
	echo Usage: $0 file1.sam [file2.sam]
	echo Wildcards are allowed.
	exit
fi

## CHECK FOR FILE EXTENSION
for file in $*; do
	if [ ! ${file##*.} == "sam" ]; then
		echo Usage: $0 file1.sam [file2.sam]
		echo Wildcards are allowed.
		exit	
	fi
done
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
### B. SAMTOOLS PROCESSING
for file in $*; do
	# Get basename of SAM file
	name=${file%\.sam}
	# Use SAMTools package to convert SAM to sorted BAM
	/import/bc2/soft/bin/samtools view -bS ${file} | /import/bc2/soft/bin/samtools sort - ${name}
	# Use SAMTools package to index BAM (BAI files generated)
	/import/bc2/soft/bin/samtools index "${name}.bam"
	# Status message
	echo "Processed file '$file'."
done
#================#
#    MAIN END    #
#================#
