#!/bin/bash

#######
### GENERAL:
### --------
### Author: 	Alexander Kanitz
### Created:	27-NOV-2012
### Modified:	06-DEC-2012
#######

#######
### FUNCTION:
### ---------
### Uses SAMTools package to convert one or more sam files to sorted, indexed bam files.
#######

#######
### ARGUMENTS:
### ----------
### 1. One or more filenames. Wildcards are accepted. File extension MUST be ".sam"
### General Notes: Compare the example in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each processed sam file, a sorted, indexed bam file and the corresponding bam.bai file are generated.
#######

#######
### USAGE:
### ------
### bash /path/to/sam2bam.sh /path/to/files/FILENAME.SAM [/path/to/more/files/FILENAME.SAM]
### Examples:
### 1. bash sam2bam.sh *.sam (processes all sam files in the current folder; script in current folder)
### 2. ./sam2bam *.sam (same as above)
#######

#######
### OVERVIEW:
### ---------
### A. FILE TESTS
### B. SAMTOOLS PROCESSING
#######

### A. FILE TESTS

## A1. CHECK FOR PRESENCE OF COMMAND LINE ARGUMENT(S)
if [ $# == 0 ]; then 
	echo Usage: $0 file1.sam [file2.sam]
	echo Wildcards are allowed.
	exit
fi
##

## A2. CHECK FOR FILE EXTENSION
for file in $*; do
	if [ ! ${file##*.} == "sam" ]; then
		echo Usage: $0 file1.sam [file2.sam]
		echo Wildcards are allowed.
		exit	
	fi
done
###

### B. SAMTOOLS PROCESSING
for file in $*; do
	name=${file%\.sam}
	# Use SAMTools package to convert SAM to sorted BAM
	/import/bc2/soft/bin/samtools view -bS ${file} | /import/bc2/soft/bin/samtools sort - ${name}
	# Use SAMTools package to index BAM (BAI files generated)
	/import/bc2/soft/bin/samtools index "${name}.bam"
done
###
