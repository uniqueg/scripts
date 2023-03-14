#!/bin/bash

### Author: Alexander Kanitz
### Created: 01-FEB-2013
### Modified: 01-FEB-2013
### Adapted from: N/A 
### Description: The script 'greps' lines containing the indicated regex from a specified subset of files in a specified folder
### Arguments: 1. regex to look for, 2. (path to) files to consider (may include wildcards), 3. --ignore-case 
### Output: File with lines from all input files containing regex
### Usage: bash ./lines_with_regex_from_files.sh ACGT ./*.fasta no ./out  

### A. TEST ARGUMENTS

## A1. CHECK FOR PRESENCE OF COMMAND LINE ARGUMENTS
if [ $# -lt 2 ] || [ $# -gt 3 ]; then 
	echo "Invalid number of arguments!"
	echo $#
	echo "Usage: $0 <regex> <path/to/files(wildcards allowed)> <optional:--ignore-case>"
	exit
fi
##

## A2. PROCESS ARGUMENTS
regex=
case_sensitive=$3
##

## A3. CHECK FOR EXISTENCE OF INPUT FILE DIRECTORY
if [ ! -d "`dirname $2`" ]; then
	echo "Input file folder does not exist!"
	echo "Usage: $0 <regex> <path/to/files(wildcards allowed)> <--ignore-case>"
	exit
fi
##

###


### B. EXECUTION 
grep $case_sensitive --no-filename $1 $2
###