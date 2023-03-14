#!/bin/bash

### Author: Alexander Kanitz
### Created: 12-FEB-2013
### Modified: 25-FEB-2013
### Adapted from: sam2bam.sh
### Function: Summarizes and converts output from C-script "generate_extended_clusters" (written by Mihaela Zavolan) to BED files indicating cluster ranges
### Arguments: 1. One or more filenames. Wildcards are accepted. File extension MUST be ".mz"
### Output: For each processed input file, a BED file indicating cluster ranges is generated
### Usage: bash ./mz_extended_clusters_to_bed.sh file1.mz [file2.mz]

### CHECK FOR PRESENCE OF COMMAND LINE ARGUMENT(S)
if [ $# == 0 ]; then 
	echo Usage: $0 file1.mz [file2.mz]
	echo Wildcards are allowed.
	exit
fi
###

### CHECK FOR FILE EXTENSION
for file in $*; do
	if [ ! ${file##*.} == "mz" ]; then
		echo Usage: $0 file1.mz [file2.mz]
		echo Wildcards are allowed.
		exit	
	fi
done
###

### B. PROCESSING
for file in $*; do
	# Get basename
	name=${file%\.mz}
	# Extract summary lines using grep, rearrange using awk
	grep ^"Cluster: " "${name}.mz"  | awk -F'[ .]' '{print $3 "\t" $4 "\t" $6 "\t" $2 "\t" ($6 - $4)}' > "${name}.bed" 
	# Screen output
	echo "File '${name}.mz' converted to '${name}.bed'."
done
###