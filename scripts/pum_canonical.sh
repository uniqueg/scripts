#!/bin/bash

### Author: Alexander Kanitz
### Created: 25-NOV-2013
### Modified: 25-NOV-2013
### Function: Parses the output of 'structure_at_miRNA_hairpin_position.pl' (Alexander Kanitz) and indicates whether the found motif in column 7 of the output contains a canonical Pumilio binding site (UGUAnAUA) by adding an additional column stating either "canonical" or "non-canonical".
### Arguments: 1. Filename
### Output: Output written to STDOUT
### Usage: bash ./pum_canonical.sh [FILE]

### CHECK FOR PRESENCE OF COMMAND LINE ARGUMENT(S)
if [ $# != 1 ]; then 
	echo Usage: $0 [FILE]
	exit
fi
###

### RUN AWK COMMAND
awk 'BEGIN { FS = "\t"; OFS = "\t" } { if ($7 ~ /..UGUA[ACGU]AUA../) print $0, "canonical"; else print $0, "non-canonical" }' $1
###

### EXIT
exit 0
###
