#!/bin/bash

### Author: Alexander Kanitz
### Created: 20-FEB-2013
### Modified: 20-FEB-2013
### Function: Lifts over hg18 assembly coordinates to assembly hg19
### Arguments: 1. path to liftOver; 2. chain file; 3. output directory; 4. filename (extension MUST be "_hg18.bed"), wildcards accepted
### Output: For each processed input file, a BED file indicating cluster ranges is generated
### Usage: bash ./mz_extended_clusters_to_bed.sh file1.mz [file2.mz]

# Pass arguments
lift_path=$1
chain=$2
out=$3
# Shift so that remaining arguments are input files
shift 4

### LIFTOVER
for file in $*; do
	# Remove extension
	name=${file%\_hg18.bed}
	# Get basename
	name=`basename $name`
	# Liftover
	"${lift_path}/liftOver" $file $chain "$out/${name}_hg19.bed" "$out/${name}_hg19_unmapped"
	# Screen output
	echo "Converted file '$file' to assembly hg19: '$out/${name}_hg19.bed'."
done
###