#!/bin/bash

### Author: Alexander Kanitz
### Created: 25-FEB-2013
### Modified: 25-FEB-2013
### Adapted from: bed_count_clusters.sh
### Function: Extracts information from overlap .bed files generated by R script "subset_overlaps_bed_bed_overlap_no.R" run on files generated by shell script "mz_extended_clusters_to_bed" after C-script "generate_extended_clusters" (written by Mihaela Zavolan); scripts works properly ONLY if processing file names of the format "..._WINSIZE.." (i.e. the second element of an underscore-delimited file name should be an integer representing the window size used for clustering) 
### Arguments: 1. output file name; 2. one or more filenames with extension ".bed" (wildcards are accepted)
### Output: One output file in tab format is generated for all files in argument 2 (one row per file); it lists the window size used for clustering (extracted from the file names!), the number of clusters containing an AsiSI motif, and the number of occurrences of the AsiSI motif (GCGATCGC) within each cluster
### Usage: bash ./mz_extended_clusters_to_bed.sh ./out_file ./in_file.bed [./in_file2.bed]

### Pass output file variable and shift
out_file=$1
shift 1
###

### Check file extension
for file in $*; do
        if [ ! ${file##*.} == "bed" ]; then
                echo Usage: $0 file1.mz [file2.mz]
                echo Wildcards are allowed.
                exit
        fi
done
###

### File processing
for file in $*; do
        # Print window size used for clustering to file 
        echo $file | awk -F'/' '{ print $NF }' | awk -v ORS="\t" -v out=$out_file -F'_' '{ print $2 >> out }'
        # Print number of clusters to file
        wc -l $file | awk -v ORS="\t" -v out=$out_file '{ print $1 >> out }'
        # Print number of AsiSI sites within each cluster to file 
        cat $file | awk -v ORS="\t" -v out=$out_file '{ print $5 >> out }'
        # Print newline character
        echo "" >> $out_file
        # Screen output
        echo "File '$file' processed."
done
###