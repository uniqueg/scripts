#!/bin/bash


# Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# Oct 22, 2015

# DESCRIPTION
# Converts output of CLIPZ tool "mRNA site extraction" to BED-like format. Reads input files from the command line. Output files will have the same basename, but the file extension '.bed'.

# USAGE
# CLIPZ_mRNA_site_extraction_to_BED.sh <FILE> (<FILE> ..)

##################
##  OPTIONS //  ##
##################

# Bash options
set -e
set -o pipefail

##################
##  // OPTIONS  ##
##################


###############
##  MAIN //  ##
###############

# Exit with usage if no arguments
if [ "$#" -eq 0 ]; then
    echo "[ERROR] No arguments given! Execution aborted."
    echo "Usage: CLIPZ_mRNA_site_extraction_to_BED.sh <FILE> (<FILE> ..)"
    exit 1
fi

# Loop over arguments...
for file in "$@"; do

    # Print status message
    echo "Converting file '$file' to BED-like format..."

    # Skip if file is not found
    if [ ! -f "$file" ]; then
        echo "[WARNING] File '$file' not found! File skipped."
        continue
    fi

    # Skip if extension is 'bed'
    ext="${file##*.}"
    if [[ $ext == "bed" ]]; then
    echo "[WARNING] File '$file' ends with extension '.bed'. File skipped."
        continue
    fi

    # Build output filename
    outFile="${file%.*}.bed"

    # Convert file to BED-like format
    awk -v OFS="\t" '{if (NR == 1) {sub(".*/", "", FILENAME); split(FILENAME, array, "_"); id = array[2]} else {if ($2 < 0) $2 = 0; if ($3 >= $6) $3 = $6 - 1; print $1, $2, $3, id":"$1"_"$2"_"$3, $9, "*", $8}}' ${file} | sort --numeric-sort --key 3,3 --stable | sort --numeric-sort --key 2,2 --stable | sort --dictionary-order --key 1,1 --stable > "${outFile}"

    # Print status message
    echo "Output written to file '$outFile'."
done

# Print status message and exit
echo "Done."
exit 0

###############
##  // MAIN  ##
###############
