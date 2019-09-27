#!/bin/bash

# Alexander Kanitz
# 04-JAN-2015

# [USAGE]
# scanAdapter.sh <FASTQ> (<NUMBER_OF_SEQUENCES> <ADAPTER_FILE>)

# [DESCRIPTION]
# The script scans a part of a specified FASTQ file (gzipped or uncompressed) for the presence of 
# adapter sequences present in a separate file. Optionally, the number of sequences to scan 
# (default: 10,000) and the location of the adapter file (default: 
# /scicore/home/zavolan/kanitz/resources/adapters/common_seq_adapters.grep) can optionally be 
# modified.

seqFile=$1
seqNumber=$((${2:-10000}*4))
adapterFile=${3:-"/scicore/home/zavolan/kanitz/RESOURCES/adapters/common_seq_adapters.grep"}

while read line; do
    echo -n "${line}: "
    grep -c "$line" <(zcat -f -- "$seqFile" | head -n $seqNumber)
done < "$adapterFile"
