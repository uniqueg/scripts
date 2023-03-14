#!/bin/bash

## Usage function
usage()
{
																															# Write proper usage function!
		echo "Generate coverage plots etc. from one or more BAM files."
}

## Initialize variables
file=""
verbose=1

## Parse options
while :
do
	case $1 in
	-h | --help | -\?)
		usage
		exit 0
		;;
	-i | --induced)
		induced=$2
		shift 2
		;;
	--induced=*)
		induced=${1#*=}	# Delete everything up till "="
		shift
		;;
	-u | --uninduced)
		uninduced=$2
		shift 2
		;;
	--uninduced=*)
		uninduced=${1#*=}	# Delete everything up till "="
		shift
		;;		
	-p | --parental)
		parental=$2
		shift 2
		;;
	--parental=*)
		parental=${1#*=}	# Delete everything up till "="
		shift
		;;
	-q | --quiet)
		verbose=0
		shift
		;;
	--) # End of all options
		shift
		break
		;;
	-*)
		echo "WARN: Unknown option (ignored): $1" >&2
		shift
		;;
	*)  # no more options. Stop while loop
		break
		;;
	esac
done


## Validate inputs
if [ ! "$induced" ]; then
	echo "[ERROR] Required option '--induced FILE' missing. See --help" >&2
	exit 1
fi

if [ ! "$uninduced" ]; then
	echo "[ERROR] Required option '--uninduced FILE' missing. See --help" >&2
	exit 1
fi

if [ ! "$parental" ]; then
	echo "[ERROR] Required option '--parental FILE' missing. See --help" >&2
	exit 1
fi




#time samtools sort -o $induced '/import/bc2/home/zavolan/kanitz/DDRNA/rna_seq/CAGE/bwa_riken/raw/ind' > '/import/bc2/home/zavolan/kanitz/DDRNA/rna_seq/CAGE/bwa_riken/sorted_indexed/induced_sorted.bam'
