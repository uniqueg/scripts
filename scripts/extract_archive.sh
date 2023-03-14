#!/bin/bash

#==================#
#   HEADER START   #
#==================#
### Name: extract_archive.sh
### Created: Feb 13, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: Bash (not POSIX compliant)
#==================#
#    HEADER END    #
#==================#

# Define variable containing usage information
script=`basename $0`
usage="bash ./$script output_folder [PATH] archive_1 (archive_2 ... archive_n)"
descr='Description:\nExtracts archives with file extensions .bz2, .gz, .tar, .tar.bz2, .tar.gz, .tbz2, .tgz and .zip into the specified (existing and writable!) directory. CAREFUL: Existing files are overwritten without prompting.\n\nRequires: Bash, GNU tar, GNU gunzip, GNU bunzip2, GNU unzip\n\nWritten on 13-FEB-2014 by Alexander Kanitz, Zavolan Group, Biozentrum, University of Basel\nVersion 1.0 (13-FEB-2014)'

## Set default exit status
exit=0

## Die if wrong number of arguments
if [[ $# -lt 2 ]]; then
	echo -e "[ERROR] Insufficient number of arguments!\n\nUsage: $usage\n\n${descr}" >&2
	exit 2
fi

## Verify that the output folder exists and is writable
if [ ! -d $1 ] ; then
	echo -e "[ERROR] Specified output folder '$1' is not a directory!\n\nUsage: $usage\n\n${descr}" >&2
	exit 2	
fi
if [ ! -w $1 ] ; then
	echo -e "[ERROR] Specified output folder '$1' is not writable!\n\nUsage: $usage\n\n${descr}" >&2
	exit 2	
fi

## Set output folder and remove argument from args array
outdir=$1
shift

## For each file...
for file in $*; do

	## If file exists...
	if [ -f "$file" ] ; then

		## ...check file extension and extract accordingly
		case "$file" in
			*.tar.bz2)	tar -C $outdir --overwrite -xjf "$file"
						;;
			*.tar.gz)	tar -C $outdir --overwrite -xzf "$file"
						;;
			*.tar)		tar -C $outdir --overwrite -xf "$file"
						;;
			*.tbz2)		tar -C $outdir --overwrite -xjf "$file"
						;;
			*.tgz)		tar -C $outdir --overwrite -xzf "$file"
						;;
			*.bz2)		base=$(basename "$file" .bz2)
						bunzip2 --quiet --force --stdout "$file" > $outdir/"$base"
						;;
			*.gz)		base=$(basename "$file" .gz)
						gunzip --quiet --force --stdout "$file" > $outdir/"$base"
						;;
			*.zip)		unzip -qq -o -d $outdir "$file"
						;;
	    	*)			echo -e "[WARNING] Unknown file extension. File '$file' could not be extracted! Continuing with other files (if any) but setting exit status to 1." >&2
						exit=1
						;;
		esac

	## ...else exit with error
	else
		echo -e "[WARNING] '$file' is not a file! Continuing with other files (if any) but setting exit status to 1." >&2
		exit=1
	fi

done

## Return exit status
exit $exit