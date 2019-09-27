#!/bin/sh

#=============#
#  HEADER //  #
#=============#
## Created: Sep 4, 2014
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: tar, gzip, bzip2, zip, samtools, fastq-dump, mktemp
#=============#
#  // HEADER  #
#=============#


#============#
#  USAGE //  #
#============#
usage()
### Returns u information for current script in a string
{
script=$(basename $0)
cat << USAGE
Usage: $script [OPTIONS] -- <FILES/ARCHIVES>

Description: Iteratively tries to extract the passed input files using different extraction tools until flat (i.e. human-readable) files are generated, then moves the flat files to a specified output directory.

Options:
	--usage | --help                        Show this screen and exit.
	--version                               Show version information and exit.
	--verbose                               Increase verbosity. Can be added multiple times.
	--temporary-directory | --tmp-dir DIR   Writable directory for storing the folder holding temporary files (default: current working directory).
	--output-directory | --out-dir DIR      Writable directory for storing extracted files/folders (default: current working directory).
	--max-files                             Maximum number of output temporary files (0 for unlimited; default: 0).

Comments:
	- Requires that all of "tar", "gzip", "bzip2", "zip", "samtools" and "fastq-dump" (SRA-Toolkit) are installed and available in your $PATH.
	- Folders resulting from extraction are currently not recursively processed but rather moved to the output directory as is; folders passed as input files are not allowed
	- If present, terminal file extension (file.EXT) are successively removed from the input filenames, regardless of the nature of the extension
	- IMPORTANT: Existing files are overwritten without prompting. When processing multiple files, filename clashes are not resolved

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 4, 2014.
Version 1.0 (Sep 5, 2014)
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
script=$(basename $0)
cat <<VERSION
$script version 1.0 (Sep 4, 2014)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 4, 2014.
VERSION
}
#============#
#  // USAGE  #
#============#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES <---#
verbose=0
tmp_dir=$PWD
out_dir=$PWD
max_files=0

#---> PARSE / ASSIGN OPTIONS <---#
while :
do
	case $1 in
	--usage | --help)
		usage
		exit 0
		;;
	--version)
		version
		exit 0
		;;
	--verbose)
		verbose=$((verbose+1))
		shift
		;;
	--tmp-dir | --temporary-directory)
		tmp_dir=$2
		shift 2
		;;
	--out-dir | --output-directory)
		out_dir=$2
		shift 2
		;;
	--max-files)
		max_files=$2
		shift 2
		;;
	--) # End of all options
		shift
		break
		;;
	-*)
		echo -e "[ERROR] Unknown option: $1\nExecution aborted." >&2
		usage
		exit 1
		;;
	*)  # no more options. Stop while loop
		break
		;;
	esac
done

#---> PARSE / ASSIGN NON-OPTIONS <---#
# Initialize variable to hold passed filenames
filenames=$@

#---> VERIFY OPTIONS <---#
# Print usage information and exit if input files is not specified
if [[ $filenames == '' ]]; then 
	echo -e "[ERROR] No input files/archives specified!\nExecution aborted." >&2
	usage
	exit 0
fi

# Verify that the specified temporary folder exists 
if [ ! -d $tmp_dir ] ; then 
        echo -e "[ERROR] Specified temporary folder '$tmp_dir' is not a directory!\nExecution aborted." >&2 
        usage 
        exit 1   
fi 
 
# Verify that the specified temporary folder is writable 
if [ ! -w $tmp_dir ] ; then 
        echo -e "[ERROR] Specified temporary folder '$tmp_dir' is not writable!\nExecution aborted." >&2 
        usage 
        exit 1   
fi 

# Verify that the specified temporary folder exists 
if [ ! -d $out_dir ] ; then 
        echo -e "[ERROR] Specified output folder '$out_dir' is not a directory!\nExecution aborted." >&2 
        usage 
        exit 1   
fi 
 
# Verify that the specified temporary folder is writable 
if [ ! -w $out_dir ] ; then 
        echo -e "[ERROR] Specified output folder '$out_dir' is not writable!\nExecution aborted." >&2 
        usage 
        exit 1   
fi 
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#
#  MAIN //  #
#===========#

#---> STATUS MESSAGE <---#
if [[ $verbose -eq 1 ]]; then echo "Starting '$0'..."; fi

#---> MAIN VARIABLES <---#


	#---> BODY <---#

	#---> Generate temporary folder, copy files to be processed to it, read names of files in it into array and determine array size <---#
	# Make temporary directory
	tmp_dir=`mktemp --directory $tmp_dir/tmp.XXXXXXXX`

	# Copy input files to temporary directory 
	cp --recursive $filenames $tmp_dir

	## Put files in temporary directory to array; skip non-files with a warning
	files=()
	for file in $tmp_dir/* ; do
		if [ -f "$file" ]
			then
				files+=("$file")
			else
				rm --recursive --force "$file"
				base=$(basename "$file")
				echo "[WARNING] '$base' is not a file. Skipped." >&2
		fi
	done

	# Determine number of files in array
	numfiles=${#files[@]}

	#---> While files in temporary directory, test whether maximum number of files exceeded, try to extract files until human readable, then move to output folder <---#
	## While the files array is not empty...
	while [ $numfiles -gt 0 ]; do

		## Test whether maximum number of files is exceeded
		if [[ $numfiles -gt $max_files ]] && [[ $max_files -ne 0 ]]; then
			echo -e "[ERROR] Total number of files in archive(s) greater than the maximum number of files specified via --max-files.\nExecution aborted." >&2
			rm -rf $tmp_dir
			exit 1
		fi

		#---> Iterate over all files <---#
		for file in "${files[@]}"; do
	   		
	   		#---> Get basename <---
	   		base=$(basename "$file")
	   		
			#---> Check whether file is human readable; move to output directory <---
			type=$(file --brief --mime-type "$file")
			if [[ $type == "text/plain" ]]; then
				mv "$file" $out_dir
				if [[ $verbose -eq 1 ]]; then echo "Moved file '$base' to output directory."; fi
				continue
			fi

			#---> Check whether file is empty; issue warning and remove <---
			type=$(file --brief --mime-type "$file")
			if [[ $type == "application/x-empty" ]]; then
				rm "$file"
				echo "[WARNING] File '$base' is empty. Skipped." >&2
				continue
			fi	

			#---> Check whether file is a directory; move to output directory <---
			if [ -d "$file" ]; then
				cp --recursive "$file" $out_dir
				rm --recursive --force "$file"
				if [[ $verbose -eq 1 ]]; then echo "Moved directory '$file' to output directory."; fi
				continue
			fi

			#---> Get file suffix and build output filename <---
			# Create random suffix if filename does not contain a suffix of format file.SUFFIX
			if [[ ! $base == *.* ]] ; then
				tmp_file=`mktemp $file.XXXXX`
				mv $file $tmp_file
				file=$tmp_file
			fi			

			# Get terminal file suffix (file.SUFFIX)
			suffix="${file##*.}"

			# Build output filename by removing suffix from filename
			out_file=$tmp_dir/$(basename "$file" .$suffix)			

			# Remove output file if it already exists
			rm --recursive --force $out_file
		
			#---> Try to extract file using different methods; if file can be extracted, remove it <---#

                        ## BAM files
                        samtools view -h "$file" > "$out_file.sam" 2> /dev/null
                        if [ $? -eq 0 ]
                                then
                                        rm --force "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'samtools'. Output filename: '$out_file.sam'"; fi
                                        continue
                                else
                                        rm --force --recursive $out_file.sam
                        fi

                        ## SRA files
                        fastq-dump --stdout "$file" > "$out_file.fq" 2> /dev/null
                        if [ $? -eq 0 ]
                                then
                                        rm --force "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'fastq-dump'. Output filename: '$out_file.fq'"; fi
                                        continue
                                else
                                        rm --force --recursive $out_file.fq
                        fi

                        ## ZIP files
                        unzip -qq -o -d $tmp_dir "$file" 2> /dev/null
                        if [ $? -eq 0 ]
                                then
                                        rm --recursive --force "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'unzip'."; fi
                                continue
                        fi

                        ## BZIP2 files
                        bzip2 --decompress --quiet --stdout "$file" > $out_file 2> /dev/null
                        if [ $? -eq 0 ]
                                then
                                        rm --force --recursive "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'bzip2'. Output filename: '$out_file'"; fi
                                        continue
                                else
                                        rm --force --recursive $out_file
                        fi

                        ## GZIP files
                        gzip --decompress --quiet --suffix $suffix --stdout "$file" > $out_file 2> /dev/null
                        if [ $? -eq 0 ]
                                then
                                        rm --force --recursive "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'gzip'. Output filename: '$out_file'"; fi
                                        continue
                                else
                                        rm --force --recursive $out_file
                        fi

                        ## TAR files
                        tar -C $tmp_dir --overwrite -xf "$file" > /dev/null 2>&1
                        if [ $? -eq 0 ]
                                then
                                        rm --recursive --force "$file"
                                        if [[ $verbose -gt 1 ]]; then echo "Extracted file '$base' with 'tar'."; fi
                                        continue
                        fi
	
		done
	
		#---> Read names of files in temporary folder into array, move directories to output folder and determine array size <---#
		## Put files in temporary directory to array; move directories to output directory
		files=()
		shopt -s nullglob
		for file in $tmp_dir/* ; do
			files+=("$file")
		done

		## Determine number of files in array
		numfiles=${#files[@]}

	done
	
        ## Remove temporary folder
        rm --recursive --force $tmp_dir

#---> STATUS MESSAGE <---#
if [[ $verbose -eq 1 ]]; then echo "Done."; fi

#---> PROGRAM EXIT <---#
exit 0

#===========#
#  // MAIN  #
#===========#


#=====================#
#  FUTURE OPTIONS //  #
#=====================#
	#	--recursive	Process files within folders
	#	--no-renaming	Do not rename files upon extraction; alternatively: only remove *matching* extensions (gz for gzip etc.)
#=====================#
#  // FUTURE OPTIONS  #
#=====================#
