#!/bin/sh

#=============#
#  HEADER //  #
#=============#
## Created: Sep 3, 2014
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: 
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
Usage: $script [OPTIONS] -- file <BED6>

Description: Sorts a BED6 file of exons by the names (dictionary order), start and stop position (numeric sort order) fields.

Options:
	--temporary-directory=DIR     Directory for saving temporary files (default: $TMP).
	--buffer=SIZE                 Buffer size (default: 2G). Recognizes suffixes % (of total memory), K (for kilobyte), M, G etc. Check GNU sort manual for details.
	--parallel=N                  Use N concurrent threads (default: 4). 
	--usage | --help              Show this screen and exit.
	--version                     Show version information and exit.

Comments:
	- Requires a recent version of GNU sort (developed and tested with version 8.9).

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 3, 2014.
Version 1.0 (Sep 3, 2014)
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
script=$(basename $0)
cat <<VERSION
$script version 1.0 (Sep 3, 2014)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 3, 2014.
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
buffer="2G"
threads=4

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
		verbose=1
		shift
		;;
	--temporary-directory)
		tmp_dir=$$1
		shift 2
		;;
	--buffer-size)
		buffer=$$1
		shift 2
		;;
	--parallel)
		threads=$$1
		shift 2
		;;
	--) # End of all options
		shift
		break
		;;
	-*)
		echo -e "[ERROR] Unknown option: $$1\nExecution aborted." >&2
		usage
		exit 1
		;;
	*)  # no more options. Stop while loop
		break
		;;
	esac
done

#---> VERIFY OPTIONS <---#
# Verify that the specified temporary folder exists
if [ ! -d $tmp_dir ] ; then
	echo -e "[ERROR] Specified output folder '$1' is not a directory!" >&2
	usage
	exit 1	
fi

# Verify that the specified temporary folder is writable
if [ ! -w $tmp_dir ] ; then
	echo -e "[ERROR] Specified output folder '$1' is not writable!" >&2
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
if [ $verbose -eq 1 ]; then echo "Starting '$0'...\n" >&2; fi

#---> MAIN VARIABLES <---#

	#---> BODY <---#

	#---> Sort file <---#
	sort --numeric-sort --field-separator=$'\t' --key=3,3 --temporary-directory=$tmp_dir --buffer-size=$buffer --parallel=$threads $1 | sort --numeric-sort --field-separator=$'\t' --key=2,2 --temporary-directory=$tmp_dir --buffer-size=$buffer --parallel=$threads | sort --dictionary-order --field-separator=$'\t' --key=4,4 --stable  --temporary-directory=$tmp_dir --buffer-size=$buffer --parallel=$threads 

#---> STATUS MESSAGE <---#
if [ $verbose -eq 1 ]; then echo "Done.\n" >&2; fi

#---> PROGRAM EXIT <---#
exit 0;

#===========#
#  // MAIN  #
#===========#
