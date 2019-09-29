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
### Returns information for current script in a string
{
script=$(basename $0)
cat << USAGE
Usage: $script [OPTIONS] -- file <BED6>

Description: Sorts a BED6 file of exons by the names (dictionary order), start and stop position (numeric sort order) fields.

Options:
	--temporary-directory DIR     Directory for saving temporary files (default: current working directory).
	--buffer-size SIZE            Buffer size (default: 2G). Recognizes suffixes % (of total memory), K (for kilobyte), M, G etc. Check GNU sort manual for details.
	--parallel N                  Use N concurrent threads (default: 4).
        --verbose                     Print additional information to STDERR during execution.
	--usage | --help              Show this screen and exit.
	--version                     Show version information and exit.

Comments:
	- Requires a recent version of GNU sort (developed and tested with version 8.9).

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 3, 2014.
Version 1.1 (Sep 16, 2014)
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
script=$(basename $0)
cat <<VERSION
$script version 1.1 (Sep 16, 2014)
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
		usage >&2
		exit 0
		;;
	--version)
		version >&2
		exit 0
		;;
	--verbose)
		verbose=1
		shift
		;;
	--temporary-directory | --tmp-dir)
		tmp_dir=$2
		shift 2
		;;
	--buffer-size | --buffer)
		buffer=$2
		shift 2
		;;
	--parallel | --threads)
		threads=$2
		shift 2
		;;
	--) # End of all options
		shift
		break
		;;
	-*)
		echo -e "[ERROR] Unknown option: $$1\nExecution aborted!\n" >&2
		usage >&2
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
	echo -e "[ERROR] Specified temporary folder '$tmp_dir' is not a directory.\nExecution aborted!\n" >&2
	usage >&2
	exit 1	
fi

# Verify that the specified temporary folder is writable
if [ ! -w $tmp_dir ] ; then
	echo -e "[ERROR] Specified temporary folder '$tmp_dir' is not writable.\nExecution aborted!\n" >&2
	usage >&2
	exit 1	
fi
## Verify that exactly one input file is specified
if [ ! $# = 1 ]; then 
        echo -e "[ERROR] Wrong number of input files! Specify exactly one BED6 file.\nExecution aborted!\n" >&2
        usage >&2
	exit 1
fi
## Verify that input file exists
for file in $*; do
        if [ ! -f $file ]; then
                echo -e "[ERROR] File '$file' not found.\nExecution aborted.\n" >&2
		usage >&2
                exit 1
        fi
done
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

