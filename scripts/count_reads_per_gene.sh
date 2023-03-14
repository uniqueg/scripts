#!/bin/sh

#=============#
#  HEADER //  #
#=============#
## Created: Oct 9, 2014
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: GNU core utilities, samtools, bedtools, perl, bam_bed12_bedtools_intersect_to_count_table.pl
#=============#
#  // HEADER  #
#=============#


#============#
#  USAGE //  #
#============#
usage()
### Returns usage information for current script in a string
{
script=$(basename $0)
cat << USAGE
Usage: $script [OPTIONS] --bam --bed12 --count_table

Description: Given a BAM alignment and a BED12 annotation/feature file, generates a count table of reads per gene.

Options:
        --bam=PATH                   Path to the BAM read alignments file (required).
        --bed12=PATH                 Path to the BED12 annotation/feature input file (required).
        --count-table=PATH           Output filename (required).
        --temporary-directory=DIR    Directory in which a folder for temporary files is to be generated (default: $TMP). The generated folder and its contents will be deleted upon completion of the script.
        --keep-temporary-directory   Do not delete temporary folder and files at the end of the script.
        --bam-sorted                 The input BAM file is already sorted by read name (QNAME field).
        --strand=STRING              Strand that shall be considered when computing overlaps. One of "sense", "antisense" and "either" (default: sense).
        --minimum-overlap=FLOAT      Fraction of a read that needs to intersect a feature for an overlap to be reported (number between 1E-9 and 1; default: 1, i.e. the whole read).
        --allow-multi-sites          Whether reads aligning to multiple loci/sites ('multimappers') and features shall be considered for counting. Reads will be split up between sites/features (may lead to non-integer values in the output count table).
        --allow-overlaps             Whether reads aligning to multiple features at a single locus/site (i.e. overlapping features) shall be considered for counting. Readss will be split up between features (may lead to non-integer values in the output count table).
        --bam-sort-threads=INT       Use N concurrent threads (default: 4).     
        --bam-sort-memory=INT        Buffer size (default: "2G"). Recognizes suffixes "%" (of total memory), "K" (for kilobyte), "M", "G" etc. Check GNU and samtools sort manuals for details.
        --samtools=PATH              Path to 'samtools' executable.
        --bedtools=PATH              Path to 'bedtools' executable.
        --perl=PATH                  Path to 'perl' executable (default: "/usr/bin/perl").
        --count-script=PATH          Path to 'bam_bed12_bedtools_intersect_to_count_table.pl' script.
        --verbose                    Write status/progress messages to STDERR.
        --usage | --help             Show this screen and exit.
        --version                    Show version information and exit.

Comments:
        - Only marginal validation of user input is performed within this wrapper script. Make sure to check all your inputs before running!

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Oct 9, 2014.
Version 1.0 (Oct 9, 2014)
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
script=$(basename $0)
cat <<VERSION
$script version 1.0 (Oct 9, 2014)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Oct 9, 2014.
VERSION
}
#============#
#  // USAGE  #
#============#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES / SET DEFAULTS <---#
tmp_dir=$TMP
keep_tmp=0
sorted=0
strand="sense"
overlap=1
allow_ms=""
allow_ol=""
threads=4
memory="2G"
samtools="/import/bc2/home/zavolan/krini/soft/samtools-0.1.19/samtools"
bedtools="/import/bc2/home/zavolan/krini/soft/bedtools2/bin/bedtools"
perl="/usr/bin/perl"
count_script="/import/bc2/home/zavolan/GROUP/ISO_BENCH/scripts/bam_bed12_bedtools_intersect_to_count_table.pl"

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
        --bam)
                bam=$2
                shift 2
                ;;
        --bed12)
                bed12=$2
                shift 2
                ;;
        --count-table)
                count_table=$2
                shift 2
                ;;
        --temporary-directory)
                tmp_dir=$2
                shift 2
                ;;
        --keep-temporary-directory)
                keep_tmp=1
                shift
                ;;
        --bam-sorted)
                sorted=1
                shift
                ;;
        --strand)
                strand=$2
                shift 2
                ;;
        --minimum-overlap)
                overlap=$2
                shift 2
                ;;
        --allow-multi-sites)
                allow_ms="--allow-multi-sites"
                shift
                ;;
        --allow-overlaps)
                allow_ol="--allow-overlaps"
                shift
                ;;
        --bam-sort-threads)
                threads=$2
                shift 2
                ;;
        --bam-sort-memory)
                memory=$2
                shift 2
                ;;
        --samtools)
                samtools=$2
                shift 2
                ;;
        --bedtools)
                bedtools=$2
                shift 2
                ;;
        --perl)
                perl=$2
                shift 2
                ;;
        --count-script)
                count_script=$2
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
# Verify that required arguments are present
if [ -z "$bam" ]; then echo -e "[ERROR] Specify the --bam option with the path to a valid BAM file as its argument!" >&2; usage; exit 1; fi
if [ -z "$bed12" ]; then echo -e "[ERROR] Specify the --bed12 option with the path to a valid BED12 file as its argument!" >&2; usage; exit 1; fi
if [ -z "$count_table" ]; then echo -e "[ERROR] Specify the --count-table option with the path to a writable location!" >&2; usage; exit 1; fi

# Verify that the specified files exist
if [ ! -f $bam ] ; then echo -e "[ERROR] Specified argument '$bam' to --bam does not exist is not a file!" >&2; usage; exit 1; fi
if [ ! -f $bed12 ] ; then echo -e "[ERROR] Specified argument '$bed12' to --bed12 does not exist or is not a file!" >&2; usage; exit 1; fi
if [ ! -f $samtools ] ; then echo -e "[ERROR] Specified argument '$samtools' to --samtools is not a file!" >&2; usage; exit 1; fi
if [ ! -f $bedtools ] ; then echo -e "[ERROR] Specified argument '$bedtools' to --bedtools is not a file!" >&2; usage; exit 1; fi
if [ ! -f $perl ] ; then echo -e "[ERROR] Specified argument '$perl' to --perl is not a file!" >&2; usage; exit 1; fi
if [ ! -f $count_script ] ; then echo -e "[ERROR] Specified argument '$count_script' to --count-script is not a file!" >&2; usage; exit 1; fi

# Verify that the specified files are readable
if [ ! -r $bam ] ; then echo -e "[ERROR] Specified argument '$bam' to --bam is not readable!" >&2; usage; exit 1; fi
if [ ! -r $bed12 ] ; then echo -e "[ERROR] Specified argument '$bed12' to --bed12 is not readable!" >&2; usage; exit 1; fi

# Verify that the specified folder exists
if [ ! -d $tmp_dir ] ; then echo -e "[ERROR] Specified argument '$tmp_dir' to --temporary-directory does not exist or is not a folder!" >&2; usage; exit 1; fi

# Verify that the specified folder is writable
if [ ! -w $tmp_dir ] ; then echo -e "[ERROR] Specified argument '$tmp_dir' to --temporary-directory is not writable!" >&2; usage; exit 1; fi

# Verify that the value of the specified string is allowed
if   [[ "$strand" == "sense" ]]; then strand="-s"
elif [[ "$strand" == "antisense" ]]; then strand="-S"
elif [[ "$strand" == "either" ]]; then strand=""
else echo -e "[ERROR] The indicated value '$strand' for the --strand option is illegal! Choose one of 'sense', 'antisense' or 'either'." >&2; usage; exit 1
fi

# Verify that the specified value is a number
if [[ $(echo $overlap | grep -Pq '^[-+]?[0-9]+\.?[0-9]*$'; echo $?) -ne 0 ]]; then echo -e "[ERROR] The indicated value '$overlap' for the --minimum-overlap option is not a valid number!" >&2; usage; exit 1; fi

# Verify that the specified number is in the allowed range
if [[ $(echo "$overlap < 0.000000001" | bc) == 1 || $(echo "$overlap > 1" | bc) == 1 ]]; then echo -e "[ERROR] The indicated value '$overlap' for the --minimum-overlap option is outside the accepted range (0.000000001 to 1)!" >&2; usage; exit 1; fi

# Verify that the specified value is an integer
if [[ ! $threads =~ ^-?[0-9]+$ ]]; then echo -e "[ERROR] The indicated value '$threads' for the --bam-sort-threads option is not an integer!" >&2; usage; exit 1; fi
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#-
#  MAIN //  #
#===========#
#---> SET FLAGS <---#
set -e

#---> MAIN VARIABLES <---#

        #---> BODY <---#

        #---> Print status message <---#
        if [ $verbose -eq 1 ]; then echo "Starting '$0'..." >&2; fi

        #---> Make temporary folder <---#
        tmp=$(mktemp -d $tmp_dir/tmp.XXXXXXXX)
        if [ $verbose -eq 1 ]; then echo "Made temporary directory '$tmp'..." >&2; fi

        #---> Sort BAM file by name <---#
        if [ $sorted -eq 1 ]; then
                bam_sorted=$bam
        else
                if [ $verbose -eq 1 ]; then echo "Sorting input BAM file by read names..." >&2; fi
		cd $tmp
                bam_sorted=$tmp/$(basename $bam .bam).name_sorted.bam
                $samtools sort -n -@ $threads -m $memory $bam $bam_sorted
		bam_sorted=$bam_sorted.bam
		cd -
        fi

        #---> Calculate overlaps <---#
        if [ $verbose -eq 1 ]; then echo "Computing BAM/BED12 overlaps..." >&2; fi
        $bedtools intersect -bed -wo $strand -f $overlap -abam $bam_sorted -b $bed12 > $tmp/overlaps

        #---> Pre-process overlaps <---#
        if [ $verbose -eq 1 ]; then echo "Processing overlaps..." >&2; fi
        cut -f1-6,16 $tmp/overlaps | cut -f1 -d "$" > $tmp/overlaps_processed

        #---> Generate count table <---#
        if [ $verbose -eq 1 ]; then echo "Counting reads per gene..." >&2; fi
        $perl $count_script --only-feature-id $allow_ms $allow_ol --overlaps $tmp/overlaps_processed --count-table $count_table

        #---> Remove temporary folder <---#
        if [ ! $keep_tmp -eq 1 ]; then
                rm -rf $tmp
                if [ $verbose -eq 1 ]; then echo "Removed temporary directory..." >&2; fi
        fi

        #---> Print status message <---#
        if [ $verbose -eq 1 ]; then echo "Done." >&2; fi

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#                                                                                                                   
