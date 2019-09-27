#!/bin/sh

#=============#
#  HEADER //  #
#=============#
## Created: Jun 25, 2015
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: GNU core utilities, Curl
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
Usage: $script [OPTIONS] --archive <URL> --download-dir <PATH> --install-dir <PATH>

Description: Installs a gzipped tarball in a desired location with 'confiugre', 'make' and 'make install'

Options:
        --archive=URL                Location (URL) of archive to install (required).
        --download-dir=PATH          Path to the root directory for the installation (default: current directory ".").
        --install-dir=PATH           Name of the installation folder (default: built from download directory and extracted folder).
        --verbose                    Write status/progress messages to STDERR.
        --usage | --help             Show this screen and exit.
        --version                    Show version information and exit.

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Jun 25, 2015.
Version 1.0.1 (June 25, 2015)
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
script=$(basename $0)
cat <<VERSION
$script version 1.0 (Jun 25, 2015)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Jun 25, 2015.
VERSION
}
#============#
#  // USAGE  #
#============#


#================#
#  // FUNCTIONS  #
#================#
topmostFile ()
### Find topmost copy of a particular file
{
    dir="$1"
    name="$2"
    find "$dir" -type f -name "$name" | awk -F'/' 'NR == 1 {line = $0; min = NF}; NR > 1 { if (NF == min) {line = ""} else if (NF < min) {line = $0; min = NF} }; END {print line}'
}

#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES / SET DEFAULTS <---#
url=""
downloadDir="$PWD"
installDir=""
verbose=0

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
        --archive)
            url=$2
            shift 2
            ;;
        --download-dir)
            downloadDir=$2
            shift 2
            ;;
        --install-dir)
            installDir=$2
            shift 2
            ;;
        --) # End of all options
            shift
            break
            ;;
        -*)
            echo -e "[ERROR] Unknown option: $${1}!\nExecution aborted.\n" >&2
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
if [ -z "$url" ]; then echo -e "[ERROR] Specify the '--archive' option with a valid argument!\nExecution aborted.\n" >&2 usage; exit 1; fi

# Verify that the specified folder is writable
if [ ! -w $downloadDir ] ; then echo -e "[ERROR] The argument '$downloadDir' to option '--download-dir' is not writable!\nExecution aborted.\n" >&2; usage; exit 1; fi
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#-
#  MAIN //  #
#===========#

#---> MAIN VARIABLES <---#

    #---> BODY <---#

    #---> Print status message <---#
    if [ $verbose -eq 1 ]; then echo "Starting installation..." >&2; fi

    #---> Download installation archive <---#
    if [ $verbose -eq 1 ]; then echo "Downloading archive..." >&2; verbosity=""; else verbosity="--silent"; fi
    mkdir --parents "$downloadDir"
    curl $verbosity --output "${downloadDir}/archive.tar.gz" "$url"
    if [ $? -ne 0 ] ; then echo -e "[ERROR] The file '$url' to option '--archive' could not be downloaded!\n" >&2; usage; exit 1; fi

    #---> Get/make installation directory <---#
    if [ $verbose -eq 1 ]; then echo "Creating installation directory..." >&2; fi
    topLevelContents=$(tar --gzip --list --file "${downloadDir}/archive.tar.gz" | sed --regexp-extended 's/\.+\///' | cut --delimiter "/" --fields 1 | uniq)
    topLevelDirs=$(echo "$topLevelContents" | wc -l)
    if [ -z "$installDir" ]; then
        if [ $topLevelDirs -gt 1 ]; then
            echo -e "[ERROR] There is more than a single top-level directory in the archive and no installation directory has been specified via '--install-dir'!\nExecution aborted.\n" >&2
            usage
            exit 1
        else
            installDir="$downloadDir"
        fi
    else
        mkdir --parents "$installDir"
    fi

    #---> Extract installation archive <---#
    if [ $verbose -eq 1 ]; then echo "Extracting archive..." >&2; verbosity="--verbose"; else verbosity=""; fi
    tar --extract --gzip $verbosity --directory "$installDir" --file "${downloadDir}/archive.tar.gz"
    if [ $topLevelDirs -eq 1  ]; then
        installDir="${installDir}/${topLevelContents}"
    fi
    mv "${downloadDir}/archive.tar.gz" "${installDir}"

    #---> Run 'configure' <---#
    if [ $verbose -eq 1 ]; then echo "Configuring..." >&2; verbosity=""; else verbosity="> /dev/null 2> /dev/null"; fi
    mkdir --parents "${installDir}/BUILD"
    cd "${installDir}/BUILD"
    config=$(topmostFile $installDir 'configure')
    $config --prefix "$installDir" $verbosity
    if [ $? -ne 0 ] ; then echo -e "[ERROR] Configuration error! Run with '--verbose' for 
    details.\nExecution aborted\n" >&2; usage; exit 1; fi

    #---> Run 'make' <---#
    if [ $verbose -eq 1 ]; then echo "Installing..." >&2; verbosity=""; else verbosity="--silent"; fi
    make $verbosity
    if [ $? -ne 0 ] ; then echo -e "[ERROR] Installation error! Run with '--verbose' for details.\nExecution aborted\n" >&2; usage; exit 1; fi
    make install $verbosity
    if [ $? -ne 0 ] ; then echo -e "[ERROR] Installation error! Run with '--verbose' for details.\nExecution aborted\n" >&2; usage; exit 1; fi
    cd -
    rm -rf "${installDir}/BUILD"

    #---> Print status message <---#
    if [ $verbose -eq 1 ]; then echo "Done." >&2; fi

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#                                                                                                                   
