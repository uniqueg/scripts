#!/usr/bin/env perl


#=====================#
#  ISSUES & IDEAS //  #
#=====================#
# TODO Include option '--outfile' to redirect output to file
# TODO Include option '--split' such that output for each identifier will be written to a separate
#      file (prefixes are left as they are, but before extension, identifier is inserted); could be
#      achieved by reading output of getResponse with Bio::SeqIO; record counting could also be
#      integrated there
# TODO Add warnings for identifers that were not found; not sure how to do this as getResponse
#      returns only an output string and no state information or anything...
# TODO The script may be problematic when requesting a lot of data, as it is kept in memory; also
#      there are options 'retstart' and 'retmax' which may restrict the max number of inputs 
#      (10000?)
#=====================#
#  // ISSUES & IDEAS  #
#=====================#


#==============#
#  PRAGMAS //  #
#==============#
use strict;
use warnings;
#==============#
#  // PRAGMAS  #
#==============#


#==============#
#  MODULES //  #
#==============#
use File::Basename;
use Getopt::Long;
use Bio::DB::EUtilities;
use Path::Class;
use POSIX;
#==============#
#  // MODULES  #
#==============#


#===============================#
#  VERSION, LICENSE & USAGE //  #
#===============================#
sub version {
### Returns version information as a string
my $scriptName = basename($0);
<<VERSION;
[VERSION INFORMATION]
    $scriptName, v1.0 (Aug 11, 2015).
    Original version from Aug 11, 2015.

[CONTACT INFORMATION]
    Alexander Kanitz <alexander.kanitz\@alumni.ethz.ch>
    Biozentrum, University of Basel
VERSION
}
#-----------------------#
sub license {
### Returns license information as a string
<<LICENSE;
Copyright (c) 2015 Alexander Kanitz

This code is released under the MIT/X11 License:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
LICENSE
}
#-----------------------#
sub usage {
### Returns usage information as a string
my $scriptName = basename($0);
my $usage = <<USAGE;
[USAGE]
    $scriptName [OPTIONS] (<IDENTIFIER> ...) > <OUTPUT_FILE>

[DESCRIPTION]
    Retrieves NCBI/GenBank nucleotide sequence records given one or more sequence identifiers (GI,
    RefSeq accession or RefSeq version).

[REQUIRED MODULES]
    File::Basename, Getopt::Long, Bio::DB:EUtilities, Bio::SeqIO,
    Path::Class, POSIX

[OPTIONS]
    --identifiers <FILE>    Read identifiers from file (one per line). If identifiers are specified
                            on the command-line as well, these are processed first. Note that it is
                            also possible to read identifiers from <STDIN> (see comments).
    --database <STRING>     NCBI database to query. One of 'nucgss', 'nucest', or
                            'nuccore' (default).
    --format <STRING>       Desired format of output. One of 'acc' (accession), 'gi' (GenBank
                            identifer), 'gss', 'est', 'fasta', 'gbwithparts' (GenBank format; only
                            required for long sequences/contigs), or 'gb' (compact GenBank format;
                            default).
    --email <EMAIL>         Email address (required by NCBI; default: 'dummy\@foo.bar'). Note that
                            the email syntax is not validated by this script.
    --verbose               Print detailed log information to STDERR.
    --version               Show version information and exit.
    --license               Show license information and exit.
    --usage | --help        Show this screen and exit.

[COMMENTS]
    (1) The script looks for identifiers on the command-line first, then in the file passed to
        option '--identifiers' (if specified). If no identifiers were found, the script will then
        (and only then!) listen to <STDIN>. Therefore, if you want to run the script as part of a
        pipe, ensure that option '--identifiers' is not set and that no plain-text identifiers are
        specified on the command-line.
USAGE
&version . "\n" . $usage;
}
#===============================#
#  // VERSION, LICENSE & USAGE  #
#===============================#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES <---#
my $idsFile = '';
my $database = 'nuccore';
my $format = 'gb';
my $mail = 'dummy@foo.bar';
my $verbose = 0;
my $version = 0;
my $license = 0;
my $usage = 0;

#---> PARSE / ASSIGN OPTIONS <---#
my $options_result = GetOptions (
    'identifiers=s' => \$idsFile,
    'database=s' => \$database,
    'format=s' => \$format,
    'email=s' => \$mail,
    'verbose' => \$verbose,
    'version' => \$version,
    'license' => \$license,
    'usage|help' => \$usage
);

#---> VERIFY OPTIONS <---#
# Print usage information and exit if option parsing was not successful
die &usage unless $options_result;

# Print usage information and exit if '--usage' or '--help' were specified
die &usage if $usage;

# Print version information and exit if '--version' was specified
die &version if $version;

# Print license information and exit if '--license' was specified
die &license if $license;

# Print usage information and die if prohibited option arguments specified
die "[ERROR] Argument to option '--database' has to be either 'nucgss', 'nucest', or 'nuccore' (default).\n",
    "[ERROR] Execution aborted!\n\n",
    &usage if $database ne 'nucgss' &&
              $database ne 'nucest' &&
              $database ne 'nuccore';
die "[ERROR] Argument to option '--format' has to be one of 'acc', 'gi', 'gss', 'est', 'fasta', 'gbwithparts', or 'gb' (default).\n",
    "[ERROR] Execution aborted!\n\n",
    &usage if $format ne 'acc' &&
              $format ne 'gi' &&
              $format ne 'fasta' &&
              $format ne 'gss' &&
              $format ne 'est' &&
              $format ne 'gbwithparts' &&
              $format ne 'gb';
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#
#  MAIN //  #
#===========#
#---> STATUS MESSAGE <---#
print STDERR "[" . &getLogDateTime() . "] Starting '$0'...\n" if $verbose;

#---> MAIN VARIABLES <---#
my $idsRef;
my $genBankRecords;

#---> BODY // <---#

    #---> Extract GenBank/NCBI identifiers <---#
    $idsRef = &getIDs($idsFile);

    #---> Handle errors <---#
    die "[ERROR] Could not open file '$idsFile'.\n",
        "[ERROR] Execution aborted!\n\n",
        &usage if $idsRef == 1;
    die "[ERROR] No identifiers were found.\n",
        "[ERROR] Execution aborted!\n\n",
        &usage if $idsRef == 2;

    #---> Retrieve GenBank entries <---#
    $genBankRecords = &retrieveGenBankNucleotideEntries($idsRef, $database, $format, $mail);

    #---> Print GenBank records <---#
    print STDOUT $genBankRecords;

#---> // BODY <---#

#---> STATUS MESSAGE <---#
print STDERR "[" . &getLogDateTime() . "] Done.\n" if $verbose;

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#


#==================#
#  SUBROUTINES //  #
#==================#
sub getIDs {
## Description: Get identifiers from command-line, dedicated file and/or <STDIN>
## Accepts: n/a (uses globals!)
## Returns: Reference to array of identifiers
## Dependencies: n/a
## Type: generic

    #---> STATUS MESSAGE <---#
    print STDERR "[" . &getLogDateTime() . "] Collecting identifiers..." . "\n" if $verbose;

    #---> PASS ARGUMENTS ---#
    my ($idsFile) = @_;

    #---> SUBROUTINE VARIABLES <---#
    my @ids;
    my $fileHandle;

    #---> BODY // <---#

        #---> Get identifiers from command-line <---#
        @ids = @ARGV;

        #---> Get identifiers from input file, if available <---#
        if ($idsFile) {

            #---> Open input filehandle <---#
            open $fileHandle, "<", $idsFile or return 1;

            #---> Traverse input file line by line <---#
            while (<$fileHandle>) {

                #---> Remove trailing newline character and add to identifiers array <---#
                chomp;
                push @ids, $_;

            }

        }

        #---> If not identifiers have been found yet, try to get from <STDIN> <---#
        if (! scalar @ids) {

            #---> Traverse <STDIN> line by line <---#
            while (<>) {

                #---> Remove trailing newline character and add to identifiers array <---#
                chomp;
                push @ids, $_;

            }

        }

        #---> Return error code if still no identifiers have been found <---#
        return 2 unless scalar @ids;

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    if (scalar @ids == 1) {
        print STDERR "[" . &getLogDateTime() . "] " . scalar @ids . " identifier found." . "\n" if $verbose;
    } else {
        print STDERR "[" . &getLogDateTime() . "] " . scalar @ids . " identifiers found." . "\n" if $verbose;
    }

    #---> RETURN VALUE <---#
    return \@ids;

}
#-----------------------#
sub retrieveGenBankNucleotideEntries {
## Description: Retrieves GenBank entries given a list of GenBank/NCBI nucleotide sequence identifiers.
## Accepts: 1. reference to array of identifiers
##          2. query database (one of 'nuccore', 'nucest', or 'nucgss')
##          3. output format (one of 'gb', 'gbwithparts', 'gss', 'est', 'fasta', 'gi', or 'acc')
##          4. email address (required by NCBI/GenBank)
## Returns: Results
## Dependencies: n/a
## Type: generic

    #---> STATUS MESSAGE <---#
    print STDERR "[" . &getLogDateTime() . "] Retrieving GenBank entries..." . "\n" if $verbose;

    #---> PASS ARGUMENTS ---#
    my ($idsRef, $database, $format, $mail) = @_;

    #---> SUBROUTINE VARIABLES <---#
    my $efetchConn;
    my $results;
    my $nRecords;

    #---> BODY // <---#

        #---> Construct efetch connector for format 'gi' <---#
        $efetchConn = Bio::DB::EUtilities->new(
            -eutil   => 'efetch',
            -id      => $idsRef,
            -db      => $database,
            -rettype => 'gi',
            -retmode => 'text',
            -email   => $mail
        );

        #---> Connect and save HTTP response for format 'gi'
        eval {
            $results = $efetchConn->get_Response->content;
            1;
        } or do {
            die "[ERROR] Server error:$@" .
                "[ERROR] Possibly no valid identifiers were passed, otherwise try again later.\n" .
                "[ERROR] Execution aborted!\n";
        };

        #---> Connect and save HTTP response to variable
        $nRecords = scalar split /\n/, $results;

        #---> If specified format is not 'gi', reconnect and save HTTP response for specified format
        $efetchConn->set_parameters(-rettype => $format);
        eval {
            $results = $efetchConn->get_Response->content;
            1;
        } or do {
            die "[ERROR] Server error:$@" .
                "[ERROR] Possibly no valid identifiers were passed, otherwise try again later.\n" .
                "[ERROR] Execution aborted!\n";
        };

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    if ($nRecords == 1) {
        print STDERR "[" . &getLogDateTime() . "] $nRecords GenBank record of format '$format' retrieved." . "\n" if $verbose;
    } else {
        print STDERR "[" . &getLogDateTime() . "] $nRecords GenBank records of format '$format' retrieved." . "\n" if $verbose;
    }

    #---> RETURN VALUE <---#
    return $results;

}
#-----------------------#
sub getLogDateTime {
## Description: Returns current date/time stamp in a format suitable for log entries.
## Accepts: 1. template string (optional; default: "%Y/%m/%d %H:%M:%S")
## Returns: Date/time string
## Dependencies: POSIX
## Type: genericHelper

    #---> REQUIRED MODULES <---#
    use POSIX;

    #---> PASSED ARGUMENTS <---#
    my $template = scalar @_ ? shift : "%Y/%m/%d %H:%M:%S";

    #---> SUBROUTINE VARIABLES <---#
    my $dateTime;

    #---> BODY // <---#

        #---> Build output filename <---#
        $dateTime = strftime( $template, localtime() );

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return $dateTime;

}
#==================#
#  // SUBROUTINES  #
#==================#
