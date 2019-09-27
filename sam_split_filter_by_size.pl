#!/usr/bin/perl


#================#
#   HEAD START   #
#================#

#---> REQUIREMENTS <---#
require Getopt::Long;
require File::Basename;

#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#---> SET VERSION <---#
my $VERSION = "1.0 (Mar 20, 2014)";

#================#
#    HEAD END    #
#================#


#===================#
#   OPTIONS START   #
#===================#

#---> HELP MESSAGE <---#
sub help {
'USAGE:
   perl ' . $0 . ' [OPTIONS] [--min INT] [--max INT] [--prefix STRING]         
                   [--header file] [--report-only] [--merge STRING] -- FILE(S)

DESCRIPTION:
   Splits or filters SAM file(s) according to read length. Either one file is   
   produced for all read lengths (read length filtering) or one file for each   
   read length (splitting). A report is printed to STDOUT.

OPTIONS:
   -h --help --usage  Show this information and die!

   --version          Show version.

   -q --quiet         Shut up.		

   --prefix STRING    Prefix for output files. The suffix "_XXX" will be        
                      appended to the prefix, where XXX is the read length      
                      [default: basename of first processed input file if       
                      reading from file(s) or "_out_file.sam" if reading from    
                      STDIN].

   --header file      Use content of FILE as header for output files. If no     
                      argument is specified, use header of first processed input
                      file. See note.

   --min INT          Minimum accepted read length [default: no minimum].

   --max INT          Maximum accepted read length [default: no maximum].

   --report-only      Only print read length counts.

   --merge STRING     Only print one output file (useful if used as read length 
                      filter). STRING is appended to --prefix. If no argument is
                      specified, "_filtered" is used as a suffix.


NOTES:
   - If input shall be read from STDIN, use the following syntax:
      perl ' . $0 . ' [OPTIONS] [--prefix STRING] [--header file] -- -
   - If the header is to be grabbed from the first processed file, it is        
      defined as the continuous stretch of lines at the very top of the first   
      processed input file that matches the regular expression "\@\w{2}\t", i.e.
      the "@" character followed by two word characters and a tab.
   - CAUTION: Input files are only marginally validated!
   - CAUTION: Read lengths of greater than 999 are not supported! 

AUTHOR:
   Alexander Kanitz, Biozentrum, University of Basel

CREATED:
   Mar 19, 2014

VERSION:
   ' . $VERSION . '
';
}

#---> VERSION MESSAGE <---#
sub version {
	$0 . ", " . $VERSION . "\n" . "Alexander Kanitz, Biozentrum, University of Basel" . "\n";
}

#---> OPTIONS VARIABLES <---#
my $show_help = '';
my $show_version = '';
my $quiet = '';
my $prefix = '';
my $header_file = 0;
my $min = 0;
my $max = 0;
my $report_only = '';
my $merge = 0;

#---> PARSE OPTIONS <---#
	my $options_result = GetOptions (
		'h|help|usage' => \$show_help,
		'version' => \$show_version,
		'quiet' => \$quiet,
		'p|prefix=s' => \$prefix,
		'h|header:s' => \$header_file,
		'min=i' => \$min,
		'max=i' => \$max,
		'report-only' => \$report_only,
		'merge:s' => \$merge,
	);

#---> VALIDATE OPTIONS <---#

	# Die if option parsing failed
	die &help() . "\n[ERROR] Option parsing failed!\nExecution aborted.\n" unless $options_result;
	
	# Show help and die if requested
	die &help() if $show_help;
	
	# Show version and die if requested
	die &version() if $show_version;
	
	# Die if no input file is given
	die &help() . "\n[ERROR] No input detected!\nExecution aborted..\n" unless @ARGV;
	
	# Die if specified header file cannot be found
	die &help() . "\n[ERROR] File '$header_file' not found!\nExecution aborted.\n" if $header_file && ! -e $header_file;

	# Die if min and/or max values are below 0 or if max is below min
	die &help() . "\n[ERROR] Argument to --min below 0!\nExecution aborted.\n" if $min < 0;
	die &help() . "\n[ERROR] Argument to --max below 0!\nExecution aborted.\n" if $max < 0;
	die &help() . "\n[ERROR] Argument to --max lesser than argument to --min !\nExecution aborted.\n" if $max < $min && $max != 0;

#---> SET SWITCHES <---#

	# Does the header need to be grabbed from the input?
	my $grab_header = $header_file eq "" ? 1 : 0;

#===================#
#    OPTIONS END    #
#===================#


#================#
#   MAIN START   #
#================#

#---> CALL MAIN SUBROUTINE <---#
&main();

#---> EXIT ORDERLY <---#
exit 0;

#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub main {

	#---> STATUS MESSAGE <---#
	print STDERR "Starting '$0'...\n" unless $quiet;

	#---> MAIN VARIABLES <---#
	my @fh;
	my $header = "";
	my $count = 0;
	my @read_lengths;

	#---> BODY <---#	

		#---> Get header if requested <---#
		if ($header_file) {
			open my $fh_header, "<", $header_file || die "Cannot open file '$header_file': $!";
				$header .= <$fh_header>;
			close $fh_header;
		}

		#---> Set default prefix if necessary <---#
		unless ($prefix) {
			# Set to basename (without extension) if reading from file
			if (defined $ARGV[0] && $ARGV[0] ne "-") {
				$prefix = basename($ARGV[0]);
			}
			# Set to "file" if reading from STDIN
			else {
				$prefix = "out_file.sam";
			}
		}
			
		#---> Process input <---#
		while (my $line = <>) {

			#---> Check for header lines <---#
			if ($line =~ m/\@\w{2}\t/) {
				# Append it to previous content of header if header needs to be grabbed
				$header .= $line if $grab_header;
				# Skip to next line
				next;
			}
			else {
				# Unset the "$grab_header" switch as soon as the first non-header line is encountered
				$grab_header = 0;
			}

			#---> Split line & extract read length <---#
			my @fields = split "\t", $line;
			my $read_length = length $fields[9];
			
			#---> Skip read if not in requested size range <---#
			next if $min && $read_length < $min;
			next if $max && $read_length > $max;
						
			#---> Initiate/increment read counter <---#
			if ( defined $read_lengths[$read_length] ) {
				$read_lengths[$read_length]++;
			}
			else {
				$read_lengths[$read_length] = 1;
			}

			#---> One output file for all filtered read lengths <---#
			if ($merge != 0 && ! $report_only) {
				# Build filename
				my $file = $merge ? $prefix . $merge : $prefix . "_filtered";
				# Open file handle
				open my $merge_fh, ">", $file || die "Cannot open file '$file': $!";
				# Print header
				print { $merge_fh } $header;
				# Print line				
				print { $merge_fh } $line;
				# Skip to next read
				next;
			}

			#---> One output file per read length <---#
			unless ($report_only) {
				# Format read length
				my $read_length_formatted = sprintf("%03d", $read_length);
				# Build filename
				my $file = $prefix . "_" . $read_length_formatted;
				# If file handle not yet open...
				unless (defined $fh[$read_length]) {
					# Open file handle
					open $fh[$read_length], ">", $file || die "Cannot open file '$file': $!";
					# Print header
					print { $fh[$read_length] } $header;
				}
				# Print line
				print { $fh[$read_length] } $line unless $report_only;
			}

			#---> Increase line count <---#
			$count++;
			
			#---> Print status message <---#
			print STDERR "[" . $count . " lines processed...]\n" if ! $quiet && ! ($count % 500000);
	
		}

		#---> Print report <---#
		for (my $read_length = 0; $read_length <= $#read_lengths; $read_length++) {
			print STDOUT $read_length . "\t" . $read_lengths[$read_length] . "\n" if defined $read_lengths[$read_length];
		}

		#---> Close output file handles <---#
		unless ($report_only) {
			#close $merge_fh if defined $merge_fh;
			for my $handle (@fh) {
				close $handle if defined $handle;
			}
		}
		
	#---> STATUS MESSAGE <---#
	print STDERR "Done.\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#
