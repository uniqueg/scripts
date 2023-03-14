#!/usr/bin/perl


#================#
#   HEAD START   #
#================#

#---> REQUIREMENTS <---#
require Getopt::Long;

#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;

#---> SET VERSION <---#
my $VERSION = "1.0 (Mar 26, 2014)";

#================#
#    HEAD END    #
#================#


#===================#
#   OPTIONS START   #
#===================#

#---> HELP MESSAGE <---#
sub help {
'USAGE:
   perl ' . $0 . ' [OPTIONS] [--delim STRING] [--field INT] [--bin INT]         
                   [--prefix STRIING] --input FILE

DESCRIPTION:
   Splits a delimited, human-readable file by the n-th different occurrence of a
   specified field.

OPTIONS:
   -h --help --usage   Show this information and die!

   -v --version        Show version.

   -q --quiet          Shut up.		

   -i --input FILE     Input file.

   -d --delim STRING   Delimiter separating fields in the input file [default:  
                       "\t"].

   -f --field INT      Reference field/column, 1-based [default: 1].

   -b --bin INT        Number of different reference field entries to bin in one
                       file [default: 100].

   -p --prefix STRING  Prefix for output files. Will be appended by "_X" where  
                       XXX is a serial number increasing by 1 with each split.  
                       [default: input filename]

NOTES:
   - PUT_NOTES_HERE!!!

AUTHOR:
   Alexander Kanitz, Biozentrum, University of Basel

CREATED:
   Mar 26, 2014

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
my $input = '';
my $delim = "\t";
my $field = 1;
my $bin = 100;
my $prefix = '';

#---> PARSE OPTIONS <---#
	my $options_result = GetOptions (
		'h|help|usage' => \$show_help,
		'v|version' => \$show_version,
		'q|quiet' => \$quiet,
		'i|input=s' => \$input,
		'd|delim=s' => \$delim,
		'f|field=i' => \$field,
		'b|bin=i' => \$bin,
		'p|prefix=s' => \$prefix
	);

#---> VALIDATE OPTIONS <---#

	# Die if option parsing failed
	die &help() . "\n[ERROR] Option parsing failed!\nExecution aborted.\n" unless $options_result;
	
	# Show help and die if requested
	die &help() if $show_help;
	
	# Show version and die if requested
	die &version() if $show_version;
	
	# Die if required options are missing
	die &help() . "\n[ERROR] No input file specified!\nExecution aborted.\n" unless $input;
	
#---> SET DEFAULTS <---#

	# Set default prefix if unset
	$prefix = $input unless $prefix;

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
	my $counter = 0;
	my $previous_ref_field_entry = "";
	my $output_fh;
	my $file_counter = 1;

	#---> BODY <---#	

		#---> Open input file handle <---#
		open my $input_fh, "<", $input or die "[ERROR] Could not open file '$input'!\n";
	
		#---> Process input <---#
		while (my $line = <$input_fh>) {

			#---> Split line <---#
			my @fields = split /$delim/, $line, $field + 1;

			#---> Handle first line: open output file handle <---#
			open $output_fh, ">", "${prefix}_${file_counter}" unless $previous_ref_field_entry;
			
			#---> Increase counter if current reference field entry is different from previous <---#
			$counter++ if $fields[$field - 1] ne $previous_ref_field_entry;
			
			#---> Close current output file handle, reset counter, increase file counter and open new one if counter reaches bin size <---#
			if ($counter > $bin) {
				close $output_fh;
				$counter = 1;
				$file_counter++;
				open $output_fh, ">", "${prefix}_${file_counter}";
			}

			#---> Print line <---#
			print { $output_fh } $line;			

			#---> Set current reference field entry as previous reference field entry <---#
			$previous_ref_field_entry = $fields[$field - 1];

		}

		#---> Close input file handle <---#
		close $input_fh;
			
	#---> STATUS MESSAGE <---#
	print STDERR "Done.\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#