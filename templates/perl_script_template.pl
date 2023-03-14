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
my $VERSION = "1.0 (Jan 01, 2099)";

#================#
#    HEAD END    #
#================#


#===================#
#   OPTIONS START   #
#===================#

#---> HELP MESSAGE <---#
sub help {
'USAGE:
   perl ' . $0 . ' [OPTIONS] [--optional] [--optional-too INT] -- REQUIRED

DESCRIPTION:
   Splits or filters SAM file(s) according to read length. Either one file is   
   produced for all read lengths (read length filtering) or one file for each   
   read length (splitting). A report is printed to STDOUT.

OPTIONS:
   -h --help --usage  Show this information and die!

   --version          Show version.

   -q --quiet         Shut up.		


NOTES:
   - PUT_NOTES_HERE!!!

AUTHOR:
   Alexander Kanitz, Biozentrum, University of Basel

CREATED:
   Jan 01, 2099

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


#---> PARSE OPTIONS <---#
	my $options_result = GetOptions (
		'h|help|usage' => \$show_help,
		'version' => \$show_version,
		'quiet' => \$quiet,

	);

#---> VALIDATE OPTIONS <---#

	# Die if option parsing failed
	die &help() . "\n[ERROR] Option parsing failed!\nExecution aborted.\n" unless $options_result;
	
	# Show help and die if requested
	die &help() if $show_help;
	
	# Show version and die if requested
	die &version() if $show_version;
	

#---> SET SWITCHES <---#


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


	#---> BODY <---#	

		#---> Get header if requested <---#


		#---> Set default prefix if necessary <---#

			
		#---> Process input <---#
		while (my $line = <>) {


	
		}

		#---> Print report <---#


		#---> Close output file handles <---#

		
	#---> STATUS MESSAGE <---#
	print STDERR "Done.\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#