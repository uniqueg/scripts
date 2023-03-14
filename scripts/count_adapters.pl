#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: count_adapters.pl
### Created: Nov 11, 2013
### Modified: Nov 11, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: Getopt::Long
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;


#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $adapters = '';
my $report = '';
my $chunk = 100000;
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'adapters=s' => \$adapters,
	'report=s' => \$report,
	'chunk=i' => \$chunk
);
die $usage_info if $usage || !$options_result;

#---> GLOBAL VARIABLES <---# 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> BODY <---#
&grep_multiple_patterns_in_multiple_files($report, $adapters, $chunk, @ARGV);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./count_adapters.pl [OPTIONS] [FILES]

Description: Count the occurrences of patterns/adapters in one or more files.

==================================================
Options:
--adapter	Indicate file with adapters (one adapter per line, no comments)
--report	Print results to file instead of standard output
--chunk		Size of the chunk to be tested (default: 100,000)
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub grep_multiple_patterns_in_multiple_files {
### Function: Uses the Perl grep function to search for multiple patterns in multiple files and prints a count report (number of occurrences per pattern per file))
### Accepts: 1. Filename of the output file (STDOUT will be used if empty string); 2. Path pointing to a file of patterns; 3. Chunk size; 4+. Filenames of files that should be searched
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $report = shift;
	my $adapters = shift;
	my $chunk = shift;
	my @files = @_;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Searching for patterns..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $fh;
	my $pattern;
	
	#---> BODY <---#

		#---> Open output filehandle <---#
		if ($report) {
			open($fh, '>', $report) or die;
		}
		else {
			open($fh, '>&', \*STDOUT) or die;
		}

		#---> Iterate over patterns <---#			
		open PATTERNS, "<", $adapters or die "[ERROR] Could not open file '$adapters'!\n";
		while ($pattern = <PATTERNS>) {
			chomp $pattern;
			print $fh "PATTERN: $pattern\n";
			
			#---> Iterate over files <---#
			for my $file (@files) {
				open FILE, "<", $file or die "[ERROR] Could not open file '$adapters'!\n";
				my $count = 0;
				while (<FILE>) {
					last if $. > $chunk;
					$count++ if $_ =~ /$pattern/;
				}
				close FILE;
				print $fh "File: " . $file . ": " . $count . "\n";
			}
			
			print $fh "---\n"
		}
	
		close PATTERNS;
		close $fh;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Patterns searched" . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return ;
}
#=======================#
#    SUBROUTINES END    #
#=======================#
