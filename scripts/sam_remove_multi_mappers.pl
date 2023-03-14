#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: sam_remove_multi_mappers.pl
### Created: Jun 19, 2013
### Modified: Jun 19, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: GetOpt::Long
#==================#
### Description: Removes from a SAM file that is sorted by read ID all lines with read IDs that occur more than once (irrespective of differences in edit distance!)
### Output: 1. SAM file with duplicates removed (original header, if any, is kept); 2. File with duplicates removed (either SAM or read IDs only in TXT file)
### Usage: perl ./sam_remove_multi_mappers.pl for information on required and optional arguments
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
my $in = '';
my $uniq = '';
my $dupl = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'd|dupl=s' => \$dupl,
	#-----------------------#
	'i|in=s' => \$in,
	'u|uniq=s' => \$uniq
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$in || !$uniq; 

#---> GLOBAL VARIABLES <---# 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'sam_remove_multi_mappers.pl'...\n\n" unless $quiet;

#---> BODY <---#
&samRemoveMultiMappers($in, $uniq, $dupl, 1, "\t", ["@..\t"]);

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
'Usage: perl ./sam_remove_multi_mappers.pl [OPTIONS] --in [FILE|SAM] --uniq [FILE|SAM]

Description: Removes from a SAM file that is sorted by read ID all lines with read IDs that occur more than once (irrespective of differences in edit distance!)

==================================================
Required arguments:
--in	Input file sorted by read ID [SAM]
--uniq	Output file containing only unique read IDs [SAM]
==================================================
Optional arguments:
--usage|help	Show this information and die
--quiet	Shut up!
--dupl	Output file containing duplicate read IDs, one per line [TXT]
';
}
#-----------------------#
sub samRemoveMultiMappers {
### Function: Removes from a tabular file separated by a specified delimiter all lines that have duplicate values for a specified field/column; the file has to be sorted by the values of that field/column
### Accepts: 1. Input file [FILE|TAB]; 2. Filename for output file with unique values in specified field/column; same format as input file [FILE|TAB]; 3. Filename for output file with duplicate values in the specified field/column; one per line; each value only once [FILE|TXT]; 4. Field/column to check for duplicates; 1-based [INT]; 5. Delimiter separating fields/columns of input file [STRING]; 6. Set of characters/strings excluding duplicate checking when at the beginning of a line (e.g. comment characters) [ARRAY REF]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($in, $uniq, $dupl, $target, $sep, $skip) = @_;
	
	#print $$skip;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Processing file '$in'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	# Initialize empty array reserved for keeping track of values of the last line (value in previous $target field and previous line)
	my @prev = ();
	# Initialize $duplicate flag to FALSE
	my $duplicate = 0;
	
	#---> BODY <---#
	
	## OPEN FILE HANDLES
	# Open file handle IN
	open IN, "$in";
	# Open file handle UNIQ
	open UNIQ, ">$uniq";
	# Open file handle DUPL
	open DUPL, ">$dupl" if $dupl;
		
	## LOOP THROUGH INPUT FILE
	# Label loop
	LINE:
	## Traverse through file handle IN line by line
	while (my $line = <IN>) {
		
		## TEST IF LINE SHALL BE SKIPPED
		## Traverse through each element $skip_str of array @$skip
		foreach my $skip_str ( @$skip ) {
			## Check whether $line starts with $skip_str
			if ( $line =~ /\A$skip_str/ ) {
				# If so, print line
				print UNIQ $line;
				# ...and proceed to next line
				next LINE;
			}
		}
		
		# SPLIT LINE BY DELIMITER $sep
		my @fields = split /$sep/, $line, $target + 1;
		
		## TEST FOR DUPLICATES IN TARGET FIELD
		## Check whether array @prev is non-empty (should evaluate to FALSE only during first occurrence of a non-skipped line!)
		if ( @prev ) {
			## If so, check whether value in $target field matches first element of @prev (i.e. is a duplicate of the value in the previous $target field!)
			if ( $fields[$target - 1] eq $prev[0] ) {
				# If so, set $duplicate flag to TRUE
				$duplicate = 1;
			}
			## If not... (i.e. value in $target field is not a duplicate of value in previous $target field)
			else {
				## Check whether the values in the previous $target fields were not duplicates (second element of @prev should then be non-zero)
				if ( !$duplicate ) {
					# If so, write value of second element of @prev (i.e. the original previous line) to file handle UNIQ
					print UNIQ $prev[1];
				}
				## If not...
				else {
					# Write value of first element of @prev (i.e. the duplicate $target field value) to file handle DUPL
					print DUPL $prev[0] . "\t" . $prev[1] . "\n" if $dupl;
					# Set $duplicate flat to FALSE
					$duplicate = 0;
				}
			}
		}
		
		# ADD VALUES OF CURRENT LINE TO ARRAY @prev
		@prev = ($fields[$target - 1], $line)
		
	}
		
	## CLOSE FILE HANDLES
	# Close file handle IN
	close IN;
	# Close file handle DUPL
	close DUPL if $dupl;
	# Close file handle UNIQ
	close UNIQ;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "File '$in' processed. Output written to '$uniq' (input file free of duplicates) and '$dupl' (duplicate values)." . "\n\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return 0;
}
#=======================#
#    SUBROUTINES END    #
#=======================#