#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: bed_score_arithmetics.pl
### Created: Dec 3, 2013
### Modified: Dec 3, 2013
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
my $op = '+';
my $outfile = '';
my $bed1 = '';
my $bed2 = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'operation=s' => \$op,
	'outfile=s' => \$outfile,
	#-----------------------#
	'bed1=s' => \$bed1,
	'bed2=s' => \$bed2
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$bed1 || !$bed2; 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print STDERR "Starting '$0'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my ($bed1_hoa_ref, $bed2_hoa_ref, $res_hoa_ref);

#---> BODY <---#
$bed1_hoa_ref = &bed_to_hoa($bed1);
$bed2_hoa_ref = &bed_to_hoa($bed2);
$res_hoa_ref = &compare_bed_scores($bed1_hoa_ref , $bed2_hoa_ref , $op);
print_hoa_values_by_alphanum_sorted_keys($outfile, $res_hoa_ref);

#---> STATUS MESSAGE <---#
print STDERR "Done.\n" unless $quiet;

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
'Usage: perl ./bed_score_arithmetics.pl [OPTIONS] --bed1 [BED] --bed2 [BED]

Description: Compares the scores of two BED files by the indicated arithmetic operator.

Notes: Compares scores for identical name fields (4th column). Warnings are issued if (1) identical names occur more than once within one file.

==================================================
Required arguments:
--bed1 [BED]	BED file 1
--bed2 [BED]	BED file 2
==================================================
Optional arguments:
--operation [STRING]	Operator used to compare --bed1 and --bed2 scores (default: "+")
--outfile [FILE]	Output filename (default: write to STDOUT)
--usage|help	Show this information and die
--quiet	Shut up!

';
}
#-----------------------#
sub bed_to_hoa {
### Function: Reads a BED file into a hash of arrays, with a specified field being the key (default: 4th/name column); a warning is issues if non-unique hash keys are found
### Accepts: 1. BED file; 2. field/column to use as hash key
### Returns: Reference to hash of arrays (hash keys = specified field; values = arrays containing the whole BED lines)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;
	my $col = (scalar @_) ? shift(@_) : 4;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$in_file'...\n" unless $quiet;
	
	#---> SUBROUTINE VARIABLES <---#
	my %HoA;	
	
	#---> BODY <---#

		#---> Open file handle <---#
		open BED, $in_file;
	
		#---> Traverse through each line of file <---#
		while (<BED>) {

			#---> Remove trailing record separator <---#
			chomp;

			#---> Split line <---#
			my @line = split /\t/, $_;

			#---> Verify correct file format <---#
			die "[ERROR] File '$in_file' does not look like a valid input file in line $.!\nExecution aborted.\n" unless defined $line[$col - 1];

			#---> Extract hash key <---#
			my $key = $line[$col - 1];

			#---> Warn if hash key not unique <---#
			warn "[WARNING] Entry '$key' not unique in line $.!\nDuplicate entry discarded." and next if exists $HoA{$key};

			#---> Add line to hash <---#
			$HoA{$key} = \@line;
		}		
	
		#---> Close file handle <---#
		close BED;

	#---> STATUS MESSAGE <---#
	print STDERR "File '$in_file' processed.\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%HoA;
}
#-----------------------#
sub compare_bed_scores {
### Function: Performs arithmetic operations on the scores of two BED files
### Accepts: 1. Hash of arrays of BED file 1; 2. Hash of arrays of BED file 2; 3. Arithmetic operator
### Returns: TBD
### Dependencies: Subroutine "bed_to_hoa" (Alexander Kanitz, 03-DEC-2013)
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($bed1, $bed2, $op) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Comparing BED scores..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @keys;
	my %tmp;
	my @uniq_keys;
	my $na = 0;
	my %out;
	
	#---> BODY <---#

		#---> Merge keys from both input hashes into a single array <---#
		push @keys, keys %$bed1, keys %$bed2;
		
		#---> Remove duplicates <---#
		@uniq_keys = grep { ! $tmp{$_}++ } @keys;

		#---> Iterate over all unique keys <---#
		foreach my $key (sort @uniq_keys) {

			#---> Initialized variable to hold the new score <---#
			my $new_score;

			#---> Verify presence of entry '$key' in first file <---#
			if ( ! exists $bed1->{$key} ) {

				#---> Set $na switch, new score to "N/A" and warn <---#
				$na = 1;
				$new_score = "N/A";
				warn "[WARNING] Entry '$key' absent in first input file.\nNew score set to 'N/A' and second entry used for writing output.\n";

			}

			#---> Verify presence of entry '$key' in second file <---#
			elsif ( ! exists $bed2->{$key} ) {

				#---> Set new score to "N/A" and warn <---#
				$new_score = "N/A";
				warn "[WARNING] Entry '$key' absent in second input file.\nNew score set to 'N/A' and first entry used for writing output.\n";

			}

			#---> Verify presence of score in first file <---#
			elsif ( ! defined $bed1->{$key}[4] ) {

				#---> Set $na switch, new score to "N/A" and warn <---#
				$na = 1;
				$new_score = "N/A";
				warn "[WARNING] Score missing for entry '$key' in first file.\nNew score set to 'N/A' and second entry used for writing output.\n";

			}

			#---> Verify presence of score in second file <---#
			elsif ( ! defined $bed2->{$key}[4] ) {

				#---> Set new score to "N/A" and warn <---#
				$new_score = "N/A";
				warn "[WARNING] Score missing for entry '$key' in second file.\nNew score set to 'N/A' and first entry used for writing output.\n";

			}

			#---> Verify presence of entry '$key' and a corresponding score in both files <---#
			else {

				#---> Extract current score from first file <---#
				my $score1 = ${$bed1->{$key}}[4];

				#---> Extract current score from second file <---#			
				my $score2 = ${$bed2->{$key}}[4];

				#---> Warn if, except for the score fields, the lines are not identical in both files <---#	
				if ( join ("|", @{$bed1->{$key}}[0..3] ) ne join ( "|", @{$bed2->{$key}}[0..3] ) ) {
					warn "[WARNING] Not all values for entry '$key' are identical in both input files! Values from first file kept.\n";
				}
	
				#---> Calculate new score <---#
				$new_score = $score1 - $score2;

			}
			
			#---> Insert new score into the entry of the first/second file, depending on $na switch <---#	
			if ($na == 1) {
				splice @{$bed2->{$key}}, 4, 1, $new_score;
				$out{$key} = \@{$bed2->{$key}};
				$na = 0;
			}
			else {
				splice @{$bed1->{$key}}, 4, 1, $new_score;
				$out{$key} = \@{$bed1->{$key}};
			}

		}
	
	#---> STATUS MESSAGE <---#
	print STDERR "BED scores compared." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%out;
}
#-----------------------#
sub print_hoa_values_by_alphanum_sorted_keys {
### Function: Writes hash values in the order of alphanumerically sorted hash keys to output file (one pair per line)
### Accepts: 1. Output file; 2. Hash reference; 3. Separator [DEFAULT: TAB]
### Returns: n/a
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS<---#
	my $out_file = shift;
	my $hash_ref = shift;
	my $sep = defined($_[0]) ? $_ : "\t";

	#---> STATUS MESSAGE <---#	
	print STDERR "Writing output to file '$out_file'...\n" if $out_file && ! $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $fh;

	#---> BODY <---#

		#---> Open output filehandle <---#
		if ($outfile) {
			open ($fh, '>', $outfile) or die "[ERROR] Could not open output file '$out_file' to write!\nExecution aborted.\n";

		}
		else {
			open ($fh, '>&', \*STDOUT) or die "[ERROR] Could not write to STDOUT!\nExecution aborted.\n";
		}

		#---> Traverse through all keys of hash referenced by '$hash_ref' <---#
		foreach my $key (sort { $a cmp $b } keys %$hash_ref) {
			
			# Print hash key TAB hash value
			print $fh join ($sep, @{$$hash_ref{$key}}) . "\n";

		}

		# Close OUT filehandle
		close $fh;
	
	#---> STATUS MESSAGE <---#
	print STDERR "File '$out_file' written.\n\n" if $out_file && ! $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#