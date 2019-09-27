#!usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: miRBase_hairpin_short_name_to_dot_bracket_structure.pl
### Created: Nov 13, 2013
### Modified: Nov 18, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: GetOpt::Long
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
my $all = '';
my $record_delim = "";
my $field = 1;
my $field_delim = "\t";
my $dot_bracket = '';
my $structures = '';
my $output = '';
my $options_result = GetOptions (
	'h|usage|help' => \$usage,
	'v|quiet|no_verbose' => \$quiet,
	'a|all' => \$all,
	'r|record-delim=s' => \$record_delim,
	'f|field=s' => \$field,
	'd|field-delim=s' => \$field_delim,
	'k|dot-bracket=s' => \$dot_bracket,
	's|structures=s' => \$structures,
	'o|output=s' => \$output
);

die $usage_info if $usage || !$options_result;

## Verify/process options/arguments
die "[ERROR] One of '--structures' and '--dot-bracket' has to be provided!\nExecution aborted.\n\n$usage_info" unless $structures || $dot_bracket;
die "[ERROR] Illegal argument to '--field'!\nExecution aborted.\n\n$usage_info" if $field <= 0;
$field-- unless $field == 0;	# Set --field to 0-based coordinates
warn "[WARNING] File '$structures' ignored, argument to '--dot-bracket' provided.\n" if $structures && $dot_bracket;
warn "[WARNING] Input file(s) ignored, option '--all' selected.\n" if $all && $ARGV[0];

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print STDERR "Starting 'miRBase_hairpin_short_name_to_dot_bracket_structure.pl'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $IDs_array_ref;
my $dot_bracket_hash_ref;

#---> BODY <---#
	
	# Process miRBase hairpin short names	
	$IDs_array_ref = &miRBase_hairpin_short_names_to_array($record_delim, $field, $field_delim) unless $all;

	# Get structures in dot-bracket format from '--dot-bracket' file (miRNA hairpin short name TAB dot-bracket structure)
	# or from '--structures' file ('miRNA.str'-like format, available at miRBase FTP site)
	$dot_bracket_hash_ref = $dot_bracket ? &read_dot_bracket_structures_to_hash($dot_bracket) : &miRBase_hairpin_structures_to_dot_bracket_hash($structures);

	# Print structures
	&print_selected_keys_or_all_hash($dot_bracket_hash_ref, $IDs_array_ref, $all, $output);

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
'Usage: perl miRBase_hairpin_short_name_to_dot_bracket_structure.pl [OPTIONS] ([INPUT FILE/S])

Description: Returns the dot-bracket structures of miRNA hairpins for miRBase miRNA hairpin short names specified via the input stream or one or more input files. Several options allow parsing of different input stream/file formats. It is also possible to select an option to return all structures, thus eliminating the requirment to provide hairpin short names. Structures have to be provided in the FASTA-like miRBase ".str" format or in dot-bracket notation in a .tab format (miRNA hairpin short name TAB dot-bracket structure). Structures are written either to standard out or a specified file in the latter file format.
At the time of writing (18-NOV-2013), the file "miRNA.str" was a suitable argument to the "--structures" option in this script. It was available at "ftp://mirbase.org/pub/mirbase/CURRENT/" and contained the structures of all miRNAs currently annotated on miRBase version 20.

Options:
--all		writes dot-bracket structures for all entries in the --structures file (STDIN and argument to --input are ignored)
--record-delim	[STRING]	indicates the delimiter separating records in the input file/stream (default: system default)
--field		[INT]	select the n-th field of the input file/stream to locate miRBase hairpin short names
--field-delim	[STRING]	indicates the delimiter separating fields in the input file/stream (default: tab)
--dot-bracket	[FILE]	subset set of dot-bracket structures from an existing file of dot-bracket structures (format: miRBase hairpin short name <TAB> dot bracket structure)
--structures	[FILE]	indicates path to "miRNA.str" file, available at miRBase FTP site; ignored if --dot-bracket is provided
--output	[FILE]	output file (default: write to STDOUT)
';
}
#-----------------------#
sub miRBase_hairpin_short_names_to_array {

	#---> PASS ARGUMENTS ---#
	my ($record_delim, $field, $field_delim) = @_;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing miRBase short names..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @IDs;
	my $orig_delim;

	#---> BODY <---#

		#---> Set record delimiter <---#
		if ($record_delim) {
			$orig_delim = $/;
			$/ = $record_delim;
		}

		#---> Push miRBase hairpin short names from input file/stream to array <---#
		while (<>) {
			chomp;
			push @IDs, ( (split /$field_delim/ )[$field] =~ m/\A(\w{3}-.*)\Z/ ) if defined ( (split /$field_delim/ )[$field] );
		}

		#---> Reset record delimiter <---#		
		$/ = $orig_delim if $record_delim;

	#---> STATUS MESSAGE <---#
	print STDERR "Processed miRBase short names." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@IDs;

}
#-----------------------#
sub miRBase_hairpin_structures_to_dot_bracket_hash {

	#---> PASS ARGUMENTS ---#
	my $structures_file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing miRNA hairpin structures..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %structures;
	my $orig_delim;

	#---> BODY <---#

		#---> Open structures file <---#
		open STR, "<", $structures_file or die "[ERROR] File '$structures_file' could not be opened! Execution aborted.\n";

		#---> Set record delimiter <---#
		$orig_delim = $/;
		$/ = ">";

		#---> Push miRBase hairpin structure records to hash <---#
		#---> Keys: miRBase hairpin short names; values: reference to array of record lines w/o ID and empty lines <---#
		while (<STR>) {

			# Remove record separator ">"
			chomp;

			# Remove empty records (first record will be empty!)
			next if /^\s*$/;

			# Split record by newline characters
			my @lines = split /\n/;

			# Remove "empty lines" (i.e. array elements consisting of only whitespace characters)
			@lines = grep /\S/, @lines;

			# Extract record ID (miRBase hairpin short name)
			my $short_name = ( split /\s+/, (shift @lines) )[0];

			## Add record to hash (warn if structure is not in proper/expected format or if miRBase hairpin short name has already been seen)
			unless (exists $structures{$short_name}) {
				$structures{$short_name} = &convert_to_dot_bracket(@lines);
				warn "[WARNING] Improper structure format for miRNA hairping '$short_name'. Dot-bracket structure set to 'N/A'.\n" and $structures{$short_name} = "N/A" if $structures{$short_name} eq "1";
			}
			else {
				warn "[WARNING] Multiple entries for miRNA hairpin '$short_name'. Only first entry kept.\n";
			}

		}

		#---> Reset record delimiter <---#		
		$/ = $orig_delim;

		#---> Close structures file <---#
		close STR;

	#---> STATUS MESSAGE <---#
	print STDERR "Processed miRNA hairpin structures." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%structures;

}
#-----------------------#
sub convert_to_dot_bracket {
	
	#---> PASS ARGUMENTS ---#
	my @structure = @_;

	#---> SUBROUTINE VARIABLES <---#
	my $length = length $structure[0];
	my $upper_arm;
	my $lower_arm;
	my $turnaround = "";
	my $dot_bracket;

	#---> BODY <---#
	
		#---> Assert that all structures have either three or five rows; pad three-row entries with spacer rows at beginning and end <---#
		if (scalar @structure == 3) {
			push @structure, " " x $length;
			unshift @structure, " " x $length;
		}
		elsif (scalar @structure != 5) {
			return 1;
		}

		#---> Assert that all rows of structure have same length <---#
		foreach (@structure) {
			return 1 unless $length == length $_;
		}

		#---> Substitute unpaired bases with dots <---#
		$structure[0] =~ tr/a-zA-Z/./;
		$structure[4] =~ tr/a-zA-Z/./;
		
		#---> Substitute paired bases with brackets <---#
		$structure[1] =~ tr/a-zA-Z/(/;
		$structure[3] =~ tr/a-zA-Z/)/;

		#---> Interleave unpaired and paired bases for each strand <---#
		$upper_arm = &interleave_two_strings( ( substr $structure[0], 0, -1 ) , ( substr $structure[1], 0, -1 ) );
		$lower_arm = &interleave_two_strings( ( substr $structure[4], 0, -1 ) , ( substr $structure[3], 0, -1 ) );

		#---> Concatenate bases in the last column ("turnaround") <---#
		$turnaround .= substr $_, -1 foreach @structure;
		return 1 if $turnaround =~ m/\|/;
		$turnaround =~ tr/a-zA-Z()/./;

		#---> Join upper arm, turnaround and reverse lower arm <---#
		$dot_bracket = $upper_arm . $turnaround . reverse $lower_arm;

		#---> Remove placeholders and spaces <---#		
		$dot_bracket =~ tr/- //d;
	
	#---> RETURN VALUE <---#
	return $dot_bracket;

}
#-----------------------#
sub interleave_two_strings {

	#---> PASS ARGUMENTS ---#
	my @first = split //, shift;
	my @second = split //, shift;

	#---> SUBROUTINE VARIABLES <---#
	my $interleaved = "";

	#---> BODY <---#

		#---> Return non-zero exit status if strings are not of same length <---#
		return 1 unless scalar @first == scalar @second;

		#---> Interleave <---#
		for (my $i = 0; $i < scalar @first; $i++) {
			$interleaved .= $first[$i] . $second[$i];
		}

	#---> RETURN VALUE <---#
	return $interleaved;

}
#-----------------------#
sub read_dot_bracket_structures_to_hash {

	#---> PASS ARGUMENTS ---#
	my $dot_bracket_file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing miRNA dot-bracket structures..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %structures;

	#---> BODY <---#

		#---> Open dot-bracket structures file <---#
		open DB, "<", $dot_bracket_file or die "[ERROR] File '$dot_bracket_file' could not be opened! Execution aborted.\n";

		#---> Push miRBase hairpin dot-bracket structure records to hash <---#
		while (<DB>) {

			# Remove trailing newline character
			chomp;

			# Split record by tab characters
			my ($short_name, $dot_bracket) = split /\t/;

			## Add record to hash (warn if ID - i.e. miRBase hairpin short name - has already been seen!)
			unless (exists $structures{$short_name}) {
				$structures{$short_name} = $dot_bracket;
			}
			else {
				warn "[WARNING] Multiple entries for miRNA hairpin '$short_name'. Only first entry kept.\n";
			}

		}

		#---> Close dot-bracket structures file <---#
		close DB;

	#---> STATUS MESSAGE <---#
	print STDERR "Processed miRNA hairpin structures." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%structures;

}
#-----------------------#
sub print_selected_keys_or_all_hash {

	#---> PASS ARGUMENTS ---#
	my $hash_ref = shift;
	my $array_ref = shift;
	my $all = shift;
	my $outfile = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Printing dot-bracket structures..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $fh;

	#---> BODY <---#

		#---> Open output filehandle <---#
		if ($outfile) {
			open ($fh, '>', $outfile) or die;
		}
		else {
			open ($fh, '>&', \*STDOUT) or die;
		}


		#---> If $all is set: Print whole hash in alphanumeric order of keys <---#
		if ($all) {
			foreach (sort { $a cmp $b } keys %$hash_ref) {
				print $fh $_ . "\t" . $$hash_ref{$_} . "\n";
			}
		}

		#---> Else: Print key -> value pairs in the order their keys appear in the array <---#
		elsif (defined $array_ref) {
			foreach (@$array_ref) {
				if (exists $$hash_ref{$_}) {
					print $fh $_ . "\t" . $$hash_ref{$_} . "\n";
				}
			}
		}

		#---> Close output filehandle <---#		
		close $fh;

	#---> STATUS MESSAGE <---#
	print STDERR "Dot-bracket structures printed." . "\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#

