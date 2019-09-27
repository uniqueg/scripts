#!usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: structure_at_miR_hairpin_position.pl
### Created: Nov 18, 2013
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
my $positions = '';
my $dot_bracket = '';
my $unpaired = '';
my $fr_unpaired = '';
my $loop_size = '';
my $loop_overlap = '';
my $dist_loop = '';
my $print_all = '';
my $output = '';
my $options_result = GetOptions (
	'h|usage|help' => \$usage,
	'v|quiet|no_verbose' => \$quiet,
	'p|positions=s' => \$positions,
	'k|dot-bracket=s' => \$dot_bracket,
	'u|unpaired' => \$unpaired,
	'f|fraction-unpaired' => \$fr_unpaired,
        's|loop-size' => \$loop_size,
        'l|loop-overlap' => \$loop_overlap,
	'd|distance-from-loop' => \$dist_loop,
	'a|all' => \$print_all,
	'o|output=s' => \$output
);

die $usage_info if $usage || !$options_result;

## Verify/process options/arguments
die "[ERROR] Required option '--positions' missing!\nExecution aborted.\n\n$usage_info" unless $positions;
die "[ERROR] Required option '--dot-bracket' missing!\nExecution aborted.\n\n$usage_info" unless $dot_bracket;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print STDERR "Starting 'structure_at_miR_hairpin_position.pl'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $pos_AoA_ref;
my $dot_bracket_hash_ref;

#---> BODY <---#
	
	# Get positions of interest
	$pos_AoA_ref = &miRNA_hairpin_positions_to_AoA($positions);

	# Get structures in dot-bracket format from '--dot-bracket' file
	$dot_bracket_hash_ref = &read_dot_bracket_structures_to_hash($dot_bracket);

	# Print structures
	&print_substructures($dot_bracket_hash_ref, $pos_AoA_ref, $output, $print_all, $unpaired, $fr_unpaired, $loop_size, $loop_overlap, $dist_loop);

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
'Usage: perl structure_at_miR_hairpin_position.pl [OPTIONS] --input [TAB] --dot-bracket [TAB]

Description: Returns the dot-bracket secondary substructure of a miRNA hairpin in a given region. Requires miRNA hairpin structures in dot-bracket notation in a .tab format (miRNA hairpin short name TAB dot-bracket structure) and a position file in a .tab format (miRNA hairpin short name as first part of ID - separated by spaces - TAB start position TAB end position TAB score TAB sequence). Additional information can be printed out (check options). Input files are generated e.g. by the scripts "miRBase_hairpin_short_name_to_dot_bracket_structure.pl" (Alexander Kanitz; structures) and ag-scan-with-PWM.pl (Andreas R. Gruber; positions). 

Options:
--positions [TAB]	Positions file. Required. See description for details.
--dot-bracket [TAB]	Structures file. Required. See description for details.
--unpaired		Also print out the number of unpaired nucleotides.
--fraction-unpaired	Also print out the fraction of unpaired nucleotides.
--loop-size		Also print out the size of the hairpin loop from first to last unpaired base.
--loop-overlap		Also print out the overlap (in nucleotides) between the motif and the loop.
--distance-from-loop	Also print out the distance (in nucleotides) between the motif center and the loop center.
--all			Print all extra information.
--output [FILE]		Output file (default: write to STDOUT).
--quiet			Do not print status messages to STDERR.
--help|usage		Print this screen.
';
}
#-----------------------#
sub miRNA_hairpin_positions_to_AoA {

	#---> PASS ARGUMENTS ---#
	my $pos_file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Extracting positions from file '$pos_file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @positions;

	#---> BODY <---#

		#---> Open positions file <---#
		open POS, "<", $pos_file or die "[ERROR] File '$pos_file' could not be opened! Execution aborted.\n";

		#---> Push miRBase hairpin short names from input file/stream to array <---#
		while (<POS>) {
			chomp;
			my @values = split /\t/;
			my @IDs = split / /, (shift @values), 3;
			unshift @values, @IDs;
			push @positions, \@values;
		}

		#---> Close positions file <---#
		close POS;

	#---> STATUS MESSAGE <---#
	print STDERR "File processed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@positions;

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
	print STDERR "Structures processed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%structures;

}
#-----------------------#
sub print_substructures {

	#---> PASS ARGUMENTS ---#
	my %str_hash = %{shift()};
	my @pos_AoA = @{shift()};
	my $outfile = shift;
	my $all = shift;
	my $unpaired = shift;
	my $fr_unpaired = shift;
	my $loop_size = shift;
	my $loop_overlap = shift;
	my $dist_loop = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Writing output..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $fh;

	#---> BODY <---#

		#---> Open output filehandle <---#
		if ($outfile) {
			open ($fh, '>', $outfile) or die "[ERROR] Could not open output file '$outfile' for writing!\nExecution aborted.\n";
		}
		else {
			open ($fh, '>&', \*STDOUT) or die "[ERROR] Could not write to STDOUT!\nExecution aborted.\n";
		}

		#---> Traverse through all positions <---#
		foreach my $pos_array_ref (@pos_AoA) {
			
			#---> Store field values in dedicated variables <---#
			my ($short_name, $hairpin_id, $long_name, $motif_start, $motif_end, $score, $seq) = @$pos_array_ref;

			#---> Extract substring from structure based on start and end coordinates <---#
			my $substructure = substr $str_hash{$short_name}, ($motif_start - 1), ($motif_end - $motif_start + 1);

			#---> Generate output string <---#
			my $out_line = join "\t", $short_name, $hairpin_id, $long_name, $motif_start, $motif_end, $score, $seq, $substructure;

			#---> Check for extra output options <---#			

				#---> I. Number and fraction of unpaired bases <---#
				if ($unpaired || $fr_unpaired || $all ) {
	
					#---> Calculate number and fraction of unpaired bases in substructure <---#
					my $unpaired_number = () = $substructure =~ /\./g;
					my $unpaired_fraction = $unpaired_number / length $substructure;
	
					#---> Append requested information to output string <---#
					$out_line .= "\t" . $unpaired_number if $unpaired || $all;
					$out_line .= "\t" . $unpaired_fraction if $fr_unpaired || $all;
	
				}

				#---> II. Loop size, overlap and distance from loop <---#
				if ($loop_size || $loop_overlap || $dist_loop || $all ) {
	
					#---> Initialize variables to hold the loop size, overlap and distance <---#
					my $size = "N/A";
					my $overlap = "N/A";
					my $distance = "N/A";

					#---> Check for presence of loops <---#
					my @loops = $str_hash{$short_name} =~ m/\(\.*\)/g;	

					#---> Assert that there is exactly one loop <---#
					if (scalar @loops == 1) {

						#---> Set loop start and end positions <---#
						my $loop_start = $-[0] + 2;	# first UNPAIRED base desired (+1) and 1-based coordinates used (+1)
						my $loop_end = $+[0] - 1;	# last UNPAIRED base desired (-2) and 1-based coordinates used (+1)

						#---> Calculate loop size <---#						
						$size = $loop_end - $loop_start + 1;

						#---> Calculate distance between loop center and motif center <---#
						$distance = ( ( $motif_end + $motif_start ) / 2 ) - ( ($loop_start + $loop_end) / 2 );
						
						#---> Check for overlap <---#
						#---> Conditions: Start or end of motif falls within loop or whole loop falls within motif <---#
						if ( &num_within_range($motif_start, $loop_start, $loop_end) || &num_within_range($motif_end, $loop_start, $loop_end) || &range_within_range($loop_start, $loop_end, $motif_start, $motif_end) ) {

							#---> Set overlap start and end positions <---#
							my $overlap_start = &num_within_range($motif_start, $loop_start, $loop_end) ? $motif_start : $loop_start;
							my $overlap_end = &num_within_range($motif_end, $loop_start, $loop_end) ? $motif_end : $loop_end;
			
							#---> Calculate overlap <---#
							$overlap = $overlap_end - $overlap_start + 1;
						}
						else {
							$overlap = 0;
						}

					}
					elsif (scalar @loops > 1) {
						warn "[WARNING] More than one loop found in dot-bracket structure of miRNA hairpin '$short_name'. Distance, loop size and loop overlap set to 'N/A'.\n";
					}
					else {
						warn "[WARNING] No loop found in dot-bracket structure of miRNA hairpin '$short_name'. Distance, loop size and loop overlap set to 'N/A'.\n";
					}
	
					#---> Append requested information to output string <---#
					$out_line .= "\t" . $size if $loop_size || $all;
					$out_line .= "\t" . $overlap if $loop_overlap || $all;
					$out_line .= "\t" . $distance if $dist_loop || $all;
	
				}

			#---> Print outputOpen output filehandle <---#
			print $fh $out_line . "\n";			
	

		}

		#---> Close output filehandle <---#		
		close $fh;

	#---> STATUS MESSAGE <---#
	print STDERR "Output written." . "\n" unless $quiet;

}
#-----------------------#
sub num_within_range {
# Returns 1 if [m,n] falls within [x,y], else returns 0; assumes x < y and swaps values if necessary

	my ($n, $x, $y) = @_;

	($x, $y) = ($y, $x) if ($x > $y);

	return 1 if $n >= $x && $n <= $y;

	return 0;

}
#-----------------------#
sub range_within_range {
# Returns 1 if [m,n] falls within [x,y], else returns 0; assumes x < y and m < n and swaps values if necessary

	my ($m, $n, $x, $y) = @_;

	($m, $n) = ($n, $m) if ($m > $n);
	($x, $y) = ($y, $x) if ($x > $y);

	return 1 if $m >= $x && $n <= $y;

	return 0;

}
#=======================#
#    SUBROUTINES END    #
#=======================#
