#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: generate_query_transcript_fasta_file.pl
### Created: Mar 6, 2013
### Modified: Mar 6, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: The script generates a set of sequences, each made up of a fixed upstream sequence fragment and a variable downstream sequence fragment, with the latter consisting of all possible, unique N-mers of a given query sequence. 
### Output: A FASTA file containing the set of sequences
### Usage: perl ./generate_query_transcript_fasta_file.pl for information on required and optional arguments
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
## Set default values for options and arguments
my $usage = '';
my $quiet = '';
my $fixed = '';
my $nt = 6;
my $query = '';
my $output = '';
## Parse options
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'fixed=s' => \$fixed,
	'nt=i' => \$nt,
	#-----------------------#
	'query=s' => \$query,
	'output=s' => \$output
);
# Die with usage information if --usage option is set or if options parsing returned FALSE
die $usage_info if $usage || !$options_result;
# Die with usage information if required arguments are not set
die $usage_info if !$query || !$output; 
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'generate_query_transcript_fasta_file.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my ($array_ref, $hash_ref);

#---> BODY <---#
# Find unique N-mers within query sequence
$array_ref = &find_unique_n_mers_in_seq($query, $nt);
# Concatenate fixed sequence to 5'-end of each of the resulting sequences
$array_ref = &cat_adapter_to_seq_array($array_ref, $fixed, 5);
# Add unique names to each sequence
$hash_ref = &name_array_elements_with_serial_number_and_add_to_hash($array_ref, length($fixed)."nt_".$nt."nt_", "prefix");
# Write to FASTA file
&print_hash_to_fasta_file($hash_ref, $output);

#---> STATUS MESSAGE <---#
print "\nDone.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current program
### Accepts: n/a
### Returns: String with usage information
### Dependencies: n/a
### Type: Specialized
'Usage: perl ./generate_query_transcript_fasta_file.pl [OPTIONS] --fixed [NUCLEOTIDE SEQUENCE | STRING | DEFAULT = EMPTY STRING ""] --query [NUCLEOTIDE SEQUENCE | STRING] --nt [N-MER SIZE | INTEGER| DEFAULT = 6] --output [FILE]
==================================================
Required arguments:
--query			Sequence to be split up into N-mers which are then concatenated to the sequence given in --fixed
--output		Output filename
==================================================
Optional arguments:
--fixed			Fixed sequence fragment present at the 5-end of each sequence in the output
--nt 			Size of N-mers generated from --query
--usage|help		Show this information and die
--quiet			Shut up!
';
}
#-----------------------#
sub find_unique_n_mers_in_seq {
### Function: Extracts all possible, unique subsequences of a specified size from a a given sequence
### Accepts: 1. Query sequence [STRING]; 2. N-mer size [INT] 
### Returns: Reference to array of subsequences
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($seq, $n) = @_;

	#--- SUBROUTINE VARIABLES ---#
	my @seqs;

	#--- BODY ---#
	# Die if query sequence is shorter than N-mer size
	die "N-mer size in '--nt' bigger than length of sequence in '--query'!\nNo output can be generated.\nDied" if length($seq) < $n;
	## Generate subsequences by looping through query sequence nucleotide by nucleotide
	for ( my $start_pos = 0 ; $start_pos <= length($seq) - $n + 1; $start_pos++ ) {
	    # Extract subsequence and push to array
	    push @seqs, substr $seq, $start_pos, $n;
	}
	my $seqs_sorted_array_ref = &unique_array_elements(@seqs);

	#---> STATUS MESSAGE <---#
	print "Extracted all possible, unique subsequences from sequence '$seq'.\n" unless $quiet;

	#--->  RETURN VALUE  <---#
	return $seqs_sorted_array_ref;
}
#-----------------------#
sub cat_adapter_to_seq_array {
### Function: Concatenate an adapter sequence to the 5'- or 3'-end of a set of sequences stored in an array
### Accepts: 1. Reference to sequence array; 2. Adapter sequence [STRING]; 3. Orientation [INT; ALLOWED = 3, 5 | DEFAULT = 5]
### Returns: Reference to array of sequences
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($seq_array_ref, $adapter, $orientation) = (@_, "");
	## Check whether orientation argument is illegal
	if ($orientation ne "3" && $orientation ne "5" ) {
		## If so, print status message and use default value
		print "Illegal value for orientation of adapter concatentation. Default value (5'-end addition) is used.\n" unless $quiet;
		$orientation = 5;		
	} 

	#--- SUBROUTINE VARIABLES ---#
	my @seqs;
	
	#--- BODY ---#
	## Traverse through each element $seq of array @$seq_array_ref
	foreach my $seq ( @$seq_array_ref ) {
	    # Push concatenated sequence to array (order depending on value in $orientation) 
		push @seqs, ($orientation == 5) ? $adapter.$seq : $seq.$adapter;
	}	
	
	#---> STATUS MESSAGE <---#
	print "Added sequence '$adapter' to the $orientation'-end of each sequence.\n" unless $quiet;	
	
	#--->  RETURN VALUE  <---#
	return \@seqs;	
}
#-----------------------#
sub name_array_elements_with_serial_number_and_add_to_hash {
### Function: Uniquely name elements of an array with serial numbers and an optional pre- or suffix; load to hash
### Accepts: 1. Reference to array; 2. Pre- or -suffix [STRING; DEFAULT = EMPTY STRING ""]; 3. Orientation [STRING; ALLOWED = "prefix", "suffix" | DEFAULT = "prefix"]
### Returns: Hash reference
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($array_ref, $fix, $orientation) = (@_, "", "");
	## Check whether orientation argument is illegal
	if ($orientation ne "prefix" && $orientation ne "suffix" ) {
		## If so, print status message and use default value
		print "Illegal value for orientation of concatentation. Default value ('prefix') is used.\n" unless $quiet;
		$orientation = "prefix";		
	} 

	#--- SUBROUTINE VARIABLES ---#
	my %hash;
	
	#--- BODY ---#
	# Determine necessary number of digits for unique numbering
	my $digits = length($#{$array_ref});
	## 
	for ( my $i = 1; $i <= $#{$array_ref}; ++$i ) {
	    # Generate unique name $id
	    my $id = ($orientation eq "prefix") ? $fix.sprintf("%0${digits}d", $i) : sprintf("%0${digits}d", $i).$fix;
		# Load hash (key: $id; value: array element)
		$hash{$id} = $$array_ref[$i-1];  
	}
	
	#---> STATUS MESSAGE <---#
	print "Added unique names.\n" unless $quiet;	
	
	#--->  RETURN VALUE  <---#
	return \%hash;	
}
#-----------------------#
sub print_hash_to_fasta_file {
### Function: Write sequence hash in FASTA format to indicated output fileUniquely name elements of an array with serial numbers and an optional pre- or suffix; load to hash
### Accepts: 1. Reference to sequence hash; 2. Output filename [STRING]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($hash_ref, $out_file) = @_;

	#--- BODY ---#
	# Open file handle OUT
	open OUT, ">$out_file";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $seq (sort {$a cmp $b} keys %{$hash_ref}) {
		# Print in fasta format
		print OUT ">" . $seq . "\n" . ${$hash_ref}{$seq} . "\n";
	}		
	# Close file handle OUT
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print "Sequences written to file '$out_file'.\n" unless $quiet;	
}
#-----------------------#
sub unique_array_elements {
### Function: Removes duplicates of array; preserves order 
### Accepts: Array
### Returns: Array reference
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my @seqs = @_;

	#--- SUBROUTINE VARIABLES ---#
	my %uniq;
	
	#--- BODY ---#	
	my @sorted = grep !$uniq{$_}++, @_;

	#--->  RETURN VALUE  <---#
	return \@sorted;
}
#=======================#
#    SUBROUTINES END    #
#=======================#