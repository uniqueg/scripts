#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: total_genome_coverage_fractions.pl
### Created: Jun 3, 2013
### Modified: Jun 3, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: For a mapped sequencing library in BED format, computes the position-wise coverage fractions for each base in the genome. The BED file has to be sorted by chromosome and position! 
### Output: One BED file for each chromosome; files contain the base-wise coverage fractions in the score columns and have the same format (half/end-open) as the input BED file
### Usage: perl ./total_genome_coverage_fractions.pl for information on required and optional arguments
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

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = 0;
my $open = 1;
my $chunk = 100000;
my $bed = '';
my $zeros = '';
my $size_fa = 1;
my $out_fa = '';
my $input = '';
my $fasta = '';
my $output = '';
my $options_result = GetOptions (
	'usage|help'	=> \$usage,
	'quiet'			=> \$quiet,
	'verbose'		=> sub { $quiet = 0 },
	'open'			=> \$open,
	'closed'		=> sub { $open = 0 },
	'chunk-size=i'	=> \$chunk,
	'bed!'			=> \$bed,
	'zeros!'		=> \$zeros,
	'fasta-seq'		=> sub { $size_fa = 0 },
	'fasta-size'	=> \$size_fa,
	'fasta-out=s'	=> \$out_fa,
	#---------------------------------------#
	'input=s'		=> \$input,
	'fasta=s'		=> \$fasta,
	'output=s'		=> \$output
);
# Show usage and exit if usage information was requested or if 'GetOptions' function failed
die &usage() if ($usage || !$options_result);
# Show usage and exit if required arguments are not set
die &usage() if (!$input || !$fasta || !$output);
## Show usage and exit if mutually exclusive options are selected
die &usage() if (); 

#---> GLOBAL VARIABLES <---# 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print STDERR "Starting 'total_genome_coverage_fractions.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $files_exist;
my $total_read_length;
my $chr_sizes_hash_ref;

#---> BODY <---#
# Validate input files
$files_exist = validate_files($input, $fasta);
# Print error message, show usage and exit if one of the supplied input files does not exist
die "The following file was not found: '$files_exist'" if $files_exist;
# Calculate total read length
$total_read_length = &total_length_bed($input, $open);
# Load FASTA file into hash
$chr_sizes_hash_ref = &fasta_to_hash($fasta);
# UNLESS provided FASTA file contains sequence lengths... 
unless ($size_fa) {
	# ...print warning message
	print STDERR "[WARNING] The conversion of the sequences of a FASTA file to their sizes/lengths is not implemented efficiently for the sake of modularity/subroutine reusability. For repeated use, transform sequences separately once, e.g. by using the Perl script 'fasta_seq_to_fasta_size.pl' by Alexander Kanitz.\n";
	# ...and replace sequences in hash with sequence lengths	
	$chr_sizes_hash_ref = &seq_hash_to_seq_sizes_hash($chr_sizes_hash_ref)
	#
}
# IF corresponding option selected, save sequence sizes to FASTA file
&hash_to_fasta_file($chr_sizes_hash_ref, $out_fa) if ($out_fa);
# Calculate coverage and print to output file
my $temp = &calc_coverage_fractions($input, $chunk, $open, $chr_sizes_hash_ref, $total_read_length, $output);

## Traverse through each element $array_ref of array @$temp
foreach my $array_ref ( @$temp ) {
	# Print array @$arra_ref
	foreach (@$array_ref) {
		print "$_\n";
	}
}

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit;
#================#
#    MAIN END    #
#================#



#===================#
#   SANDBOX START   #
#===================#

	###########
	## TO DO ##
	###########
	# 
	# POPULATE ARRAYS		take care of overhead
	#
	# WRITE OUTPUT SUB 
	#  
	# INCLUDE OPTIONS: 3. INPUT SORTED; 4. INPUT: SAM/BAM; 5. OUTPUT DIR 
	#
	# DIE IF SEQUENCE SIZE IS NOT IN HASH
	#
	###########	

#===================#
#    SANDBOX END    #
#===================#



#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current program
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./total_genome_coverage_fractions.pl [OPTIONS] --input [FILE] --fasta [FILE] --output [FILE]
==================================================
Required arguments:
--input			Input file in sorted(!) BED format [FILE]
--fasta			Filename of FASTA file of sequences or sizes (see --fasta-seq and --fasta-size) [FILE]; ID lines have to match sequence names that appear in --input (e.g. ">chr1" for "chr1")
--output		Name for output file [FILE]
==================================================
Optional arguments:
--usage|--help		Show this information and die
--verbose		Print status messages [DEFAULT]
--quiet			Shut up!
--open			BED file in half/end-open format [DEFAULT]
--closed		BED file in end-closed format (i.e. end position is included in sequence!)
--chunk-size		Number of positions to be written to the output file at a time [INTEGER | DEFAULT: 100,000]
--bed|--nobed		Output is written in BED format (--bed) or a simple TAB format of the following structure: "POSITION/TAB/SCORE"(--nobed) [DEFAULT: --bed] 
--zeros|--nozeros	Lines with zeros are included|excluded in the output [DEFAULT: --zeros]
--fasta-size		FASTA file contains sequence lengths [DEFAULT]
--fasta-seq		FASTA file contains sequences (conversion to length required; memory-intensive for large sequences/genomes!)
--fasta-out		If indicated with argument, saves the sequence sizes to FASTA file [FILE | DEFAULT: empty string; inactive]
';
}
#-----------------------#
sub validate_files {
### Function: Assures the existence of files
### Accepts: 1./2./... Filenames of files to be checked for their existence
### Returns: '0' if all files exist, else the name of the first found non-existing file 
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my @filenames = @_;

	#---> STATUS MESSAGE <---#
	print STDERR "Verifiying the existence of specified input files...\n" unless $quiet;
	
	#---> BODY <---#
	## Traverse through each element $file of array @filenames
	foreach my $file ( @filenames ) {
    	## Check file existence
		unless (-e $file) {
 			# Return filename $file if corresponding file does not exist
 			return $file;
 		} 
	}
	
	#---> STATUS MESSAGE <---#
	print STDERR "The existence of the specified input files has been verified!\n\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return 0;
}
#-----------------------#
sub total_length_bed {
### Function: Calculates the sum of the lengths of all regions in a BED file.
### Accepts: 1. Name of input BED file [FILE]; 2. Boolean (0 or 1) indicating whether the input BED file is in half/end-open format (default: 1)
### Returns: Sum of region lengths [INTEGER]
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;
	my $open = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Counting total length of all regions in BED file '$in_file'...\n" unless $quiet;			

	#---> SUBROUTINE VARIABLES <---#
	my $sum = 0;
	
	#---> BODY <---#
	# Open file handle BED
	open BED, "$in_file";
	
	## Traverse through file handle BED line by line
	while (<BED>) {
		# Split line by tabs
		my ($chr, $start, $end) = split /\t/, $_;
		# Add current region length to sum of lengths
		$sum += $end - $start + 1 - $open;	
	}
		
	# Close file handle BED
	close BED;
	
	#---> STATUS MESSAGE <---#
	print STDERR "BED file '$in_file' processed! Total length of all regions: $sum\n\n" unless $quiet;		
	
	#---> RETURN VALUE <---#
	return $sum;
}
#-----------------------#
sub fasta_to_hash {
### Function: Read FASTA file into hash; multiple sequence lines per entry, if present, are concatenated; no alphabet checking is performed
### Accepts: 1. FASTA filename
### Returns: Reference to sequence hash
### Dependencies: Subroutine 'check_key_value_pair'
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	

	#---> STATUS MESSAGE <---#
	print STDERR "Processing FASTA file '$in_file'...\n" unless $quiet;			

	#---> SUBROUTINE VARIABLES <---#
	my (%seq_hash, $id, $last_id_line);
	my $seq = "";
	
	#---> BODY <---#
	# Open file handle FA
	open FA, "$in_file";	
	#-----------------------#
	## Traverse through file handle FA line by line
	while (my $line = <FA>) {
		# Remove trailing newline character "\n"
		chomp $line;
		## Check whether current line is identifier line
		if ( $line =~ s/^>// ) {
			## Check whether identifier $id is defined and sequence $seq is not empty (the very first case!); if so... 
			if (defined $id && $seq ne "") {
				# Test whether key or key/value sequence pair is duplicate and generate unique $id (saved in $$result_array_ref[1]) if necessary
				my ($err_code, $new_id) = &check_key_value_pair($id, $seq, \%seq_hash);
				# Write warning message to STDERR if error code in $err_code is 1
				print STDERR "[WARNING] Sequence entry starting in line $last_id_line is a duplicate. Only one instance kept." . "\n" if $err_code == 1;
				# Write warning message to STDERR if error code in $err_code is 2				
				print STDERR "[WARNING] Identifier in line $last_id_line of '$in_file' not unique. Adding suffix." . "\n" if $err_code == 2;
				# Assign sequence $seq to hash key $new_id
				$seq_hash{$new_id} = $seq;				
			}
			# Set $id to current line, $last_id_line to current line number in $in_file and $seq to empty string ""	
			($id, $last_id_line, $seq) = ($line, $., "");
		}
		## Else...
		else {
		    # Append current line to existing sequence string $seq
		    $seq = $seq.$line;
		}
	}
	#-----------------------#
	## Add last entry to hash
	# Test whether key or key/value sequence pair is duplicate and generate unique $id (saved in $$result_array_ref[1]) if necessary
	my ($err_code, $new_id) = &check_key_value_pair($id, $seq, \%seq_hash);
	# Write warning message to STDERR if error code in $err_code is 1
	print STDERR "[WARNING] Sequence entry starting in line $last_id_line is a duplicate. Only one instance kept." . "\n" if $err_code == 1;
	# Write warning message to STDERR if error code in $err_code is 2				
	print STDERR "[WARNING] Identifier in line $last_id_line of '$in_file' not unique. Adding suffix." . "\n" if $err_code == 2;
	# Assign sequence $seq to hash key $new_id
	$seq_hash{$new_id} = $seq;		
	#-----------------------#
	# Close file handle FA
	close FA;
		
	#---> STATUS MESSAGE <---#
	print STDERR "FASTA File $in_file read!\n\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return \%seq_hash;
}
#-----------------------#
sub check_key_value_pair {
### Function: Tests whether a specified key already exists in a hash; if so, and if the values match (i.e. entry is a duplicate), prints warning message to STDERR, but only one instance is kept; if so, and if the values do not match, prints warning message to STDERR and appends (via an optional specified linker) to key
### Accepts: 1. Hash key to test; 2. Hash value to test; 3. Hash reference; 4. Linker [STRING | DEFAULT = "_"]
### Returns: 1. Error code (0: original key unique ('normal' case); 1: key/value pair is a duplicate; 2: original key not unique, linker and serial number appended); 2. Unique key
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($key, $value, $hash_ref, $linker) = (@_, "_");
	
	#---> BODY <---#
	## Check whether key with current identifier $key exists in hash %$hash_ref
	if ( defined $$hash_ref{$key} ) {
		## If so, check whether $key/$value pair is a duplicate (i.e. hash value of current identifier $key matches current $value)
		if ( $$hash_ref{$key} eq $value ) {
			# If so, return error code 1 and original key $key
			return (1, $key);
		}
		## Else...
		else {
			# Try to add suffix by initializing a serial number $suf_nr... 
			my $suf_nr = 1;
			# ...and increasing it as long as a hash key made of the current serial number $suf_nr appended to $key via linker $linker exists
			$suf_nr++ while defined $$hash_ref{$key . "_" . $suf_nr};
			# Append latest serial number $suf_nr to $key via linker $linker
			$key .= $linker . $suf_nr;
			# Return error code 2 and new, unique key $key
			return (2, $key);
		}
	}
	# Else...
	else {
		# Return error code 0 and original key $key
		return (0, $key);	
	}	
}
#-----------------------#
sub seq_hash_to_seq_sizes_hash {
### Function: Calculates sizes of sequence hash
### Accepts: 1. Sequence hash (keys = IDs, values = sequences)
### Returns: Sequence size hash (keys = IDs, values = sequence lengths)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my %seq_hash = %{shift()};	

	#---> STATUS MESSAGE <---#	
	print STDERR "Calculating sequence sizes...\n" unless $quiet;
	
	#---> SUBROUTINE VARIABLES <---#
	
	
	#---> BODY <---#
	## Traverse through sequence hash sorted by keys (i.e. sequence IDs); for each...
	foreach my $seq_id (sort keys %seq_hash) {
		# ...replace hash value (i.e. sequence with sequence length)
		$seq_hash{$seq_id} = length $seq_hash{$seq_id};
	}	
	
	#---> STATUS MESSAGE <---#
	print STDERR "Calculated sequence sizes!\n\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return \%seq_hash;
}
#-----------------------#
sub hash_to_fasta_file {
### Function: Write sequence hash in FASTA format to indicated output fileUniquely name elements of an array with serial numbers and an optional pre- or suffix; load to hash
### Accepts: 1. Reference to sequence hash; 2. Output filename [STRING]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($hash_ref, $out_file) = @_;

	#---> STATUS MESSAGE <---#
	print STDERR "Writing to file '$out_file'...\n" unless $quiet;	

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
	print STDERR "File '$out_file' written.\n\n" unless $quiet;	
}
#-----------------------#
sub calc_coverage_fractions {
### Function: Calculates the coverage fractions for each base on each strand from the specified input file in BED format for all the positions between positions 0 and the chromosome limits specified in the chromosome sizes hash; coverage fractions are calculated as the fractions of the coverage (number of reads spanning a given position) and the total read length passed as an argument to the function.
### Accepts: 1. Input BED file [FILE]; 2. Number of positions after which to print; 3. chromosome sizes hash with chromosome names (hash keys) that are compatible with those in 1. [HASH REFERENCE]; 4. total read length (i.e. the sum of the lengths of all regions in 1.) [INTEGER]; 5. output filename [FILE]
### Returns: Output file in BED format giving the coverage fractions in the score column
### Dependencies: ?????????????????
### Type: ????????????????
	#---> PASS ARGUMENTS ---#
	my ($bed, $chunk, $open, $chr_sizes_hash_ref, $total_len, $out_file) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$bed'...\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $chr_now;																# holds current chromosome
	my $chunk_last = 0;															# holds the maximum chromosome position of the previous chunk; is initialized to 0 for the first chunk
	my $chunk_now;																# holds the maximum chromosome position of the current chunk
	my $chr_size_limit;															# holds size limit of the current chromosome
	my (@cov_fract_plus, @cov_fract_minus, @cov_fract_star);					# arrays to hold the coverage fractions for the plus, minus and combined/undefined strand for the current chunk
	my @array_refs = (\@cov_fract_plus, \@cov_fract_minus, \@cov_fract_star);	# holds references to coverage fraction arrays
	my $incr = 1/$total_len;													# holds value to bed added to coverage fractions for each read covering a given position 
	
	#---> BODY <---#

	# Open BED file handle
	open BED, "$bed";
	## Crawl through BED file line by line
	while (<BED>) {	

		print $_;
		
		# Remove trailing new line characters
		chomp;
		# Split line into relevant coordinates by tabs
		my ($chr, $start, $end, $name, $score, $strand) = split /\t/;

		## IF current chromosome is defined...
		if (defined $chr_now) {
			## Test whether it is different from new chromosome; IF so...
			if ($chr_now ne $chr) {
				
print "NEW CHR!\n";
				
				# Write output
#!				&write_to_file();
				# Write rest of chromosome
#!				&write_rest_to_file();
				# Re-initialize array
				&initialize_arrays(0, ($chunk_now - $chunk_last + 1), \@array_refs);
				# Set the size limit for the new chromosome (subtract 1 because not 0-based)
				$chr_size_limit = $$chr_sizes_hash_ref{$chr} - 1;
				# Write current to previous chunk
				$chunk_last = $chunk_now;
				# Reset maximum chromosome position of current chunk (either to chromosome size limit or $chunk - 1)
				$chunk_now = ($chr_size_limit < $chunk - 1) ? $chr_size_limit : $chunk - 1;
				# Set current chromosome
				$chr_now = $chr;
			}
		}
		## ELSE (only in first iteration!)...
		else {

print "FIRST INSTANCE!\n";

			# Set the size limit for the new chromosome (subtract 1 because not 0-based)
			$chr_size_limit = $$chr_sizes_hash_ref{$chr};
			# Set current chromosome
			$chr_now = $chr;
			# Set maximum chromosome position of current chunk (either to chromosome size limit or $chunk - 1)
			$chunk_now = ($chr_size_limit < $chunk - 1) ? $chr_size_limit : $chunk - 1;
		}
		
		## Initialize arrays UNLESS they are already initialized (only in first iteration!)
		unless (@cov_fract_minus || @cov_fract_plus || @cov_fract_star) {

print "FIRST INSTANCE! JUST ONCE TOO!\n";

			&initialize_arrays(0, $chunk, \@array_refs);
		}

		## Check IF end of current read (subtract $open to account for potential half-/end-open format of BED) lies beyond chromosome size limit 
		if ($end - $open > $chr_size_limit) {
			# Print warning message
			print STDERR "[WARNING] Region '$name' on chromosome '$chr_now' is out of bounds. Coverage fractions for '$chr_now' will be printed between position 1 (0 in BED file) and position $chr_size_limit (as indicated by the provided FASTA file) only!\n\n";
			# Write output
#!			&write_to_file();
			# Re-initialize array
			&initialize_arrays(0, ($chunk_now - $chunk_last + 1), \@array_refs);
			# Write current to previous chunk
			$chunk_last = $chunk_now;
			# Reset maximum chromosome position of current chunk (either to chromosome size limit or to current value + $chunk - 1)
			$chunk_now = ($chr_size_limit < $chunk_now + $chunk - 1) ? $chr_size_limit : $chunk_now + $chunk - 1;
			# Skip remainder of iteration
			next;			
		}
		
		## Check IF current read is completely beyond size limit of current chunk
		if ($start > $chunk_now) {

print "WHEN DOES THIS HAPPEN?\nCHUNK NOW: $chunk_now\n";
			# Write output
#!			&write_to_file();
			# Re-initialize array
			&initialize_arrays(0, ($chunk_now - $chunk_last + 1), \@array_refs);
			# Write current to previous chunk
			$chunk_last = $chunk_now;
			# Reset maximum chromosome position of current chunk (either to chromosome size limit or to current value + $chunk - 1)
			$chunk_now = ($chr_size_limit < $chunk_now + $chunk - 1) ? $chr_size_limit : $chunk_now + $chunk - 1;			
		}

print $strand;

		## Populate coverage fraction arrays
		## Check IF current region is on the '+' strand		
		if ($strand eq "+") {
			## If so, for each position covered by the current region...
			for my $index ($start .. $end) {
				# ...increase corresponding array element value by $incr
				$cov_fract_plus[$index] += $incr;
				print "a\n";
			}
		}
		## Check ELSIF current region is on the '-' strand	
		elsif ($strand eq "-") {
			## For each position on the '+' strand covered by the current region...
			for my $index ($start .. $end) {
				# ...increase corresponding array element value by $incr
				$cov_fract_minus[$index] += $incr;				
			}			
		}
		## Check ELSIF current region has unrecognized strand information		
		elsif ($strand ne "*" || ! defined($strand)) {
			# Print warning message
			print STDERR "[WARNING] Region '$name' has unrecognized strand information '$strand'. Only '+', '-', '*' and undefined/empty (assumed '*') are recognized. Region skipped!\n\n";
			next;
		}
		## For each position on the '+', '-' or '*'/empty strand covered by the current region...
		for my $index ($start .. $end) {
			# ...increase corresponding array element value by $incr
			$cov_fract_star[$index] += $incr;			
		}	

## Print array @cov_fract_star
#foreach (@cov_fract_star) {
#	print "$_\n";
#}


	}




	#---> STATUS MESSAGE <---#
	print STDERR "File '$bed' processed!\n\n" unless $quiet;

	#---> RETURN VALUE <---#
	return \@array_refs;

}
#-----------------------#
sub initialize_arrays {
### Function: Populates a specified number of elements of one or more arrays with a specified variable (e.g. '0')  
### Accepts: 1. Value to populate the array(s) with; 2. number of array elements that shall be populated; 3. Reference to array of references to arrays that shall be populated
### Returns: Reference to array of references to populated arrays
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $value = shift;
	my $times = shift;
	my $array_refs = shift;
	
	#---> BODY <---#
	## Traverse through each element $array_ref of array @array_refs
	foreach my $array_ref ( @$array_refs ) {
	    # Empty array
	    @$array_ref = ();
		## Initialize array...
		for my $pos (0 .. ($times - 1)) {
			## ...by pushing $value to the first $times elements
			push @$array_ref, $value;
		}		
	}	
	
	#---> RETURN VALUE <---#
	return $array_refs;
}
#=======================#
#    SUBROUTINES END    #
#=======================#