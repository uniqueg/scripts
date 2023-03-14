## PRINT HASH OF ARRAYS
foreach my $qu (keys %query) {
	print "$qu\t@{$query{$qu}}\n";
}
foreach my $tar (keys %target) {
	print "$tar\t@{$target{$tar}}\n";
}

my @array4 = ( "id1", "start1", "stop1", "name1", "id2", "start2", "stop2", "name2", "id3", "start3", "stop3", "name3" );
my @array3 = ( "id1", "start1", "stop1", "id2", "start2", "stop2", "id3", "start3", "stop3" );
my @array2 = ( "start1", "stop1", "start2", "stop2", "start3", "stop3" );

my @a = @{&everyNthElementOfArray(\@array4, 4)};

print "\n@a";

sub bedToHashOfArrays {
###
###
###
	#
	my $bed_file = shift;
	#
	my %HoA;
	#
	open BED, $bed_file;
	#
	while (<BED>) {
		chomp;
		my @line = split /\t/, $_, 4;
		unless (exists $HoA{$line[0]}) {
			$HoA{$line[0]} = [ $line[3], @line[1..2] ];
		}
		else {
			push @{ $HoA{$line[0]} }, $line[3], @line[1..2];
		}
	}
	#
	close BED;
	#
	return \%HoA;
}

sub everyNthElementOfArray {			### AoAs!!!
	my @array = @{shift()};#( "id1", "start1", "stop1", "id2", "start2", "stop2", "id3", "start3", "stop3" );
	my $no = shift;
	my @out;
	for (my $i = $no; $i >= 1; $i--) {
		push @out, @array[grep {!(($_ + $i) % $no)} 0..$#array];
	}
	print "@out\n";
	return \@out;
}


#foreach my $key (keys %query) {
#	print "$key\t@{$query{$key}{'rest'}}\n";
#	print "$key\t$query{$key}{'rest'}[0]\n";
#}
#my %target = %$target_HoHoA_ref;
#foreach my $key (keys %target) {
#	print "$key\t@{$target{$key}{'start'}}\n";
#	print "$key\t$target{$key}{'rest'}[0]\n";
#}

sub fastaToHash {
### Pass file in fasta format and load into hash
### Keys: ID lines (minus leading '>')
### Values: Sequences (multiple lines are concatenated, separating newline characters - if present - are removed) 
	# Pass arguments
	my $file = shift;
	# Screen output
	print "Processing fasta file \"$file\".\n";
	# Declare initialize variables
	my %hash;
	my $seq = "";
	my $id;
	# Open fasta file handle
	open FA, "$file";
	## Crawl through fasta file line by line
	while (<FA>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block if identifier line (contains leading '>')
		if ($line =~ s/^>//) {
			# Unless empty string, pass sequence ($seq) to hash (key = $id)
			$hash{$id} = $seq if $seq ne "";
			# Set identifier variable ($id)
			$id = $line;
			# Empty sequence string
			$seq = "";
		}
		## Else concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	# Add last entry to hash
	$hash{$id} = $seq;
	# Close fasta file handle
	close FA;
	# Screen output
	print "File \"$file\" processed.\n";
	# Return hash reference
	return \%hash;
}	


sub lookupID {
	# Pass arguments
	my $line = shift;
	my %id_hash = %{shift()};
	# Split line by tabs
	my @entries = split /\t/, $line; 
	# Change all letters to uppercase
	$entries[0] = uc $entries[0];
	# Obtain ID from separate file
	$entries[0] = $id_hash{$entries[0]}; 
	# Re-join line by tab
	$line = join "\t", @entries;
	# Return line
	return "$line\n"; 
}

sub fastaToHashRev {
### Pass file in fasta format and load into hash
### Keys: ID lines (minus leading '>')
### Values: Sequences (multiple lines are concatenated, separating newline characters - if present - are removed) 
	# Pass arguments
	my $file = shift;
	# Screen output
	print "Processing file \"$file\".\n";
	# Declare initialize variables
	my %hash;
	my $seq = "";
	my $id;
	# Open fasta file handle
	open FA, "$file";
	## Crawl through fasta file line by line
	while (<FA>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block if identifier line (contains leading '>')
		if ($line =~ s/^>//) {
			# Unless empty string, pass sequence ($seq) to hash (key = $id)
			$hash{$seq} = $id if $seq ne "";
			# Set identifier variable ($id)
			$id = $line;
			# Empty sequence string
			$seq = "";
		}
		## Else concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	# Add last entry to hash
	$hash{$seq} = $id;
	# Close fasta file handle
	close FA;
	# Screen output
	print "File processed.\n";
	# Return hash reference
	return \%hash;
}




#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	04-JAN-2013
### Modified:	10-JAN-2013
#######

#######
### FUNCTION:
### ---------
### 
#######

#######
### ARGUMENTS:
### ----------
### 1. 
### 2. 
### 3. 
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### 
#######

#######
### USAGE:
### ------
### perl /path/to/wm_in_fasta.pl <path/to/wm_input_file> </path/to/fasta_file> <top_hits/integer> <positional_separator/char> <nucleotide_separator/char> <path/to/output_file>
### Examples:
### 1. perl 
### 2. perl  
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Read weight matrices
### C. Read fasta file
### D. Slide window and calculate scores
### E. Clean-up
#######

### A. Pre-requisites
use warnings;
use strict;
use Benchmark;
my $wm_file = shift;
my $fa_file = shift;
my $no_hits = 1000;#shift;
my $pos_sep = "/";#shift;
my $nt_sep = ",";#shift;
my $out_file = "out_test";#shift;
# Print time
my @t = localtime(time);
print "STARTED: $t[2]:$t[1]:$t[0]\n";
###

### B. Read weight matrices
# Open weight matrix file handle
open WM, $wm_file;
# Open output file handle
open OUT, ">$out_file";
# Traverse through weight matrices one by one
while (<WM>) {
	# Assign read line to identifier variable $id_line
	my $id_line = $_;
	print OUT $id_line;
	# Read line containing sequence
	my $wm_line = <WM>;
	## Remove trailing new line characters from input lines
	chomp ($id_line);
	chomp ($wm_line);
	# Load 2-dimensional array (sequence position, nucleotide probability)
	my @wm_plus = &loadAoAfromString($wm_line, $pos_sep, $nt_sep);
	# Compute complement of weight matrix
	my @wm_minus = @{&wmComplement(\@wm_plus)};

	### C. Read fasta file
	# Read fasta file handle
	open FA, $fa_file;
	# Traverse through fasta file 
	while (<FA>) {
		# Remove trailing new line character
		chomp;
		# Assign read line to chromosome variable $chr
		my $chr = $_;
		# Remove leading ">" character
		$chr = substr($chr,1);
		# Read line containing sequence
		my $seq_line = <FA>;
		# Remove trailing new line character
		chomp($seq_line);
		# Transform lowercase to uppercase and split into array
		my @seq = split //, uc($seq_line);
		# Index chromosome (no requirement for conditional statements in sliding window loop)
		foreach (@seq) {
			if		($_ eq "A")	{$_ = 0}
			elsif	($_ eq "C")	{$_ = 1}
			elsif	($_ eq "G")	{$_ = 2}
			elsif	($_ eq "T")	{$_ = 3}
			else				{$_ = 4}
		}
			
		### D. Slide window and calculate scores
		for (my $start = 0; $start <= (@seq - @wm_plus); $start++) {
			my $end = $start + $#wm_plus;
			my @window = @seq[$start..$end];
			my $score = &wmLogScoreSum(\@window, \@wm_plus);
			if ($score ne "FALSE") {
				print OUT $chr . "\t" . ($start + 1) . "\t" . ($start + @wm_plus) . "\t" . "+" . "\t" . $score . "\n";
			} 
			$score = &wmLogScoreSum(\@window, \@wm_minus);
			if ($score ne "FALSE") {
				print OUT $chr . "\t" . ($start + 1) . "\t" . ($start + @wm_plus) . "\t" . "-" . "\t" . $score . "\n";
			} 		 
		}
		print "$id_line, $chr processed.\n"
	}
	# Close fasta file handle
	close FA;
}
# Close output file handle
close OUT;
# Close weight matrix file handle
close WM;

### E. Clean-up
# Print time
@t = localtime(time);
print "FINISHED: $t[2]:$t[1]:$t[0]\n";
# Print status message
print "Results written to file: $out_file\n";
# Exit program
exit;
###

### F. Subroutines

## E1. Load array of arrays from string with two separators 
sub loadAoAfromString {
	my $string = shift;
	my $sep1 = shift;
	my $sep2 = shift;
	my @tmp = split $sep1, $string;
	my @AoA;
	foreach (@tmp) {
		my @logs = split $sep2, $_;
		# Calculate the logs of each probability
		foreach (@logs) {
			$_ = log($_);
		}
		# Add value for other letters than ACGT (e.g. N; arbitrary value: ln of 0.25)
		push @logs, -1.386294361;						
		push @AoA, [ @logs ];
	}
	return @AoA;
}
##

## E2. Compute complement of weight matrix				################# DO REVERSE!!!
sub wmComplement {
	my $wm_ref = shift;
	my @wm = @$wm_ref;
	my @wm_comp = map { [@$_] } @wm;
	for (my $i = 0; $i <= $#wm_comp; $i++) {
		my $tmp = $wm_comp[$i][0];
		$wm_comp[$i][0] = $wm_comp[$i][3];
		$wm_comp[$i][3] = $tmp;
		$tmp = $wm_comp[$i][1];
		$wm_comp[$i][1] = $wm_comp[$i][2];
		$wm_comp[$i][2] = $tmp;
	}
	return \@wm_comp;
}
##

## E3. Compute log likelihood of weight matrix in window
sub wmLogScoreSum {
	my $win_ref = shift;
	my $wm_ref = shift;
	my $score = 0;
	my @win = @$win_ref;
	for (my $pos = 0; $pos < @win; $pos++) {
		$score += ${${$wm_ref}[$pos]}[$win[$pos]];
		if ($score < -18) {
			return "FALSE"
		}
	}
	return $score;
}
###

### TO DO
# only keep top hits ($no_hits) -> use hash, see notes ???
# add MIN SCORE to parameters
# remove all unnecessary stuff
# Sequences must be in 1 line!



#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	04-JAN-2013
### Modified:	10-JAN-2013
#######

#######
### FUNCTION:
### ---------
### 
#######

#######
### ARGUMENTS:
### ----------
### 1. 
### 2. 
### 3. 
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### 
#######

#######
### USAGE:
### ------
### perl /path/to/wm_in_fasta.pl <path/to/wm_input_file> </path/to/fasta_file> <top_hits/integer> <positional_separator/char> <nucleotide_separator/char> <path/to/output_file>
### Examples:
### 1. perl 
### 2. perl  
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Read weight matrices
### C. Read fasta file
### D. Slide window and calculate scores
### E. Clean-up
#######

### A. Pre-requisites
use warnings;
use strict;
## Print time
my @t = localtime(time);
print "STARTED: $t[2]:$t[1]:$t[0]\n";
## Parameter initialization
my $wm_file = "wm_test.fa";#shift;
my $fa_file = "seq_test.fa";#shift;
my $no_hits = 1000;#shift;
my $def_prob = "0.25";#shift;
my $pos_sep = "/";#shift;
my $nt_sep = ",";#shift;
my $out_file = "out_test_2";#shift;
# Calculate ln of default probability
$def_prob = log($def_prob);
###

### B. Read fasta file
# Open fasta file handle
open FA, $fa_file;
# Initialize sequence identifier @id_seq
my $id_seq;
# Initialize sequence array @seq
my @seq;
# Crawl through genome fasta file line by line
while (<FA>) {
	# Remove trailing new line character from input line
	chomp;
	## If line starts with ">" (identifier line), extract identifier $id_seq and empty sequence array @seq; else, transform all letters to uppercase, push to array and replace with numbers (see below) 
	if (substr($_, 0, 1) eq ">") {
		if ( @seq) {
			
		}
		$id_seq = substr $_, 1;
		undef(@seq);
	}
	else {
		# Transform lowercase to uppercase and push to array
		push @seq, split //, uc($_);
		## Index sequence (i.e. replace letters with numbers corresponding to wm positons; no requirement for conditional statements in sliding window loop)
		foreach (@seq) {
			if		($_ eq "A")	{$_ = 0}
			elsif	($_ eq "C")	{$_ = 1}
			elsif	($_ eq "G")	{$_ = 2}
			elsif	($_ eq "T")	{$_ = 3}
			else				{$_ = 4}
		}
	}
	
	### C. Read weight matrices
	# Open weight matrix file handle
	open WM, $wm_file;
	# Open output file handle
	open OUT, ">$out_file";
	# Traverse through weight matrices one by one
	while (<WM>) {
		# Assign read line to identifier variable $id_wm
		my $id_wm = $_;
		print OUT $id_wm;
		# Read line containing sequence
		my $wm_line = <WM>;
		## Remove trailing new line characters from input lines
		chomp ($id_wm);
		chomp ($wm_line);
		# Load 2-dimensional array (sequence position, nucleotide probability)
		my @wm_plus = &loadAoAfromString($wm_line, $pos_sep, $nt_sep, $def_prob);
		# Compute complement of weight matrix
		my @wm_minus = @{&wmComplement(\@wm_plus)};
		# Calculate length of wm array (starting from 0, i.e. sequence length - 1; to be used in window sliding)
		my $wm_len = $#wm_plus;
	
		### D. Slide window and calculate scores
		for (my $start = 0; $start <= (@seq - @wm_plus); $start++) {
			my @test = @seq[$start..($start + $wm_len)];
			my $score = &wmLogScoreSum(\@test, \@wm_plus);
			if ($score ne "FALSE") {
				print OUT $id_seq . "\t" . ($start + 1) . "\t" . ($start + @wm_plus) . "\t" . "+" . "\t" . $score . "\n";
			} 
			$score = &wmLogScoreSum(\@test, \@wm_minus);
			if ($score ne "FALSE") {
				print OUT $id_seq . "\t" . ($start + 1) . "\t" . ($start + @wm_plus) . "\t" . "-" . "\t" . $score . "\n";
			} 		 
		}
		print "$id_wm, $id_seq processed.\n"
	}
	# Close weight matrix file handle
	close WM;	
}
# Close output file handle
close OUT;
# Close fasta file handle
close FA;

### E. Clean-up
## Print time
@t = localtime(time);
print "FINISHED: $t[2]:$t[1]:$t[0]\n";
# Print status message
print "Results written to file: $out_file\n";
# Exit program
exit;
###

### F. Subroutines

## E1. Load array of arrays from string with two separators; transform wm probabilites to log space 
sub loadAoAfromString {
	my $string = shift;
	my $sep1 = shift;
	my $sep2 = shift;
	my $default = shift;
	my @tmp = split $sep1, $string;
	my @AoA;
	foreach (@tmp) {
		my @logs = split $sep2, $_;
		# Calculate the logs of each probability
		foreach (@logs) {
			$_ = log($_);
		}
		## Add value for other letters than ACGT (e.g. N; arbitrary value: ln of 0.25)
		push @logs, $default;						
		push @AoA, [ @logs ];
	}
	return @AoA;
}
##

## E2. Compute complement of weight matrix				################# DO REVERSE!!!
sub wmComplement {
	my $wm_ref = shift;
	my @wm = @$wm_ref;
	my @wm_comp = map { [@$_] } @wm;
	for (my $i = 0; $i <= $#wm_comp; $i++) {
		my $tmp = $wm_comp[$i][0];
		$wm_comp[$i][0] = $wm_comp[$i][3];
		$wm_comp[$i][3] = $tmp;
		$tmp = $wm_comp[$i][1];
		$wm_comp[$i][1] = $wm_comp[$i][2];
		$wm_comp[$i][2] = $tmp;
	}
	return \@wm_comp;
}
##

## E3. Compute log likelihood of weight matrix in window
sub wmLogScoreSum {
	my $win_ref = shift;
	my $wm_ref = shift;
	my @win = @$win_ref;
	my $score = 0;
	for (my $pos = 0; $pos < @win; $pos++) {
		$score += ${${$wm_ref}[$pos]}[$win[$pos]];
		if ($score < -30) {
			return "FALSE"
		}
	}
	return $score;
}
###

### TO DO
# only keep top hits ($no_hits) -> use hash (key = score, value = id), see notes ???
# add MIN SCORE to parameters
# remove all unnecessary stuff
# Sequences must be in 1 line!

#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	11-JAN-2013
### Modified:	15-JAN-2013
#######

#######
### FUNCTION:
### ---------
### Calculates the absolute distance for each region of a query BED file to the closest region in a SORTED target BED file
#######

#######
### ARGUMENTS:
### ----------
### 1. Query BED file name
### 2. Target BED file name
### 3. Genome sizes file name
### 4. Output file name
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### A tab-separated file containing the IDs (4th column!) of the query BED file and the closest distance to any of the regions in the target BED file
#######

#######
### USAGE:
### ------
### perl /path/to/distance_bed_bed.pl <path/to/query_bed_file> <path/to/target_bed_file> <path/to/genome_sizes_fasta_file> <path/to/output_file>
### Example: perl ./distance_bed_bed.pl ./query.bed ./target.bed ./hg19_sizes.fa ./out
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Main program
### C. Subroutines
#######


### A. Pre-requisites
## Pragmas
use warnings;
use strict;
## Command-line arguments / initialization
my $query_file = shift;
my $target_file = shift;
my $gen_file = shift;
my $out_file = shift;
###


### B. Main program 
# Load query BED file to hash of arrays
my $query_HoHoA_ref = &bedToHashOfHashesOfArrays($query_file, 4, 1, "start", "stop", "rest");
# Load target BED file to hash of arrays
my $target_HoHoA_ref = &bedToHashOfHashesOfArrays($target_file, 4, 1, "start", "stop", "rest");
# Create chromosome sizes lookup hash
my $chr_size_hash_ref = &fastaToHash($gen_file);
# Calculate distances
&bedDistance($query_HoHoA_ref, $target_HoHoA_ref, $chr_size_hash_ref);
# Exit program
exit;
###


### C. Subroutines

sub fastaToHash {
### Pass file in fasta format and load into hash
### Keys: ID lines (minus leading '>')
### Values: Sequences (multiple lines are concatenated, separating newline characters - if present - are removed) 
	# Pass arguments
	my $file = shift;
	# Declare initialize variables
	my %hash;
	my $seq = "";
	my $id;
	# Open fasta file handle
	open FA, "$file";
	## Crawl through fasta file line by line
	while (<FA>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block if identifier line (contains leading '>')
		if ($line =~ s/^>//) {
			# Unless empty string, pass sequence ($seq) to hash (key = $id)
			$hash{$id} = $seq if $seq ne "";
			# Set identifier variable ($id)
			$id = $line;
			# Empty sequence string
			$seq = "";
		}
		## Else concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	# Add last entry to hash
	$hash{$id} = $seq;
	# Close fasta file handle
	close FA;
	# Screen output
	print "File $file processed.\n";
	# Return hash reference
	return \%hash;
}

sub bedToHashOfHashesOfArrays {
### Accepts file in bed format (or other tab-separated file)
### Returns a hash of hashes of arrays containing all column values (addressable by column name) grouped by a common column
### The total number of columns, the number of the grouping column (= key for the outer hash), and the names of the remaining columns (= keys for the inner hashes) are passed when calling (e.g. 3, 1, "start", "stop", "strand")
### Values of the non-grouping columns are used to populate the values of the inner hashes (i.e. arrays)
	## Pass arguments
	my $bed_file = shift;
	my $no = shift;
	my $first = shift;
	## Pass hash key arguments
	my @keys;
	foreach (@_) {
		push @keys, $_;
	}
	# Declare/initialize variables
	my %HoHoA;
	# Open bed file handle
	open BED, $bed_file;
	## Crawl through bed file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_, $no;
		# Extract lookup/key value for outer hash from @line array
		my $value = splice @line, ($first - 1), 1;
		## IF outer hash entry with key $value does not exist, create both outer and inner hash
		unless (exists $HoHoA{$value}) {
			## Add value for each key of inner hash
			foreach (@keys) {
    			$HoHoA{$value}{$_} = [ shift @line ];
			}
		}
		## ELSE push to values (i.e. arrays) of inner hash
		else {
			## For each key of inner hash, push values to existing ones (array!)
			foreach (@keys) {
    			push @{ $HoHoA{$value}{$_} } , shift @line;
			}
		}
	}
	# Close bed file handle
	close BED;
	# Print status message
	print "File $bed_file processed.\n";
	# Return bed file hash of hashes of arrays
	return \%HoHoA;
}

sub bedDistance {
###
###
###
	## Pass arguments
	my %query = %{shift()};
	my %target = %{shift()};
	my %chr_size = %{shift()};
	# Print status message
	print "Calculating distances...\n";
	# Open output file handle
	open OUT, ">$out_file";
	## Crawl through chromosomes (= outer hash of HoHoA) one by one
	foreach my $chr (keys %query) {
		## Repeat the following as long as there are values in HoHoA array "start"
		while (@{$query{$chr}{"start"}}) {
			## Declare/initialize distance variables ($min_distance gets chromosome size as initial value))
			my $dist;
			my $min_dist = $chr_size{$chr};
			## Populate query-related variables by shifting from HoHoA arrays
			my $start_query = shift @{$query{$chr}{"start"}};
			my $stop_query = shift @{$query{$chr}{"stop"}};
			my $mid_query = ($start_query + $stop_query) / 2;
			my $id_query = shift @{$query{$chr}{"rest"}};
			## Declare target-related variables
			my ($start_target, $stop_target, $mid_target, $id_target);
			## Test whether current chromosome (= outer hash of HoHoA) exists in target HoHoA; else skip chromosome and print nothing!
			if (exists $target{$chr}) {
				## Crawl through values in HoHoA arrays (i.e. regions in target bed file)
				for (my $i = 0; $i < @{$target{$chr}{"start"}}; $i++) {	
					## Populate target-related variables by shifting from HoHoA arrays
					$start_target = ${$target{$chr}{"start"}}[$i];
					$stop_target = ${$target{$chr}{"stop"}}[$i];
					$mid_target = ($start_target + $stop_target) / 2;
					# Calculate absolute distance between mid-point of current query and target ranges
					$dist = abs($mid_query - $mid_target);
					# Proceed to next query range if distance increases compared to minimum distance (i.e. last distance)
					last if $dist > $min_dist;
					# Set minimum distance to current distance
					$min_dist = $dist;
				}
				# Print ID of query and minimum distance to output file
				print OUT "$id_query\t$min_dist\n";
			}
		}
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Minimum distances written to file: $out_file\n";
}


sub rmArrayDuplicates {
	# Pass argument
	@array = @{shift()};
	# Map array to hash
	my %hash = map { $_, 1 } @array;
    # Re-load array with hash keys
    @array = keys %hash;
    # Return array reference
    return \@array;
}







































	# Pass argument
	my %hash = %{shift()};
	# Obtain chromosome names
	my @array = keys %hash;
	## Declare variables
	my (@num, @non_num, @num_sorted, @non_num_sorted, @array_sorted);
	my %hash_sorted;
	## Split numeric from non-numeric chromosomes
	foreach my $chr (@array) {
		## Test IF chromosome is numeric
		if ($chr =~ /^chr\d+/) {
			# Push to array @num
			push @num, $chr;
		}
		## ELSE
		else {
			# Push to array @non_num
			push @non_num, $chr;
		}		
	}
	# Sort numerical chromosome names
	@num_sorted = 	map { $_->[0] }						# 4. Return (sorted) capture group [0], i.e. input string / original chromosome name 
					sort { $a->[1] <=> $b->[1] }		# 3. Numerical sort according to capture group [1]
					sort { $a->[2] cmp $b->[2] }		# 2. ASCII-sort according to capture group [2]
					map { [$_, /chr(\d+)(.*)/] }		# 1. Extract digit(s) following 'chr' and rest of name into capture groups [1] and [2] 
					@num;
	# Sort non-numerical chromosome names
	@non_num_sorted = 	map { $_->[0] }								# 4. Return (sorted) capture group [0], i.e. input string / original chromosome name 
						sort { 										# 3. Sort according to capture group [1]:						
							if ($a->[1] eq 'X') { return -1 }		## 3a. 'X' first
							elsif ($b->[1] eq 'X') { return 1 }
							elsif ($a->[1] eq 'Y') { return -1 }	## 3b. 'Y' second
							elsif ($b->[1] eq 'Y') { return 1 }
							elsif ($a->[1] eq 'M') { return 1 }		## 3c. 'M' last
							elsif ($b->[1] eq 'M') { return -1 }
							else { return $a->[1] cmp $b->[1] } 	# 3d. Rest ASCII-sorted
						}
						sort { $a->[2] cmp $b->[2] }				# 2. ASCII-sort according to capture group [2]	
						map { [$_, /chr(.)(.*)/] }					# 1. Extract letter following 'chr' and rest of name into capture groups [1] and [2]
						@non_num;
	# Add sorted non-numerical to sorted numerical array 
	@array_sorted = (@num_sorted, @non_num_sorted); 
	# Build sorted hash by crawling through each element of @array_sorted...
	foreach my $chr (@array_sorted) {
		# ...and assigning the corresponding value of the unsorted hash %hash
		$hash_sorted{$chr} = $hash{$chr};
	}
	
	foreach my $chr (keys %hash_sorted) {
	# ...and assigning the corresponding value of the unsorted hash %hash
		print "$chr\t$hash_sorted{$chr}\n";
	}
	
	# Return reference to sorted hash
	return \%hash_sorted;
	
sub bedToHashOfSortedArraysOfArrays {    ########?!?!?! KILL???
### Accepts file in bed format
### Returns a hash (key = chr) of arrays (positions) of arrays (start, end positions); the outer arrays are sorted by ascending start positions
### Arrays are sorted by size in ascending order
	## Arguments / declarations
	my $bed_file = shift;
	my %HoAoA;
	# Open bed file handle
	open BED, $bed_file;
	## Crawl through bed file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into relevant coordinates by tabs
		my ($chr, $start, $end) = split /\t/;
		# Use current chromosome value as hash key to array of positions (outer array); push pair of start/end positions (inner array) to outer array
		push @{$HoAoA{$chr}}, [$start, $end];
	}
	# Close bed file handle
	close BED;
	## Sort each outer array of %HoAoA by first, then second elements in inner arrays (number sort)
	foreach my $AoA_ref (values %HoAoA) {
		my @sortedAoA = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$AoA_ref};
		$AoA_ref = \@sortedAoA;
	}
	# Print status message
	print "File $bed_file processed.\n";
	# Return bed file hash of arrays of arrays reference
	return \%HoAoA;
}