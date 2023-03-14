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
my $wm_file = shift;
my $wm_len = 18;#shift;
my $fa_file = shift;
my $min_score = 18;#shift;
my $no_hits = 1000;#shift;
my $pos_sep = "/";#shift;
my $nt_sep = ",";#shift;
my $out_file = "out_test_2";#shift;
# Variable declaration
my %wm;
my $start;
###

### B. Read weight matrices
# Open weight matrix file handle
open WM, $wm_file;
# Declare weight matrix hash %wm

# Traverse through weight matrices one by one
while (<WM>) {
	# Remove trailing new line characters from input lines
	chomp;
	# Assign read line to identifier variable $id_wm
	my $id_wm = $_;
	# Read line containing nucleotide probabilities
	my $wm_line = <WM>;
	# Remove trailing new line characters from input lines
	chomp $wm_line;
	# Load into hash (keys: ids, values: probabilities)
	$wm{$id_wm} = $wm_line;
}
# Close weight matrix file handle
close WM;	
# Load 2-dimensional array (sequence position, nucleotide probability)
for my $key (keys %wm) {
	$wm{$key} = &loadAoAfromString($wm{$key}, $pos_sep, $nt_sep);
	# Compute complement of weight matrix
	$wm{"${key}_comp"} = &wmComplement($wm{$key})
}

### C. Read fasta file
# Open output file handle
open OUT, ">$out_file";
# Open fasta file handle
open FA, $fa_file;
# Declare sequence identifier @id_seq
my $id_seq;
# Declare sequence array @seq
my @seq;
# Crawl through genome fasta file line by line
while (<FA>) {
	# Remove trailing new line character from input line
	chomp;
	## If line starts with ">" (identifier line), calculate score for previous sequence, extract identifier $id_seq and empty sequence array @seq; else, transform all letters to uppercase, push to array and replace with numbers (see below) 
	if (substr($_, 0, 1) eq ">") {
		&calcWinScore if @seq;
		$id_seq = substr $_, 1;
		@seq = ();
	}
	else {
		# Transform lowercase to uppercase and push to array
		push @seq, split(//, uc($_));
		## Index sequence (i.e. replace letters with numbers corresponding to wm positons; no requirement for conditional statements in sliding window loop)
		foreach (@seq) {
			if		($_ eq "A")	{$_ = 0}
			elsif	($_ eq "C")	{$_ = 1}
			elsif	($_ eq "G")	{$_ = 2}
			elsif	($_ eq "T")	{$_ = 3}
			else				{$_ = 4}
		}
	}
}
# Calculate scores for last sequence
&calcWinScore;
# Close fasta file handle
close FA;
# Close output file handle
close OUT;

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
		# Add value for other letters than ACGT (e.g. N; arbitrary value: ln of 0.25)
		push @logs, -1.386294361;						
		push @AoA, [ @logs ];
	}
	
	return \@AoA;
}
##

## E2. Compute complement of weight matrix		################# DO REVERSE!!!
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

### E3. Slide window and compute log likelihoods of weight matrices in window
sub calcWinScore {
	for ($start = 0; $start <= (@seq - $wm_len); $start++) {
		my @win = @seq[$start..($start + $wm_len - 1)];
		for my $key (keys %wm) {
			my $score = 0;
			for (my $pos = 0; $pos < $wm_len; $pos++) {
				$score += ${${$wm{$key}}[$pos]}[$win[$pos]];
			}
			if ($score >= $min_score) {
				print OUT $key . "\t" . $id_seq . "\t" . ($start + 1) . "\t" . ($start + $wm_len) . "\t" . $score . "\n";
			}
		}
	}
}

### TO DO
# only keep top hits ($no_hits) -> use hash (key = score, value = id), see notes ???
# REVERSE complement!
# remove all unnecessary stuff