#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	04-JAN-2013
### Modified:	08-JAN-2013
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
my $no_hits = shift;
my $pos_sep = shift;
my $nt_sep = shift;
my $out_file = shift;
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

## E2. Compute complement of weight matrix
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