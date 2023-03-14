#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 12-DEC-2012 / Modified: 18-DEC-2012
### Description: Reduce number of sequences in "mappingtests.fa" (created by Andreas Gruber) used for mapper benchmarking. A given maximum number of sequences for each category is kept.
### Usage: perl /path/to/compare_tsv.pl </path/to/reference_file/mappingtests.tab> <path/to/fasta/file/mappingtests.fa> <max_number_per_category>

### A. Pre-requisites
use warnings;
use strict;
my $usage = 'Usage: perl /path/to/compare_tsv.pl </path/to/reference_file> <path/to/fasta/file> <max_number_per_category>\n';
my $ref_file = shift or die "$usage";
my $fa_file = shift or die "$usage";
my $max_cat = shift or die "$usage";
###

### B. Reference file processing
#!# Open file, split line, obtain category by concatenating length and type
#!# Generate hash (category => running number)
#!# Print sequence identifier if running number does not exceed treshold indicated in command line ("$max_cat")
my %cat;
my @ids;
open (REF, $ref_file);
while (<REF>) {
	chomp;
	my @cols = split("\t");
	my $id = shift(@cols);
	my $cat = $cols[4] . ':' . length($cols[5]);
	# Checks if hash entry exists for this category; if so, running number is increased, else it is set to 1
	if (exists $cat{$cat}) {
		$cat{$cat}++;	
	}
	else {
		$cat{$cat} = 1;
	}
	# Checks whether running number is below or equal to treshold; if so, sequence id with leading ">" is pushed to array '@ids', else the current line is discarded 
	push (@ids, ">$id") if ($cat{$cat} <= $max_cat);
}
print "Processed file $ref_file.\n";
close (REF);
###

### C. Fasta file processing
#!# 
my %hash = map { $_ => 1 } @ids;
open (FA, $fa_file);
open (OUT, ">$fa_file.$max_cat.out");
while (<FA>) {
	chomp;
	my $id = $_;
	my $seq = <FA>;
	print OUT "$id\n$seq" if (exists $hash{$id});
}
print "Processed file $fa_file.\n";
close (OUT);
close (FA);
###

### D. Clean-up
exit;
###