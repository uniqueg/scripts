#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 14-JAN-2013/ Modified: 15-JAN-2013
### Description: The program creates and submits jobs for the perl script pat_search.pl
### Usage: perl /path/to/wrapper_pat_search.pl

## Pragmas
use strict;
use warnings;

# Initialize directory variables
my $dir_gen = "/import/bc2/data/databases/UCSC/hg19";
my $dir_pat = "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2";

## Load chromosome files
# Open directory
opendir CHR, $dir_gen;
# Read fasta files from directory
my @chromosomes= map {s/\.[^.]+$//;$_} grep {/\.fa$/} readdir CHR;
# Close directory
closedir CHR;

## Load pattern files
# Open directory
opendir PAT, $dir_pat;
# Read fasta files from directory
my @patterns= map {s/\.[^.]+$//;$_} grep {/\.fa$/} readdir PAT;
# Close directory
closedir PAT;

## Traverse through @chromsomes array
foreach (@chromosomes) {
	# Initialize variable $chr
	my $chr = $_;
	## Traverse through @patterns array
	foreach (@patterns) {
		# Initialize variable $pat
		my $pat = $_;
		# Open output file
		open OUT, ">pat_search_${chr}_${pat}.job";
		## Print shell script for job
		print OUT '#!/bin/bash' . "\n";
		print OUT '#$ -S /bin/bash' . "\n";
		print OUT '#$ -P project_zavolan' . "\n";
		print OUT '#$ -q fs_long' . "\n";
		print OUT '#$ -j y' . "\n";
		print OUT '#$ -cwd' . "\n";
		print OUT "\n";
		print OUT "time perl pat_search.pl ${dir_pat}/${pat}.fa ${dir_gen}/${chr}.fa 0 AsiSI_perm_${chr}_${pat}.bed\n";
		# Close output file
		close OUT;
		# Submit job
		system("qsub", "pat_search_${chr}_${pat}.job");
	}
}