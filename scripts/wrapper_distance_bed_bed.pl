#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 14-JAN-2013/ Modified: 15-JAN-2013
### Description: The program creates and submits jobs for the perl script distance_bed_bed.pl   
### Usage: perl /path/to/wrapper_distance_bed_bed.pl

## Pragmas
use strict;
use warnings;

# Initialize directory variable $dir
my $dir = "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/bed";

# Open directory
opendir DIR, $dir;
# Read fasta files from directory
my @files= map {s/\.[^.]+$//;$_} grep {/\.bed$/} readdir DIR;
# Close directory
closedir DIR;

## INDUCED
# Traverse through @files array
foreach (@files) {
	# Open output file
	open OUT, ">/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/distance_bed_bed_${_}_induced.job";
	## Print shell script for job
	print OUT '#!/bin/bash' . "\n";
	print OUT '#$ -S /bin/bash' . "\n";
	print OUT '#$ -P project_zavolan' . "\n";
	print OUT '#$ -q fs_short' . "\n";
	print OUT '#$ -o /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/' . "\n";
	print OUT '#$ -j y' . "\n";
	print OUT '#$ -cwd' . "\n";
	print OUT "\n";
	print OUT "time perl /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/distance_bed_bed.pl $dir/$_.bed /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/output/Alex_yh2ax_induced/OUTPUT/yH2Ax_ind_FgBg-peakmerge/allpeaks /import/bc2/home/zavolan/kanitz/big/hg19_sizes.fa /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/tab/induced/distance_${_}_induced.tab\n";
	# Close output file
	close OUT;
	# Submit job
	system("qsub", "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/distance_bed_bed_${_}_induced.job");
}

## UNINDUCED
# Traverse through @files array
foreach (@files) {
	# Open output file
	open OUT, ">/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/distance_bed_bed_${_}_uninduced.job";
	## Print shell script for job
	print OUT '#!/bin/bash' . "\n";
	print OUT '#$ -S /bin/bash' . "\n";
	print OUT '#$ -P project_zavolan' . "\n";
	print OUT '#$ -q fs_short' . "\n";
	print OUT '#$ -o /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/' . "\n";
	print OUT '#$ -j y' . "\n";
	print OUT '#$ -cwd' . "\n";
	print OUT "\n";
	print OUT "time perl /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/distance_bed_bed.pl $dir/$_.bed /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/output/Alex_yh2ax/OUTPUT/yH2Ax_FgBg-peakmerge/allpeaks /import/bc2/home/zavolan/kanitz/big/hg19_sizes.fa /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/tab/uninduced/distance_${_}_uninduced.tab\n";
	# Close output file
	close OUT;
	# Submit job
	system("qsub", "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/pat_search_2/jobs/distance_bed_bed/distance_bed_bed_${_}_uninduced.job");
}