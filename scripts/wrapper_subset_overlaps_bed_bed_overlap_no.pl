#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 25-FEB-2013/ Modified: 25-FEB-2013
### Description: The program creates and submits jobs for the R script subset_overlaps_bed_bed_overlap_no.R for the overlap analysis of yH2Ax ChIP peak clusters generated with various window sizes and AsiSI motifs
### Usage: perl /path/to/wrapper_subset_overlaps_bed_bed_overlap_no.pl

## Pragmas
use strict;
use warnings;

# Initialize directory variable $dir
my $dir = '/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/clusters/output/bed';

# Open directory
opendir DIR, $dir;
## Read files from directory
my @files = grep {/^clusters.*bed$/} readdir DIR;
# Close directory
closedir DIR;

# Traverse through @files array
for my $file (@files) {
	# Split filename
	my @filename = split /[_|\.]/, $file;
	# Build job filename
	my $job = "subset_overlaps_bed_bed_overlap_no_" . $filename[0] . "_" . $filename[1] . "_" . $filename[4] . ".job";
	# Build output filename
	my $out = $filename[0] . "_" . $filename[1] . "_with_AsiSI_sites_" . $filename[4] . "." . $filename[5];
	# Open output file
	open OUT, ">/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/jobs/$job";
	## Print shell script for job
	print OUT '#!/bin/bash' . "\n";
	print OUT '#$ -S /bin/bash' . "\n";
	print OUT '#$ -P project_zavolan' . "\n";
	print OUT '#$ -q fs_long' . "\n";
	print OUT '#$ -j y' . "\n";
	print OUT '#$ -cwd' . "\n";
	print OUT '#$ -o /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/jobs' . "\n";
	print OUT "\n";	
	print OUT '## GENERAL INFO' . "\n";
	print OUT '# Author: Alexander Kanitz' . "\n";
	print OUT '# Created: 25-FEB-2013' . "\n";
	print OUT '# Modified: 25-FEB-2013' . "\n";
	print OUT "\n";
	print OUT '## STATUS MESSAGE' . "\n";
	print OUT 'echo "JOB STARTED."' . "\n";
	print OUT "\n";
	print OUT '## EXPORT' . "\n";
	print OUT '# Export R library path' . "\n";
	print OUT 'export R_LIBS_USER=/import/bc2/home/zavolan/kanitz/R/R-2.15.2/library' . "\n";
	print OUT "\n";
	print OUT '## VARIABLES' . "\n";
	print OUT '# Path to Rscript' . "\n";
	print OUT "Rscript='/import/bc2/home/zavolan/kanitz/R/R-2.15.2/bin/Rscript'" . "\n";
	print OUT '# Path to script file' . "\n";
	print OUT "s_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/scripts'" . "\n";
	print OUT '# Path to query files' . "\n";
	print OUT "qu_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/clusters/output/bed'" . "\n";
	print OUT '# Path to subject file' . "\n";
	print OUT "sb_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/input'" . "\n";	
	print OUT '# Genome variable' . "\n";
	print OUT 'gen="hg19"' . "\n";
	print OUT '# Minimum overlap' . "\n";
	print OUT 'min_ol=1' . "\n";
	print OUT '# Minimum score' . "\n";
	print OUT 'min_sc=0' . "\n";
	print OUT '# Path to output files' . "\n";
	print OUT "o_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/output'" . "\n";
	print OUT "\n";
	print OUT '## COMMAND' . "\n";
	print OUT 'time $Rscript $s_dir/subset_overlaps_bed_bed_overlap_no.R $qu_dir/' . $file . ' $sb_dir/AsiSI_orig.bed $gen $min_ol $min_sc $o_dir/' . $out . "\n";
	print OUT "\n";
	print OUT '## STATUS MESSAGE' . "\n";
	print OUT 'echo "JOB COMPLETED."' . "\n";
	# Close output file
	close OUT;
	# Submit job
	system("qsub", "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/jobs/$job");
}