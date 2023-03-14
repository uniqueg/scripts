#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 20-FEB-2013/ Modified: 20-FEB-2013
### Description: The program creates and submits jobs for the R script subset_overlaps_bed_bed_avg_score.R for the overlap analysis of the methylomes of two different ESC lines
### Usage: perl /path/to/wrapper_subset_overlaps_bed_bed_avg_score.pl

## Pragmas
use strict;
use warnings;

# Initialize directory variable $dir
my $dir = '/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/output/hg19';

# Open directory
opendir DIR, $dir;
## Read files from directory
my @files_h1 = grep {/^mc_h1/} readdir DIR;
# Close directory
closedir DIR;

# Open directory
opendir DIR, $dir;
## Read files from directory
my @files_i90 = grep {/^mc_i90/} readdir DIR;
# Close directory
closedir DIR;

# Sort file arrays
@files_h1 = sort @files_h1;
@files_i90 = sort @files_i90;

## Extract chromosome names
# Declare variable for array that holds chromosome names
my @chr;
# Traverse through each element of files array @files_h1
foreach my $file (@files_h1) {
	# Split filename of type "mc_h1_?_hg19.bed" by underscore
	my @name = split /_/, $file;
	# Add 3rd element to chromosome array
	push @chr, "chr" . $name[2]; 
}
# Sort chromosome array
@chr = sort @chr;

# Traverse through @files array
for my $i (0..$#chr) {
	# Open output file
	open OUT, ">/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/jobs/overlaps_cell_lines/subset_overlaps_bed_bed_avg_score_$chr[$i].job";
	## Print shell script for job
	print OUT '#!/bin/bash' . "\n";
	print OUT '#$ -S /bin/bash' . "\n";
	print OUT '#$ -P project_zavolan' . "\n";
	print OUT '#$ -q fs_long' . "\n";
	print OUT '#$ -j y' . "\n";
	print OUT '#$ -cwd' . "\n";
	print OUT '#$ -o /import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/jobs/overlaps_cell_lines' . "\n";
	print OUT "\n";	
	print OUT '## GENERAL INFO' . "\n";
	print OUT '# Author: Alexander Kanitz' . "\n";
	print OUT '# Created: 20-FEB-2013' . "\n";
	print OUT '# Modified: 20-FEB-2013' . "\n";
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
	print OUT "s_dir='/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/scripts'" . "\n";
	print OUT '# Path to input files' . "\n";
	print OUT "i_dir='/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/output/hg19'" . "\n";
	print OUT '# Path to output files' . "\n";
	print OUT "o_dir='/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/output'" . "\n";
	print OUT "\n";
	print OUT '## COMMAND' . "\n";
	print OUT 'time $Rscript $s_dir/subset_overlaps_bed_bed_avg_score.R $i_dir/' . $files_h1[$i] . ' $i_dir/' . $files_i90[$i] . ' hg19 1 0 $o_dir/overlap_mc_h1_i90_' . $chr[$i] . '.bed' . "\n";
	print OUT "\n";
	print OUT '## STATUS MESSAGE' . "\n";
	print OUT 'echo "JOB COMPLETED."' . "\n";
	# Close output file
	close OUT;
	# Submit job
	system("qsub", "/import/bc2/home/zavolan/kanitz/DDRNA/pat_search/methylation/Lister_Nature/jobs/overlaps_cell_lines/subset_overlaps_bed_bed_avg_score_$chr[$i].job");
}