#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 08-FEB-2013/ Modified: 08-FEB-2013
### Description: The program creates and submits jobs for the C script generate_extended_clusters   
### Usage: perl ./wrapper_generate_extended_clusters.pl

## Pragmas
use strict;
use warnings;

# Initialize directory variable $dir
my $dir = "/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/merge_regions/jobs/generate_extended_clusters";

## Traverse through window sizes 1000 to 1,000,000 (step size = 1000)
for (my $size = 1000; $size <= 3000000; $size = $size + 1000) {
        # Open output file
        open OUT, ">${dir}/generate_extended_clusters_${size}.job";
        ## Print shell script for job
        print OUT '#!/bin/bash' . "\n";
        print OUT '#$ -S /bin/bash' . "\n";
        print OUT '#$ -P project_zavolan' . "\n";
        print OUT '#$ -q fs_short' . "\n";
        print OUT '#$ -j y' . "\n";
        print OUT '#$ -cwd' . "\n";
        print OUT '#$ -o /import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/merge_regions/jobs/generate_extended_clusters' . "\n";
        print OUT "\n";
        print OUT '## GENERAL INFO' . "\n";
        print OUT '# Author: Alexander Kanitz' . "\n";
        print OUT '# Created: 08-FEB-2013' . "\n";
        print OUT '# Modified: 08-FEB-2013' . "\n";
        print OUT "\n";
        print OUT '## STATUS MESSAGE' . "\n";
        print OUT 'echo "JOB STARTED."' . "\n";
        print OUT "\n";
        print OUT '## VARIABLES' . "\n";
        print OUT '# Path to script file' . "\n";
        print OUT "s_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/merge_regions/scripts'" . "\n";
        print OUT '# Path to input files' . "\n";
        print OUT "i_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/merge_regions/input'" . "\n";
        print OUT '# Window size' . "\n";
        print OUT "win=${size}" . "\n";
        print OUT '# Path to output file' . "\n";
        print OUT "o_dir='/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/merge_regions/output/mz'" . "\n";
        print OUT "\n";
        print OUT '## COMMANDS' . "\n";
        print OUT 'time $s_dir/generate_extended_clusters $i_dir/z_cutoff_induced.mz $win > $o_dir/clusters_${win}_z_cutoff_induced.mz' . "\n";
        print OUT 'time $s_dir/generate_extended_clusters $i_dir/z_cutoff_uninduced.mz $win > $o_dir/clusters_${win}_z_cutoff_uninduced.mz' . "\n";
        print OUT "\n";
        print OUT '## STATUS MESSAGE' . "\n";
        print OUT 'echo "JOB COMPLETED."' . "\n";
        # Close output file
        close OUT;
        # Submit job
        system("qsub", "${dir}/generate_extended_clusters_${size}.job");
}