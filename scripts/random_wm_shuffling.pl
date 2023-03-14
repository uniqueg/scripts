#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 04-JAN-2013/ Modified: 10-JAN-2013
### Description: The program generates random permutations of a given weight matrix (sequence positions separated by an indicated character) while keeping the probabilities for each position (and thus sequence composition). It concludes if an indicated (maximum) number of permutations is reached. A fourth parameter allows the replacement of 0 to a user-defined value. A fifth parameter ('repetition factor'; integer equal to or greater than 0) defines how to deal with multiple occurences of the same permutation. If set to zero, the indicated number of desired weight matrices will surely be reached, but may contain repetitions. If set to an integer > 0, the parameter serves as a factor which - when multiplied with the indicated number of desired variants - sets a limit to the maximum number of random matrices that are generated in order to reach the desired number of DIFFERENT permutations (e.g. desired permutations: 200, repetition factor: 3, maximum generated matrices: 600). When this limit is reached, the program concludes even if the desired number of permutations is not. Thus, if most or all permutations of a given weight matrix are desired, a high number should be chosen for the repetition factor. However, in order to find ALL permutations of a given weight matrix, it is best to infer them using a rational approach, as a random approach a) will be very slow, b) may not capture all permutations even with a very high repetition factor, and c) can only guarantee that all weight matrices are found if the expected number of permutations is known (can then be compared with the 'serial number' of the last found permutation in the output). An output file with the indicated name is produced.    
### Usage: perl /path/to/random_wm_shuffling.pl <weight_matrix/string> <separator/character> <number_of_variants/integer greater than 0> <null replacement/float> <repetition_factor/integer equal to or greater than 0> <path/to/output/file.fa>
### Usage: The weight matrix must be provided in the first line of a text file (all other lines are ignored). Sequence positions must be separated by the indicated separator character and nucleotide probabilities/fractions must use a dot '.' as a decimal separator. 
### Example: perl ./random_wm_shuffling.pl ./wm / 23 0 10 results.fa (writes up to 23, i.e. different permutations of the indicated weight matrix to results.fa; using a high repetition factor will increase the chances that really all permutations are captured; zeros are left as is)
### Example: perl ./random_wm_shuffling.pl ./wm / 100 0.01 0 results.fa (writes 100 permutations of the indicated weight matrix to results.fa; as there are only 24 variations possible, most will constitute repetitions; 0s are set to 0.01) 

### A. Pre-requisites
## Pragmas
use warnings;
use strict;
## Command-line arguments
my $in_file = shift;
my $sep = shift;
my $variants = shift;
my $null_repl = shift;
my $repetitions = shift;
my $out_file = shift;
## Variable declarations/initializations
my $in_wm;
my @out_wm;
my @wm;
my $tmp_wm;
my $mltp_occ = 0;
###

### B. Generate permutations
# Open input file handle
open (IN, "$in_file");
# Read weight matrix from first line of $in_file
$in_wm = <IN>;
# Remove trailing newline character
chomp($in_wm);
# Replace zeros: 1. '0.0' followed by no significant digit 2. '0' at beginning not followed by '.'; 3. '0' not preceded by digit and not preceded or followed by '.'; 4. '0' at end not preceded by digit or '.'
$in_wm =~ s/0\.(0)+(?!\d)|^0(?!\.)|(?<!(\d|\.))0(?!\.)|(?<!(\d|\.))0$/$null_repl/g;
# Close input file handle
close (IN);
# Add input matrix $in_wm to @out_wm results array
push(@out_wm, $in_wm);
# Split matrix $in_wm into @wm array of characters (i.e. nucleotides)
@wm = split(/$sep/, $in_wm);
# Generate shuffled matrices until $variants count is reached and the number of multiple occurences is below the product of the desired number of $variants and the indicated repetition factor $repetitions
while ($#out_wm < $variants && $mltp_occ <= $repetitions * $variants) {
  # Permutes array @wm in place
  &fisher_yates_shuffle( \@wm );
  # Join permutated matrix to string using separator $sep
  $tmp_wm = join($sep,@wm);
  ## Unless permutations are not required to be different from each other, check for existence of permutation in @out_wm array; if found, increase $mltp_occ and next iteration of loop immediately
  unless ($repetitions == 0) {
  	if (grep {$_ eq $tmp_wm} @out_wm) {
  		$mltp_occ++;
  		next;
   	}
  }
  # Else push permutated matrix $tmp_wm to array
  push(@out_wm, $tmp_wm);	
}

### C. Print output
# Open output file handle
open (OUT, ">$out_file");
# Print original matrix $in_wm
print OUT ">original_wm\n$in_wm\n";
# Calculate number of maximum leading zeroes for ID formatting
my $digits = length($#out_wm);
## Print resulting permutations
for (my $i=1; $i <= $#out_wm; ++$i) {
    my $id = sprintf("%0${digits}d", $i);
    print OUT ">perm_wm_$id\n$out_wm[$i]\n";      
}
# Close output file handle
close (OUT);
###

### D. Clean-up
# Print status message
print "Results written to $out_file.\n";
# Exit program
exit;
###

### E. Subroutines

## E1. Generate random permutation
# from Perl Cookbook, Chapter 4 (http://docstore.mik.ua/orelly/perl/cookbook/ch04_18.htm)
# fisher_yates_shuffle( \@array ) : generate a random permutation of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
###