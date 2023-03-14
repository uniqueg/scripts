#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 03-JAN-2013/ Modified: 10-JAN-2013
### Description: The program generates random permutations of a given sequence while keeping sequence composition. It concludes if an indicated (maximum) number of permutations is reached. A third parameter ('repetition factor'; integer equal to or greater than 0) defines how to deal with multiple occurences of the same permutation. If set to zero, the indicated number of desired sequences will surely be reached, but may contain repetitions. If set to an integer > 0, the parameter serves as a factor which - when multiplied with the indicated number of desired variants - sets a limit to the maximum number of random sequences that are generated in order to reach the desired number of DIFFERENT permutations (e.g. desired permutations: 200, repetition factor: 3, maximum generated sequences: 600). When this limit is reached, the program concludes even if the desired number of permutations is not. Thus, if most or all permutations of a given sequence are desired, a high number should be chosen for the repetition factor. However, in order to find ALL permutations of a given sequence, it is best to infer them using a rational approach, as a random approach a) will be very slow, b) may not capture all permutations even with a very high repetition factor, and c) can only guarantee that all sequences are found if the expected number of permutations is known (can then be compared with the 'serial number' of the last found permutation in the output). An output file with the indicated name is produced.    
### Usage: perl /path/to/random_sequence_shuffling.pl <path/to/seq_file> <number_of_variants/integer greater than 0> <repetition_factor/integer equal to or greater than 0> <path/to/output/file.fa>
### Usage: The weight matrix must be provided in the first line of a text file (all other lines are ignored).
### Example: perl ./random_sequence_shuffling.pl ./seq 23 10 results.fa (writes up to 23, i.e. ALL different permutations of the indicated sequence to results.fa; using a high repetition factor will increase the chances that really all permutations are captured)
### Example: perl ./random_sequence_shuffling.pl ./seq 100 0 results.fa (writes 100 permutations of the indicated sequence to results.fa; as there are only 24 variations possible, most will constitute repetitions)

### A. Pre-requisites
## Pragmas
use warnings;
use strict;
## Arguments
my $in_file = shift;
my $variants = shift;
my $repetitions = shift;
my $out_file = shift;
## Variable declarations/initializations
my @out_seq;
my @seq;
my $tmp_seq;
my $mltp_occ = 0;
###

### B. Generate permutations
# Open input file handle
open (IN, "$in_file");
# Read weight matrix from first line of $in_file
my $in_seq = <IN>;
# Remove trailing newline character
chomp($in_seq);
# Close input file handle
close (IN);
# Add input sequence $in_seq to @out_seq results array
push(@out_seq, $in_seq);
# Split sequence $in_seq into @seq array of characters (i.e. nucleotides)
@seq = split(//, $in_seq);
# Generate shuffled sequences until $variants count is reached and the number of multiple occurences is below the product of the desired number of $variants and the indicated repetition factor $repetitions
while ($#out_seq < $variants && $mltp_occ <= $repetitions * $variants) {
  # Permutes array @seq in place
  &fisher_yates_shuffle( \@seq );
  # Join permutated sequence to string
  $tmp_seq = join("",@seq);
  ## Unless permutations are not required to be different from each other, check for existence of permutation in @out_seq array; if found, increase $mltp_occ and next iteration of loop immediately
  unless ($repetitions == 0) {
  	if (grep {$_ eq $tmp_seq} @out_seq) {
  		$mltp_occ++;
  		next;
   	}
  }
  # Else push permutated sequence $tmp_seq to array
  push(@out_seq, $tmp_seq);	
}

### C. Print output
# Open output file handle
open (OUT, ">$out_file");
# Print original sequence $in_seq
print OUT ">original_seq\n$in_seq\n";
# Calculate number of maximum leading zeroes for ID formatting
my $digits = length($#out_seq);
## Print resulting permutations
for (my $i=1; $i <= $#out_seq; ++$i) {
    my $id = sprintf("%0${digits}d", $i);
    print OUT ">perm_seq_$id\n$out_seq[$i]\n";      
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