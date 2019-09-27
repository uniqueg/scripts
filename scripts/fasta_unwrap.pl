#!/usr/bin/perl

## Alexander Kanitz, Biozentrum, University of Basel
## 15-AUG-2014 

use strict;
use warnings;

# Set input field separator to ">"
$/ = ">";

# Discard first occurence of input field separator (empty record)
<>;

# Traverse file record by record
while (<>) {

	# Remove input field separator
	chomp;
	# Extract ID and sequence lines
	my ($id, @seq) = split /\n/;
	# Join sequence lines
	my $seq = join "", @seq;
	# Print record
	print STDOUT ">", $id, "\n", $seq, "\n";

}
