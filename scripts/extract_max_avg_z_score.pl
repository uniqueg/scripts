# Works on output files of 'bedtools window' ran on two bed files A (primary) and B
# Prints out the average and max scores of each set of B regions that overlap each of the primary regions from file A (e.g. average z-score in window around each region in A)
# Output is in a BED-like tab-separated format: chr/start/end/id of primary region (file A)/average/max score

my ($chr, $start, $end, $id);
my $sum = 0;
my $count = 0;
my $max = 0;

## Traverse over each line in STDIN
while (<>) {
	# Split line by TAB
	my @line = split "\t", $_;
	# If line contains new primary region...
	if ($line[3] ne $id) {
		# Print metrics for previous region to STDOUT
		print "$chr\t$start\t$end\t$id\t" . $sum/$count . "\t$max\n" if defined $id;
		# Reset metrics
		($sum, $count, $max) = 0;
		# ...and define new region
		($chr, $start, $end, $id) = ($line[0], $line[1], $line[2], $line[3]);
	} 
	
	else {
		## Update metrics (sum, count and max score)
		$sum += $line[10];
		$count++;		
		if ($line[10] > $max) {
			$max = $line[10];
		}
	}
}
# Print metrics for last region to STDOUT
print "$chr\t$start\t$end\t$id\t" . $sum/$count . "\t$max\n" if defined $id;
