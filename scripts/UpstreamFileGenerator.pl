use strict;
use warnings;

my $in='/import/bc2/home/zavolan/stamelaki/test';
my $out='/import/bc2/home/zavolan/stamelaki/outUpstream.txt';

open my $IN, '<', $in or die $!;
open my $OUT,'+>', $out or die $!;

while ( <$IN> ) {
        chomp; # strip off the record separator
        my @line = split;
	my @name= split '_', $line[0];
	
#print $line[0], "\t";
	my($start,$end);
	if ($name[2] eq "+") {
		$start=$name[3]-200;
                $end=$name[3]-1;

#		print $OUT $row;
#		print $row;
	}
	if ($name[2] eq "-"){
		$start=$name[4]+1;
		$end=$name[4]+200;
#		 print $OUT $row;
#                print $row;
			
	}
	my $row=$name[1]."\t".$start."\t".$end."\t".$name[2]."\t".$line[0]."_up200"."\n";
        print $OUT $row;

#	print $row;
	
}
close ($IN);
close ($OUT);
 
