use strict;
use warnings;
use Data::Dumper;

my $in='/import/bc2/home/zavolan/stamelaki/test';
my $out='/import/bc2/home/zavolan/stamelaki/outDownstream.txt';

open my $IN, '<', $in or die $!;
open my $OUT,'+>', $out or die $!;

while ( <$IN> ) {
        chomp; # strip off the record separator
        my @line = split;
	my @name = split '_',$line[0];
	#print Dumper(@name);
        #my $row= $name[0]."\t".$name[1]."\t";
#print $line[0], "\t";
	#my $coords;
	my ($start,$end);
        if ($name[2] eq "+") {
                $start=$name[4]+1;
		$end=$name[4]+200;
                #$coords=$start."\t".$end;
#               print $OUT $row;
#               print $row;
        	#my $row=$name[1]."\t".$start."\t".$end."\t".$name[2]."\t".$line[0]."_down200"."\n";
	}
        if ($name[2] eq "-"){
	        $end=$name[3]-1;
		$start=$name[3]-200;

               # $coords=$x."\t".$name[4];
#                print $OUT $row;
#                print $row;
		#my $row=$name[1]."\t".$coords."\t".$name[2]."\t".$line[0]."_down200"."\n";
        }
	my $row=$name[1]."\t".$start."\t".$end."\t".$name[2]."\t".$line[0]."_down200"."\n";
        print $OUT $row;
#       print $row;

}
close ($IN);
close ($OUT);


