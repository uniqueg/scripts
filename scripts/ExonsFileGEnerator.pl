use strict;
use warnings;
use Data::Dumper;

my $in='/import/bc2/home/zavolan/stamelaki/MotEvoFiles/OUT12_B.txt';
my $out='/import/bc2/home/zavolan/stamelaki/MotEvoFiles/outExons.txt';

open my $IN, '<', $in or die $!;
open my $OUT,'+>', $out or die $!;

while ( <$IN> ) {
        chomp; # strip off the record separator
        my @line = split;
        my @name = split '_',$line[0];
	my $row=$name[1]."\t".$name[3]."\t".$name[4]."\t".$name[2]."\t".$line[0]."\n";
        print $OUT $row;
}
close ($IN);
close ($OUT);

