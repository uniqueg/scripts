use strict;
use warnings;

my %hash;

while ( <> ) {
        chomp; # strip off the record separator
        my @line = split;
        my $key = join "\t", @line[0..4];

        if (defined $hash{$key}) {

                $hash{$key}= $hash{$key} . ':' . $line[5];
        } else {
                $hash{$key}= $line[5];
        }
}

foreach my $key (sort keys %hash) {

        my @coor = split "\t", $key;
        my $ex_id = join '_', $coor[4], $coor[0], $coor[3], $coor[1], $coor[2];
        print $ex_id . "\t" . $hash{$key} . "\n";
}

