#!/usr/bin/perl

# Alexander Kanitz
# 11-DEC-2013

# Inserts outlinks for expressioin plots for each ENSEMBL ID into the html shared with Keisuke Kaji

use strict;
use warnings;

$/ = "<tr>";

my $head = <>;
print $head;

my $tab_head = <>;
print $tab_head;

while (my $line = <>) {

	my @elements = split /\<\/?td\>/, $line;
	
	my $regex = quotemeta $elements[3];

        my $ENSEMBL = $elements[3];
	$ENSEMBL =~ s/^c\(\\(.*)\\"\)"$/$1/;

	my @ENSEMBL = split /\\", \\"/, $ENSEMBL;

	foreach my $ID (@ENSEMBL) {

		$ID = "<a href='PlotsPopulationsEleni/" . $ID . ".png'>" . $ID . "</a>";

	}

	my $replace = join "<br>", @ENSEMBL;

	$line =~ s/$regex/$replace/;

	print $line;
	
}
