#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#35.5	XI_190201	1	Cover
my $file1=shift;#bed
my %mismatch;
my %Cover;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @mid=split /\t/,$line;
	if ($line=~/Mis/) {
		$mismatch{int($mid[0])}+=$mid[2];
	}elsif ($line=~/Cover/) {
		$Cover{int($mid[0])}+=$mid[2];
	}
	
}
close IN;

foreach my $BQ(sort keys %Cover) {
	if (exists $Cover{$BQ}) {
		print "$BQ\t$mismatch{$BQ}\t$Cover{$BQ}\t".($mismatch{$BQ}/$Cover{$BQ})."\n";
	}else {
		print "$BQ\t0\t$Cover{$BQ}\t0\n";
	}
}