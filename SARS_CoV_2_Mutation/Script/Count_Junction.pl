#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2     S100007198L1C007R0570111342     65      28254
my $file=shift;#file name
my %JCR;
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
		my @read=split /\t/,$line;
		$JCR{"$read[2]\t$read[3]"}+=1;
}
close IN;

foreach my $key(sort keys %JCR) {
	print "$key\t$JCR{$key}\n";
}