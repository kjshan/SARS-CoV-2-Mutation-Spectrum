#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

my %BQ_mismatch;
my %BQ_match;
my $file1=shift;#pileup file
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	unless ($line=~/Barcode/) {
		chomp $line;
		my @SNP=split /\t/,$line;
		if ($SNP[3]>=15 &&
			$SNP[2] !=4402 && 
			$SNP[2] !=5062 && 
			$SNP[2] !=8782 && 
			$SNP[2] !=28144 ) {
			$BQ_mismatch{$SNP[0]}+=$SNP[5];
			$BQ_match{$SNP[0]}+=$SNP[4];
		}
	}
}
close IN;

foreach my $BQ(keys %BQ_match) {
	print $BQ."\t".$BQ_mismatch{$BQ}/$BQ_match{$BQ}."\n";
}