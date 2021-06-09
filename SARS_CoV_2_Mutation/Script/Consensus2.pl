#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
#use POSIX;
#AverBQ  Barcode Pos     Dis     non_PCR Match
#29.5    1002-28709      28724   15      2       Match

my %BQ_mismatch;
my %BQ_match;
my $file1=shift;#pileup file
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	unless ($line=~/Barcode/) {
		chomp $line;
		my @SNP=split /\t/,$line;
		if ($SNP[3]>=15 && $SNP[4]>=2  &&
			$SNP[2] !=4402 && 
			$SNP[2] !=5062 && 
			$SNP[2] !=8782 && 
			$SNP[2] !=28144 ) {
			
			if ($SNP[5] eq "Match") {
				$BQ_match{int($SNP[0])}+=1;
			}else {$BQ_mismatch{int($SNP[0])}+=1;}
		}
	}
}
close IN;

foreach my $BQ(keys %BQ_match) {
	print $BQ."\t".$BQ_mismatch{$BQ}/($BQ_mismatch{$BQ}+$BQ_match{$BQ})."\n";
}