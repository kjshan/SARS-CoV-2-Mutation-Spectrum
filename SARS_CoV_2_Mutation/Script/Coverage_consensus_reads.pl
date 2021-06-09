#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

my %count;
my %barcode;
my $file1=shift;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
		chomp $line;
		my @SNP=split /\t/,$line;
		if ($SNP[3]>=15 && $SNP[4]>=2 &&
			$SNP[1] !=4402 && 
			$SNP[1] !=5062 && 
			$SNP[1] !=8782 && 
			$SNP[1] !=28144 ) {
			$count{uc($SNP[2])}+=1;
			#$barcode{uc($SNP[2])}.=;
		}
}
close IN;

foreach my $base(keys %count) {
	print $base."\t".$count{$base}."\n";
}