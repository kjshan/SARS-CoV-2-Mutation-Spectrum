#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


#II-266921-II-266970-V300081983L2C004R0331321697 99      II      266921  255     150M    =       266970  199     GTGTGTTCGTGCCTATCCAATTCTACAACGTAGTAAGATGGAGAAAACGTGTACGAATTAGTCAACAAAGCTTTCAGTTTTTTGCTATATGGGCAGCCAGTCTTGCTAAATACAATCATGGGTGA
my $file2=shift;
my $header;
my %barcode;
my %NUM;
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	if ($line=~/^@/) {
		$header.=$line;
	}else{
		my @read=split /\:/,$line;
		$NUM{"$read[0]:$read[1]:$read[2]:$read[3]"}+=1;
		$barcode{"$read[0]:$read[1]:$read[2]:$read[3]"}.=$NUM{"$read[0]:$read[1]:$read[2]:$read[3]"}."-".$line;
		
	}
}
close IN;

foreach my $barcode(keys %barcode) {
	print "$barcode.$NUM{$barcode}.sam"."\n";
	open (OUT,">$barcode.$NUM{$barcode}.sam")||die;
		print OUT $header;
		print OUT $barcode{$barcode};
	close OUT;
}