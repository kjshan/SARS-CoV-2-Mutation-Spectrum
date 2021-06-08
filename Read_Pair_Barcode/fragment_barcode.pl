#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#Sam file
#S100007198L1C004R0320729945     163     NC_045512.2     1       255     22S75M3S        =       1       28289   CTTATATCATTTGTATTTCCACATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACAAA    FAAAA@AFAA77BA@@A@FF?

#V300081983L2C001R0060930119     99      XII     453360  1       150M    =       453405  195     GAAAACTATTCCTTCCTGTGGATTTTCACGGGCCGTCACAAGCGCACCGGAGCCAGCAAAGGTGCTGGCCTCTTCCAGCCATAAGACCCCATCTCCGGATAAACCAATTCCGGGGTGATAAGCTGTTAAGAAGAAAAGATAACTCCTCCC  EEFFFFFFFFFFFFGFFGDFFFFFFFFFFFFFFGFEFFFFFFFFFFFGFF6FGFFFFFFFGFBGFFFFFFFFFFGFFFFGFFGGFDFEFFGFGFFFFGFDGEFGGGFFFFGFFFF6GDFDDFFGGFFFDFEGFFGF6GEFFDEGGGGGFG  NH:i:3  HI:i:1  AS:i:298        nM:i:0

my $file=shift;#file name
my %pair;
my %pair_read;
my %barcode;

open(IN, $file) or die ("can not open $file\n");

while (my $line=<IN>) {
	if ($line=~/^@/) {
		print $line;
	}else{
		my @read=split /\t/,$line;
		if ($read[5] !~/S|H|P|I|D/i && $read[11] eq "NH:i:1") {
			$pair{$read[0]}+=1;
			$pair_read{$read[0]}.=$line;
			$barcode{$read[0]}.="$read[2]:$read[3]:";
		}
	}
}
close IN;


foreach my $read(keys %pair) {
	if ($pair{$read}==2) {
		my @mid=split /\n/,$pair_read{$read};
		print $barcode{$read}.$mid[0]."\n";
		print $barcode{$read}.$mid[1]."\n";
	}
}