#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#9284-19426      1       9231:19465:S100007198L1C006R0290432564,

my %UMI;
my %read;
my $file=shift;#file name
open(IN, $file) or die ("can not open $file\n");
#28M2D18M28190N54M
while (my $line=<IN>) {
	chomp $line;
	my @read=split /[\t|,]/,$line;
	for (my $i=2;$i<=$#read ;$i++) {
		 $read{$read[$i]}+=1;
	}
}
close IN;

open(IN, $file) or die ("can not open $file\n");
#28M2D18M28190N54M
while (my $line=<IN>) {
	chomp $line;
	my @read=split /[\t|,]/,$line;
	my $count=0;
	my $reads=();
	
	for (my $i=2;$i<=$#read ;$i++) {
		 if ($read{$read[$i]}==1) {
			 $reads.=$read[$i].",";
			 $count+=1;
		 }
	}
	if ($count>=2) {
		print $read[0]."\t".$count."\t$reads";
		print "\n";
	}

}
close IN;
