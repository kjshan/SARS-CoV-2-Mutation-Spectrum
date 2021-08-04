#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


my $file=shift;#total sequence
my $ref=shift;
$/="\n>";
my %seq;
my @num=split /,/,$ref;
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
	my @seq=split /\n/,$line;
	$seq[0]=~s/\>//;
	$seq{$seq[0]}=$seq[1];
}
close IN1;

foreach my $name(@num) {
	print ">$name\n$seq{$name}\n";
}