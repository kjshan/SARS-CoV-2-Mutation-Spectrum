#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


my $file=shift;#total sequence
my $ref=shift;
$/="\n>";
my $num=1;
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
	my @seq=split /\n/,$line;
	$seq[0]=~s/\>//;
	$num+=1;
	open (OUT,">$num")||die;
		print OUT ">$seq[0]\n$seq[1]\n";
	close OUT;
	my $newname=$num.".fas";
	`cat $ref ./$num > $newname `;
	`rm ./$num `;
}
close IN1;
