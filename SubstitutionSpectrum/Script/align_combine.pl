#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


my $file=shift;#total sequence
$/="\n>";
my $ref=shift;
my %gap;
my $align;
my %num;
my @file1=split /\//,$file;
my @file2=split /\.fas\.out/,$file1[-1];
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
	my @seq=split /\n/,$line;
	$seq[0]=~s/\>//;
	my $align;
	for (my $i=1;$i<=$#seq ;$i++) {
		$align.=$seq[$i];
	}
	my @base=split //,$align;
	if ($seq[0] =~/$ref/) {
		for (my $j=0;$j<=$#base ;$j++) {
			if ($base[$j] eq "-") {
				$gap{$j}=1;
			}
		}
	}else {
			print ">$file2[0]\n";
			for (my $j=0;$j<=$#base ;$j++) {
				unless (exists $gap{$j}) {
					print $base[$j];
					
				}
			}
			print "\n";
	}
}
close IN1;

