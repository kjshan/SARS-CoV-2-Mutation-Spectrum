#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#I       77      T       1       .       H       S
#1       14415   1       T       T:1:1

my $file3=shift;
my $cutoff=shift;
my %base;
my %Chr_Site;
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	chomp $line;
	my @SNP=split /\s+/,$line;
	if ($SNP[2]>=$cutoff) {
		$Chr_Site{$SNP[0]}{$SNP[1]}="$SNP[3]\t$SNP[2]";
	}
}
close IN1;

#1       14363   29806   -
print "Strand\tChr\tSite\tBase\tcount\n";
my $file1=shift;#"/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Reference/human/Gene.bed";
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	my @site=split /\s+/,$line;
	foreach my $key2 (keys %{$Chr_Site{$site[0]}}) {
		if ( $site[1]<= $key2 && $key2<= $site[2]) {
			my @mid=split /\t/,$Chr_Site{$site[0]}{$key2};
			my $mid=$mid[0];
			if ($site[3] eq "-") {
				if ($mid=~/A/i) {
					$mid="T";
				}elsif ($mid=~/T/i) {
					$mid="A";
				}elsif ($mid=~/C/i) {
					$mid="G";
				}elsif ($mid=~/G/i) {
					$mid="C";
				}
			}
			print "$site[3]\t$site[0]\t$key2\t$mid\t$mid[1]\n";
		}
	}
}
close IN;

