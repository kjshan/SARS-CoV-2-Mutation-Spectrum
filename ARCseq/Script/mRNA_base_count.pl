#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

my %Ref;
my %Chr_Site;
#L-BC-La 236     7       C       C:6:0.857142857142857   A:1:0.142857142857143
my $file3=shift;

open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	chomp $line;
	my @SNP=split /\s+/,$line;
	if ($SNP[2]>=100) {
		$Chr_Site{$SNP[0]}{$SNP[1]}="$SNP[3]\t$SNP[2]";
	}
}
close IN1;

#1       14363   29806   -
my %uniq;
my %count;
my $file1=shift;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	my @site=split /\s+/,$line;
	foreach my $key2 (keys %{$Chr_Site{$site[0]}}) {
		if ( $site[1]<= $key2 && $key2<= $site[2] ) {
			if ($site[3] eq "+") {
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				$uniq{"$site[0]\t$key2"}= $Chr_Site{$site[0]}{$key2};
			}elsif ($site[3] eq "-") {
				my $SNP=$Chr_Site{$site[0]}{$key2};
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				#$SNP[0]\t$SNP[1]\t$SNP[2]\t$SNP[3] > ".uc($num[0])."\t$num[1]\t$num[2]
					my @mid=split /\t/,$SNP;
					if (uc($mid[0])=~/A/i) {
						$mid[0]="T";
					}elsif (uc($mid[0])=~/T/i) {
						$mid[0]="A";
					}elsif (uc($mid[0])=~/C/i) {
						$mid[0]="G";
					}elsif (uc($mid[0])=~/G/i) {
						$mid[0]="C";
					}
					$uniq{"$site[0]\t$key2"}= "$mid[0]\t$mid[1]";
			}
		}
	}
}
close IN;

my %Base_count;
foreach my $site(sort keys %count) {
	if ($count{$site}==1) {
		my @mid1=split /\t/,$site;
		my $Base=$uniq{$site};
		my @mid=split /\t/,$Base;
		my $key=$mid1[0]."\t".uc($mid[0]);
		$Base_count{$key}+=$mid[1];
	}
}

foreach my $base(keys %Base_count) {
	print $base."\t".$Base_count{$base}."\n";
}