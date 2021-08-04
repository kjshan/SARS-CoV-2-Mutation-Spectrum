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
	if ($#SNP>=5) {
		for (my $i=4; $i<=$#SNP ;$i++) {
			if ($SNP[$i]=~/:/) {
				my @num=split /:/,$SNP[$i];
				unless (uc($SNP[3]) eq uc($num[0])) {
					 if ($SNP[3]!~/N/i && $num[0]!~/N/i  ) {
						 $Chr_Site{$SNP[0]}{$SNP[1]}.="$SNP[0]\t$SNP[1]\t$SNP[2]\t$SNP[3] > ".uc($num[0])."\t$num[1]\t$num[2]\n";
					 }
				}else{ 
					if ($SNP[3]!~/N/i && $num[0]!~/N/i ) {
						$Ref{$SNP[0]}{$SNP[1]}="$num[1]";
					}
				}
			}
		}
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
		if ( $site[1]<= $key2 && $key2<= $site[2] && exists $Ref{$site[0]}{$key2}) {
			if ($site[3] eq "+") {
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				my @total_SNP=split /\n/,$Chr_Site{$site[0]}{$key2};
				foreach my $SNP(@total_SNP) {
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t".$SNP."\t".$Ref{$site[0]}{$key2}."\n";
				}
			}elsif ($site[3] eq "-") {
				my @total_SNP=split /\n/,$Chr_Site{$site[0]}{$key2};
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				#$SNP[0]\t$SNP[1]\t$SNP[2]\t$SNP[3] > ".uc($num[0])."\t$num[1]\t$num[2]
				foreach my $SNP (@total_SNP) {
					my @mid=split /\s+/,$SNP;
					if (uc($mid[5])=~/A/i) {
						$mid[5]="T";
					}elsif (uc($mid[5])=~/T/i) {
						$mid[5]="A";
					}elsif (uc($mid[5])=~/C/i) {
						$mid[5]="G";
					}elsif (uc($mid[5])=~/G/i) {
						$mid[5]="C";
					}
					if (uc($mid[3])=~/A/i) {
						$mid[3]="T";
					}elsif (uc($mid[3])=~/T/i) {
						$mid[3]="A";
					}elsif (uc($mid[3])=~/C/i) {
						$mid[3]="G";
					}elsif (uc($mid[3])=~/G/i) {
						$mid[3]="C";
					}
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t$mid[0]\t$mid[1]\t$mid[2]\t$mid[3] > $mid[5]\t$mid[6]\t$mid[7]\t".$Ref{$site[0]}{$key2}."\n" ;
				}
			}
		}
	}
}
close IN;


foreach my $site(sort keys %count) {
	if ($count{$site}==1) {
		print $uniq{$site};
	}
	
}