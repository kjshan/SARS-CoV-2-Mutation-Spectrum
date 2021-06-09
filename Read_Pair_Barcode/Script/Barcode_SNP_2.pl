#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#II_414287	C > T	1	0.04	0.96	25

my %Chr_Site;
my $file=shift;#SNP file
open (SNP,$file) or die ("can not open $file\n");
while (my $line1=<SNP>) {
	chomp $line1;
	my @mid=split /\t/,$line1;
	if ($mid[0]!~/20S_RNA_narnavirus/) {
		my @site=split /\_/,$mid[0];
		$Chr_Site{$site[0]}{$site[1]}.="$line1\n";
	}else{
		my @site=split /\_/,$mid[0];
		$Chr_Site{"20S_RNA_narnavirus"}{$site[$#site]}.="$line1\n";
	} 
}
close SNP;



#KE148125.1      2278    2379    +
my %uniq;
my %count;
my $file1=shift;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	my @site=split /\s+/,$line;
	foreach my $key2 (keys %{$Chr_Site{$site[0]}}) {
		if ( $site[1]<= $key2 && $key2<= $site[2]) {
			if ($site[3] eq "+") {
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				my @total_SNP=split /\n/,$Chr_Site{$site[0]}{$key2};
				foreach my $SNP(@total_SNP) {
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t".$SNP."\n";
				}
			}elsif ($site[3] eq "-") {
				my @total_SNP=split /\n/,$Chr_Site{$site[0]}{$key2};
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				foreach my $SNP (@total_SNP) {
					my @mid=split /\s+/,$SNP;
					if (uc($mid[1])=~/A/i) {
						$mid[1]="T";
					}elsif (uc($mid[1])=~/T/i) {
						$mid[1]="A";
					}elsif (uc($mid[1])=~/C/i) {
						$mid[1]="G";
					}elsif (uc($mid[1])=~/G/i) {
						$mid[1]="C";
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
					##print "Chr\tSite\tSNP\tread_number\tread_fraction\treads\tRef_reads_number\tTotal_reads_number\tFragment_number\tFragment_ref\n";
					#II_414287	C > T	1	0.04	0.96
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t$mid[0]\t$mid[1] > $mid[3]\t$mid[4]\t$mid[5]\t$mid[6]\t$mid[7]\n" ;
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