#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#L-A     1682    T       144

my %Chr_Site;
my $file=shift;#SNP file
open (SNP,$file) or die ("can not open $file\n");
while (my $line1=<SNP>) {
	chomp $line1;
	my @site=split /\t/,$line1;
	$Chr_Site{$site[0]}{$site[1]}.="$line1\n";
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
				my $SNP=$Chr_Site{$site[0]}{$key2};
				$uniq{"$site[0]\t$key2"}= "$site[3]\t".$SNP."\n";
			}elsif ($site[3] eq "-") {
				my $SNP=$Chr_Site{$site[0]}{$key2};
				my @mid=split /\s+/,$SNP;
				if (uc($mid[2])=~/A/i) {
					$mid[2]="T";
				}elsif (uc($mid[2])=~/T/i) {
					$mid[2]="A";
				}elsif (uc($mid[2])=~/C/i) {
					$mid[2]="G";
				}elsif (uc($mid[2])=~/G/i) {
					$mid[2]="C";
				}
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				##print "Chr\tSite\tSNP\tread_number\tread_fraction\treads\tRef_reads_number\tTotal_reads_number\tFragment_number\tFragment_ref\n";
				#II_414287	C > T	1	0.04	0.96
				$uniq{"$site[0]\t$key2"}.= "$site[3]\t$mid[0]\t$mid[1]\t$mid[2]\t$mid[3]\n" ;
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