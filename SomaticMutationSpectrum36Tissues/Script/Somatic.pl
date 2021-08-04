#!/usr/bin/perl -w
use strict;
use warnings;

#chr     pos     ref     alt     context coverage        alt_count       tissue  sample_id       subject_id
#chr1    25451781        G       A       TTGGG   57      8       Adipose_Subcutaneous    03105   129
my %SNP;
my $file=shift;#call SNP fas
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
		my @ref=split /\s+/,$line;
		$ref[0]=~s/chr//g;
		#print $ref[0]."\n";
		$SNP{$ref[0]}{$ref[1]}.="$ref[2] > $ref[3]\t$ref[7]\n";
}
close IN1;



#1       14363   29806   -
my %uniq;
my %count;
my $bed="/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Reference/human/mRNA.bed";#call SNP fas
open(IN, $bed) or die ("can not open $bed\n");
while (my $line=<IN>) {
	chomp $line;
	my @site=split /\s+/,$line;
	foreach my $key2 (keys %{$SNP{$site[0]}}) {
		if ( $site[1]<= $key2 && $key2<= $site[2]) {
			if ($site[3] eq "+") {
				my @mid=split /\n/,$SNP{$site[0]}{$key2};
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				foreach my $SNP(@mid) {
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t$site[0]\t$key2\t$site[1]-$site[2]\t".$SNP."\n";
				}
			}elsif ($site[3] eq "-") {
				my @mid1=split /\n/,$SNP{$site[0]}{$key2};
				$count{"$site[0]\t$key2"}+=1;#discard SNPs exist in multiple gene
				foreach my $SNP(@mid1) {
					my @mid=split /\s+/,$SNP;
					if (uc($mid[0])=~/A/i) {
						$mid[0]="T";
					}elsif (uc($mid[0])=~/T/i) {
						$mid[0]="A";
					}elsif (uc($mid[0])=~/C/i) {
						$mid[0]="G";
					}elsif (uc($mid[0])=~/G/i) {
						$mid[0]="C";
					}
					if (uc($mid[2])=~/A/i) {
						$mid[2]="T";
					}elsif (uc($mid[2])=~/T/i) {
						$mid[2]="A";
					}elsif (uc($mid[2])=~/C/i) {
						$mid[2]="G";
					}elsif (uc($mid[2])=~/G/i) {
						$mid[2]="C";
					}
					$uniq{"$site[0]\t$key2"}.= "$site[3]\t$site[0]\t$key2\t$site[1]-$site[2]\t$mid[0] > $mid[2]\t$mid[3]\n" ;
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