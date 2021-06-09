#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2	29310	T	4	.+1T.+1T,+1t,+1t	BBBB	91,91,46,46	1:S100007198L1C007R0440892555,1:S100007198L1C007R0440892555,2:S100007198L1C007R0440892555,2:S100007198L1C007R0440892555

my $file1=shift;#pileup file
my @file=split /\./,$file1;
my @JC=split /\-/,$file[0];#start 0 and end 1
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	my %mutation=();
	my %UMI_read=();
	my %UMI_read_POS=();

	my %UMI_count2;
	my @SNP=split /\t/,$line;
	  #UMI UMI_read_pair_number
	unless ($SNP[4]=~/[\>|\<]/) {
		  $SNP[4]=~s/\^.//g;
		  $SNP[4]=~s/\$//g;
		  my @read=split /,/,$SNP[7];
		 #$SNP[4]=~ tr/a-z/A-z/;
		  my $count=0;
		while ($SNP[4]=~/[\-|\+][0-9]+[A|T|C|G|a|t|c|g]+/gi) {#|\d
			$mutation{uc($&)}+=1;#SNP count
			my $mid11=uc($&);
			my $mid33=uc($&);
			$mid11=~s/[A|T|C|G|a|t|c|g]+//gi;
			my $mid22=$mid11;
			$mid11=~s/[\-|\+]//gi;
			my $indel_len=$mid11+length($mid22);
			$count+=$indel_len;
			
			my $site;
			$site=(pos $SNP[4])-$count-1;
			my $mid1=$read[($site)];#read
			$UMI_read{uc($mid33)}.="$mid1,";

			my @read_POS=split /,/,$SNP[6];
			my $mid2=$read_POS[($site)];#read
			$UMI_read_POS{uc($mid33)}.="$mid2,";
		}

		#read pair number
		my %hash=();
		foreach my $read_name(@read) {
			if ($read_name=~/:/) {
				my @real=split /:/,$read_name;
				#print $real[1]."\n";
				$hash{$real[1]}+=1;
			}
		}
		my $len=keys %hash;
		my $cover=$#read+1;

		my $UMI_Dis;
		my $left=$SNP[1]-$JC[0];
		my $right=$SNP[1]-$JC[1];
		if (abs($left)<abs($right)) {
			$UMI_Dis=0-abs($left);
		}else{$UMI_Dis=abs($right);}

		foreach my $base(keys %mutation) {
			my $base1=$base;
			if ($base=~/\-/) {
				$base=~s/[\-|\+][0-9]+//gi;
				print "$file[0]\t$file[1]\t$cover\t$len\t$SNP[1]\tDeletion\t$base\t";
				print -length($base)."\t";
			}elsif ($base=~/\+/) {
				$base=~s/[\-|\+][0-9]+//gi;
				print "$file[0]\t$file[1]\t$cover\t$len\t$SNP[1]\tInsertion\t$base\t";
				print length($base)."\t";
			}
			
				print "$UMI_Dis\t$UMI_read_POS{$base1}\t$UMI_read{$base1}\n"; 
		}
	}
	  
}
close IN;
