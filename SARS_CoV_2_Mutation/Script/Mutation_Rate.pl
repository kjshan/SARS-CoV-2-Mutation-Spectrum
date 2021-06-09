#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


#UMI     UMI_read_pair_number    read
#1000-28254      6       920:28354:S100007198L1C001R0060138612,832:975:S100007198L1C004R0520434463,740:977:S100007198L1C008R0460373726,960:28337:S100007198L1C001R0020730201,991:28357:S100007198L1C003R0010708388,920:28354:S100007198L1C008R0600769019,

my $file2="/Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/STAR/Junction_read/SCV2_Unique_UMI.txt";#bed
my %UMI;

my $file1=shift;#pileup file
my @file=split /\./,$file1;
my @JC=split /\-/,$file[0];#start 0 and end 1


#print $file[0]."pp\n";
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	
	if ( $read[0] eq "$file[0]") {#&& $read[1] <= 100
		#print $read[0]."pp\n";
		my @read1=split /,/,$read[2];
		foreach my $readname(@read1) {
			if ($readname=~/:/) {
				my @read2=split /:/,$readname;
				$UMI{"1:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
				$UMI{"2:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
				#print "1:$read2[2]\n";
			}
			
		}
	}
}
close IN;



#NC_045512.2	27548	A	1	^~.	0	1	1:S100007198L1C007R0290862550
my %consensus;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	my %hash=();
	my $len=0;
	my @SNP=split /\t/,$line;
	
	unless ($SNP[4]=~/[\>|\<]/) {
		if ($SNP[3]>1 && $SNP[1] !=4402 && $SNP[1] !=5062 && $SNP[1] !=8782 && $SNP[1] !=28144) {
		  $SNP[4]=~s/\^.//g;
		  $SNP[4]=~s/\$//g;
		    if ($SNP[4] !~/[\.|\,][A|T|C|G]/gi && $SNP[4]!~/[A|T|C|G][\.|,]/gi  && 
				$SNP[4]!~/A[T|C|G]/gi && $SNP[4]!~/T[A|C|G]/gi && 
				$SNP[4]!~/C[T|A|G]/gi && $SNP[4]!~/G[A|C|T]/gi) {
		  
				#if (($SNP[1]<= $JC[0]-15 || $SNP[1]>= $JC[1]+15)) {#|\d
					my @read_name=split /,/,$SNP[$#SNP];
					for (my $i=0;$i<=$#read_name ;$i++) {
						if ($read_name[$i]=~/S100/) {
							my $PCR=$UMI{$read_name[$i]};
							#print $PCR."\n";
							$hash{$PCR}+=1;
						}
						
					}
				$len=keys %hash;
				
				#if ($len>=2) {
					my $left=$SNP[1]-$JC[0];
					my $right=$SNP[1]-$JC[1];
					my $UMI_Dis;
					if (abs($left)<abs($right)) {
						$UMI_Dis=abs($left);
					}else{$UMI_Dis=abs($right);}
					print "$file[0]\t$SNP[1]\t".uc($SNP[2])."\t$UMI_Dis\t$len\n";
					
				#}
				
			}
		}
	}

	  
}
close IN;

