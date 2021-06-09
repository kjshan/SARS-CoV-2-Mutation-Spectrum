#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


#UMI     UMI_read_pair_number    read
#1000-28254      6       920:28354:S100007198L1C001R0060138612,832:975:S100007198L1C004R0520434463,740:977:S100007198L1C008R0460373726,960:28337:S100007198L1C001R0020730201,991:28357:S100007198L1C003R0010708388,920:28354:S100007198L1C008R0600769019,

my $file1="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/JCR_indel.txt2.gt1";#"/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/SCV2_Unique_UMI_gt1.txt";#bed
my %UMI;

open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	#if ($read[1] !~/UMI/ && $read[1] >= 2 && $read[1] <= 20) {#&& $read[1] <= 100
		my @read1=split /,/,$read[2];
		foreach my $readname(@read1) {
			if ($readname=~/:/) {
				my @read2=split /:/,$readname;
				$UMI{"1:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
				$UMI{"2:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
			}
		}
	#}
}
close IN;


my $file2="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/UMI_indel.information";#"/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/UMI2-20.information";#bed
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	chomp $line;
	my @read=split /\t/,$line;
	print "$line\t";
		my %hash=();
		my $Read_count=0;
		my @read_name=split /,/,$read[$#read];
		for (my $i=0;$i<=$#read_name ;$i++) {
			if ($read_name[$i]=~/S100/) {
				my $PCR=$UMI{$read_name[$i]};
				print "$PCR:$read_name[$i],";
				#PCR duplication
				$Read_count+=1;
				$hash{$PCR}+=1;
			}
		}
		
		my $len=keys %hash;
		print "\t$len\t$Read_count\n";
}
close IN;
