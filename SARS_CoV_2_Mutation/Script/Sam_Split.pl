#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#Sam file
my $header;
my %read_pair;
my %count=();
my $file1=shift;#"/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/UMI_2_20.sort.sam";
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	if ($line=~/^\@/) {
		$header.=$line;
	}else {
		my @read=split /\t/,$line;
		$count{$read[0]}+=1;
		$read_pair{$read[0]}.="$count{$read[0]}:$line";
	}
}
close IN;


#UMI     UMI_read_pair_number    read
#1000-28254      6       920:28354:S100007198L1C001R0060138612,832:975:S100007198L1C004R0520434463,740:977:S100007198L1C008R0460373726,960:28337:S100007198L1C001R0020730201,991:28357:S100007198L1C003R0010708388,920:28354:S100007198L1C008R0600769019,

my $file2=shift;#"/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/SCV2_Unique_UMI_gt1.txt";
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	#if ($read[1] >= 2 && $read[1] <= 20) {
		open (OUT,">".$read[0].".$read[1].sam")||die;
			print OUT $header;
			my @read1=split /,/,$read[2];
			foreach my $readname(@read1) {
			if ($readname=~/:/) {
				my @read2=split /:/,$readname;
				print OUT $read_pair{$read2[2]};
			}
		}
		close OUT;
		#`samtools mpileup --max-depth 0 --reference /Dell/Dell13/shankj/projects/Cov/NC_045512.2.fna  --output-QNAME -Q 0 -B --adjust-MQ 0 --output-BP /Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/UMI_2_20.sort.bam  -o  /Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/UMI_2_20.sort.mpileup`;
	#}
}
close IN;
#parallel  -j 50 "samtools sort  {} -o ./{}.bam" :::: list.txt
#ls *.bam > list.txt
#parallel  -j 50 "samtools mpileup  --reference /Dell/Dell13/shankj/projects/Cov/NC_045512.2.fna  --output-QNAME -Q 0 -B --adjust-MQ 0 --output-BP  {} -o  ./{}.mpileup" :::: list.txt

#parallel  -j 50 "perl /Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/Information.pl {} " :::: list.txt > ../UMI2-20.information

