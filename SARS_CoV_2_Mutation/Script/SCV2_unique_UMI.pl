#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
#S100007198L1C009R0450431713     163     NC_045512.2     37      255     28M28190N72M    =       28354   28417   CAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACAAACTAAAATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGG FA@F@CDFIACF@F7@1@@IABF@FCAF@FBBBAFIABFBABFABBABCDAF7FA@A7CIEAFFIIABABCFBIFIAABAI7AFFFFCF?BD@7?BAB6I    NH:i:1  HI:i:1  AS:i:192        nM:i:0
my $file2=shift;#"/Dell/Dell9/shankj/Cov/RTmutation/STAR/Junction_read/SCV2.Unique.Perl.JCR.sam";
my %read1;
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	#chomp $line;
	unless ($line=~/^@/) {
		my @read=split /\t/,$line;
		$read1{$read[0]}.="$read[3]:";
	}
}
close IN;
#NC_045512.2     S100007198L1C002R0250205842     66      25380
my $file=shift;#"/Dell/Dell9/shankj/Cov/RTmutation/STAR/Junction_read/SCV2.Unique.JCR.bed";
my %read;
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
		my @read=split /\t/,$line;
		$read{"$read[2]-$read[3]"}.="$read1{$read[1]}$read[1],";
}
close IN;


#38	27676	2
my $file1=shift;#"/Dell/Dell9/shankj/Cov/RTmutation/STAR/Junction_read/SCV2.Unique.JCR.bed.count";#bed
print "UMI\tUMI_read_pair_number\tread\n";
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
		my @read=split /\t/,$line;
		if (exists $read{"$read[0]-$read[1]"}) {
			print "$read[0]-$read[1]\t$read[2]\t".$read{"$read[0]-$read[1]"}."\n";
		}
}
close IN;

