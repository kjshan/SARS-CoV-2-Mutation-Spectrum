#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2     8       2593    G       A:40:0.0154261473197069 C:32:0.0123409178557655 T:28:0.0107983031237948
#NC_045512.2     9       3845    T       A:15:0.00390117035110533        C:2:0.000520156046814044        G:11:0.00286085825747724

my $file1=shift;#Total Mutation
my %mutation;
my %read;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	chomp $line;
	unless ($line=~/Position/) {
		my @read=split /\t/,$line;
		my $mid=$read[2];
			foreach my $readname(@read) {
			if ($readname=~/:/) {
				my @mut=split /:/,$readname;
				$mid-=$mut[1];
				#print "$read[1]\t$mid\n";
			}
		}
		#print "$read[1]\t$mid;\n";
		foreach my $readname(@read) {
			if ($readname=~/:/) {
				my @mut=split /:/,$readname;
				$mutation{$read[1]}{"$read[3] > $mut[0]"}="$mut[1]\t$mid";#pos mut -> read_NO:Total read
			}
		}
	}
}
close IN;


#Position        Ref > Alt       Alt_UMI_number  Alt_read_number Mutation_UMI:Mutation_UMI_Read_W/O_Duplication:Mutation_UMI_Read:Mutation_UMI_Total_Read:Mutation_UMI_same_Read
#8       G > T   1       1       65-28254:1:1:5734950
#8       G > A   2       3       65-28254:2:2:5734950    66-25380:1:1:496750
#9       T > G   2       2       66-25380:1:1:496750     65-28254:1:1:5734950

my $file2=shift;#Junction Mutation

print "Pos\tSNP\tUMI\tUMI_reads\tUMI_Alt_no_PCR_reads\tUMI_Alt_reads\tUMI_ref_reads\tAll_Alt_reads\tAll_Ref_reads\tDis\n";
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	chomp $line;
	unless ($line=~/UMI_reads/) {
		my @read=split /\t/,$line;
		foreach my $UMI(@read) {
			if ($UMI=~/:/) {
				my @mut=split /:/,$UMI;
				my @dis=split /\-/,$mut[0];
				#距离 splice site
				if ($mut[1] !~/UMI/ && $mut[1] >=1) {#至少有2个不同起始位置的reads cover
					my $mutation=$mutation{$read[0]}{$read[1]};#pos mut -> read_NO:Total read
					#if ($#mut<5) {
						#$mut[4]=0;
					#}
					if (((abs($dis[0]-$read[0])+1) > (abs($dis[1]-$read[0])+1)))  {
						print "$read[0]\t$read[1]\t$mut[0]\t$mut[3]\t$mut[1]\t$mut[2]\t$mut[4]\t$mutation\t".(abs($dis[1]-$read[0])+1)."\n";
					}else {print "$read[0]\t$read[1]\t$mut[0]\t$mut[3]\t$mut[1]\t$mut[2]\t$mut[4]\t$mutation\t".(abs($dis[0]-$read[0])+1)."\n";}
				}
			}
		}
	}
}
close IN;

