#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#Chr     SNP     Ctrl    Stress
#mRNA    A > T   5       5

my %Ctrl;
my %Stress;
my $file1=shift;#Somatic
open(IN1, $file1) or die ("can not open $file1\n");
my %SNP;
while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	$Ctrl{$SNP[0]}{$SNP[1]}=$SNP[2];
	$Stress{$SNP[0]}{$SNP[1]}=$SNP[3];
}
close IN1;


#Chr     SNP     SNP_number      Base_number     Mutation_Freq
#L-BC-La G > C   1       149829  6.6742753405549e-06

my $file2=shift;#Ctrl2
open (OUT,">$file2"."somatic_minus")||die;

open(IN1, $file2) or die ("can not open $file2\n");
while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	unless (exists $Ctrl{$SNP[0]}{$SNP[1]}) {
		print OUT $line;
	}else{
		if ($line=~/Chr/) {
			print OUT $line;
		}else{
		print OUT "$SNP[0]\t$SNP[1]\t".($SNP[2]-$Ctrl{$SNP[0]}{$SNP[1]})."\t$SNP[3]\t".(($SNP[2]-$Ctrl{$SNP[0]}{$SNP[1]})/$SNP[3])."\n";}
	}
}
close IN1;
close OUT;

my $file3=shift;#Stress
open (OUT,">$file3"."somatic_minus")||die;

open(IN1, $file3) or die ("can not open $file3\n");

while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	unless (exists $Stress{$SNP[0]}{$SNP[1]}) {
		print OUT $line;
	}else{
		if ($line=~/Chr/) {
			print OUT $line;
		}else{
		print OUT "$SNP[0]\t$SNP[1]\t".($SNP[2]-$Stress{$SNP[0]}{$SNP[1]})."\t$SNP[3]\t".(($SNP[2]-$Stress{$SNP[0]}{$SNP[1]})/$SNP[3])."\n";}
	}
}
close IN1;
close OUT;