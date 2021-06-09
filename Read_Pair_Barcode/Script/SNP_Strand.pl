#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#+ L-A     1682    T       144
my %count_base;
my $file=shift;#./Coverage.default
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	my @Chr=split /\s+/,$line;
	if ($line=~/\d/ && $Chr[4]>=100) {
		
		if ($Chr[1]=~/L-A/  || $Chr[1]=~/20S/ || $Chr[1]=~/23S/ || 
			$Chr[1]=~/L-BC-La/  || $Chr[1]=~/L-BC-2/ ) {
			$count_base{$Chr[1]}{$Chr[3]}+=$Chr[4];
		}elsif ($Chr[1] !~/"Mito"/) {
			$count_base{"mRNA"}{$Chr[3]}+=$Chr[4];
		}
	}
}
close IN;


#-       I_139527        T > C   2       0.0198019801980198      0.98019801980198        101
my $file1=shift;#Ctrl2.Consensuse.Mismatch.Strand
open(IN1, $file1) or die ("can not open $file1\n");
my %SNP;
while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	if ($SNP[4]<=0.01 && $SNP[5]>0.9) {
		my @Chr=split /_/,$SNP[1];
		if ( $Chr[0]=~/L-A/  || $Chr[0]=~/20S/ || $Chr[0]=~/23S/ || 
			$Chr[0]=~/L-BC-La/  || $Chr[0]=~/L-BC-2/ ) {
			$SNP{$Chr[0]}{$SNP[2]}+=$SNP[3];
		}elsif ($Chr[0] !~/"Mito"/) {
			$SNP{"mRNA"}{$SNP[2]}+=$SNP[3];
		}
		
	}
}
close IN1;

print "Chr\tSNP\tSNP_number\tBase_number\tMutation_Freq\n";
foreach my $Chr(keys %SNP) {
	foreach my $SNP(keys %{$SNP{$Chr}}) {
		my @ref=split /\s+/,$SNP;
		print $Chr."\t".$SNP."\t".$SNP{$Chr}{$SNP}."\t".$count_base{$Chr}{$ref[0]}."\t".($SNP{$Chr}{$SNP}/$count_base{$Chr}{$ref[0]})."\n";
	}
}
