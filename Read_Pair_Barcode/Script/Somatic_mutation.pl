#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


#-       I_139527        T > C   2       0.0198019801980198      0.98019801980198        101
my $file1=shift;#Ctrl2
open(IN1, $file1) or die ("can not open $file1\n");
my %SNP;
while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	if ($SNP[4]<=0.01 && $SNP[5]>0.9) {
		$SNP{$SNP[1]}{$SNP[2]}+=$SNP[3];
	}
}
close IN1;


my $file2=shift;#Ctrl2
open(IN1, $file2) or die ("can not open $file2\n");
my %SNP_somatic;
my %SNP_somatic_Ctrl;
while (my $line=<IN1>) {
	my @SNP=split /\t/,$line;
	if ($SNP[4]<=0.01 && $SNP[5]>0.9 && exists $SNP{$SNP[1]}{$SNP[2]}) {
		my @Chr=split /_/,$SNP[1];
		if ( $Chr[0]=~/L-A/  || $Chr[0]=~/20S/ || $Chr[0]=~/23S/ || 
			$Chr[0]=~/L-BC-La/  || $Chr[0]=~/L-BC-2/ ) {
			$SNP_somatic{$Chr[0]}{$SNP[2]}+=$SNP[3];
			$SNP_somatic_Ctrl{$Chr[0]}{$SNP[2]}+=$SNP{$SNP[1]}{$SNP[2]};
		}elsif ($Chr[0] !~/"Mito"/) {
			$SNP_somatic{"mRNA"}{$SNP[2]}+=$SNP[3];
			$SNP_somatic_Ctrl{"mRNA"}{$SNP[2]}+=$SNP{$SNP[1]}{$SNP[2]};
		}
	}
}
close IN1;



print "Chr\tSNP\tCtrl\tStress\n";
foreach my $Chr(keys %SNP_somatic) {
	foreach my $SNP(keys %{$SNP_somatic{$Chr}}) {
		print $Chr."\t".$SNP."\t".$SNP_somatic_Ctrl{$Chr}{$SNP}."\t".$SNP_somatic{$Chr}{$SNP}."\n";
	}
}
