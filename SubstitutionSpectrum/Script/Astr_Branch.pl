#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
my $Species=shift;
#S1	Bat_MF593268.1
my $file1="SeqCodes";#
my %name;
open(IN1, $file1) or die ("can not open $file1\n");
while (my $line=<IN1>) {
	chomp $line;
	my @seq=split /\s+/,$line;
	$name{$seq[0]}=$seq[1];
}
close IN1;

#S2 N1
my $file2="New";
my %C_P;
open(IN1, $file2) or die ("can not open $file2\n");
while (my $line=<IN1>) {
	if ($line!~/^name/ && $line!~/^\#/ && $line !~/root/ && $line=~/\d/) {
		chomp $line;
		my @seq=split /\s+/,$line;
		$C_P{$seq[0]}=$seq[1];
	}
}
close IN1;

#S1 \n seq \n
my %seq;
my $file="seq.joint.txt";#call SNP fas
$/="\n>";
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\>//;
	my @seq=split /\n/,$line;
	$seq{$seq[0]}=$seq[1];
}
close IN;

open (OUT,">SNP")||die;
my %SNP;
my %G_base;
my %C_base;
foreach my $seq(sort keys %C_P) {
	my $ref=$C_P{$seq};
	my $ref_seq=$seq{$ref};
	my @ref_base=split //,$ref_seq;
	my $branch=$seq;
	if (exists $name{$seq}) {
		$branch=$name{$seq};
	}
	my $branch_seq=$seq{$seq};
	my @branch_base=split //,$branch_seq;
	for (my $j=0;$j<=$#ref_base;$j++) {
			if ($ref_base[$j]=~/[A,T,C,G]{1}/) {
				if (uc($ref_base[$j]) eq "G") {
					$G_base{"$ref:$branch"}+=1;
				}elsif (uc($ref_base[$j]) eq "C") {
					$C_base{"$ref:$branch"}+=1;
				}
				if ($branch_base[$j]=~/[A,T,C,G]{1}/) {
					unless ($ref_base[$j] eq $branch_base[$j]) {
						print OUT ($j+1);
						print OUT "\t$ref_base[$j]$branch_base[$j]\t$ref:$branch\n";
						$SNP{"$ref:$branch"}{"$ref_base[$j]$branch_base[$j]"}+=1;
					}
				}
			}
	}
}
close OUT;
print "Branch\tSNP\tCount\tG\tC\tSpecies\n";
foreach my $key1 (sort keys %SNP) {
	$SNP{$key1}{"GT"}+=0;
	$SNP{$key1}{"CA"}+=0;
	$SNP{$key1}{"CT"}+=0;
	$SNP{$key1}{"GA"}+=0;

	$SNP{$key1}{"GC"}+=0;
	$SNP{$key1}{"CG"}+=0;
	$SNP{$key1}{"AT"}+=0;
	$SNP{$key1}{"TA"}+=0;

	$SNP{$key1}{"TG"}+=0;
	$SNP{$key1}{"TC"}+=0;
	$SNP{$key1}{"AG"}+=0;
	$SNP{$key1}{"AC"}+=0;
     foreach my $key2 (keys %{$SNP{$key1}}) {
        print $key1."\t".$key2."\t".$SNP{$key1}{$key2}."\t".$G_base{$key1}."\t".$C_base{$key1}."\t$Species\n";
    }
}