#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
#perl RemoveRefGap.pl *fasta reference_name

my $in=shift;#MSA file; fasta 
my $ref=shift;
my %gap;
#>reference
#-------ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATC--TCTTGTAGATCTGT-TCTCTAAACGAACTTTAAAA
#>seq1
#---------------------------------------CAACCAACTTTCGATCAATCTTGTAGATCTGTATCTCTAAACGAACTTTAAAA
#>seq2
#-------ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATC--TCTTGTAGATCTGTATCTCTAAA--AACTTTAAAA
$/="\n>";
open (IN,$in);
my $count=0;
while(my $line=<IN>){
	$line=~s/\>//;
	my @line=split /\n/,$line;
	if ($line[0] eq $ref) {
		while ($line[1]=~/[\-]/gi) {
			my $site=pos $line[1];
			$gap{($site-1)}+=1;
		}
	}
}


open (IN,$in);
while(my $line=<IN>){
	$line=~s/\>//;
	my @line=split /\n/,$line;
	print ">$line[0]\n";
	my @base=split //,$line[1];
	for (my $j=0;$j<=$#base;$j++) {
		unless (exists $gap{$j}) {
			print $base[$j];
		}
	}
	print "\n";
}
close IN;