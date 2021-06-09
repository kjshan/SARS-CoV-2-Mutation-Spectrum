#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


my $file=shift;#call SNP fas
my $refname=uc(shift);
$/="\n>";
my @ref;
my %SNP;
my %gap;
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
	chomp $line;
	$line=~s/\>//;
	my $line1=uc($line);
	my @seq=split /\n/,$line1;
	$seq[1]=~s/\s+//g;
	#print "$seq[0]";
	#unless ($seq[0]=~/EMBOSS/) {
		unless ($seq[0]=~/$refname/) {
			my @base=split //,$seq[1];
			$SNP{$seq[0]}=\@base;
			for (my $i=0;$i<=$#base ;$i++) {
				if ($base[$i] eq "-") {
					$gap{$seq[0]}{$i}+=1;
				}
			}
			#print ${$SNP{$seq[0]}}[1]."\t$seq[1]\n";
		}else {@ref=split //,$seq[1];}
	#}
}
close IN1;
#print @ref;
foreach my $seq (keys %SNP) {
	for (my $j=0;$j<=$#ref ;$j++) {
		#unless ($ref[$j] eq "-") {
			#print $j."\n";
			unless (exists $gap{$seq}{$j}) {
				#print "$seq\t$j\t".${$SNP{$seq}}[$j]."\n";
				unless (${$SNP{$seq}}[$j] eq $ref[$j]) {
					if ($ref[$j]=~/[A,T,C,G]{1}/ && ${$SNP{$seq}}[$j]=~/[A,T,C,G]{1}/) {
						print ($j+1);
						print "\t$ref[$j] > ${$SNP{$seq}}[$j]\t$seq\n";
					}
				}
				#unless (${$SNP{$seq}}[$j]=~/[A,T,C,G]{1}/ || ${$SNP{$seq}}[$j]=~/[A,T,C,G]{1}/) {
					#print $j."\n";
				#}
			}
		#}
	}
}