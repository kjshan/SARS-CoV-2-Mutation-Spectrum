#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#65      25382   S100007198L1C008R0320875124     163     NC_045512.2     30      255     36M25315N64M    =       25496   25571   AACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTATGGATTTGTTTATGAGAATCTTCACAATTGGAACTGTAACTTTGAAGCAAG >:FF?=CI?=:??AHI;A;AI=A7DAI?@:AI:>IDGA@;AHE@=CDA8D@@ADADIDDD9AI>+;(/I?D836AA/AIIA?FDIAA+EAAA?>3.0;=H    NH:i:1  HI:i:1  AS:i:180        nM:i:1
#64      25381   S100007198L1C008R0080009028     163     NC_045512.2     25      255     19M1D21M25315N60M       =       25496   25571   TAACAAACCAACCAACTTTGATCTCTTGTAGATCTGTTCTCTAAACGAACTTATGGATTTGTTTATGAGAATCTTCACAATTGGAACTGTAACTTTGAAG D@@FA@@FGAAEG@@FDDDI?D?DFDAID<I:DBDI@@HA@DAA?>I?AGADA@IA8AADIDDDBD<AI@0DID?7@3A@DA<H@>?AHA?>2*DDC;4-    NH:i:1  HI:i:1  AS:i:190        nM:i:1

my %UMI;
my %read;
my $file=shift;#file name
open(IN, $file) or die ("can not open $file\n");
#28M2D18M28190N54M
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	unless ($line =~/^@/) {
		if ( $read[5] =~/N/i) {
			my @migar=split /N/,$read[5];
			my $start=$read[3];
			my $len=$migar[0];
			my $end=0;
			my @match=split /M/,$len;
			for (my $i=0;$i<$#match;$i++) {
				my $match=$match[$i];
				unless ($match=~/D|I/) {
					$start+=$match;
				}elsif ($match=~/I/) {
					my @mid=split /I/,$match;
					$start+=$mid[1];
				}elsif ($match=~/D/) {
					my @mid=split /D/,$match;
					$start+=$mid[1];
					$start+=$mid[0];
				}
				
			}
			unless ($match[$#match]=~/D|I/) {
				$end=$start+$match[$#match];
			}else{
				my @mid=split /[D,I]/,$match[$#match];
				$end=$start+$mid[1];
			}
			my $UMI1=($start);
			my $UMI2=($end-1);
			$UMI{"$UMI1-$UMI2"}.=$read[0].",";#read -> UMI
			
		}
		$read{$read[0]}.="$read[3]:";
	}
}
close IN;

my %UMI_count;
my %UMI_count2;
foreach my $UMI(sort keys %UMI) {
	my @reads=split /,/,$UMI{$UMI};
	#if ($#reads>=1) {
		my %ha=();
		foreach my $read(grep{++$ha{$_}<2} @reads) {
			$UMI_count{$UMI}.=$read{$read}.$read.",";
			$UMI_count2{$UMI}+=1;
		}
	#}
}
foreach my $UMI(sort keys %UMI_count) {
	print "$UMI\t$UMI_count2{$UMI}\t$UMI_count{$UMI}\n";
}