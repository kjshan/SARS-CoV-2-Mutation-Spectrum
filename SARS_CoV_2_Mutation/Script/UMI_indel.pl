#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#Sam file
#S100007198L1C004R0320729945     163     NC_045512.2     1       255     22S75M3S        =       1       28289   CTTATATCATTTGTATTTCCACATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACAAA    FAAAA@AFAA77BA@@A@FF?

my $file="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_Unique_reads/SCV2.Unique.mapped.sam.Junction.read.sam";#file name
my %Junction;
my %Junction_time;
my %Junction1;
open(IN, $file) or die ("can not open $file\n");

while (my $line=<IN>) {
		my @read=split /\t/,$line;
		if ($read[5] !~/S|H|P/i) {
			$Junction{$read[0]}.=$line;
			$Junction_time{$read[0]}+=1;
			unless (exists $Junction1{$read[0]}) {
				if ( $read[5] =~/I|D/i ) {
					if (exists $Junction_time{$read[0]} && $Junction_time{$read[0]}==2) {
						print  $Junction{$read[0]};
						#print $line;
					}elsif (exists $Junction_time{$read[0]} && $Junction_time{$read[0]}==1) {
						#print $line; 
						$Junction1{$read[0]}.=$line;
					}
				}
			}else{
				print  $Junction{$read[0]};
			}
		}
}
close IN;