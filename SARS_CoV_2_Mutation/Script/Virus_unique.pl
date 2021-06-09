#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#S100007198L1C001R0050794690     163     NC_045512.2     27761   255     66M34S  =       27783   122     TGAACTTTCATTAATTAACTACTATTTGTNATTTTTACCCTTACTCCTATACANTGTTTTAATTATCCTTATATACTTTTNTTACANACTTAAACTCCAA    A+AACAD?FAAA@AAA/B4A/'A?AAC)B!?@AAAAA+F'AB0>C1%BAA8>>!A*AAAA@AA?AD(%7C@A+BA37AAA!(A5:A!7)AA2?A%B87AB    NH:i:1  HI:i:1  AS:i:146        nM:i:8

my $file1=shift;#human
my $file3=shift;#virus read

my %read;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	#chomp $line;
	unless ($line=~/^@/) {
		my @name=split /\t/,$line;
		unless ($name[2] eq "*") {
			$read{$name[0]}.=$line;
		}
	}
}
close IN;


my %virus;
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	unless ($line=~/^@/) {
		my @name=split /\t/,$line;#
		unless ($name[2] eq "*") {
				unless (exists $read{$name[0]}) {
						if ( $name[11] eq "NH:i:1") {
							print  $line;
						}
				}
			}
	}else {print $line;}
}
close IN1;
