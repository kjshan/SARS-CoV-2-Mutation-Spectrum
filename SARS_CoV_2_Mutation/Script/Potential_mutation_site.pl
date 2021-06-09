#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#1000-28515      3       958:28649:S100007198L1C008R0210064379,971:28648:S100007198L1C002R0340492848,951:28615:S100007198L1C008R0440693383,
#1000-29306      2       853:960:S100007198L1C005R0320636011,706:957:S100007198L1C007R0170680541,

my $file1=shift;#bed
my %read;
my %UMI_infact_count;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	if ($read[1] >= 2 && $read[1] <= 20 && $read[1] !~/UMI_read_pair_number/) {
		$UMI_infact_count{$read[0]}=$read[1];
		my @read1=split /,/,$read[2];
		foreach my $readname(@read1) {
			if ($readname=~/:/) {
				my @read2=split /:/,$readname;
				$read{$read2[2]}="$read[0]:$read2[0]:$read2[1]";#read name -> UMI:1:2
			}
		}
	}
}
close IN;





my %mutation;
#NC_045512.2     8       G       3       ^~.^~.^~.       III     1,1,1   S100007198L1C007R0580859955:69:4579,S100007198L1C008R0340265535:68:29053,S100007198L1C003R0450748162:66:29375
#opendir(DIR,"/Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/STAR/Junction_read/Mpileup_Split/")||die;
#my @files=grep (/^x/,readdir DIR);

my $file3=shift;#mpileup
#print "Site\tRef\n";
#foreach my $file3 (@files) {
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
		chomp $line;
		my %hash;
		my $site=0;
		my @name=split /\t/,$line;
		$name[4]=~s/\^.//g;
		$name[4]=~s/\$//g;
		#my %UMI_count= ();
		if ($name[4]=~/A|T|C|G|\.|\,/i) {
			my @read=split /,/,$name[7];
			while ($name[4]=~/[A|T|C|G|\.|\,]{1}/gi) {#|\d
				$site=pos $name[4];
				my $mid1=$read[($site-1)];#read
				if (exists $read{$mid1}) {
					my @UMI=split /:/,$read{$mid1};
					my $mid1111="$UMI[1]:$UMI[2]";
					$hash{$mid1111}+=1;
					my $len=keys %hash;
					if ($len >=2) {
						print "$name[1]\t$name[2]\n";
						last;
					}
				}
				
			}
		}
}
close IN1;
#}
#closedir DIR;