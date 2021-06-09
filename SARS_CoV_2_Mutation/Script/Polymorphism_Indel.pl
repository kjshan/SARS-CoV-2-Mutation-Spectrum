#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#Chr	Position	Read	Ref	Alt1	Alt2	Alt3
#NC_045512.2	1	261	A
#NC_045512.2	2	379	T
#NC_045512.2	3	491	T	A:7:0.0142566191446029	C:3:0.00610997963340122	G:6:0.0122199592668024
#NC_045512.2	4	718	A	C:9:0.0125348189415042	G:22:0.0306406685236769	T:38:0.052924791086351

my $file1="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/bin/SCV2.Unique.summarize";#bed
my %POS;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	unless ($read[2] =~ /Position/) {
		$POS{$read[1]}=$read[2];
	}
}
close IN;


#70      70      26      Deletion        AC      -2      8,      S100007198L1C001R0450234805,

my $file2="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/UMI_indel.information.background.default";#bed
my %Indel;
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	my @num=split /,/,$read[7];
	$POS{$read[2]}+=$#num;
	$Indel{$read[2]}{"$read[4]\t$read[5]"}+=$#num;
	#print $read[2]."\t"."$read[4]\t$read[5]\n";
}
close IN;


#1273-29234	3	3	3	1241	Deletion	GTG	-3	-32	63,51,18,	2:S100007198L1C007R0320891052,1:S100007198L1C005R0090614797,1:S100007198L1C003R0010393465,	1129:1179:2:S100007198L1C007R0320891052,1191:29265:1:S100007198L1C005R0090614797,1224:29266:1:S100007198L1C003R0010393465,	3	3
#1273-5946	2	2	2	1241	Deletion	GTG	-3	-32	44,34,	2:S100007198L1C001R0060022787,2:S100007198L1C002R0160519481,	1077:1198:2:S100007198L1C001R0060022787,1108:1208:2:S100007198L1C002R0160519481,	2	2
#read[13] SNP_read_number
my $file3="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/Consensuse.Indel.txt";#bed

open(IN, $file3) or die ("can not open $file3\n");
while (my $line=<IN>) {
	chomp $line;
	my @read=split /\t/,$line;
	my @num=split /,/,$read[7];
	$POS{$read[4]}+=$read[13];# SNP_read_number
	$Indel{$read[4]}{"$read[6]\t$read[7]"}+=$read[13];
	#print "$read[4]\t$read[6]\t$read[7]\t$read[13]\n";
}


open(IN, $file3) or die ("can not open $file3\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	my @num=split /,/,$read[7];
	print $Indel{$read[4]}{"$read[6]\t$read[7]"}/$POS{$read[4]}."\t".$line;
}
close IN;
