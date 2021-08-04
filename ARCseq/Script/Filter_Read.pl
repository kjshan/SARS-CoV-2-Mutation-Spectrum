#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#@HD     VN:1.4
#@SQ     SN:NC_045512.2  LN:29903
#@PG     ID:STAR PN:STAR VN:2.7.3a       CL:STAR   --runThreadN 1000   --genomeDir /Dell/Dell9/shankj/Cov/RTmutation/SARS_CoV_2/STAR_index/   --readFilesIn /Dell/Dell9/shankj/Cov/RTmutation/Vero_SCV2_1.fq   /Dell/Dell9/shankj/Cov/RTmutation/Vero_SCV2_2.fq      --outFileNamePrefix /Dell/Dell9/shankj/Cov/RTmutation/STAR/SCV2   --outSJfilterCountUniqueMin 1   1   1   1      --outSJfilterCountTotalMin 1   1   1   1      --outSJfilterOverhangMin 12   12   12   12      --outSJfilterDistToOtherSJmin 0   0   0   0      --outFilterType BySJout   --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   --outFilterMismatchNoverReadLmax 0.04   --scoreGapNoncan -4   --scoreGapATAC -4   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJstitchMismatchNmax -1   -1   -1   -1      --chimScoreJunctionNonGTAG 0   --chimOutType WithinBAM   HardClip
#@CO     user command line: STAR --runThreadN 1000 --outFilterMultimapNmax 20 --genomeDir /Dell/Dell9/shankj/Cov/RTmutation/SARS_CoV_2/STAR_index/ --readFilesIn /Dell/Dell9/shankj/Cov/RTmutation/Vero_SCV2_1.fq /Dell/Dell9/shankj/Cov/RTmutation/Vero_SCV2_2.fq --outFileNamePrefix /Dell/Dell9/shankj/Cov/RTmutation/STAR/SCV2 --outFilterType BySJout --alignSJoverhangMin 8 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimOutType WithinBAM HardClip --chimScoreJunctionNonGTAG 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
#S100007198L1C001R0050074372     163     NC_045512.2     28403   255     100M    =       28598   295     GGTTTACCCAATAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAAGGCGTACCAATTAACACCA    IIA@A@I=F?@B@9A@>AIFIBAD@II@=38GGI?AB=I?8@H?><A8IG2@@II@>G@GFDD@@;A?@<@@5.9I.?<@AG72HD;H/@@D5A?C@?(@    NH:i:1  HI:i:1  AS:i:196        nM:i:1
#S100007198L1C001R0050074372     83      NC_045512.2     28598   255     100M    =       28403   -295    TATTTCTACTACCTAGGAACTGGGCCAGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAATACAC    @A9=AI?DIAAIIBAFI@AFAFIIICDIDDFI8II@IA?FIIBAAIFAIIAADIDDAFDIII=>@IDADAFII:@FIDAI;I?FFFAIII?AI@DBDIDI    NH:i:1  HI:i:1  AS:i:196        nM:i:1
#S100007198L1C001R0050074380     163     NC_045512.2     26850   255     10M200N10M    =       26918   167     TGGTCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTCTGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCTGAGCTG    B4=DIBBAIBBAA<@-BBBIABBIAABIC<7CFB@@+D@B?AI>I6BA:,<BDBBAA+C@B8(B4BI)6D<AIAB0BBBE@,BBIB?0AB8*32.B2%A.    NH:i:1  HI:i:1  AS:i:195        nM:i:1

my $file=shift;#file name
my @file=split /\//,$file;

#my $out1=shift;
#open (OUT1,">$out1")||die;
open(IN, $file) or die ("can not open $file[-1]\n");

while (my $line=<IN>) {
	#chomp $line;
	if ($line=~/^@/) {
		print  $line;
	}else{
		my @read=split /\t/,$line;
		if ($read[2] !~/\*/ && 
			$read[5] !~/I|D|N|S|H|P/i && 
			$read[11] eq "NH:i:1" && 
			$read[4] > 0 ) {
			#&& $line !~/NM:i:2[0-9]/i && $line !~/NM:i:1[0-9]/i && ($line =~/NM:i:1/i || $line =~/NM:i:0/i || $line =~/NM:i:2/i)
			#foreach my $mid (@read) {
				#if ($mid=~/NM:i:/) {
					#my @mismatch=split /NM:i:/,$mid;
					#if ($mismatch[1]<=5) {
						print "$read[3]~".($read[3]+length($read[9])-1)."~$line";
						
					#}
				#}
			#}
			
		}
	}
}
close IN;
