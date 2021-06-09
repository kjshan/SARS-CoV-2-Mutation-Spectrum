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

#my $out1=$file[-1]."_Junction.read";
#open (OUT1,">$out1")||die;
my %JCR;
open(IN, $file) or die ("can not open $file[-1]\n");
my %match;
while (my $line=<IN>) {
	#chomp $line;
	unless ($line=~/^@/) {
		my @read=split /\t/,$line;
		#if ($read[5] !~/S|I|D|H|P/i) {
			if ($read[5] =~/N/) {
				my @site1=split /M/,$read[5];
				my $start=($read[3]+$site1[0]);
					my @site=split /N|M/;
					my $mid=0;
					foreach my $len(@site) {
						unless ($len=~/D|I/) {
							$mid+=$len;
						}else {
							
						}
					}
					my $end=$mid-1;
					$JCR{"$read[2]\t$read[0]\t$start\t$end\n"}+=1;
				#}
			}else {
				$match{$read[0]}+=1;
			}
		#}
	}
}
close IN;
foreach my $key(keys %JCR) {
	if ($JCR{$key}==2) {
		print $key;
	}else {
		my @read=split /\t/,$key;
		if ($match{$read[1]}==1) {
			print $key;
		}
	}
}