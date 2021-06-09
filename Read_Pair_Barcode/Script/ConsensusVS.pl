#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#XII:282988:XII:282989.10.sam

my $file1=shift;#bed
my %UMI_count1;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @UMI=split /\./,$line;
	$UMI_count1{$UMI[0]}=$UMI[1];
}
close IN;


my $file3=shift;#mpileup
open (OUT,">$file3.Consensus")||die;
open (OUT1,">$file3.Unconsensus")||die;
#L-BC-La 17      c       13      ..^~.^~.^~.^~.^~.^~.^~.^~.^~.^~.^~.     FEEEECEEEEEEE   10,2,1,1,1,1,1,1,1,1,1,1,1      L-BC-La:8:L-BC-La:56:V350014000L3C005R0200477440,L-BC-La:16:L-BC-La:82:V350014000L3C001R0290916070,L-BC-La:17:L-BC-La:69:V350014000L3C001R0150208226,L-BC-La:17:L-BC-La:159:V350014000L3C001R0020295100,L-BC-La:17:L-BC-La:93:V350014000L3C006R0020407205,L-BC-La:17:L-BC-La:374:V350014000L3C002R0480134999,L-BC-La:17:L-BC-La:85:V350014000L3C004R0050639294,L-BC-La:17:L-BC-La:34:V350014000L3C002R0300644641,L-BC-La:17:L-BC-La:89:V350014000L3C003R0090529831,L-BC-La:17:L-BC-La:175:V350014000L3C003R0270405748,L-BC-La:17:L-BC-La:117:V350014000L3C002R0460366032,L-BC-La:17:L-BC-La:160:V350014000L3C001R0500464871,L-BC-La:17:L-BC-La:138:V350014000L3C003R0311129640
#print OUT "Chr_Position\tRef > Alt\tAlt_number\tAlt_fraction\tRef_fraction\tCover\n";
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
		chomp $line;
		my @name=split /\t/,$line;
		if ($name[3] >= 2 && $name[0]!~/Mito/) {
			my $site=0;
			my $gap=0;my $num=0;
			my $mid=0;
			my %mutation=();
			my %same=();
			my %SNP_count=();

			$name[4]=~s/\^.//g;
			$name[4]=~s/\$//g;
			$name[2]=uc($name[2]);
			#reads
			my @read=split /,/,$name[7];
			#BQ
			my %BQ_mis=();
			my %BQ_same=();
			my @BQ=split //,$name[5];
			my %un_Consensus_SNP=();
			my %un_Consensus_Cover=();

			while ($name[4]=~/[A|T|C|G]{1}/gi) {#|\d
				$site=pos $name[4];
				my $mid1=$read[($site-1)];#read
				my @UMI=split /:/,$mid1;#$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]
				my $base=uc($&);
				$UMI[0]=~s/\d+\-//;$UMI[0]=~s/\d+\-//;
				#print "$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]"."\n";
				$mutation{$base}{"$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]"}+=1;
				$BQ_mis{$base}{"$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]"}+=(ord($BQ[($site-1)])-33);

				$un_Consensus_SNP{(ord($BQ[($site-1)])-33)}+=1;
				$un_Consensus_Cover{(ord($BQ[($site-1)])-33)}+=1;
			}
			while ($name[4]=~/[\.|\,]{1}/gi) {#|\d
				$site=pos $name[4];
				my $mid1=$read[($site-1)];#read
				my @UMI=split /:/,$mid1;#$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]
				$UMI[0]=~s/\d+\-//;$UMI[0]=~s/\d+\-//;
				#print $UMI[0]."\n";
				$same{"$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]"}+=1;
				$BQ_same{"$UMI[0]:$UMI[1]:$UMI[2]:$UMI[3]"}+=(ord($BQ[($site-1)])-33);
				$un_Consensus_Cover{(ord($BQ[($site-1)])-33)}+=1;
			}
			my $cover=0;
			my $same=0;
			my %AverBQ_Cover;
			foreach my $UMI(keys %same) {
				if ($same{$UMI}==$UMI_count1{$UMI}) {
					$AverBQ_Cover{int($BQ_same{$UMI}/$same{$UMI})}+=1;
				}
			}
			foreach my $mutation(keys %mutation) {
				#print $mutation."\n";;
				foreach my $UMI_name (keys %{$mutation{$mutation}}) {
					if ($mutation{$mutation}{$UMI_name}==$UMI_count1{$UMI_name}) {
						$AverBQ_Cover{int($BQ_mis{$mutation}{$UMI_name}/$mutation{$mutation}{$UMI_name})}+=1;
						$SNP_count{int($BQ_mis{$mutation}{$UMI_name}/$mutation{$mutation}{$UMI_name})}+=1;
					}
					
				}
				
			}
			foreach my $BQ(keys %AverBQ_Cover) {
				print OUT "$BQ\t$name[0]_$name[1]\t$AverBQ_Cover{$BQ}\tCover\n";
			}
			foreach my $BQ(keys %SNP_count) {
				print OUT "$BQ\t$name[0]_$name[1]\t$SNP_count{$BQ}\tMis\n";
			}

			foreach my $BQ(keys %un_Consensus_SNP) {
				print OUT1 "$BQ\t$name[0]_$name[1]\t".$un_Consensus_SNP{$BQ}."\tMis\n";
			}
			foreach my $BQ(keys %un_Consensus_Cover) {
				print OUT1 "$BQ\t$name[0]_$name[1]\t".$un_Consensus_Cover{$BQ}."\tCover\n";
			}
		}
}
close IN1;

