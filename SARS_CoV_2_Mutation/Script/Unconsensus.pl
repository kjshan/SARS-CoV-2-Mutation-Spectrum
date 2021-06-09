#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2     8       G       3       ^~.^~.^~.       III    
my $file3=shift;#pileup file
my @file=split /\./,$file3;
my @JC=split /\-/,$file[0];#start 0 and end 1
if ($file[1]>1) {
my $count=0;
#print "BQ\tBarcode\tSNP\tPos\tDis\tCovered_reads_number\tMismatch_reads_number\n";
print "BQ\tBarcode\tPos\tDis\tMatch_reads_number\tMismatch_reads_number\n";
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	#$count+=1;
	#if ($count !=4402 && 
		#$count !=5062 && 
		#$count !=8782 && 
		#$count !=28144 && 
		#($count<=$JC[0]-15 || $count>=$JC[1]+15)) {
		chomp $line;
		my @name=split /\t/,$line;
	
	if ($name[1] !=4402 && 
		$name[1] !=5062 && 
		$name[1] !=8782 && 
		$name[1] !=28144 && 
		($name[1]<=$JC[0]-15 || $name[1]>=$JC[1]+15)) {
		#unless ($name[4]=~/[\>|\<]/) {
			my $left=$name[1]-$JC[0];
			my $right=$name[1]-$JC[1];
			my $UMI_Dis;
			my %cover=();
			if (abs($left)<abs($right)) {
				$UMI_Dis=abs($left);
			}else{$UMI_Dis=abs($right);}

				$name[2]=uc($name[2]);
				my %mutation=();
				my %mis=();my %same=();
				my $same_read1=();

				$name[4]=~s/\^.//g;
				$name[4]=~s/\$//g;
				my @BQ=split //,$name[5];
				while ($name[4]=~/[A|T|C|G]{1}/gi) {#|\d
					my $site=pos $name[4];

					my $SNP11=$name[2]." > ".uc($&);#print "$SNP11\n";
					my $BQ1=(ord($BQ[($site-1)])-33);
					$mutation{$BQ1}{$SNP11}+=1;#SNP count
					$mis{$BQ1}+=1;#$cover{$BQ1}+=1;

				}
				while ($name[4]=~/[\.|\,]{1}/gi) {#|\d
					my $site=pos $name[4];
					my $BQ1=(ord($BQ[($site-1)])-33);
					$same{$BQ1}+=1;#$cover{$BQ1}+=1;
				}

				foreach my $BQ(keys %same) {
					if (exists $mis{$BQ}) {
						#foreach my $SNP(keys %{$mutation{$BQ}}) {
							#print "$BQ\t$SNP\t$file[0]\t$name[1]\t$UMI_Dis\t".$cover{$BQ}."\t".$mutation{$BQ}{$SNP}."\n";
							print "$BQ\t$file[0]\t$name[1]\t$UMI_Dis\t".$same{$BQ}."\t".$mis{$BQ}."\n";
						#}
					}else {
						#print "$BQ\tNULL\t$file[0]\t$name[1]\t$UMI_Dis\t".$cover{$BQ}."\t0\n";
						print "$BQ\t$file[0]\t$name[1]\t$UMI_Dis\t".$same{$BQ}."\t0\n";
					}
				}

				foreach my $BQ(keys %mis) {
					unless (exists $same{$BQ}) {
						print "$BQ\t$file[0]\t$name[1]\t$UMI_Dis\t0\t$mis{$BQ}\n";
						#foreach my $SNP(keys %{$mutation{$BQ}}) {
							#print "$BQ\t$SNP\t$file[0]\t$name[1]\t$UMI_Dis\t".$cover{$BQ}."\t".$mutation{$BQ}{$SNP}."\n";
						#}
					}
				}
		#}
	}
}
close IN1;

}

