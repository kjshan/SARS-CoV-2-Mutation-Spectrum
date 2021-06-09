#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
my $file3=shift;#pileup file
my @file=split /\./,$file3;
my @JC=split /\-/,$file[0];#start 0 and end 1

if ($file[1]>1) {
my $file2="/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_UMI2_20/SCV2_Unique_UMI_gt1.txt";#bed
my %UMI;
#print $file[0]."pp\n";
open(IN, $file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	
	if ( $read[0] eq "$file[0]") {#&& $read[1] <= 100
		#print $read[0]."pp\n";
		my @read1=split /,/,$read[2];
		foreach my $readname(@read1) {
			if ($readname=~/:/) {
				my @read2=split /:/,$readname;
				$UMI{"1:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
				$UMI{"2:$read2[2]"}="$read2[0]:$read2[1]";#read name -> 1:2
				#print "1:$read2[2]\n";
			}
			
		}
	}
}
close IN;


#NC_045512.2     8       G       3       ^~.^~.^~.       III    



#print "BQ\tBarcode\tSNP\tPos\tDis\tCovered_reads_number\tMismatch_reads_number\n";
print "AverBQ\tBarcode\tPos\tDis\tnon_PCR\tMatch\n";
my $count=0;
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	chomp $line;
	my @name=split /\t/,$line;

	if ($name[1] !=4402 && 
		$name[1] !=5062 && 
		$name[1] !=8782 && 
		$name[1] !=28144 && 
		($name[1]<=$JC[0]-15 || $name[1]>=$JC[1]+15)) {
		#unless ($name[4]=~/[\>|\<]/) {
			my $BQ=0;
			$name[4]=~s/\^.//g;
			$name[4]=~s/\$//g;

				#Dis
				my $left=$name[1]-$JC[0];
				my $right=$name[1]-$JC[1];
				my $UMI_Dis;
				my %cover=();
				if (abs($left)<abs($right)) {
					$UMI_Dis=abs($left);
				}else{$UMI_Dis=abs($right);}
			if ($UMI_Dis>=15) {
			
			if (($name[4]=~/[\.|,]{$name[3]}/gi )) {#|\d
				my %BQ=();
				my @BQ=split //,$name[5];
				foreach my $BQ1(@BQ) {
					$BQ+=(ord($BQ1)-33);
				}

				my @read_name=split /,/,$name[$#name];my %hash;
				#my $name;
				for (my $i=0;$i<=$#read_name ;$i++) {
					if ($read_name[$i]=~/S100/) {
						my $PCR=$UMI{$read_name[$i]};
						#$name.="$PCR:$read_name[$i],";
						$hash{$PCR}+=1;
					}
					
				}
				my $len=keys %hash;
				if ($len>=2) {
					print $BQ/($#BQ+1)."\t$file[0]\t$name[1]\t$UMI_Dis\t$len\tMatch\n";#\t$name
				}
			}elsif ( $name[4]=~/[A]{$name[3]}/gi || $name[4]=~/[T]{$name[3]}/gi ||  $name[4]=~/[C]{$name[3]}/gi || $name[4]=~/[G]{$name[3]}/gi ) {
				my %BQ=();
				my @BQ=split //,$name[5];
				foreach my $BQ1(@BQ) {
					$BQ+=(ord($BQ1)-33);
				}

				my @read_name=split /,/,$name[$#name];my %hash;
					#my $name;
					for (my $i=0;$i<=$#read_name ;$i++) {
					if ($read_name[$i]=~/S100/) {
						my $PCR=$UMI{$read_name[$i]};
						#$name.="$PCR:$read_name[$i],";
						$hash{$PCR}+=1;
					}
					
					}
				my $len=keys %hash;
				if ($len>=2) {
				print $BQ/($#BQ+1)."\t$file[0]\t$name[1]\t$UMI_Dis\t$len\tMismatch\n";
				}#\t$name
			}
			}

	}
}
close IN1;

}

