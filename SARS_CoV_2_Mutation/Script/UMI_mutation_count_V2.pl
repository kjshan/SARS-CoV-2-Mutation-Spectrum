#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#UMI     UMI_read_pair_number    read
#1000-28254      6       920:28354:S100007198L1C001R0060138612,832:975:S100007198L1C004R0520434463,740:977:S100007198L1C008R0460373726,960:28337:S100007198L1C001R0020730201,991:28357:S100007198L1C003R0010708388,920:28354:S100007198L1C008R0600769019,

my $file1=shift;#bed
my %read;
my %UMI_infact_count;
open(IN, $file1) or die ("can not open $file1\n");
while (my $line=<IN>) {
	my @read=split /\t/,$line;
	if ($read[1] !~/UMI/ && $read[1] >= 2 ) {#&& $read[1] <= 100
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

#open (OUT,">UMI_read.summarise")||die;
my %mutation;
#NC_045512.2     8       G       3       ^~.^~.^~.       III     1,1,1   S100007198L1C007R0580859955:69:4579,S100007198L1C008R0340265535:68:29053,S100007198L1C003R0450748162:66:29375
my $file3=shift;#mpileup
print "Position\tRef > Alt\tAlt_UMI_number\tAlt_read_number\tMutation_UMI:Mutation_UMI_Read_W/O_Duplication:Mutation_UMI_Read:Mutation_UMI_Total_Read:Mutation_UMI_same_Read\n";
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
		chomp $line;
		my $site=0;
		my $gap=0;my $num=0;
		my $mid=0;
		my @name=split /\t/,$line;
		if ($name[3] > 0) {
			$name[4]=~s/\^.//g;
			$name[4]=~s/\$//g;
			my %mutation = ();
			my %UMI_read = ();
			my %match_read = ();
			#my %UMI_count= ();
			if ($name[4]=~/A|T|C|G/i) {
				my @read=split /,/,$name[7];
				while ($name[4]=~/[A|T|C|G]{1}/gi) {#|\d
					$site=pos $name[4];
					my $mid1=$read[($site-1)];#read
					if (exists $read{$mid1}) {
							$mutation{uc($&)}+=1;#SNP count
							my @UMI=split /:/,$read{$mid1};#UMI:920:28354
							#$UMI_count{"$UMI[1]-$UMI[2]"}+=1;#该位点有多少个UMI cover
							#$Mutation_UMI_count{uc($mismatch[$i])}{"$UMI[1]-$UMI[2]"}+=1;#每个UMI对应多少个突变的reads
							my $mid1111="$UMI[1]:$UMI[2]";
							$UMI_read{uc($&)}{$UMI[0]}.="$mid1111,";#该突变该UMI对应的reads pair
							#$UMI_read{uc($mismatch[$i])}{"$UMI[0]:$UMI[1]:$UMI[2]"}+=1;#该突变该UMI中非PCR duplication的reads数目
							#print $site."\t$mid1\t$name[2] > $&\t$name[1]\n";
					}
				}
				while ($name[4]=~/[\.|\,]{1}/gi) {#|\d
					$site=pos $name[4];
					my $mid1=$read[($site-1)];#read
					if (exists $read{$mid1}) {
						my @UMI=split /:/,$read{$mid1};
						$match_read{$UMI[0]}+=1
					}
				}
				foreach my $mutation(keys %UMI_read) {
					#my $len1=keys %UMI_count;#UMI_count
					my $len3=keys %{$UMI_read{$mutation}};
					#Position\tRead_count\tRef > Alt\tAlt_UMI_number\tAlt_read_number
					print  "$name[1]\t"."$name[2] > $mutation\t$len3\t$mutation{$mutation}";#Position\tRead_Number\tUMI_count\tRef > Alt\t
					foreach my $UMI_name (keys %{$UMI_read{$mutation}}) {
						my %hash=();
						my $len4=0;
						my $read123=$UMI_read{$mutation}{$UMI_name};
						my @read13=split /,/,$read123;
						foreach my $key(@read13) {
							if ($key=~/:/) {
								$hash{$key}+=1;
								$len4+=1;
							}
						}
						my $len2=keys %hash;
						if (exists $match_read{$UMI_name}) {
							print  "\t$UMI_name:$len2:".$len4.":".$UMI_infact_count{$UMI_name}.":".$match_read{$UMI_name};
						}else{
							print  "\t$UMI_name:$len2:".$len4.":".$UMI_infact_count{$UMI_name}.":0";
						}
					}
					print  "\n";
				}
			}
		}
}
close IN1;
