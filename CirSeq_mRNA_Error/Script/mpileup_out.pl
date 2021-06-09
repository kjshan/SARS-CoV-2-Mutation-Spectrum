#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2     8       G       3       ^~.^~.^~.       III    
my $file3=shift;
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	chomp $line;
	my @name=split /\t/,$line;
	my %mutation=();
	my $mis=0;my $same=0;my $site=();
	my $same_read1=();
	
	$name[4]=~s/\^.//g;
	$name[4]=~s/\$//g;
	while ($name[4]=~/[A|T|C|G]{1}/gi) {#|\d
		$mis+=1;
		$mutation{uc($&)}+=1;#SNP count
	}
	while ($name[4]=~/[\.|\,]{1}/gi) {#|\d
		$same+=1;
	}
	$name[2]=uc($name[2]);
		
	if ($same==0) {
		print "$name[0]\t$name[1]\t".($same+$mis)."\t$name[2]\t$name[2]:0:0";
		#         chr        Pos        Number       Ref            Ref£ºNo.£ºFraction£ºReads   
	}else{print "$name[0]\t$name[1]\t".($same+$mis)."\t$name[2]\t$name[2]:$same:".($same/($same+$mis));}
	foreach my $key(sort keys %mutation) {
		#        SNP£ºNo.£ºFraction£ºReads
		print "\t$key:$mutation{$key}:".($mutation{$key}/($same+$mis));
	}
	print "\n";

}
close IN1;