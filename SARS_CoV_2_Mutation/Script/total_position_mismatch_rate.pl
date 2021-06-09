#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2	1044	176349	C	A:176:0.000998020969781513	G:48:0.00027218753721314	T:71:0.000402610732127769
my $file3=shift;#/Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/STAR/Junction_read/20200531/SCV2.Unique.summarize
my $total=0;
my %mutation;
my %base;
my $mis;
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	unless ($line=~/Position/) {
		chomp $line;
		my @mid=split /\t/,$line;
		if ($mid[1] != 28144 && $mid[1] !=8782 && $mid[1] !=5062 && $mid[1] !=4402) {
			$total=$mid[2];
			$mis=0;
			#$base{$mid[3]}+=$mid[2];
			for (my $i=4;$i<=$#mid;$i++) {
				my @frac=split /:/,$mid[$i];
					$mis+=$frac[1];
					$mutation{$mid[3]}{$frac[0]}+=$frac[1];
			}
		}
		print "$mid[1]\t$mis\t$total\n";
	}
}
close IN1;



#print "Error_Type\tMismatch_read_Number\tTotal_read_Number\tMismatch_read_Number/Total_read_Number\n";

#foreach my $ref (sort keys %mutation) {#首先对key1进行排序
    
#     foreach my $mut (keys %{$mutation{$ref}}) {

#        print "$ref > $mut"."\t".$mutation{$ref}{$mut}."\t$base{$ref}\t".($mutation{$ref}{$mut}/$base{$ref})."\n";
#    }
#}
#print "Sum\t$mis\t$total\t".($mis/$total)."\n";
