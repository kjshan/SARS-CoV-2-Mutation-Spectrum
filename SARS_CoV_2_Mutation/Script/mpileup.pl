#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#NC_045512.2     4       A       718     ....T..................T..........C.....C......T..T......................
my $file3=shift;#mpileup
my %mutation;
my $num;
my $mid;

print "Chr\tPosition\tRead\tRef\tAlt1\tAlt2\tAlt3\n";
open(IN1, $file3) or die ("can not open $file3\n");
while (my $line=<IN1>) {
	#if ($line=~/^NC_045512\.2/) {
		chomp $line;
		my $mis=0;my $same=0;
		my @name=split /\t/,$line;		
		$name[4]=~s/\^.//g;
		$name[4]=~s/\$//g;
		my @mismatch=split //,$name[4];
		for (my $i=0; $i<=$#mismatch;$i++) {
			if ($mismatch[$i]=~/A|T|C|G/i) {
				$mutation{uc($mismatch[$i])}+=1;
				$mis+=1;
			}elsif ($mismatch[$i]=~/\+|\-/) {
					for (my $j=$i+1;$j<=$#mismatch;$j++) {
						if ($mismatch[$j]=~/[0-9]/) {
							$num+=1;
						}else{$j+=$#mismatch;}
					}#print "\n$i\t$num\t$mismatch[$i]\t$mismatch[$i+1]\t$mismatch[$i+2]\t$mismatch[$i+3]\n";
					for (my $f=$num;$f > 0;$f--) {
						$mid+=$mismatch[$i+$num-$f+1] * (10 ** ($f-1));
						#if ($mismatch[$i+$num-$f+1] !~/0-9/) {
								#print "\n$mismatch[$i]\t$mismatch[$i+1]\t$mismatch[$i+2]\t$mismatch[$i+3]\n";
						#}
						#print "\n$f\t$mismatch[$i+$num-$f+1]\t".(10 ** ($f-1))."\t$mid\n";
					}
				    $i+=($mid+$num);
					$num=0;$mid=0;
			}elsif ($mismatch[$i]=~/\.|\,/) {
				$same+=1;
			}
		}
		if (($same+$mis) > 0) {
			print "$name[0]\t$name[1]\t".($same+$mis)."\t$name[2]";
			foreach my $key(sort keys %mutation) {
				print "\t$key:$mutation{$key}:".($mutation{$key}/($same+$mis));
			}
		}elsif (($same+$mis) == 0) {
			print "$name[0]\t$name[1]\t0\t$name[2]";
			
		}
		print "\n";
	#}
	%mutation = ();
}
close IN1;






