#!usr/bin/perl
use strict;

my $file=shift;#file name
my %split;
$/="\n>";
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	$line=~s/\>//;
	my @mid=split /\n/,$line;
	for (my $i=1;$i<=$#mid ;$i++) {
		$split{$i}.=">$mid[0]\n$mid[$i]\n";
	}
	
}
close IN;

foreach my $num(keys %split) {
	open (OUT,">$num".".fas")||die;
		print OUT  $split{$num};
	close OUT;
}