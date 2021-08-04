#!usr/bin/perl
use strict;
my $file=shift;
open INPUT,"$file";
open OUTPUT,">result.txt";
my(@input,$input,$output,$len,$star);
@input=<INPUT>;
foreach $input(@input){
   chomp($input);
if($input=~/>/){
   print OUTPUT "$input\n";
}else{
$len=length($input);
      for($star=0;$star<$len;$star+=70){
          $output=substr($input,$star,70);
          print OUTPUT "$output\n";
          } 
     }
}