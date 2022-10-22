#!/usr/bin/perl
use strict;

my $inputfile=shift;
open SEQ, $inputfile || die "$!";
my $count = 0;
my $n = 0;
while(<SEQ>) {
   chomp;
   next if(/^#/);
   next if(/^$/);
   my @temp = split(/\t/);
   if($temp[1]!~/^\d+$/ or @temp < 3) {
     $n++;
   }
   $count++;
   if($count > 2000) {
      last;
   }
}
close SEQ;
print "Format error in the upload sequences. Only files in VCF format are supported." if $n != 0;
