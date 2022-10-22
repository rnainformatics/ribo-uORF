#!/usr/bin/perl -w

use strict;

my $readfile = $ARGV[0];
my $cleanFaFile = $ARGV[1];
my $outdir = $ARGV[2];
 
my %reads;
open(IN, "samtools view -F 4 $readfile |");
while (<IN>) {
	chomp;
	my @tem = split(/\t/,$_);
	$reads{$tem[0]} = 1;
}
close IN;


open CHART, "$cleanFaFile" ||die $!;
open OUT,">$outdir/temp/clean.uniqueMapping.fa" or die $!;
$/=">";<CHART>;$/="\n";
     while(<CHART>){
       chomp;
       my $seqname = $_;
       $/=">";
       my $seq=<CHART>;
       chomp($seq);
	   $seq =~ s/\n//;
       if(exists $reads{$seqname}){
           print OUT ">$seqname\n$seq\n";
       }
       $/="\n";
}
close CHART; 
close OUT;