#!/usr/bin/perl -w
use strict;

my $minlen = $ARGV[0];
my $maxlen = $ARGV[1];
my $file1 = $ARGV[2];
my $jobid = $ARGV[3];

my $outdir = "./data/results/$jobid";
my $resultTemp = "./data/results/$jobid/temp";

if ($file1 =~ /.gz$/) {
    open(CHART, "gzip -dc $file1 |") || die "cannot open $file1\n";
} else {
    open(CHART, "<$file1") || die "cannot open $file1\n";
}


my %lengthstat;
my $totalStat = 0;
open OUT, ">$resultTemp/filter.fa" ||die $!;
my %seqs;
$/=">";<CHART>;$/="\n";
     while(<CHART>){
       chomp;
       my $seqname = $_;
       $/=">";
       my $seq=<CHART>;
       chomp($seq);
	   $seq =~ s/\s//g;
	   $seq =~ s/\n//g;
	   my $seqlen = length($seq);
       if($seqlen >= $minlen and $seqlen <= $maxlen){
           print OUT ">$seqname\n$seq\n";
       }
	   next if $seqlen > 50;
	   if(not exists $lengthstat{$seqlen}){
	       $lengthstat{$seqlen} = 1;
		   $totalStat++;
	   } else {
	       $lengthstat{$seqlen}++;
		   $totalStat++;
	   }
       $/="\n";
}
close CHART; 
close OUT;

open OUT, ">$outdir/raw_length_stat.txt" ||die $!;
print OUT "seqlength\tcount\n";
foreach my $key (sort keys %lengthstat) {
     print OUT $key,"\t",$lengthstat{$key}/$totalStat,"\n";
}
close OUT;