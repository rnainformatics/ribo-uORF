#!/usr/bin/perl -w
use strict;

my %control;
my %uniqStat;
my %totalStat;

my $resultdir = $ARGV[0];
my $file1 = "$resultdir/temp/filter.fa ";
if ($file1 =~ /.gz$/) {
    open(CHART, "gzip -dc $file1 |") || die "cannot open $file1\n";
} else {
    open(CHART, "<$file1") || die "cannot open $file1\n";
}

my @seqTypes = ("rRNA","tRNA");
my %excludeReads;
open IN, "$resultdir/temp/rRNA.txt" || die "$!";
while(<IN>){
   chomp;
   my @temp = split(/\t/,$_);
   if(not exists $control{$temp[0]}) {
	   $uniqStat{"rRNA"} += 1;
	   if($temp[0] =~ /_[xX](\d+)/){
		  $totalStat{"rRNA"} += $1;
	   }
	   $control{$temp[0]} = 1;
   }
   $excludeReads{$temp[0]} = 1;
}
close IN;

open IN, "$resultdir/temp/tRNA.txt" || die "$!";
while(<IN>){
   chomp;
   my @temp = split(/\t/,$_);
   if (not exists $control{$temp[0]}) {
	   $uniqStat{"tRNA"} += 1;
	   if($temp[0] =~ /_[xX](\d+)/){
		  $totalStat{"tRNA"} += $1;
	   }
	   $control{$temp[0]} = 1;
   }
   $excludeReads{$temp[0]} = 1;
}
close IN;

if(-e "$resultdir/temp/snRNA.txt") {
	open IN, "$resultdir/temp/snRNA.txt" || die "$!";
	while(<IN>){
	   chomp;
	   my @temp = split(/\t/,$_);
	   if (not exists $control{$temp[0]}) {
		   $uniqStat{"snRNA"} += 1;
		   if($temp[0] =~ /_[xX](\d+)/){
			  $totalStat{"snRNA"} += $1;
		   }
		   $control{$temp[0]} = 1;
	   }
	   $excludeReads{$temp[0]} = 1;
	}
	close IN;
	push @seqTypes, "snRNA";
}

open OUT, ">$resultdir/temp/clean.fa" ||die $!;
my %seqs;
$/=">";<CHART>;$/="\n";
while(<CHART>){
       chomp;
       my $seqname = $_;
       $/=">";
       my $seq=<CHART>;
       chomp($seq);
	   $seq =~ s/\n//;
       if(not exists $excludeReads{$seqname}){
		   $uniqStat{"Clean"} += 1;
		   if($seqname =~ /_[xX](\d+)/){
			  $totalStat{"Clean"} += $1;
		   }
           print OUT ">$seqname\n$seq\n";
       }
       $/="\n";
}
close CHART; 
close OUT;


push @seqTypes, "Clean";
open OUT, ">$resultdir/cleanStat.txt" ||die $!;
open CSV, ">$resultdir/cleanStat.csv" ||die $!;
print OUT "Type\tUnique\tTotal\n";
print CSV "Type,Unique,Total\n";
foreach my $t (@seqTypes){
	my $uniqCount;
	my $totalCount;
	if(exists $uniqStat{$t}) {
	   $uniqCount = $uniqStat{$t};
	} else {
	   $uniqCount = 0;
	}
	if(exists $totalStat{$t}) {
	   $totalCount = $totalStat{$t};
	} else {
	   $totalCount = 0;
	}
	print OUT "$t\t$uniqCount\t$totalCount\n";
	print CSV "$t,$uniqCount,$totalCount\n"; 
}
close OUT;
close CSV;
