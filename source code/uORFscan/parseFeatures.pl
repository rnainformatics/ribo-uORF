#!/usr/bin/perl -w

use strict;

my $resultdir = $ARGV[0];
my $cleanFaFile = $ARGV[1];
my $deDup = $ARGV[2];

my $countFile = "$resultdir/temp/uniqueMap.read.txt";
my $dupReadFile = "$resultdir/temp/dupReads.txt";

my %dupReads;
if(-e $dupReadFile) {
	open IN,$dupReadFile or die $!;
	while(<IN>){
		chomp;
		$dupReads{$_} =1;
	}
	close IN;
}


open IN,$countFile or die $!;
my %stat;
my %reads;
while(<IN>){
	chomp;
	next if(/^#/);
	my @tem = split(/\t/,$_);
	$reads{$tem[0]} = 1;
}
close IN;

open CHART, "$cleanFaFile" ||die $!;
open OUT,">$resultdir/temp/clean.uniqueMapping.fa" or die $!;
$/=">";<CHART>;$/="\n";

my $uniqReadNum = 0;
my $totalReadNum = 0;
my $uniqRedupMapNum = 0;
my $totalRedupMapNum = 0;
my $uniqMapNum = 0;
my $totalMapNum = 0;
my $uniqNonMapNum;
my $totalNonMapNum;
my $uniqNonRedupMapNum;
my $totalNonRedupMapNum;

while(<CHART>){
       chomp;
       my $seqname = $_;
	   my @nameTem = split(/\t/,$seqname);
	   my @nameTem2 = split(/_[xX]/,$seqname);
	   $uniqReadNum++;
	   $totalReadNum += $nameTem2[1];
       $/=">";
       my $seq=<CHART>;
       chomp($seq);
	   $seq =~ s/\n//;
	   if(exists $reads{$seqname}) {
		   $uniqMapNum++;
		   $totalMapNum += $nameTem2[1];
		   if (exists $dupReads{$seqname}) {
				$uniqRedupMapNum++;
				$totalRedupMapNum += $nameTem2[1];
		   }
		   if($deDup ne "N") {
				if (not exists $dupReads{$seqname}) {
					#my @tem = split(/_[xX]/,$seqname);
					#$tem[1] = 1;
					#$seqname = $tem[0]."_x".$tem[1];
					print OUT ">$seqname\n$seq\n";
				}
		   } else {
				print OUT ">$seqname\n$seq\n";
		   }
	   }
       $/="\n";
}
close CHART; 
close OUT;

$uniqNonMapNum = $uniqReadNum - $uniqMapNum;
$totalNonMapNum = $totalReadNum - $totalMapNum;
$uniqNonRedupMapNum = $uniqMapNum - $uniqRedupMapNum;
$totalNonRedupMapNum = $totalMapNum - $totalRedupMapNum;

open OUT, ">$resultdir/mapRedupStat.txt" ||die $!;
open CSV, ">$resultdir/mapRedupStat.csv" ||die $!;
print OUT "Type\tUnique\tTotal\n";
if($deDup ne "N") {
	print OUT "Duplicates_mapped\t$uniqRedupMapNum\t$totalRedupMapNum\n"; 
	print OUT "De_duplicates_mapped\t$uniqNonRedupMapNum\t$totalNonRedupMapNum\n"; 
	print OUT "Non_mapped\t$uniqNonMapNum\t$totalNonMapNum\n";
} else {
	print OUT "Mapped\t$uniqMapNum\t$totalMapNum\n"; 
	print OUT "Non_mapped\t$uniqNonMapNum\t$totalNonMapNum\n";
}
close OUT;
print CSV "Type,Unique,Total\n";
if($deDup ne "N") {
	print CSV "Duplicates_mapped,$uniqRedupMapNum,$totalRedupMapNum\n"; 
	print CSV "De_duplicates_mapped,$uniqNonRedupMapNum,$totalNonRedupMapNum\n"; 
	print CSV "Non_mapped,$uniqNonMapNum,$totalNonMapNum\n";
} else {
	print CSV "Mapped,$uniqMapNum,$totalMapNum\n"; 
	print CSV "Non_mapped,$uniqNonMapNum,$totalNonMapNum\n";
}
close CSV;