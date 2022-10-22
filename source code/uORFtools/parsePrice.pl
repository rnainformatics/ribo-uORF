#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use JSON qw(to_json);

my $rawDataDir = $ARGV[0];
my $species = $ARGV[1];

my %orfinfors;
my %samples;
my %oktemp;
open OUT, ">$rawDataDir/uorf.df.txt" || die "$!";
my $infile = "$rawDataDir/hsv.orfs.tsv";
next if not -e $infile;
my @fileinfor = stat ($infile);
my $size = $fileinfor[7];
if($size > 0) {
	open CHART,"$infile" or die $!;
	<CHART>;
	my %control;
	while(<CHART>){
			chomp;
			next if(/^#/);
			my @tem2 = split(/\t/,$_);
			my $txname = $tem2[1];
			$txname =~ s/_.*//;
			if (not defined $tem2[5]) {
				next;
			}
			next if($tem2[5] ne "uoORF" and $tem2[5] ne "uORF");
			if($tem2[2] =~ /^(\w+)([+-]):(\d+)-.*-(\d+)$/) {
					my $chr = $1;
					my $strand = $2;
					my $start = $3;
					my $end = $4;
					$start = $start+1;
					my $score = $tem2[6]."\t".$tem2[7]."\t".$tem2[8]."\t".$tem2[10];
					#print OUT "chr".$chr."\t".$start."\t".$start."\t".$tem2[1]."_start\t".$tem2[5]."\t".$strand."\t".$tem2[4]."\n";
					#print OUT "chr".$chr."\t".$end."\t".$end."\t".$tem2[1]."_end\t".$tem2[5]."\t".$strand."\t".$tem2[4]."\n";
					print OUT "chr".$chr."\t".$start."\t".$end."\t".$tem2[1]."\t".$tem2[5]."\t".$strand."\t".$tem2[4]."\t".$tem2[9]."\t".$tem2[8]."\n";
					
					#$oktemp{$tem2[1]} =  "chr".$chr."\t".$end."\t".$end."\t".$tem2[1]."\t".$tem2[5]."\t".$strand."\t".$tem2[4]."\t".$score;
			}
			if(exists $control{$txname}) {
				next;
			}else {
				$control{$txname} = 1;
			}
	}
	close CHART;
}
close OUT;


=pod
system("Rscript program/mapTotranscript.R $species $rawDataDir");
open OUT2, ">$rawDataDir/uorf.df.$species.forsql.txt" || die "$!";
open CHART,"$rawDataDir/uorf.df.rev.txt" or die $!;
<CHART>;
my %uorf2txLoc;
while(<CHART>){
	chomp;
	my @tem2 = split(/\t/,$_);
	my @tem = split(/\|/,$tem2[4]);
	next if $tem2[2]-$tem2[1] < 10;
	push @{$uorf2txLoc{$tem2[4]}},$tem2[4].":".$tem2[5].":".$tem2[1].":".$tem2[2].":".$tem2[3].":".$tem2[6];
	if(exists $oktemp{$tem2[4]}) {
		print OUT2 $oktemp{$tem2[4]}."\t".$tem2[1]."\t".$tem2[2]."\n";
	}
}
close CHART;
close OUT2;
=cut

