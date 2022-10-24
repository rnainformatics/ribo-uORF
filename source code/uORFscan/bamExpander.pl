#!/usr/bin/perl -w

use strict;

my $inputdir = $ARGV[0];
my $dedup = $ARGV[1];

my $genomeBamfile = $inputdir."/genomeAligned.sortedByCoord.out.bam";

#open(IN, "samtools view -H $genomeBamfile |");
#while (<IN>) {
#	print $_;
#}
#close IN;

my %forDeDup;
my $prevChr="";
my $prevStart="";
my $prevLength="";
my $prevStrand="";
my $prevUMI="";
my %dupReads;
open(IN, "samtools view -h -F 4 $genomeBamfile |");
open(OUT, " > $inputdir/temp/dupReads.txt") || die "$!";
while (<IN>) {
	chomp;
	if(/^@/) {
		print $_,"\n";
		next;
	}
	my @tem = split(/\t/,$_);
	my $strand;
	my $readLength = get_length_from_CIGAR($tem[5]);
	if($tem[1] & 16) {
	   $strand = "-";
	} else {
	   $strand = "+";
	}
	
	if($dedup ne "N") {
		my @tem2 = split(/:/,$tem[0]);
				  
		if(($tem[2] eq $prevChr) and ($tem[3] == $prevStart) and ($readLength == $prevLength) and ($strand eq $prevStrand) and ($tem2[0] eq $prevUMI)) {
			print OUT $tem[0],"\n";
			next;
		}
		print $_."\n";
		$prevChr = $tem[2];
		$prevStart = $tem[3];
		$prevLength = $readLength;
		$prevStrand = $strand;
		$prevUMI = $tem2[0];
	} else {
		#next if($tem[5] =~ /\d+S/);
		my @tem2 = split(/_x/,$tem[0]);
		if($tem2[1] == 1) {
			print $_."\n";
		} else {
			my @rec = @tem;
			for(my $i=1; $i<=$tem2[1]; $i++) {
				$rec[0] = $tem[0]."_".$i;
				print join("\t",@rec)."\n";
			}
		}
	}
}
close IN;
close OUT;

sub get_length_from_CIGAR {
    my($cigar) = @_;
    # length of read is sum of M/I/S/=/X operations, per the SAM spec
    my $read_length = 0;
    while($cigar =~ /(\d+)[MIS=X]/g) {
	$read_length += $1;
    }
    return $read_length;
}