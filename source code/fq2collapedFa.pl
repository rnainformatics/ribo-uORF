#!/usr/bin/perl  -w


use strict;
use warnings;
use Getopt::Long;

my $version = "1.1";
my $options = get_options($version);


my $fastqInput = $$options{'i'};

if ($fastqInput =~ /.gz$/) {
    open(CHART, "gzip -dc $fastqInput |") || die "cannot open $fastqInput\n";
} else {
    open(CHART, "<$fastqInput") || die "cannot open $fastqInput\n";
}

open OUT, ">$$options{'o'}" ||die $!;
my %seqs;
while(<CHART>) {
   chomp;
   my $seqname = $_;
   my $seq=<CHART>;
   my $strand=<CHART>;
   my $qual=<CHART>;
   $seq=~s/\s+//g;
   $seq=~s/\n+//g;
   if(not exists $seqs{$seq}) {
     $seqs{$seq} = 1;
   } else {
     $seqs{$seq}++;
   }
}
close CHART;

my $n=0;
for my $key (reverse sort {$seqs{$a}<=>$seqs{$b}} keys %seqs) {
    $n++;
    print OUT ">seq$n"."_x$seqs{$key}\n$key\n";
}
close OUT;



sub usage_message {
    my $version_num = shift;
    my $usage_message = "\nConvert fastq to collapsed fasta:
    Usage: perl $0 -i <fq> -o <fa>
    <fq> : fastq file without adaptor (.fq .fastq).
    <fa> : output fasta file (.fa .fasta).\n";
    return $usage_message;
}

sub get_options {
    my($v_num) = shift;
    my %options = ();
    GetOptions(\%options,
	       'help',
	       'i=s',
	       'o=s'
    );
    unless(%options) {
	my $usage_message = usage_message($v_num);
	die "$usage_message\n";
    }
    return \%options;
}
