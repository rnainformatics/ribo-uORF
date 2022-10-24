use strict;

if ($#ARGV < 2) {
	print "##Annotate potential ORF"."\n";
	print "usage: perl ORFcandidate.pl -g genomeSequenceFile -t transcriptomeFile -o outputDir [-s startCodon] [-l orfLengthCutoff]"."\n";
	print "-g genomeSequenceFile: input the genome assembly file in fasta format;"."\n";
	print "-t transcriptomeFile: reference transcriptome annotation file in genePred format;"."\n";
	print "-o outputDir: output directory;"."\n";
	print "-s startCodon [optional]: start codon types, default: ATG;CTG;GTG;TTG;ACG;"."\n";
	print "-l orfLengthCutoff [optional]: cutoff of minimum candidate ORF length, default: 6."."\n";
	exit;
}

use Getopt::Std;
use Data::Dumper;

### get arguments ###
my %args; 
getopt("gtosl",\%args);
my $genome=$args{g}; 
if (! $genome) {
	print "No reference genome file"."\n";
    exit;
}
my $genepred=$args{t}; 
if (! $genepred) {
	print "No reference transcriptome annotation file"."\n";
    exit;
}
my $outdir=$args{o}; 
if (! $outdir) {
	print "No output directory"."\n";
    exit;
}
my $scodon=$args{s}; 
my $lenoff=$args{l};

#($scodon="ATG/CTG/GTG/TTG/ACG/AAG/AGG/ATC/ATA/ATT") if (!$scodon);

($scodon="ATG/CTG/GTG/TTG/ACG") if (!$scodon);

($lenoff=6) if (! $lenoff);

my @s=split /\//, $scodon;
my %sc;
for (my $i=0; $i<=$#s; $i++) {
	$sc{$s[$i]}=1;
}

### read genome file
my %seq;
my $id;
open (SE, "$genome");
while (<SE>) {
	chomp;
	if ($_ =~ /^>/) {
		my @a=split /\s+/, (substr ($_, 1));
		$id=$a[0];
	} else {
		$seq{$id}.=$_;
	}
}
close SE;

### output candidate ORF
open (IN, "$genepred");
my $line;
open (OUT2, ">$outdir/all.uORF.extented.nt.fa");
open (OUT3, ">$outdir/all.uORF.extented.aa.fa");
open (OUT5, ">$outdir/all.uORF.extented.nt.all.fa");
open (OUT6, ">$outdir/all.uORF.extented.aa.all.fa");
open (OUT7, ">$outdir/all.uORF.genepred.txt");
open (OUT9, ">$outdir/all.uORF.utr5.bed");
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my @val;
	my $out;
	my $len=0;
	my $loc1=0;
	my $loc2=0;
	my $chr = $s1[1];
	my $k=$s1[0].":".$s1[1].":".$s1[2];
	for (my $i=0; $i<=$#s2; $i++) {
		$out.=substr($seq{$s1[1]}, $s2[$i], ($s3[$i]-$s2[$i]));
		$len+=$s3[$i]-$s2[$i];
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			push @val, $j;
			if ($j == $s1[5]) {
				$loc1=$#val+1;
			}
			if ($j == $s1[6]) {
				$loc2=$#val+1;
			}
		}
		
	}
	if($loc2==0) {
	    $loc2=$#val+1;
	}
	push @val, $s3[$#s2];
	#print Dumper(@val);
	$out=uc($out);
	if ($s1[2] eq '-') {
		$out=reverse($out);
		$out =~ tr/ACGT/TGCA/;
		@val = reverse (@val);
		my $a=$loc1;
		if($len == $loc2) {
			$loc1=1;
		} else {
			$loc1=$len-$loc2+2;
		}
		$loc2=$len-$a+2;
	}
	if ($s1[5] == $s1[6]) {
		$loc1=0;
		$loc2=0;
	}
	open (OUT1, ">>$outdir/$chr.candidateORF.genepred.txt");
	my $ra=0;
	for (my $i=0; $i<=($len-2); $i++) {
		my $b=substr ($out, $i, 3);
		if (exists ($sc{$b})) {
			my $ind=0; 
			for (my $j=$i; $j<=$len && $ind==0; $j+=3) {
				my $d=substr ($out, $j, 3);
				if ($d eq 'TAG' || $d eq 'TAA' || $d eq 'TGA') {
					if (($j-$i+3) >= $lenoff) {
						my $seq=substr($out, $i, $j-$i+3);
						my $seq2 = "";
						$ra++;
						my $type;
						my $loc3=$i+1;
						my $loc4=$j+4;
						if ($loc1==0) {
							$type="noncoding";
						} elsif ($loc3==$loc1 && $loc4==$loc2) {
							$type="canonical";
						} elsif ($loc3==$loc1 && $loc4==$loc2+1 && $loc4 == $len+1) {
							$type="canonical";
						} elsif ($loc3 < $loc1 && $loc4 <= $loc1) {
							$type="uORF";
							$seq2=$seq;
						} elsif ($loc3 < $loc1 && $loc4 > $loc1 && $loc4 < $loc2) {
							if(($loc1-1-$i) % 3 == 0) {
								$seq2=substr($out, $i, $loc1-$i-1);
								$type="inframe.overlap.uORF";
								print "inframe.overlap.uORF ".$k."|".$ra."|".$len.":".($i+1).":".($j+4)."|".$type."|".$b."|".$loc1."-".$loc2."-".$loc3."-".$loc4."-".$_."\n";
							} else {
								$type="outframe.overlap.uORF";
								$seq2 = $seq;
							}
						} elsif ($loc3 >= $loc1 && $loc4 < $loc2) {
							$type="internal";
						} elsif ($loc3 > $loc1 && $loc3 < $loc2 && $loc4 > $loc2) { 
							$type="external";
						} elsif ($loc3 > $loc1 && $loc4 == $loc2) {
							$type="truncation";
						} elsif ($loc3 < $loc1 && $loc4 == $loc2) {
							$type="extension";
							if(($loc1-1-$i) % 3 == 0) {
								$seq2=substr($out, $i, $loc1-1-$i);
							} else {
								#warn $k."|".$ra."|".$len.":".($i+1).":".($j+4)."|".$type."|".$b."|".$loc1."-".$loc2."-".$loc3."-".$loc4."-".$_;
							}
						} elsif ($loc3 >= $loc2) {
							$type="polycistronic";
						} elsif ($loc3 <= $loc1 && $loc4 > $loc2) {
							$type="readthrough";
						} elsif ($loc3 == $loc1 && $loc4 != $loc2) {
							$type="seqerror";
						} else {
							$type="other";
						}
						#my $id=$k."|".$ra."|".$len.":".($i+1).":".($j+4)."|".$type."|".$b;
						my $id=$k."|".$ra."|".$len.":".($i+1).":".($j+4)."|".$type."|".$b."|".$loc1."-".$loc2."-".$loc3."-".$loc4;
						
						my $out2=$id."\t".$s1[1]."\t".$s1[2]."\t".$s1[3]."\t".$s1[4];
						my $out3 = $s1[1];
						if ($s1[2] eq '+') {
						    #print $val[$i],"\n";
							$out2.="\t".$val[$i]."\t".$val[$j+3];
						} else {
						    my $end = $val[$i];
							$out2.="\t".$val[$j+2]."\t".$end; 
						} 
						$out2 .= "\t".$s1[7]."\t".$s1[8]."\t".$s1[9];
						print OUT1 $out2."\n";
						if($type =~ /overlap.uORF/ or $type eq "extension") {
							    if ($s1[2] eq '+') {
									my $end1 = $s1[5]-1;
									$out3.="\t".$val[$i]."\t".$end1;
									#my $len2 = $j+4-$i-1
									#$out3.="\t".$s1[5]-$len2."\t".$s1[5];
								} else {
									my $start2 = $s1[6]+1;
									my $end2 = $val[$i];
									$out3.="\t".$start2."\t".$end2; 
								} 
								$out3 .= "\t".$id."\t".$s1[2]."\t".$s1[7];
								print OUT9 $out3."\n";
						} elsif($type eq "uORF") {
							    if ($s1[2] eq '+') {
									$out3.="\t".$val[$i]."\t".$val[$j+3];
								} else {
									my $end3 = $val[$i];
									$out3.="\t".$val[$j+2]."\t".$end3; 
								} 
								$out3 .= "\t".$id."\t".$s1[2]."\t".$s1[7];
								print OUT9 $out3."\n";
						}
						
						if($seq2 ne "") {
							my $aaseq = &translateToProtein($seq2);
        					print OUT5 ">".$id."\n".$seq2."\n"; 
							print OUT6 ">".$id."\n".$aaseq."\n"; 
							if(length($seq2) > 30) {
								print OUT2 ">".$id."\n".$seq2."\n"; 
								print OUT3 ">".$id."\n".$aaseq."\n"; 
							}
						}
						my ($gsi,$gei);
						my @newStarts;
						my @newEnds;
						my $gstart = $s1[5]+$i+1;
						my $gend = $s1[5]+$j+1;
						for (my $t=0; $t<=$#s2; $t++) {
							if($gstart >= $s2[$t] & $gstart < $s3[$t]) {
								#
								#$gsi = $t;
								push @newStarts,$gstart;
								push @newEnds,$s3[$t];
							}
							if($s2[$t] > $gstart & $s3[$t] < $gend) {
								push @newStarts,$s2[$t];
								push @newEnds,$s3[$t];
							}
							if($gend >= $s2[$t] & $gend < $s3[$t]) {
								push @newStarts,$s2[$t];
								push @newEnds,$gend;
							}
						}
						print OUT7 $s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$gstart."\t".$gend."\t".$gstart."\t".$gend."\t".join(",",@newStarts).",\t".join(",",@newEnds).",\n";
						
					} 
					$ind=1;					
				}
			}
		}
	}
	close OUT1;
}
close OUT2;
close OUT3;
close OUT5;
close OUT6;
close OUT7;
close OUT9;
close IN;


sub translateToProtein {
	my $inseq = shift;
	$inseq =~ tr/T/U/;

	my @codon;
	my $j = 0;
	for (my $i = 0 ; $i < length ($inseq); $i = $i+3) {
		@codon[$j] = substr( $inseq, $i, 3);
		$j++;
	}
	my(%genetic_code_table) = (
		'UUU' => 'F',    
		'UUC' => 'F',
		'UUA' => 'L',
		'UUG' => 'L',
		'CUU' => 'L',
		'CUC' => 'L',
		'CUA' => 'L',
		'CUG' => 'L',
		'AUU' => 'I',
		'AUC' => 'I',
		'AUA' => 'I',
		'AUG' => 'M',
		'GUU' => 'V',
		'GUC' => 'V',
		'GUA' => 'V',
		'GUG' => 'V',
		'UCU' => 'S',
		'UCC' => 'S',
		'UCA' => 'S',
		'UCG' => 'S',
		'CCU' => 'P',
		'CCC' => 'P',
		'CCA' => 'P',
		'CCG' => 'P',
		'ACU' => 'T',
		'ACC' => 'T',
		'ACA' => 'T',
		'ACG' => 'T',
		'GCU' => 'A',
		'GCC' => 'A',
		'GCA' => 'A',
		'GCG' => 'A',
		'UAU' => 'Y',
		'UAC' => 'Y',
		'UAA' => '',
		'UAG' => '',
		'CAU' => 'H',
		'CAC' => 'H',
		'CAA' => 'Q',
		'CAG' => 'Q',
		'AAU' => 'N',
		'AAC' => 'N',
		'AAA' => 'K',
		'AAG' => 'K',
		'GAU' => 'D',
		'GAC' => 'D',
		'GAA' => 'E',
		'GAG' => 'E',
		'UGU' => 'C',
		'UGC' => 'C',
		'UGA' => '',
		'UGG' => 'W',
		'CGU' => 'R',
		'CGC' => 'R',
		'CGA' => 'R',
		'CGG' => 'R',
		'AGU' => 'S',
		'AGC' => 'S',
		'AGA' => 'R',
		'AGG' => 'R',
		'GGU' => 'G',
		'GGC' => 'G',
		'GGA' => 'G',
		'GGG' => 'G',
	);
	my @result;
	for (my $j = 0 ; $j < length ($inseq)/3; $j++) {
		if( exists ($genetic_code_table{@codon[$j]}) ) {
			@result[$j] = $genetic_code_table{@codon[$j]};
		}
	}
	return(join("",@result));
}
