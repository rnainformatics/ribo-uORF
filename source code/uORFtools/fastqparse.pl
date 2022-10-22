#!/usr/bin/perl
use strict;

use File::Basename;
use List::MoreUtils qw( minmax );
use LWP::Simple;
use Getopt::Std;
use FindBin qw($Bin);
use vars qw($opt_i $opt_d $opt_p $opt_g $opt_m $opt_n $opt_a);
getopts('i:d:p:m:n:a:');

my $fileid       = $opt_i;
my $jobid        = $opt_d;
my $species      = $opt_p;
my $minlength    = $opt_m;
my $maxlength    = $opt_n;
my $adapterParameters = $opt_a;


my ($adapter3, $adapter5, $adaptermismatch, $adapterminlen, $rm5p, $rm3p) = split(/:/,$adapterParameters);

my $uploadPath = "./data/uploads/fq/$jobid";
my $resultPath = "./data/uploads/fq/";
if(not -d $uploadPath){
   system("mkdir $uploadPath");
}
my $logPath = "./data/results/$jobid";

unless($opt_i or $opt_d) {
    die "Usage perl fastqparse.pl <id> <jobid> <species> <minlen> <maxlen>\n";
}

my $file1;
my $objsrrid;
my $intype;
if($fileid =~ /^\d+$/){
    $file1 = $resultPath."/".$fileid.".filter.fq";
	$intype = "fastq";
} elsif($fileid =~ /\/.*$/i) {
	$intype = "link";
	my $filename = $fileid;
	$filename =~ s/.*\///;
	if($fileid =~ /bioinformatics\.sc\.cn/) {
		get_log_status("Uploading data","Uploading");
	} else {
		get_log_status("Downloading from link from $fileid","Downloading");
	}
	print("getstore($fileid, $uploadPath/$filename)\n");
	my $rc = getstore($fileid, "$uploadPath/$filename");
	if (is_error($rc)) {
	  system("wget $fileid -P $uploadPath");
	} 
	if(not -e "$uploadPath/$filename")  {
		get_log_status("Cannot download the link $fileid, please check the link.","Error");
		die "Download <$fileid> failed";
	}
	   
    if($filename =~/.gz$/) {
		system("zcat  $uploadPath/$filename > $resultPath/$jobid.upload.tem.fa");
    } elsif($filename =~/.zip$/){
		system ("unzip -p $uploadPath/$filename > $resultPath/$jobid.upload.tem.fa");
    } else {
		system("mv $uploadPath/$filename $resultPath/$jobid.upload.tem.fa");
    }
	checkfa($minlength,$maxlength,"$resultPath/$jobid.upload.tem.fa","$resultPath/$jobid.upload.fa");
	unlink("$uploadPath/$filename");
	exit;
}

get_log_status("Quality type checking","Processing");
my $qualtype=&qualityCheck($file1);

my $ada3string = "";
if($adapter3 ne "NA") {
	$ada3string .= "-a $adapter3";
}
#if($adapter5 ne "NA") {
#	$ada3string .= " -a2 $adapter5";
#}

my $cleanedFastqFile = "";
if($ada3string ne "") {
	if($rm5p > 0) {
		$ada3string .= " --clip_R1 $rm5p";
	}
	if($rm3p > 0) {
		$ada3string .= " --three_prime_clip_R1 $rm3p";
	}
	if($adapterminlen eq "NA") {
		$adapterminlen = 6;
	}
	if($adaptermismatch eq "NA") {
		$adaptermismatch = 0.1;
	}
	my $qualityFilterCommand = "trim_galore $file1 --phred$qualtype $ada3string -length $minlength --max_length $maxlength -e $adaptermismatch -j 4 --max_n 1 -q 20 --stringency $adapterminlen -o $uploadPath > $uploadPath/$jobid.report.txt";		
	get_log_status("Adapter cleaning and low quality filtering","Processing");
	system($qualityFilterCommand);
	my $filename = basename($file1);
	$filename =~ s/\.fastq//;
	$cleanedFastqFile = "$uploadPath/$filename"."_trimmed.fq";
} else {
	$cleanedFastqFile = "$file1";
}


open CHART,  $cleanedFastqFile||die $!;
my %seqs;
my $seqNum=0;
while(<CHART>) {
   my $seqname = $_;
   my $seq = <CHART>;
   my $dire = <CHART>;
   my $qual=<CHART>;
   chomp($seq);
   $seq=~s/\s+//g;
   $seq=~s/\n+//g;
   my $polyA=0;
   if($seq =~ /(A{6,})/i){
       $polyA = 1 if(length($1) > 8);
   }
   if($polyA == 0) {
     if(not exists $seqs{$seq}) {
       $seqs{$seq} = 1;
     } else {
       $seqs{$seq}++;
     }
   }
   $seqNum++;
}
close CHART;


open OUT, ">$resultPath/$jobid.upload.fa" || die $!;
my $n=0;
my %lengthStat;
for my $seq (reverse sort {$seqs{$a}<=>$seqs{$b}} keys %seqs) {
	$n++;
	if(not exists $lengthStat{length($seq)}) {
		$lengthStat{length($seq)} = 1;
	} else {
		$lengthStat{length($seq)}++;
	}
	print OUT ">seq$n"."_x$seqs{$seq}\n$seq\n";
}
close OUT;
#system("rm $uploadPath/$jobid.ada.fastq");
system("rm $cleanedFastqFile");
#if(-e $file1){
#	system("rm $file1");
#}

if (scalar(keys %lengthStat) < 5) {
	  get_log_status("Adaptor information may be wrong! Please check!","Error");
	  system("rm $resultPath/$jobid.upload.fa");
	  die "Adaptor information is wrong!";
}

if ($seqNum < 100000) {
	  get_log_status("There are less than 100k clean sequences which are in the defined length range. Please check and make sure you have inputted the correct adapter information","Error");
	  unlink("$resultPath/$jobid.upload.fa");
	  die "Not enough clean sequences!";
}

sub qualityCheck {
    my $inputfile = shift;
	if ($inputfile =~ /.gz$/) {
		open(CHART, "gzip -dc $inputfile |") || die "cannot open $inputfile\n";
	} else {
		open(CHART, "<$inputfile") || die "cannot open $inputfile\n";
	}
	my $cnt=0;
	my ($min, $max); # global min and max values
	my $limit = 1000;
	while (my $id = <CHART>) {
			$id =~ m/^@/ || die "expected @ not found in line 1!\n";
			my $seq = <CHART>;
			my $sep = <CHART>;
			$sep =~ m/^\+/ || die "expected + not found in line 3!\n";
			my $qual = <CHART>;
			chomp($qual);
			$cnt++;
			$cnt>=$limit && last;
			my @chars = split("", $qual);
			my @nums = sort { $a <=> $b } (map { unpack("C*", $_ )} @chars);
			if ($cnt==1) {
					($min, $max) = minmax @nums;
			} else {
					my ($lmin, $lmax) = minmax @nums; # local values for this read
					$lmin<$min ? $min=$lmin : $min=$min;
					$lmax>$max ? $max=$lmax : $max=$max;
			}
	}
	close CHART;

	my %diag;
	my %comment=(
					'Sanger' => 'Phred+33,  Q[33; 73],  (0, 40)',
					'Solexa' => 'Solexa+64, Q[59; 104], (-5, 40)',
					'Illumina 1.3+' => 'Phred+64,  Q[64; 104], (0, 40)',
					'Illumina 1.5+' => 'Phred+64,  Q[66; 104], (3, 40), with 0=N/A, 1=N/A, 2=Read Segment Quality Control Indicator',
					'Illumina 1.8+' => 'Phred+33,  Q[33; 74],  (0, 41)',
	);

	if ($min < 33) { 
		warn "No known encodings with chars < 33";
		$diag{'33'}='x';
	} elsif ($min < 64) {
		$diag{'33'}='x';
	} elsif ($min <= 126) {
		$diag{'64'}='x';
	} else {
		warn "No known encodings with chars > 126";
		$diag{'64'}='x';
	}
	
	my $qualtype;
	my @types = keys %diag;
	if(@types == 2) {
	   print join("\t",@types),"\n";
	   $qualtype = "33"; 
	   print "Quality values warning. found two type of quality\n";
	} else {
		$qualtype = $types[0];
	}
	print "Phred".$qualtype,"\n";
	return $qualtype;
}


sub get_log_status()
{
    my $info = $_[0];
	my $status = $_[1];
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon +=1;
	my $logdate=sprintf("$info\t%04d-%02d-%02d %02d:%02d:%02d\n",$year,$mon,$mday,$hour,$min,$sec);
	#print($logdate);
	system "echo '{\"type\":\"singlecase\",\"species\":\"$species\",\"intype\":\"$intype\",\"current\":\"$status\",\"message\":\"$info\",\"percent\":\"3\"}' > $logPath/status.json";
	system "echo '$logdate' >> $logPath/log.txt";
}

sub checkfa()
{
	my $minlength=shift;
	my $maxlength=shift;
	my $inputfile=shift;
	my $uploadfile=shift;
	open SEQ, $inputfile || die "$!";
	open OUT, ">$uploadfile" || die "$!";
	my $n = 0;
	$/=">";<SEQ>;$/="\n";
	while(<SEQ>) {
	   chomp;
	   my $seqhead = $_;
	   $/=">";
	   my $seq=<SEQ>;
	   chomp($seq);
	   $seq =~ s/\s+$//;
	   if(length($seq) >= $minlength and length($seq) <= $maxlength) {
		   if($seqhead =~ /(\S+?_x\d+)/) {
			  print OUT ">$seqhead\n$seq\n";
			  $n++;
		   }
	   }
	   $/="\n";
	}
	close SEQ;
	close OUT;
	unlink($inputfile);
	if($n == 0){
	   &get_log_status("When link is provided, only file in collapsed fasta format is supported (fasta head like >A<span style='color:red'>_x</span>B, A is seq ID while B is frequency)","Error");
	   unlink($uploadfile);
	} elsif($n < 10000){
	   &get_log_status("There are less than 10k RPF sequences. Please check and make sure you have enough data","Error");
	   unlink($uploadfile);
	} else {
		if($fileid =~ /bioinformatics\.sc\.cn/) {
			&get_log_status("Data uploaded","Uploading");
		} else {
			&get_log_status("Data downloaded","Downloading");
		}	   
	}
}

