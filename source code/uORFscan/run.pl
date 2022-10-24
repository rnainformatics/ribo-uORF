#!/usr/bin/perl -w

use JSON::XS;
use Data::Dumper;

my ($species, $minlen, $maxlen, $mismatch, $maxMultiMapping, $detectOffset, $deDup, $ORFpvalue, $jobid, $uploadfile, $ncpu);
($species, $minlen, $maxlen, $mismatch, $maxMultiMapping, $detectOffset, $deDup,  $ORFpvalue, $jobid, $uploadfile)=@ARGV;

my $samplenames = "NA";

my $dbbasepath =   "./db/";
my $genomeindex = "$dbbasepath/STAR_index/$species";
my $genomefa = "$dbbasepath/genome/$species.fa";
my $gtfindex = "$dbbasepath/annotation/$species.annotation.gtf";
my $ribocodeLib = "$dbbasepath/ribocode/$species";
my $genepredindex = "$dbbasepath/ribORF_uORF/$species";
my $mRNAindex = "$dbbasepath/mRNA/$species.txdb.fa";
my $rRNAindex = "$dbbasepath/rRNA/$species.rRNA.fa";
my $tRNAindex = "$dbbasepath/tRNA/$species.tRNA.fa";
my $snRNAindex = "$dbbasepath/snRNA/$species.snRNA.fa";

$minlen //= 26;
$maxlen //= 34;
$mismatch //= 2;
$detectOffset //= 1;
$deDup //= "N";
$ncpu //= 5;

if(not -e $uploadfile) {
	print "No input data!";
	exit;
}
system  "mkdir -p ./data/results/$jobid/temp" if ! -d "./data/results/$jobid/temp";
$resultdir = "./data/results/$jobid";
$resultTemp = "./data/results/$jobid/temp";

if($samplenames eq "NA") {
   $samplenames = $jobid;
}

#system "echo '' > $resultdir/log.txt";

&get_log_status("Data importing and pre-processing","data import",5);
system("perl ./program/filterLength.pl $minlen $maxlen $uploadfile $jobid $deDup");
&get_log_status("RPF mapping to tRNA and rRNA using Bowtie v1.2.2","RPF mapping",10);
system("bowtie -f -p $ncpu -v $mismatch --un $resultTemp/norRNA.fa $rRNAindex $resultTemp/filter.fa > $resultTemp/rRNA.txt");
system("bowtie -f -p $ncpu -v $mismatch --un $resultTemp/norRNAtRNA.fa $tRNAindex $resultTemp/norRNA.fa > $resultTemp/tRNA.txt");
if(-e $snRNAindex) {
	system("bowtie -f -p $ncpu -v $mismatch --un $resultTemp/norRNAtRNAsnRNA.fa $snRNAindex $resultTemp/norRNAtRNA.fa > $resultTemp/snRNA.txt");
}

##filter
system("perl ./program/filter_stat.pl $resultdir");

my $nseq = 0;
open SEQ, "$resultTemp/clean.fa" || die "$!";
$/=">";<SEQ>;$/="\n";
while(<SEQ>) {
   chomp;
   my $seqhead = $_;
   $/=">";
   my $seq=<SEQ>;
   $nseq++;
   $/="\n";
}
close SEQ;
if($nseq < 10000){
   &get_log_status("There are less than 10k clean RPF sequences. Please check and make sure you have enough data","Error",10);
   exit;
}

&get_log_status("RPF mapping to genome using STAR v2.7.3a","RPF mapping",10);
system("STAR --outFilterType BySJout --runThreadN $ncpu --outFilterMismatchNmax $mismatch --genomeDir $genomeindex --readFilesIn $resultTemp/clean.fa  --quantMode TranscriptomeSAM GeneCounts --outSAMattributes MD NH --alignEndsType EndToEnd --outFilterMultimapNmax $maxMultiMapping --outFileNamePrefix $resultdir/genome --outSAMtype BAM SortedByCoordinate");

&get_log_status("Extract end to end mapped reads","RPF mapping",20);
system("perl  program/bamExpander.pl $resultdir $deDup|samtools view -bS -o $resultdir/genomeFqAligned.sortedByCoord.out.bam -");
system("samtools index $resultdir/genomeFqAligned.sortedByCoord.out.bam");

#system("ribotish quality -b $resultdir/genomeFqAligned.sortedByCoord.out.bam -g $gtfindex -o $resultdir/ribotish_qual -f $resultdir/ribotish_qual.pdf -p $ncpu -l 26,34");

#count expression RPF count
#&get_log_status("RPF counting and expression profiling using featureCounts in Subread package v1.6.3 ","Ribosome profiling",25);
#system("featureCounts -T $ncpu -t exon -g gene_id -a $gtfindex -R CORE --Rpath $resultdir -o $resultdir/featurecount.unique.txt $resultdir/genomeAligned.sortedByCoord.out.bam");

system("samtools view -F4 $resultdir/genomeAligned.sortedByCoord.out.bam|awk '{print \$1}'|uniq > $resultTemp/uniqueMap.read.txt");

system("perl ./program/parseFeatures.pl $resultdir $resultTemp/clean.fa $deDup");
open MAPRATE, "$resultdir/mapRedupStat.txt" || die "$!";
my $mappingRate;
if($deDup ne "N") {
	my ($du_mapped, $non_du_mapped, $non_mapped);
	while(<MAPRATE>){
		chomp;
		my @tem =split(/\t/);
		if($tem[0] eq "Duplicates_mapped") {
			$du_mapped = $tem[2];
		} elsif($tem[0] eq "De_duplicates_mapped") {
			$non_du_mapped = $tem[2];
		} elsif($tem[0] eq "Non_mapped") {
			$non_mapped = $tem[2];
		}
	}
	$mappingRate = sprintf("%.2f",(($du_mapped+$non_du_mapped)/($non_mapped+$du_mapped+$non_du_mapped))*100);
} else {
	my ($mapped, $non_mapped);
	while(<MAPRATE>){
		chomp;
		my @tem =split(/\t/);
		if($tem[0] eq "Mapped") {
			$mapped = $tem[2];
		}elsif($tem[0] eq "Non_mapped") {
			$non_mapped = $tem[2];
		}
	}
	$mappingRate = sprintf("%.2f",($mapped/($mapped+$non_mapped))*100);
}
close MAPRATE;
if($mappingRate < 5) {
	&get_log_status("The genome mapping rate is only <font color='red'>$mappingRate%</font>! Please check! Maybe you selected the wrong species.","Error",25);
	exit;
}
#system("perl ./program/filterUniqueFa.pl $resultdir/genomeAligned.sortedByCoord.out.bam $resultTemp/clean.fa $resultdir");


&get_log_status("RPF mapping to transcriptome using using Bowtie v1.2.2","RPF profiling",30);
system("bowtie2 -f -p $ncpu -a -x $mRNAindex -U $resultTemp/clean.uniqueMapping.fa | samtools view -bS > $resultTemp/mRNA.bam"); 
system("samtools sort -o $resultdir/mRNA.sort.bam $resultTemp/mRNA.bam");
#system("samtools view -F4 $resultTemp/mRNA.bam|awk '{print \$1}'|uniq > $resultTemp/mRNA.read.txt");
system("samtools index $resultdir/mRNA.sort.bam");

##p-site infer
my $offsetString = "26,27,28,29,30,31 12,12,12,12,12,12";
if($detectOffset == 1) {
	&get_log_status("Inferring P-site using plastid v0.4.8","RPF profiling",30);
	my $psiteInferTool = "plastid";
	if($psiteInferTool eq "ribotish") {
	  #system("/share/software/Python-2.7.15/bin/ribotish quality -b $resultdir/genomeFqAligned.sortedByCoord.out.bam -g $gtfindex -o $resultTemp/qual -p 12 -l $minlen,$maxlen"); 
	  open IN, "$resultdir/genomeFqAligned.sortedByCoord.out.bam.para.py" || die "$!";
	  while(<IN>) {
		chomp;
		s/offdict = //;
		$offsetJson = $_;
	  }
	  close IN;
	  $offsetJson =~ s/(\d+):/'$1':/g;
	  $offsetJson =~ s/'/"/g;
	  my $offset= decode_json($offsetJson);
	  open OUT, ">$resultdir/offsets.conf.txt" || die "$!";
	  print OUT "length\tp_offset\n";
	  my @lengths;
	  my @offsets;
	  foreach my $key (sort keys %{$offset}) {
	    next if($key !~ /^\d+$/);
		push @lengths,$key;
		my $os = $$offset{$key};
		if($os < 10 or $os > 15) {
			$os = 12;
		}
		push @offsets,$os;
		print OUT "$key\t$os\n";
	  }
	  $offsetString = join(",",@lengths)." ".join(",",@offsets);
	  close OUT;
	} elsif($psiteInferTool eq "plastid") {
      system("psite ./db/annotation/$species.orfs_rois.txt $resultdir/psite --min_length $minlen --max_length $maxlen --require_upstream --count_files  $resultdir/genomeFqAligned.sortedByCoord.out.bam");
	  if(-e "$resultdir/psite_p_offsets.txt") {										
	  open IN, "$resultdir/psite_p_offsets.txt" || die "$!";
	  open OUT, ">$resultdir/offsets.conf.txt" || die "$!";
	  print OUT "length\tp_offset\n";
	  my @lengths;
	  my @offsets;
	  while(<IN>) {
		chomp;
		if(/^\d+/) {
			my @temp = split(/\t/);
			push @lengths,$temp[0];
			if($temp[1] < 9) {
				$temp[1] = 9;
			}
			if($temp[1] > 15) {
				$temp[1] = 15;
			}
			push @offsets,$temp[1];
			print OUT $temp[0]."\t".$temp[1]."\n";
		}
	  }
	  close IN; 
	  close OUT;
	  $offsetString = join(",",@lengths)." ".join(",",@offsets);
	  }
	}
} elsif(not -e "$resultdir/offsets.conf.txt") {
	system("cp ./db/offsets.conf.txt $resultdir/");
}
 

print STDERR "ORF finding ...\n";
&get_log_status("Active translated ORF calling","Annotating",30);

my $orftool = "ribORF";
if($orftool eq "Price") {
	system("/share/software/Price_1.0.3/gedi -e Price -reads $resultdir/genomeFqAligned.sortedByCoord.out.bam -genomic ./db/price/$species.oml -prefix $resultdir/hsv -progress -plot");
	system("perl ./program/parsePrice.pl $resultdir $species");
} elsif($orftool eq "ribotish") {
	system("ribotish predict -b $resultdir/genomeFqAligned.sortedByCoord.out.bam -g $gtfindex -f $genomefa -o $resultdir/orflist.txt -p $ncpu");
} elsif($orftool eq "ribORF") {
	#print("perl ./program/ribORF.parrel.pl -f $resultdir/genomeAligned.sortedByCoord.out.bam -s $resultdir/offsets.conf.txt -d $genepredindex -o $resultdir\n");
	system("perl ./program/ribORF.parrel.pl -f $resultdir/genomeAligned.sortedByCoord.out.bam -s $resultdir/offsets.conf.txt -d $genepredindex -o $resultdir");
}
#./data/results/Je3Re096PwUXpzY1
#system("Rscript ./program/parseORF.R $species $jobid $orftool");


my $nrecord = 0;
open SEQ, "$resultdir/repre.valid.pred.pvalue.parameters.txt" || die "$!";
<SEQ>;
while(<SEQ>) {
   chomp; 
   $nrecord++;
}
close SEQ;
if($nrecord < 2){
   &get_log_status("No uORF detected in your data","Error",70);
   exit;
}

print STDERR "Ribo metagene analysis ...\n";
&get_log_status("Ribo metagene analysis","Annotating", 70);
system("Rscript ./program/ribo-meta_web_single.R $species $jobid $samplenames");

system("Rscript ./program/get_plot_highcharts.R $jobid");
system("perl ./program/loadMoreMake.pl $jobid");
#system("chown -R www-data:www-data $resultdir");

#system("rm $resultTemp -fr");
#if(-e $jobid."_psites.hd5") {
#	system("rm $jobid"."_psites.hd5");
#}
#system("rm $resultTemp -fr");
&get_log_status("Job completed. Thanks for using Ribo-uORF!","Job complete",100);



sub get_log_status()
{
    my $info = $_[0];
	my $status = $_[1];
	my $percent = $_[2];
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon +=1;
	my $logdate=sprintf("%04d-%02d-%02d %02d:%02d:%02d\n",$year,$mon,$mday,$hour,$min,$sec);
	$logdate = "$info\t$logdate";
	print($logdate);
	system "echo '{\"type\":\"singlecase\",\"species\":\"$species\",\"sample\":\"$samplenames\",\"current\":\"$status\",\"percent\":\"$percent\",\"message\":\"$info\"}' > $resultdir/status.json";
	system "echo '$logdate' >> $resultdir/log.txt";
}
